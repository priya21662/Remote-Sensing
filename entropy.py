import numpy as np
from osgeo import gdal
import matplotlib.pyplot as plt
def read_bin(file_path):
    ds = gdal.Open(file_path)
    band = ds.GetRasterBand(1)
    return band.ReadAsArray().astype(np.complex64)
def read_complex_bin(file_path, rows, cols):
    """Reads a binary file of interleaved float32 real/imaginary values (off-diagonals like T12)."""
    # np.complex64 automatically pairs two 32-bit floats into one complex number
    data = np.fromfile(file_path, dtype=np.complex64)
    return data.reshape((rows, cols))

def write_polsar_bin(data, output_path, reference_ds):
    rows, cols = data.shape
    driver = gdal.GetDriverByName("ENVI")
    out_ds = driver.Create(output_path, cols, rows, 1, gdal.GDT_Float32)
    out_ds.SetGeoTransform(reference_ds.GetGeoTransform())
    out_ds.SetProjection(reference_ds.GetProjection())
    # Write data and clean up
    band = out_ds.GetRasterBand(1)
    band.WriteArray(data)
    band.FlushCache()
    # Setting to None closes the file properly
    out_ds = None
    print(f"Successfully saved: {output_path}")


dataset=gdal.Open(r"C:\Users\supri\OneDrive\Desktop\RS Internship\sample_data\full_pol\T3\T11.bin", gdal.GA_ReadOnly)
rows=dataset.RasterYSize
cols=dataset.RasterXSize
# Load Diagonal Elements
t11 = read_bin(r"C:\Users\supri\OneDrive\Desktop\RS Internship\sample_data\full_pol\T3\T11.bin")
t22 = read_bin(r"C:\Users\supri\OneDrive\Desktop\RS Internship\sample_data\full_pol\T3\T22.bin")
t33 = read_bin(r"C:\Users\supri\OneDrive\Desktop\RS Internship\sample_data\full_pol\T3\T33.bin")

# Load Off-Diagonal Elements (Ensure you combine Real and Imaginary)
t12 = read_bin(r"C:\Users\supri\OneDrive\Desktop\RS Internship\sample_data\full_pol\T3\T12_real.bin") + 1j * read_bin(r"C:\Users\supri\OneDrive\Desktop\RS Internship\sample_data\full_pol\T3\T12_imag.bin")
t13 = read_bin(r"C:\Users\supri\OneDrive\Desktop\RS Internship\sample_data\full_pol\T3\T13_real.bin") + 1j * read_bin(r"C:\Users\supri\OneDrive\Desktop\RS Internship\sample_data\full_pol\T3\T13_imag.bin")
t23 = read_bin(r"C:\Users\supri\OneDrive\Desktop\RS Internship\sample_data\full_pol\T3\T23_real.bin") + 1j * read_bin(r"C:\Users\supri\OneDrive\Desktop\RS Internship\sample_data\full_pol\T3\T23_imag.bin")

# Step 1: Calculate the Trace (Span)
trace = (t11 + t22 + t33).real

# Step 2: Calculate the Determinant for a 3x3 Hermitian Matrix
T3 = np.zeros((rows, cols, 3, 3), dtype=np.complex64)

T3[:, :, 0, 0] = t11
T3[:, :, 0, 1] = t12
T3[:, :, 0, 2] = t13

T3[:, :, 1, 0] = np.conj(t12)
T3[:, :, 1, 1] = t22
T3[:, :, 1, 2] = t23

T3[:, :, 2, 0] = np.conj(t13)
T3[:, :, 2, 1] = np.conj(t23)
T3[:, :, 2, 2] = t33

det_t = np.linalg.det(T3).real
# Step 3: Degree of Polarization
# We use np.clip to avoid issues with negative values due to floating point errors
dop = np.sqrt(np.clip(1 - (27 * det_t) / (trace**3 + 1e-9), 0, 1))
print(T3.shape)
path_dop=r"C:\Users\supri\OneDrive\Desktop\RS Internship\sample_data\full_pol\dop_act.bin"
write_polsar_bin(dop,path_dop,dataset)


#Eigen Values
eigvals = np.linalg.eigvals(T3)
print(eigvals.shape)
# Sort eigenvalues in descending order along the last axis
eigvals = np.sort(eigvals, axis=-1)[:, :, ::-1]

output_file = r"C:\Users\supri\OneDrive\Desktop\RS Internship\sample_data\full_pol\eigval.bin"
driver = gdal.GetDriverByName("ENVI")
# Read metadata
geotransform = dataset.GetGeoTransform()
projection = dataset.GetProjection()
# 1. Change the number of bands from 1 to 3
out_ds = driver.Create(output_file, cols, rows, 3, gdal.GDT_Float32)
# Copy spatial metadata
out_ds.SetGeoTransform(geotransform)
out_ds.SetProjection(projection)
# 2. Write each eigenvalue to its own band
for i in range(3):
    band = out_ds.GetRasterBand(i + 1)
    band.WriteArray(eigvals[:, :, i])
out_ds.FlushCache()
out_ds = None


#Entropy
eigvals = np.real(eigvals)
eigvals = np.sort(eigvals, axis=-1)[..., ::-1]
total = np.sum(eigvals, axis=-1)
p1 = eigvals[:,:,0] / total
p2 = eigvals[:,:,1] / total
p3 = eigvals[:,:,2] / total
eps = 1e-9
H = -(p1*np.log(p1+eps) + 
      p2*np.log(p2+eps) + 
      p3*np.log(p3+eps)) / np.log(3)
path=r"C:\Users\supri\OneDrive\Desktop\RS Internship\sample_data\full_pol\entropy.bin"
write_polsar_bin(H,path,dataset)
dataset=None