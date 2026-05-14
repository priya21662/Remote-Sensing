from osgeo import gdal
import matplotlib.pyplot as plt
import numpy as np
def read (path):
    ds=gdal.Open(path,gdal.GA_ReadOnly)
    band=ds.GetRasterBand(1)
    arr=band.ReadAsArray()
    return arr

#Converting bin files to numpy arrays
t11 = read(r"C:\Users\supri\OneDrive\Desktop\RS Internship\sample_data\full_pol\T3\T11.bin")
t22 = read(r"C:\Users\supri\OneDrive\Desktop\RS Internship\sample_data\full_pol\T3\T22.bin")
t33 = read(r"C:\Users\supri\OneDrive\Desktop\RS Internship\sample_data\full_pol\T3\T33.bin")
sum_arr=t11+t22+t33 #Trace of T3
print(sum_arr.shape)

# Read metadata
ds=gdal.Open(r"C:\Users\supri\OneDrive\Desktop\RS Internship\sample_data\full_pol\T3\T11.bin",gdal.GA_ReadOnly)
geotransform = ds.GetGeoTransform()
projection = ds.GetProjection()
rows, cols = sum_arr.shape
# Output file
output_file = r"C:\Users\supri\OneDrive\Desktop\RS Internship\sample_data\trace.bin"
# Create ENVI binary
driver = gdal.GetDriverByName("ENVI")
out_ds = driver.Create(output_file, cols, rows, 1, gdal.GDT_Float32)
# Copy spatial metadata
out_ds.SetGeoTransform(geotransform)
out_ds.SetProjection(projection)
# Write array
out_ds.GetRasterBand(1).WriteArray(sum_arr)
# Save
out_ds.FlushCache()
# Close files
ds=None
out_ds = None


