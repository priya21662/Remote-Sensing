from osgeo import gdal
import numpy as np
def read_bin(file_path):
    data=gdal.Open(file_path,gdal.GA_ReadOnly)
    band=data.GetRasterBand(1)
    arr=band.ReadAsArray()
    return arr
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

dop_path=r"C:\Users\supri\OneDrive\Desktop\RS Internship\sample_data\full_pol\dop_act.bin"
dop_fp_path=r"C:\Users\supri\OneDrive\Desktop\RS Internship\sample_data\full_pol\T3\dop_fp.tif"
dop_cal=read_bin(dop_path)
dop_fp=read_bin(dop_fp_path)
diff=dop_fp-dop_cal
# Create a mask where the condition is true
#positive_values = diff[diff > 0]
pos_val=diff[diff>-1*np.e-6]
neg_val=diff[diff<1*np.e-6]
ref_ds=gdal.Open(r"C:\Users\supri\OneDrive\Desktop\RS Internship\sample_data\full_pol\dop_act.bin",gdal.GA_ReadOnly)
print(f"Number of positive pixels: {len(pos_val)}")
print(f"Number of negative pixels: {len(neg_val)}")
path=r"C:\Users\supri\OneDrive\Desktop\RS Internship\sample_data\full_pol\diff_dop.bin"
write_polsar_bin(diff,path,ref_ds)
ref_ds=None
data=None

