from osgeo import gdal
import numpy as np
from pathlib import Path
import os 
def read_bin(file_path):
    ds=gdal.Open(file_path,gdal.GA_ReadOnly)
    band=ds.GetRasterBand(1)
    arr=band.ReadAsArray()
    ds=None
    return arr
def write_polsar_bin(data, output_path,t3_folder):
    rows, cols = data.shape
    driver = gdal.GetDriverByName("ENVI")
    ref_path = os.path.join(t3_folder, "T11.bin")
    reference_ds = gdal.Open(ref_path)
    out_ds = driver.Create(output_path, cols, rows, 1, gdal.GDT_Float32)
    out_ds.SetGeoTransform(reference_ds.GetGeoTransform())
    out_ds.SetProjection(reference_ds.GetProjection())
    band = out_ds.GetRasterBand(1)
    band.WriteArray(data)
    band.FlushCache()
    out_ds = None
    reference_ds=None
    print(f"Successfully saved: {output_path}")

def chunks(file_path):
    folder=Path(file_path)
    t11 = read_bin(os.path.join(folder,"T11.bin"))
    t22 = read_bin(os.path.join(folder,"T22.bin"))
    t33 = read_bin(os.path.join(folder,"T33.bin"))
    t12 = read_bin(os.path.join(folder,"T12_real.bin")) +1j* read_bin(os.path.join(folder,"T12_imag.bin"))
    t13 = read_bin(os.path.join(folder,"T13_real.bin")) +1j* read_bin(os.path.join(folder,"T13_imag.bin"))
    t23 = read_bin(os.path.join(folder,"T23_real.bin")) +1j* read_bin(os.path.join(folder,"T23_imag.bin"))
    rows,columns=t11.shape
    T3=np.zeros((rows,columns,3,3),dtype=np.complex64)
    T3[:, :, 0, 0] = t11
    T3[:, :, 0, 1] = t12
    T3[:, :, 0, 2] = t13

    T3[:, :, 1, 0] = np.conj(t12)
    T3[:, :, 1, 1] = t22
    T3[:, :, 1, 2] = t23

    T3[:, :, 2, 0] = np.conj(t13)
    T3[:, :, 2, 1] = np.conj(t23)
    T3[:, :, 2, 2] = t33
    return T3

def dop(arr):
    det=np.linalg.det(arr)
    trace=np.linalg.trace(arr)
    m1= np.sqrt(np.clip(1 - (27 * det) / (trace**3 + 1e-9), 0, 1))
    return m1

folder=r"C:\Users\supri\OneDrive\Desktop\RS Internship\sample_data\full_pol\T3"
t3=chunks(folder)
m=dop(t3)

output=r"C:\Users\supri\OneDrive\Desktop\RS Internship\calc_data\dop_calc.bin" 
write_polsar_bin(m,output,folder)

    

