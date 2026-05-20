import numpy as np
from osgeo import gdal
import os
def read_bin(folder):
    data=gdal.Open(folder,gdal.GA_ReadOnly)
    band=data.GetRasterBand(1)
    arr=band.ReadAsArray()
    data=None
    return arr

def write_bin(data,output,t3_folder):
    rows, cols = data.shape
    driver = gdal.GetDriverByName("ENVI")
    ref_path = os.path.join(t3_folder, "T11.bin")
    reference_ds = gdal.Open(ref_path)
    out_ds = driver.Create(output, cols, rows, 1, gdal.GDT_Float32)
    out_ds.SetGeoTransform(reference_ds.GetGeoTransform())
    out_ds.SetProjection(reference_ds.GetProjection())
    band = out_ds.GetRasterBand(1)
    band.WriteArray(data)
    band.FlushCache()
    out_ds = None
    reference_ds=None
    print(f"Successfully saved: {output}")

def chunks(folder_path):
    eps=1e-10
   
    # 1. Try to find dimensions 
    t11_path = os.path.join(folder_path, 'T11.bin')
    t11=read_bin(t11_path)
    rows,cols=(t11).shape
    

    # 0:T11, 1:T12_re, 2:T12_im, 3:T13_re, 4:T13_im, 5:T22, 6:T23_re, 7:T23_im, 8:T33
    file_map = {
        0: 'T11.bin', 1: 'T12_real.bin', 2: 'T12_imag.bin',
        3: 'T13_real.bin', 4: 'T13_imag.bin', 5: 'T22.bin',
        6: 'T23_real.bin', 7: 'T23_imag.bin', 8: 'T33.bin'
    }
    data = {i: read_bin(os.path.join(folder_path,name)) for i, name in file_map.items()}
    t3=np.zeros((rows,cols,3,3),dtype=np.complex64)

    t3[:,:, 0, 0] = data[0] + eps
    t3[:,:, 1, 1] = data[5] + eps
    t3[:,:, 2, 2] = data[8] + eps

    # Upper Triangle (Complex)
    # Logic: M[0][1][0] = eps + avg[1]; M[0][1][1] = eps + avg[2];
    t3[:,:, 0, 1] = (data[1] + eps) + 1j * data[2]
    t3[:,:, 0, 2] = (data[3] + eps) + 1j * data[4]
    t3[:,:, 1, 2] = (data[6] + eps) + 1j * data[7]

    # Lower Triangle (Conjugate Transpose)
    # Logic: M[1][0][0] = M[0][1][0]; M[1][0][1] = -M[0][1][1];
    t3[:,:, 1, 0] = np.conj(t3[:,:, 0, 1])
    t3[:,:, 2, 0] = np.conj(t3[:,:, 0, 2])
    t3[:,:, 2, 1] = np.conj(t3[:,:, 1, 2])

    return t3

def reyni(T3,alpha):
    eps=1e-10
    eigvals=np.linalg.eigvalsh(T3) 
    eigvals= np.maximum(eigvals,eps)
    #eigvals=np.sort(eigvals, axis=-1)[:, :, ::-1]
    span=np.sum(eigvals,axis=-1,keepdims=True)
    p=eigvals/(eps+span)

    # Handle the alpha = 1 limit (Shannon Entropy)
    
    if np.isclose(alpha, 1.0):
        reyni_index = -np.sum(p * np.log(p), axis=-1) / np.log(3)
    else:
        
        reyni_index = np.log(np.sum(p**alpha, axis=-1)) / (np.log(3) * (1 - alpha))
    return reyni_index

folder=r"C:\Users\supri\OneDrive\Desktop\RS Internship\sample_data\full_pol\T3"
T3=chunks(folder)
reyni_index3=reyni(T3,3)
reyni_index4=reyni(T3,4)
output1=r"C:\Users\supri\OneDrive\Desktop\RS Internship\calc_data\reyni3.bin"
output2=r"C:\Users\supri\OneDrive\Desktop\RS Internship\calc_data\reyni4.bin"
write_bin(reyni_index3,output1,folder)
write_bin(reyni_index4,output2,folder)