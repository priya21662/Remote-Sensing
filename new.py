from osgeo import gdal
import matplotlib.pyplot as plt
import numpy as np

# 1. Open the file
# Note: For raw binary files without headers, you may need to use NumPy
file_path =r"C:\Users\supri\OneDrive\Desktop\RS Internship\sample_data\full_pol\T3\T33.bin"
dataset = gdal.Open(file_path, gdal.GA_ReadOnly)


# 2. Extract a specific band (starting from index 1)
band = dataset.GetRasterBand(1)

# 3. Read the band data as a NumPy array
image_array = band.ReadAsArray()

# 5. Display the image using Matplotlib
plt.figure(figsize=(10, 8))
plt.imshow(image_array, cmap='gray')  # Use 'gray' for single-band images
plt.colorbar(label='Pixel Value')
plt.title(f"Binary File Visualization: {file_path}")
plt.show()

# Clean up
dataset = None
