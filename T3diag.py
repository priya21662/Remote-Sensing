from osgeo import gdal
import matplotlib.pyplot as plt
import numpy as np

# File paths
file_paths = {
    "T11": r"C:\Users\supri\OneDrive\Desktop\RS Internship\sample_data\full_pol\T3\T11.bin",
    "T22": r"C:\Users\supri\OneDrive\Desktop\RS Internship\sample_data\full_pol\T3\T22.bin",
    "T33": r"C:\Users\supri\OneDrive\Desktop\RS Internship\sample_data\full_pol\T3\T33.bin"
}

# Create figure with 3 subplots
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# Loop through files
for i, (title, file_path) in enumerate(file_paths.items()):
    
    dataset = gdal.Open(file_path, gdal.GA_ReadOnly)

    if dataset is None:
        print(f"Could not open {file_path}")
        continue

    # Read first band
    band = dataset.GetRasterBand(1)
    image_array = band.ReadAsArray()

    # Display image
    im = axes[i].imshow(image_array, cmap='gray')
    axes[i].set_title(title)
    axes[i].set_xlabel("Columns")
    axes[i].set_ylabel("Rows")

    # Add colorbar for each subplot
    plt.colorbar(im, ax=axes[i], fraction=0.046, pad=0.04)

    dataset = None

# Adjust spacing
plt.tight_layout()

# Show all images
plt.show()