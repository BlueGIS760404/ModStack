# Extracting bands 1 and 2 for one file using OSGeo4W
gdal_translate "HDF4_EOS:EOS_GRID:MOD09GA.hdf:MODIS_Grid_500m_2D:sur_refl_b01_1" band1.tif
gdal_translate "HDF4_EOS:EOS_GRID:MOD09GA.hdf:MODIS_Grid_500m_2D:sur_refl_b02_1" band2.tif

# Extracting bands 1 and 2 for all files in a directory -----------------------------------------------------------------------------------------------------------
import os
import subprocess

# === Configure your input and output directories ===
input_dir = r'C:\Users\Reza'  # Path to HDF files
output_dir = r'C:\Users\Reza\bands_1_2'  # Output GeoTIFFs

os.makedirs(output_dir, exist_ok=True)

# === Define the subdatasets you want to extract ===
layers = {
    'band1': ('MODIS_Grid_500m_2D', 'sur_refl_b01_1'),
    'band2': ('MODIS_Grid_500m_2D', 'sur_refl_b02_1'),
    'state_1km': ('MODIS_Grid_1km_2D', 'state_1km_1')  # 1km QA layer
}

# === Loop through all HDF files ===
for filename in os.listdir(input_dir):
    if filename.endswith('.hdf') and 'MOD09GA' in filename:
        base_name = os.path.splitext(filename)[0]
        input_path = os.path.join(input_dir, filename)
        print(f"\n📂 Processing {filename}")

        for key, (grid, subdataset_name) in layers.items():
            output_filename = f"{base_name}_{key}.tif"
            output_path = os.path.join(output_dir, output_filename)

            subdataset = f'HDF4_EOS:EOS_GRID:{input_path}:{grid}:{subdataset_name}'
            cmd = ['gdal_translate', subdataset, output_path]

            print(f"➡️  Extracting {key} → {output_filename}")
            result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

            if result.returncode != 0:
                print(f"❌ Error extracting {key}:\n{result.stderr.strip()}")
            else:
                print(f"✅ Saved: {output_filename}")

print("\n🎉 Done extracting Band 1, Band 2, and state_1km for all files.")

# Calculating NDVI ------------------------------------------------------------------------------------------------------------------------------------------------
import os
import numpy as np
import rasterio
from rasterio.enums import Resampling

def calculate_ndvi_with_relaxed_mask(red_path, nir_path, qa_path, output_path):
    """
    Calculate NDVI with relaxed masking: exclude only QA pixels with bits 0-1 == 11,
    and scale to [-2000, 10000]. Output is int16 GeoTIFF with nodata = -3000.
    """

    # Read red band
    with rasterio.open(red_path) as red_src:
        red = red_src.read(1).astype('float32')
        profile = red_src.profile.copy()
        red_shape = red.shape

    # Read NIR band
    with rasterio.open(nir_path) as nir_src:
        nir = nir_src.read(1).astype('float32')

    # Read QA band and resample to match 500m bands (state_1km is 1km resolution)
    with rasterio.open(qa_path) as qa_src:
        qa = qa_src.read(1, out_shape=red_shape, resampling=Resampling.nearest)

    # Define nodata value
    nodata_val = -3000

    # QA bitmask filtering — only exclude pixels with bits 0-1 == 11 (invalid)
    cloud_bits = qa & 0b11
    valid_qa_mask = cloud_bits != 0b11

    # NDVI calculation: (NIR - Red) / (NIR + Red)
    denominator = nir + red
    with np.errstate(divide='ignore', invalid='ignore'):
        ndvi = np.where(
            (denominator != 0) & valid_qa_mask,
            (nir - red) / denominator,
            np.nan
        )

    # Clip NDVI to [-1, 1] to avoid outliers
    ndvi = np.clip(ndvi, -1.0, 1.0)

    # Scale NDVI to [-2000, 10000] linearly
    ndvi_scaled = ((ndvi + 1) / 2) * 12000 - 2000  # full range: 12000

    # Diagnostics
    total = ndvi_scaled.size
    nodata_pixels = np.sum(np.isnan(ndvi_scaled))
    valid_pixels = total - nodata_pixels
    cloudy_pixels = np.sum(cloud_bits == 0b11)

    print(f"📊 Valid NDVI pixels: {valid_pixels}/{total} ({(valid_pixels / total) * 100:.2f}%)")
    print(f"☁️ Cloudy pixels excluded (bits 11): {cloudy_pixels} ({(cloudy_pixels / total) * 100:.2f}%)")

    # Replace NaNs with nodata and convert to int16
    ndvi_scaled = np.where(np.isnan(ndvi_scaled), nodata_val, ndvi_scaled).astype(np.int16)

    # Update profile for output GeoTIFF
    profile.update(
        dtype=rasterio.int16,
        count=1,
        nodata=nodata_val,
        compress='lzw'
    )

    # Save output
    with rasterio.open(output_path, 'w', **profile) as dst:
        dst.write(ndvi_scaled, 1)

    print(f"✅ NDVI saved: {os.path.basename(output_path)}\n")

if __name__ == "__main__":
    input_dir = r'C:\Users\Reza\bands_1_2'
    output_dir = r'C:\Users\Reza\bands_1_2\ndvi'

    os.makedirs(output_dir, exist_ok=True)

    for fname in os.listdir(input_dir):
        if fname.endswith('_band1.tif') and 'MOD09GA' in fname:
            band1_path = os.path.join(input_dir, fname)
            band2_path = os.path.join(input_dir, fname.replace('_band1.tif', '_band2.tif'))
            qa_path = os.path.join(input_dir, fname.replace('_band1.tif', '_state_1km.tif'))

            if os.path.exists(band2_path) and os.path.exists(qa_path):
                ndvi_path = os.path.join(output_dir, fname.replace('_band1.tif', '_ndvi.tif'))
                print(f"🔄 Calculating NDVI for {fname}")
                calculate_ndvi_with_relaxed_mask(band1_path, band2_path, qa_path, ndvi_path)
            else:
                print(f"❌ Missing band2 or QA for {fname}")
