# Extracting bands 1 and 2 for one file using OSGeo4W
gdal_translate "HDF4_EOS:EOS_GRID:MOD09GA.hdf:MODIS_Grid_500m_2D:sur_refl_b01_1" band1.tif
gdal_translate "HDF4_EOS:EOS_GRID:MOD09GA.hdf:MODIS_Grid_500m_2D:sur_refl_b02_1" band2.tif

# Extracting bands 1 and 2 for all files in a directory -----------------------------------------------------------------------------------------------------------
import os
import subprocess

# Directory where the input HDF files are located
input_dir = 'D:/Module11/PySEBAL_data/Satellite_data'  # you can change this to your folder path
# Directory where the output TIFF files will be saved
output_dir = 'osgeo4w'

# Make sure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Band info to extract
bands = {
    'band1': 'sur_refl_b01_1',
    'band2': 'sur_refl_b02_1'
}

# Find all HDF files containing "MOD09GA"
for filename in os.listdir(input_dir):
    if filename.endswith('.hdf') and 'MOD09GA' in filename:
        base_name = os.path.splitext(filename)[0]
        input_path = os.path.join(input_dir, filename)
        
        for band_key, band_name in bands.items():
            output_filename = f"{base_name}_{band_key}.tif"
            output_path = os.path.join(output_dir, output_filename)
            
            gdal_path = "gdal_translate"  # Make sure gdal_translate is in your PATH
            # Construct the GDAL subdataset string
            subdataset = f'HDF4_EOS:EOS_GRID:{input_path}:MODIS_Grid_500m_2D:{band_name}'
            
            cmd = [gdal_path, subdataset, output_path]
            print("Running command:", ' '.join(cmd))
            subprocess.run(cmd, check=True)

print("Done extracting bands.")

# Calculating NDVI ------------------------------------------------------------------------------------------------------------------------------------------------
import os
import numpy as np
import rasterio
import matplotlib.pyplot as plt

def calculate_ndvi(red_band_path, nir_band_path, output_path):
    """
    Calculate NDVI from Red and NIR bands, scale to [0,1], and save as float32 GeoTIFF.
    """
    with rasterio.open(red_band_path) as red_src:
        red = red_src.read(1).astype('float32')
        profile = red_src.profile.copy()

    with rasterio.open(nir_band_path) as nir_src:
        nir = nir_src.read(1).astype('float32')

        if red_src.crs != nir_src.crs:
            print(f"⚠️ CRS mismatch between {os.path.basename(red_band_path)} and {os.path.basename(nir_band_path)}.")
        if red_src.transform != nir_src.transform:
            print(f"⚠️ Transform mismatch between {os.path.basename(red_band_path)} and {os.path.basename(nir_band_path)}.")

    denominator = nir + red
    valid_mask = (denominator != 0) & (~np.isinf(red)) & (~np.isinf(nir)) & (~np.isnan(red)) & (~np.isnan(nir))

    ndvi = np.full_like(red, -2.0, dtype=np.float32)
    ndvi[valid_mask] = (nir[valid_mask] - red[valid_mask]) / denominator[valid_mask]
    ndvi = np.clip(ndvi, -1.0, 1.0)
    ndvi_scaled = (ndvi + 1) / 2
    ndvi_scaled[~valid_mask] = 0.0

    profile.update(
        dtype=rasterio.float32,
        count=1,
        nodata=0.0,
        compress='lzw'
    )

    # ✅ Use correct output file path here
    with rasterio.open(output_path, 'w', **profile) as dst:
        dst.write(ndvi_scaled, 1)

    print(f"✅ NDVI saved: {os.path.basename(output_path)}")

if __name__ == "__main__":
    input_dir = 'D:/Module11/PySEBAL_data/Satellite_data/bands_1_2'
    output_dir = 'D:/Module11/PySEBAL_data/Satellite_data/ndvi'

    os.makedirs(output_dir, exist_ok=True)  # Ensure output directory exists

    for fname in os.listdir(input_dir):
        if fname.endswith('_band1.tif') and fname.startswith('MOD09GA'):
            band1_path = os.path.join(input_dir, fname)
            band2_fname = fname.replace('_band1.tif', '_band2.tif')
            band2_path = os.path.join(input_dir, band2_fname)

            if os.path.exists(band2_path):
                ndvi_fname = fname.replace('_band1.tif', '_ndvi.tif')
                ndvi_path = os.path.join(output_dir, ndvi_fname)  # ✅ Save to output folder
                print(f"🔄 Calculating NDVI for {fname}")
                calculate_ndvi(band1_path, band2_path, ndvi_path)
            else:
                print(f"❌ Band2 file not found for {fname}")
