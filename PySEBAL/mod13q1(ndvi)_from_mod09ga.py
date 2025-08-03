# OSGeo4W
gdal_translate "HDF4_EOS:EOS_GRID:MOD09GA.hdf:MODIS_Grid_500m_2D:sur_refl_b01_1" band1.tif
gdal_translate "HDF4_EOS:EOS_GRID:MOD09GA.hdf:MODIS_Grid_500m_2D:sur_refl_b02_1" band2.tif


# Extract bands 1 and 2
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


