import numpy as np
import rasterio
from rasterio.plot import reshape_as_image
import glob
import os

# -----------------------------
# Helper: Load and average rasters
# -----------------------------
def load_and_average(file_list):
    arrays = []
    ref_raster = None

    for file in file_list:
        with rasterio.open(file) as src:
            array = src.read(1).astype(np.float32)
            array[array < 0] = np.nan  # mask NoData
            arrays.append(array)
            if ref_raster is None:
                ref_raster = src

    if not arrays:
        raise ValueError(f"No files found: {file_list}")

    return np.nanmean(arrays, axis=0), ref_raster

# -----------------------------
# Saxton & Rawls (2006) PTFs
# -----------------------------
def compute_soil_moisture_params(sand, clay, om, bd):
    om = np.clip(om, 0, 20)  # Limit organic matter
    fc = (0.2576 + 0.0036 * clay + 0.002 * om - 0.0009 * sand).clip(0, 0.6)
    wp = (0.026 + 0.005 * clay + 0.015 * om).clip(0, 0.3)
    porosity = 1 - (bd / 2.65)
    theta_s = porosity.clip(0.3, 0.6)
    theta_r = (0.01 + 0.02 * clay / 100).clip(0.01, 0.08)
    return theta_s, theta_r, fc, wp

# -----------------------------
# Paths: update accordingly
# -----------------------------
folder = r"D:\Module11\PySEBAL_data\Soil"
output_dir = r"D:\Module11\PySEBAL_data\Soil\soil_variables"
os.makedirs(output_dir, exist_ok=True)

# Depths
depths_top = ['0-5', '5-15', '15-30']
depths_sub = ['30-60', '60-100']

def get_files(var, depths):
    return sorted([f for f in glob.glob(os.path.join(folder, f"{var}_*.tif")) if any(d in f for d in depths)])

# -----------------------------
# Load and average soil layers
# -----------------------------

# Topsoil
clay_top, ref = load_and_average(get_files('clay', depths_top))
sand_top, _ = load_and_average(get_files('sand', depths_top))
ocd_top_raw, _ = load_and_average(get_files('ocd', depths_top))
ocd_top = ocd_top_raw / 10  # g/kg to %
bdod_top_raw, _ = load_and_average(get_files('bdod', depths_top))
bdod_top = bdod_top_raw / 100  # kg/m³ to g/cm³

# Subsoil
clay_sub, _ = load_and_average(get_files('clay', depths_sub))
sand_sub, _ = load_and_average(get_files('sand', depths_sub))
ocd_sub_raw, _ = load_and_average(get_files('ocd', depths_sub))
ocd_sub = ocd_sub_raw / 10
bdod_sub_raw, _ = load_and_average(get_files('bdod', depths_sub))
bdod_sub = bdod_sub_raw * 10 / 1000

# -----------------------------
# Compute soil moisture parameters
# -----------------------------
theta_s_top, theta_r_top, fc_top, wp_top = compute_soil_moisture_params(sand_top, clay_top, ocd_top, bdod_top)
theta_s_sub, theta_r_sub, fc_sub, wp_sub = compute_soil_moisture_params(sand_sub, clay_sub, ocd_sub, bdod_sub)

# -----------------------------
# Save output rasters
# -----------------------------
def save_raster(data, ref_raster, out_path):
    meta = ref_raster.meta.copy()
    meta.update(dtype='float32', count=1, compress='lzw')
    with rasterio.open(out_path, 'w', **meta) as dst:
        dst.write(data.astype(np.float32), 1)

save_raster(theta_s_top, ref, os.path.join(output_dir, "theta_s.tif"))
save_raster(theta_s_sub, ref, os.path.join(output_dir, "theta_s_subsoil.tif"))
save_raster(theta_r_top, ref, os.path.join(output_dir, "theta_r.tif"))
save_raster(theta_r_sub, ref, os.path.join(output_dir, "theta_r_subsoil.tif"))
save_raster(fc_top, ref, os.path.join(output_dir, "field_capacity.tif"))
save_raster(wp_top, ref, os.path.join(output_dir, "wilting_point.tif"))

print("✅ All soil parameters calculated and saved.")
