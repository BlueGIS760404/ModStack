himport os
os.environ['XARRAY_NO_DASK'] = '1'
os.environ['PROJ_LIB'] = r'C:\Users\Reza\anaconda3\envs\gldas_env\Library\share\proj'

import xarray as xr
import numpy as np
import rasterio
from rasterio.mask import mask
import geopandas as gpd
import shapely
from datetime import datetime, timedelta
from glob import glob

# === Paths ===
gldas_raw = r'C:\Users\Reza\Module11\PySEBAL_data\Meteo\gldas_raw'
gldas_tiffs = r'C:\Users\Reza\Module11\PySEBAL_data\Meteo\gldas_tiffs'
os.makedirs(gldas_tiffs, exist_ok=True)

# === Load shapefile ===
mississippi = gpd.read_file(r"C:\Users\Reza\Module11\PySEBAL_data\Meteo\mississippi.shp")
mississippi = mississippi.to_crs("EPSG:4326")
mississippi = mississippi[~mississippi.geometry.is_empty & mississippi.geometry.is_valid]
mississippi.geometry = mississippi.geometry.buffer(0)
geojson_geom = [shapely.geometry.mapping(geom) for geom in mississippi.geometry]
shapely_geom = mississippi.geometry.unary_union

# === Utilities ===
def compute_rh(tair, qair, pressure=101325):
    tair = np.asarray(tair)
    qair = np.asarray(qair)
    rh = np.full_like(tair, np.nan)
    valid_mask = (~np.isnan(tair)) & (~np.isnan(qair)) & (qair > 0) & (tair > 100)
    T = tair[valid_mask] - 273.15
    es = 6.112 * np.exp((17.67 * T) / (T + 243.5)) * 100
    e = qair[valid_mask] * pressure / (0.622 + 0.378 * qair[valid_mask])
    rh[valid_mask] = np.clip(100 * e / es, 0, 100)
    return rh

def k_to_c(temp_k):
    return temp_k - 273.15

# === Check if all files are valid ===
def validate_files(files):
    for f in files:
        try:
            with xr.open_dataset(f) as _:
                pass
        except Exception as e:
            print(f"  ❌ File corrupted or unreadable: {os.path.basename(f)} — {e}")
            return False
    return True

# === Main processor ===
def process_day(date):
    date_str = date.strftime("%Y%m%d")
    files = sorted(glob(os.path.join(gldas_raw, f"GLDAS_NOAH025_3H.A{date_str}*.nc4")))
    print(f"📅 {date_str} — Found {len(files)} files")

    if len(files) != 8 or not validate_files(files):
        print(f"❌ Skipping {date_str} due to missing or invalid files")
        return

    tair_stack, qair_stack, wind_stack, swdown_stack = [], [], [], []

    for f in files:
        with xr.open_dataset(f) as ds:
            tair = ds['Tair_f_inst'].squeeze().load()
            qair = ds['Qair_f_inst'].squeeze().load()
            wind = ds['Wind_f_inst'].squeeze().load()
            swdown = ds['SWdown_f_tavg'].squeeze().load()
            swdown = swdown.where(swdown >= 0, 0)

            tair_stack.append(tair)
            qair_stack.append(qair)
            wind_stack.append(wind)
            swdown_stack.append(swdown)

    tair_stack = xr.concat(tair_stack, dim='time')
    qair_stack = xr.concat(qair_stack, dim='time')
    wind_stack = xr.concat(wind_stack, dim='time')
    swdown_stack = xr.concat(swdown_stack, dim='time')

    inst_idx = 5
    temp_inst = k_to_c(tair_stack[inst_idx])
    rh_inst = compute_rh(tair_stack[inst_idx].values, qair_stack[inst_idx].values)
    wind_inst = wind_stack[inst_idx]
    rs_inst = swdown_stack[inst_idx]

    temp_24 = k_to_c(tair_stack.mean(dim='time'))
    rh_24 = compute_rh(tair_stack.mean(dim='time').values, qair_stack.mean(dim='time').values)
    wind_24 = wind_stack.mean(dim='time')
    rs_24 = swdown_stack.mean(dim='time')

    vars_to_save = {
        'Temp_inst': temp_inst,
        'Temp_24': temp_24,
        'RH_inst': xr.DataArray(rh_inst, dims=temp_inst.dims, coords=temp_inst.coords),
        'RH_24': xr.DataArray(rh_24, dims=temp_24.dims, coords=temp_24.coords),
        'Wind_inst': wind_inst,
        'Wind_24': wind_24,
        'Rs_inst': rs_inst,
        'Rs_24': rs_24
    }

    for var_name, data in vars_to_save.items():
        arr = np.flipud(data.values.astype('float32'))
        arr = np.where(np.isnan(arr), -9999, arr)

        transform = rasterio.transform.from_origin(west=-180, north=90, xsize=0.25, ysize=0.25)
        temp_file = os.path.join(gldas_tiffs, f"{var_name}_{date_str}_raw.tif")

        with rasterio.open(
            temp_file, 'w', driver='GTiff',
            height=arr.shape[0], width=arr.shape[1],
            count=1, dtype='float32',
            crs="EPSG:4326", transform=transform,
            nodata=-9999
        ) as dst:
            dst.write(arr, 1)

        try:
            with rasterio.open(temp_file) as src:
                raster_bounds = shapely.geometry.box(*src.bounds)

                if not raster_bounds.intersects(shapely_geom):
                    print(f"  ⚠️ No intersection for {var_name}, saving full extent")
                    final_path = os.path.join(gldas_tiffs, f"{var_name}_{date_str}_FULL.tif")
                    os.rename(temp_file, final_path)
                    continue

                out_image, out_transform = mask(
                    src, geojson_geom, crop=True,
                    filled=False, nodata=-9999, all_touched=True
                )

                out_meta = src.meta.copy()
                out_meta.update({
                    "height": out_image.shape[1],
                    "width": out_image.shape[2],
                    "transform": out_transform,
                    "nodata": -9999
                })

                final_path = os.path.join(gldas_tiffs, f"{var_name}_{date_str}.tif")
                with rasterio.open(final_path, "w", **out_meta) as dest:
                    dest.write(out_image)

            valid_mask = (out_image != -9999)
            if np.any(valid_mask):
                valid_values = out_image[valid_mask]
                if 'Temp' in var_name:
                    valid_values = valid_values - 273.15 if np.nanmean(valid_values) > 100 else valid_values
                print(f"  📊 {var_name} → min: {np.min(valid_values):.2f}, max: {np.max(valid_values):.2f}, mean: {np.mean(valid_values):.2f}")
            else:
                print(f"  ⚠️ No valid data for {var_name} after clipping")
            print(f"    ✅ Saved: {final_path}")

        except Exception as e:
            print(f"  ❌ Failed to process {var_name} — {str(e)}")
            fallback_path = os.path.join(gldas_tiffs, f"{var_name}_{date_str}_FULL.tif")
            os.rename(temp_file, fallback_path)

        if os.path.exists(temp_file):
            os.remove(temp_file)

    print(f"✅ Done: {date_str}")

# === Loop through dates ===
start_date = datetime(2024, 1, 1)
end_date = datetime(2024, 12, 31)
current = start_date

while current <= end_date:
    process_day(current)
    current += timedelta(days=1)
