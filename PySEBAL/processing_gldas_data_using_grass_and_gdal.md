### Install OSGeo4W and confirmation
# Step 1: Follow steps mentioned in the link below:
https://pysebal.readthedocs.io/en/latest/installation.html

# Step 2: Open OSGeo4W Shell and enter the following
# Install following packages
pip3 install openpyxl joblib
pip3 install grass_session

# Step 3: Confirm the installation
# Specify drive
D:
# Change to the directory with SEBAL code
cd PySEBAL_dev\SEBAL
# open python
python
# import one of the PySEBAL Script
import pysebal_py3
# If there are no errors, the installation is successful
# To exit from python
exit()





### Extracting meteo variables from GLDAS
Repeat the following steps for all NC files, including:
	- GLDAS_NOAH025_3H.A20240312.0000.021.nc4
	- GLDAS_NOAH025_3H.A20240312.0300.021.nc4,
	- GLDAS_NOAH025_3H.A20240312.0600.021.nc4,
	- GLDAS_NOAH025_3H.A20240312.0900.021.nc4,
	- GLDAS_NOAH025_3H.A20240312.1200.021.nc4,
	- GLDAS_NOAH025_3H.A20240312.1500.021.nc4,
	- GLDAS_NOAH025_3H.A20240312.1800.021.nc4,
	- GLDAS_NOAH025_3H.A20240312.2100.021.nc4
	
# Step 1: Specify drive
D:
Change to the directory with meteo file (GLDAS)
cd D:\Module11\PySEBAL_data\Meteo\GLDAS_NOAH025_3H.A20240312.0000.021

# Step 2: See the metadata of GLDAS file
gdalinfo GLDAS_NOAH025_3H.A20240312.0000.021.nc4

# Step 3: Extracting meteo variables from GLDAS for the whole world
gdal_translate -a_nodata -9999 NETCDF:"GLDAS_NOAH025_3H.A20240312.0000.021.nc4":Qair_f_inst GLDAS_NOAH025_3H_20240312_0000_Qair.tif
gdal_translate -a_nodata -9999 NETCDF:"GLDAS_NOAH025_3H.A20240312.0000.021.nc4":Psurf_f_inst GLDAS_NOAH025_3H_20240312_0000_Psurf.tif
gdal_translate -a_nodata -9999 NETCDF:"GLDAS_NOAH025_3H.A20240312.0000.021.nc4":Tair_f_inst GLDAS_NOAH025_3H_20240312_0000_Tair.tif
gdal_translate -a_nodata -9999 NETCDF:"GLDAS_NOAH025_3H.A20240312.0000.021.nc4":Wind_f_inst GLDAS_NOAH025_3H_20240312_0000_Wind.tif
gdal_translate -a_nodata -9999 NETCDF:"GLDAS_NOAH025_3H.A20240312.0000.021.nc4":SWdown_f_tavg GLDAS_NOAH025_3H_20240312_0000_SWdown.tif

# Step 4: Clip to Mississippi shapefile
gdalwarp -overwrite -s_srs EPSG:4326 -cutline mississippi_boundary.shp -crop_to_cutline -dstnodata -9999 GLDAS_NOAH025_3H_20240312_0000_Qair.tif GLDAS_NOAH025_3H_20240312_0000_Qair_clipped.tif
gdalwarp -overwrite -s_srs EPSG:4326 -cutline mississippi_boundary.shp -crop_to_cutline -dstnodata -9999 GLDAS_NOAH025_3H_20240312_0000_Psurf.tif GLDAS_NOAH025_3H_20240312_0000_Psurf_clipped.tif
gdalwarp -overwrite -s_srs EPSG:4326 -cutline mississippi_boundary.shp -crop_to_cutline -dstnodata -9999 GLDAS_NOAH025_3H_20240312_0000_Tair.tif GLDAS_NOAH025_3H_20240312_0000_Tair_clipped.tif
gdalwarp -overwrite -s_srs EPSG:4326 -cutline mississippi_boundary.shp -crop_to_cutline -dstnodata -9999 GLDAS_NOAH025_3H_20240312_0000_Wind.tif GLDAS_NOAH025_3H_20240312_0000_Wind_clipped.tif
gdalwarp -overwrite -s_srs EPSG:4326 -cutline mississippi_boundary.shp -crop_to_cutline -dstnodata -9999 GLDAS_NOAH025_3H_20240312_0000_SWdown.tif GLDAS_NOAH025_3H_20240312_0000_SWdown_clipped.tif

# Step 5: Remove the initial tif files and rename clipped tif files  





### Import created maps into GRASS GIS
# Step 1: Open GRASS GIS GUI
grass84 --gui

# Step 2: import these tif files into GRASS GIS
r.import input="GLDAS_NOAH025_3H_20240312_0000_Qair.tif" output="GLDAS_NOAH025_3H_20240312_0000_Qair" -o --overwrite
r.import input="GLDAS_NOAH025_3H_20240312_0000_Psurf.tif" output="GLDAS_NOAH025_3H_20240312_0000_Psurf" -o --overwrite
r.import input="GLDAS_NOAH025_3H_20240312_0000_Tair.tif" output="GLDAS_NOAH025_3H_20240312_0000_Tair" -o --overwrite
r.import input="GLDAS_NOAH025_3H_20240312_0000_Wind.tif" output="GLDAS_NOAH025_3H_20240312_0000_Wind" -o --overwrite
r.import input="GLDAS_NOAH025_3H_20240312_0000_SWdown.tif" output="GLDAS_NOAH025_3H_20240312_0000_SWdown" -o --overwrite

# Step 3: Seting the computational region in GRASS GIS
g.region res=0.25 -a
g.region -p





### Perform conversions and export the layers
# Air temperature convert from Kelvin to Celsius
r.mapcalc "GLDAS_NOAH025_3H_20240312_0000_Tair_deg = GLDAS_NOAH025_3H_20240312_0000_Tair - 273.15" --o
# Export the Air temperature (Celsius) map
r.out.gdal input=GLDAS_NOAH025_3H_20240312_0000_Tair_deg output="D:/Module11/PySEBAL_data/Meteo/GLDAS_NOAH025_3H.A20240312.0000.021/GLDAS_NOAH025_3H_20240312_0000_Tair_deg.tif" format=GTiff createopt="COMPRESS=DEFLATE" --overwrite
# Pressure convert from pa to mb
r.mapcalc "GLDAS_NOAH025_3H_20240312_0000_Psurf_mb = GLDAS_NOAH025_3H_20240312_0000_Psurf / 100" --o
# Export the Pressure (mb) map
r.out.gdal input=GLDAS_NOAH025_3H_20240312_0000_Psurf_mb output="D:/Module11/PySEBAL_data/Meteo/GLDAS_NOAH025_3H.A20240312.0000.021/GLDAS_NOAH025_3H_20240312_0000_Psurf_mb.tif" format=GTiff createopt="COMPRESS=DEFLATE" --overwrite
# Saturation vapour pressure
r.mapcalc "es = 6.112 * exp((17.67 * GLDAS_NOAH025_3H_20240312_0000_Tair_deg) / (GLDAS_NOAH025_3H_20240312_0000_Tair_deg + 243.5))" --o
# Export the Saturation vapour pressure map
r.out.gdal input=es output="D:/Module11/PySEBAL_data/Meteo/GLDAS_NOAH025_3H.A20240312.0000.021/es.tif" format=GTiff createopt="COMPRESS=DEFLATE" --overwrite
# vapour pressure
r.mapcalc "e = (GLDAS_NOAH025_3H_20240312_0000_Qair * GLDAS_NOAH025_3H_20240312_0000_Psurf_mb) / (0.378 * GLDAS_NOAH025_3H_20240312_0000_Qair + 0.622)" --o
# Export the vapour pressure map
r.out.gdal input=e output="D:/Module11/PySEBAL_data/Meteo/GLDAS_NOAH025_3H.A20240312.0000.021/e.tif" format=GTiff createopt="COMPRESS=DEFLATE" --overwrite
# Calculate relative humidity
r.mapcalc "GLDAS_NOAH025_3H_20240312_0000_Rh1 = (e / es) * 100" --o
# Export the relative humidity map
r.out.gdal input=GLDAS_NOAH025_3H_20240312_0000_Rh1 output="D:/Module11/PySEBAL_data/Meteo/GLDAS_NOAH025_3H.A20240312.0000.021/GLDAS_NOAH025_3H_20240312_0000_Rh1.tif" format=GTiff createopt="COMPRESS=DEFLATE" --overwrite
# Remove outliers
r.mapcalc "GLDAS_NOAH025_3H_20240312_0000_Rh = float(if(GLDAS_NOAH025_3H_20240312_0000_Rh1 > 100, 100, if(GLDAS_NOAH025_3H_20240312_0000_Rh1 < 0, 0, GLDAS_NOAH025_3H_20240312_0000_Rh1)))" --o
# Export the relative humidity_no_outliers map
r.out.gdal input=GLDAS_NOAH025_3H_20240312_0000_Rh output="D:/Module11/PySEBAL_data/Meteo/GLDAS_NOAH025_3H.A20240312.0000.021/GLDAS_NOAH025_3H_20240312_0000_Rh.tif" format=GTiff createopt="COMPRESS=DEFLATE" --overwrite





### Create instantaneous averages
# Specify drive
D:
Change to the directory with meteo file (GLDAS)
cd D:\Module11\PySEBAL_data\Meteo\GLDAS_NOAH025_3H.A20240312.0900.021

# Open GRASS GIS GUI
grass84 --gui

# Importing GLDAS_NOAH025_3H_20240312_0900 layers (Tair_deg, SWdown, Wind and Rh) into GRASS GIS (as Landsat acquisition time is around 8:30 GMT) 
r.import input="GLDAS_NOAH025_3H_20240312_0900_Tair_deg.tif" output="GLDAS_NOAH025_3H_20240312_0900_Tair_deg" -o --overwrite
r.import input="GLDAS_NOAH025_3H_20240312_0900_SWdown.tif" output="GLDAS_NOAH025_3H_20240312_0900_SWdown" -o --overwrite
r.import input="GLDAS_NOAH025_3H_20240312_0900_Wind.tif" output="GLDAS_NOAH025_3H_20240312_0900_Wind" -o --overwrite
r.import input="GLDAS_NOAH025_3H_20240312_0900_Rh.tif" output="GLDAS_NOAH025_3H_20240312_0900_Rh" -o --overwrite

# Air temperature instantaneous
r.mapcalc "GLDAS_NOAH025_3H_20240312_Tair_inst = GLDAS_NOAH025_3H_20240312_0900_Tair_deg"
# Export the air temperature instantaneous map
r.out.gdal input=GLDAS_NOAH025_3H_20240312_Tair_inst output="D:\Module11\PySEBAL_data\Meteo\Instantaneous_Averages/GLDAS_NOAH025_3H_20240312_Tair_inst.tif" format=GTiff createopt="COMPRESS=DEFLATE" --overwrite
# Shortwave radiation instantaneous
r.mapcalc "GLDAS_NOAH025_3H_20240312_SWdown_inst = GLDAS_NOAH025_3H_20240312_0900_SWdown"
# Export the shortwave radiation instantaneous map
r.out.gdal input=GLDAS_NOAH025_3H_20240312_SWdown_inst output="D:\Module11\PySEBAL_data\Meteo\Instantaneous_Averages/GLDAS_NOAH025_3H_20240312_SWdown_inst.tif" format=GTiff createopt="COMPRESS=DEFLATE" --overwrite
# Wind speed instantaneous
r.mapcalc "GLDAS_NOAH025_3H_20240312_Wind_inst = GLDAS_NOAH025_3H_20240312_0900_Wind"
# Export the wind speed instantaneous map
r.out.gdal input=GLDAS_NOAH025_3H_20240312_Wind_inst output="D:\Module11\PySEBAL_data\Meteo\Instantaneous_Averages/GLDAS_NOAH025_3H_20240312_Wind_inst.tif" format=GTiff createopt="COMPRESS=DEFLATE" --overwrite
## Relative humidity instantaneous
r.mapcalc "GLDAS_NOAH025_3H_20240312_Rh_inst = GLDAS_NOAH025_3H_20240312_0900_Rh"
# Export the relative humidity instantaneous map
r.out.gdal input=GLDAS_NOAH025_3H_20240312_Rh_inst output="D:\Module11\PySEBAL_data\Meteo\Instantaneous_Averages/GLDAS_NOAH025_3H_20240312_Rh_inst.tif" format=GTiff createopt="COMPRESS=DEFLATE" --overwrite





### Create daily averages
Repeat the following steps for all 4 variables, including:
	- GLDAS_NOAH025_3H_20240312_Tair_deg
	- GLDAS_NOAH025_3H_20240312_Wind,
	- GLDAS_NOAH025_3H_20240312_SWdown,
	- GLDAS_NOAH025_3H_20240312_Rh

# Specify drive
D:
Change to the directory with meteo file (GLDAS)
cd D:\Module11\PySEBAL_data\Meteo\Daily_Averages\GLDAS_NOAH025_3H_20180606_Tair_deg

# Open GRASS GIS GUI
grass84 --gui

# Importing 8 tif layers extracted from each GLDAS .nc4 file into GRASS GIS
r.import input="GLDAS_NOAH025_3H_20240312_0000_Tair_deg.tif" output="GLDAS_NOAH025_3H_20240312_0000_Tair_deg" -o --overwrite
r.import input="GLDAS_NOAH025_3H_20240312_0300_Tair_deg.tif" output="GLDAS_NOAH025_3H_20240312_0300_Tair_deg" -o --overwrite
r.import input="GLDAS_NOAH025_3H_20240312_0600_Tair_deg.tif" output="GLDAS_NOAH025_3H_20240312_0600_Tair_deg" -o --overwrite
r.import input="GLDAS_NOAH025_3H_20240312_0900_Tair_deg.tif" output="GLDAS_NOAH025_3H_20240312_0900_Tair_deg" -o --overwrite
r.import input="GLDAS_NOAH025_3H_20240312_1200_Tair_deg.tif" output="GLDAS_NOAH025_3H_20240312_1200_Tair_deg" -o --overwrite
r.import input="GLDAS_NOAH025_3H_20240312_1500_Tair_deg.tif" output="GLDAS_NOAH025_3H_20240312_1500_Tair_deg" -o --overwrite
r.import input="GLDAS_NOAH025_3H_20240312_1800_Tair_deg.tif" output="GLDAS_NOAH025_3H_20240312_1800_Tair_deg" -o --overwrite
r.import input="GLDAS_NOAH025_3H_20240312_2100_Tair_deg.tif" output="GLDAS_NOAH025_3H_20240312_2100_Tair_deg" -o --overwrite

# Averaging the values for all 8 layers
# Air temperature daily average
r.series input=GLDAS_NOAH025_3H_20240312_0000_Tair_deg,GLDAS_NOAH025_3H_20240312_0300_Tair_deg,GLDAS_NOAH025_3H_20240312_0600_Tair_deg,GLDAS_NOAH025_3H_20240312_0900_Tair_deg,GLDAS_NOAH025_3H_20240312_1200_Tair_deg,GLDAS_NOAH025_3H_20240312_1500_Tair_deg,GLDAS_NOAH025_3H_20240312_1800_Tair_deg,GLDAS_NOAH025_3H_20240312_2100_Tair_deg output=GLDAS_NOAH025_3H_20240312_Tair_deg_daily method=average
# Export the air temperature daily average map
r.out.gdal input=GLDAS_NOAH025_3H_20240312_Tair_deg_daily output="D:\Module11\PySEBAL_data\Meteo\Daily_Averages/GLDAS_NOAH025_3H_20240312_Tair_deg_daily.tif" format=GTiff createopt="COMPRESS=DEFLATE" --overwrite

# Resample the instantaneous and daily averaged layers
# Save the following commands inside a bat file named "resample_instantaneous.bat"
@echo off
REM --- Set paths and variables
SET "MAPSET=C:\Users\Reza\Documents\grassdata\world_latlong_wgs84\PERMANENT"
SET "DATA_DIR=D:\Module11\PySEBAL_data\Meteo\Instantaneous_Averages"
SET "DATE=20240312"

REM --- Use the full path to grass84.bat
SET "GRASS_CMD=C:\OSGeo4W\bin\grass84.bat"

REM --- Check if GRASS command exists
if not exist "%GRASS_CMD%" (
    echo Error: GRASS command not found: %GRASS_CMD%
    pause
    exit /b 1
)

REM --- Check if directories exist
if not exist "%MAPSET%" (
    echo Error: MAPSET directory does not exist: %MAPSET%
    pause
    exit /b 1
)
if not exist "%DATA_DIR%" (
    echo Error: DATA_DIR directory does not exist: %DATA_DIR%
    pause
    exit /b 1
)

REM --- Check if input files exist
if not exist "%DATA_DIR%\GLDAS_NOAH025_3H_%DATE%_Tair_deg_inst.tif" (
    echo Error: Input file does not exist: %DATA_DIR%\GLDAS_NOAH025_3H_%DATE%_Tair_deg_inst.tif
    pause
    exit /b 1
)
if not exist "%DATA_DIR%\GLDAS_NOAH025_3H_%DATE%_Wind_inst.tif" (
    echo Error: Input file does not exist: %DATA_DIR%\GLDAS_NOAH025_3H_%DATE%_Wind_inst.tif
    pause
    exit /b 1
)
if not exist "%DATA_DIR%\GLDAS_NOAH025_3H_%DATE%_SWdown_inst.tif" (
    echo Error: Input file does not exist: %DATA_DIR%\GLDAS_NOAH025_3H_%DATE%_SWdown_inst.tif
    pause
    exit /b 1
)
if not exist "%DATA_DIR%\GLDAS_NOAH025_3H_%DATE%_Rh_inst.tif" (
    echo Error: Input file does not exist: %DATA_DIR%\GLDAS_NOAH025_3H_%DATE%_Rh_inst.tif
    pause
    exit /b 1
)

REM --- Import raster maps
echo Importing raster maps...
call "%GRASS_CMD%" "%MAPSET%" --exec r.import input="%DATA_DIR%\GLDAS_NOAH025_3H_%DATE%_Tair_deg_inst.tif" output="GLDAS_NOAH025_3H_%DATE%_Tair_deg_inst" --overwrite -l --verbose
if errorlevel 1 (
    echo Failed to import Tair_deg_inst
    pause
    exit /b 1
)

call "%GRASS_CMD%" "%MAPSET%" --exec r.import input="%DATA_DIR%\GLDAS_NOAH025_3H_%DATE%_Wind_inst.tif" output="GLDAS_NOAH025_3H_%DATE%_Wind_inst" --overwrite -l --verbose
if errorlevel 1 (
    echo Failed to import Wind_inst
    pause
    exit /b 1
)

call "%GRASS_CMD%" "%MAPSET%" --exec r.import input="%DATA_DIR%\GLDAS_NOAH025_3H_%DATE%_SWdown_inst.tif" output="GLDAS_NOAH025_3H_%DATE%_SWdown_inst" --overwrite -l --verbose
if errorlevel 1 (
    echo Failed to import SWdown_inst
    pause
    exit /b 1
)

call "%GRASS_CMD%" "%MAPSET%" --exec r.import input="%DATA_DIR%\GLDAS_NOAH025_3H_%DATE%_Rh_inst.tif" output="GLDAS_NOAH025_3H_%DATE%_Rh_inst" --overwrite -l --verbose
if errorlevel 1 (
    echo Failed to import Rh_inst
    pause
    exit /b 1
)

REM --- Check raster bounds
echo Checking raster bounds...
call "%GRASS_CMD%" "%MAPSET%" --exec r.info map=GLDAS_NOAH025_3H_%DATE%_Tair_deg_inst
if errorlevel 1 (
    echo Failed to retrieve raster info for Tair_deg_inst
    pause
    exit /b 1
)

REM --- Set computational region with explicit numerical bounds
echo Setting computational region...
call "%GRASS_CMD%" "%MAPSET%" --exec g.region north=90 south=-90 east=180 west=-180 nsres=0.25 ewres=0.25 -p
if errorlevel 1 (
    echo Failed to set computational region
    pause
    exit /b 1
)

REM --- Verify computational region
echo Verifying computational region...
call "%GRASS_CMD%" "%MAPSET%" --exec g.region -p
if errorlevel 1 (
    echo Failed to verify computational region
    pause
    exit /b 1
)

REM --- Resample maps using r.resamp.interp
echo Resampling maps...
call "%GRASS_CMD%" "%MAPSET%" --exec r.resamp.interp input="GLDAS_NOAH025_3H_%DATE%_Tair_deg_inst" output="GLDAS_NOAH025_3H_%DATE%_Tair_deg_inst_interp" method=bicubic --overwrite
if errorlevel 1 (
    echo Failed to resample Tair_deg_inst
    pause
    exit /b 1
)

call "%GRASS_CMD%" "%MAPSET%" --exec r.resamp.interp input="GLDAS_NOAH025_3H_%DATE%_Wind_inst" output="GLDAS_NOAH025_3H_%DATE%_Wind_inst_interp" method=bicubic --overwrite
if errorlevel 1 (
    echo Failed to resample Wind_inst
    pause
    exit /b 1
)

call "%GRASS_CMD%" "%MAPSET%" --exec r.resamp.interp input="GLDAS_NOAH025_3H_%DATE%_SWdown_inst" output="GLDAS_NOAH025_3H_%DATE%_SWdown_inst_interp" method=bicubic --overwrite
if errorlevel 1 (
    echo Failed to resample SWdown_inst
    pause
    exit /b 1
)

call "%GRASS_CMD%" "%MAPSET%" --exec r.resamp.interp input="GLDAS_NOAH025_3H_%DATE%_Rh_inst" output="GLDAS_NOAH025_3H_%DATE%_Rh_inst_interp" method=bicubic --overwrite
if errorlevel 1 (
    echo Failed to resample Rh_inst
    pause
    exit /b 1
)

REM --- Export maps
echo Exporting resampled maps...
call "%GRASS_CMD%" "%MAPSET%" --exec r.out.gdal input="GLDAS_NOAH025_3H_%DATE%_Tair_deg_inst_interp" output="%DATA_DIR%\GLDAS_NOAH025_3H_%DATE%_Tair_deg_inst_interp.tif" format=GTiff createopt="COMPRESS=DEFLATE" --overwrite
if errorlevel 1 (
    echo Failed to export Tair_deg_inst_interp
    pause
    exit /b 1
)

call "%GRASS_CMD%" "%MAPSET%" --exec r.out.gdal input="GLDAS_NOAH025_3H_%DATE%_Wind_inst_interp" output="%DATA_DIR%\GLDAS_NOAH025_3H_%DATE%_Wind_inst_interp.tif" format=GTiff createopt="COMPRESS=DEFLATE" --overwrite
if errorlevel 1 (
    echo Failed to export Wind_inst_interp
    pause
    exit /b 1
)

call "%GRASS_CMD%" "%MAPSET%" --exec r.out.gdal input="GLDAS_NOAH025_3H_%DATE%_SWdown_inst_interp" output="%DATA_DIR%\GLDAS_NOAH025_3H_%DATE%_SWdown_inst_interp.tif" format=GTiff createopt="COMPRESS=DEFLATE" --overwrite
if errorlevel 1 (
    echo Failed to export SWdown_inst_interp
    pause
    exit /b 1
)

call "%GRASS_CMD%" "%MAPSET%" --exec r.out.gdal input="GLDAS_NOAH025_3H_%DATE%_Rh_inst_interp" output="%DATA_DIR%\GLDAS_NOAH025_3H_%DATE%_Rh_inst_interp.tif" format=GTiff createopt="COMPRESS=DEFLATE" --overwrite
if errorlevel 1 (
    echo Failed to export Rh_inst_interp
    pause
    exit /b 1
)

echo All operations completed successfully.
pause
