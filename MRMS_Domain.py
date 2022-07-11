#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os,math,pyproj,subprocess,fiona
import rasterio
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import xarray as xr
import rioxarray as rxr
import geopandas as gpd
from glob import glob
from netCDF4 import Dataset
from functools import partial
from pyproj import CRS
from pyproj import Transformer
from shapely.geometry import Point, mapping
from shapely.ops import transform

#os.environ['GDAL_DATA'] = os.environ['CONDA_PREFIX'] + r'\Library\share\gdal'
#os.environ['PROJ_LIB'] = os.environ['CONDA_PREFIX'] + r'\Library\share'


# In[ ]:


## input params
start_yr = 2015
end_yr = 2021
city = 'CITY NAME'
lon = 
lat = 
buf_dis1 = 60 # radius of the round buffer, in km
buf_dis2 = 40
buf_dis3 = 20
ang_range = 90


# In[ ]:


def geodesic_point_buffer(lat, lon, km):
    # Azimuthal equidistant projection
    aeqd_proj = f"+proj=aeqd +lat_0={lat} +lon_0={lon} +x_0=0 +y_0=0"
    transformer = Transformer.from_crs("EPSG:4326", aeqd_proj, always_xy=True)
    point = Point(lon, lat)
    point_aeqd = transform(transformer.transform, point)
    circle_aeqd = point_aeqd.buffer(km * 1000)
    return mapping(transform(partial(transformer.transform, direction="INVERSE"), circle_aeqd))


# In[ ]:


# check wind data
## MERRA2 data (M2T1NXSLV) for computing wind direction

for year in range(start_yr,end_yr+1):
    # create directories
    if not os.path.exists('./MRMS_PrecipData/'+str(year)):
        os.mkdir('./MRMS_PrecipData/'+str(year))
    
    filename_wind = './MRMS_PrecipData/'+str(year)+'/wind_drct_'+city+'.txt'
    if not os.path.exists(filename_wind):
        wind_file = sorted(glob('./MRMS_PrecipData/'+str(year)+'/wind/*.nc4'))
        wind_data = Dataset(wind_file[0])

        idx_y = (np.abs(lat-wind_data.variables['lat'][:])).argmin()
        idx_x = (np.abs(lon-wind_data.variables['lon'][:])).argmin()

        try:
            os.remove(filename_wind)  # remove the file to store wind directions if it exists
        except OSError:
            pass

        for i in range(len(wind_file)):
            wind_data = xr.open_dataset(wind_file[i])
            for j in range(len(wind_data.time)):
                name = wind_data.time[j].values
                name = str(name)[:13]
                wnd_deg = np.rad2deg(math.atan2(wind_data.V850.values[j,idx_y,idx_x],wind_data.U850.values[j,idx_y,idx_x]))+180
                wnd_vel = np.sqrt(wind_data.V850.values[j,idx_y,idx_x]**2+wind_data.U850.values[j,idx_y,idx_x]**2)
                with open(filename_wind,'a') as f:
                    f.write('%s %s %s\n' % (name,round(wnd_deg,2),round(wnd_vel,2)))


# In[ ]:


for year in range(start_yr,end_yr+1):
    # create directories
    if not os.path.exists('./MRMS_PrecipData/'+str(year)):
        os.mkdir('./MRMS_PrecipData/'+str(year))
    
    filename_wind = './MRMS_PrecipData/'+str(year)+'/wind_drct_'+city+'.txt'
    
    # read wind direction
    wind = pd.read_csv(filename_wind,delimiter = ' ',header = None)
    wind.columns = ['time','direction','speed']
    wind = wind.astype({'time':'str'})
    
    # QPE precipitation data
    nc_file = glob('./MRMS_PrecipData/'+str(year)+'/nc_file/'+str(year)+'0[6-8][0-3][0-9]*.nc')
    nc_file = sorted(nc_file)
    
    # stats output initialization
    header = "time mean median min max std sum percentile(25) percentile(75) percentile(10) percentile(20) percentile(30) percentile(40) xgrids(90) xgrids(95)"
    
    buffer1 = geodesic_point_buffer(lat,lon,buf_dis1) # 60km buffer
    buffer2 = geodesic_point_buffer(lat,lon,buf_dis2)
    buffer3 = geodesic_point_buffer(lat,lon,buf_dis3)
    
    # domain stats
    if not os.path.exists('./MRMS_PrecipData/'+str(year)+'/Results/domain stat'):
        os.mkdir('./MRMS_PrecipData/'+str(year)+'/Results/domain stat/Test')
    
    filename_domain1 = './MRMS_PrecipData/'+str(year)+'/Results/domain stat/'+city+str(year)+'_'+str(buf_dis1)+'.txt'
    filename_domain2 = './MRMS_PrecipData/'+str(year)+'/Results/domain stat/'+city+str(year)+'_'+str(buf_dis2)+'.txt'
    filename_domain3 = './MRMS_PrecipData/'+str(year)+'/Results/domain stat/'+city+str(year)+'_'+str(buf_dis3)+'.txt'
    header = "time domain_avg domain_std"
    try:
        os.remove(filename_domain1)  # remove the file if it exists
    except OSError:
        pass
    try:  # remove the file if it exists
        os.remove(filename_domain2)
    except OSError:
        pass
    try:
        os.remove(filename_domain3)
    except OSError:
        pass 
    
    with open(filename_domain1,'a') as f1:
        f1.write(header + '\n')
        f1.close
    with open(filename_domain2,'a') as f2:
        f2.write(header + '\n')
        f2.close
    with open(filename_domain3,'a') as f3:
        f3.write(header + '\n')
        f3.close
    
    for item in nc_file:
    ##extract time
        time = item.split('/')[-1]
        time = time[:4]+'-'+time[4:6]+'-'+time[6:8]+'T'+time[9:11]
        print('%s is being processed...'%time)
    
    ##get wind direction value
        wind_val = wind[wind['time']==time]['direction'].values
    
    ## open circular buffer
        nc_data = xr.open_dataset(item)
        raw_buffer1 = nc_data[list(nc_data.keys())[0]].rio.set_crs('WGS84').rio.clip([buffer1])
        raw_buffer2 = nc_data[list(nc_data.keys())[0]].rio.set_crs('WGS84').rio.clip([buffer2])
        raw_buffer3 = nc_data[list(nc_data.keys())[0]].rio.set_crs('WGS84').rio.clip([buffer3])
        #print(raw_buffer.mean(skipna=True).values)
        
        with open(filename_domain1,'a') as f1:
            f1.write('%s %s %s\n' % (time,raw_buffer1.mean(skipna=True).values,raw_buffer1.std(skipna=True).values))
        with open(filename_domain2,'a') as f2:
            f2.write('%s %s %s\n' % (time,raw_buffer2.mean(skipna=True).values,raw_buffer2.std(skipna=True).values))
        with open(filename_domain3,'a') as f3:
            f3.write('%s %s %s\n' % (time,raw_buffer2.mean(skipna=True).values,raw_buffer3.std(skipna=True).values))

