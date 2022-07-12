#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os,math
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
#from windrose import WindroseAxes


# In[ ]:


# Drawing km circle around what we're interested in
def geodesic_point_buffer(lat, lon, km):
    # Azimuthal equidistant projection
    aeqd_proj = f"+proj=aeqd +lat_0={lat} +lon_0={lon} +x_0=0 +y_0=0"
    transformer = Transformer.from_crs("EPSG:4326", aeqd_proj, always_xy=True)
    point = Point(lon, lat)
    point_aeqd = transform(transformer.transform, point)
    circle_aeqd = point_aeqd.buffer(km * 1000)
    return mapping(transform(partial(transformer.transform, direction="INVERSE"), circle_aeqd))


# ## spatial map distinguished by wind speed

# In[ ]:


start_yr = 2015
end_yr = 2021
ls_city = ['Atlanta','Austin','Charlotte','Cincinnati','Columbus','Dallas','Houston','Indianapolis','KC','Louisville','Memphis','Miami','Minneapolis','Nashville','NYC','OKC','Omaha','Orlando',
           'Phoenix','Pittsburgh','Richmond','SanAntonio','StLouis','Tucson','Philadelphia','DC','Baltimore']
# 
ls_lon = [-84.3880,-97.7431,-80.8431,-84.5120,-82.9988,-96.7970,-95.3698,-86.1581,-94.5786,-85.7585,-90.0490,-80.1918,-93.2650,-86.7816,-74.0060,-97.5164,-95.9345,
          -81.3789,-112.0740,-79.9959,-77.4360,-98.4936,-90.1994,-110.9747,-75.1652,-77.0369,-76.6122]
#
ls_lat = [33.7490,30.2672,35.2271,39.1031,39.9612,32.7767,29.7604,39.7684,39.0997,38.2527,35.1495,25.7617,44.9778,36.1627,40.7128,35.4676,41.2565,28.5384,
          33.4484,40.4406,37.5407,29.4241,38.6270,32.2226,39.9526,38.9072,39.2904]
#
buf_dis = 60

os.chdir('/glade/scratch/YOUR_USERNAME/MRMS_PrecipData/')

wnd_stat = pd.read_csv('DmntWndStats.csv',sep=' ')

for city,lon,lat in zip(ls_city,ls_lon,ls_lat):
    
    ## get quantile values for wind speed
    wnd_q1 = wnd_stat[wnd_stat['City']==city]['Q1'].values[0]
    wnd_q3 = wnd_stat[wnd_stat['City']==city]['Q3'].values[0]
    
    wind = pd.DataFrame() # store wind data for all years
    
    ## QPE precipitation data
    ls_nc = []
    for year in range(start_yr,end_yr+1):
        nc_file = glob(str(year)+'/nc_file/'+str(year)+'0[6-8][0-3][0-9]*.nc')
        nc_file = sorted(nc_file)
        ls_nc = ls_nc + nc_file
        
        # get hourly wind data
        file_wind = str(year)+'/wind_drct_'+city+'.txt'
        tmp_wind = pd.read_csv(file_wind,sep = ' ',header = None)
        wind = pd.concat([wind,tmp_wind],ignore_index = True)
        
    wind.columns = ['time','direction','speed']
    wind = wind.astype({'time':'str'})
         
    buffer = geodesic_point_buffer(lat,lon,buf_dis) # create buffer
    
    res_h = 0
    res_l = 0
    for item in ls_nc:
        
        '''extract time'''
        time = item.split('/')[-1]
        time = time[:4]+'-'+time[4:6]+'-'+time[6:8]+'T'+time[9:11]
        
        if wind[wind['time']==time]['speed'].values[0] <= wnd_q1:
            print('%s is being processed (Q1)...'%time)
        
            nc_data = xr.open_dataset(item)
            raw_buffer = nc_data[list(nc_data.keys())[0]].rio.set_crs('WGS84').rio.clip([buffer])
            res_l = res_l + raw_buffer[0]
            
        if wind[wind['time']==time]['speed'].values[0] >= wnd_q3:
            print('%s is being processed (Q3)...'%time)
        
            nc_data = xr.open_dataset(item)
            raw_buffer = nc_data[list(nc_data.keys())[0]].rio.set_crs('WGS84').rio.clip([buffer])
            res_h = res_h + raw_buffer[0]
    res_l.rio.to_raster('Wnd_Quantile/'+city+'_Precip60Q1.tif')
    res_h.rio.to_raster('Wnd_Quantile/'+city+'_Precip60Q3.tif')

