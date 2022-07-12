#!/usr/bin/env python
# coding: utf-8

# In[ ]:


## normalize by hourly data

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


# In[ ]:


os.chdir('/glade/scratch/YOUR_USERNAME/MRMS_PrecipData/')

wnd_stat = pd.read_csv('DmntWndStats.csv',sep=' ')

for city,lon,lat in zip(ls_city[24:27],ls_lon[24:27],ls_lat[24:27]):
    
    ## get quantile values for wind speed
    wnd_q1 = wnd_stat[wnd_stat['City']==city]['Q1'].values[0]
    wnd_q2 = wnd_stat[wnd_stat['City']==city]['Median'].values[0]
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
    
    res_q1 = 0
    res_q2 = 0
    res_q3 = 0
    res_q4 = 0
    
    count_q1 = 0
    count_q2 = 0
    count_q3 = 0
    count_q4 = 0
    
    for item in ls_nc:
        
        '''extract time'''
        time = item.split('/')[-1]
        time = time[:4]+'-'+time[4:6]+'-'+time[6:8]+'T'+time[9:11]
        
        nc_data = xr.open_dataset(item)
        raw_buffer = nc_data[list(nc_data.keys())[0]].rio.set_crs('WGS84').rio.clip([buffer])
        raw_buffer = raw_buffer[0]
        
        raw_mean = raw_buffer.mean().values
        raw_std = raw_buffer.std().values
        
        if raw_mean >= 0.1:
            
            raw_buffer = (raw_buffer-raw_mean)/raw_std # normalization
            
            if wind[wind['time']==time]['speed'].values[0] <= wnd_q1:
                print('%s is being processed (Q1)...'%time)
                res_q1 = res_q1+raw_buffer
                count_q1 = count_q1 + 1
            if wnd_q1 < wind[wind['time']==time]['speed'].values[0] < wnd_q2:
                print('%s is being processed (Q1-Q2)...'%time)
                res_q2 = res_q2+raw_buffer
                count_q2 = count_q2 + 1
            if wnd_q2 < wind[wind['time']==time]['speed'].values[0] < wnd_q3:
                print('%s is being processed (Q2-Q3)...'%time)
                res_q3 = res_q3+raw_buffer
                count_q3 = count_q3 + 1
            if wind[wind['time']==time]['speed'].values[0] >= wnd_q3:
                print('%s is being processed (Q3)...'%time)
                res_q4 = res_q4+raw_buffer
                count_q4 = count_q4 + 1
                
    res_q1 = res_q1/count_q1
    res_q2 = res_q2/count_q2
    res_q3 = res_q3/count_q3
    res_q4 = res_q4/count_q4
    
    res_q1.rio.to_raster('Precip_map/Wnd_Quantile/Hourly_Normalization/'+city+'_Precip60Q1.tif')
    res_q2.rio.to_raster('Precip_map/Wnd_Quantile/Hourly_Normalization/'+city+'_Precip60Q1_Q2.tif')
    res_q3.rio.to_raster('Precip_map/Wnd_Quantile/Hourly_Normalization/'+city+'_Precip60Q2_Q3.tif')
    res_q4.rio.to_raster('Precip_map/Wnd_Quantile/Hourly_Normalization/'+city+'_Precip60Q3.tif')

