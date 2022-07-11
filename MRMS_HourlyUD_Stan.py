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
#import rioxarray as rxr
import geopandas as gpd
from glob import glob
from netCDF4 import Dataset
from functools import partial
from pyproj import CRS
from pyproj import Transformer
from shapely.geometry import Point, mapping
from shapely.ops import transform


# In[ ]:


## input params
start_yr = 2015
end_yr = 2021
city = 'CITY NAME'
lon = -84.3880
lat = 33.7490
buf_dis = 60 # radius of the round buffer, in km
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


for year in range(start_yr,end_yr+1):
    # create directories
    if not os.path.exists(str(year)):
        os.mkdir(str(year))
    
    filename_wind = '../'+str(year)+'/wind_drct_'+city+'.txt'
    filename_upwind = str(year)+'/upwind_'+city+str(buf_dis)+'.txt'
    filename_downwind = str(year)+'/downwind_'+city+str(buf_dis)+'.txt'
    
    # read wind direction
    wind = pd.read_csv(filename_wind,delimiter = ' ',header = None)
    wind.columns = ['time','direction','speed']
    wind = wind.astype({'time':'str'})
    
    # QPE precipitation data
    nc_file = glob('../'+str(year)+'/nc_file/'+str(year)+'0[6-8][0-3][0-9]*.nc')
    nc_file = sorted(nc_file)
    
    # stats output initialization
    header = "time mean median min max std sum percentile(25) percentile(75) percentile(10) percentile(20) percentile(30) percentile(40) xgrids(90) xgrids(95)"
    try:
        os.remove(filename_upwind)  # remove the file if it exists
    except OSError:
        pass
 
    try:
        os.remove(filename_downwind)
    except OSError:
        pass   
    
    with open(filename_upwind,'a') as fu:
        fu.write(header + '\n')
        fu.close

    with open(filename_downwind,'a') as fd:
        fd.write(header + '\n')
        fd.close
    
    buffer = geodesic_point_buffer(lat,lon,buf_dis) # create buffer
    
    for item in nc_file:
        nc_data = xr.open_dataset(item)
        raw_buffer = nc_data[list(nc_data.keys())[0]].rio.set_crs('WGS84').rio.clip([buffer])
        if raw_buffer.mean()<=0:
            continue
        raw_buffer = (raw_buffer-raw_buffer.mean())/raw_buffer.std() # standardization
        
    ##extract time
        time = item.split('/')[-1]
        time = time[:4]+'-'+time[4:6]+'-'+time[6:8]+'T'+time[9:11]
        print('%s is being processed...'%time)
    
    ##get wind direction value
        wind_val = wind[wind['time']==time]['direction'].values
    
    ##create array of angle
        center_x = int(np.floor(raw_buffer.longitude.size/2))
        center_y = int(np.floor(raw_buffer.latitude.size/2))
        xmat = np.array([[x-raw_buffer.longitude[center_x].values for x in raw_buffer.longitude.values]]*len(raw_buffer.latitude.values))
        ymat = np.array([[y-raw_buffer.latitude[center_y].values for y in raw_buffer.latitude.values]]*len(raw_buffer.longitude.values)).T
        arr_ang = np.rad2deg(np.arctan2(ymat,xmat))
        arr_ang = np.array(arr_ang<0)*360+arr_ang
    
        ang_buffer = xr.DataArray(
             data = arr_ang,
             dims = ['latitude','longitude'], # control the x,y layout when plotting
             coords = dict(
               longitude = (['longitude'],raw_buffer.longitude.values),
               latitude = (['latitude'],raw_buffer.latitude.values),
               )
            )
    
    ## define upwind and downwind areas, now 90 deg sector each'''
    # upwind: wind_val-ang_range/2 to wind_val+ang_range/2; downwind:wind_val-ang_range/2+180 to wind_val+ang_range/2+180
        upwnd_low = wind_val-ang_range/2
        upwnd_high = wind_val+ang_range/2
        dnwnd_low = upwnd_low+180
        dnwnd_high = upwnd_high+180
        if upwnd_low >=0 and upwnd_high <360:
            upwnd_masked = raw_buffer.where((ang_buffer>=upwnd_low)&(ang_buffer<=upwnd_high))
            if dnwnd_high < 360:
                dnwnd_masked = raw_buffer.where((ang_buffer>=dnwnd_low)&(ang_buffer<=dnwnd_high))
            else: #dnwnd_high >= 360:
                if dnwnd_low < 360:
                    dnwnd_masked = raw_buffer.where((ang_buffer>=dnwnd_low)|(ang_buffer<=dnwnd_high%360))
                else:
                    dnwnd_masked = raw_buffer.where((ang_buffer>=dnwnd_low%360)&(ang_buffer<=dnwnd_high%360))
        if upwnd_low >= 0 and upwnd_high >=360:
            upwnd_masked = raw_buffer.where((ang_buffer>=upwnd_low)|(ang_buffer<=upwnd_high%360))
            dnwnd_masked = raw_buffer.where((ang_buffer>=dnwnd_low%360)&(ang_buffer<=dnwnd_high%360))
        if upwnd_low < 0 and upwnd_high > 0:
            upwnd_masked = raw_buffer.where((ang_buffer>=upwnd_low+360)|(ang_buffer<=upwnd_high))
            dnwnd_masked = raw_buffer.where((ang_buffer>=dnwnd_low)&(ang_buffer<=dnwnd_high))
    
     ##############################################
     ### threshold for extreme rainfall ###
    
     ### filtering ###

    #############################################
    
        stats_upwnd = [upwnd_masked.mean(skipna=True).values,
                   upwnd_masked.median(skipna=True).values,
                   upwnd_masked.min(skipna=True).values,
                   upwnd_masked.max(skipna=True).values,
                   upwnd_masked.std(skipna=True).values,
                   upwnd_masked.sum(skipna=True).values,
                   upwnd_masked.quantile(0.25,skipna=True).values,
                   upwnd_masked.quantile(0.75,skipna=True).values,
                   upwnd_masked.quantile(0.1,skipna=True).values,
                   upwnd_masked.quantile(0.2,skipna=True).values,
                   upwnd_masked.quantile(0.3,skipna=True).values,
                   upwnd_masked.quantile(0.4,skipna=True).values]
        with open(filename_upwind,'a') as fu:
                fu.write('%s %s %s %s %s %s %s %s %s %s %s %s %s\n' % 
                         (time,
                          stats_upwnd[0],
                          stats_upwnd[1],
                          stats_upwnd[2],
                          stats_upwnd[3],
                          stats_upwnd[4],
                          stats_upwnd[5],
                          stats_upwnd[6],
                          stats_upwnd[7],
                          stats_upwnd[8],
                          stats_upwnd[9],
                          stats_upwnd[10],
                          stats_upwnd[11]))
                fu.write('\n')
        stats_dnwnd = [dnwnd_masked.mean(skipna=True).values,
                   dnwnd_masked.median(skipna=True).values,
                   dnwnd_masked.min(skipna=True).values,
                   dnwnd_masked.max(skipna=True).values,
                   dnwnd_masked.std(skipna=True).values,
                   dnwnd_masked.sum(skipna=True).values,
                   dnwnd_masked.quantile(0.25,skipna=True).values,
                   dnwnd_masked.quantile(0.75,skipna=True).values,
                   dnwnd_masked.quantile(0.1,skipna=True).values,
                   dnwnd_masked.quantile(0.2,skipna=True).values,
                   dnwnd_masked.quantile(0.3,skipna=True).values,
                   dnwnd_masked.quantile(0.4,skipna=True).values]
        with open(filename_downwind,'a') as fd:
                fd.write('%s %s %s %s %s %s %s %s %s %s %s %s %s\n' % 
                         (time,
                          stats_dnwnd[0],
                          stats_dnwnd[1],
                          stats_dnwnd[2],
                          stats_dnwnd[3],
                          stats_dnwnd[4],
                          stats_dnwnd[5],
                          stats_dnwnd[6],
                          stats_dnwnd[7],
                          stats_dnwnd[8],
                          stats_dnwnd[9],
                          stats_dnwnd[10],
                          stats_dnwnd[11]))
                fd.write('\n')

