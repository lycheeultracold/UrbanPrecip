#!/usr/bin/env python
# coding: utf-8

# # Domain maps of precipitation accumulation at each wind direction hours
# ## i.e. if 2015-06-01 and 2015-06-02 are two hours with the same wind direction, then add them together

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
from windrose import WindroseAxes
from mpl_toolkits.axes_grid1 import make_axes_locatable


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
    
def drct2val(x):
    # within list: [upwind lower, upwind upper, downwind lower, downwind upper]
    if x == 'E':
        return [337.5,22.5,157.5,202.5]
    if x == 'NE':
        return [22.5,67.5,202.5,247.5]
    if x == 'N':
        return [67.5,112.5,247.5,292.5]
    if x == 'NW':
        return [112.5,157.5,292.5,337.5]
    if x == 'W':
        return [157.5,202.5,337.5,22.5]
    if x == 'SW':
        return [202.5,247.5,22.5,67.5]
    if x == 'S':
        return [247.5,292.5,67.5,112.5]
    if x == 'SE':
        return [292.5,337.5,112.5,157.5]


# In[ ]:


# add together precip map at hours with same direction (wind speed: all)

start_yr = 2015
end_yr = 2021
ls_city = ['Atlanta','Austin','Baltimore','Charlotte','Cincinnati','Columbus','Dallas','DC','Houston','Indianapolis','KC','Louisville','Memphis','Miami','Minneapolis','Nashville','NYC','OKC','Omaha','Orlando',
           'Phoenix','Pittsburgh','Richmond','SanAntonio','StLouis','Tucson','Philadelphia']
# 
ls_lon = [-84.3880,-97.7431,-76.6122,-80.8431,-84.5120,-82.9988,-96.7970,-77.0369,-95.3698,-86.1581,-94.5786,-85.7585,-90.0490,-80.1918,-93.2650,-86.7816,-74.0060,-97.5164,-95.9345,
          -81.3789,-112.0740,-79.9959,-77.4360,-98.4936,-90.1994,-110.9747,-75.1652]
#
ls_lat = [33.7490,30.2672,39.2904,35.2271,39.1031,39.9612,32.7767,38.9072,29.7604,39.7684,39.0997,38.2527,35.1495,25.7617,44.9778,36.1627,40.7128,35.4676,41.2565,28.5384,
          33.4484,40.4406,37.5407,29.4241,38.6270,32.2226,39.9526]

all_drcts = ['E','NE','N','NW','W','SW','S','SE'] 



for city in ls_city:
    
    wind = pd.DataFrame() # store wind data for all years
    for year in range(start_yr,end_yr+1):  
        ## get hourly wind data
        file_wind = '../MRMS_PrecipData/'+str(year)+'/wind_drct_'+city+'.txt'
        tmp_wind = pd.read_csv(file_wind,sep = ' ',header = None)
        wind = pd.concat([wind,tmp_wind],ignore_index = True)
    
    wind.columns = ['time','direction','speed']
    wind['zone'] = pd.cut(x=wind['direction'],bins=[0]+list(np.arange(22.5,338,45))+[360],right=False,
      labels=['E','NE','N','NW','W','SW','S','SE','E'],ordered=False)
    
    ls_name = glob('Results/Hourly Map/'+city+'/*.tif')
    ls_name = sorted(ls_name)
    
    root = os.getcwd()
    path = '/Results/Directional Hourly Map/'+city
    new_path = root+path
    if not os.path.exists(new_path):
        os.mkdir(new_path)
    
    for item in all_drcts:
        wind_sel = wind[wind['zone']==item]
        matches = [x for x in ls_name if x.split('/')[-1].split('.')[0].split('_')[1] in wind_sel['time'].values] 
        maps = 0
        for x in matches:
            tmp = rxr.open_rasterio(x)[0]
            tmp = tmp.where(tmp<999)
            maps = maps+tmp
        maps.rio.to_raster(new_path+'/'+city+'Map_'+item+'.tif')

