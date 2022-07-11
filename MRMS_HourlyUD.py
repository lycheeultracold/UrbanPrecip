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
#from windrose import WindroseAxes

#os.environ['GDAL_DATA'] = os.environ['CONDA_PREFIX'] + r'\Library\share\gdal'
#os.environ['PROJ_LIB'] = os.environ['CONDA_PREFIX'] + r'\Library\share'


# In[ ]:


## input params
start_yr = 2015
end_yr = 2021
city = 'CITY NAME'
lon = # LONGITUDE 
lat = # LATITUDE
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
    filename_upwind = './MRMS_PrecipData/'+str(year)+'/Results/upwind_'+city+str(buf_dis)+'.txt'
    filename_downwind = './MRMS_PrecipData/'+str(year)+'/Results/downwind_'+city+str(buf_dis)+'.txt'
    
    # read wind direction
    wind = pd.read_csv(filename_wind,delimiter = ' ',header = None)
    wind.columns = ['time','direction','speed']
    wind = wind.astype({'time':'str'})
    
    # QPE precipitation data
    nc_file = glob('./MRMS_PrecipData/'+str(year)+'/nc_file/'+str(year)+'0[6-8][0-3][0-9]*.nc')
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
    
    # domain stats
    if not os.path.exists('./MRMS_PrecipData/'+str(year)+'/Results/domain stat'):
        os.mkdir('./MRMS_PrecipData/'+str(year)+'/Results/domain stat')
    
    filename_domain = './MRMS_PrecipData/'+str(year)+'/Results/domain stat/'+city+str(year)+'.txt'
    header = "time domain_avg domain_std"
    try:
        os.remove(filename_domain)  # remove the file if it exists
    except OSError:
        pass  
    
    with open(filename_domain,'a') as f:
        f.write(header + '\n')
        f.close
    
    for item in nc_file:
    ##extract time
        time = item.split('/')[-1]
        time = time[:4]+'-'+time[4:6]+'-'+time[6:8]+'T'+time[9:11]
        print('%s is being processed...'%time)
    
    ##get wind direction value
        wind_val = wind[wind['time']==time]['direction'].values
    
    ## open 60km circle buffer
        nc_data = xr.open_dataset(item)
        raw_buffer = nc_data[list(nc_data.keys())[0]].rio.set_crs('WGS84').rio.clip([buffer])
        #print(raw_buffer.mean(skipna=True).values)
        
        with open(filename_domain,'a') as f:
            f.write('%s %s %s\n' % (time,raw_buffer.mean(skipna=True).values,raw_buffer.std(skipna=True).values))
    
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
        
        # remove drizzels
        upwnd_thres = upwnd_masked.where(upwnd_masked>0)
        dnwnd_thres = dnwnd_masked.where(dnwnd_masked>0)
    
     ##############################################
     ### threshold for extreme rainfall ###
        upwnd_xrain_90 = upwnd_thres.quantile(0.9,skipna=True).values
        dnwnd_xrain_90 = dnwnd_thres.quantile(0.9,skipna=True).values
        upwnd_xrain_95 = upwnd_thres.quantile(0.95,skipna=True).values
        dnwnd_xrain_95 = dnwnd_thres.quantile(0.95,skipna=True).values
    
     ### filtering ###
        upwnd_xrain90 = upwnd_thres.where(upwnd_masked>upwnd_xrain_90)
        dnwnd_xrain90 = dnwnd_thres.where(dnwnd_masked>dnwnd_xrain_90)
        upwnd_xrain95 = upwnd_thres.where(upwnd_masked>upwnd_xrain_95)
        dnwnd_xrain95 = dnwnd_thres.where(dnwnd_masked>dnwnd_xrain_95)
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
                   upwnd_masked.quantile(0.4,skipna=True).values,
                   upwnd_xrain90.count().values,
                   upwnd_xrain95.count().values]
        with open(filename_upwind,'a') as fu:
                fu.write('%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n' % 
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
                          stats_upwnd[11],
                          stats_upwnd[12],
                          stats_upwnd[13],))
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
                   dnwnd_masked.quantile(0.4,skipna=True).values,
                   dnwnd_xrain90.count().values,
                   dnwnd_xrain95.count().values]
        with open(filename_downwind,'a') as fd:
                fd.write('%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n' % 
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
                          stats_dnwnd[11],
                          stats_dnwnd[12],
                          stats_dnwnd[13],))
                fd.write('\n')

