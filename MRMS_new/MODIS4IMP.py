#!/usr/bin/env python
# coding: utf-8

# In[1]:


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
from scipy.optimize import curve_fit
from sympy import symbols, solve


# In[2]:


# Drawing km circle around what we're interested in
def geodesic_point_buffer(lat, lon, km):
    # Azimuthal equidistant projection
    aeqd_proj = f"+proj=aeqd +lat_0={lat} +lon_0={lon} +x_0=0 +y_0=0"
    transformer = Transformer.from_crs("EPSG:4326", aeqd_proj, always_xy=True)
    point = Point(lon, lat)
    point_aeqd = transform(transformer.transform, point)
    circle_aeqd = point_aeqd.buffer(km * 1000)
    return mapping(transform(partial(transformer.transform, direction="INVERSE"), circle_aeqd))
    
def drct2val(x): # angle range is 45
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
    
def drct2val_90(x): # angle range is 90
    # within list: [upwind lower, upwind upper, downwind lower, downwind upper]
    if x == 'E':
        return [315,45,135,225]
    if x == 'NE':
        return [0,90,180,270]
    if x == 'N':
        return [45,135,225,315]
    if x == 'NW':
        return [90,180,270,360]
    if x == 'W':
        return [135,225,315,45]
    if x == 'SW':
        return [180,270,0,90]
    if x == 'S':
        return [225,315,45,135]
    if x == 'SE':
        return [270,360,90,180]
    
def com_row(df1,df2,filt_col,ext_col,val): # get common rows of two dataframes
    # filt_col -- the column used to set conditions for filtering, e.g.'mean'
    # val -- conditional value as the filter
    # ext_col -- the column used to extract common rows, e.g.'time'
    df1 = df1[df1[filt_col]>=val]
    df2 = df2[df2[filt_col]>=val]
    intersect = set(df1[ext_col])&set(df2[ext_col])
    df1 = df1[df1[ext_col].isin(intersect)]
    df2 = df2[df2[ext_col].isin(intersect)]
    return df1,df2

def ang_mat(basemat):
    ## create angle matrix 
    ## basemat need to be in the form opened by rioxarray.open_rasterio 
    
    ## be careful with the dimension name in any opened tif, they need to be consistent!
    center_x = int(np.floor(basemat.x.size/2))
    center_y = int(np.floor(basemat.y.size/2))
    xmat = np.array([[x-basemat.x[center_x].values for x in basemat.x.values]]*len(basemat.y.values))
    ymat = np.array([[y-basemat.y[center_y].values for y in basemat.y.values]]*len(basemat.x.values)).T
    arr_ang = np.rad2deg(np.arctan2(ymat,xmat))
    arr_ang = np.array(arr_ang<0)*360+arr_ang
        
    ang_buffer = xr.DataArray(
             data = arr_ang,
             dims = ['y','x'], # control the x,y layout when plotting
             coords = dict(
               x = (['x'],basemat.x.values),
               y = (['y'],basemat.y.values),
               )
            )
    return ang_buffer

# function of calculating imp
def drct_px(city,lon,lat,ang_range):    
    
    file_land = glob('/glade/scratch/USERNAME/Modis/*.tif')
    file_land = file_land[0]
    #year = file_land.split('_')[-2][3:7]

    ras = rxr.open_rasterio(file_land)
    buffer = geodesic_point_buffer(lat,lon,buf_dis) # create buffer
    raw_buffer = ras.rio.set_crs('WGS84').rio.clip([buffer])
    raw_buffer = raw_buffer[0]
    raw_buffer = raw_buffer.where(raw_buffer < 255)
    px_sum = raw_buffer.count().values
    px_imp = raw_buffer.where(raw_buffer == 13).count().values
    frac_imp = px_imp/px_sum
    frac_imp = round(frac_imp,4)

    center_x = int(np.floor(raw_buffer.x.size/2))
    center_y = int(np.floor(raw_buffer.y.size/2))
    xmat = np.array([[x-raw_buffer.x[center_x-1].values for x in raw_buffer.x.values]]*len(raw_buffer.y.values))
    ymat = np.array([[y-raw_buffer.y[center_y-1].values for y in raw_buffer.y.values]]*len(raw_buffer.x.values)).T
    arr_ang = np.rad2deg(np.arctan2(ymat,xmat))
    arr_ang = np.array(arr_ang<0)*360+arr_ang
    
    ang_buffer = xr.DataArray(
            data = arr_ang,
            dims = ['y','x'], # control the x,y layout when plotting
            coords = dict(
            x = (['x'],raw_buffer.x.values),
            y = (['y'],raw_buffer.y.values),
            )
        )

    tmp_all_up = [city]
    tmp_imp_up = [city]
    tmp_all_dn = [city]
    tmp_imp_dn = [city]
    wind_val = [0,45,90,135,180,225,270,315]
    for i in wind_val:
        bnd_low_up = i-ang_range/2
        bnd_high_up = i+ang_range/2
        bnd_low_dn = i-ang_range/2+180
        bnd_high_dn = i+ang_range/2+180
    
        ## upwind
        if (bnd_low_up>=0) & (bnd_high_up<=360):
            masked_up = raw_buffer.where((ang_buffer>=bnd_low_up)&(ang_buffer<=bnd_high_up))
        if (bnd_low_up<0) & (bnd_high_up>0):
            masked_up = raw_buffer.where((ang_buffer>=bnd_low_up+360)|(ang_buffer<=bnd_high_up))
        ## downwind
        if (bnd_low_dn%360) < (bnd_high_dn%360):
            masked_dn = raw_buffer.where((ang_buffer>=bnd_low_dn%360)&(ang_buffer<=bnd_high_dn%360))
        if (bnd_low_dn%360) > (bnd_high_dn%360):
            masked_dn = raw_buffer.where((ang_buffer>=bnd_low_dn%360)|(ang_buffer<=bnd_high_dn%360))
        
        px_sum_up = masked_up.count().values
        px_imp_up = masked_up.where(masked_up == 13).count().values
        tmp_all_up.append(px_sum_up)
        tmp_imp_up.append(px_imp_up)
        px_sum_dn = masked_dn.count().values
        px_imp_dn = masked_dn.where(masked_dn == 13).count().values
        tmp_all_dn.append(px_sum_dn)
        tmp_imp_dn.append(px_imp_dn)
    tmp_all_up = [tmp_all_up]
    tmp_imp_up = [tmp_imp_up]
    tmp_all_dn = [tmp_all_dn]
    tmp_imp_dn = [tmp_imp_dn]
    return tmp_all_up,tmp_imp_up,tmp_all_dn,tmp_imp_dn,[[city,frac_imp]]

def rev_drct(x):
    drct = ['E','NE','N','NW','W','SW','S','SE']
    idx = drct.index(x)
    return drct[(idx+4)%len(drct)]


# In[3]:


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

ang_range = 90

df_city = pd.DataFrame({'Name':ls_city,'lon':ls_lon,'lat':ls_lat})
df_city = df_city.sort_values('Name',ignore_index=True)
print(df_city.head(6))


# In[ ]:


if not os.path.exists('/glade/scratch/USERNAME/MRMS_new/Results/Modis'):
    os.mkdir('/glade/scratch/USERNAME/MRMS_new/Results/Modis')
    
os.chdir('/glade/scratch/USERNAME/MRMS_new/Results/Modis')

# calculate number of imp pixels and total pixels with different radius
for buf_dis in np.arange(10,70,10):
    px_all_up = pd.DataFrame()
    px_imp_up = pd.DataFrame()
    px_all_dn = pd.DataFrame()
    px_imp_dn = pd.DataFrame()
    imp_domain = pd.DataFrame()
    for city,lon,lat in zip(ls_city,ls_lon,ls_lat):
        res = drct_px(city,lon,lat,90)
        px_all_up = px_all_up.append(res[0],ignore_index=True)
        px_imp_up = px_imp_up.append(res[1],ignore_index=True)
        px_all_dn = px_all_dn.append(res[2],ignore_index=True)
        px_imp_dn = px_imp_dn.append(res[3],ignore_index=True)
        imp_domain = imp_domain.append(res[4],ignore_index=True) # overall fraction of imp (non-directional)
    col_name = ['City','E','NE','N','NW','W','SW','S','SE']
    px_all_up.columns = col_name
    px_imp_up.columns = col_name
    px_all_dn.columns = col_name
    px_imp_dn.columns = col_name
    imp_domain.columns = ['City','imp_fraction']
    
    px_all_up.to_csv('Drct_all_up'+str(buf_dis).zfill(2)+'.csv',sep = ' ')
    px_imp_up.to_csv('Drct_imp_up'+str(buf_dis).zfill(2)+'.csv',sep = ' ')
    px_all_dn.to_csv('Drct_all_dn'+str(buf_dis).zfill(2)+'.csv',sep = ' ')
    px_imp_dn.to_csv('Drct_imp_dn'+str(buf_dis).zfill(2)+'.csv',sep = ' ')
    imp_domain.to_csv('Domain_imp'+str(buf_dis).zfill(2)+'.csv',sep = ' ')


# In[24]:


# directional and ring-based imp calculation

os.chdir('/glade/scratch/USERNAME/MRMS_new/Results/Modis')

file_allpx = glob("Drct_all_up*")
file_imppx = glob("Drct_imp_up*")

df_drct = pd.read_csv('../8PrevailDrcts.csv',sep=' ')

for city in df_city['Name']:
    df_res = pd.DataFrame()
    
    for i in range(1,9):
        drct = df_drct[df_drct['City']==city][str(i)].values[0]
        
        ls_all_up = []
        ls_imp_up = []
        ls_all_dn = []
        ls_imp_dn = []
        for f_all,f_imp in zip(sorted(file_allpx),sorted(file_imppx)):
            
            df_all = pd.read_csv(f_all,sep = ' ')
            df_imp = pd.read_csv(f_imp,sep = ' ')
            
            ls_all_up.append(df_all[df_all['City']==city][drct].values[0])
            ls_imp_up.append(df_imp[df_imp['City']==city][drct].values[0])
            ls_all_dn.append(df_all[df_all['City']==city][rev_drct(drct)].values[0])
            ls_imp_dn.append(df_imp[df_imp['City']==city][rev_drct(drct)].values[0])
            
        imp_rate_up = [ls_imp_up[0]/ls_all_up[0]]
        
        for k in range(1,len(file_allpx)):
            imp_rate_up.append((ls_imp_up[k]-ls_imp_up[k-1])/(ls_all_up[k]-ls_all_up[k-1]))
        imp_rate_up = imp_rate_up[::-1]
        imp_rate_dn = [ls_imp_dn[0]/ls_all_dn[0]]
        for k in range(1,len(file_allpx)):
            imp_rate_dn.append((ls_imp_dn[k]-ls_imp_dn[k-1])/(ls_all_dn[k]-ls_all_dn[k-1]))
            
        imp_rate = imp_rate_up+imp_rate_dn
        imp_rate = [round(i,4) for i in imp_rate]
        imp_rate = [0 if x<0 else x for x in imp_rate]
        imp_rate = [1 if x>1 else x for x in imp_rate]
        df_res = df_res.append(pd.DataFrame([[drct]+imp_rate]), ignore_index=True)
    colnames = ['Direction','-60_-50','-50_-40','-40_-30','-30_-20','-20_-10','-10_0',
                '0_10','10_20','20_30','30_40','40_50','50_60']
    df_res.columns = colnames
    df_res.to_csv('Ring_Result/'+city+'_RingImp.csv',sep=' ')

