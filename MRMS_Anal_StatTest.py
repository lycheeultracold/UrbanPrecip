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
#from netCDF4 import Dataset
from functools import partial
from pyproj import CRS
from pyproj import Transformer
from shapely.geometry import Point, mapping
from shapely.ops import transform
#from windrose import WindroseAxes
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy import stats

##### EXAMPLE #####
## All directions

os.chdir('/glade/scratch/USERNAME/MRMS_PrecipData/')
start_yr = 2015
end_yr = 2021
ls_city = ['Atlanta','Austin','Charlotte','Cincinnati','Columbus','Dallas','Houston','Indianapolis','KC','Louisville',
           'Memphis','Miami','Minneapolis','Nashville','NYC','OKC','Omaha','Orlando','Phoenix','Pittsburgh',
           'Richmond','SanAntonio','StLouis','Tucson','Philadelphia','DC','Baltimore']
ls_lon = [-84.3880,-97.7431,-80.8431,-84.5120,-82.9988,-96.7970,-95.3698,-86.1581,-94.5786,-85.7585,-90.0490,
       -80.1918,-93.2650,-86.7816,-74.0060,-97.5164,-95.9345,-81.3789,-112.0740,-79.9959,-77.4360,-98.4936,-90.1994,
          -110.9747,-75.1652,-77.0369,-76.6122]
ls_lat = [33.7490,30.2672,35.2271,39.1031,39.9612,32.7767,29.7604,39.7684,39.0997,38.2527,35.1495,
       25.7617,44.9778,36.1627,40.7128,35.4676,41.2565,28.5384,33.4484,40.4406,37.5407,29.4241,38.6270,32.2226,39.9526,
         38.9072,39.2904]
buf_dis = 60 # radius of the round buffer, in km

ls_ks = []
ls_u = []
ls_size = []
for city in ls_city:
    df_upwnd = pd.DataFrame()
    df_dnwnd = pd.DataFrame()
    dat_wind = pd.DataFrame()
    df_domain = pd.DataFrame()
    for year in range(start_yr,end_yr+1):
        os.chdir('/glade/scratch/USERNAME/MRMS_PrecipData/'+str(year))
        filename_upwind = 'Results/upwind_'+city+str(buf_dis)+'.txt'
        filename_downwind = 'Results/downwind_'+city+str(buf_dis)+'.txt'
        filename_wind = 'wind_drct_'+city+'.txt'
        dat_precip_up = pd.read_csv(filename_upwind,delimiter=' ')
        dat_precip_dn = pd.read_csv(filename_downwind,delimiter=' ')
        wind = pd.read_csv(filename_wind,delimiter = ' ',header = None)
        df_upwnd = pd.concat([df_upwnd,dat_precip_up],ignore_index=True)
        df_dnwnd = pd.concat([df_dnwnd,dat_precip_dn],ignore_index=True)
        dat_wind = pd.concat([dat_wind,wind],ignore_index=True)
        #open domain stats
        file_domain = 'Results/domain stat/'+city+str(year)+'_'+str(buf_dis)+'.txt'
        domain = pd.read_csv(file_domain,sep=' ')
        df_domain = pd.concat([df_domain,domain],ignore_index=True)
    #print(df_domain)    
    dat_wind.columns = ['time','direction','speed']
    dat_wind['zone'] = dat_wind['direction'].map(wnd_dmnt)
    zone_count = dat_wind['zone'].value_counts()
    
    colnames = ["time","mean","median","min","max","std","sum","percentile(25)","percentile(75)","percentile(10)",
            "percentile(20)","percentile(30)","percentile(40)","xgrids(90)","xgrids(95)"]
    df_upwnd.columns = colnames
    df_upwnd['region'] = 'upwind'
    df_dnwnd.columns = colnames
    df_dnwnd['region'] = 'downwind'

    no_dup_up = pd.concat([dat_wind.time,df_upwnd.time]).drop_duplicates(keep=False)
    dat_wind = dat_wind.drop(no_dup_up.index.values)
    df_upwnd = dat_wind.merge(df_upwnd,how='outer',on='time')
    df_dnwnd = dat_wind.merge(df_dnwnd,how='outer',on='time')

    # filter
    df_upwnd = df_upwnd[0.1<=df_upwnd['mean']]
    df_dnwnd = df_dnwnd[0.1<=df_dnwnd['mean']]

    intersect = set(df_upwnd['time'])&set(df_dnwnd['time'])
    df_upwnd = df_upwnd[df_upwnd['time'].isin(intersect)]
    df_dnwnd = df_dnwnd[df_dnwnd['time'].isin(intersect)]
    df_domain = df_domain[df_domain['time'].isin(intersect)]
    df_upwnd = df_upwnd.merge(df_domain,how='outer',on='time')
    df_dnwnd = df_dnwnd.merge(df_domain,how='outer',on='time')

    df_upwnd['norm_mean'] = (df_upwnd['mean']-df_upwnd['domain_avg'])/df_upwnd['domain_std']
    df_dnwnd['norm_mean'] = (df_dnwnd['mean']-df_dnwnd['domain_avg'])/df_dnwnd['domain_std']
    
    ks_mean = stats.ks_2samp(df_upwnd['norm_mean'],df_dnwnd['norm_mean'])[1]
    ls_ks.append(round(ks_mean,2))
    u_mean = stats.mannwhitneyu(df_upwnd['norm_mean'],df_dnwnd['norm_mean'])[1]
    ls_u.append(round(u_mean,2))
    
    #sample size
    ls_size.append(df_upwnd.shape[0])
    
print(ls_ks)
print(ls_u)

res_p = pd.DataFrame(
    {
        'City':ls_city,
        'Longitude':ls_lon,
        'Latitude':ls_lat,
        'KS-P':ls_ks,
        'U-P':ls_u,
        'Sample size':ls_size
    } 
)
os.chdir('')
res_p.to_csv('AllDrct_SigLevel'+str(buf_dis)+'.csv',sep=' ')