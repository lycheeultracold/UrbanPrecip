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
from mpl_toolkits.axes_grid1 import make_axes_locatable

def geodesic_point_buffer(lat, lon, km):
    # Azimuthal equidistant projection
    aeqd_proj = f"+proj=aeqd +lat_0={lat} +lon_0={lon} +x_0=0 +y_0=0"
    transformer = Transformer.from_crs("EPSG:4326", aeqd_proj, always_xy=True)
    point = Point(lon, lat)
    point_aeqd = transform(transformer.transform, point)
    circle_aeqd = point_aeqd.buffer(km * 1000)
    return mapping(transform(partial(transformer.transform, direction="INVERSE"), circle_aeqd))

## input params
year = 2015
city = 'CITY NAME'
lon = 
lat = 
buf_dis = 60 # radius of the round buffer, in km
os.chdir('/glade/scratch/USERNAME/MRMS_PrecipData/'+str(year))

path = 'Results'
if not os.path.exists(path):
    os.mkdir(path)

filename_wind = 'wind_drct_'+city+'.txt'
filename_upwind = path+'/upwind_'+city+str(buf_dis)+'.txt'
filename_downwind = path+'/downwind_'+city+str(buf_dis)+'.txt'

## QPE precipitation data
nc_file = glob('nc_file/'+str(year)+'0[6-8][0-3][0-9]*.nc')
nc_file = sorted(nc_file)

from scipy import ndimage

file_cen = path+'/'+city+'_cen'+'.txt'
header = "Time X Y Loc"
try:
    os.remove(file_cen)
except OSError:
    pass

# read wind direction
wind = pd.read_csv(filename_wind,delimiter = ' ',header = None)
wind.columns = ['time','direction','speed']
wind = wind.astype({'time':'str'})
    
# create file for storing centroids
with open(file_cen,'a') as f:
    f.write(header + '\n')
    f.close

buffer = geodesic_point_buffer(lat,lon,buf_dis) # create buffer
#cen_ls = [] # used to store coordinates of the center
for item in nc_file:
    time = item.split('/')[1]
    time = time[:4]+'-'+time[4:6]+'-'+time[6:8]+'T'+time[9:11]
    print('%s is being processed...'%time)
    

    '''get wind direction value'''
    wind_val = wind[wind['time']==time]['direction'].values
    wind_val = wind_val%360 # make 0<=wind direction<360
    
    nc_data = xr.open_dataset(item)
    raw_buffer = nc_data['GaugeCorrQPE01H_0mabovemeansealevel'].rio.set_crs('WGS84').rio.clip([buffer])
    raw_buffer = raw_buffer.fillna(0)
    cen = ndimage.measurements.center_of_mass(raw_buffer.values)
    
    center_x = raw_buffer.longitude.size/2
    center_y = raw_buffer.latitude.size/2
    
    ang = np.rad2deg(np.arctan2(cen[1]-center_y,cen[2]-center_x))
    ang = ang if ang>=0 else ang+360
    
    '''decide if the centroid is in upwind or downwind region or neither'''
    '''define upwind and downwind areas, now 90 deg sector each'''
    ang_range = 90 # supposedly <= 90
    # upwind: wind_val-ang_range/2 to wind_val+ang_range/2; downwind:wind_val-ang_range/2+180 to wind_val+ang_range/2+180
    upwnd_low = wind_val-ang_range/2
    upwnd_high = wind_val+ang_range/2
    dnwnd_low = upwnd_low+180
    dnwnd_high = upwnd_high+180
    if upwnd_low >=0 and upwnd_high <360:
        if dnwnd_high < 360:
            if upwnd_low<=ang<=upwnd_high:
                region = 'upwind'
            elif dnwnd_low<=ang<=dnwnd_high:
                region = 'downwind'
            else:
                region = 'neither'
        else: #dnwnd_high >= 360:
            if dnwnd_low < 360:
                if upwnd_low<=ang<=upwnd_high:
                    region = 'upwind'
                elif (dnwnd_low<=ang) | (ang<=dnwnd_high%360):
                    region = 'downwind'
                else:
                    region = 'neither'
            else:
                if upwnd_low<=ang<=upwnd_high:
                    region = 'upwind'
                elif dnwnd_low%360 <= ang <= dnwnd_high%360:
                    region = 'downwind'
                else:
                    region = 'neither'
    if upwnd_low >= 0 and upwnd_high >=360:
        if (upwnd_low<=ang) | (upwnd_high%360>=ang):
            region = 'upwind'
        elif dnwnd_low%360 <= ang <= dnwnd_high%360:
            region = 'downwind'
        else:
            region = 'neither'
    if upwnd_low < 0 and upwnd_high > 0:
        if (upwnd_low<=ang-360) | (upwnd_high>=ang):
            region = 'upwind'
        elif dnwnd_low <= ang <= dnwnd_high:
            region = 'downwind'
        else:
            region = 'neither'
    
    with open(file_cen,'a') as f:
                f.write('%s %s %s %s\n' % 
                         (time,
                          cen[2],
                          cen[1],
                          region))
                f.write('\n')

## open generated results for plotting
file_cen = path+'/'+city+'_cen'+'.txt'
centroids = pd.read_csv(file_cen,delimiter = ' ')
cen_upwnd = centroids.loc[centroids['Loc']=='upwind']
cen_dnwnd = centroids.loc[centroids['Loc']=='downwind']
cen_null = centroids.loc[centroids['Loc']=='neither']