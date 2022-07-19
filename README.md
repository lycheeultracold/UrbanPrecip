# UrbanPrecip

-------

*Scripts for reference to do urban precipitation analysis.*

* MCD12Q1.006_LC_Type2_doy2018001_aid0001.tif: Land use data
* MRMS_Anal_CenOfMass.py: calculate center of mass
* MRMS_Anal_StatTest.py: for stat tests
* MRMS_Domain.py: domain avg and std for precip at each hour
* MRMS_HourlyUD.py: hourly precip metrics at upwind and downwind regions
* MRMS_HourlyUD_Stan.py: standardized results of above
* MRMS_Spatial.py: to generate spatial maps with all hourly maps overlaid
* MRMS_Spatial_WndQuantile.py: spatial maps above distinguished by wind speed quantile values
* MRMS_SpNorm_WndQuantile.py: standardized results of above

* MODIS_LULC.ipynb: mostly for analysis while LULC is involved. check comments within cells.
* Wind_related.ipynb: wind related analysis, e.g. speed, direction. check comments.
* ring_increment.ipynb: overall ring-based results

* shell files: for preprocessing use


**MRMS_V2** folder:  changes applied after changing the way to distinguish hours with different wind directions.

Results/Modis/
Drct_all_dn10: 
Drct_all_up10: 
Total number of grids in quarter circle of radius 10 in the down (dn) and upwind diretions, for eight different wind directions.

Drct_imp_dn10:
Drct_imp_up10:
Total number of impervious grids in quarter circle of radius 10 in the down (dn) and upwind diretions, for eight different wind directions.

Domain_imp10:
Impervious fraction of circular region of radiums 10 (# of impervious grids/ # of grids in circular region) 

Results/Modis/Ring_Result

Atlanta_RingImp:
Impervious fraction from r to r+dr ring.
First row (distance to the city center) -ve upwind; +ve downwind

Results/RingPrecip

Atlanta_RingPrecip
Standardized precipitation 
[standardized using 60 km radius spatial mean and std, i.e., statstics in Results/Directional Hourly Map/Atlanta_Domain60.xsl]
First row (distance to the city center) -ve upwind; +ve downwind

Results/Directional Hourly Map/Atlanta

AtlantaMap_WD*: WD* = E, N, NE, NW, S, SE, SW,W 
Circular region of 60 km radius for raw precip 
summed over total number of hours with corresponding WD [need to get the number of hours with that WD]

AtlantaMap10_WD*: 
Circular region of 10 km radius for raw precip 
summed over total number of hours with corresponding WD 

Atlanta_Domain60, xsl: 
Spatial mean and Std of raw precip 
summed over total number of hours with corresponding WD in circular region of 60 km radius. i.e., sptial mean and std for AtlantaMap_WD*

/Results/
8PrevioDrcts.xsl: Wind direction in decreasing number of hours in each city
DEI_AllDrct90.xsl: DEI computed from quarter circle in eight directions
DEI_Dom_NonDom.xsl: DEI for dominant and non-dominant wind directions
