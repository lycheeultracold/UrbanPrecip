# UrbanPrecip

-------

*Scripts for reference to do urban precipitation analysis.*

**MRMS_V2** folder:  changes applied after changing the way to distinguish hours with different wind directions.
(codes outside MRMS_V2 folder are mostly for reference use)

* MCD12Q1.006_LC_Type2_doy2018001_aid0001.tif: Land use data
* MRMS_Anal_CenOfMass.py: calculate center of mass
* MRMS_Anal_StatTest.py: for stat tests
* MRMS_Domain.py: domain avg and std for precip at each hour
* MRMS_HourlyUD.py: hourly precip metrics at upwind and downwind regions
* MRMS_HourlyUD_Stan.py: standardized results of above
* MRMS_Spatial.py: to generate spatial maps with all hourly maps overlaid
* MRMS_Spatial_WndQuantile.py: spatial maps above distinguished by wind speed quantile values
* MRMS_SpNorm_WndQuantile.py: standardized results of above
