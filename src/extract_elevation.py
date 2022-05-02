# -*- coding: utf-8 -*-
"""
Created on Mon May  2 16:33:27 2022

@author: jullienn
"""
#Import packages
#import rasterio
#from rasterio.plot import show
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pandas as pd
from os import listdir
from os.path import isfile, join
import pickle
#from pysheds.grid import Grid
import pdb
import numpy as np
from pyproj import Transformer
import matplotlib.gridspec as gridspec
import scipy.io
from osgeo import gdal
import geopandas as gpd  # Requires the pyshp package

from matplotlib.colors import ListedColormap, BoundaryNorm
from shapely.geometry import Point, Polygon
from matplotlib.patches import Patch
import cartopy.crs as ccrs
from matplotlib.lines import Line2D

import seaborn as sns
sns.set_theme(style="whitegrid")
from scalebar import scale_bar
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


#Set fontsize plot
plt.rcParams.update({'font.size': 10})

create_elevation_dictionaries='FALSE'
#pdb.set_trace()

### -------------------------- Load GrIS DEM ----------------------------- ###
#https://towardsdatascience.com/reading-and-visualizing-geotiff-images-with-python-8dcca7a74510
import rasterio
from rasterio.plot import show

path_GrIS_DEM = r'C:/Users/jullienn/switchdrive/Private/research/backup_Aglaja/working_environment/greenland_topo_data/elevations/greenland_dem_mosaic_100m_v3.0.tif'
GrIS_DEM = rasterio.open(path_GrIS_DEM)

#fig, axs = plt.subplots()
#show(img,ax=axs)
### -------------------------- Load GrIS DEM ----------------------------- ###

### -------------------------- Load shapefiles --------------------------- ###
#from https://gis.stackexchange.com/questions/113799/how-to-read-a-shapefile-in-python
path_IceBridgeArea_Shape='C:/Users/jullienn/switchdrive/Private/research/backup_Aglaja/working_environment/greenland_topo_data/IceBridge Area Shapefiles/IceBridge Area Shapefiles/'
IceBridgeArea_Shape=gpd.read_file(path_IceBridgeArea_Shape+'IceBridgeArea_Shape.shp')

path_regional_masks='C:/Users/jullienn/switchdrive/Private/research/backup_Aglaja/working_environment/greenland_topo_data/masks_for_2002_2003_calculations'

NW_icecap_greenland_mask=gpd.read_file(path_regional_masks+'/NW_icecap_greenland_mask_3413.shp')
NW_north_greenland_mask=gpd.read_file(path_regional_masks+'/NW_north_greenland_mask_3413.shp')
NW_west_greenland_mask=gpd.read_file(path_regional_masks+'/NW_west_greenland_mask_3413.shp')
SW_lower_greenland_mask=gpd.read_file(path_regional_masks+'/SW_lower_greenland_mask_3413.shp')
SW_middle_greenland_mask=gpd.read_file(path_regional_masks+'/SW_middle_greenland_mask_3413.shp')
SW_upper_greenland_mask=gpd.read_file(path_regional_masks+'/SW_upper_greenland_mask_3413.shp')

#Load Rignot et al., 2016 Greenland drainage bassins
path_rignotetal2016_GrIS_drainage_bassins='C:/Users/jullienn/switchdrive/Private/research/backup_Aglaja/working_environment/greenland_topo_data/GRE_Basins_IMBIE2_v1.3/'
GrIS_drainage_bassins=gpd.read_file(path_rignotetal2016_GrIS_drainage_bassins+'GRE_Basins_IMBIE2_v1.3_EPSG_3413.shp',rows=slice(51,57,1)) #the regions are the last rows of the shapefile

#Extract indiv regions and create related indiv shapefiles
NO_rignotetal=GrIS_drainage_bassins[GrIS_drainage_bassins.SUBREGION1=='NO']
NE_rignotetal=GrIS_drainage_bassins[GrIS_drainage_bassins.SUBREGION1=='NE']
SE_rignotetal=GrIS_drainage_bassins[GrIS_drainage_bassins.SUBREGION1=='SE']
SW_rignotetal=GrIS_drainage_bassins[GrIS_drainage_bassins.SUBREGION1=='SW']
CW_rignotetal=GrIS_drainage_bassins[GrIS_drainage_bassins.SUBREGION1=='CW']
NW_rignotetal=GrIS_drainage_bassins[GrIS_drainage_bassins.SUBREGION1=='NW']

#Load Rignot et al., 2016 GrIS mask
path_rignotetal2016_GrIS='C:/Users/jullienn/switchdrive/Private/research/backup_Aglaja/working_environment/greenland_topo_data/GRE_IceSheet_IMBIE2/GRE_IceSheet_IMBIE2/'
GrIS_rignotetal2016=gpd.read_file(path_rignotetal2016_GrIS+'GRE_IceSheet_IMBIE2_v1_EPSG3413.shp',rows=slice(1,2,1)) #the regions are the last rows of the shapefile
GrIS_mask=GrIS_rignotetal2016[GrIS_rignotetal2016.SUBREGION1=='ICE_SHEET']
### -------------------------- Load shapefiles --------------------------- ###



################### Load 2002-2003 ice lenses location ##################
path_data='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/icelens_identification'
  
#Open the file and read it
f_icelens_flightlines = open(path_data+'/metadata_coord_icelens_2002_2003_26022020', "rb")
icelens_2002_3_flightlines = pickle.load(f_icelens_flightlines)
f_icelens_flightlines.close()

lat_icelens=[]
lon_icelens=[]
colorcode_icelens=[]
Track_name=[]

for year in list(icelens_2002_3_flightlines.keys()):
    for days in list(icelens_2002_3_flightlines[year].keys()):
        for indiv_file in list(icelens_2002_3_flightlines[year][days].keys()):
            print(indiv_file)
            if (indiv_file[0:7]=='quality'):
                print('Quality file, continue')
                continue
            elif (not(bool(icelens_2002_3_flightlines[year][days][indiv_file]))):
                print('No ice lens, continue')
                continue
            else:
                lat_icelens=np.append(lat_icelens,icelens_2002_3_flightlines[year][days][indiv_file][0])
                lon_icelens=np.append(lon_icelens,icelens_2002_3_flightlines[year][days][indiv_file][1])
                colorcode_icelens=np.append(colorcode_icelens,icelens_2002_3_flightlines[year][days][indiv_file][2])
                #Create an empty vector of strings
                Track_name=np.append(Track_name,[indiv_file for x in range(0,len(icelens_2002_3_flightlines[year][days][indiv_file][0]))])

#Create a dataframe out of it
df_2002_2003=pd.DataFrame(lat_icelens, columns =['lat_3413'])
df_2002_2003['lon_3413']=lon_icelens
df_2002_2003['colorcode_icelens']=colorcode_icelens
df_2002_2003['Track_name']=Track_name
################### Load 2002-2003 ice lenses location ##################

################### Load 2010-2018 ice slabs location ##################
#Load the data
filename_20102018='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2010_2018/excel_spatial_aggreation_and_other/final_excel/high_estimate/Ice_Layer_Output_Thicknesses_2010_2018_jullienetal2021_high_estimate.csv'
df_20102018 = pd.read_csv(filename_20102018, sep=",", decimal='.')

#Transform the coordinated from WGS84 to EPSG:3413
#Example from: https://pyproj4.github.io/pyproj/stable/examples.html
transformer = Transformer.from_crs("EPSG:4326", "EPSG:3413", always_xy=True)
points=transformer.transform(np.array(df_20102018.lon),np.array(df_20102018.lat))

lon_3413_20102018=points[0]
lat_3413_20102018=points[1]
################### Load 2010-2018 ice slabs location ##################

#######################################################################
###          Inland expansion of iceslabs from 2002 to 2018         ###
#######################################################################

# compare min and max of lat/lon of the track with respect to shapefile

#Store lat/lon 3413
df_20102018['lat_3413']=lat_3413_20102018
df_20102018['lon_3413']=lon_3413_20102018

#Initialise the elevation and shapefile belonging column
df_20102018['key_shp']=np.nan
df_20102018['elevation']=np.nan
df_20102018['year']=np.nan

#I. Load regional shapefile I have created on QGIS
#Done before the if statement

#II. Do the intersection between the mask and 2010-2018 data and keep only the matching one

#This part of code is from 'refine_location_2017_2018.py'
#Loop over all data point to check whether it belongs to one of the four shapefile
for i in range(0,lon_3413_20102018.size):
    #select the point i
    single_point=Point(lon_3413_20102018[i],lat_3413_20102018[i])
    
    #Do the identification between the point i and the regional shapefiles
    #From: https://automating-gis-processes.github.io/CSC18/lessons/L4/point-in-polygon.html
    check_NO_rignotetal=np.asarray(NO_rignotetal.contains(single_point)).astype(int)
    check_NE_rignotetal=np.asarray(NE_rignotetal.contains(single_point)).astype(int)
    check_SE_rignotetal=np.asarray(SE_rignotetal.contains(single_point)).astype(int)
    check_SW_rignotetal=np.asarray(SW_rignotetal.contains(single_point)).astype(int)
    check_CW_rignotetal=np.asarray(CW_rignotetal.contains(single_point)).astype(int)
    check_NW_rignotetal=np.asarray(NW_rignotetal.contains(single_point)).astype(int)

    #Associated the point of interest to its regional shapefile in data_iceslabs
    if (np.sum(check_NO_rignotetal)>0):
        df_20102018['key_shp'].iloc[i]='NO'
    elif (np.sum(check_NE_rignotetal)>0):
        df_20102018['key_shp'].iloc[i]='NE'
    elif (np.sum(check_SE_rignotetal)>0):
        df_20102018['key_shp'].iloc[i]='SE'
    elif (np.sum(check_SW_rignotetal)>0):
        df_20102018['key_shp'].iloc[i]='SW'
    elif (np.sum(check_CW_rignotetal)>0):
        df_20102018['key_shp'].iloc[i]='CW'
    elif (np.sum(check_NW_rignotetal)>0):
        df_20102018['key_shp'].iloc[i]='NW'
    else:
        df_20102018['key_shp'].iloc[i]='Out'
    
    #Add the year
    df_20102018['year'].iloc[i]=int(df_20102018['Track_name'].iloc[i][0:4])
    
    #Calcul elevation
    if (np.isnan(df_20102018['lon_3413'].iloc[i])):
        continue
    
    #This is from https://gis.stackexchange.com/questions/190423/getting-pixel-values-at-single-point-using-rasterio
    for val in GrIS_DEM.sample([(df_20102018['lon_3413'].iloc[i], df_20102018['lat_3413'].iloc[i])]): 
        #Calculate the corresponding elevation
        df_20102018['elevation'].iloc[i]=val
    
    #Monitor the process
    print(i/lon_3413_20102018.size*100,'%')

#Save the dictionary into a picke file
filename_tosave='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2010_2018/excel_spatial_aggreation_and_other/df_20102018_with_elevation_high_estimate_rignotetalregions'
outfile= open(filename_tosave, "wb" )
pickle.dump(df_20102018,outfile)
outfile.close()

pdb.set_trace()

#III. Do the intersection between the mask and 2002-2003 data

#Initialise the shapefile belonging column
df_2002_2003['key_shp']=np.nan
df_2002_2003['elevation']=np.nan
df_2002_2003['year']=np.nan

#This part of code is from 'refine_location_2017_2018.py'
#Loop over all data point to check whether it belongs to one of the four shapefile
for i in range(0,len(df_2002_2003)):
    #select the point i
    single_point=Point(df_2002_2003['lon_3413'].iloc[i],df_2002_2003['lat_3413'].iloc[i])
    
    #Do the identification between the point i and the regional shapefiles
    #From: https://automating-gis-processes.github.io/CSC18/lessons/L4/point-in-polygon.html
    check_NO_rignotetal=np.asarray(NO_rignotetal.contains(single_point)).astype(int)
    check_NE_rignotetal=np.asarray(NE_rignotetal.contains(single_point)).astype(int)
    check_SE_rignotetal=np.asarray(SE_rignotetal.contains(single_point)).astype(int)
    check_SW_rignotetal=np.asarray(SW_rignotetal.contains(single_point)).astype(int)
    check_CW_rignotetal=np.asarray(CW_rignotetal.contains(single_point)).astype(int)
    check_NW_rignotetal=np.asarray(NW_rignotetal.contains(single_point)).astype(int)
    
    #Associated the point of interest to its regional shapefile in data_iceslabs
    if (np.sum(check_NO_rignotetal)>0):
        df_2002_2003['key_shp'].iloc[i]='NO'
    elif (np.sum(check_NE_rignotetal)>0):
        df_2002_2003['key_shp'].iloc[i]='NE'
    elif (np.sum(check_SE_rignotetal)>0):
        df_2002_2003['key_shp'].iloc[i]='SE'
    elif (np.sum(check_SW_rignotetal)>0):
        df_2002_2003['key_shp'].iloc[i]='SW'
    elif (np.sum(check_CW_rignotetal)>0):
        df_2002_2003['key_shp'].iloc[i]='CW'
    elif (np.sum(check_NW_rignotetal)>0):
        df_2002_2003['key_shp'].iloc[i]='NW'
    else:
        df_2002_2003['key_shp'].iloc[i]='Out'
    
    #Add the year
    if (df_2002_2003['Track_name'].iloc[i][6:8] == '02'):
        year_to_write=2002
    elif (df_2002_2003['Track_name'].iloc[i][6:8] == '03'):
        year_to_write=2003
    else:
        print('Year not known, error')
        break
    df_2002_2003['year'].iloc[i]=year_to_write
    
    #Calcul elevation
    if (np.isnan(df_2002_2003['lon_3413'].iloc[i])):
        continue
    
    #This is from https://gis.stackexchange.com/questions/190423/getting-pixel-values-at-single-point-using-rasterio
    for val in GrIS_DEM.sample([(df_2002_2003['lon_3413'].iloc[i], df_2002_2003['lat_3413'].iloc[i])]): 
        #Calculate the corresponding elevation
        df_2002_2003['elevation'].iloc[i]=val
    
    #Monitor the process
    print(i/len(df_2002_2003)*100,'%')

#Only work with green slabs
df_2002_2003_green=df_2002_2003[df_2002_2003['colorcode_icelens']==1]

#Save the dictionary into a picke file
filename_tosave='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2010_2018/excel_spatial_aggreation_and_other/df_2002_2003_with_elevation_rignotetalregions'
outfile= open(filename_tosave, "wb" )
pickle.dump(df_2002_2003,outfile)
outfile.close()
