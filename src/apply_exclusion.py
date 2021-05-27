# -*- coding: utf-8 -*-
"""
Created on Wed May 26 10:23:18 2021

@author: JullienN
"""

#In this code, I perform the matching of exclusion lat/lon with data using a 
#certain radius around the exclusion, and I perform the transcription of
#lat/lon exclusion into pixel exclusion

##############################################################################
############### Define function for discrete colorbar display ###############
##############################################################################

def discrete_cmap(N, base_cmap=None):
    """Create an N-bin discrete colormap from the specified input map"""
    #This piece of code is from: https://gist.github.com/jakevdp/91077b0cae40f8f8244a
    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:
    
    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)
##############################################################################
############### Define function for discrete colorbar display ###############
##############################################################################



##############################################################################
################## Define function for distance calculation ##################
##############################################################################
#This function comes from https://stackoverflow.com/questions/4913349/haversine-formula-in-python-bearing-and-distance-between-two-gps-points/4913653#4913653
def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a))
    r = 6371 # Radius of earth in kilometers. Use 3956 for miles
    return c * r
##############################################################################
################## Define function for distance calculation ##################
##############################################################################

import pickle
import pdb
import h5py
import geopandas as gpd
import descartes
from osgeo import gdal
import matplotlib.pyplot as plt
from pysheds.grid import Grid
import matplotlib.colors as mcolors
import numpy as np
from pyproj import Transformer
from shapely.geometry import Point, Polygon
import matplotlib.patches as mpatches
from math import radians, cos, sin, asin, sqrt


#Let's first start with a test data sample to make the matching work
#Looking at QGIS project, I choose the following sequence of traces:
#   - 20170422_01_168_170
#   - 20140424_01_002_004 
#   - 20130405_01_165_167 
#   - 20120418_01_129_131 
#   - 20110419_01_008_010 
#   - 20100508_01_114_115

### --------------- Open and read the exclusion dictionnary -------------- ###
path_exclusion_dict='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/Exclusion_folder/'
#Open the file and read it
f_dict = open(path_exclusion_dict+'exclusiondic_2010_2014', "rb")
exclusion_dict = pickle.load(f_dict)
f_dict.close()
### --------------- Open and read the exclusion dictionnary -------------- ###

### ------------------- Plot all the exclusion lat/lon ------------------- ###
### --------------------------- Load shapefile --------------------------- ###
#from https://gis.stackexchange.com/questions/113799/how-to-read-a-shapefile-in-python
path_IceBridgeArea_Shape='C:/Users/jullienn/switchdrive/Private/research/backup_Aglaja/working_environment/greenland_topo_data/IceBridge Area Shapefiles/IceBridge Area Shapefiles/IceBridgeArea_Shape.shp'
IceBridgeArea_Shape = gpd.read_file(path_IceBridgeArea_Shape)
### --------------------------- Load shapefile --------------------------- ###

### ----------------------------- Load the DEM --------------------------- ###
# This part is from plot_2002_2003.py'
filename_raster = "C:/Users/jullienn/switchdrive/Private/research/backup_Aglaja/working_environment/greenland_topo_data/elevations/greenland_dem_mosaic_100m_v3.0.tif" #path to raster
grid = Grid.from_raster(filename_raster,data_name='dem')
#Minnor slicing on borders to enhance colorbars
elevDem=grid.dem[:-1,:-1]              
#Scale the colormap
divnorm = mcolors.DivergingNorm(vmin=0, vcenter=1250, vmax=2500)

#Extract elevation from DEM to associated with coordinates. This piece of code
#is from https://gis.stackexchange.com/questions/221292/retrieve-pixel-value-with-geographic-coordinate-as-input-with-gdal
driver = gdal.GetDriverByName('GTiff')

dataset_dem = gdal.Open(filename_raster)
band = dataset_dem.GetRasterBand(1)

cols = dataset_dem.RasterXSize
rows = dataset_dem.RasterYSize

transform_elev = dataset_dem.GetGeoTransform()

xOrigin = transform_elev[0]
yOrigin = transform_elev[3]
pixelWidth = transform_elev[1]
pixelHeight = -transform_elev[5]

data_dem = band.ReadAsArray(0, 0, cols, rows)

#Define ther zero for longitude:
#Where lon==0 is in may09_03_15:
    #lon[886]=23.53372773084396 and lon[887]=-40.08804568537925
    #lat[886]=-3120053.856912824, lat[887]=-3120048.666364133
avg_lon_zero=(23.53372773084396+-40.08804568537925)/2
index_lon_zero=int((avg_lon_zero-xOrigin) / pixelWidth)
### ----------------------------- Load the DEM --------------------------- ###

### --------------------- Display the DEM and shapefile ------------------ ###
# This part is from plot_2002_2003.py'
fig, (ax1) = plt.subplots(1, 1)#, gridspec_kw={'width_ratios': [1, 3]})
fig.suptitle('Exclusion overview')
cb=ax1.imshow(elevDem, extent=grid.extent,cmap=discrete_cmap(10,'cubehelix_r'),norm=divnorm,alpha=0.5)
IceBridgeArea_Shape.plot(ax=ax1)
### --------------------- Display the DEM and shapefile ------------------ ###

### ---------------------------- Plot exclusion -------------------------- ###
transformer = Transformer.from_crs("EPSG:4326", "EPSG:3413", always_xy=True)

#Create an empty dictionnary for appending all the data of a specific year together
yearly_exclusion = {k: {} for k in list(['2010','2011','2012','2013','2014'])}

for indiv_year in yearly_exclusion.keys():
    yearly_exclusion[indiv_year]={k: [] for k in list(['lat_append','lon_append'])}

#Visualize the exclusion
for exclusion_indiv in exclusion_dict.keys():
    
    if (exclusion_dict[exclusion_indiv]=={}):
        #No exclusion, continue
        print(exclusion_indiv,' is empty, continue')
        continue
    else:
        #Exclusion lat/lon to recover
        points_excl=transformer.transform(np.array(exclusion_dict[exclusion_indiv]['lon_exclusion']),np.array(exclusion_dict[exclusion_indiv]['lat_exclusion']))
        
        lon_3413_excl=points_excl[0]
        lat_3413_excl=points_excl[1]
        
        #Recover the year to associate a colour code
        if (exclusion_indiv[0:4] == '2010'):
            color_code='#fef0d9'
        elif (exclusion_indiv[0:4] == '2011'):
            color_code='#fdcc8a'
        elif (exclusion_indiv[0:4] == '2012'):
            color_code='#fc8d59'
        elif (exclusion_indiv[0:4] == '2013'):
            color_code='#e34a33'
        elif (exclusion_indiv[0:4] == '2014'):
            color_code='#b30000'
        
        #Plot lat/lon exclusion
        ax1.scatter(lon_3413_excl,lat_3413_excl,c=color_code,s=0.3)
        
        #Take all the data of that year that have already been appended
        lat_append=yearly_exclusion[exclusion_indiv[0:4]]['lat_append']
        lon_append=yearly_exclusion[exclusion_indiv[0:4]]['lon_append']
        
        #Append all the data together as a function of the year. Take the lat/lon as WGS_84, not EPSG_3413!
        lat_append=np.append(lat_append,exclusion_dict[exclusion_indiv]['lat_exclusion'])
        lon_append=np.append(lon_append,exclusion_dict[exclusion_indiv]['lon_exclusion'])
        
        #Fill in the yearly_exclusion dictionnary
        yearly_exclusion[exclusion_indiv[0:4]]['lat_append']=lat_append
        yearly_exclusion[exclusion_indiv[0:4]]['lon_append']=lon_append       

#Legend preparation
patch2010 = mpatches.Patch(color='#fef0d9', label='2010')
patch2011 = mpatches.Patch(color='#fdcc8a', label='2011')
patch2012 = mpatches.Patch(color='#fc8d59', label='2012')
patch2013 = mpatches.Patch(color='#e34a33', label='2013')
patch2014 = mpatches.Patch(color='#b30000', label='2014')

plt.legend(handles=[patch2010,patch2011,patch2012,patch2013,patch2014])
plt.show()
pdb.set_trace()
### ---------------------------- Plot exclusion -------------------------- ###

### The idea would be to loop through any date in the exclusion_dict to check
# whether there are any lat/lon exclusion available to apply. If yes,
# convert into pixel, if not, leave blank

path_2017_data='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/2017_Greenland_P3/CSARP_qlook'

#folder_day
folder_day='20170422_01'
onlyfiles=list(['Data_20170422_01_168.mat','Data_20170422_01_169.mat','Data_20170422_01_170.mat'])

for indiv_file in onlyfiles:
    #Update path
    pdb.set_trace()
    path_2017_indiv=path_2017_data+'/'+folder_day
    #1. Open 2017 data
    #Open the file and read it
    
    with h5py.File(path_2017_indiv+'/'+indiv_file, 'r') as f:
        f.keys()
        #Select radar echogram
        lat=f['Latitude'][:].transpose() #2017 data should be transposed
        lon=f['Longitude'][:].transpose() #2017 data should be transposed
    
    ### --------------------------- Plot lat/lon ----------------------------- ###
    #Transform the coordinated from WGS84 to EPSG:3413
    #Example from: https://pyproj4.github.io/pyproj/stable/examples.html
    transformer = Transformer.from_crs("EPSG:4326", "EPSG:3413", always_xy=True)
    points=transformer.transform(np.array(lon),np.array(lat))
    
    lon_3413=points[0]
    lat_3413=points[1]
    
    ax1.scatter(lon_3413,lat_3413,c='lime',s=0.1)    
    ax1.set_xlim(-120000,-20000)
    ax1.set_ylim(-2.5e6,-2.4e6)
    plt.show()
    
    #1. Calculate the distance between the point of interest and the exclusion lat/lon
    for i in range(0,lat.size):
        single_lat=lat[0,i]
        single_lon=lon[0,i]
        
        pdb.set_trace()
        
        #Calculate the great circle distance between ith
        #and jth element. This haversine function does not take into
        # account the variation of the earth radius as a function of the
        #latitude. How to use: haversine(lon1, lat1, lon2, lat2)
        dist_pt=[]
        dist_pt=haversine(single_lon, single_lat, df_agg['lon'][j], df_agg['lat'][j])

#Calculation if pixel is done    
#2. If the distance if less than 50m away (to modify?), then take the closest point
#3. Convert the exclusion in pixel
#4. Go through the exclusion pixels and determine whether this is a case x-, -x or x-y
#5. Log the corresponding date with exclusion in pixel.
# RQ: NOTE THAT THE PIXEL LOCATION CONSIDERS TRACEBEGIN_TRACEEND

print('Done plotting exclusion with data of interest')

#1. if exclusion available in less than 50m away, go take the closest exclusion
#available
"""
    ### --------------------------- Plot lat/lon ----------------------------- ###
    
    #2. Loop over exclusion_dict to retreive lat/lon exclusion
    #visualize the exclusion?
    for exclusion_indiv in exclusion_dict.keys():
        if (exclusion_dict[exclusion_indiv]=={}):
            #No exclusion, continue
            print(exclusion_indiv,' is empty, continue')
            continue
        else:
            #Exclusion lat/lon to recover
            #Plot lat/lon exclusion
            
            #exclusion_dict[exclusion_indiv]['lon_exclusion']
            #exclusion_dict[exclusion_indiv]['lat_exclusion']
            
            #Extract lat/lon exclusion
            print('Extracting exclusion')
            
    
    print('Done plotting exclusion lat/lon')
    plt.show()
    pdb.set_trace()

    




#3. Save the corresponding
"""