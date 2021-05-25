# -*- coding: utf-8 -*-
"""
Created on Mon May 24 14:58:51 2021

@author: JullienN
"""

# This code helps to fine the selection of data for 2017 and 2018 dataset
# in order to process only data where iceslabs can be expected.
# To do so, the coordinates of any track is checked whether it belongs to
# the shapefile 'IceBridgeArea_Shape' from MacFerrin et al., 2019. If yes,
# then we keep we associated track number. If not, it is checked it the
# associated track number does indeed not contain ice slabs. If positive,
# then the track can be deleted


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

import geopandas as gpd
import descartes
from osgeo import gdal
import matplotlib.pyplot as plt
from pysheds.grid import Grid
import matplotlib.colors as mcolors
import numpy as np
import os
import h5py
from shapely.geometry import Point, Polygon

from os import listdir
import os.path
from os.path import isfile, join

import pdb

import shapefile as shp  # Requires the pyshp package
import matplotlib.pyplot as plt
from pyproj import Transformer

#1. Load shapefile
### --------------------------- Load shapefile --------------------------- ###
#from https://gis.stackexchange.com/questions/113799/how-to-read-a-shapefile-in-python
path_IceBridgeArea_Shape='C:/Users/jullienn/switchdrive/Private/research/backup_Aglaja/working_environment/greenland_topo_data/IceBridge Area Shapefiles/IceBridge Area Shapefiles/IceBridgeArea_Shape.shp'
IceBridgeArea_Shape = gpd.read_file(path_IceBridgeArea_Shape)
### --------------------------- Load shapefile --------------------------- ###

#2. Load the DEM
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

#3. Plot the DEM and the ice slabs shalefile
### ---------------------------- Display the DEM ------------------------- ###
# This part is from plot_2002_2003.py'
fig, (ax1) = plt.subplots(1, 1)#, gridspec_kw={'width_ratios': [1, 3]})
fig.suptitle('Iceslabs area overview')
cb=ax1.imshow(elevDem, extent=grid.extent,cmap=discrete_cmap(10,'cubehelix_r'),norm=divnorm,alpha=0.5)
### ---------------------------- Display the DEM ------------------------- ###

### ------------------------- Display the shapefile ---------------------- ###

IceBridgeArea_Shape.plot(ax=ax1)
plt.show()
### ------------------------- Display the shapefile ---------------------- ###

#### --------------------------- Plot shapefile --------------------------- ###
##Not usefull anymore
## this part of code is from https://gis.stackexchange.com/questions/131716/plot-shapefile-with-matplotlib
#sf = shp.Reader(path_IceBridgeArea_Shape)
#for shape in sf.shapeRecords():
#    x = [i[0] for i in shape.shape.points[:]]
#    y = [i[1] for i in shape.shape.points[:]]
#    ax1.plot(x,y,'blue','linewidth',0.1)
#plt.show()
### --------------------------- Plot shapefile --------------------------- ###

#4. Load individual 2017 data
path_2017='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/2017_Greenland_P3/CSARP_qlook'
os.chdir(path_2017)

#Prepare the list where the filenames where data is within the iceslabs shapefile
#interest will be saved
data_intersect=[]

# Read the folders of data
folders_track = [ f.name for f in os.scandir(path_2017) if f.is_dir() ]

for folder in folders_track:
    #Update the path to get into the folder of interest
    path_2017_indiv_track=path_2017+'/'+folder
    # Read the files of this specific day
    onlyfiles = [f for f in listdir(path_2017_indiv_track) if isfile(join(path_2017_indiv_track, f))]
    
    #Loop over any file
    for indiv_file in onlyfiles:
        print(indiv_file)
        #Open the file of interest
        with h5py.File(path_2017_indiv_track+'/'+indiv_file, 'r') as f:
            f.keys()
            #Load lat and lon
            lat=f['Latitude'][:].transpose() #2017 data should be transposed
            lon=f['Longitude'][:].transpose() #2017 data should be transposed
        
        ### --------------------------- Plot lat/lon ----------------------------- ###
        #Transform the coordinated from WGS84 to EPSG:3413
        #Example from: https://pyproj4.github.io/pyproj/stable/examples.html
        transformer = Transformer.from_crs("EPSG:4326", "EPSG:3413", always_xy=True)
        points=transformer.transform(np.array(lon),np.array(lat))
        
        lon_3413=points[0]
        lat_3413=points[1]
        
        #ax1.scatter(lon_3413,lat_3413,c='k',s=0.1)
        #ax1.set_xlim(-533400,-532400)
        #ax1.set_ylim(-1335350,-1335450)
        #plt.show()
        ### --------------------------- Plot lat/lon ----------------------------- ###
        
        #Loop over any data point to check whether it belongs to the ice slabs
        #shape area or not
        for i in range(0,lon_3413.size):
            
            single_point=Point(lon_3413[0,i],lat_3413[0,i])
            
            #5. Is data in shapefile? allow a ~100km radius maybe?
            #From: https://automating-gis-processes.github.io/CSC18/lessons/L4/point-in-polygon.html
            check_point=np.asarray(IceBridgeArea_Shape.contains(single_point)).astype(int)
            
            if (np.sum(check_point)>0):
                
                #if yes, the studied point is contained in at least within 1 polygon
                #We can break out from the loop, it is enough, we keep this date as
                #potentially holding iceslabs
                
                #6. Create the list of data to be kept
                data_intersect=np.append(data_intersect,indiv_file[5:20])
                
                print('Trace is within ice slabs area, break')
                break
            
            ##Display the point of interest for checking if this technique works.
            ##Yes, it does work!
            #ax1.scatter(lon_3413[0,i],lat_3413[0,i],c='r',s=0.1)
            #plt.show()
            
            #nb 5 <=> id=6 is TRUE. It works!

print('End of processing')



