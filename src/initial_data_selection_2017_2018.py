# -*- coding: utf-8 -*-
"""
Created on Thu Nov 25 18:53:55 2021

@author: jullienn

This code is inspired from refine_location_2017_2018.py
"""


# This code helps to fine the selection of data for 2017 and 2018 dataset
# in order to process only data where iceslabs can be expected.
# To do so, the coordinates of any track is checked whether it belongs to
# the shapefile 'IceBridgeArea_Shape' from MacFerrin et al., 2019. If yes,
# then we keep we associated track number. Once this process is done, I should
# carefully check the first and last tracks of any continuous trace contained
# within this shape to make sure I do not miss ice slabs below and above this
# area. Then I can delete the remaining tracks.

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
from pysheds.grid import Grid
import matplotlib.colors as mcolors
from osgeo import gdal

import descartes
import matplotlib.pyplot as plt

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
import pandas as pd

year_to_refine='2017'

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

#3. Load the 2017-2018 coordinates
### -------------- Load the 2017 and 2018 coordinates -------------------- ###
path_coordinates='C:/Users/jullienn/Documents/working_environment/Extended_Greenland_iceslabs/data/2017_2018_coordinates/'
coordinates_2017=pd.read_csv(path_coordinates+'2017_Greenland_P3.csv',dtype={'FRAME':str})
coordinates_2018=pd.read_csv(path_coordinates+'2018_Greenland_P3.csv',dtype={'FRAME':str})
### -------------- Load the 2017 and 2018 coordinates -------------------- ###

#Add a column of str of FRAMES for data selection
coordinates_2017['dates']=[x[:10] for x in coordinates_2017['FRAME']]
coordinates_2018['dates']=[x[:10] for x in coordinates_2018['FRAME']]
coordinates_2017['files']=[x for x in coordinates_2017['FRAME']]
coordinates_2018['files']=[x for x in coordinates_2018['FRAME']]

#Identify unique dates in 2017 and 2018 data
list_2017_indiv=list(np.unique([x[:10] for x in coordinates_2017['FRAME']]))
list_2018_indiv=list(np.unique([x[:10] for x in coordinates_2018['FRAME']]))

#Append the list
list_2017_2018=list_2017_indiv+list_2018_indiv

#Define path where figures need to be saved
path_figname='C:/Users/jullienn/Documents/working_environment/Extended_Greenland_iceslabs/data/2017_2018_coordinates/figures_selection_datatobeprocessed/'

#Prepare the list where the filenames where data is within the iceslabs shapefile
#interest will be saved
data_intersect=[]

count=0
#5. Loop through individual day of data collection
for indiv_date in list_2017_2018:
    print(count/len(list_2017_2018)*100,' %')
    print('   Perform identification')
    
    #a. Retrieve the coordinates of that date, as vector
    if (indiv_date[0:4]=='2017'):
        df_indiv=coordinates_2017[coordinates_2017['dates']==indiv_date]
    elif (indiv_date[0:4]=='2018'):
        df_indiv=coordinates_2018[coordinates_2018['dates']==indiv_date]
    else:
        print('Year not known')
        break
    
    ### --------------------------- Plot lat/lon ----------------------------- ###
    #Transform the coordinated from WGS84 to EPSG:3413
    #Example from: https://pyproj4.github.io/pyproj/stable/examples.html
    transformer = Transformer.from_crs("EPSG:4326", "EPSG:3413", always_xy=True)
    points=transformer.transform(np.array(df_indiv.LON),np.array(np.array(df_indiv.LAT)))
    
    lon_3413=points[0]
    lat_3413=points[1]
    
    #b. Display coordinates above GrIS DEM and ice slabs shapefile
    fig, (ax1) = plt.subplots(1, 1)#, gridspec_kw={'width_ratios': [1, 3]})
    fig.suptitle(indiv_date)
    cb=ax1.imshow(elevDem, extent=grid.extent,cmap=discrete_cmap(10,'cubehelix_r'),norm=divnorm,alpha=0.5)
    IceBridgeArea_Shape.plot(ax=ax1)
    
    ax1.scatter(lon_3413,lat_3413,c='gray',s=0.1)
    #Display the start of the trace
    ax1.scatter(lon_3413[0],lat_3413[0],c='red',s=10)
    #Zoom over the trace
    ax1.set_xlim(np.min(lon_3413)-1000,np.max(lon_3413)+1000)
    ax1.set_ylim(np.min(lat_3413)-1000,np.max(lat_3413)+1000)
    ### --------------------------- Plot lat/lon ----------------------------- ###
    
    #c. Loop through the indiv files
    for indiv_file in np.unique(df_indiv['files']):
        # Retrieve the coordinates of that indiv file, as vector
        df_indiv_file=df_indiv[df_indiv['files']==indiv_file]
        
        #Transform the coordinates from WGS84 to EPSG:3413
        #Example from: https://pyproj4.github.io/pyproj/stable/examples.html
        transformer = Transformer.from_crs("EPSG:4326", "EPSG:3413", always_xy=True)
        points_indiv=transformer.transform(np.array(df_indiv_file.LON),np.array(np.array(df_indiv_file.LAT)))
        
        lon_3413_indiv=points_indiv[0]
        lat_3413_indiv=points_indiv[1]
        
        #d. Do the intersection with ice slabs shapefile
        #Loop over any data point to check whether it belongs to the ice slabs
        #shape area or not
        for i in range(0,lon_3413_indiv.size):
            single_point=Point(lon_3413_indiv[i],lat_3413_indiv[i])
            #From: https://automating-gis-processes.github.io/CSC18/lessons/L4/point-in-polygon.html
            check_point=np.asarray(IceBridgeArea_Shape.contains(single_point)).astype(int)

            if (np.sum(check_point)>0):
                
                #if yes, the studied point is contained in at least within 1 polygon
                #We can break out from the loop, it is enough, we keep this date as
                #potentially holding iceslabs
                
                #6. Create the list of data to be kept
                data_intersect=np.append(data_intersect,indiv_file[0:10]+'_'+indiv_file[10:13])
                ax1.scatter(lon_3413_indiv,lat_3413_indiv,c='lime',s=0.1)
                break
        
        #Display the start and end of trace, and number of the trace
        ax1.scatter(lon_3413_indiv[0],lat_3413_indiv[0],marker="|",s=0.75,c='black')
        ax1.scatter(lon_3413_indiv[-1],lat_3413_indiv[-1],marker="|",s=0.75,c='black')
        ax1.text(np.median(lon_3413_indiv),np.median(lat_3413_indiv),int(indiv_file[10:13]),fontsize=5)
    
    print('   Saving figure')
    #Define the figname
    figName=path_figname+indiv_date+'_data_selection.png'
    # Get the current figure like in MATLAB
    fig = plt.gcf()
    plt.show() # show it here (important, if done before you will get blank picture)
    fig.set_size_inches((8.9, 5), forward=False)
    #fig.savefig(figName, dpi=500) # Change is over here
    plt.close()
    count=count+1
'''    
#Save the data
with open(path_figname+"intial_data_selection_20172018.txt", 'w') as output:
    for row in data_intersect:
        output.write(str(row) + '\n')
'''        
print('End of processing')



