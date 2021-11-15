# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 16:07:38 2021

@author: jullienn
"""

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

def calcul_elevation(lon,lat,data_dem,yOrigin,pixelHeight,pixelWidth,index_lon_zero):
    
    if (np.isnan(lon) or np.isnan(lat)):
        #elev_all=np.append(elev_all,np.nan)
        elevation=np.nan
    else:
        #The origin is top left corner!!
        #y will always be negative
        row = int((yOrigin - lat ) / pixelHeight)
        if (lon<0):
            # if x negative
            col = index_lon_zero-int((-lon-0) / pixelWidth)
        elif (lon>0):
            # if x positive
            col = index_lon_zero+int((lon-0) / pixelWidth)
        #Read the elevation
        elevation=data_dem[row][col]
    
    return elevation


#Import packages
import rasterio
from rasterio.plot import show
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pandas as pd
from os import listdir
from os.path import isfile, join
import pickle
from pysheds.grid import Grid
import pdb
import numpy as np
from pyproj import Transformer
import matplotlib.gridspec as gridspec
import scipy.io
from osgeo import gdal
import geopandas as gpd  # Requires the pyshp package

from matplotlib.colors import ListedColormap, BoundaryNorm

create_elevation_dictionaries='FALSE'
#pdb.set_trace()

########################## Load GrIS elevation ##########################
#Open the DEM
grid = Grid.from_raster("C:/Users/jullienn/switchdrive/Private/research/backup_Aglaja/working_environment/greenland_topo_data/elevations/greenland_dem_mosaic_100m_v3.0.tif",data_name='dem')
#Minnor slicing on borders to enhance colorbars
elevDem=grid.dem[:-1,:-1]              
#Scale the colormap
divnorm = mcolors.DivergingNorm(vmin=0, vcenter=1250, vmax=2500)
########################## Load GrIS elevation ##########################


############################ Load DEM information ############################
#Extract elevation from DEM to associated with coordinates. This piece of code
#is from https://gis.stackexchange.com/questions/221292/retrieve-pixel-value-with-geographic-coordinate-as-input-with-gdal
driver = gdal.GetDriverByName('GTiff')
filename_raster = "C:/Users/jullienn/switchdrive/Private/research/backup_Aglaja/working_environment/greenland_topo_data/elevations/greenland_dem_mosaic_100m_v3.0.tif" #path to raster

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
############################ Load DEM information ############################

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
### -------------------------- Load shapefiles --------------------------- ###
#pdb.set_trace()

if (create_elevation_dictionaries == 'TRUE'):
    
    ################# Load 2002-2003 flightlines coordinates ################
    path_data='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/icelens_identification'
    
    #Open the file and read it
    f_flightlines = open(path_data+'/metadata_coord_2002_2003', "rb")
    all_2002_3_flightlines = pickle.load(f_flightlines)
    f_flightlines.close()
    ################# Load 2002-2003 flightlines coordinates ################
    
    lat_all=[]
    lon_all=[]
    #elev_all=[]
    
    elevation_dictionnary = {k: {} for k in list(['2002','2003'])}
    
    for year in list(all_2002_3_flightlines.keys()):
        
        elevation_dictionnary[year]={k: {} for k in list(all_2002_3_flightlines[year].keys())}
        
        for days in list(all_2002_3_flightlines[year].keys()):
            
            elevation_dictionnary[year][days]={k: {} for k in list(all_2002_3_flightlines[year][days].keys())}
            
            for indiv_file in list(all_2002_3_flightlines[year][days].keys()):
                if (indiv_file[0:7]=='quality'):
                    continue
                else:
                    print(indiv_file)
                    lat_all=np.append(lat_all,all_2002_3_flightlines[year][days][indiv_file][0])
                    lon_all=np.append(lon_all,all_2002_3_flightlines[year][days][indiv_file][1])
                    #Extract the elevation:
                    lat_elev=[]
                    lon_elev=[]
                    if (days=='jun04'):
                        lat_elev=np.transpose(all_2002_3_flightlines[year][days][indiv_file][0])
                        lon_elev=np.transpose(all_2002_3_flightlines[year][days][indiv_file][1])
                    else:
                        lat_elev=all_2002_3_flightlines[year][days][indiv_file][0]
                        lon_elev=all_2002_3_flightlines[year][days][indiv_file][1]
                    
                    latlon_tuple=[]
                    latlon_tuple=list(zip(lon_elev,lat_elev))
                    
                    elev_indiv_file=[]
                    for indiv_coord in latlon_tuple:
                        if (np.isnan(indiv_coord[0]) or np.isnan(indiv_coord[1])):
                            #elev_all=np.append(elev_all,np.nan)
                            elev_indiv_file=np.append(elev_indiv_file,np.nan)
                        else:
                            #The origin is top left corner!!
                            #y will always be negative
                            row = int((yOrigin - indiv_coord[1] ) / pixelHeight)
                            if (indiv_coord[0]<0):
                                # if x negative
                                col = index_lon_zero-int((-indiv_coord[0]-0) / pixelWidth)
                            elif (indiv_coord[0]>0):
                                # if x positive
                                col = index_lon_zero+int((indiv_coord[0]-0) / pixelWidth)
                            #Read the elevation
                            #elev_all=np.append(elev_all,data_dem[row][col])
                            elev_indiv_file=np.append(elev_indiv_file,data_dem[row][col])
                    
                    #Store data into the dictionnary
                    elevation_dictionnary[year][days][indiv_file]=elev_indiv_file
                    
    ################# Load 2002-2003 flightlines coordinates ################
    
    ################### Load 2002-2003 ice lenses location ##################
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
    filename_20102018='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2010_2018/excel_spatial_aggreation_and_other/final_excel/Ice_Layer_Output_Thicknesses_2010_2018_jullienetal2021.csv'
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
    
    #Prepare plot
    fig, (ax1) = plt.subplots(1, 1)#, gridspec_kw={'width_ratios': [1, 3]})
    fig.suptitle('Iceslabs area overview')
    
    #Display DEM
    cb1=ax1.imshow(elevDem, extent=grid.extent,cmap=discrete_cmap(10,'cubehelix_r'),alpha=0.5,norm=divnorm)
    cbar1=fig.colorbar(cb1, ax=[ax1], location='left')
    cbar1.set_label('Elevation [m]')
    
    #Display the shapefile
    IceBridgeArea_Shape.plot(ax=ax1)
    
    #Plot all the 2002-2003 icelenses according to their condifence color
    #1. Red
    ax1.scatter(df_2002_2003[df_2002_2003['colorcode_icelens']==-1]['lon_3413'], df_2002_2003[df_2002_2003['colorcode_icelens']==-1]['lat_3413'],s=1,facecolors='#c9662c', edgecolors='none')
    #2. Orange
    ax1.scatter(df_2002_2003[df_2002_2003['colorcode_icelens']==0]['lon_3413'], df_2002_2003[df_2002_2003['colorcode_icelens']==0]['lat_3413'],s=1,facecolors='#fed976', edgecolors='none')
    #3. Green
    ax1.scatter(df_2002_2003[df_2002_2003['colorcode_icelens']==1]['lon_3413'], df_2002_2003[df_2002_2003['colorcode_icelens']==1]['lat_3413'],s=1,facecolors='#238b45', edgecolors='none')
    ##Purple
    #ax1.scatter(lon_icelens[colorcode_icelens==2], lat_icelens[colorcode_icelens==2],s=1,facecolors='purple', edgecolors='none')
    
    #Correct zoom
    ax1.set_xlim(-650000,900000)
    ax1.set_ylim(-3360000,-650000)
    
    plt.show()
    pdb.set_trace()
    # compare min and max of lat/lon of the track with respect to shapefile
    from shapely.geometry import Point, Polygon
    
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
        check_NW_icecap_greenland=np.asarray(NW_icecap_greenland_mask.contains(single_point)).astype(int)
        check_NW_north_greenland=np.asarray(NW_north_greenland_mask.contains(single_point)).astype(int)
        check_NW_west_greenland=np.asarray(NW_west_greenland_mask.contains(single_point)).astype(int)
        check_SW_lower_greenland=np.asarray(SW_lower_greenland_mask.contains(single_point)).astype(int)
        check_SW_middle_greenland=np.asarray(SW_middle_greenland_mask.contains(single_point)).astype(int)
        check_SW_upper_greenland=np.asarray(SW_upper_greenland_mask.contains(single_point)).astype(int)
    
        #Associated the point of interest to its regional shapefile in data_iceslabs
        if (np.sum(check_NW_icecap_greenland)>0):
            df_20102018['key_shp'][i]='NW_icecap'
        elif (np.sum(check_NW_north_greenland)>0):
            df_20102018['key_shp'][i]='NW_north'
        elif (np.sum(check_NW_west_greenland)>0):
            df_20102018['key_shp'][i]='NW_west'
        elif (np.sum(check_SW_lower_greenland)>0):
            df_20102018['key_shp'][i]='SW_lower'
        elif (np.sum(check_SW_middle_greenland)>0):
            df_20102018['key_shp'][i]='SW_middle'
        elif (np.sum(check_SW_upper_greenland)>0):
            df_20102018['key_shp'][i]='SW_upper'
        else:
            df_20102018['key_shp'][i]='Out'
        
        #Calculate the corresponding elevation
        df_20102018['elevation'][i]=calcul_elevation(df_20102018['lon_3413'][i],df_20102018['lat_3413'][i],data_dem,yOrigin,pixelHeight,pixelWidth,index_lon_zero)
        #Add the year
        df_20102018['year'][i]=int(df_20102018['Track_name'][i][0:4])
    
        #Monitor the process
        print(i/lon_3413_20102018.size*100,'%')
    
    
    #Save the dictionary into a picke file
    filename_tosave='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2010_2018/excel_spatial_aggreation_and_other/df_20102018_with_elevation'
    outfile= open(filename_tosave, "wb" )
    pickle.dump(df_20102018,outfile)
    outfile.close()
    
    
    #Display the keys
    fig, (ax1) = plt.subplots(1, 1)#, gridspec_kw={'width_ratios': [1, 3]})
    fig.suptitle('Iceslabs keys')
    
    #Display DEM
    cb1=ax1.imshow(elevDem, extent=grid.extent,cmap=discrete_cmap(10,'cubehelix_r'),alpha=0.5,norm=divnorm)
    cbar1=fig.colorbar(cb1, ax=[ax1], location='left')
    cbar1.set_label('Elevation [m]')
    
    #Display the shapefile
    NW_icecap_greenland_mask.plot(ax=ax1)
    NW_north_greenland_mask.plot(ax=ax1)
    NW_west_greenland_mask.plot(ax=ax1)
    SW_lower_greenland_mask.plot(ax=ax1)
    SW_middle_greenland_mask.plot(ax=ax1)
    SW_upper_greenland_mask.plot(ax=ax1)
    
    #Display the data as a function of their belonging keys
    ax1.scatter(df_20102018['lon_3413'][df_20102018['key_shp']=='NW_icecap'],df_20102018['lat_3413'][df_20102018['key_shp']=='NW_icecap'],facecolors='orange')
    ax1.scatter(df_20102018['lon_3413'][df_20102018['key_shp']=='NW_west'],df_20102018['lat_3413'][df_20102018['key_shp']=='NW_west'],facecolors='blue')
    ax1.scatter(df_20102018['lon_3413'][df_20102018['key_shp']=='NW_north'],df_20102018['lat_3413'][df_20102018['key_shp']=='NW_north'],facecolors='purple')
    ax1.scatter(df_20102018['lon_3413'][df_20102018['key_shp']=='SW_lower'],df_20102018['lat_3413'][df_20102018['key_shp']=='SW_lower'],facecolors='red')
    ax1.scatter(df_20102018['lon_3413'][df_20102018['key_shp']=='SW_middle'],df_20102018['lat_3413'][df_20102018['key_shp']=='SW_middle'],facecolors='green')
    ax1.scatter(df_20102018['lon_3413'][df_20102018['key_shp']=='SW_upper'],df_20102018['lat_3413'][df_20102018['key_shp']=='SW_upper'],facecolors='k')
    
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
        check_NW_icecap_greenland=np.asarray(NW_icecap_greenland_mask.contains(single_point)).astype(int)
        check_NW_north_greenland=np.asarray(NW_north_greenland_mask.contains(single_point)).astype(int)
        check_NW_west_greenland=np.asarray(NW_west_greenland_mask.contains(single_point)).astype(int)
        check_SW_lower_greenland=np.asarray(SW_lower_greenland_mask.contains(single_point)).astype(int)
        check_SW_middle_greenland=np.asarray(SW_middle_greenland_mask.contains(single_point)).astype(int)
        check_SW_upper_greenland=np.asarray(SW_upper_greenland_mask.contains(single_point)).astype(int)
    
        #Associated the point of interest to its regional shapefile in data_iceslabs
        if (np.sum(check_NW_icecap_greenland)>0):
            df_2002_2003['key_shp'].iloc[i]='NW_icecap'
        elif (np.sum(check_NW_north_greenland)>0):
            df_2002_2003['key_shp'].iloc[i]='NW_north'
        elif (np.sum(check_NW_west_greenland)>0):
            df_2002_2003['key_shp'].iloc[i]='NW_west'
        elif (np.sum(check_SW_lower_greenland)>0):
            df_2002_2003['key_shp'].iloc[i]='SW_lower'
        elif (np.sum(check_SW_middle_greenland)>0):
            df_2002_2003['key_shp'].iloc[i]='SW_middle'
        elif (np.sum(check_SW_upper_greenland)>0):
            df_2002_2003['key_shp'].iloc[i]='SW_upper'
        else:
            df_2002_2003['key_shp'].iloc[i]='Out'
        
        #Calculate the corresponding elevation
        df_2002_2003['elevation'].iloc[i]=calcul_elevation(df_2002_2003['lon_3413'].iloc[i],df_2002_2003['lat_3413'].iloc[i],data_dem,yOrigin,pixelHeight,pixelWidth,index_lon_zero)
        
        #Add the year
        if (df_2002_2003['Track_name'][i][6:8] == '02'):
            year_to_write=2002
        elif (df_2002_2003['Track_name'][i][6:8] == '03'):
            year_to_write=2003
        else:
            print('Year not known, error')
            break
        
        df_2002_2003['year'][i]=year_to_write
        
        #Monitor the process
        print(i/len(df_2002_2003)*100,'%')
    
    #Only work with green slabs
    df_2002_2003_green=df_2002_2003[df_2002_2003['colorcode_icelens']==1]
    
    #Display the data as a function of their belonging keys
    ax1.scatter(df_2002_2003_green['lon_3413'][df_2002_2003_green['key_shp']=='NW_icecap'],df_2002_2003_green['lat_3413'][df_2002_2003_green['key_shp']=='NW_icecap'],facecolors='brown')
    ax1.scatter(df_2002_2003_green['lon_3413'][df_2002_2003_green['key_shp']=='NW_west'],df_2002_2003_green['lat_3413'][df_2002_2003_green['key_shp']=='NW_west'],facecolors='cyan')
    ax1.scatter(df_2002_2003_green['lon_3413'][df_2002_2003_green['key_shp']=='NW_north'],df_2002_2003_green['lat_3413'][df_2002_2003_green['key_shp']=='NW_north'],facecolors='pink')
    ax1.scatter(df_2002_2003_green['lon_3413'][df_2002_2003_green['key_shp']=='SW_lower'],df_2002_2003_green['lat_3413'][df_2002_2003_green['key_shp']=='SW_lower'],facecolors='yellow')
    ax1.scatter(df_2002_2003_green['lon_3413'][df_2002_2003_green['key_shp']=='SW_middle'],df_2002_2003_green['lat_3413'][df_2002_2003_green['key_shp']=='SW_middle'],facecolors='olive')
    ax1.scatter(df_2002_2003_green['lon_3413'][df_2002_2003_green['key_shp']=='SW_upper'],df_2002_2003_green['lat_3413'][df_2002_2003_green['key_shp']=='SW_upper'],facecolors='gray')
    
    plt.show()
    
    #Save the dictionary into a picke file
    filename_tosave='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2010_2018/excel_spatial_aggreation_and_other/df_2002_2003_with_elevation'
    outfile= open(filename_tosave, "wb" )
    pickle.dump(df_2002_2003,outfile)
    outfile.close()
    
    pdb.set_trace()
else:
    #Dictionnaries have already been created, load them
    path_df_with_elevation='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2010_2018/excel_spatial_aggreation_and_other/' 
    
    #Load 2002-2003
    f_20022003 = open(path_df_with_elevation+'df_2002_2003_with_elevation', "rb")
    df_2002_2003 = pickle.load(f_20022003)
    f_20022003.close()
    
    #Load 2010-2018
    f_20102018 = open(path_df_with_elevation+'df_20102018_with_elevation', "rb")
    df_2010_2018 = pickle.load(f_20102018)
    f_20102018.close()

#IV. From here on, work with the different periods separated by strong melting summers.
#    Work thus with 2002-2003 VS 2010 VS 2011-2012 VS 2013-2014 VS 2017-2018
#    Select the absolute low and absolute high of 2002-2003, 2010-2014 and 2017-2018

#Let's create ~10km latitudinal (resp. longitudinal) slices for SW Greenland (resp. NW Greenland)
#and calculate the low and high end in each slice for elevation difference:

#1. Create the latitudinal (resp. longitudinal) slices
############ Is this grid to change?? This is about 10km width but not even whether north or south!
lat_slices=np.linspace(-2800000,-1490000,int((np.abs(-2800000)-np.abs(-1490000))/10000))
lon_slices=np.linspace(-600000,650000,int((np.abs(650000)+np.abs(-600000))/10000))

#2. Select and store all the data belonging to the lon/lat slices in a dictionnary.
#### ------------------------- 2010-2018 -------------------------------- ####
#   Retreive and store min and max elevation of each slice in a dataframe
#   ----- Latitudinal slices

#Create a dictionnary where to store slices information
dict_lat_slice={}

#Create a dictionnary to store np arrays storing slices min and max elevation for each region
dict_lat_slices_summary={k: {} for k in list(df_2010_2018['key_shp'].unique())}

#loop over the regions, create the room for each time period in each region
for region in list(df_2010_2018['key_shp'].unique()):
    dict_lat_slices_summary[region]={k: {} for k in list(['2010','2011-2012','2013-2014','2017-2018'])}
    
    for time_period in dict_lat_slices_summary[region].keys():
        #Fill the dict_lat_slices_summary dictionnary with zeros
        dict_lat_slices_summary[region][time_period]=np.zeros((len(lat_slices),2))*np.nan

#Loop over each boundary of lat slices and store dataset related to slices
for i in range(1,len(lat_slices)):
    
    #Identify low and higher end of the slice
    low_bound=lat_slices[i-1]
    high_bound=lat_slices[i]
    
    #Select all the data belonging to this slice
    ind_slice=np.logical_and(np.array(df_2010_2018['lat_3413']>=low_bound),np.array(df_2010_2018['lat_3413']<high_bound))
    df_slice=df_2010_2018[ind_slice]
    
    #Store the associated df
    dict_lat_slice[str(int(lat_slices[i-1]))+' to '+str(int(lat_slices[i]))]=df_slice   
    
    #Loop over the regions present in df_slice
    for region in list(df_slice['key_shp'].unique()):
        #Select only the data belonging to this region
        df_region=df_slice[df_slice['key_shp']==region]
        
        #Loop over the different time periods (2010, 2011-2012, 2013-2014, 2017-2018)
        for time_period in list(['2010','2011-2012','2013-2014','2017-2018']):
            if (time_period == '2010'):
                df_region_period=df_region[df_region['year']==2010]
            elif (time_period == '2011-2012'):
                df_region_period=df_region[(df_region['year']>=2011) & (df_region['year']<=2012)]
            elif (time_period == '2013-2014'):
                df_region_period=df_region[(df_region['year']>=2013) & (df_region['year']<=2014)]
            elif (time_period == '2017-2018'):
                df_region_period=df_region[(df_region['year']>=2017) & (df_region['year']<=2018)]
            else:
                print('Time period not known, break')
                break
            #Identify min and max of each region and store them into a dataframe
            #Retreive the stored array
            array_region_indiv=dict_lat_slices_summary[region][time_period]
            #Store min and max of this regional slice
            array_region_indiv[i,0]=np.min(df_region_period['elevation'])
            array_region_indiv[i,1]=np.max(df_region_period['elevation'])
            #Store again data into dict_lat_slices_summary
            dict_lat_slices_summary[region][time_period]=array_region_indiv
            
#   ----- Longitudinal slices
#Create a dictionnary where to store slices information
dict_lon_slice={}

#Create a dictionnary to store np arrays storing slices min and max elevation for each region
dict_lon_slices_summary={k: {} for k in list(df_2010_2018['key_shp'].unique())}

#loop over the regions, create the room for each time period in each region
for region in list(df_2010_2018['key_shp'].unique()):
    dict_lon_slices_summary[region]={k: {} for k in list(['2010','2011-2012','2013-2014','2017-2018'])}
    
    for time_period in dict_lon_slices_summary[region].keys():
        #Fill the dict_lon_slices_summary dictionnary with zeros
        dict_lon_slices_summary[region][time_period]=np.zeros((len(lon_slices),2))*np.nan

#Loop over each boundary of lon slices and store dataset related to slices
for i in range(1,len(lon_slices)):
    
    #Identify low and higher end of the slice
    low_bound=lon_slices[i-1]
    high_bound=lon_slices[i]
    
    #Select all the data belonging to this slice
    ind_slice=np.logical_and(np.array(df_2010_2018['lon_3413']>=low_bound),np.array(df_2010_2018['lon_3413']<high_bound))
    df_slice=df_2010_2018[ind_slice]
    
    #Store the associated df
    dict_lon_slice[str(int(lon_slices[i-1]))+' to '+str(int(lon_slices[i]))]=df_slice   
    
    #Loop over the regions present in df_slice
    for region in list(df_slice['key_shp'].unique()):
        #Select only the data belonging to this region
        df_region=df_slice[df_slice['key_shp']==region]
        
        #Loop over the different time periods (2010, 2011-2012, 2013-2014, 2017-2018)
        for time_period in list(['2010','2011-2012','2013-2014','2017-2018']):
            if (time_period == '2010'):
                df_region_period=df_region[df_region['year']==2010]
            elif (time_period == '2011-2012'):
                df_region_period=df_region[(df_region['year']>=2011) & (df_region['year']<=2012)]
            elif (time_period == '2013-2014'):
                df_region_period=df_region[(df_region['year']>=2013) & (df_region['year']<=2014)]
            elif (time_period == '2017-2018'):
                df_region_period=df_region[(df_region['year']>=2017) & (df_region['year']<=2018)]
            else:
                print('Time period not known, break')
                break
            #Identify min and max of each region and store them into a dataframe
            #Retreive the stored array
            array_region_indiv=dict_lon_slices_summary[region][time_period]
            #Store min and max of this regional slice
            array_region_indiv[i,0]=np.min(df_region_period['elevation'])
            array_region_indiv[i,1]=np.max(df_region_period['elevation'])
            #Store again data into dict_lat_slices_summary
            dict_lon_slices_summary[region][time_period]=array_region_indiv

#### ------------------------- 2010-2018 -------------------------------- ####

pdb.set_trace()
#3. Associate each slice to its belonging region.
#   Not needed! Already present in dataframes!

#4. Calculate the average minimum and maximum of each region among the slices

#5. Flag the more or less perpendicularly crossing 2002-2003 flight lines and exclude the one not crossing
flag_low=['jun04_02proc_4.mat','jun04_02proc_36.mat','jun04_02proc_52.mat','jun04_02proc_53.mat',
      'may09_03_0_aggregated','may09_03_1_aggregated','may09_03_30_aggregated',
      'may09_03_37_aggregated','may11_03_8_aggregated','may11_03_12_aggregated',
      'may11_03_13_aggregated','may11_03_16_aggregated','may11_03_20_aggregated',
      'may11_03_21_aggregated','may11_03_38_aggregated','may11_03_39_aggregated',
      'may12_03_1_aggregated','may12_03_2_aggregated','may12_03_11_aggregated',
      'may12_03_15_aggregated','may12_03_36_aggregated','may13_03_30_aggregated',
      'may14_03_1_aggregated','may14_03_2_aggregated','may14_03_7_aggregated',
      'may14_03_8_aggregated','may14_03_20_aggregated','may14_03_21_aggregated',
      'may15_03_0_aggregated','may15_03_2_aggregated','may15_03_4_aggregated',
      'may15_03_9_aggregated','may18_02_27_aggregated']

flag_high=['jun04_02proc_4.mat','jun04_02proc_36.mat','jun04_02proc_52.mat','jun04_02proc_53.mat',
      'may09_03_0_aggregated','may09_03_1_aggregated','may09_03_30_aggregated',
      'may09_03_37_aggregated','may11_03_20_aggregated','may11_03_21_aggregated',
      'may11_03_37_aggregated','may11_03_38_aggregated','may12_03_1_aggregated',
      'may12_03_2_aggregated','may12_03_11_aggregated','may12_03_36_aggregated',
      'may13_03_30_aggregated','may14_03_1_aggregated','may14_03_2_aggregated',
      'may14_03_7_aggregated','may14_03_20_aggregated','may14_03_21_aggregated',
      'may15_03_2_aggregated','may15_03_4_aggregated','may15_03_9_aggregated',
      'may18_02_27_aggregated']

unique_flags=np.unique(np.append(flag_low,flag_high))

#6. Take the absolute min and max of all 2002-2003 ice slabs in a specific region
#A suite of 2002-2003 traces do not belong to different region, which ease coding

#Here are the traces. For consecutive ones, the ice slabs range elevation is distributed through consecutive traces
traces=[['jun04_02proc_4.mat'],
        ['jun04_02proc_36.mat'],
        ['jun04_02proc_52.mat','jun04_02proc_53.mat'],
        ['may09_03_0_aggregated','may09_03_1_aggregated'],
        ['may09_03_30_aggregated'],
        ['may09_03_37_aggregated'],
        ['may11_03_8_aggregated'],
        ['may11_03_12_aggregated','may11_03_13_aggregated'],
        ['may11_03_16_aggregated'],
        ['may11_03_20_aggregated','may11_03_21_aggregated'],
        ['may11_03_37_aggregated','may11_03_38_aggregated','may11_03_39_aggregated'],
        ['may12_03_1_aggregated','may12_03_2_aggregated'],
        ['may12_03_11_aggregated'],
        ['may12_03_15_aggregated'],
        ['may12_03_36_aggregated'],
        ['may13_03_30_aggregated'],
        ['may14_03_1_aggregated','may14_03_2_aggregated'],
        ['may14_03_7_aggregated','may14_03_8_aggregated'],
        ['may14_03_20_aggregated','may14_03_21_aggregated'],
        ['may15_03_0_aggregated'],
        ['may15_03_2_aggregated'],
        ['may15_03_4_aggregated'],
        ['may15_03_9_aggregated'],
        ['may18_02_27_aggregated']]

list_traces=[item for sublist in traces for item in sublist]

#Loop over the traces, check the flags and populate a low end and high end vector where applicable.
#If consecutive traces, consider the suite of traces!

#Create the dictionary
dict_summary_2002_2003={k: {} for k in list(df_2002_2003['key_shp'].unique())}

#Fill the dict_summary_2002_2003 dictionnary with a np.nan
for region in list(df_2002_2003['key_shp'].unique()):
    dict_summary_2002_2003[region]=np.zeros((len(traces),2))*np.nan

count=0
#Loop over the traces
for trace in traces:
    
    #Check whether we are dealing with single or consecutive traces
    if(len(trace)>1):
        #We are dealing with consecutive traces
        #Select the data related to the first trace
        data_trace=df_2002_2003[df_2002_2003['Track_name']==trace[0]]
        
        #loop over the traces and append data to each other, do not take the first one
        for indiv_trace in list(trace[1:]):
            #Select all the data related to this trace
            data_trace=data_trace.append(df_2002_2003[df_2002_2003['Track_name']==indiv_trace])
            
    else:
        #We are dealing with individual traces
        #Select all the data related to this trace
        data_trace=df_2002_2003[df_2002_2003['Track_name']==trace[0]]

    #Now my data_trace datasets are ready to be worked with
    #Keep only green ice slabs
    data_trace=data_trace[data_trace['colorcode_icelens']==1]
    
    if (len(data_trace)<1):
        #No green ice slabs, continue
        continue
    else:
        #Identify the region
        region=list(np.unique(data_trace['key_shp']))
        
        #Retreive the stored array
        array_region_indiv=dict_summary_2002_2003[region[0]]
        
        #Check the flag: shall we store data?
        if trace[0] in list(flag_low):
            #Store min in array_region_indiv
            array_region_indiv[count,0]=np.min(data_trace['elevation'])
        
        if trace[0] in list(flag_high):
            #Store max in array_region_indiv
            array_region_indiv[count,1]=np.max(data_trace['elevation'])
        
        #Update count
        count=count+1
        #Store again data into dict_lat_slices_summary
        dict_summary_2002_2003[region[0]]=array_region_indiv

#7. Do the elevation difference and eventually the corresponding distance calculation in each region
pdb.set_trace()
#Create a dictionnary where to store relevant information
dict_summary={k: {} for k in list(df_2010_2018['key_shp'].unique())}

#Loop over the regions
for region in list(df_2010_2018['key_shp'].unique()):
    
    #Continue building the dictionnary
    dict_summary[region]={k: {} for k in list(['2002-2003','2010','2011-2012','2013-2014','2017-2018'])}
    
    #Loop over the 5 time periods
    
    for time_period in list(['2002-2003','2010','2011-2012','2013-2014','2017-2018']):
        dict_summary[region][time_period]={k: {} for k in list(['min_elev','max_elev'])}
        
        #Take the average of low and high elevation where ice slabs have been
        #identified in this region, no matter the year in this specific time
        #period, and store relevant information
        
        if (time_period=='2002-2003'):
            #Retreive the corresponding matrix where data are stored
            dict_temp=dict_summary_2002_2003[region]
        else:
            #The dictionnary to select is different whether we are in north or south greenland
            if (region in list(['NW_north','NW_west','NW_icecap'])):
                dict_temp=dict_lon_slices_summary[region][time_period]
            else:
                dict_temp=dict_lat_slices_summary[region][time_period]
        
        #Calculate and store averages
        dict_summary[region][time_period]['min_elev']=np.nanmean(dict_temp[:,0])
        dict_summary[region][time_period]['max_elev']=np.nanmean(dict_temp[:,1])

pdb.set_trace()
#Plot the inland expansion as a graph

#Display the keys
fig, axs = plt.subplots(2, 3)#, gridspec_kw={'width_ratios': [1, 3]})
fig.suptitle('Iceslabs inland progression')

axs = axs.ravel()

#count for subplot
i=0
#Loop over the region
for region in list(dict_summary.keys()):
    if (region == 'Out'):
        continue
    #Create empty vectors
    low_end=np.zeros(1)
    high_end=np.zeros(1)

    for time_period in list(dict_summary[region].keys()):
        low_end=np.append(low_end,dict_summary[region][time_period]['min_elev'])
        high_end=np.append(high_end,dict_summary[region][time_period]['max_elev'])
    
    #Remove zeros from low_end and high_end vectors
    low_end = low_end[~(low_end==0)]
    high_end = high_end[~(high_end==0)]
    
    #Plot the low end and high end of each region
    axs[i].plot(np.linspace(0,2,len(low_end)),low_end,label='Low end')
    axs[i].plot(np.linspace(0,2,len(high_end)),high_end,label='High end')
    
    #Set title
    axs[i].title.set_text(region)
    
    #Set x tick
    axs[i].set_xticks(np.linspace(0,2,len(high_end)))
    axs[i].set_xticklabels(list(dict_summary[region].keys()))
    
    axs[i].set_xlim(0,2)
    
    axs[i].grid()
    #Update count
    i=i+1
    
plt.legend()
plt.show()
#######################################################################
###          Inland expansion of iceslabs from 2002 to 2018         ###
#######################################################################

#Try the violin plot - Do not require any latitudinal and longitudinal averaging!

import seaborn as sns

sns.set_theme(style="whitegrid")

#Set the year for plotting
df_2002_2003_green['year']=["2002-2003" for x in range(len(df_2002_2003_green))]
df_MacFerrin['year']=["2010-2014" for x in range(len(df_MacFerrin))]
df_2017_2018['year']=["2017-2018" for x in range(len(df_2017_2018))]

#Append all the dataframes together
df_all=df_2002_2003_green
df_all=df_all.append(df_MacFerrin)
df_all=df_all.append(df_2017_2018)

#Prepare plot
fig, axs = plt.subplots(2, 3)#, gridspec_kw={'width_ratios': [1, 3]})
fig.suptitle('Iceslabs inland progression')

axs = axs.ravel()
i=0

for region in df_all['key_shp'].unique():
    if (region == 'Out'):
        continue
    
    #Add 2010-2017 and 2017-2018!
    sns.violinplot(ax=axs[i], data=df_all[df_all['key_shp']==region], x="year", y="elevation",
               inner="quart", linewidth=1)
    sns.despine(left=True)

    #Set title
    axs[i].title.set_text(region)
    
    axs[i].grid()
    #Update count
    i=i+1
                        
    '''
    catplot is noce but did not managed to to subplots with it
    sns.catplot(data=df_2002_2003_green[df_2002_2003_green['key_shp']==region], kind="violin", x="year", y="elevation", hue="colorcode_icelens",ax = axs[i])
    '''
    






