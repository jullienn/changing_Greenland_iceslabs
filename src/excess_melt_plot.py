# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 14:29:11 2021

@author: JullienN

Read and plot the excess melt calculation made by Andrew Tedstone
"""

#Import librairies

#Import packages

import matplotlib.pyplot as plt
import rasterio as rio
import geopandas as gpd
from os import listdir
import os
from os.path import isfile, join
import numpy as np
import xarray as xr
import pickle
import pdb
import pandas as pd
import pyproj
import scipy.io


#Define path of the working directory
global_path='C:/Users/jullienn/Documents/working_environment/'

#Define the desired year to plot
desired_year='2006'
generate_raw_excess_melt='FALSE'
generate_excess_melt_traces='TRUE'

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

##########################################################################
###                      Load Greenland DEM and contours               ###
##########################################################################

dem_filename = 'C:/Users/jullienn/Documents/working_environment/greenland_topo_data/elevations/greenland_dem_mosaic_100m_v3.0.tif'
contours_filename = 'C:/Users/jullienn/Documents/working_environment/greenland_topo_data/elevations/greenland_contours_100m_v3.0.shp'

with rio.open(dem_filename) as fh:
    dem = fh.read()
    print(fh)
    dem_bounds = fh.bounds
    dem_crs = fh.crs
contours = gpd.read_file(contours_filename)
# Make sure that the coordinate reference system matches the image
contours = contours.to_crs(dem_crs)
         
##########################################################################
###                      Load Greenland DEM and contours               ###
##########################################################################


##########################################################################
###                          Load excess melt data   	               ###
##########################################################################
#Define path
path='C:/Users/jullienn/Documents/working_environment/excess_melt/'

#Load the data
data_path = path+'excess_melt_mbyear.nc'
DS = xr.open_dataset(data_path)

#Extract coordinates
lat_M_e=DS.x.data
lon_M_e=DS.y.data

#Load melt data
melt_data= DS.M_e
##########################################################################
###                          Load excess melt data   	               ###
##########################################################################

if (generate_excess_melt_traces=='TRUE'):
    #Create the figures of old dataset traces on top of annual excess melt figures
    
    #Define the working environment for loading the data
    path= 'C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data'
    os.chdir(path) # relative path: scripts dir is under Lab
    
    # Read the years of data
    folder_years = [ f.name for f in os.scandir(path) if f.is_dir() ]
    
    for folder_year in folder_years:
        if (folder_year in list(['2002','2003'])):
            print('Treating the year',folder_year)
            
            #Go into the yearly folders 
            folder_year_name=path+'/'+folder_year
            os.chdir(folder_year_name)
    
            # Read the days of this specific year
            folder_days = [ f.name for f in os.scandir(folder_year_name) if f.is_dir() ]
            
            for folder_day in folder_days:
                
                print('Now in year',folder_year,'day',folder_day)
                
                #Go into the daily folders 
                folder_day_name=folder_year_name+'/'+folder_day
                os.chdir(folder_day_name)
                
                # Read the files of this specific day
                onlyfiles = [f for f in listdir(folder_day_name) if isfile(join(folder_day_name, f))]
                #pdb.set_trace()
                for indiv_file in onlyfiles:
                    print('Treating file',indiv_file)
    
                    #pdb.set_trace()
                    #If indiv_file is the quality file, continue
                    if (indiv_file[0:7]==('quality')):
                        #pdb.set_trace()
                        continue
                    
                    if (folder_day=='jun04'):
                        
                        fdata= scipy.io.loadmat(folder_day_name+'/'+indiv_file)
                        #Select radar echogram and corresponding lat/lon
                        radar_echo=fdata['data']
                        lat=fdata['latitude']
                        lon=fdata['longitude']
                        
                        #Transform lat and lon variables into series to be able
                        #to create the pandas dataframes:
                        lat=pd.Series(lat.flatten())
                        lon=pd.Series(lon.flatten())
    
                    else:
                        #Open the file and read it
                        f_agg = open(folder_day_name+'/'+indiv_file, "rb")
                        data = pickle.load(f_agg)
                        f_agg.close()
                                        
                        #Select radar echogram and corresponding lat/lon
                        radar_echo=data['radar_echogram']
                        
                        latlontime=data['latlontime']
                        lat=latlontime['lat_gps']
                        lon=latlontime['lon_gps']

                    #II. Plot radar echogram localisation
                    #II.a. Plot radar track
                    #II.a.1. Reproject the track from WGS 84 to EPSG 3413
                    #Some index have lat and lon equal to 0 because of jumps in data aggregation.
                    #Replace these 0 by NaNs
                    
                    if (not(folder_day=='jun04')):
                        lat.replace(0, np.nan, inplace=True)
                        lon.replace(0, np.nan, inplace=True)
                    
                    #Transform the longitudes. The longitudes are ~46 whereas they should be ~-46! So add a '-' in front of lon
                    lon=-lon
                    
                    ######################################################
                    ###     Convert lat/lon into MAR's projection      ###
                    #This piece of code is adapted from Andrew
                    
                    #1. Create a dataframe of traces with lat/lon in EPSG:4326
                    df_latlon=pd.DataFrame({'lon':pd.Series(lon),
                                            'lat':pd.Series(lat)})
                    
                    #2. Create the geodataframe
                    df_latlon_gdf = gpd.GeoDataFrame(df_latlon, geometry=gpd.points_from_xy(df_latlon.lon, df_latlon.lat), crs=4326)
                    
                    #3. Convert the data from EPSG:4326 to MAR's projection
                    mar_proj4 = '+proj=stere +lat_0=70.5 +lon_0=-40 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'
                    df_latlon_gdf = df_latlon_gdf.to_crs(mar_proj4)
                    
                    ###     Convert lat/lon into MAR's projection      ###
                    ######################################################

                    ######################################################                        
                    ###           Load excess melt data                ###
                    
                    #Select the data associated with the wanted year (desired_year)
                    melt_year = melt_data.sel(time=desired_year)
                    melt_year_np = melt_year.values
                    
                    melt_year_plot=np.asarray(melt_year_np)
                    melt_year_plot=melt_year_plot[0,0,:,:,0]
                    
                    ###           Load excess melt data                ###
                    ######################################################                        


                    #Plot excess melt data with traces over it
                    plt.rcParams.update({'font.size': 20})
                    plt.figure(figsize=(48,40))
                    ax = plt.subplot(111)
                    dem_extent = (dem_bounds[0], dem_bounds[2], dem_bounds[1], dem_bounds[3])
                    #plot excess melt
                    plt.imshow(np.squeeze(np.flipud(melt_year_plot)), extent=dem_extent,cmap=discrete_cmap(5,'hot_r'))
                    plt.colorbar(label='Excess melt [mm w.e./year]')
                    plt.clim(0,1000)
                    ax.grid()
                    
                    #plot tracks
                    df_latlon_gdf.plot(ax=ax, marker='o', markersize=1, zorder=45, color='blue')
                    
                    #Display begining of tracks
                    begining_traces=df_latlon_gdf.loc[0:10]
                    begining_traces.plot(ax=ax, marker='o', markersize=1, zorder=45, color='magenta')
                    
                    #contours.plot(ax=ax, edgecolor='black')
                    plt.title(indiv_file.replace("_aggregated","")+' with excess melt year: '+desired_year)
                    
                    plt.show()
                    pdb.set_trace()

 
                    #if (folder_day=='jun04'):
                    #    ax1.set_xlim(lon_3413[0,0]-500000, lon_3413[0,0]+500000)
                    #    ax1.set_ylim(lat_3413[0,0]-500000, lat_3413[0,0]+500000)
                    #else:
                    #    ax1.set_xlim(lon_3413[0]-500000, lon_3413[0]+500000)
                    #    ax1.set_ylim(lat_3413[0]-500000, lat_3413[0]+500000)


                    ##II.a.4. Plot the tracks
                    #ax.scatter(lon_3413, lat_3413,s=0.1)
                    #
                    #if (folder_day=='jun04'):
                    #    ax1.scatter(lon_3413[0,0],lat_3413[0,0],c='m',s=0.1) #Plot the start in green
                    #else:
                    #    ax1.scatter(lon_3413[0],lat_3413[0],c='m',s=0.1)
                    
                    #ax1.grid()



                    
                    ##Create the figure name
                    #fig_name=[]
                    #fig_name='C:/Users/jullienn/Documents/working_environment/excess_melt/figures_excess_melt/excess_melt_'+wanted_year+'.png'
                    
                    ##Save the figure
                    #plt.savefig(fig_name)
                    #plt.clf()




if (generate_raw_excess_melt=='TRUE'):
    #Generate and save the raw annual excess melt figures
    for year in list(np.arange(1990,2020)):
        
        #Define the year
        wanted_year=str(year)
        
        #Select the data associated with the wanted year
        melt_year = melt_data.sel(time=wanted_year)
        melt_year_np = melt_year.values
        
        melt_year_plot=np.asarray(melt_year_np)
        melt_year_plot=melt_year_plot[0,0,:,:,0]
        
        #Plot dem and contours elevation
        plt.rcParams.update({'font.size': 20})
        plt.figure(figsize=(48,40))
        ax = plt.subplot(111)
        dem_extent = (dem_bounds[0], dem_bounds[2], dem_bounds[1], dem_bounds[3])
        plt.imshow(np.squeeze(np.flipud(melt_year_plot)), extent=dem_extent,cmap=discrete_cmap(5,'hot_r'))
        
        plt.colorbar(label='Excess melt [mm w.e./year]')
        plt.clim(0,1000)
        ax.grid()
        #contours.plot(ax=ax, edgecolor='black')
        plt.title('Excess melt plot, year: '+wanted_year)
    
        #Create the figure name
        fig_name=[]
        fig_name='C:/Users/jullienn/Documents/working_environment/excess_melt/figures_excess_melt/excess_melt_'+wanted_year+'.png'
                            
        #Save the figure
        plt.savefig(fig_name)
        plt.clf()








