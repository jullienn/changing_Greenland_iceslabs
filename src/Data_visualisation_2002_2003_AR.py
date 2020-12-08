# -*- coding: utf-8 -*-
"""
Created on Sat Dec  5 11:34:11 2020

@author: Jullien Nicolas
"""

#Import packages

import scipy.io
import rasterio
from rasterio.plot import show
from matplotlib import pyplot
import numpy as np
import h5py
import matplotlib.colors as mcolors
import pandas as pd
from os import listdir
from os.path import isfile, join
import pdb
import pickle
import os.path
import os

from pysheds.grid import Grid

import osgeo.ogr as ogr
import osgeo.osr as osr

from pyproj import Transformer

##############################################################################
############################## Define variables ##############################
##############################################################################
dt = 2.034489716724874e-09 #Timestep for 2002/2003 traces
t0 = 0; # Unknown so set to zero

#Compute the speed (Modified Robin speed):
# self.C / (1.0 + (coefficient*density_kg_m3/1000.0))
v= 299792458 / (1.0 + (0.734*0.873/1000.0))
##############################################################################
############################## Define variables ##############################
##############################################################################

#Open the DEM
grid = Grid.from_raster("C:/Users/jullienn/Documents/working_environment/greenland_topo_data/elevations/greenland_dem_mosaic_100m_v3.0.tif",data_name='dem')
#Minnor slicing on borders to enhance colorbars
elevDem=grid.dem[:-1,:-1]              
#Scale the colormap
divnorm = mcolors.DivergingNorm(vmin=0, vcenter=1250, vmax=2500)
                             
#Define the working environment
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
            
            if (folder_day=='jun04'):
                print('Folder',folder_day,'. Not aggregated files, treat this one at the end. Continue ...')
                continue
            
            print('Now in year',folder_year,'day',folder_day)
            
            #Go into the daily folders 
            folder_day_name=folder_year_name+'/'+folder_day
            os.chdir(folder_day_name)
            
            # Read the files of this specific day
            onlyfiles = [f for f in listdir(folder_day_name) if isfile(join(folder_day_name, f))]
            #pdb.set_trace()
            for indiv_file in onlyfiles:
                print('Treating file',indiv_file)
                
                #If figure have already been generated, continue
                filename_to_check='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/2002_2003_data_localisation/'+indiv_file+'.png'
                if (os.path.isfile(filename_to_check)):
                    print('Figure already existent, move on to the next date')
                    continue
                
                #Open the file and read it
                f_agg = open(folder_day_name+'/'+indiv_file, "rb")
                data = pickle.load(f_agg)
                f_agg.close()
                
                #Select radar echogram and corresponding lat/lon
                radar_echo=data['radar_echogram']
                
                latlontime=data['latlontime']
                lat=latlontime['lat_gps']
                lon=latlontime['lon_gps']
                
                ###############################################################
                ##I. Process and plot radar echogram
                #
                ##Pick the surface!
                #
                ##Select the first 30m of radar echogram
                ##1. Compute the vertical resolution
                ##a. Time computation according to John Paden's email.
                #Nt = radar_echo.shape[0]
                #Time = t0 + dt*np.arange(1,Nt+1)
                ##b. Calculate the depth:
                ##self.SAMPLE_DEPTHS = self.radar_speed_m_s * self.SAMPLE_TIMES / 2.0
                #depths = v * Time / 2.0
                #
                ##2. Select the first 100 meters.
                #depths_100=depths[depths <= 100]
                #radar_echo_100=radar_echo[depths <= 100]
                #
                ##Plot the first 100m of radar echogram and lat/lon on map
                #
                ##ticks_yplot=np.around(np.linspace(0, 1400, 432))
                ##ticks_yplot=ticks_yplot.astype(int)
                ##labels_yplot=depths[ticks_yplot]
                #
                ##Should I use flipud or not??? np.flipud()
                #
                ##Change the size of the figure
                ##pyplot.rcParams["figure.figsize"]=30,30
                ##Plot the data
                #pyplot.figure()
                #color_map=pyplot.pcolor(radar_echo,cmap=pyplot.get_cmap('gray'))#,norm=divnorm)
                #pyplot.ylabel('Depth [m]')
                #pyplot.xlabel('Horizontal distance')
                ##pyplot.yticks(ticks=ticks_yplot,labels=labels_yplot)
                ##pyplot.xticks(fontsize=20)
                ##pyplot.yticks(fontsize=20)
                ##pyplot.ylim(0, 200)
                #pyplot.title('Radar echogram complete')
                #cbar=pyplot.colorbar()
                #cbar.set_label('Signal strength')
                #pyplot.show()
                #pdb.set_trace()
                ##Save the figure
                
                ###############################################################
                #II. Process and plot radar echogram localisation
                
                #Plot dem
                pyplot.figure(figsize=(48,40))
                pyplot.imshow(elevDem, extent=grid.extent,cmap='hot_r',norm=divnorm)
                pyplot.colorbar(label='Elevation [m]')
                pyplot.grid()
                
                #Plot radar track
                #1. reproject the track from WGS 84 to EPSG 3413             
               
                #if (indiv_file=='may18_02_10_aggregated'):
                #    pdb.set_trace()
                #    
                #if (indiv_file=='may18_02_18_aggregated'):
                #    pdb.set_trace()
                #
                #if (indiv_file=='may18_02_30_aggregated'):
                #    pdb.set_trace()
                #    
                #if (indiv_file=='may18_02_26_aggregated'):
                #    pdb.set_trace()

                #Some index have lat and lon equal to 0 because of jumps in data aggregation.
                #Replace these 0 by NaNs
                lat.replace(0, np.nan, inplace=True)
                lon.replace(0, np.nan, inplace=True)
                
                #Transform the longitudes. The longitudes are ~46 whereas they should be ~-46! So add a '-' in front of lon
                lon=-lon

                #pdb.set_trace()
                
                #Example from: https://pyproj4.github.io/pyproj/stable/examples.html
                transformer = Transformer.from_crs("EPSG:4326", "EPSG:3413", always_xy=True)
                points=transformer.transform(np.array(lon),np.array(lat))
                
                lon_3413=points[0]
                lat_3413=points[1]

                ###############################################################
                ############# Begin from IceBridgeGPR_Manager_v2.py ###########
                #'''Return the coordinates in North Polar Stereo projection, in (eastings, northings)'''
                #points = list(zip(lon, lat))
                ## 2. Convert all points to Polar Stereo grid projection
                #gcsSR = osr.SpatialReference()
                #gcsSR.SetWellKnownGeogCS("WGS84") # WGS84 geographical coordinates
                #npsSR = osr.SpatialReference()
                #npsSR.ImportFromEPSG(3413) # NSIDC North Pole Stereographic
                #GCS_TO_NPS = osr.CoordinateTransformation(gcsSR, npsSR)
        
                ## Transform to North Polar Stereo
                #nps_points = np.array(GCS_TO_NPS.TransformPoints(points))
                ## Cut off the zero "height" dimension
                #eastings  = np.array([p[0] for p in nps_points])
                #northings = np.array([p[1] for p in nps_points])
                ############# End from IceBridgeGPR_Manager_v2.py #############            
                ###############################################################
                
                #2. plot the tracks
                pyplot.scatter(lon_3413, lat_3413)
                pyplot.grid()
                #pyplot.show()
                
                #Create the figure name
                fig_name=[]
                fig_name='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/2002_2003_data_localisation/'+indiv_file+'.png'
                
                #Save the figure
                pyplot.savefig(fig_name)
                
                
                
    else:
        print('Folder',folder_year,', continue ...')
        continue

#Open files




#Visualise the data

            #Store everything into one dictionnary (matrix and vectors of data)
            #if (folder_day == 'jun04'):
            #    dic_file_being_read=[]
            #    dic_file_being_read = { "trace_id" : indiv_file.replace(".mat",""),
            #         "radar_echogram" : file_being_read['data'],
            #         "latlontime" : df_final}
            #else:
            #    dic_file_being_read=[]
            #    dic_file_being_read = { "trace_id" : indiv_file.replace(".mat",""),
            #         "radar_echogram" : file_being_read['filtfin'],
            #         "latlontime" : df_final}