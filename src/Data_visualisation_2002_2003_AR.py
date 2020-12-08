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


##############################################################################
############# Define kernel function for surface identification ##############
##############################################################################
#_gaussian function taken from IceBridgeGPR_Manager_v2.py
# Define a quick guassian function to scale the cutoff mask above
def _gaussian(x,mu,sigma):
    return numpy.exp(-numpy.power((x-mu)/sigma, 2.)/2.)

#This function have been taken from 'IceBridgeGPR_Manager_v2.py
def kernel_function(traces_input,suggested_pixel):
    traces = traces_input
    traces = np.log10(traces)
    # We do not have the original indicies to use as a starter
    
    # 3) Perform surface pick crawling threshold behavior mask (assume a step-change analysis [goes from weak->strong at surface], and continuity of surface in original file.)
    # Create a step-change mask to optimze where the returns transition from "dark" to "bright"
    MASK_RADIUS = 50
    vertical_span_mask = numpy.empty([MASK_RADIUS*2,], dtype=numpy.float)
    vertical_span_mask[:MASK_RADIUS] = -1.0
    vertical_span_mask[MASK_RADIUS:] = +3.0
    
    vertical_span_mask = vertical_span_mask * _gaussian(numpy.arange(vertical_span_mask.shape[0]),mu=(MASK_RADIUS-5),sigma=(float(MASK_RADIUS)/3.0))
    
    # Expand the shape to handle array broadcasting below
    vertical_span_mask.shape = vertical_span_mask.shape[0], 1
    
    # This is the vertical window size of the extent of the search.  Should be bigger than any jump from one surface pixel to the next.
    MASK_SEARCH_RADIUS = 150
    
    improved_indices = numpy.empty(traces.shape[1], dtype='int64')
    # CHECK THAT traces.shape[1] CORRESPONDS TO HORIZONTAL DISTANCE!
    
    # Start at the left with the hand-picked "suggested surface pick" in the ICEBRIDGE_SURFACE_PICK_SUGGESTIONS_FILE as starting point
    
    last_best_index = suggested_pixel
     
    pdb.set_trace()
    # A template graph to use, just have to add in the center vertical index at each point and go from there.
    search_indices_template = numpy.sum(numpy.indices((vertical_span_mask.shape[0], 2*MASK_SEARCH_RADIUS)),axis=0) - MASK_SEARCH_RADIUS - MASK_RADIUS
    for i in range(traces.shape[1]):
        # Create an array of indices spanning the top-to-bottom of the MASK_SEARCH_RADIUS, and fanning out MASK_RADIUS above and below that point.
        search_indices = search_indices_template + last_best_index
        # Handle overflow indices if below zero or above max (shouldn't generally happen)... just assign to the top or bottom pixel
        search_indices[search_indices < 0] = 0
        search_indices[search_indices >= traces.shape[0]] = traces.shape[0]-1
        
        bestfit_sum = numpy.sum(traces[:,i][search_indices] * vertical_span_mask, axis=0)
        
        assert bestfit_sum.shape[0] == 2*MASK_SEARCH_RADIUS
        
        # Get the best fit (with the highest value from the transformation fit)
        last_best_index = search_indices[MASK_RADIUS,numpy.argmax(bestfit_sum)]
        improved_indices[i] = last_best_index
        
        #If there are pixels with particularly strong echo that are being erroneously
        #picked up as the surface, then I could do that:
        ###### Erase most the little "jump" artifacts in the surface picker.
        ##### improved_indices = self._get_rid_of_false_surface_jumps(improved_indices)
        
        #I do not use any mask so I think I shouldn't need to use that:
        ###### Must re-expand the surface indices to account for masked values (filled w/ nan)
        ##### improved_indices_expanded = self._refill_array(improved_indices, surface_maskname)
        
    return improved_indices
    
##############################################################################
############# Define kernel function for surface identification ##############
##############################################################################

#Open the DEM
grid = Grid.from_raster("C:/Users/jullienn/Documents/working_environment/greenland_topo_data/elevations/greenland_dem_mosaic_100m_v3.0.tif",data_name='dem')
#Minnor slicing on borders to enhance colorbars
elevDem=grid.dem[:-1,:-1]              
#Scale the colormap
divnorm = mcolors.DivergingNorm(vmin=0, vcenter=1250, vmax=2500)

#Open, read and close the file of suggested surface picks
f = open('C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/Exclusion_folder/txt/SURFACE_STARTING_PICKS_Suggestions_2002_2003.txt','r')
lines = [line.strip() for line in f.readlines() if len(line.strip()) > 0]
f.close()
                             
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
                
                #If indiv_file is the quality file, continue
                if (indiv_file[0:7]==('quality')):
                    #pdb.set_trace()
                    continue
                
                ##If figure have already been generated, continue
                #filename_to_check='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/2002_2003_data_localisation/'+indiv_file+'.png'
                #if (os.path.isfile(filename_to_check)):
                #    print('Figure already existent, move on to the next date')
                #    continue
                
                pdb.set_trace()
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
                #I. Process and plot radar echogram
                
                #Pick the surface!
                #There is no 'Surface' variable in 2002/2003 dataset such as 2010/2014 datset. I have to overcome this issue.
                
                # Load the suggested pixel for the specific date
                for date_pix in lines:
                    pdb.set_trace()
                    if (date_pix[0:10]==str(indiv_file.replace("_aggregated",""))):
                        suggested_pixel=int(date_pix[12:16])
                
                #index_pos = lines[0][0:10].index(str(indiv_file.replace("_aggregated","")))
                #lines[:][0:10].index(str('may18_02_4'))
                

                #if indiv_file in lines:
                #    lines
                #    date_and_pixel=lines[count_ite]
                    
                
                
                #Call the kernel_function to compute the surface
                surface_indices=kernel_function(radar_echo, suggested_pixel)
                
                # Get our slice
                #radar_slice = self._return_radar_slice_given_surface(traces,
                                                             #surface_indices,
                                                             #meters_cutoff_above=meters_cutoff_above,
                                                             #meters_cutoff_below=meters_cutoff_below)
                
                #Select the first 30m of radar echogram
                #1. Compute the vertical resolution
                #a. Time computation according to John Paden's email.
                Nt = radar_echo.shape[0]
                Time = t0 + dt*np.arange(1,Nt+1)
                #b. Calculate the depth:
                #self.SAMPLE_DEPTHS = self.radar_speed_m_s * self.SAMPLE_TIMES / 2.0
                depths = v * Time / 2.0
                
                #2. Select the first 100 meters.
                #depths_100=depths[depths <= 100]
                #radar_echo_100=radar_echo[depths <= 100]
                
                #Plot the first 100m of radar echogram and lat/lon on map
                pdb.set_trace()
                
                #ticks_yplot=np.around(np.linspace(0, 1400, 432))
                #ticks_yplot=ticks_yplot.astype(int)
                #labels_yplot=depths[ticks_yplot]
                
                ##############################################################
                ############### Begin explanations on pcolor #################
                
                #Explainations on how pcolor works. I convinced myself doing a small example that I plotted.
                #Further explanations: https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.pcolor.html#:~:text=pcolor()%20displays%20all%20columns,Similarly%20for%20the%20rows.
                #Example: I have a matrix Z that I want to plot
                
                #Z=[10,3,24,70,                            40,5,48,22
                #    2,6,87,21,      ----pcolor(Z)--->     2,6,87,21
                #    40,5,48,22]                           10,3,24,70
                
                #I must use np.flipud() to display the data from top to bottom.
                
                ################ End explanations on pcolor ##################
                ##############################################################
                
                #Change the size of the figure
                #pyplot.rcParams["figure.figsize"]=30,30
                #Plot the data
                pyplot.figure()
                color_map=pyplot.pcolor(np.flipud(radar_echo),cmap=pyplot.get_cmap('gray'))#,norm=divnorm)
                pyplot.ylabel('Depth [m]')
                pyplot.xlabel('Horizontal distance')
                #pyplot.yticks(ticks=ticks_yplot,labels=labels_yplot)
                #pyplot.xticks(fontsize=20)
                #pyplot.yticks(fontsize=20)
                #pyplot.ylim(0, 200)
                pyplot.title('Radar echogram complete: '+indiv_file[0:10])
                cbar=pyplot.colorbar()
                cbar.set_label('Signal strength')
                pyplot.show()
                pdb.set_trace()
                #Save the figure
                
                ###############################################################
                #II. Process and plot radar echogram localisation
                
                ##Plot dem
                #pyplot.figure(figsize=(48,40))
                #pyplot.imshow(elevDem, extent=grid.extent,cmap='hot_r',norm=divnorm)
                #pyplot.colorbar(label='Elevation [m]')
                #pyplot.grid()
                #
                ##Plot radar track
                ##1. reproject the track from WGS 84 to EPSG 3413             
                #
                ##Some index have lat and lon equal to 0 because of jumps in data aggregation.
                ##Replace these 0 by NaNs
                #lat.replace(0, np.nan, inplace=True)
                #lon.replace(0, np.nan, inplace=True)
                #
                ##Transform the longitudes. The longitudes are ~46 whereas they should be ~-46! So add a '-' in front of lon
                #lon=-lon
                #
                ##Transform the coordinated from WGS84 to EPSG:3413
                ##Example from: https://pyproj4.github.io/pyproj/stable/examples.html
                #transformer = Transformer.from_crs("EPSG:4326", "EPSG:3413", always_xy=True)
                #points=transformer.transform(np.array(lon),np.array(lat))
                #
                #lon_3413=points[0]
                #lat_3413=points[1]
                #
                ##2. plot the tracks
                #pyplot.scatter(lon_3413, lat_3413)
                #pyplot.grid()
                ##pyplot.show()
                #
                ##Create the figure name
                #fig_name=[]
                #fig_name='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/2002_2003_data_localisation/'+indiv_file+'.png'
                #
                ##Save the figure
                #pyplot.savefig(fig_name)
                #pyplot.clf()
 
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