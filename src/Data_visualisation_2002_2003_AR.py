# -*- coding: utf-8 -*-
"""
Created on Sat Dec  5 11:34:11 2020

@author: Jullien Nicolas
"""

#Import packages

import scipy.io
#import rasterio
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
            pdb.set_trace()
            for indiv_file in onlyfiles:
                print('Treating file',indiv_file)
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
                
                #Pick the surface!
                
                #Select the first 30m of radar echogram
                #1. Compute the vertical resolution
                #a. Time computation according to John Paden's email.
                Nt = radar_echo.shape[0]
                Time = t0 + dt*np.arange(1,Nt+1)
                #b. Calculate the depth:
                #self.SAMPLE_DEPTHS = self.radar_speed_m_s * self.SAMPLE_TIMES / 2.0
                depths = v * Time / 2.0
                
                #2. Select the first 100 meters.
                depths_100=depths[depths <= 100]
                radar_echo_100=radar_echo[depths <= 100]
                
                #Plot the first 100m of radar echogram and lat/lon on map
                
                
                #ticks_yplot=np.around(np.linspace(0, 1400, 432))
                #ticks_yplot=ticks_yplot.astype(int)
                #labels_yplot=depths[ticks_yplot]
                
                #Change the size of the figure
                #pyplot.rcParams["figure.figsize"]=30,30
                #Plot the data
                pyplot.figure()
                color_map=pyplot.pcolor(radar_echo_100,cmap=pyplot.get_cmap('gray'))#,norm=divnorm)
                pyplot.ylabel('Depth [m]')
                pyplot.xlabel('Horizontal distance')
                #pyplot.yticks(ticks=ticks_yplot,labels=labels_yplot)
                #pyplot.xticks(fontsize=20)
                #pyplot.yticks(fontsize=20)
                #pyplot.ylim(0, 200)
                pyplot.title('Radar echogram')
                cbar=pyplot.colorbar()
                cbar.set_label('Signal strength')
                pyplot.show()
                
                #Save the figure
                
        
        
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