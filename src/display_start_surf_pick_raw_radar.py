# -*- coding: utf-8 -*-
"""
Created on Sat Dec 16 10:30 2020

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

import matplotlib.gridspec as gridspec

##############################################################################
############################## Define variables ##############################
##############################################################################
dt = 2.034489716724874e-09 #Timestep for 2002/2003 traces
t0 = 0; # Unknown so set to zero

#Compute the speed (Modified Robin speed):
# self.C / (1.0 + (coefficient*density_kg_m3/1000.0))
v= 299792458 / (1.0 + (0.734*0.873/1000.0))

raw_radar_echograms='TRUE'

#Choose the year I want to diplay
year_display='2002'

##############################################################################
############################## Define variables ##############################
##############################################################################

##############################################################################
################### Define function for ice lenses logging ###################
##############################################################################
#This function if adapted from https://stackoverflow.com/questions/37363755/python-mouse-click-coordinates-as-simply-as-possible
def onclick(event):
    #This functions print and save the x and y coordinates in pixels!
    print(event.xdata, event.ydata)
    #Fill in the file to log on the information
    #filename_flog='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/2002_2003_radar_slice/flog_icelenses_alldates.txt'
    #f_log = open(filename_flog, "a")
    #f_log.write(str(round(event.xdata,2))+','+str(round(event.ydata,2))+'\n')
    #f_log.close() #Close the quality assessment file when we’re done!
##############################################################################
################### Define function for ice lenses logging ###################
##############################################################################

#Open, read and close the file of dates that require surface pick improvement for 2002/2003
f = open('C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/2002_2003_radar_raw_echogram/dates_for_surf_pick_start.txt','r')
dates_surf = [line.strip() for line in f.readlines() if len(line.strip()) > 0]
f.close()

#Open, read and close the file of dates that require surface pick improvement for 2017
f = open('C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/2002_2003_radar_raw_echogram/dates_for_surf_pick_2017.txt','r')
dates_surf_2017 = [line.strip() for line in f.readlines() if len(line.strip()) > 0]
f.close()

#Define the working environment
path= 'C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data'
os.chdir(path) # relative path: scripts dir is under Lab

# Read the years of data
folder_years = [ f.name for f in os.scandir(path) if f.is_dir() ]
#pdb.set_trace()
for folder_year in folder_years:
    
    if (year_display in list(['2002','2003'])):
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
                    
                    #If indiv_file is the quality file, continue
                    if (indiv_file[0:7]==('quality')):
                        #pdb.set_trace()
                        continue
                    
                    #If files does not deblong to 'dates_surf', the improvement of
                    #the start surf pick is not neccessary, continue
                    #pdb.set_trace()
                    
                    if (folder_day=='jun04'):
                        if (not(indiv_file.replace(".mat","") in dates_surf)):
                            print('No need to improve start surf pick of',indiv_file.replace(".mat",""))
                            continue
                    else:
                        if (not(indiv_file.replace("_aggregated","") in dates_surf)):
                            print('No need to improve start surf pick of',indiv_file.replace("_aggregated",""))
                            continue
                    
                    if (folder_day=='jun04'):
                        print('Folder is',folder_day,', Special treatment')
                        
                        fdata= scipy.io.loadmat(folder_day_name+'/'+indiv_file)
                        #Select radar echogram and corresponding lat/lon
                        radar_echo=fdata['data']
                        lat=fdata['latitude']
                        lon=fdata['longitude']
                        #pdb.set_trace()
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
                    
                    #Select the first 30m of radar echogram
                    #1. Compute the vertical resolution
                    #a. Time computation according to John Paden's email.
                    Nt = radar_echo.shape[0]
                    Time = t0 + dt*np.arange(1,Nt+1)
                    #b. Calculate the depth:
                    #self.SAMPLE_DEPTHS = self.radar_speed_m_s * self.SAMPLE_TIMES / 2.0
                    depths = v * Time / 2.0
                    
                    #If raw_radar_echograms is set to 'TRUE', then plot the raw
                    #radar echogram of that date and save it
                    if (raw_radar_echograms=='TRUE'):
                        ##If file have already been created, continue
                        #filename_to_check='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/2002_2003_radar_raw_echogram/'+indiv_file+'.png'
                        #if (os.path.isfile(filename_to_check)):
                        #    print('Figure already existent, move on to the next date')
                        #    continue
                        
                        #Generate the pick for vertical distance display
                        ticks_yplot=np.arange(0,radar_echo.shape[0],200)
                        
                        #Plot the raw radar echogram
                        fig=pyplot.figure(figsize=(48,40))
                        
                        #Change label font
                        pyplot.rcParams.update({'font.size': 40})
                        color_map=pyplot.pcolor(radar_echo[:,0:100],cmap=pyplot.get_cmap('gray'))#,norm=divnorm)
                        pyplot.gca().invert_yaxis() #Invert the y axis = avoid using flipud.
                        pyplot.yticks(ticks=ticks_yplot,labels=(np.round(depths[ticks_yplot])))
                        pyplot.ylabel('Depth [m]')
                        pyplot.xlabel('Horizontal distance')
                        pyplot.title('Raw radar echogram, first 100 horizontal pixels: '+indiv_file.replace("_aggregated",""))
                        cbar=pyplot.colorbar()
                        cbar.set_label('Signal strength')
                        
                        fig.canvas.mpl_connect('button_press_event', onclick)
                        pyplot.show()
    
                        pdb.set_trace()
                        
                        continue
                    
        else:
            print('Folder',folder_year,', continue ...')
            continue
                    
    if(year_display=='2017'):
        if (folder_year=='2017_Greenland_P3'):
            print('Treating the year',folder_year)
            #pdb.set_trace()
            
            #Go into the yearly folders 
            folder_year_name=path+'/'+folder_year+'/CSARP_qlook'
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
                    
                    #If files does not belong to 'dates_surf_2017', we do not need
                    #to know the start surf pick, continue
                    #pdb.set_trace()
                    if (not(indiv_file.replace(".mat","") in dates_surf_2017)):
                        print('No need to improve start surf pick of',indiv_file)
                        continue
                    
                    #Open the file and read it
                    with h5py.File(folder_day_name+'/'+indiv_file, 'r') as f:
                        f.keys()
                        #Select radar echogram
                        radar_echo=f['Data'][:].transpose() #2017 data should be transposed
                    
                    #If raw_radar_echograms is set to 'TRUE', then plot the raw
                    #radar echogram of that date and save it
                    if (raw_radar_echograms=='TRUE'):
                        
                        #Generate the pick for vertical distance display
                        ticks_yplot=np.arange(0,radar_echo.shape[0],200)
                        
                        #Plot the raw radar echogram
                        fig=pyplot.figure(figsize=(48,40))
                        
                        #Change label font
                        pyplot.rcParams.update({'font.size': 20})
                        color_map=pyplot.pcolor(np.log10(radar_echo[:,0:100]),cmap=pyplot.get_cmap('gray'))#,norm=divnorm)
                        pyplot.gca().invert_yaxis() #Invert the y axis = avoid using flipud.
                        #pyplot.yticks(ticks=ticks_yplot,labels=(np.round(depths[ticks_yplot])))
                        pyplot.ylabel('Depth [pixel]')
                        pyplot.xlabel('Horizontal distance')
                        pyplot.title('Log10(raw radar echogram), first 100 horizontal pixels: '+indiv_file.replace(".mat",""))
                        cbar=pyplot.colorbar()
                        cbar.set_label('Signal strength')
                        
                        fig.canvas.mpl_connect('button_press_event', onclick)
                        pyplot.show()
    
                        pdb.set_trace()
                        
                        
                        continue

        else:
            print('Folder',folder_year,', continue ...')
            continue