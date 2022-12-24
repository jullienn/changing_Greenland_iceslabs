# -*- coding: utf-8 -*-
"""
Created on Sat Dec 16 10:30 2020

@author: Jullien Nicolas
"""

#Import packages

import scipy.io
import rasterio
from rasterio.plot import show
import matplotlib.pyplot as plt
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
import cartopy.crs as ccrs
import geopandas as gpd
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
year_display='2012'

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


#### ---------- Investigation NO and NW 2010-12 data availability --------- ###

'''
20120510_01_046 OK -
20120510_01_075 OK -
20120510_01_062 OK -
20120510_01_030 OK -
20120510_01_089 OK -
20120330_01_007 OK -
20110511_01_166 OK -
20110511_01_109 OK -
20110511_01_096 OK -
'''

###################### From Tedstone et al., 2022 #####################
#from plot_map_decadal_change.py
# Define the CartoPy CRS object.
crs = ccrs.NorthPolarStereo(central_longitude=-45., true_scale_latitude=70.)
# This can be converted into a `proj4` string/dict compatible with GeoPandas
crs_proj4 = crs.proj4_init
###################### From Tedstone et al., 2022 #####################

#Transform coordinates from WGS84 to EPSG:3413
transformer = Transformer.from_crs("EPSG:4326", "EPSG:3413", always_xy=True)

### -------------------------- Load shapefiles --------------------------- ###
#Load Rignot et al., 2016 Greenland drainage bassins
path_rignotetal2016_GrIS_drainage_bassins='C:/Users/jullienn/switchdrive/Private/research/backup_Aglaja/working_environment/greenland_topo_data/GRE_Basins_IMBIE2_v1.3/'
GrIS_drainage_bassins=gpd.read_file(path_rignotetal2016_GrIS_drainage_bassins+'GRE_Basins_IMBIE2_v1.3_EPSG_3413.shp') #the regions are the last rows of the shapefile

path_df_with_elevation='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2002_2018/' 
#Load corresponding shapefile
iceslabs_jullien_highend_20102018=gpd.read_file(path_df_with_elevation+'shapefiles/iceslabs_jullien_highend_20102018.shp') 

#Prepare plot
fig = plt.figure()
fig.set_size_inches(14, 10) # set figure's size manually to your full screen (32x18), this is from https://stackoverflow.com/questions/32428193/saving-matplotlib-graphs-to-image-as-full-screen
#projection set up from https://stackoverflow.com/questions/33942233/how-do-i-change-matplotlibs-subplot-projection-of-an-existing-axis
ax1 = plt.subplot(projection=crs)

GrIS_drainage_bassins.plot(ax=ax1,color='none', edgecolor='black',linewidth=0.5)
iceslabs_jullien_highend_20102018.plot(ax=ax1,color='red', edgecolor='black',linewidth=0.5)


path_data='C:/Users/jullienn/Downloads/'
path_save='C:/Users/jullienn/switchdrive/Private/research/RT1/manuscript/revisions/v4/investigation_NW_Greenland/'

if(year_display=='2012'):

    #onlyfiles=['Data_img_01_20120511_01_031.mat','Data_img_01_20120511_01_032.mat','Data_img_01_20120511_01_038.mat','Data_img_01_20120511_01_039.mat','Data_img_01_20120511_01_040.mat']
    
    #onlyfiles=['Data_img_01_20120510_01_045.mat','Data_img_01_20120510_01_046.mat','Data_img_01_20120510_01_047.mat']
    #onlyfiles=['Data_img_01_20120510_01_074.mat','Data_img_01_20120510_01_075.mat','Data_img_01_20120510_01_076.mat']
    #onlyfiles=['Data_img_01_20120510_01_061.mat','Data_img_01_20120510_01_062.mat','Data_img_01_20120510_01_063.mat']
    #onlyfiles=['Data_img_01_20120510_01_029.mat','Data_img_01_20120510_01_030.mat','Data_img_01_20120510_01_031.mat']
    #onlyfiles=['Data_img_01_20120510_01_088.mat','Data_img_01_20120510_01_089.mat','Data_img_01_20120510_01_090.mat']
    #onlyfiles=['Data_img_01_20120330_01_006.mat','Data_img_01_20120330_01_007.mat','Data_img_01_20120330_01_008.mat']
    #onlyfiles=['Data_20110511_01_165.mat','Data_20110511_01_166.mat','Data_20110511_01_167.mat']
    #onlyfiles=['Data_20110511_01_108.mat','Data_20110511_01_109.mat','Data_20110511_01_110.mat']
    onlyfiles=['Data_20110511_01_095.mat','Data_20110511_01_096.mat','Data_20110511_01_097.mat']



    for indiv_file in onlyfiles:
        print('Treating file',indiv_file)

        
        #Open the file and read it
        fdata= scipy.io.loadmat(path_data+indiv_file)
        #Select radar echogram and corresponding lat/lon
        radar_echo=fdata['Data']
        
        lat_data=fdata['Latitude']
        lon_data=fdata['Longitude']
        points=transformer.transform(np.asarray(lon_data),np.asarray(lat_data))
        lon_3413=points[0]
        lat_3413=points[1]
        
        #Display data loc on map
        ax1.scatter(lon_3413,lat_3413,s=0.1)
        ax1.text(np.median(lon_3413),np.median(lat_3413),indiv_file[-7:-4],size=20)

        time=fdata['Time']
        surface=fdata['Surface']
        
        depths = v * time / 2.0
        depths=depths+np.abs(depths[0])
        depths_20m=depths[depths<=20]
        
        index=np.zeros((surface.shape[0],surface.shape[1]))
        radar_20m=np.zeros((len(depths_20m),surface.shape[1]))
        radar_20m[:]=np.nan
        
        for i in range(0,surface.shape[1]):
            index[0,i]=np.where(surface[0,i]==time)[0][0]
            radar_20m[:,i]=radar_echo[index[0][i].astype(int):index[0][i].astype(int)+len(depths_20m),i]

        #If raw_radar_echograms is set to 'TRUE', then plot the raw
        #radar echogram of that date and save it
        if (raw_radar_echograms=='TRUE'):
            '''
            #Generate the pick for vertical distance display
            ticks_yplot=np.arange(0,radar_30m.shape[0],200)
            '''

            fig = plt.figure()
            ax2 = plt.subplot()

            color_map=ax2.pcolor(np.log10(radar_20m),cmap=plt.get_cmap('gray'))#,norm=divnorm)
            plt.gca().invert_yaxis() #Invert the y axis = avoid using flipud.
            #pyplot.yticks(ticks=ticks_yplot,labels=(np.round(depths[ticks_yplot])))

            plt.title(indiv_file.replace(".mat",""))
            ax2.set_aspect(2)
            
            #Maximize plot size - This is from Fig1.py from Grenland ice slabs expansion and thickening paper.
            figManager = plt.get_current_fig_manager()
            figManager.window.showMaximized()
            
            plt.show()

            #Save figure
            plt.savefig(path_save+indiv_file.replace(".mat",""),dpi=300,bbox_inches='tight')
            #bbox_inches is from https://stackoverflow.com/questions/32428193/saving-matplotlib-graphs-to-image-as-full-screen)

#### ---------- Investigation NO and NW 2010-12 data availability --------- ###
pdb.set_trace()

plt.savefig(path_save+'map_'+indiv_file.replace(".mat",""),dpi=300,bbox_inches='tight')


##Open, read and close the file of dates that require surface pick improvement for 2002/2003
#f = open('C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/2002_2003_radar_raw_echogram/dates_for_surf_pick_start.txt','r')
#dates_surf = [line.strip() for line in f.readlines() if len(line.strip()) > 0]
#f.close()

#Open, read and close the file of dates that require surface pick improvement for 2017
f = open('C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/dates_for_surf_pick_2017.txt','r')
dates_surf_2017 = [line.strip() for line in f.readlines() if len(line.strip()) > 0]
f.close()

f = open('C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/dates_for_surf_pick_2018.txt','r')
dates_surf_2018 = [line.strip() for line in f.readlines() if len(line.strip()) > 0]
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
                    
                    ##If files does not belong to 'dates_surf_2017', we do not need
                    ##to know the start surf pick, continue
                    ##pdb.set_trace()
                    #if (not(indiv_file.replace(".mat","") in dates_surf_2017)):
                    #    print('No need to improve start surf pick of',indiv_file)
                    #    continue
                    
                    if (not(indiv_file.replace(".mat","")=='Data_20170417_01_170')):
                        continue
                    
                    #Open the file and read it
                    with h5py.File(folder_day_name+'/'+indiv_file, 'r') as f:
                        f.keys()
                        #Select radar echogram
                        radar_echo=f['Data'][:].transpose() #2017 data should be transposed
                    
                    #If raw_radar_echograms is set to 'TRUE', then plot the raw
                    #radar echogram of that date and save it
                    if (raw_radar_echograms=='TRUE'):
                        pdb.set_trace()
                        #Generate the pick for vertical distance display
                        ticks_yplot=np.arange(0,radar_echo.shape[0],200)
                        
                        #Plot the raw radar echogram
                        fig=pyplot.figure(figsize=(48,40))
                        
                        #Change label font
                        pyplot.rcParams.update({'font.size': 10})
                        color_map=pyplot.pcolor(np.log10(radar_echo[500:2000,:]),cmap=pyplot.get_cmap('gray'))#,norm=divnorm)
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

    if(year_display=='2018'):
        if (folder_year=='2018_Greenland_P3'):
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
                    
                    ##If files does not belong to 'dates_surf_2017', we do not need
                    ##to know the start surf pick, continue
                    ##pdb.set_trace()
                    #if (not(indiv_file.replace(".mat","") in dates_surf_2018)):
                    #    print('No need to improve start surf pick of',indiv_file)
                    #    continue
                    
                    if (not(indiv_file.replace(".mat","")=='Data_20180405_01_161')):
                        continue
                    
                    #Open the file and read it
                    with h5py.File(folder_day_name+'/'+indiv_file, 'r') as f:
                        
                        f.keys()
                        #Select radar echogram
                        radar_echo=f['Data'][:].transpose() #2017 data should be transposed
                        
                        ### Begin investigation roll and pitch of aircraft
                        roll=f['Roll'][:]
                        heading=f['Heading'][:]
                        
                        fig, (ax1, ax2) = pyplot.subplots(2, 1)#, gridspec_kw={'width_ratios': [1, 3]})
                        fig.suptitle('Investigation roll')
                        ax1.plot(np.arange(0,roll.size),roll)
                        ax1.set_title('Roll')
                        ax2.plot(np.arange(0,heading.size),heading)
                        ax2.set_title('Heading')
                        pyplot.show()
                        
                        pdb.set_trace()
                        
                        ### End investigation roll and pitch of aircraft
                    
                    #If raw_radar_echograms is set to 'TRUE', then plot the raw
                    #radar echogram of that date and save it
                    if (raw_radar_echograms=='TRUE'):
                        
                        #Generate the pick for vertical distance display
                        ticks_yplot=np.arange(0,radar_echo.shape[0],200)
                        
                        #Plot the raw radar echogram
                        fig=pyplot.figure(figsize=(48,40))
                        
                        #Change label font
                        pyplot.rcParams.update({'font.size': 10})
                        color_map=pyplot.pcolor(np.log10(radar_echo[500:2000,:]),cmap=pyplot.get_cmap('gray'))#,norm=divnorm)
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

    if(year_display=='2011'):
        if (folder_year=='2011_Greenland_P3'):
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
                    
                    ##If files does not belong to 'dates_surf_2017', we do not need
                    ##to know the start surf pick, continue
                    #if (not(indiv_file.replace(".mat","") == 'Data_20110329_01_018')):
                    #    print('No need to improve start surf pick of',indiv_file)
                    #    continue
                    pdb.set_trace()
                    
                    #Open the file and read it
                    fdata= scipy.io.loadmat(folder_day_name+'/'+indiv_file)
                    #Select radar echogram and corresponding lat/lon
                    radar_echo=fdata['Data']
                    
                    #If raw_radar_echograms is set to 'TRUE', then plot the raw
                    #radar echogram of that date and save it
                    if (raw_radar_echograms=='TRUE'):
                        
                        #Generate the pick for vertical distance display
                        ticks_yplot=np.arange(0,radar_echo.shape[0],200)
                        
                        #Plot the raw radar echogram
                        fig=pyplot.figure(figsize=(48,40))
                        
                        #Change label font
                        pyplot.rcParams.update({'font.size': 20})
                        color_map=pyplot.pcolor(np.log10(radar_echo),cmap=pyplot.get_cmap('gray'))#,norm=divnorm)
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


