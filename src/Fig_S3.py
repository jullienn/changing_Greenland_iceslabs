# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 17:24:59 2022

@author: jullienn
"""

def compute_distances(eastings,northings):
    #This function is from plot_2002_2003.py, which was originally taken from MacFerrin et al., 2019
    '''Compute the distance (in m here, not km as written originally) of the traces in the file.'''
    # C = sqrt(A^2  + B^2)
    distances = np.power(np.power((eastings[1:] - eastings[:-1]),2) + np.power((northings[1:] - northings[:-1]),2), 0.5)

    #Calculate the cumsum of the distances
    cumsum_distances=np.nancumsum(distances)
    #Seeting the first value of the cumsum to be zero as it is the origin
    return_cumsum_distances=np.zeros(eastings.shape[0])
    return_cumsum_distances[1:eastings.shape[0]]=cumsum_distances

    return return_cumsum_distances

#Import libraries
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import h5py
import pandas as pd
import pdb
import pickle
from pyproj import Transformer
import matplotlib.gridspec as gridspec
import matplotlib.image as mpimg

general_path='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/'
path_new_method='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2010_2018/'

#Choose the date to illustrate the process
chosen_trace='20100507_01_008_010'

#1. Import original radar trace
path_data_open=general_path+'data/'+chosen_trace[0:4]+'_Greenland_P3/CSARP_qlook/'+chosen_trace[0:11]+'/'

lat_appended=[]
lon_appended=[]
radar_appended=[]

#Define all the individual files of the chosen trace
for indiv_nb in range(int(chosen_trace[12:15]),int(chosen_trace[16:19])+1):
    if (indiv_nb<10):
        indiv_filename=chosen_trace[0:12]+'00'+str(indiv_nb)
    elif ((indiv_nb>=10) and (indiv_nb<100)):
        indiv_filename=chosen_trace[0:12]+'0'+str(indiv_nb)
    else:
        indiv_filename=chosen_trace[0:12]+str(indiv_nb)
    print(indiv_filename)

    if (int(chosen_trace[0:4])>=2014):
        fdata_filename = h5py.File(path_data_open+'Data_'+indiv_filename+'.mat')
        lat_variable=fdata_filename['Latitude'][:,:]
        lon_variable=fdata_filename['Longitude'][:,:]
        time_variable=fdata_filename['Time'][:,:]
        radar_variable=fdata_filename['Data'][:,:]
    else:
        fdata_filename = scipy.io.loadmat(path_data_open+'Data_'+indiv_filename+'.mat')
        lat_variable = fdata_filename['Latitude']
        lon_variable = fdata_filename['Longitude']
        time_variable = fdata_filename['Time']
        radar_variable = fdata_filename['Data']
    
    #Append data
    lat_appended=np.append(lat_appended,lat_variable)
    lon_appended=np.append(lon_appended,lon_variable)
    
    if(indiv_nb==int(chosen_trace[12:15])):
        radar_appended=radar_variable
    else:
        radar_appended=np.append(radar_appended,radar_variable,axis=1)

#Compute the speed (Modified Robin speed):
# self.C / (1.0 + (coefficient*density_kg_m3/1000.0))
v= 299792458 / (1.0 + (0.734*0.873/1000.0))

#calculate depth
depths = v * time_variable / 2.0
#Reset depths to 0
depths=depths-depths[0]
#Identify index where depth > 20 m
ind_lower_20m=np.where(depths<=20)[0]

#Define transformer for coordinates transform from "EPSG:4326" to "EPSG:3413"
transformer = Transformer.from_crs("EPSG:4326", "EPSG:3413", always_xy=True)

#Transform the coordinated from WGS84 to EPSG:3413
#Example from: https://pyproj4.github.io/pyproj/stable/examples.html
points=transformer.transform(np.array(lon_appended),np.array(lat_appended))
lon_3413=points[0]
lat_3413=points[1]

#Calculate distances
distances_with_start_transect=compute_distances(lon_3413,lat_3413)

#2. Show after surface picking - Let's do it manually
#Load surface pick
path_surfacepick=general_path+'data/exported/Surface_Indices_Picklefiles/'
filename_surfacepick=chosen_trace+'_SURFACE.pickle'
f_surfacepick = open(path_surfacepick+filename_surfacepick, "rb")
surfacepick_file = pickle.load(f_surfacepick)
f_surfacepick.close()

#Define index limits for 30m deep
#For 2010-2011, 30m deep index is from 0 to 128 included. From 2012 to 2018, index is from 0:60 included
if (chosen_trace[0:4] in list(['2010','2011'])):
    ind_30m=surfacepick_file+128
    for_raw=428
elif (chosen_trace[0:4] in list(['2012','2013','2014','2017','2018'])):
    ind_30m=surfacepick_file+60
    for_raw=201
else:
    print('Error, year not known')

#set empty slice for radar_30m and raw_data
radar_30m=np.zeros((np.unique(np.abs(surfacepick_file-ind_30m))[0],len(surfacepick_file)))
radar_30m[:]=np.nan

raw_data=np.zeros((for_raw,len(surfacepick_file)))
raw_data[:]=np.nan

#Select radar data
for i in range(0,len(surfacepick_file)):
    radar_30m[:,i]=radar_appended[surfacepick_file[i]:ind_30m[i],i]
    raw_data[:,i]=radar_appended[surfacepick_file[0]-250:surfacepick_file[0]-250+for_raw,i]

#3. After lakes and other exclusions - We need to use the file _SURFACE_SLICE_100M.pickle
path_lakes_excl=general_path+'data/exported/Surface_Slice_100m_Picklefiles/'
filename_lakes_excl=chosen_trace+'_SURFACE_SLICE_100M.pickle'
f_lakes_excl = open(path_lakes_excl+filename_lakes_excl, "rb")
lakes_excl_file = pickle.load(f_lakes_excl)
f_lakes_excl.close()
#Select only the first 20m!
lakes_excl_file_20m=lakes_excl_file[ind_lower_20m,:]

#4. After roll correction
path_roll_corrected=general_path+'data/exported/Roll_Corrected_Picklefiles/'
filename_roll_corrected=chosen_trace+'_ROLL_CORRECTED.pickle'
f_roll_corrected = open(path_roll_corrected+filename_roll_corrected, "rb")
roll_corrected_file = pickle.load(f_roll_corrected)
f_roll_corrected.close()
#Select only the first 20m!
roll_corrected_20m=roll_corrected_file[ind_lower_20m,:]

#5. After surface removal
#GENERATE IT!!!


#6. After depth correction
'''
path_depth_corrected=path_new_method+'ii_out_from_iceslabs_processing_jullien.py/pickles/'
'''
path_depth_corrected=general_path+'/inspection_april2022/indiv_prob/'
filename_depth_corrected=chosen_trace+'_Depth_Corrected_surf_removal.pickle'
f_depth_corrected = open(path_depth_corrected+filename_depth_corrected, "rb")
depth_corrected_file = pickle.load(f_depth_corrected)
f_depth_corrected.close()
#Select only the first 20m!
depth_corrected_20m=depth_corrected_file[ind_lower_20m,:]

#7. After ice slabs likelihood computation
'''
path_likelihood=path_new_method+'/iii_out_from_probabilistic_iceslabs.py/before_DF_exclusion/pickles/'
'''
path_likelihood=general_path+'/inspection_april2022/probability_iceslabs/before_DF_appliance/pickles/'
filename_likelihood=chosen_trace+'_probability_iceslabs_presence.pickle'
f_likelihood = open(path_likelihood+filename_likelihood, "rb")
likelihood_file = pickle.load(f_likelihood)
f_likelihood.close()

#8. After dry firn exclusion, final ice slab product, with average columnal ice
#   likelihood displayed as colorshade at the surface??
'''
path_likelihood_after_DF=path_new_method+'/iii_out_from_probabilistic_iceslabs.py/after_DF_exclusion/pickles/'
'''
path_likelihood_after_DF=general_path+'/inspection_april2022/probability_iceslabs/after_DF_appliance/pickles/'
filename_likelihood_after_DF=chosen_trace+'_probability_iceslabs_presence_after_DF.pickle'
f_likelihood_after_DF = open(path_likelihood_after_DF+filename_likelihood_after_DF, "rb")
likelihood_file_after_DF = pickle.load(f_likelihood_after_DF)
f_likelihood_after_DF.close()

print('ATTENTION, b IS RESCALED AS WENT THROUGH THE RESCALING PROCEDURE WHILE EXPORTING')
print('=> EITHER EXTRACT THE NON RESCALED SLICE OR RESCALE ALL THE REMAINING!')
print('NEED TO SHOW INFLUENCE, THUS CONSIDER FIXING COLORBAR TO ALL SUBPLOTS')
#Create the figure
fig = plt.figure(figsize=(24,10))
gs = gridspec.GridSpec(24, 10)
gs.update(wspace=0.001)
ax1 = plt.subplot(gs[0:3, 0:10])
ax2 = plt.subplot(gs[3:6, 0:10])
ax3 = plt.subplot(gs[6:9, 0:10])
ax4 = plt.subplot(gs[9:12, 0:10])
ax5 = plt.subplot(gs[12:15, 0:10])
ax6 = plt.subplot(gs[15:18, 0:10])
ax7 = plt.subplot(gs[18:21, 0:10])
ax8 = plt.subplot(gs[21:24, 0:10])

#ax1 is raw radar
#ax1 = radar_appended
ax1.pcolor(distances_with_start_transect, depths[0:raw_data.shape[0]], np.log10(raw_data),cmap=plt.get_cmap('gray'),zorder=-2)#,norm=divnorm)
ax1.invert_yaxis() #Invert the y axis = avoid using flipud.    
ax1.set_ylim(100,0)
#ax2.setp(ax2.get_xticklabels(), visible=False)
#ax2.set_yticklabels(['0','10',''])
ax1.text(0.01, 0.75,'a',ha='center', va='center', transform=ax1.transAxes,fontsize=15,zorder=10,weight='bold',color='white')#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot

#ax2 is after surface picking
ax2.pcolor(distances_with_start_transect, depths[ind_lower_20m], np.log10(radar_30m[ind_lower_20m,:]),cmap=plt.get_cmap('gray'),zorder=-2)#,norm=divnorm)
ax2.invert_yaxis() #Invert the y axis = avoid using flipud.    
ax2.set_ylim(20,0)
#ax2.setp(ax2.get_xticklabels(), visible=False)
#ax2.set_yticklabels(['0','10',''])
ax2.text(0.01, 0.75,'b',ha='center', va='center', transform=ax2.transAxes,fontsize=15,zorder=10,weight='bold',color='white')#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot

#ax3 is after appliance of lakes and other exclusions
ax3.pcolor(distances_with_start_transect, depths[ind_lower_20m], np.log10(lakes_excl_file_20m),cmap=plt.get_cmap('gray'),zorder=-2)#,norm=divnorm)
ax3.invert_yaxis() #Invert the y axis = avoid using flipud.    
ax3.set_ylim(20,0)
#ax3.setp(ax3.get_xticklabels(), visible=False)
#ax3.set_yticklabels(['0','10',''])
ax3.text(0.01, 0.75,'c',ha='center', va='center', transform=ax3.transAxes,fontsize=15,zorder=10,weight='bold',color='white')#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot

#ax4 is after appliance of roll correction
ax4.pcolor(distances_with_start_transect, depths[ind_lower_20m], roll_corrected_20m,cmap=plt.get_cmap('gray'),zorder=-2)#,norm=divnorm)
ax4.invert_yaxis() #Invert the y axis = avoid using flipud.    
ax4.set_ylim(20,0)
#ax4.setp(ax4.get_xticklabels(), visible=False)
#ax4.set_yticklabels(['0','10',''])
ax4.text(0.01, 0.75,'d',ha='center', va='center', transform=ax4.transAxes,fontsize=15,zorder=10,weight='bold',color='white')#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot

#ax5 is after surface removal

#ax6 is after depth correction
ax6.pcolor(distances_with_start_transect, depths[ind_lower_20m], depth_corrected_20m,cmap=plt.get_cmap('gray'),zorder=-2)#,norm=divnorm)
ax6.invert_yaxis() #Invert the y axis = avoid using flipud.    
ax6.set_ylim(20,0)
#ax6.setp(ax6.get_xticklabels(), visible=False)
#ax6.set_yticklabels(['0','10',''])
ax6.text(0.01, 0.75,'f',ha='center', va='center', transform=ax6.transAxes,fontsize=15,zorder=10,weight='bold',color='white')#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot

#ax7 is after ice layer likelihood computation
ax7.pcolor(distances_with_start_transect, depths[ind_lower_20m], likelihood_file,cmap=plt.get_cmap('gray_r'),zorder=-2)#,norm=divnorm)
ax7.invert_yaxis() #Invert the y axis = avoid using flipud.    
ax7.set_ylim(20,0)
#ax7.setp(ax7.get_xticklabels(), visible=False)
#ax7.set_yticklabels(['0','10',''])
ax7.text(0.01, 0.75,'g',ha='center', va='center', transform=ax7.transAxes,fontsize=15,zorder=10,weight='bold',color='white')#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot

#ax8 is after appliance of dry firn exclusions appliance - final ice slabs product
ax8.pcolor(distances_with_start_transect, depths[ind_lower_20m], likelihood_file_after_DF,cmap=plt.get_cmap('gray_r'),zorder=-2)#,norm=divnorm)
ax8.invert_yaxis() #Invert the y axis = avoid using flipud.    
ax8.set_ylim(20,0)
#ax8.setp(ax8.get_xticklabels(), visible=False)
#ax8.set_yticklabels(['0','10',''])
ax8.text(0.01, 0.75,'h',ha='center', va='center', transform=ax8.transAxes,fontsize=15,zorder=10,weight='bold',color='white')#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot

plt.show()