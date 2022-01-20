# -*- coding: utf-8 -*-
"""
Created on Thu Jan 20 10:41:00 2022

@author: JullienN
"""
    
import pickle
import scipy.io
import numpy as np
import pdb
import h5py
import matplotlib.pyplot as plt
import pandas as pd
import geopandas as gpd
import matplotlib.gridspec as gridspec
from pyproj import Transformer
import seaborn as sns
sns.set_theme(style="whitegrid")
import matplotlib.image as mpimg


### -------------------------- Load shapefiles --------------------------- ###
#Load Rignot et al., 2016 Greenland drainage bassins
path_rignotetal2016_GrIS_drainage_bassins='C:/Users/jullienn/switchdrive/Private/research/backup_Aglaja/working_environment/greenland_topo_data/GRE_Basins_IMBIE2_v1.3/'
GrIS_drainage_bassins=gpd.read_file(path_rignotetal2016_GrIS_drainage_bassins+'GRE_Basins_IMBIE2_v1.3.shp',rows=slice(51,57,1)) #the regions are the last rows of the shapefile

#Extract indiv regions and create related indiv shapefiles
SW_rignotetal=GrIS_drainage_bassins[GrIS_drainage_bassins.SUBREGION1=='SW']
CW_rignotetal=GrIS_drainage_bassins[GrIS_drainage_bassins.SUBREGION1=='CW']
### -------------------------- Load shapefiles --------------------------- ###
panel_a={2010:'empty',
      2011:['Data_20110516_01_009.mat','Data_20110516_01_010.mat'],
      2012:'empty',
      2013:['Data_20130402_01_008.mat'],
      2014:'empty',
      2017:['Data_20170412_01_075.mat','Data_20170412_01_076.mat'],
      2018:'empty'}

panel_b={2010:['Data_20100517_02_001.mat','Data_20100517_02_002.mat'],
      2011:['Data_20110502_01_171.mat'],
      2012:['Data_20120516_01_002.mat'],
      2013:['Data_20130419_01_004.mat','Data_20130419_01_005.mat'],
      2014:['Data_20140507_03_007.mat','Data_20140507_03_008.mat'], #test with 20140514_02_087_089 and 20140515_02_173_175 also
      2017:['Data_20170417_01_171.mat','Data_20170417_01_172.mat','Data_20170417_01_173.mat','Data_20170417_01_174.mat'],
      2018:'empty'}

panel_c={2010:['Data_20100508_01_114.mat','Data_20100508_01_115.mat'],
      2011:['Data_20110419_01_008.mat','Data_20110419_01_009.mat','Data_20110419_01_010.mat'],
      2012:['Data_20120418_01_129.mat','Data_20120418_01_130.mat','Data_20120418_01_131.mat'],
      2013:['Data_20130405_01_165.mat','Data_20130405_01_166.mat','Data_20130405_01_167.mat'],
      2014:['Data_20140424_01_002.mat','Data_20140424_01_003.mat','Data_20140424_01_004.mat'],
      2017:['Data_20170422_01_168.mat','Data_20170422_01_169.mat','Data_20170422_01_170.mat','Data_20170422_01_171.mat'],
      2018:['Data_20180427_01_170.mat','Data_20180427_01_171.mat','Data_20180427_01_172.mat']}

#This one is collocated with FS1, 2, 3.
panel_d={2010:'empty',
      2011:'empty',
      2012:['Data_20120423_01_137.mat','Data_20120423_01_138.mat'],
      2013:['Data_20130409_01_010.mat','Data_20130409_01_011.mat','Data_20130409_01_012.mat'],
      2014:'empty',
      2017:'empty',
      2018:['Data_20180421_01_004.mat','Data_20180421_01_005.mat','Data_20180421_01_006.mat','Data_20180421_01_007.mat']}

panel_e={2010:['Data_20100513_01_001.mat','Data_20100513_01_002.mat'],
      2011:['Data_20110411_01_116.mat','Data_20110411_01_117.mat','Data_20110411_01_118.mat'],
      2012:['Data_20120428_01_125.mat','Data_20120428_01_126.mat'],
      2013:'empty',
      2014:['Data_20140408_11_024.mat','Data_20140408_11_025.mat','Data_20140408_11_026.mat'],
      2017:['Data_20170508_02_165.mat','Data_20170508_02_166.mat','Data_20170508_02_167.mat','Data_20170508_02_168.mat','Data_20170508_02_169.mat','Data_20170508_02_170.mat','Data_20170508_02_171.mat'],
      2018:'empty'}

#Define the panel to study
investigation_year=panel_b

fig = plt.figure()
gs = gridspec.GridSpec(30, 10)
gs.update(wspace=0.1)
#gs.update(wspace=0.001)
ax1 = plt.subplot(gs[0:4, 0:10])
ax2 = plt.subplot(gs[4:8, 0:10])
ax3 = plt.subplot(gs[8:12, 0:10])
ax4 = plt.subplot(gs[12:16, 0:10])
ax5 = plt.subplot(gs[16:20, 0:10])
ax6 = plt.subplot(gs[20:24, 0:10])
ax7 = plt.subplot(gs[24:28, 0:10])

#Fig_prob
fig_prob = plt.figure()
gs = gridspec.GridSpec(30, 10)
gs.update(wspace=0.1)
#gs.update(wspace=0.001)
ax1_prob = plt.subplot(gs[0:4, 0:10])
ax2_prob = plt.subplot(gs[4:8, 0:10])
ax3_prob = plt.subplot(gs[8:12, 0:10])
ax4_prob = plt.subplot(gs[12:16, 0:10])
ax5_prob = plt.subplot(gs[16:20, 0:10])
ax6_prob = plt.subplot(gs[20:24, 0:10])
ax7_prob = plt.subplot(gs[24:28, 0:10])

#Define paths
path_data_MacFerrin='C:/Users/jullienn/switchdrive/Private/research/RT1/MacFerrin_FigShare/IceBridge_Accumulation_Radar/IceBridge_Accumulation_Radar/Images_PostProcessed_Greyscale/'
path_data_Jullien='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2010_2018/i_out_from_IceBridgeGPR_Manager_v2.py/pickles_and_images/'
path_data_Jullien_probability='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2010_2018/iii_out_from_probabilistic_iceslabs.py/images/'
path_data='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/'

#Compute the speed (Modified Robin speed):
# self.C / (1.0 + (coefficient*density_kg_m3/1000.0))
v= 299792458

dataframe={}

#lon_for_min_max
lon_for_min_max=[]

for single_year in investigation_year.keys():
    print(single_year)
        
    #If no data, continue
    if (investigation_year[single_year]=='empty'):
        print('No data for year '+str(single_year)+', continue')
        continue
    
    ###1. Load the images:
    start_date_track=investigation_year[single_year][0]
    end_date_track=investigation_year[single_year][-1]
    date_track=start_date_track[5:20]+'_'+end_date_track[17:20]
    
    #pat_to_use=path_data_MacFerrin+date_track+'_1'
    pat_to_use=path_data_Jullien+date_track+'_X'
    
    filename_image=pat_to_use+'DEPTHCORRECT_AFTER.png'
    #Load image
    img = mpimg.imread(filename_image)
    
    #Load probability
    path_prob=path_data_Jullien_probability+date_track+'_probability_iceslabs_presence.png'
    #Load probability image
    img_prob = mpimg.imread(path_prob)
    
    ###2. Load the latitude and longitude
    lat_appended=[]
    lon_appended=[]
    
    for indiv_file_load in investigation_year[single_year]:
        print(indiv_file_load)
        
        #Create the path
        path_raw_data=path_data+str(single_year)+'_Greenland_P3/CSARP_qlook/'+indiv_file_load[5:16]+'/'

        # If 20180421_01_007, different path
        if (indiv_file_load=='Data_20180421_01_007.mat'):
            path_raw_data='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/2018_Greenland_P3/'
        
        #Load data
        if (single_year>=2014):
            
            fdata_filename = h5py.File(path_raw_data+indiv_file_load)
            lat_filename=fdata_filename['Latitude'][:,:]
            lon_filename=fdata_filename['Longitude'][:,:]
            time_filename=fdata_filename['Time'][:,:]
            
        else:
            fdata_filename = scipy.io.loadmat(path_raw_data+indiv_file_load)
            lat_filename = fdata_filename['Latitude']
            lon_filename = fdata_filename['Longitude']
            time_filename = fdata_filename['Time']
            
        #Append data
        lat_appended=np.append(lat_appended,lat_filename)
        lon_appended=np.append(lon_appended,lon_filename)
        
    #Check whether the data are acquired ascending or descending elevation wise.
    #I choose the ascending format. For the data that are descending, reverse them
    #To have ascending data, the longitude should increase
    #(-48 = low elevations, -46 = higher elevations)
    
    if (np.sum(np.diff(lon_appended))<0):
        #It is going toward lower elevations, thus flip left-right
        #(or up-down) all the data!
        
        lat_appended=np.flipud(lat_appended)
        lon_appended=np.flipud(lon_appended)
        img=np.fliplr(img)
        img_prob=np.fliplr(img_prob)
    
    #Calculate the depth from the time
    #########################################################################
    # From plot_2002_2003.py - BEGIN
    #########################################################################
    depth_check = v * time_filename / 2.0
    
    #If 2014, transpose the vector
    if (str(single_year)>='2014'):
        depth_check=np.transpose(depth_check)
    
    #Reset times to zero! This is from IceBridgeGPR_Manager_v2.py
    if (depth_check[10]<0):
        #depth_check[10] so that I am sure that the whole vector is negative and
        #not the first as can be for some date were the proccessing is working
        depth_check=depth_check+abs(depth_check[0])
        depth = depth_check
    else:
        depth = depth_check
    
    if (str(single_year) in list(['2011','2012','2014','2017','2018'])):
        if (depth_check[10]>1):
            #depth_check[10] so that I am sure that the whole vector is largely positive and
            #not the first as can be for some date were the proccessing is working
            depth_check=depth_check-abs(depth_check[0])
            depth = depth_check
        
    #Store reunited lat/lon, slice output and mask in a dictionnary:
    dataframe[str(single_year)]={'lat_appended':lat_appended,
                                 'lon_appended':lon_appended,
                                 'depth':depth,
                                 'img':img,
                                 'img_prob':img_prob}
    
    #Append for lon min/max computation
    lon_for_min_max=np.append(lon_for_min_max,lon_appended)
    
#Define min and max
x_start=np.min(lon_for_min_max)
x_end=np.max(lon_for_min_max)

for single_year in investigation_year.keys():
    
    if (investigation_year[single_year]=='empty'):
        continue
    print(single_year)
    
    if (single_year==2010):
        ax_plot=ax1
        ax_prob=ax1_prob
        color_toplot="#4dac26"
        ax_plot.set_xticklabels([])
    elif (single_year==2011):
        ax_plot=ax2
        ax_prob=ax2_prob
        color_toplot="#0571b0"
        ax_plot.set_xticklabels([])
    elif (single_year==2012):
        ax_plot=ax3
        ax_prob=ax3_prob
        color_toplot="#e31a1c"
        ax_plot.set_xticklabels([])
    elif (single_year==2013):
        ax_plot=ax4
        ax_prob=ax4_prob
        color_toplot="#e31a1c"
        ax_plot.set_xticklabels([])
    elif (single_year==2014):
        ax_plot=ax5
        ax_prob=ax5_prob
        color_toplot="#2171b5"
        ax_plot.set_xticklabels([])
    elif (single_year==2017):
        ax_plot=ax6
        ax_prob=ax6_prob
        color_toplot="#2171b5"
        ax_plot.set_xticklabels([])
    elif (single_year==2018):
        ax_plot=ax7
        ax_prob=ax7_prob
        color_toplot="#e31a1c"
        ax_plot.set_xlabel('Longitude [°]')
        #Activate ticks xlabel
        ax_plot.xaxis.tick_bottom()
    else:
        print('year not know')
    
    #Load data
    X=dataframe[str(single_year)]['lon_appended']
    Y=np.arange(0,30,30/dataframe[str(single_year)]['img'].shape[0])
    Y_prob=np.arange(0,20,20/dataframe[str(single_year)]['img_prob'].shape[0])
    C=dataframe[str(single_year)]['img']
    C_prob=dataframe[str(single_year)]['img_prob']
    
    #plot data
    cb=ax_plot.pcolor(X, Y, C,cmap=plt.get_cmap('gray'),zorder=-1)#,norm=divnorm)
    ax_plot.invert_yaxis() #Invert the y axis = avoid using flipud.  
    
    ax_prob.pcolor(X, Y_prob, C_prob,cmap=plt.get_cmap('gray'),zorder=-1)#,norm=divnorm)
    ax_prob.invert_yaxis() #Invert the y axis = avoid using flipud.  
    
    #Set xlim
    ax_plot.set_xlim(x_start,x_end)
    ax_prob.set_xlim(x_start,x_end)
    
    #Set ylim
    ax_plot.set_ylim(20,0)

    #Activate ticks ylabel
    ax_plot.yaxis.tick_left()
    ax_prob.yaxis.tick_left()
    
    '''
    # instantiate a second axes that shares the same x-axis. This is from https://stackoverflow.com/questions/13369888/matplotlib-y-axis-label-on-right-side
    ax2 = ax_plot.twinx()
    ax2.set_yticklabels([])
    ax2.set_yticks([])
    ax2_prob = ax_prob.twinx()
    ax2_prob.set_yticklabels([])
    ax2_prob.set_yticks([])    
    
    #Add year on radargram
    ax2.set_ylabel(str(single_year), color=color_toplot, weight='bold')#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot
    ax2_prob.set_ylabel(str(single_year), color=color_toplot, weight='bold')#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot
    '''

ax1.set_ylabel('Depth [m]')
ax1_prob.set_ylabel('Depth [m]')

figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()

plt.show()