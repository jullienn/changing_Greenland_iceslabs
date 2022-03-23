# -*- coding: utf-8 -*-
"""
Created on Sun Nov 28 12:56:32 2021

@author: jullienn

Code adapted from lenses_thickening_visualisation.py
"""

def plot_thickness(dictionnary_case_study,dataframe,df_2010_2018_elevation,axt,my_pal):
    #This function is adapted from plot_thickness_evolution from fig_2_paper_iceslabs.py
    
    ###########################################################################
    ###                    Display thickness evolution                      ###
    ###########################################################################
    #Define empty dictionnary for elevation slice definition
    df_for_elev=pd.DataFrame(columns=list(df_2010_2018_elevation.keys()))
    
    #Loop over the years
    for year in dictionnary_case_study.keys():
        if (dictionnary_case_study[year] == 'empty'):
            continue  
        #Select data for the trace
        df_for_elev_temp=df_2010_2018_elevation[df_2010_2018_elevation['Track_name']==dictionnary_case_study[year][0][5:20]+'_'+dictionnary_case_study[year][-1][17:20]]
        #Append data to each other
        df_for_elev=df_for_elev.append(df_for_elev_temp)

    #Create an empty df_sampling
    df_sampling=pd.DataFrame(columns=['Track_name','time_period','low_bound', 'high_bound', 'bound_nb', 'mean', 'median', 'q025', 'q075','stddev'])

    #Elev divide every 1m.
    elev_bin_desired=1
    elev_divide=np.arange(np.floor(np.min(df_for_elev['elevation'])),np.floor(np.max(df_for_elev['elevation']))+1+elev_bin_desired,elev_bin_desired)
    
    #Define window size for smoothing
    winsize=10

    #Loop over the different time periods (2010, 2011-2012, 2013-2014, 2017-2018)
    for indiv_year in range(2012,2019):
        print(indiv_year)
        #Get data for that specific time period
        df_trace_year=df_for_elev[df_for_elev['year']==indiv_year]
        
        if (len(df_trace_year)==0):
            #No data in this time period, continue
            print('   No data for this year, continue')
            continue
        else:            
            #Sort df_trace_year from low to high elevations
            df_trace_year_sorted=df_trace_year.sort_values(by=['elevation'])
            
            #Set bound_nb to 0
            bound_nb=0
            #Loop over the lon divide
            for i in range(1,len(elev_divide)):
                
                #Identify low and higher end of the slice
                low_bound=elev_divide[i-1]
                high_bound=elev_divide[i]
        
                #Select all the data belonging to this elev slice
                ind_slice=np.logical_and(np.array(df_trace_year_sorted['elevation']>=low_bound),np.array(df_trace_year_sorted['elevation']<high_bound))
                df_select=df_trace_year_sorted[ind_slice]
                
                if (len(df_select)==0):
                    continue
                
                #Fill in dictionnary
                df_temp=pd.DataFrame(columns=['Track_name','time_period','low_bound', 'high_bound', 'bound_nb', 'mean', 'median', 'q025', 'q075','stddev'])
                df_temp['Track_name']=df_select['Track_name'].unique()
                df_temp['time_period']=indiv_year
                df_temp['low_bound']=low_bound
                df_temp['high_bound']=high_bound
                df_temp['bound_nb']=str(bound_nb)
                df_temp['mean']=np.nanmean(df_select['20m_ice_content_m'])
                df_temp['median']=np.nanmedian(df_select['20m_ice_content_m'])
                df_temp['q025']=np.nanquantile(df_select['20m_ice_content_m'],0.25)
                df_temp['q075']=np.nanquantile(df_select['20m_ice_content_m'],0.75)
                df_temp['stddev']=np.nanstd(df_select['20m_ice_content_m'])
                
                #Append dictionnary
                df_sampling=df_sampling.append(df_temp)
                
                #Update bound_nb
                bound_nb=bound_nb+1
    
    for time_period in range(2012,2019):
    
        if (len(df_sampling[df_sampling['time_period']==time_period])==0):
            #Empty time period, continue
            continue
        else:
            df_plot=df_sampling[df_sampling['time_period']==time_period]
            
            #Rolling window, size = 10
            df_plot['rolling_10_median']=df_plot.rolling(winsize, win_type=None,center=True).quantile(quantile=0.5)['median'] 
            df_plot['rolling_10_q025']=df_plot.rolling(winsize, win_type=None,center=True).quantile(quantile=0.5)['q025'] 
            df_plot['rolling_10_q075']=df_plot.rolling(winsize, win_type=None,center=True).quantile(quantile=0.5)['q075'] 
            
            # Plot the median
            axt.plot(df_plot["low_bound"],df_plot["rolling_10_median"],color=my_pal[time_period])
            #Display IQR
            axt.fill_between(df_plot['low_bound'], df_plot['rolling_10_q025'], df_plot['rolling_10_q075'], alpha=0.3,color=my_pal[time_period])
            
            #Display the median where outside of average window range
            #Create array_fill_start and _end for filling at the start and at the end
            array_fill_start=np.zeros(6,)
            array_fill_start[:]=np.nan
            array_fill_start[0:5]=np.asarray(df_plot["median"].iloc[0:int(winsize/2)])
            array_fill_start[-1]=np.asarray((df_plot['rolling_10_median'].iloc[int(winsize/2)]))
            
            array_fill_end=np.zeros(5,)
            array_fill_end[:]=np.nan
            array_fill_end[0]=np.asarray((df_plot['rolling_10_median'].iloc[int(len(df_plot)-winsize/2)]))
            array_fill_end[1:5]=np.asarray(df_plot["median"].iloc[int(len(df_plot)-winsize/2+1):len(df_plot)])
            
            #Display
            axt.plot(df_plot["low_bound"].iloc[0:int(winsize/2)+1],array_fill_start,alpha=0.5,color=my_pal[time_period])
            axt.plot(df_plot["low_bound"].iloc[int(len(df_plot)-winsize/2):len(df_plot)],array_fill_end,alpha=0.5,color=my_pal[time_period])
    
    ###########################################################################
    ###                    Display thickness evolution                      ###
    ###########################################################################
    
    
    ###########################################################################
    ###                           Display radargrams                        ###
    ###########################################################################
    #Loop over the years
    for year in dictionnary_case_study.keys():
        #Pick up the corresponding datetrack
        date_track=dictionnary_case_study[year][0][5:20]+'_'+dictionnary_case_study[year][-1][17:20]
        
        #Reset depths to 0
        dataframe[str(year)]['depth']=dataframe[str(year)]['depth']-dataframe[str(year)]['depth'][0]
        
        #Select radar slice
        depth_corrected_file=dataframe[str(year)]['radar']
        
        #Identify index where time < 20 m
        ind_lower_20m=np.where(dataframe[str(year)]['depth']<20)[0]
        depth_corrected_20m=depth_corrected_file[ind_lower_20m,:]
        
        #Identify axis for plotting
        if (year==2010):
            ax_plotting=ax1r
            ax1r.set_xlabel('Longitude [°]')
            label_for_map=u'\u03B1'
            #Adapt xticklabels
            ax_plotting.set_yticklabels(['0', '10', ''])
            ax_plotting.set_xticklabels([])
            ax_plotting.set_yticks([0,10,])
            #Define pannel label
            casestudy_nb='a'
        elif (year==2011):
            ax_plotting=ax2r
            ax2r.set_xlabel('Latitude [°]')
            ax2r.set_ylabel('Depth [m]')
            label_for_map=u'\u03B2'
            #Activate ticks xlabel
            ax_plotting.xaxis.tick_bottom()
            #Define pannel label
            casestudy_nb='b'
        elif (year==2012):
            ax_plotting=ax3r
            ax_plotting.axvline(x=-47.11,zorder=1,linestyle='--',color='k')
            ax_plotting.axvline(x=-47.02,zorder=1,linestyle='--',color='k')
            label_for_map=u'\u03B3'
            #Adapt xticklabels
            ax_plotting.set_yticklabels(['0', '10', ''])
            ax_plotting.set_xticklabels([])
            ax_plotting.set_yticks([0,10,])
            #Define pannel label
            casestudy_nb='c'
        elif (year==2013):
            ax_plotting=ax4r
            ax4r.set_xlabel('Depth [m]')
            ax_plotting.axvline(x=-47.11,zorder=1,linestyle='--',color='k')
            ax_plotting.axvline(x=-47.02,zorder=1,linestyle='--',color='k')
            label_for_map=u'\u03B3'
            #Adapt xticklabels
            ax_plotting.set_yticklabels(['0', '10', ''])
            ax_plotting.set_xticklabels([])
            ax_plotting.set_yticks([0,10,])
            #Define pannel label
            casestudy_nb='d'
        elif (year==2014):
            ax_plotting=ax5r
            label_for_map=u'\u03B4'
            #Adapt xticklabels, from https://stackoverflow.com/questions/43673884/change-x-axis-ticks-to-custom-strings
            ax_plotting.set_yticklabels(['0', '10', ''])
            ax_plotting.set_xticklabels([])
            ax_plotting.set_yticks([0,10,])
            #Define pannel label
            casestudy_nb='e'
        elif (year==2017):
            ax_plotting=ax6r
            label_for_map=u'\u03B4'
            #Adapt xticklabels
            ax_plotting.set_yticklabels(['0', '10', ''])
            ax_plotting.set_xticklabels([])
            ax_plotting.set_yticks([0,10,])
            #Define pannel label
            casestudy_nb='f'
        elif (year==2018):
            ax_plotting=ax7r
            ax_plotting.axvline(x=-47.11,zorder=1,linestyle='--',color='k')
            ax_plotting.axvline(x=-47.02,zorder=1,linestyle='--',color='k')
            label_for_map=u'\u03B3'
            #Activate ticks xlabel
            ax_plotting.xaxis.tick_bottom()
            #Define pannel label
            casestudy_nb='g'
        else:
            print('Year not existing')
        
        #Select x vector
        if (year==2011):
            X=dataframe[str(year)]['lat_appended']
        else:
            X=dataframe[str(year)]['lon_appended']
        
        #Select y vector and radargram
        Y=np.arange(0,100,100/dataframe[str(year)]['radar'].shape[0])
        C=dataframe[str(year)]['radar']
                
        cb=ax_plotting.pcolor(X, Y, C,cmap=plt.get_cmap('gray'),zorder=-1)#,norm=divnorm)
        ax_plotting.invert_yaxis() #Invert the y axis = avoid using flipud.    
        ax_plotting.set_ylim(20,0)
        
        ###########################################################################
        ###                           Display radargrams                        ###
        ########################################################################### 
        
        ###########################################################################
        ###                       Display data localisation                     ###
        ###########################################################################
        
        #Create lat/lon vectors for display
        lon_plot=dataframe[str(year)]['lon_appended'][dataframe[str(year)]['mask']]
        lat_plot=dataframe[str(year)]['lat_appended'][dataframe[str(year)]['mask']]
        
        lon3413_plot=dataframe[str(year)]['lon_3413'][dataframe[str(year)]['mask']]
        lat3413_plot=dataframe[str(year)]['lat_3413'][dataframe[str(year)]['mask']]
        
        if (year==2011):
            ax_plotting.set_xlim(66.8707,67.2)
            ind_map=np.logical_and(lat_plot>=66.8707,lat_plot<=67.2)
            #Display KAN_U
            ax_plotting.scatter(67.000425,1,s=10,c='r')
        else:
            ax_plotting.set_xlim(-47.5,-46.66)
            ind_map=np.logical_and(lon_plot>=-47.5,lon_plot<=-46.66)
            #Display KAN_U
            ax_plotting.scatter(-47.030473,1,s=10,c='r')
        
        #display loc on map
        ax8map.scatter(lon3413_plot[ind_map],lat3413_plot[ind_map],c='k',s=0.5)

        #Add year on radargram
        ax_plotting.text(0.96, 0.90,str(year)+', '+label_for_map, color=my_pal[year],zorder=10, ha='center', va='center', transform=ax_plotting.transAxes, weight='bold')#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot
        
        #Add pannel label
        ax_plotting.text(0.01, 0.875,casestudy_nb,ha='center', va='center', transform=ax_plotting.transAxes,weight='bold',fontsize=20,zorder=10)#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot
        
        #Activate ticks ylabel
        ax_plotting.yaxis.tick_left()

        ###########################################################################
        ###                       Display data localisation                     ###
        ###########################################################################
    
    return np.min(df_for_elev['elevation']),np.max(df_for_elev['elevation'])


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
import cartopy.crs as ccrs
from pyproj import Transformer

#Define palette as a function of time
#This is from https://www.python-graph-gallery.com/33-control-colors-of-boxplot-seaborn
my_pal = {2010: "k", 2011: "k", 2012: "#1a9850", 2013: "#542788", 2014: "#2166ac", 2017:"#bf812d",2018:"#b2182b"}

### -------------------------- Load shapefiles --------------------------- ###
#Load Rignot et al., 2016 Greenland drainage bassins
path_rignotetal2016_GrIS_drainage_bassins='C:/Users/jullienn/switchdrive/Private/research/backup_Aglaja/working_environment/greenland_topo_data/GRE_Basins_IMBIE2_v1.3/'
GrIS_drainage_bassins=gpd.read_file(path_rignotetal2016_GrIS_drainage_bassins+'GRE_Basins_IMBIE2_v1.3_EPSG_3413.shp',rows=slice(51,57,1)) #the regions are the last rows of the shapefile

#Extract indiv regions and create related indiv shapefiles
SW_rignotetal=GrIS_drainage_bassins[GrIS_drainage_bassins.SUBREGION1=='SW']
CW_rignotetal=GrIS_drainage_bassins[GrIS_drainage_bassins.SUBREGION1=='CW']
### -------------------------- Load shapefiles --------------------------- ###

#Best continuity between overlap through 2012-2018
investigation_year={2010:['Data_20100515_01_007.mat','Data_20100515_01_008.mat','Data_20100515_01_009.mat'],
                    2011:['Data_20110408_01_087.mat','Data_20110408_01_088.mat','Data_20110408_01_089.mat',
                          'Data_20110408_01_090.mat','Data_20110408_01_091.mat','Data_20110408_01_092.mat',
                          'Data_20110408_01_093.mat','Data_20110408_01_094.mat','Data_20110408_01_095.mat',
                          'Data_20110408_01_096.mat','Data_20110408_01_097.mat','Data_20110408_01_098.mat',
                          'Data_20110408_01_099.mat','Data_20110408_01_100.mat','Data_20110408_01_101.mat',
                          'Data_20110408_01_102.mat','Data_20110408_01_103.mat'],
                    2012:['Data_20120423_01_137.mat','Data_20120423_01_138.mat'],
                    2013:['Data_20130409_01_010.mat','Data_20130409_01_011.mat','Data_20130409_01_012.mat'],
                    2014:['Data_20140416_05_035.mat','Data_20140416_05_036.mat','Data_20140416_05_037.mat'],
                    2017:['Data_20170502_01_171.mat','Data_20170502_01_172.mat','Data_20170502_01_173.mat'],
                    2018:['Data_20180421_01_004.mat','Data_20180421_01_005.mat','Data_20180421_01_006.mat','Data_20180421_01_007.mat']}
'''
#Best absolute closest point for each year, without any overlapping consideration
20180426_01_004_006
20170506_01_010_012
20170508_02_011_013
20140408_04_001_003
20130408_01_009_021
20120423_01_137_138
20110408_01_087_103
20100515_01_007_009

#Best trade off between continuity and overlap
20180426_01_004_006
20170508_02_011_013
20140408_04_001_003
20130409_01_010_012
20120423_01_137_138
2011 does not matter where
2010 does not matter where

#grouping overlapping
#1
20120423_01_137_138
20130409_01_010_012
20180421_01_004_007

#2
20140416_05_035_037
20170502_01_171_173

#3 closer to KAN_u compared to #2
20140408_04_001_003
20170508_02_011_013

#4
20180426_01_004_006
'''
'''
#Closer to KAN_U
investigation_year={2010:'empty',
                    2011:'empty',
                    2012:'empty',
                    2013:['20130408_01_009_021'],
                    2014:['20140408_04_001_003'],
                    2017:['20170508_02_011_013'],
                    2018:['20180426_01_004_006']}
'''
#Define paths
path_data='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/'
path_depth_corrected='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2010_2018/ii_out_from_iceslabs_processing_jullien.py/pickles/'
path_mask='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2010_2018/i_out_from_IceBridgeGPR_Manager_v2.py/pickles_and_images/Boolean_Array_Picklefiles/'
path_probabilistic='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2010_2018/iii_out_from_probabilistic_iceslabs.py/pickles/'

#Define transformer for coordinates transform from "EPSG:4326" to "EPSG:3413"
transformer = Transformer.from_crs("EPSG:4326", "EPSG:3413", always_xy=True)

#Compute the speed (Modified Robin speed):
# self.C / (1.0 + (coefficient*density_kg_m3/1000.0))
v= 299792458 / (1.0 + (0.734*0.873/1000.0))

dataframe={}

for single_year in investigation_year.keys():
    print(single_year)
        
    #If no data, continue
    if (investigation_year[single_year]=='empty'):
        print('No data for year '+str(single_year)+', continue')
        continue
    
    ###1. Load the depth_corrected and probabilistic files:
    start_date_track=investigation_year[single_year][0]
    end_date_track=investigation_year[single_year][-1]
    date_track=start_date_track[5:20]+'_'+end_date_track[17:20]
    
    '''
    if (single_year==2011):
        date_track='20110408_01_087_103'
    if (single_year==2018):
        date_track='20180421_01_004_007'
    #pdb.set_trace()
    
    '''
    
    '''
    #Define filename of depth corrected and probabilistic data
    if (date_track not in list_trace_failed):
        filename_depth_corrected=date_track+'_Depth_Corrected_surf_removal.pickle'
    else:
        filename_depth_corrected=date_track+'_Depth_Corrected_surf_removal_rescaled.pickle'
    '''
    
    filename_depth_corrected=date_track+'_Depth_Corrected_surf_removal.pickle'
    filename_mask=date_track+'_mask.pickle'
    filename_probabilistic=date_track+'_probability_iceslabs_presence.pickle'
    
    #Open files
    f_depth_corrected = open(path_depth_corrected+filename_depth_corrected, "rb")
    radar = pickle.load(f_depth_corrected)
    f_depth_corrected.close()
    
    f_mask = open(path_mask+filename_mask, "rb")
    mask = pickle.load(f_mask)
    f_mask.close()
    
    f_probabilistic = open(path_probabilistic+filename_probabilistic, "rb")
    probabilistic_file = pickle.load(f_probabilistic)
    f_probabilistic.close() 
    
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
        radar=np.fliplr(radar)
        mask=np.flipud(mask)
        probabilistic_file=np.fliplr(probabilistic_file)
    
    #Transform the coordinated from WGS84 to EPSG:3413
    #Example from: https://pyproj4.github.io/pyproj/stable/examples.html
    points=transformer.transform(np.array(lon_appended),np.array(lat_appended))
    lon_3413=points[0]
    lat_3413=points[1]
    
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
                                 'lat_3413':lat_3413,
                                 'lon_3413':lon_3413,
                                 'depth':depth,
                                 'radar':radar,
                                 'mask':mask,
                                 'probabilistic':probabilistic_file}


#Load all 2010-2018 data without spatial aggregation
path='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2010_2018/excel_spatial_aggreation_and_other/final_excel/prob00/'
df_2010_2018_csv = pd.read_csv(path+'Ice_Layer_Output_Thicknesses_2010_2018_jullienetal2021_prob00.csv',delimiter=',',decimal='.')
#Transform the coordinated from WGS84 to EPSG:3413
transformer = Transformer.from_crs("EPSG:4326", "EPSG:3413", always_xy=True)
points=transformer.transform(np.asarray(df_2010_2018_csv["lon"]),np.asarray(df_2010_2018_csv["lat"]))

#Store lat/lon in 3413
df_2010_2018_csv['lon_3413']=points[0]
df_2010_2018_csv['lat_3413']=points[1]

#Load 2010-2018 elevation dataset
path_df_with_elevation='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2010_2018/excel_spatial_aggreation_and_other/old_elevation_computation/' 
f_20102018 = open(path_df_with_elevation+'df_20102018_with_elevation_prob00_rignotetalregions', "rb")
df_2010_2018_elevation = pickle.load(f_20102018)
f_20102018.close()

###################### From Tedstone et al., 2022 #####################
#from plot_map_decadal_change.py
# Define the CartoPy CRS object.
crs = ccrs.NorthPolarStereo(central_longitude=-45., true_scale_latitude=70.)
# This can be converted into a `proj4` string/dict compatible with GeoPandas
crs_proj4 = crs.proj4_init
###################### From Tedstone et al., 2022 #####################

#Plot
fig = plt.figure()
gs = gridspec.GridSpec(35, 20)
gs.update(wspace=0.1)
gs.update(hspace=0.1)

ax1r = plt.subplot(gs[0:4, 0:10])
ax2r = plt.subplot(gs[0:4, 10:20])
ax3r = plt.subplot(gs[4:8, 0:10])
ax4r = plt.subplot(gs[8:12, 0:10])
ax5r = plt.subplot(gs[12:16, 0:10])
ax6r = plt.subplot(gs[16:20, 0:10])
ax7r = plt.subplot(gs[20:24, 0:10])
ax11t = plt.subplot(gs[27:35, 0:10])

ax8map = plt.subplot(gs[26:35, 10:20],projection=crs)
ax10m = plt.subplot(gs[7:20, 10:20])

#Display GrIS drainage bassins on the map subplot
SW_rignotetal.plot(ax=ax8map,color='white', edgecolor='black',linewidth=0.5) 
CW_rignotetal.plot(ax=ax8map,color='white', edgecolor='black',linewidth=0.5) 

#Plot thickness change for that case study on axis ax11t, display the radargrams, map and shallowest and deepest slab
min_elev,max_elev=plot_thickness(investigation_year,dataframe,df_2010_2018_elevation,ax11t,my_pal)

#Finalize axis ax2r
ax2r.yaxis.set_label_position("right")
ax2r.yaxis.tick_right()

#Finalize axis ax11t
ax11t.set_xlim(min_elev,max_elev)
ax11t.set_xlabel('Elevation [m]')
ax11t.set_ylabel('Column ice thickness [m]')
#Activate ticks label
ax11t.xaxis.tick_bottom()
ax11t.yaxis.tick_left()

#Add vertical lines where the analysed section is
ax11t.axvline(x=1882,zorder=1,linestyle='--',color='k') #from QGIS, reading along 2014-2017 flightline
ax11t.axvline(x=1852,zorder=1,linestyle='--',color='k') #from QGIS, reading along 2014-2017 flightline
ax11t.set_xlim(1735.68640136718750,1983.92419433593750)
#xmin is from 20170502_01_171_173 in df_for_elev that is the closest from -47.5
#xmax is from 20140416_05_035_037 in df_for_elev that is the closest from -46.66
#Note that 2014 and 2017 are perfectly overlapping.
ax11t.scatter(1879,15.8,s=10,c='r')
#Add pannel label
ax11t.text(0.01, 0.875,'h',ha='center', va='center', transform=ax11t.transAxes,weight='bold',fontsize=20)#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot

#Finalize radargrams plot
ax7r.set_xlabel('Longitude [°]')
ax4r.set_ylabel('Depth [m]')

#Finalize map plot
#Display flightlines correpsondance with year
ax8map.text(-100400,-2505000,u'\u03B2')
ax8map.text(-94240,-2554000,u'\u03B1')
ax8map.text(-79280,-2522000,u'\u03B3')
ax8map.text(-79280,-2533000,s=u'\u03B4')
#Show KAN_U
#Show pannel numbers on the map
ax8map.scatter(-89205.404,-2522571.489,s=15,c='#b2182b',label='KAN_U',zorder=10)
#Add pannel label
ax8map.text(-114400,-2505000,'j',ha='center', va='center',fontsize=20)#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot

###################### From Tedstone et al., 2022 #####################
#from plot_map_decadal_change.py
# x0, x1, y0, y1
ax8map.set_extent([-114500, -70280, -2556000, -2495000], crs=crs)
gl=ax8map.gridlines(draw_labels=True, xlocs=[-47, -47.5], ylocs=[67], x_inline=False, y_inline=False,linewidth=0.5)
#Customize lat labels
gl.ylabels_right = False
gl.xlabels_bottom = False
ax8map.axis('off')
#ax8map.legend(loc='upper right')
###################### From Tedstone et al., 2022 #####################

figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()

#Load KAN_U data
path_KAN_U_data='C:/Users/jullienn/switchdrive/Private/research/RT1/KAN_U_data/'
df_KAN_U_csv = pd.read_csv(path_KAN_U_data+'KAN_U_hourly_v03_csv.csv',sep=';',decimal=',',header=0,na_values=-999)

#1. compute hourly melt

#Melt equation: M=SWd-SWup+LWdown-LWup+QH+QL W/m2. 1W=1J/s
df_melt=df_KAN_U_csv
df_melt['melt']=df_KAN_U_csv['ShortwaveRadiationDown_Cor(W/m2)']-df_KAN_U_csv['ShortwaveRadiationUp_Cor(W/m2)']+df_KAN_U_csv['LongwaveRadiationDown(W/m2)']-df_KAN_U_csv['LongwaveRadiationUp(W/m2)']+df_KAN_U_csv['SensibleHeatFlux(W/m2)']+df_KAN_U_csv['LatentHeatFlux(W/m2)']

#2. compute hourly NRJ [J]

#Transform flux into energy: 1W = 1J/s => W*s=J
df_melt['NRJ']=df_melt['melt']*60*60

#3. discard negative melt
df_melt_positive=df_melt
df_melt_positive[df_melt_positive['NRJ']<0]=np.nan

#4. Keep only summer data

#Select only from 1st May to 30th September
df_melt_summer=df_melt_positive[(df_melt_positive['MonthOfYear']>=5) & (df_melt_positive['MonthOfYear']<=9)]

#Keep only data from 2009 to 2017
df_melt_summer_time_period=df_melt_summer[(df_melt_summer['Year']>=2009) & (df_melt_summer['Year']<=2017)]

#Transform data from J to kJ => /1000
df_melt_summer_time_period['NRJ_kJ']=df_melt_summer_time_period['NRJ']/1000

#Set the year as integer for data display
df_melt_summer_time_period.Year=df_melt_summer_time_period.Year.astype(int)

#5. show total cumulative melt
ax = sns.barplot(x="Year", y="NRJ_kJ", data=df_melt_summer_time_period,palette=['lightgrey'],ax=ax10m,estimator=sum,ci=None)

#Set y tick to the right
ax10m.yaxis.set_label_position("right")
ax10m.yaxis.tick_right()
ax10m.set_ylabel('Melt energy availaibility [kJ]')
ax10m.set_xlabel('Year')
#Activate ticks xlabel
ax10m.xaxis.tick_bottom()
#Add pannel label
ax10m.text(0.01, 0.875,'i',ha='center', va='center', transform=ax10m.transAxes,weight='bold',fontsize=20)#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot

'''
ax10m.legend_.remove()
plt.show()
'''
######################## Assess the data coverage #############################

#Create an empty dataframe for data coverage
year_coverage=[]
month_coverage=[]
day_coverage=[]
data_coverage=[]

for indiv_year in df_melt['Year'].unique():
    #Select indiv year
    df_melt_year=df_melt[df_melt['Year']==indiv_year]
    
    for indiv_month in df_melt_year['MonthOfYear'].unique():
        #Select indiv month
        df_melt_year_month=df_melt_year[df_melt_year['MonthOfYear']==indiv_month]
    
        for indiv_day in df_melt_year_month['DayOfMonth'].unique():
            #Select indiv day
            df_melt_year_month_day=df_melt_year_month[df_melt_year_month['DayOfMonth']==indiv_day]
            
            #Store number of non-NaN data I have for this day (complete = 24)
            indiv_day_data=np.asarray(df_melt_year_month_day['melt'])
            
            #Get rid of NaNs
            indiv_day_data_nonNaNs=indiv_day_data[~np.isnan(indiv_day_data)]
            
            #Store associated number of melt data for this day
            year_coverage=np.append(year_coverage,indiv_year)
            month_coverage=np.append(month_coverage,indiv_month)
            day_coverage=np.append(day_coverage,indiv_day)
            data_coverage=np.append(data_coverage,len(indiv_day_data_nonNaNs))

#Create a coverage dataframe
df_coverage=pd.DataFrame(year_coverage, columns=['Year'])
df_coverage['Month']=month_coverage
df_coverage['Day']=day_coverage
df_coverage['Data_coverage']=data_coverage

pdb.set_trace()

######################## Assess the data coverage #############################


#Save figure
plt.savefig('C:/Users/jullienn/switchdrive/Private/research/RT1/figures/fig3/v2/fig3.png',dpi=1000)

pdb.set_trace()

indiv_file_load='Data_20120423_01_137.mat'
path_raw_data=path_data+str(2012)+'_Greenland_P3/CSARP_qlook/'+'20120423_01'+'/'

fdata_filename = scipy.io.loadmat(path_raw_data+indiv_file_load)
lon_filename = fdata_filename['Longitude']


