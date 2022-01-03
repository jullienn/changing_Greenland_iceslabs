# -*- coding: utf-8 -*-
"""
Created on Sun Nov 28 12:56:32 2021

@author: jullienn

Code adapted from lenses_thickening_visualisation.py
"""

def plot_thickness(dictionnary_case_study,df_2010_2018_csv,axt):
    #This function is adapted from plot_thickness_evolution from fig_2_paper_iceslabs.py
    
    #Define empty dictionnary for longitudinal slice definition
    df_for_lon=pd.DataFrame(columns=list(df_2010_2018_csv.keys()))
    
    #Loop over the years
    for year in dictionnary_case_study.keys():
        if (str(year) in list(['2010','2011'])):
            print('Not on the transect, continue')
            continue
        #Select data for the trace
        df_for_lon_temp=df_2010_2018_csv[df_2010_2018_csv['Track_name']==dictionnary_case_study[year][0][5:20]+'_'+dictionnary_case_study[year][-1][17:20]]
        #Append data to each other
        df_for_lon=df_for_lon.append(df_for_lon_temp)
    
    #Desired number of longitudinal slices
    desired_nb=20
    
    #Create empty dataframe for storing data
    df_sampling=pd.DataFrame(columns=['Track_name','year','low_bound', 'high_bound', 'bound_nb', 'mean', 'stddev', '20m_ice_content_m'])
    
    #Loop over the years
    for year in dictionnary_case_study.keys():
        if (str(year) in list(['2010','2011'])):
            print('Not on the transect, continue')
            continue
        
        #Select data for the trace
        df_trace=df_2010_2018_csv[df_2010_2018_csv['Track_name']==dictionnary_case_study[year][0][5:20]+'_'+dictionnary_case_study[year][-1][17:20]]

        #Define the longitudinal sampling THIS WORKS ONLY FOR NEGATIVE LON SO FAR!!!!
        #lon_divide=np.arange(np.floor(np.min(df_for_lon['lon_3413'])),(np.floor(np.max(df_for_lon['lon_3413']))+1)+(np.abs(np.floor(np.min(df_for_lon['lon_3413'])))-np.abs(np.floor(np.max(df_for_lon['lon_3413']))+1))/desired_nb,(np.abs(np.floor(np.min(df_for_lon['lon_3413'])))-np.abs(np.floor(np.max(df_for_lon['lon_3413']))+1))/desired_nb)
        
        #Lon divide every 4km. I have compared between 2500, 3000, 4000 and 5000m. Best trade off between visualisation and overaggrgation seems to be 4000m
        km_bin_desired=4000
        lon_divide=np.arange(-110240,np.floor(np.max(df_for_lon['lon_3413'])).astype(int)+1+km_bin_desired,km_bin_desired)
        
        #Set bound_nb to 0
        bound_nb=0
        #Loop over the lon divide
        for i in range(1,len(lon_divide)):
            
            #Identify low and higher end of the slice
            low_bound=lon_divide[i-1]
            high_bound=lon_divide[i]
    
            #Select all the data belonging to this lon slice
            ind_slice=np.logical_and(np.array(df_trace['lon_3413']>=low_bound),np.array(df_trace['lon_3413']<high_bound))
            df_select=df_trace[ind_slice]
            
            #Fill in dictionnary
            df_temp=pd.DataFrame(columns=['Track_name','year','low_bound', 'high_bound', 'bound_nb', 'mean', 'stddev', '20m_ice_content_m'])
            df_temp['20m_ice_content_m']=np.asarray(df_select['20m_ice_content_m'])
            df_temp['Track_name']=np.asarray([df_select['Track_name'].unique()]*len(df_select))
            df_temp['year']=np.asarray([year]*len(df_select))
            df_temp['low_bound']=np.asarray([str(low_bound)]*len(df_select))
            df_temp['high_bound']=np.asarray([str(high_bound)]*len(df_select))
            df_temp['bound_nb']=np.asarray([str(bound_nb)]*len(df_select))
            df_temp['mean']=np.asarray([np.nanmean(df_select['20m_ice_content_m'])]*len(df_select))
            df_temp['stddev']=np.asarray([np.nanstd(df_select['20m_ice_content_m'])]*len(df_select))
            
            #Append dictionnary
            df_sampling=df_sampling.append(df_temp)
            
            #Update bound_nb
            bound_nb=bound_nb+1
    
    #Set order to display data
    order_plot=np.arange(np.min(np.asarray(df_sampling['bound_nb']).astype(int)),np.max(np.asarray(df_sampling['bound_nb']).astype(int))+1)
    #Define palette plot
    #This is from https://www.python-graph-gallery.com/33-control-colors-of-boxplot-seaborn
    my_pal = {2010: "#ffffcc", 2011: "#d9f0a3", 2012:"#addd8e", 2013:"#78c679", 2014:"#41ab5d", 2017:"#238443" ,2018:"#005a32"}
    
    #plot thickness data
    sns.boxplot(x="bound_nb", y="20m_ice_content_m", hue="year",data=df_sampling, palette=my_pal, ax=axt,order=order_plot.astype(str))
    
    #Set y tick to the right
    axt.yaxis.set_label_position("right")
    axt.yaxis.tick_right()
    #axt.set_xticklabels([])
    
    #Get rid of legend
    axt.legend_.remove()
    axt.set_xlabel('')
    axt.set_ylabel('')


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


### -------------------------- Load shapefiles --------------------------- ###
#Load Rignot et al., 2016 Greenland drainage bassins
path_rignotetal2016_GrIS_drainage_bassins='C:/Users/jullienn/switchdrive/Private/research/backup_Aglaja/working_environment/greenland_topo_data/GRE_Basins_IMBIE2_v1.3/'
GrIS_drainage_bassins=gpd.read_file(path_rignotetal2016_GrIS_drainage_bassins+'GRE_Basins_IMBIE2_v1.3.shp',rows=slice(51,57,1)) #the regions are the last rows of the shapefile

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
path_df_with_elevation='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2010_2018/excel_spatial_aggreation_and_other/' 
f_20102018 = open(path_df_with_elevation+'df_20102018_with_elevation_prob00_rignotetalregions', "rb")
df_2010_2018_elevation = pickle.load(f_20102018)
f_20102018.close()

#Plot
fig = plt.figure()
gs = gridspec.GridSpec(35, 20)
gs.update(wspace=0.1)
#gs.update(wspace=0.001)

ax1r = plt.subplot(gs[0:5, 0:10])
ax2r = plt.subplot(gs[5:10, 0:10])
ax3r = plt.subplot(gs[10:15, 0:10])
ax4r = plt.subplot(gs[15:20, 0:10])
ax5r = plt.subplot(gs[20:25, 0:10])
ax6r = plt.subplot(gs[25:30, 0:10])
ax7r = plt.subplot(gs[30:35, 0:10])

ax8map = plt.subplot(gs[0:5, 10:20])
ax9l = plt.subplot(gs[5:15, 10:20])
ax10m = plt.subplot(gs[15:25, 10:20])
ax11t = plt.subplot(gs[25:35, 10:20])

#Display GrIS drainage bassins on the map subplot
SW_rignotetal.plot(ax=ax8map,color='white', edgecolor='black') 
CW_rignotetal.plot(ax=ax8map,color='white', edgecolor='black') 

#Plot thickness change for that case study on axis ax11t
plot_thickness(investigation_year,df_2010_2018_csv,ax11t)

#Finalize axis ax11t
ax11t.set_xticklabels(np.arange(0,10*4,4))
ax11t.set_xlabel('Longitude [km]')
ax11t.set_ylabel('Ice slabs thickness [m]')

pdb.set_trace()

#Display the data        
for single_year in investigation_year.keys():
    
    print(str(single_year))
    
    #If no data, continue
    if (investigation_year[single_year]=='empty'):
        print('No data for year '+str(single_year)+', continue')
        continue
    
    start_date_track=investigation_year[single_year][0]
    end_date_track=investigation_year[single_year][-1]
    date_track=start_date_track[5:20]+'_'+end_date_track[17:20]
    
    '''
    #Calculate distances (in m)
    distances=compute_distances(dataframe[str(single_year)]['lon_appended'],dataframe[str(single_year)]['lat_appended'])
    '''
        
    #Reset depths to 0
    dataframe[str(single_year)]['depth']=dataframe[str(single_year)]['depth']-dataframe[str(single_year)]['depth'][0]
    
    #Select radar slice
    depth_corrected_file=dataframe[str(single_year)]['radar']
    
    #Identify index where time < 20 m
    ind_lower_20m=np.where(dataframe[str(single_year)]['depth']<20)[0]
    depth_corrected_20m=depth_corrected_file[ind_lower_20m,:]
    
    #Identify axis for plotting
    if (single_year==2010):
        ax_plotting=ax1r
        color_toplot="#ffffcc"
    elif (single_year==2011):
        ax_plotting=ax2r
        color_toplot="#d9f0a3"
    elif (single_year==2012):
        ax_plotting=ax3r
        color_toplot="#addd8e"
    elif (single_year==2013):
        ax_plotting=ax4r
        color_toplot="#78c679"
    elif (single_year==2014):
        ax_plotting=ax5r
        color_toplot="#41ab5d"
    elif (single_year==2017):
        ax_plotting=ax6r
        color_toplot="#238443"
    elif (single_year==2018):
        ax_plotting=ax7r
        color_toplot="#005a32"
    else:
        print('Year not existing')
    
    #ax_plotting.set_aspect(4)    
    X=dataframe[str(single_year)]['lon_appended']
    Y=np.arange(0,100,100/dataframe[str(single_year)]['radar'].shape[0])
    C=dataframe[str(single_year)]['radar']
            
    cb=ax_plotting.pcolor(X, Y, C,cmap=plt.get_cmap('gray'))#,norm=divnorm)
    ax_plotting.invert_yaxis() #Invert the y axis = avoid using flipud.    
    ax_plotting.set_ylim(20,0)
    plt.show()
    ax_plotting.set_xlim(-48,-46.5655)
    
    #Display radar trace on map with mask applied on data
    ax8map.scatter(dataframe[str(single_year)]['lon_appended'][dataframe[str(single_year)]['mask']],dataframe[str(single_year)]['lat_appended'][dataframe[str(single_year)]['mask']],c=color_toplot,s=0.1)
    
    ##########################################################################
    ###                        Extract ice content                         ###
    ##########################################################################

    #Extract ice content
    indiv_probability_slice=dataframe[str(single_year)]['probabilistic']
    
    #Compute depth_delta_m
    depth_delta_m = np.mean(dataframe[str(single_year)]['depth'][1:] - dataframe[str(single_year)]['depth'][:-1])
    
    #Let's transform the probabilistic ice slabs into an ice content
    #We must derive a low end and high end of ice slabs likelihood
    #for low end: slabs identified in 19 quantiles out of 19 => likelihood = 19/19=1
    #for high end: slabs identified in 1 quantile out of 19 => likelihood = 1/19 = 0.05263
    index_prob=indiv_probability_slice>=0.1579 # >3/19
    
    #Create slice full of nans
    slice_for_calculation=np.zeros((indiv_probability_slice.shape[0],indiv_probability_slice.shape[1]))
    
    #fill in slice_for_calculation by ones where likelihood >= 0.5
    slice_for_calculation[index_prob]=1

    # Number of pixels times the thickness of each pixel
    ice_content_m = np.sum(slice_for_calculation, axis=0) * depth_delta_m
    
    #Moving window to average results
    ice_content_m_avg= np.round(np.convolve(ice_content_m, np.ones(50)/50, mode='same'))
    
    #Store total ice content
    dataframe[str(single_year)]['ice_content']=ice_content_m_avg
    
    ##########################################################################
    ###                        Extract ice content                         ###
    ##########################################################################   
    
    if (str(single_year) in list(['2010','2011'])):
        continue
    else:
        ax9l.plot(X,ice_content_m_avg,label=str(single_year),color=color_toplot)
        plt.legend()
        ax9l.set_xlim(-47.427,-46.5655)
        plt.show()
 
figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()
plt.show()

#Finalize map plot
ax8map.set_xlim(-48.25,-46.10)
ax8map.set_ylim(66.8,67)





pdb.set_trace()

indiv_file_load='Data_20120423_01_137.mat'
path_raw_data=path_data+str(2012)+'_Greenland_P3/CSARP_qlook/'+'20120423_01'+'/'

fdata_filename = scipy.io.loadmat(path_raw_data+indiv_file_load)
lon_filename = fdata_filename['Longitude']








