# -*- coding: utf-8 -*-
"""
Created on Sun Dec 19 12:14:06 2021

@author: jullienn
"""
def plot_thickness_evolution(dictionnary_case_study,df_2010_2018_csv):
    
    #Desired number of slices
    desired_nb=20
    
    #Create empty dataframe for storing data
    df_sampling=pd.DataFrame(columns=['Track_name','year','low_bound', 'high_bound', 'bound_nb', 'mean', 'stddev', '20m_ice_content_m'])
    
    #Loop over the years
    for year in dictionnary_case_study.keys():
        if (dictionnary_case_study[year] == 'empty'):
            continue
        
        #Select data for the trace
        df_trace=df_2010_2018_csv[df_2010_2018_csv['Track_name']==dictionnary_case_study[year][0][5:20]+'_'+dictionnary_case_study[year][-1][17:20]]
        
        #Define the longitudinal sampling THIS WORKS ONLY FOR NEGATIVE LON SO FAR!!!!
        lon_divide=np.arange(np.floor(np.min(df_trace['lon_3413'])),(np.floor(np.max(df_trace['lon_3413']))+1)+(np.abs(np.floor(np.min(df_trace['lon_3413'])))-np.abs(np.floor(np.max(df_trace['lon_3413']))+1))/desired_nb,(np.abs(np.floor(np.min(df_trace['lon_3413'])))-np.abs(np.floor(np.max(df_trace['lon_3413']))+1))/desired_nb)
        
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
    
    #plot data
    fig, ax = plt.subplots()
    fig.suptitle('Case study')
    ax = sns.boxplot(x="bound_nb", y="20m_ice_content_m", hue="year",
                     data=df_sampling, palette="Set3")
    
    plt.show()
    
    return

def plot_thickness_high_end(df_2010_2018,df_recent,df_old,elevDem,grid,slice_lon_summary,lat_slices,list_high_end):
    
    #Plot differences
    diff_to_plot=df_recent-df_old
    pos_diff_to_plot=diff_to_plot
    pos_diff_to_plot[pos_diff_to_plot.avg_20m_icecontent<0]=np.nan
    
    fig, (ax1) = plt.subplots(1, 1)#, gridspec_kw={'width_ratios': [1, 3]})
    fig.suptitle('Spatial aggregation, positive difference '+str(np.unique(df_recent['year'])[0])+'-'+str(np.unique(df_old['year'])[0]))
    
    # Make the norm for difference plotting
    #divnorm_diff = mcolors.DivergingNorm(vmin=0, vcenter=5, vmax=10)
    '''
    #Display 2017 and 2018 data
    plt.scatter(df_2010_2018_csv[df_2010_2018_csv.Track_name.str[:4]=='2017']['lon_3413'],df_2010_2018_csv[df_2010_2018_csv.Track_name.str[:4]=='2017']['lat_3413'],s=0.1,color='#737373')
    plt.scatter(df_2010_2018_csv[df_2010_2018_csv.Track_name.str[:4]=='2018']['lon_3413'],df_2010_2018_csv[df_2010_2018_csv.Track_name.str[:4]=='2018']['lat_3413'],s=0.1,color='#737373',label='2017-2018 ice slabs')
    
    #Display 2011 and 2012 data
    plt.scatter(df_2010_2018_csv[df_2010_2018_csv.Track_name.str[:4]=='2011']['lon_3413'],df_2010_2018_csv[df_2010_2018_csv.Track_name.str[:4]=='2011']['lat_3413'],s=0.1,color='#969696')
    plt.scatter(df_2010_2018_csv[df_2010_2018_csv.Track_name.str[:4]=='2012']['lon_3413'],df_2010_2018_csv[df_2010_2018_csv.Track_name.str[:4]=='2012']['lat_3413'],s=0.1,color='#969696',label='2011-2012 ice slabs')
    '''
    
    #Display the difference between 2011 and 2010 if aggregated data
    sc= ax1.scatter(df_recent['avg_lon_3413'],df_recent['avg_lat_3413'],c=pos_diff_to_plot['avg_20m_icecontent'],cmap=discrete_cmap(10,'Blues'))#,norm=divnorm_diff)
    cbar=fig.colorbar(sc)
    cbar.set_label('Difference in iceslabs thickness', fontsize=15)
    
    #Add upper iceslabs
    if ('2002-2003' in list_high_end):
        ax1.step(slice_lon_summary[:,0],lat_slices,color='#fee5d9',label='2002-2003')
    if ('2010' in list_high_end):
        ax1.step(slice_lon_summary[:,1],lat_slices,color='#fcae91',label='2010')
    if ('2011-2012' in list_high_end):
        ax1.step(slice_lon_summary[:,2],lat_slices,color='#fb6a4a',label='2011-2012')
    if ('2013-2014' in list_high_end):
        ax1.step(slice_lon_summary[:,3],lat_slices,color='#de2d26',label='2013-2014')
    if ('2017-2018' in list_high_end):
        ax1.step(slice_lon_summary[:,4],lat_slices,color='#a50f15',label='2017-2018')
    
    plt.legend()
    
    ax1.set_xlim(-240000,-65000)
    ax1.set_ylim(-2650000,-2250000)
    
    #Allows to open plot in full size directly
    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()
    
    plt.show()
    


###     This is from iceslabs_20102018_thickening_analysis.py       ###

#Import librairies
import datetime
from scipy import spatial
import pandas as pd
from pyproj import Transformer
import numpy as np
import pdb
import matplotlib.pyplot as plt
import geopandas as gpd
import seaborn as sns
sns.set_theme(style="whitegrid")

### -------------------------- Load shapefiles --------------------------- ###
path_regional_masks='C:/Users/jullienn/switchdrive/Private/research/backup_Aglaja/working_environment/greenland_topo_data/masks_for_2002_2003_calculations'

NW_icecap_greenland_mask=gpd.read_file(path_regional_masks+'/NW_icecap_greenland_mask_3413.shp')
NW_north_greenland_mask=gpd.read_file(path_regional_masks+'/NW_north_greenland_mask_3413.shp')
NW_west_greenland_mask=gpd.read_file(path_regional_masks+'/NW_west_greenland_mask_3413.shp')
SW_lower_greenland_mask=gpd.read_file(path_regional_masks+'/SW_lower_greenland_mask_3413.shp')
SW_middle_greenland_mask=gpd.read_file(path_regional_masks+'/SW_middle_greenland_mask_3413.shp')
SW_upper_greenland_mask=gpd.read_file(path_regional_masks+'/SW_upper_greenland_mask_3413.shp')

#Load Rignot et al., 2016 Greenland drainage bassins
path_rignotetal2016_GrIS_drainage_bassins='C:/Users/jullienn/switchdrive/Private/research/backup_Aglaja/working_environment/greenland_topo_data/GRE_Basins_IMBIE2_v1.3/'
GrIS_drainage_bassins=gpd.read_file(path_rignotetal2016_GrIS_drainage_bassins+'GRE_Basins_IMBIE2_v1.3_EPSG_3413.shp',rows=slice(51,57,1)) #the regions are the last rows of the shapefile

#Extract indiv regions and create related indiv shapefiles
NO_rignotetal=GrIS_drainage_bassins[GrIS_drainage_bassins.SUBREGION1=='NO']
NE_rignotetal=GrIS_drainage_bassins[GrIS_drainage_bassins.SUBREGION1=='NE']
SE_rignotetal=GrIS_drainage_bassins[GrIS_drainage_bassins.SUBREGION1=='SE']
SW_rignotetal=GrIS_drainage_bassins[GrIS_drainage_bassins.SUBREGION1=='SW']
CW_rignotetal=GrIS_drainage_bassins[GrIS_drainage_bassins.SUBREGION1=='CW']
NW_rignotetal=GrIS_drainage_bassins[GrIS_drainage_bassins.SUBREGION1=='NW']
### -------------------------- Load shapefiles --------------------------- ###

#Load the spatial aggregated data. All the points within a radius of 100m are averaged
path='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2010_2018/excel_spatial_aggreation_and_other/final_excel/prob00/'
df_2010_2018_spatially_aggregated = pd.read_csv(path+'jullien_etal_20102018_spatial_aggregation_grid_1000_prob00.csv',delimiter=';',decimal=',')

#Load all 2010-2018 data without spatial aggregation
df_2010_2018_csv = pd.read_csv(path+'Ice_Layer_Output_Thicknesses_2010_2018_jullienetal2021_prob00.csv',delimiter=',',decimal='.')
#Transform the coordinated from WGS84 to EPSG:3413
transformer = Transformer.from_crs("EPSG:4326", "EPSG:3413", always_xy=True)
points=transformer.transform(np.asarray(df_2010_2018_csv["lon"]),np.asarray(df_2010_2018_csv["lat"]))

#Store lat/lon in 3413
df_2010_2018_csv['lon_3413']=points[0]
df_2010_2018_csv['lat_3413']=points[1]

#Loop over the keys and create 1 dataframe per year. Where no data for this particular year, store a nan
#Create empty arrays
nan_array=np.zeros((1,7))
nan_array[:]=np.nan

array_2010=np.zeros((1,7))
array_2010[:]=np.nan

array_2011=np.zeros((1,7))
array_2011[:]=np.nan

array_2012=np.zeros((1,7))
array_2012[:]=np.nan

array_2013=np.zeros((1,7))
array_2013[:]=np.nan

array_2014=np.zeros((1,7))
array_2014[:]=np.nan

array_2017=np.zeros((1,7))
array_2017[:]=np.nan

array_2018=np.zeros((1,7))
array_2018[:]=np.nan

#Loop over the keys
for indiv_key in np.unique(df_2010_2018_spatially_aggregated.key):
    
    #Select data where key
    df_key=df_2010_2018_spatially_aggregated[df_2010_2018_spatially_aggregated.key==indiv_key]
    
    # ----- 2010
    if (not(2010 in np.asarray(df_key.year))):
        #No data for this year, store NaNs
        array_2010=np.append(array_2010,nan_array,axis=0)
    else:
        #There are data for this year, store them
        array_2010=np.append(array_2010,np.asarray(df_key[df_key.year==2010]),axis=0)
    # ----- 2010
    
    # ----- 2011
    if (not(2011 in np.asarray(df_key.year))):
        #No data for this year, store NaNs
        array_2011=np.append(array_2011,nan_array,axis=0)
    else:
        #There are data for this year, store them
        array_2011=np.append(array_2011,np.asarray(df_key[df_key.year==2011]),axis=0)
    # ----- 2011
    
    # ----- 2012
    if (not(2012 in np.asarray(df_key.year))):
        #No data for this year, store NaNs
        array_2012=np.append(array_2012,nan_array,axis=0)
    else:
        #There are data for this year, store them
        array_2012=np.append(array_2012,np.asarray(df_key[df_key.year==2012]),axis=0)
    # ----- 2012
    
    # ----- 2013
    if (not(2013 in np.asarray(df_key.year))):
        #No data for this year, store NaNs
        array_2013=np.append(array_2013,nan_array,axis=0)
    else:
        #There are data for this year, store them
        array_2013=np.append(array_2013,np.asarray(df_key[df_key.year==2013]),axis=0)
    # ----- 2013
    
    # ----- 2014
    if (not(2014 in np.asarray(df_key.year))):
        #No data for this year, store NaNs
        array_2014=np.append(array_2014,nan_array,axis=0)
    else:
        #There are data for this year, store them
        array_2014=np.append(array_2014,np.asarray(df_key[df_key.year==2014]),axis=0)
    # ----- 2014
    
    # ----- 2017
    if (not(2017 in np.asarray(df_key.year))):
        #No data for this year, store NaNs
        array_2017=np.append(array_2017,nan_array,axis=0)
    else:
        #There are data for this year, store them
        array_2017=np.append(array_2017,np.asarray(df_key[df_key.year==2017]),axis=0)
    # ----- 2017
    
    # ----- 2018
    if (not(2018 in np.asarray(df_key.year))):
        #No data for this year, store NaNs
        array_2018=np.append(array_2018,nan_array,axis=0)
    else:
        #There are data for this year, store them
        array_2018=np.append(array_2018,np.asarray(df_key[df_key.year==2018]),axis=0)
    # ----- 2018
    
    print(indiv_key/len(np.unique(df_2010_2018_spatially_aggregated.key))*100,' %')

#Delete the first line of all the array_year because NaNs
array_2010=np.delete(array_2010,0,0)
array_2011=np.delete(array_2011,0,0)
array_2012=np.delete(array_2012,0,0)
array_2013=np.delete(array_2013,0,0)
array_2014=np.delete(array_2014,0,0)
array_2017=np.delete(array_2017,0,0)
array_2018=np.delete(array_2018,0,0)

#Store array as dataframes
df_spatially_aggregated_2010=pd.DataFrame(data=array_2010,
                                          columns=['index', 'avg_20m_icecontent', 'std_20m_icecontent', 'avg_lat_3413','avg_lon_3413', 'year', 'key'])
df_spatially_aggregated_2011=pd.DataFrame(data=array_2011,
                                          columns=['index', 'avg_20m_icecontent', 'std_20m_icecontent', 'avg_lat_3413','avg_lon_3413', 'year', 'key'])
df_spatially_aggregated_2012=pd.DataFrame(data=array_2012,
                                          columns=['index', 'avg_20m_icecontent', 'std_20m_icecontent', 'avg_lat_3413','avg_lon_3413', 'year', 'key'])
df_spatially_aggregated_2013=pd.DataFrame(data=array_2013,
                                          columns=['index', 'avg_20m_icecontent', 'std_20m_icecontent', 'avg_lat_3413','avg_lon_3413', 'year', 'key'])
df_spatially_aggregated_2014=pd.DataFrame(data=array_2014,
                                          columns=['index', 'avg_20m_icecontent', 'std_20m_icecontent', 'avg_lat_3413','avg_lon_3413', 'year', 'key'])
df_spatially_aggregated_2017=pd.DataFrame(data=array_2017,
                                          columns=['index', 'avg_20m_icecontent', 'std_20m_icecontent', 'avg_lat_3413','avg_lon_3413', 'year', 'key'])
df_spatially_aggregated_2018=pd.DataFrame(data=array_2018,
                                          columns=['index', 'avg_20m_icecontent', 'std_20m_icecontent', 'avg_lat_3413','avg_lon_3413', 'year', 'key'])


list_high_end=list(['2002-2003','2010','2011-2012','2013-2014','2017-2018'])


#Plot 2010, 2011, 2012, 2013, 2014 ,2017 2018, select overlapping case study: use clean and clear ice slabs tramsects

#CHOOSE LOC1, loc 2, loc 3, loc 6, 15, 17, 23

loc1={2010:['Data_20100507_01_008.mat','Data_20100507_01_009.mat','Data_20100507_01_010.mat'],
      2011:['Data_20110426_01_009.mat','Data_20110426_01_010.mat','Data_20110426_01_011.mat'],
      2012:'empty',
      2013:'empty',
      2014:['Data_20140421_01_009.mat','Data_20140421_01_010.mat','Data_20140421_01_011.mat','Data_20140421_01_012.mat','Data_20140421_01_013.mat'],
      2017:['Data_20170424_01_008.mat','Data_20170424_01_009.mat','Data_20170424_01_010.mat','Data_20170424_01_011.mat','Data_20170424_01_012.mat','Data_20170424_01_013.mat','Data_20170424_01_014.mat'],
      2018:'empty'}

loc2={2010:['Data_20100513_01_001.mat','Data_20100513_01_002.mat'],
      2011:['Data_20110411_01_116.mat','Data_20110411_01_117.mat','Data_20110411_01_118.mat'],
      2012:['Data_20120428_01_125.mat','Data_20120428_01_126.mat'],
      2013:'empty',
      2014:['Data_20140408_11_024.mat','Data_20140408_11_025.mat','Data_20140408_11_026.mat'],
      2017:['Data_20170508_02_165.mat','Data_20170508_02_166.mat','Data_20170508_02_167.mat','Data_20170508_02_168.mat','Data_20170508_02_169.mat','Data_20170508_02_170.mat','Data_20170508_02_171.mat'],
      2018:'empty'}

#This one is collocated with FS1, 2, 3.
loc3={2010:'empty',
      2011:'empty',
      2012:['Data_20120423_01_137.mat','Data_20120423_01_138.mat'],
      2013:['Data_20130409_01_010.mat','Data_20130409_01_011.mat','Data_20130409_01_012.mat'],
      2014:'empty',
      2017:'empty',
      2018:['Data_20180421_01_004.mat','Data_20180421_01_005.mat','Data_20180421_01_006.mat','Data_20180421_01_007.mat']}


fig, (ax1) = plt.subplots(1, 1)#, gridspec_kw={'width_ratios': [1, 3]})

#Display GrIS drainage bassins
NO_rignotetal.plot(ax=ax1,color='white', edgecolor='black')
NE_rignotetal.plot(ax=ax1,color='white', edgecolor='black') 
SE_rignotetal.plot(ax=ax1,color='white', edgecolor='black') 
SW_rignotetal.plot(ax=ax1,color='white', edgecolor='black') 
CW_rignotetal.plot(ax=ax1,color='white', edgecolor='black') 
NW_rignotetal.plot(ax=ax1,color='white', edgecolor='black') 

plt.scatter(df_spatially_aggregated_2010['avg_lon_3413'],df_spatially_aggregated_2010['avg_lat_3413'],c=df_spatially_aggregated_2010['avg_20m_icecontent'],s=0.2)

#Display
plt.scatter(df_2010_2018_csv[df_2010_2018_csv['Track_name']==loc1[2010][0][5:20]+'_'+loc1[2010][2][17:20]]['lon_3413'],
            df_2010_2018_csv[df_2010_2018_csv['Track_name']==loc1[2010][0][5:20]+'_'+loc1[2010][2][17:20]]['lat_3413'],
            s=0.1,color='#737373')
plt.scatter(df_2010_2018_csv[df_2010_2018_csv['Track_name']==loc1[2011][0][5:20]+'_'+loc1[2011][2][17:20]]['lon_3413'],
            df_2010_2018_csv[df_2010_2018_csv['Track_name']==loc1[2011][0][5:20]+'_'+loc1[2011][2][17:20]]['lat_3413'],
            s=0.1,color='#737373')
plt.scatter(df_2010_2018_csv[df_2010_2018_csv['Track_name']==loc1[2014][0][5:20]+'_'+loc1[2014][2][17:20]]['lon_3413'],
            df_2010_2018_csv[df_2010_2018_csv['Track_name']==loc1[2014][0][5:20]+'_'+loc1[2014][2][17:20]]['lat_3413'],
            s=0.1,color='#737373')
plt.scatter(df_2010_2018_csv[df_2010_2018_csv['Track_name']==loc1[2017][0][5:20]+'_'+loc1[2017][2][17:20]]['lon_3413'],
            df_2010_2018_csv[df_2010_2018_csv['Track_name']==loc1[2017][0][5:20]+'_'+loc1[2017][2][17:20]]['lat_3413'],
            s=0.1,color='#737373')


plot_thickness_evolution(loc1,df_2010_2018_csv)



pdb.set_trace()


plot_thickness_high_end(df_2010_2018,df_spatially_aggregated_2017,df_spatially_aggregated_2010,elevDem,grid,slice_lon_summary,lat_slices,list_high_end)

'''
list_high_end=list(['2002-2003','2011-2012','2017-2018'])
plot_thickness_high_end(df_2010_2018,df_spatially_aggregated_2017,df_spatially_aggregated_2011,elevDem,grid,slice_lon_summary,lat_slices,list_high_end)

list_high_end=list(['2002-2003','2011-2012','2017-2018'])
plot_thickness_high_end(df_2010_2018,df_spatially_aggregated_2017,df_spatially_aggregated_2012,elevDem,grid,slice_lon_summary,lat_slices,list_high_end)

list_high_end=list(['2002-2003','2013-2014','2017-2018'])
plot_thickness_high_end(df_2010_2018,df_spatially_aggregated_2017,df_spatially_aggregated_2013,elevDem,grid,slice_lon_summary,lat_slices,list_high_end)

list_high_end=list(['2002-2003','2013-2014','2017-2018'])
plot_thickness_high_end(df_2010_2018,df_spatially_aggregated_2017,df_spatially_aggregated_2014,elevDem,grid,slice_lon_summary,lat_slices,list_high_end)

list_high_end=list(['2002-2003','2013-2014','2011-2012'])
plot_thickness_high_end(df_2010_2018,df_spatially_aggregated_2013,df_spatially_aggregated_2012,elevDem,grid,slice_lon_summary,lat_slices,list_high_end)

list_high_end=list(['2002-2003','2011-2012','2010'])
plot_thickness_high_end(df_2010_2018,df_spatially_aggregated_2011,df_spatially_aggregated_2010,elevDem,grid,slice_lon_summary,lat_slices,list_high_end)
'''
