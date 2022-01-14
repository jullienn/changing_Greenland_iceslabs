# -*- coding: utf-8 -*-
"""
Created on Sun Dec 19 12:14:06 2021

@author: jullienn
"""
def plot_thickness_evolution(dictionnary_case_study,df_2010_2018_csv,df_2010_2018_elevation,ax1,axt,custom_angle,offset_x,offset_y,casestudy_nb):
    
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
                
        #Display data
        ax1.scatter(df_2010_2018_elevation[df_2010_2018_elevation['Track_name']==dictionnary_case_study[year][0][5:20]+'_'+dictionnary_case_study[year][-1][17:20]]['lon_3413'],
                    df_2010_2018_elevation[df_2010_2018_elevation['Track_name']==dictionnary_case_study[year][0][5:20]+'_'+dictionnary_case_study[year][-1][17:20]]['lat_3413'],
                    s=0.1,color='#737373')
            
    #Display rectangle around data    
    x=(np.min(df_for_elev.lon_3413)-offset_x)
    y=(np.min(df_for_elev.lat_3413)-offset_y)
    width=10000
    height=np.sqrt(np.power(abs(np.min(df_for_elev.lon_3413))-abs(np.max(df_for_elev.lon_3413)),2)+np.power(abs(np.min(df_for_elev.lat_3413))-abs(np.max(df_for_elev.lat_3413)),2))+2*offset_x
    #This is from https://stackoverflow.com/questions/37435369/matplotlib-how-to-draw-a-rectangle-on-image
    # Create a Rectangle patch
    rect = patches.Rectangle((x,y),width,height, angle=custom_angle, linewidth=1, edgecolor='blue', facecolor='none')
    # Add the patch to the Axes
    ax1.add_patch(rect)
    
    if (casestudy_nb=='a'):
        #Add number of case study on fig localisation
        ax1.text(x+30000,y-40000,casestudy_nb,color='r',weight='bold',fontsize=12)
    elif (casestudy_nb=='c'):
        #Add number of case study on fig localisation
        ax1.text(x-15000,y+30000,casestudy_nb,color='r',weight='bold',fontsize=12)
    else:
        #Add number of case study on fig localisation
        ax1.text(x-35000,y-15000,casestudy_nb,color='r',weight='bold',fontsize=12)

    #Add number of case study on fig localisation    
    axt.text(0.005, 0.875,casestudy_nb, ha='center', va='center', transform=axt.transAxes,weight='bold',fontsize=20)#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot
    
    #Define palette for time periods
    #This is from https://www.python-graph-gallery.com/33-control-colors-of-boxplot-seaborn
    my_pal = {'2010': "#fdd49e", '2011-2012': "#fc8d59", '2013-2014':"#d7301f",'2017-2018':"#7f0000"}
    
    #Create an empty df_sampling
    df_sampling=pd.DataFrame(columns=['Track_name','time_period','low_bound', 'high_bound', 'bound_nb', 'mean', 'median', 'q025', 'q075','stddev'])

    #Elev divide every 1m.
    elev_bin_desired=1
    elev_divide=np.arange(np.floor(np.min(df_for_elev['elevation'])),np.floor(np.max(df_for_elev['elevation']))+1+elev_bin_desired,elev_bin_desired)
    
    #Define window size for smoothing
    winsize=10
    
    #Loop over the different time periods (2010, 2011-2012, 2013-2014, 2017-2018)
    for time_period in list(['2010','2011-2012','2013-2014','2017-2018']):
        
        #Get data for that specific time period
        if (time_period == '2010'):
            df_trace_year=df_for_elev[df_for_elev['year']==2010]
        elif (time_period == '2011-2012'):
            df_trace_year=df_for_elev[(df_for_elev['year']>=2011) & (df_for_elev['year']<=2012)]
        elif (time_period == '2013-2014'):
            df_trace_year=df_for_elev[(df_for_elev['year']>=2013) & (df_for_elev['year']<=2014)]
        elif (time_period == '2017-2018'):
            df_trace_year=df_for_elev[(df_for_elev['year']>=2017) & (df_for_elev['year']<=2018)]
        else:
            print('Time period not known, break')
            break
        
        if (len(df_trace_year)==0):
            #No data in this time period, continue
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
                df_temp['time_period']=time_period
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
        
    for time_period in list(['2010','2011-2012','2013-2014','2017-2018']):
        
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

    #Get rid of legend
    #axt.legend_.remove()
    axt.set_xlabel('')
    axt.set_ylabel('')
    plt.show()
    
    print('End plotting fig 2')
    return

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
import pickle
import matplotlib.patches as patches
import matplotlib as mpl
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
import seaborn as sns
sns.set_theme(style="whitegrid")
from scipy import signal
import cartopy.crs as ccrs

#Set fontsize plot
plt.rcParams.update({'font.size': 10})

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

path='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2010_2018/excel_spatial_aggreation_and_other/final_excel/high_estimate/'

#Load all 2010-2018 data without spatial aggregation
df_2010_2018_csv = pd.read_csv(path+'Ice_Layer_Output_Thicknesses_2010_2018_jullienetal2021_high_estimate.csv',delimiter=',',decimal='.')
#Transform the coordinated from WGS84 to EPSG:3413
transformer = Transformer.from_crs("EPSG:4326", "EPSG:3413", always_xy=True)
points=transformer.transform(np.asarray(df_2010_2018_csv["lon"]),np.asarray(df_2010_2018_csv["lat"]))

#Store lat/lon in 3413
df_2010_2018_csv['lon_3413']=points[0]
df_2010_2018_csv['lat_3413']=points[1]

'''
#Load the spatial aggregated data. All the points within a radius of 100m are averaged
df_2010_2018_spatially_aggregated = pd.read_csv(path+'jullien_etal_20102018_spatial_aggregation_grid_1000_prob00.csv',delimiter=';',decimal=',')

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
'''

#Plot 2010, 2011, 2012, 2013, 2014 ,2017 2018, select overlapping case study: use clean and clear ice slabs tramsects

#CHOOSE LOC1, loc 2, loc 3, loc 15, loc 23

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

'''
#in 2017 overestimation of ice content
loc4={2010:['Data_20100512_04_073.mat','Data_20100512_04_074.mat'],
      2011:'empty',
      2012:'empty',
      2013:'empty',
      2014:'empty',
      2017:['Data_20170421_01_171.mat','Data_20170421_01_172.mat','Data_20170421_01_173.mat','Data_20170421_01_174.mat'],
      2018:['Data_20180425_01_166.mat','Data_20180425_01_167.mat','Data_20180425_01_168.mat','Data_20180425_01_169.mat']}
'''
'''
#Only 2 years in SW, do not display it
loc5={2010:'empty',
      2011:'empty',
      2012:['Data_20120412_01_095.mat'],
      2013:'empty',
      2014:'empty',
      2017:'empty',
      2018:['Data_20180421_01_174.mat','Data_20180421_01_175.mat','Data_20180421_01_176.mat','Data_20180421_01_177.mat']}
'''
loc6={2010:'empty',
      2011:['Data_20110516_01_009.mat','Data_20110516_01_010.mat'],
      2012:'empty',
      2013:['Data_20130402_01_008.mat'],
      2014:'empty',
      2017:['Data_20170412_01_075.mat','Data_20170412_01_076.mat'],
      2018:'empty'}

'''
#loc7 is in NE but is curved
loc7={2010:'empty',
      2011:'empty',
      2012:'empty',
      2013:'empty',
      2014:['Data_20140508_03_019.mat','Data_20140508_03_020.mat','Data_20140508_03_021.mat','Data_20140508_03_022.mat','Data_20140508_03_023.mat','Data_20140508_03_024.mat'],
      2017:['Data_20170328_01_095.mat','Data_20170328_01_096.mat','Data_20170328_01_097.mat','Data_20170328_01_098.mat','Data_20170328_01_099.mat','Data_20170328_01_100.mat','Data_20170328_01_101.mat'],
      2018:'empty'}
'''

loc8={2010:['Data_20100517_02_001.mat','Data_20100517_02_002.mat'],
      2011:['Data_20110502_01_171.mat'],
      2012:['Data_20120516_01_002.mat'],
      2013:['Data_20130419_01_004.mat','Data_20130419_01_005.mat'],
      2014:['Data_20140507_03_007.mat','Data_20140507_03_008.mat'], #test with 20140514_02_087_089 and 20140515_02_173_175 also
      2017:['Data_20170417_01_171.mat','Data_20170417_01_172.mat','Data_20170417_01_173.mat','Data_20170417_01_174.mat'],
      2018:'empty'}

'''
20100517_02_001_002, 20100519_01_005_005
20110509_01_177_177, 20110502_01_171_171
20120516_01_002_002, 20120330_01_124_125, 20120516_01_115_115
20130419_01_004_005
20140507_03_007_008, 20140514_02_087_089, #20140519_02_002_004 diverging at the start. do not consider 20140429_02_160_161
20170417_01_171_174
'''

loc9={2010:['Data_20100508_01_114.mat','Data_20100508_01_115.mat'],
      2011:['Data_20110419_01_008.mat','Data_20110419_01_009.mat','Data_20110419_01_010.mat'],
      2012:['Data_20120418_01_129.mat','Data_20120418_01_130.mat','Data_20120418_01_131.mat'],
      2013:['Data_20130405_01_165.mat','Data_20130405_01_166.mat','Data_20130405_01_167.mat'],
      2014:['Data_20140424_01_002.mat','Data_20140424_01_003.mat','Data_20140424_01_004.mat'],
      2017:['Data_20170422_01_168.mat','Data_20170422_01_169.mat','Data_20170422_01_170.mat','Data_20170422_01_171.mat'],
      2018:['Data_20180427_01_170.mat','Data_20180427_01_171.mat','Data_20180427_01_172.mat']}


###################### From Tedstone et al., 2022 #####################
#from plot_map_decadal_change.py
# Define the CartoPy CRS object.
crs = ccrs.NorthPolarStereo(central_longitude=-45., true_scale_latitude=70.)

# This can be converted into a `proj4` string/dict compatible with GeoPandas
crs_proj4 = crs.proj4_init
###################### From Tedstone et al., 2022 #####################
        
fig = plt.figure(figsize=(32,48))
gs = gridspec.GridSpec(40, 20)
gs.update(wspace=0.001)
#gs.update(wspace=0.001)
#projection set up from https://stackoverflow.com/questions/33942233/how-do-i-change-matplotlibs-subplot-projection-of-an-existing-axis
ax1 = plt.subplot(gs[0:25, 0:2],projection=crs)
ax_legend = plt.subplot(gs[35:40, 0:2])

ax2t = plt.subplot(gs[0:5, 3:20])
ax3t = plt.subplot(gs[7:12, 3:20])
ax4t = plt.subplot(gs[14:19, 3:20])
ax5t = plt.subplot(gs[21:26, 3:20])
ax6t = plt.subplot(gs[28:33, 3:20])
ax7t = plt.subplot(gs[35:40, 3:20])

ax1.set_facecolor('white')

#Display GrIS drainage bassins
NO_rignotetal.plot(ax=ax1,color='white', edgecolor='black',linewidth=0.5)
NE_rignotetal.plot(ax=ax1,color='white', edgecolor='black',linewidth=0.5) 
SE_rignotetal.plot(ax=ax1,color='white', edgecolor='black',linewidth=0.5) 
SW_rignotetal.plot(ax=ax1,color='white', edgecolor='black',linewidth=0.5) 
CW_rignotetal.plot(ax=ax1,color='white', edgecolor='black',linewidth=0.5) 
NW_rignotetal.plot(ax=ax1,color='white', edgecolor='black',linewidth=0.5) 

#plt.scatter(df_spatially_aggregated_2010['avg_lon_3413'],df_spatially_aggregated_2010['avg_lat_3413'],c=df_spatially_aggregated_2010['avg_20m_icecontent'],s=0.2)

#Load 2010-2018 elevation dataset
path_df_with_elevation='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2010_2018/excel_spatial_aggreation_and_other/' 
f_20102018 = open(path_df_with_elevation+'df_20102018_with_elevation_high_estimate_rignotetalregions', "rb")
df_2010_2018_elevation = pickle.load(f_20102018)
f_20102018.close()


#Plot data
plot_thickness_evolution(loc6,df_2010_2018_csv,df_2010_2018_elevation,ax1,ax2t,custom_angle=-120,offset_x=7000,offset_y=-18000,casestudy_nb='a')

plot_thickness_evolution(loc8,df_2010_2018_csv,df_2010_2018_elevation,ax1,ax3t,custom_angle=-90,offset_x=10000,offset_y=-5000,casestudy_nb='b')

plot_thickness_evolution(loc1,df_2010_2018_csv,df_2010_2018_elevation,ax1,ax4t,custom_angle=-52,offset_x=10000,offset_y=1000,casestudy_nb='c')

plot_thickness_evolution(loc9,df_2010_2018_csv,df_2010_2018_elevation,ax1,ax5t,custom_angle=-90,offset_x=10000,offset_y=-5000,casestudy_nb='d')

plot_thickness_evolution(loc3,df_2010_2018_csv,df_2010_2018_elevation,ax1,ax6t,custom_angle=-90,offset_x=10000,offset_y=-5000,casestudy_nb='e')

plot_thickness_evolution(loc2,df_2010_2018_csv,df_2010_2018_elevation,ax1,ax7t,custom_angle=-90,offset_x=10000,offset_y=-5000,casestudy_nb='f')


###################### From Tedstone et al., 2022 #####################
#from plot_map_decadal_change.py
# x0, x1, y0, y1
ax1.set_extent([-580000, -44000, -2650000, -1290000], crs=crs)
gl=ax1.gridlines(draw_labels=True, xlocs=[-50, -35], ylocs=[70, 75], x_inline=False, y_inline=False,linewidth=0.5)
#Customize lat labels
gl.ylabels_right = False
gl.xlabels_bottom = False
ax1.axis('off')
###################### From Tedstone et al., 2022 #####################

#panels a-b share axis, panels c-d-e-f share axis
ax2t.set_xlim(1130,1440)
ax3t.set_xlim(1130,1440)
ax4t.set_xlim(1600,2080)
ax5t.set_xlim(1600,2080)
ax6t.set_xlim(1600,2080)
ax7t.set_xlim(1600,2080)

#Display distance as Elevation [m]
ax7t.set_ylabel('Column ice thickness [m]')
ax7t.set_xlabel('Elevation [m]')

#Custom legend myself
legend_elements = [Patch(facecolor='#fdd49e',label='2010'),
                   Patch(facecolor='#fc8d59',label='2011-2012'),
                   Patch(facecolor='#d7301f',label='2013-2014'),
                   Patch(facecolor='#7f0000',label='2017-2018')]

ax_legend.legend(handles=legend_elements)
plt.legend()

#Get rid of axis in legend axis
ax_legend.axis('off')
ax_legend.set_title('Legend')
plt.show()
ax7t.legend_.remove()

#Save the figure
plt.savefig('C:/Users/jullienn/switchdrive/Private/research/RT1/figures/fig2/v2/fig2.png',dpi=500)


'''
#NW but requires additional coding + turning trace
plot_thickness_evolution(loc7,df_2010_2018_csv,df_2010_2018_elevation,ax1,custom_angle=-90,offset_x=10000,offset_y=-5000)
'''