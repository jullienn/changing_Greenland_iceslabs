# -*- coding: utf-8 -*-
"""
Created on Sun Dec 19 12:14:06 2021

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

def compute_distances_pyproj(eastings,northings):
    #This is from https://pyproj4.github.io/pyproj/stable/api/geod.html#pyproj.Geod.line_length
    from pyproj import Geod
    geod = Geod(ellps="WGS84")
    
    distances=[]
    pdb.set_trace()
    
    for line_length in geod.line_lengths(northings, eastings):
        distances=np.append(distances,line_length)
        
    return return_cumsum_distances

def plot_thickness_evolution(dictionnary_case_study,df_2010_2018_csv,df_2010_2018_elevation,ax1,axt,custom_angle,offset_x,offset_y,casestudy_nb):
    
    #Define empty dictionnary for elevation slice definition
    df_for_elev=pd.DataFrame(columns=list(df_2010_2018_elevation.keys()))
    
    #Loop over the years
    for year in dictionnary_case_study.keys():
        if (dictionnary_case_study[year] == 'empty'):
            continue  
        #Select data for the trace
        df_for_elev_temp=df_2010_2018_elevation[df_2010_2018_elevation['Track_name']==dictionnary_case_study[year][0][5:20]+'_'+dictionnary_case_study[year][-1][17:20]]
        
        #If panel, b, then start only at lon=-530298
        if (casestudy_nb=='b'):
            #Do not keep where lon_3413 < -530298 because not monotoneously elevation increase
            df_for_elev_temp=df_for_elev_temp[df_for_elev_temp['lon_3413']>=-530298]
        
        #Append data to each other
        df_for_elev=df_for_elev.append(df_for_elev_temp)
                
        #Display data
        ax1.scatter(df_for_elev_temp['lon_3413'],
                    df_for_elev_temp['lat_3413'],
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
    
    #pdb.set_trace()
    #Add number of case study on fig localisation    
    axt.text(0.99, 0.825,casestudy_nb, ha='center', va='center', transform=axt.transAxes,fontsize=15)#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot
    
    #Define palette for time periods
    #This is from https://www.python-graph-gallery.com/33-control-colors-of-boxplot-seaborn
    my_pal = {'2010': "#fdd49e", '2011-2012': "#fc8d59", '2013-2014':"#d7301f",'2017-2018':"#7f0000"}
        
    #Create an empty df_sampling
    df_sampling=pd.DataFrame(columns=['Track_name','time_period','low_bound', 'high_bound', 'bound_nb', 'mean', 'median', 'q025', 'q075','stddev','rolling_10_median_scatter'])
        
    #Sort df_for_elev from low to high longitude (from west to east)
    df_for_elev_sorted=df_for_elev.sort_values(by=['lon_3413'])
    
    #Create a nan array for storing distances
    df_for_elev_sorted['distances']=np.nan
    
    #Store coordinates of the bounds of the transect
    bounds_transect=np.array([[df_for_elev_sorted.iloc[0]['lon_3413'], df_for_elev_sorted.iloc[0]['lat_3413']],
                              [df_for_elev_sorted.iloc[-1]['lon_3413'], df_for_elev_sorted.iloc[-1]['lat_3413']]])
    
    #Compute distance between the westernmost and easternmost point
    bounds_distances=compute_distances(bounds_transect[:,0],bounds_transect[:,1])
        
    #grizou= compute_distances_pyproj(np.asarray(df_for_elev_sorted['lon']),np.asarray(df_for_elev_sorted['lat']))
    
    #Distance divide every 300m.
    dist_bin_desired=300
    dist_divide=np.arange(np.floor(np.min(bounds_distances)),np.floor(np.max(bounds_distances))+1+dist_bin_desired,dist_bin_desired)
        
    #Define window size for smoothing
    winsize=3
        
    #Calculate distance for every single year
    for indiv_year in np.array([2010,2011,2012,2013,2014,2017,2018]):
        #Extract the indexes of the corresponding year
        ind_indiv_year=np.where(df_for_elev_sorted['year']==indiv_year)
        #Select the corresponding year
        df_trace_year_sorted_for_dist=df_for_elev_sorted.iloc[ind_indiv_year]
        
        if (len(df_trace_year_sorted_for_dist)>0):            
            #Calculate the distance compared to start of transect
            #a. Add the start of the transect for calculating distances
            coordinates_df=[np.append(bounds_transect[0,0],np.asarray(df_trace_year_sorted_for_dist['lon_3413'])),
                            np.append(bounds_transect[0,1],np.asarray(df_trace_year_sorted_for_dist['lat_3413']))]
            #b. Calculate the distances
            distances_with_start_transect=compute_distances(coordinates_df[0],coordinates_df[1])
            #c. Store the distances
            df_for_elev_sorted['distances'].iloc[ind_indiv_year]=distances_with_start_transect[1:]
        
    #Define empty list
    app_time_period=[]
    app_low_bound=[]
    app_high_bound=[]
    app_bound_nb=[]
    #app_Track_name=[]
    app_mean=[]
    app_median=[]
    app_q025=[]
    app_q075=[]
    app_stddev=[]
    
    #Loop over the different time periods (2010, 2011-2012, 2013-2014, 2017-2018)
    for time_period in list(['2010','2011-2012','2013-2014','2017-2018']):
        
        #Get data for that specific time period
        if (time_period == '2010'):
            df_trace_year_sorted=df_for_elev_sorted[df_for_elev_sorted['year']==2010]
        elif (time_period == '2011-2012'):
            df_trace_year_sorted=df_for_elev_sorted[(df_for_elev_sorted['year']>=2011) & (df_for_elev_sorted['year']<=2012)]
        elif (time_period == '2013-2014'):
            df_trace_year_sorted=df_for_elev_sorted[(df_for_elev_sorted['year']>=2013) & (df_for_elev_sorted['year']<=2014)]
        elif (time_period == '2017-2018'):
            df_trace_year_sorted=df_for_elev_sorted[(df_for_elev_sorted['year']>=2017) & (df_for_elev_sorted['year']<=2018)]
        else:
            print('Time period not known, break')
            break
          
        if (len(df_trace_year_sorted)==0):
            #No data in this time period, continue
            continue
        else:
            #Set bound_nb to 0
            bound_nb=0
            #Loop over the dist divide
            for i in range(1,len(dist_divide)):
                
                #Identify low and higher end of the slice
                low_bound=dist_divide[i-1]
                high_bound=dist_divide[i]
        
                #Select all the data belonging to this elev slice
                ind_slice=np.logical_and(np.array(df_trace_year_sorted['distances']>=low_bound),np.array(df_trace_year_sorted['distances']<high_bound))
                df_select=df_trace_year_sorted[ind_slice]
                
                #Append data to each other - general info                
                app_time_period=np.append(app_time_period,np.asarray(time_period))
                app_low_bound=np.append(app_low_bound,low_bound)
                app_high_bound=np.append(app_high_bound,high_bound)
                app_bound_nb=np.append(app_bound_nb,str(bound_nb))
                
                if (len(df_select)==0):
                    #Append data to each other - data
                    #app_Track_name=np.append(app_Track_name,np.nan)
                    app_mean=np.append(app_mean,np.nan)
                    app_median=np.append(app_median,np.nan)
                    app_q025=np.append(app_q025,np.nan)
                    app_q075=np.append(app_q075,np.nan)
                    app_stddev=np.append(app_stddev,np.nan)
                else:
                    #Append data to each other -data
                    #app_Track_name=np.append(app_Track_name,np.asarray(df_select['Track_name'].unique()))
                    app_mean=np.append(app_mean,df_select["20m_ice_content_m"].mean())
                    app_median=np.append(app_median,df_select["20m_ice_content_m"].median())
                    app_q025=np.append(app_q025,df_select["20m_ice_content_m"].quantile(q=0.25))
                    app_q075=np.append(app_q075,df_select["20m_ice_content_m"].quantile(q=0.75))
                    app_stddev=np.append(app_stddev,df_select["20m_ice_content_m"].std())

                #Update bound_nb
                bound_nb=bound_nb+1
    
    #Create df_sampling which is the dataframe reuniting all the appended lists
    df_sampling = pd.DataFrame(data={'time_period': app_time_period})
    df_sampling['low_bound']=app_low_bound
    df_sampling['high_bound']=app_high_bound
    df_sampling['bound_nb']=app_bound_nb
    #df_sampling['Track_name']=app_Track_name
    df_sampling['mean']=app_mean
    df_sampling['median']=app_median
    df_sampling['q025']=app_q025
    df_sampling['q075']=app_q075
    df_sampling['stddev']=app_stddev
    df_sampling['rolling_10_median_scatter']=[np.nan]*len(app_mean)
    
    for time_period in list(['2010','2011-2012','2013-2014','2017-2018']):
        if (len(df_sampling[df_sampling['time_period']==time_period])==0):
            #Empty time period, continue
            continue
        else:
            df_plot=df_sampling[df_sampling['time_period']==time_period]
            
            #Rolling window, size = winsize
            df_plot['rolling_10_median']=df_plot.rolling(winsize, win_type=None,center=True).quantile(quantile=0.5)['median'] 
            df_plot['rolling_10_q025']=df_plot.rolling(winsize, win_type=None,center=True).quantile(quantile=0.5)['q025'] 
            df_plot['rolling_10_q075']=df_plot.rolling(winsize, win_type=None,center=True).quantile(quantile=0.5)['q075'] 
            
            #Where window rolling shows NaN because of NaN surrounding, add raw data
            for i in range(0,len(df_plot)):
                if (np.isnan(np.asarray(df_plot['rolling_10_median'].iloc[i]))):
                    df_plot['rolling_10_median_scatter'].iloc[i]=df_plot['median'].iloc[i]
                    df_plot['rolling_10_q025'].iloc[i]=np.nan#df_plot['q025'].iloc[i]
                    df_plot['rolling_10_q075'].iloc[i]=np.nan#df_plot['q075'].iloc[i]
                
                if (i>0):
                    if (np.isnan(np.asarray(df_plot['rolling_10_median'].iloc[i-1])) and not(np.isnan(np.asarray(df_plot['rolling_10_median'].iloc[i])))):
                        df_plot['rolling_10_median_scatter'].iloc[i]=df_plot['rolling_10_median'].iloc[i]
                        
                if (i<(len(df_plot)-1)):
                    if (np.isnan(np.asarray(df_plot['rolling_10_median'].iloc[i+1])) and not(np.isnan(np.asarray(df_plot['rolling_10_median'].iloc[i])))):
                        df_plot['rolling_10_median_scatter'].iloc[i]=df_plot['rolling_10_median'].iloc[i]
                
            # Plot the median
            axt.plot(df_plot["low_bound"],df_plot["rolling_10_median"],color=my_pal[time_period])
            #Display IQR
            axt.fill_between(df_plot['low_bound'], df_plot['rolling_10_q025'], df_plot['rolling_10_q075'], alpha=0.3,color=my_pal[time_period])
            #Display raw data where moving window do not display
            axt.plot(df_plot['low_bound'], df_plot['rolling_10_median_scatter'],color=my_pal[time_period],alpha=0.5)
            axt.scatter(df_plot['low_bound'], df_plot['rolling_10_median_scatter'],c=my_pal[time_period],marker='.',s=0.5,alpha=1)
            
    #Get rid of legend
    #axt.legend_.remove()
    axt.set_xlabel('')
    axt.set_ylabel('')
    
    #Activate ticks x and y label
    axt.yaxis.tick_left()
    axt.xaxis.tick_bottom()
    
    #Place yticks on the right hand side
    axt.yaxis.tick_right()
        
    #4. Display elevation
    #Store the xticks for the distance
    xtick_distance=axt.get_xticks()
    #Set the xticks
    axt.set_xticks(xtick_distance)
    
    #Find closest corresponding elevation
    #This is from https://stackoverflow.com/questions/11244514/modify-tick-label-text
    elevation_display=[np.nan]*len(xtick_distance)
    count=0
    for indiv_dist in xtick_distance:
        if (indiv_dist<0):
            elevation_display[count]=''
        else:
            #Extract index where distance is minimal
            index_closest=np.argmin(np.abs(np.abs(df_for_elev_sorted['distances'])-np.abs(indiv_dist)))
            #If minimum distance is higher than 1km, store nan. If not, store corresponding elevation
            if (np.abs(np.abs(df_for_elev_sorted['distances'])-np.abs(indiv_dist)).iloc[index_closest] > 1000):
                elevation_display[count]=''
            else:
                elevation_display[count]=np.round(df_for_elev_sorted.iloc[index_closest]['elevation']).astype(int)
            
        count=count+1
        
    #Display elevation on the top xticklabels
    #This is from https://stackoverflow.com/questions/19884335/matplotlib-top-bottom-ticks-different "Zaus' reply"
    ax_t = axt.secondary_xaxis('top')
    ax_t.set_xticks(xtick_distance)
    ax_t.set_xticklabels(elevation_display)
    
    #Display bottom xtick in km instead of m
    axt.set_xticklabels((xtick_distance/1000).astype(int))
    
    #Modify spacing between xticklabels and xticks
    axt.tick_params(pad=1.2)
    ax_t.tick_params(pad=1.2)

    #Set xlims
    #axt.set_xlim(0,70000)
    
    '''
    # Hide grid lines, from https://stackoverflow.com/questions/45148704/how-to-hide-axes-and-gridlines-in-matplotlib-python
    axt.grid(False)
    '''
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
import matplotlib.ticker as mticker


#Set fontsize plot
plt.rcParams.update({'font.size': 10})

### -------------------------- Load shapefiles --------------------------- ###
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

path='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2002_2018/final_excel/high_estimate/'

#Load all 2010-2018 data without spatial aggregation
df_2010_2018_csv = pd.read_csv(path+'Ice_Layer_Output_Thicknesses_2010_2018_jullienetal2021_high_estimate.csv',delimiter=',',decimal='.')
#Transform the coordinated from WGS84 to EPSG:3413
transformer = Transformer.from_crs("EPSG:4326", "EPSG:3413", always_xy=True)
points=transformer.transform(np.asarray(df_2010_2018_csv["lon"]),np.asarray(df_2010_2018_csv["lat"]))

#Store lat/lon in 3413
df_2010_2018_csv['lon_3413']=points[0]
df_2010_2018_csv['lat_3413']=points[1]

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


#in 2017 overestimation of ice content
loc4={2010:['Data_20100512_04_073.mat','Data_20100512_04_074.mat'],
      2011:'empty',
      2012:'empty',
      2013:'empty',
      2014:'empty',
      2017:['Data_20170421_01_171.mat','Data_20170421_01_172.mat','Data_20170421_01_173.mat','Data_20170421_01_174.mat'],
      2018:['Data_20180425_01_166.mat','Data_20180425_01_167.mat','Data_20180425_01_168.mat','Data_20180425_01_169.mat']}

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

'''
in the NO
loc10={2010:'empty',
       2011:'empty',
       2012:'empty',
       2013:['Data_20130420_08_045.mat','Data_20130420_08_046.mat','Data_20130420_08_047.mat','Data_20130420_08_048.mat'],
       2014:['Data_20140519_08_066.mat','Data_20140519_08_067.mat','Data_20140519_08_068.mat','Data_20140519_08_069.mat'],
       2017:['Data_20170413_01_126.mat','Data_20170413_01_127.mat','Data_20170413_01_128.mat',
             'Data_20170413_01_129.mat','Data_20170413_01_130.mat','Data_20170413_01_131.mat',
             'Data_20170413_01_132.mat','Data_20170413_01_133.mat','Data_20170413_01_134.mat'],
       2018:'empty'}

loc11={2010:'empty',
       2011:'empty',
       2012:['Data_20120511_01_059.mat'],
       2013:'empty',
       2014:'empty',
       2017:['Data_20170417_01_104.mat','Data_20170417_01_105.mat','Data_20170417_01_106.mat'],
       2018:'empty'}

'''


pkpas={2010:['Data_20100512_04_073.mat','Data_20100512_04_074.mat'],
                    2011:'empty',
                    2012:'empty',
                    2013:'empty',
                    2014:'empty',
                    2017:'empty',
                    2018:['Data_20180425_01_166.mat','Data_20180425_01_167.mat','Data_20180425_01_168.mat','Data_20180425_01_169.mat']}

pkpas={2010:'empty',
                    2011:'empty',
                    2012:'empty',
                    2013:['Data_20130405_01_011.mat','Data_20130405_01_012.mat','Data_20130405_01_013.mat'],
                    2014:['Data_20140424_03_046.mat','Data_20140424_03_047.mat','Data_20140424_03_048.mat'],
                    2017:['Data_20170422_01_007.mat','Data_20170422_01_008.mat','Data_20170422_01_009.mat'],
                    2018:'empty'}


pkpas={2010:['Data_20100512_04_073.mat','Data_20100512_04_074.mat'],
                    2011:'empty',
                    2012:'empty',
                    2013:'empty',
                    2014:'empty',
                    2017:['Data_20170421_01_171.mat','Data_20170421_01_172.mat','Data_20170421_01_173.mat','Data_20170421_01_174.mat'],
                    2018:['Data_20180425_01_166.mat','Data_20180425_01_167.mat','Data_20180425_01_168.mat','Data_20180425_01_169.mat']}

plutotbien={2010:'empty',
            2011:'empty',
            2012:['Data_20120418_01_005.mat','Data_20120418_01_006.mat','Data_20120418_01_007.mat'],
            2013:'empty',
            2014:'empty',
            2017:['Data_20170505_02_181.mat','Data_20170505_02_182.mat','Data_20170505_02_183.mat'],
            2018:'empty'}



###################### From Tedstone et al., 2022 #####################
#from plot_map_decadal_change.py
# Define the CartoPy CRS object.
crs = ccrs.NorthPolarStereo(central_longitude=-45., true_scale_latitude=70.)

# This can be converted into a `proj4` string/dict compatible with GeoPandas
crs_proj4 = crs.proj4_init
###################### From Tedstone et al., 2022 #####################
        
fig = plt.figure(figsize=(32,48))
gs = gridspec.GridSpec(39, 20)
gs.update(wspace=0.001)
#gs.update(wspace=0.001)
#projection set up from https://stackoverflow.com/questions/33942233/how-do-i-change-matplotlibs-subplot-projection-of-an-existing-axis
ax1 = plt.subplot(gs[0:25, 0:2],projection=crs)
ax_legend = plt.subplot(gs[34:39, 0:2])

ax2t = plt.subplot(gs[0:4, 3:20])
ax3t = plt.subplot(gs[7:11, 3:20])
ax4t = plt.subplot(gs[14:18, 3:20])
ax5t = plt.subplot(gs[21:25, 3:20])
ax6t = plt.subplot(gs[28:32, 3:20])
ax7t = plt.subplot(gs[35:39, 3:20])

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
path_df_with_elevation='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2002_2018/final_excel/high_estimate/' 
f_20102018 = open(path_df_with_elevation+'df_20102018_with_elevation_high_estimate_rignotetalregions', "rb")
df_2010_2018_elevation = pickle.load(f_20102018)
f_20102018.close()

#Plot data
plot_thickness_evolution(loc6,df_2010_2018_csv,df_2010_2018_elevation,ax1,ax2t,custom_angle=-120,offset_x=7000,offset_y=-18000,casestudy_nb='a')

plot_thickness_evolution(loc8,df_2010_2018_csv,df_2010_2018_elevation,ax1,ax3t,custom_angle=-90,offset_x=10000,offset_y=-5000,casestudy_nb='b')
#previousl b was loc8
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

#Display distance as Elevation [m]
ax5t.yaxis.set_label_position("right")
ax5t.set_ylabel('Column ice thickness [m]')
ax7t.set_xlabel('Distance [km]')
ax2t.xaxis.set_label_position("top")
ax2t.set_xlabel('Elevation [m]')

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

pdb.set_trace()
#Save the figure
plt.savefig('C:/Users/jullienn/switchdrive/Private/research/RT1/figures/fig2/v6/fig2.png',dpi=500)


'''
#NW but requires additional coding + turning trace
plot_thickness_evolution(loc7,df_2010_2018_csv,df_2010_2018_elevation,ax1,custom_angle=-90,offset_x=10000,offset_y=-5000)
'''