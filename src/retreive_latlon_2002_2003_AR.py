# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 15:35:11 2020

@author: Nicolas Jullien
"""
import scipy.io
import rasterio
from matplotlib import pyplot
import numpy as np
import h5py
import matplotlib.colors as mcolors
import pandas as pd

# 1. Let's start with 2003 dataset
    # 1.a. Let's starts with data acquired on May 13th, 2003.
    
    #Load the master file of May 13, 2003
    may13_03_gpslatlontime= scipy.io.loadmat('C://Users//Nicolas Jullien//Documents//PhD//iceslabs_processing//iceslabs_MacFerrin//data//2003//may13//may13_03_gpslatlontime.mat')
    may13_03_gpslatlontime['seconds']
    may13_03_gpslatlontime['useconds']
    may13_03_gpslatlontime['lat_gps']
    may13_03_gpslatlontime['lon_gps']
    may13_03_gpslatlontime['time_gps']
    
    #Create the dataframe        
    df_gps_may13_03=pd.DataFrame({'seconds':pd.Series(np.ndarray.flatten(np.transpose(may13_03_gpslatlontime['seconds']))),
    'useconds':pd.Series(np.ndarray.flatten(np.transpose(may13_03_gpslatlontime['useconds']))),
    'lat_gps':pd.Series(np.ndarray.flatten(np.transpose(may13_03_gpslatlontime['lat_gps']))),
    'lon_gps':pd.Series(np.ndarray.flatten(np.transpose(may13_03_gpslatlontime['lon_gps']))),
    'time_gps':pd.Series(np.ndarray.flatten(np.transpose(may13_03_gpslatlontime['time_gps'])))})
    
    #Set the seconds column to be the index of the dataframe
    df_gps_may13_03=df_gps_may13_03.set_index('seconds')
        
    #Load an individual file
    may13_03_0= scipy.io.loadmat('C://Users//Nicolas Jullien//Documents//PhD//iceslabs_processing//iceslabs_MacFerrin//data//2003//may13//may13_03_0.mat')
    
    #We have to round the timearr to the lower integer
    int_may13_03_0 = np.floor(may13_03_0['timearr'])
    #Transform decimals into integer
    int_may13_03_0=int_may13_03_0.astype(int)
    
    #Create the dataframe
    df_may13_03_0_timearr=pd.DataFrame({'timearr':pd.Series(np.ndarray.flatten(np.transpose(int_may13_03_0))),
                                       'index_vector':pd.Series(np.arange(0,int_may13_03_0.size,1)),
                                       'timearr_trace':pd.Series(np.ndarray.flatten(np.transpose(int_may13_03_0)))})
 
    #Set the timearr column to be the index of the dataframe
    df_may13_03_0_timearr=df_may13_03_0_timearr.set_index('timearr')
    
    #Rename the column 'timearr_trace' to 'timearr'
    df_may13_03_0_timearr.columns = ['index_vector', 'timearr']
    #Make the correspondance between timearr and seconds and join datasets
    result_join=df_may13_03_0_timearr.join(df_gps_may13_03, lsuffix='_time_arr', rsuffix='_gps_latlontime')
    
    # Be careful about the remove duplicates with respect to the time_gps. This
    # is the only solution I found to indeed remove duplictes but maybe check
    # on other datasets to be sure that it does the right thing!
    result_join_without_duplicates=result_join.drop_duplicates(subset=['time_gps'])
    
    #Store everything into one dictionnary (matrix and vectors of data)
    
    may13_03_trace0 = { "trace_id" : 'may13_03_trace0',
         "radar_echogram" : may13_03_0['filtfin'],
         "latlontime" : result_join_without_duplicates }
    
    
#load individual data
may13_03_0= scipy.io.loadmat('C://Users//Nicolas Jullien//Documents//PhD//iceslabs_processing//iceslabs_MacFerrin//data//2003//may13//may13_03_0.mat')
may13_03_1= scipy.io.loadmat('C://Users//Nicolas Jullien//Documents//PhD//iceslabs_processing//iceslabs_MacFerrin//data//2003//may13//may13_03_1.mat')
may13_03_2= scipy.io.loadmat('C://Users//Nicolas Jullien//Documents//PhD//iceslabs_processing//iceslabs_MacFerrin//data//2003//may13//may13_03_2.mat')
may13_03_3= scipy.io.loadmat('C://Users//Nicolas Jullien//Documents//PhD//iceslabs_processing//iceslabs_MacFerrin//data//2003//may13//may13_03_3.mat')
may13_03_4= scipy.io.loadmat('C://Users//Nicolas Jullien//Documents//PhD//iceslabs_processing//iceslabs_MacFerrin//data//2003//may13//may13_03_4.mat')
may13_03_5= scipy.io.loadmat('C://Users//Nicolas Jullien//Documents//PhD//iceslabs_processing//iceslabs_MacFerrin//data//2003//may13//may13_03_5.mat')
may13_03_6= scipy.io.loadmat('C://Users//Nicolas Jullien//Documents//PhD//iceslabs_processing//iceslabs_MacFerrin//data//2003//may13//may13_03_6.mat')
may13_03_7= scipy.io.loadmat('C://Users//Nicolas Jullien//Documents//PhD//iceslabs_processing//iceslabs_MacFerrin//data//2003//may13//may13_03_7.mat')
may13_03_8= scipy.io.loadmat('C://Users//Nicolas Jullien//Documents//PhD//iceslabs_processing//iceslabs_MacFerrin//data//2003//may13//may13_03_8.mat')
may13_03_9= scipy.io.loadmat('C://Users//Nicolas Jullien//Documents//PhD//iceslabs_processing//iceslabs_MacFerrin//data//2003//may13//may13_03_9.mat')
may13_03_10= scipy.io.loadmat('C://Users//Nicolas Jullien//Documents//PhD//iceslabs_processing//iceslabs_MacFerrin//data//2003//may13//may13_03_10.mat')
may13_03_11= scipy.io.loadmat('C://Users//Nicolas Jullien//Documents//PhD//iceslabs_processing//iceslabs_MacFerrin//data//2003//may13//may13_03_11.mat')
may13_03_12= scipy.io.loadmat('C://Users//Nicolas Jullien//Documents//PhD//iceslabs_processing//iceslabs_MacFerrin//data//2003//may13//may13_03_12.mat')
may13_03_13= scipy.io.loadmat('C://Users//Nicolas Jullien//Documents//PhD//iceslabs_processing//iceslabs_MacFerrin//data//2003//may13//may13_03_13.mat')
may13_03_14= scipy.io.loadmat('C://Users//Nicolas Jullien//Documents//PhD//iceslabs_processing//iceslabs_MacFerrin//data//2003//may13//may13_03_14.mat')
may13_03_15= scipy.io.loadmat('C://Users//Nicolas Jullien//Documents//PhD//iceslabs_processing//iceslabs_MacFerrin//data//2003//may13//may13_03_15.mat')
may13_03_16= scipy.io.loadmat('C://Users//Nicolas Jullien//Documents//PhD//iceslabs_processing//iceslabs_MacFerrin//data//2003//may13//may13_03_16.mat')
may13_03_17= scipy.io.loadmat('C://Users//Nicolas Jullien//Documents//PhD//iceslabs_processing//iceslabs_MacFerrin//data//2003//may13//may13_03_17.mat')
may13_03_18= scipy.io.loadmat('C://Users//Nicolas Jullien//Documents//PhD//iceslabs_processing//iceslabs_MacFerrin//data//2003//may13//may13_03_18.mat')
may13_03_19= scipy.io.loadmat('C://Users//Nicolas Jullien//Documents//PhD//iceslabs_processing//iceslabs_MacFerrin//data//2003//may13//may13_03_19.mat')
may13_03_20= scipy.io.loadmat('C://Users//Nicolas Jullien//Documents//PhD//iceslabs_processing//iceslabs_MacFerrin//data//2003//may13//may13_03_20.mat')
may13_03_21= scipy.io.loadmat('C://Users//Nicolas Jullien//Documents//PhD//iceslabs_processing//iceslabs_MacFerrin//data//2003//may13//may13_03_21.mat')
may13_03_22= scipy.io.loadmat('C://Users//Nicolas Jullien//Documents//PhD//iceslabs_processing//iceslabs_MacFerrin//data//2003//may13//may13_03_22.mat')
may13_03_23= scipy.io.loadmat('C://Users//Nicolas Jullien//Documents//PhD//iceslabs_processing//iceslabs_MacFerrin//data//2003//may13//may13_03_23.mat')
may13_03_24= scipy.io.loadmat('C://Users//Nicolas Jullien//Documents//PhD//iceslabs_processing//iceslabs_MacFerrin//data//2003//may13//may13_03_24.mat')
may13_03_25= scipy.io.loadmat('C://Users//Nicolas Jullien//Documents//PhD//iceslabs_processing//iceslabs_MacFerrin//data//2003//may13//may13_03_25.mat')
may13_03_26= scipy.io.loadmat('C://Users//Nicolas Jullien//Documents//PhD//iceslabs_processing//iceslabs_MacFerrin//data//2003//may13//may13_03_26.mat')
may13_03_27= scipy.io.loadmat('C://Users//Nicolas Jullien//Documents//PhD//iceslabs_processing//iceslabs_MacFerrin//data//2003//may13//may13_03_27.mat')
may13_03_28= scipy.io.loadmat('C://Users//Nicolas Jullien//Documents//PhD//iceslabs_processing//iceslabs_MacFerrin//data//2003//may13//may13_03_28.mat')
may13_03_29= scipy.io.loadmat('C://Users//Nicolas Jullien//Documents//PhD//iceslabs_processing//iceslabs_MacFerrin//data//2003//may13//may13_03_29.mat')
may13_03_30= scipy.io.loadmat('C://Users//Nicolas Jullien//Documents//PhD//iceslabs_processing//iceslabs_MacFerrin//data//2003//may13//may13_03_30.mat')

#store the filtfin to obtain the length of the depth
may13_03_0_filtfin= may13_03_0['filtfin']
may13_03_1_filtfin= may13_03_1['filtfin']
may13_03_2_filtfin= may13_03_2['filtfin']
may13_03_3_filtfin= may13_03_3['filtfin']
may13_03_4_filtfin= may13_03_4['filtfin']
may13_03_5_filtfin= may13_03_5['filtfin']
may13_03_6_filtfin= may13_03_6['filtfin']
may13_03_7_filtfin= may13_03_7['filtfin']
may13_03_8_filtfin= may13_03_8['filtfin']
may13_03_9_filtfin= may13_03_9['filtfin']
may13_03_10_filtfin= may13_03_10['filtfin']
may13_03_11_filtfin= may13_03_11['filtfin']
may13_03_12_filtfin= may13_03_12['filtfin']
may13_03_13_filtfin= may13_03_13['filtfin']
may13_03_14_filtfin= may13_03_14['filtfin']
may13_03_15_filtfin= may13_03_15['filtfin']
may13_03_16_filtfin= may13_03_16['filtfin']
may13_03_17_filtfin= may13_03_17['filtfin']
may13_03_18_filtfin= may13_03_18['filtfin']
may13_03_19_filtfin= may13_03_19['filtfin']
may13_03_20_filtfin= may13_03_20['filtfin']
may13_03_21_filtfin= may13_03_21['filtfin']
may13_03_22_filtfin= may13_03_22['filtfin']
may13_03_23_filtfin= may13_03_23['filtfin']
may13_03_24_filtfin= may13_03_24['filtfin']
may13_03_25_filtfin= may13_03_25['filtfin']
may13_03_26_filtfin= may13_03_26['filtfin']
may13_03_27_filtfin= may13_03_27['filtfin']
may13_03_28_filtfin= may13_03_28['filtfin']
may13_03_29_filtfin= may13_03_29['filtfin']
may13_03_30_filtfin= may13_03_30['filtfin']

