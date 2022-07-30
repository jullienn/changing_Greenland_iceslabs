# -*- coding: utf-8 -*-
"""
Created on Fri Dec 10 11:05:30 2021

@author: JullienN
"""

import scipy.io
import numpy as np
import pandas as pd
from os import listdir
from os.path import isfile, join
import pdb
import os.path
import os

year_to_download='2010'

#Create empty dataframe
coordinates_2010=pd.DataFrame(columns=['LAT','LON','FRAME'])

########################## Download 2010-2014 AR data #########################
#This is from Download_AR_data.py
#Code from: https://gist.github.com/nasrulhazim/cfd5f01e3b261b09d54f721cc1a7c50d
if (year_to_download=='2010'):
    
    from ftplib import FTP
    from datetime import datetime
    
    print('Initialisation ...')
    #Set data we want to download
    download_mat='TRUE'
    
    #Load the data that have already been downloaded
    #f = open('C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/Exclusion_folder/datetrack_20102018.txt','r')
    f = open('/flash/jullienn/data/threshold_processing/datetrack_20102018.txt','r')
    lines = [line.strip() for line in f.readlines() if len(line.strip()) > 0]
    f.close()
        
    data_already_downloaded=[]
    for indiv_trace in lines:
        
        if (not(indiv_trace[0:4] == '2010')):
            print('Not 2010, continue')
            continue
        
        #Create a vector from begin_trace_nb to end_trace_nb
        vect_traces_nb=np.arange(int(indiv_trace[12:15]),int(indiv_trace[16:19])+1,1)
        
        #Create the right suffixe for imediate check later on
        vect_str=[]
        for trace_nb_indiv in list(vect_traces_nb):
            if (trace_nb_indiv<10):
                trace_nb_indiv_str='00'+str(trace_nb_indiv)
            elif ((trace_nb_indiv>=10) and (trace_nb_indiv<100)):
                trace_nb_indiv_str='0'+str(trace_nb_indiv)
            else:
                trace_nb_indiv_str=str(trace_nb_indiv)
            #Store str trace nb                        
            vect_str=np.append(vect_str,trace_nb_indiv_str)
            #Store the complete filenumber
            data_already_downloaded=np.append(data_already_downloaded,indiv_trace[0:11]+'_'+trace_nb_indiv_str)
    
    print('Logging on the cresis website ...')
    
    start = datetime.now()
    ftp = FTP('data.cresis.ku.edu')
    ftp.login()
    
    path='/data/accum/'
    ftp.cwd(path)
    
    # Get folders_years name
    folders_years = ftp.nlst()
    
    #folders I want to download
    list_download=['2010_Greenland_P3']
    
    print('Starting to download ...')
    for folder_year in folders_years:
        
        if (folder_year in list(list_download)):
            print('    Downloading ' + folder_year)
            
            #Go to folder year
            folder_year_name=[]
                
            if(download_mat=='TRUE'):
                folder_year_name=path + folder_year + '/CSARP_qlook/'
                print('    Downloading the .mat files ...')
            
            #Go to folder CSARP_standard
            ftp.cwd(folder_year_name)
            
            # For this particular year, get folders name
            folders=[]
            folders = ftp.nlst()
            #Go to the folder
            
            #Loop over the folders, and download all the data, except the ones already downloaded
            for folder in folders:
                
                print('        '+folder+' folder')
                folder_name=[]
                folder_name=folder_year_name + folder + '/'
                ftp.cwd(folder_name)
            
                # Get all Files in this folder
                files=[]
                files = ftp.nlst()
                
                if (download_mat=='TRUE'):
                    #Define the paths
                    #path_data_there='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/' + folder_year + '/' + 'CSARP_qlook/' + folder
                    #path_to_save='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/' + '2010_Greenland_P3_flightlines/' + 'CSARP_qlook/' + folder
                    
                    path_data_there='/flash/jullienn/data/threshold_processing/' + folder_year + '/' + 'CSARP_qlook/' + folder
                    path_to_save='/flash/jullienn/flightlines/data/2010_Greenland_P3_flightlines/' + 'CSARP_qlook/' + folder
                    
                    #Create the directory to store the data
                    #this is from: https://thispointer.com/how-to-create-a-directory
                    #-in-python/#:~:text=Python%27s%20OS%20module%20provides%20an%20
                    #another%20function%20to%20create%20a%20directories%20i.e.&text=.
                    #makedirs(path)-,os.,mkdir%20%2Dp%20command%20in%20linux.
                    os.mkdir(path_to_save)
                    print('            Created the directory'+path_to_save)
                    
                    for file in files:
                        
                        #Create empty temporaire dataframe
                        temp_coordinates_2010=pd.DataFrame(columns=['LAT','LON','FRAME'])
                        if (file[5:20] in list(data_already_downloaded)):
                            print(file[5:20],'             already downloaded, load lat/lon')                            
                            #Load data
                            fdata_filename = scipy.io.loadmat(path_data_there+'/Data_'+file[5:20])
                        
                            #Create temporaire df
                            temp_coordinates_2010['LAT']=np.ndarray.flatten(fdata_filename['Latitude'])
                            temp_coordinates_2010['LON']=np.ndarray.flatten(fdata_filename['Longitude'])
                            temp_coordinates_2010['FRAME']=[file[5:20]]*fdata_filename['Latitude'].size
                            
                            #Store data
                            coordinates_2010=coordinates_2010.append(temp_coordinates_2010)
                        else:
                            #Download .mat files
                            print("                Downloading... " + file)
                            ftp.retrbinary("RETR " + file ,open(path_to_save + '/' + file, 'wb').write)
                            
                            #Load data
                            fdata_filename = scipy.io.loadmat(path_to_save+'/'+file)
                        
                            #Create temporaire df
                            temp_coordinates_2010['LAT']=np.ndarray.flatten(fdata_filename['Latitude'])
                            temp_coordinates_2010['LON']=np.ndarray.flatten(fdata_filename['Longitude'])
                            temp_coordinates_2010['FRAME']=[file[5:20]]*fdata_filename['Latitude'].size
                            
                            #Store data
                            coordinates_2010=coordinates_2010.append(temp_coordinates_2010)
                            
                            #Delete the file
                            #os.remove(path_to_save+'/Data_'+file[5:20]+'.mat')
                        
        else:
            print('Not 2010, continue')
            continue
            
    ftp.close()
    end = datetime.now()
    diff = end - start
    print('Aggregation of 2010 coordinates done for ' + str(diff.seconds) + 's')

#Save the lat/lon dataframe as csv file
#coordinates_2010.to_csv('C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/flightlines/2010_Greenland_P3.csv')
coordinates_2010.to_csv('/flash/jullienn/flightlines/data/2010_Greenland_P3.csv')

print('--- End of processing ---')
########################## Download 2010-2014 AR data #########################