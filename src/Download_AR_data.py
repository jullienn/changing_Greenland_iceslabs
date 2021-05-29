# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 10:10:39 2020

@author: JullienN
"""

import scipy.io
from matplotlib import pyplot
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
import csv

year_to_download='2018'#'2010_2014'

############################ Download old AR data #############################
if (year_to_download=='old'):
    #Code from: https://gist.github.com/nasrulhazim/cfd5f01e3b261b09d54f721cc1a7c50d
    
    from ftplib import FTP
    from datetime import datetime
    
    
    start = datetime.now()
    ftp = FTP('data.cresis.ku.edu')
    ftp.login()
    
    #path='/data/accum/old_format/2003/'
    path='/data/accum/old_format/'
    ftp.cwd(path)
    
    # Get folders_years name
    folders_years = ftp.nlst()
    
    for folder_year in folders_years:
        
        if (folder_year == '2006'):
            print('Break because year 2006')
            break
        elif (folder_year == '2009'):
            print('Break because year 2009')
            break
        elif (folder_year == '2011'):
            print('Break because year 2011')
            break
        else:
            print('This is either the year 2002 or 2003')
            folder_year_name=[]
            folder_year_name=path + folder_year + '/'
            ftp.cwd(folder_year_name)
            
            # For this particular year, get folders name
            folders=[]
            folders = ftp.nlst()
            #Go to the folder
            
            for folder in folders:
                folder_name=[]
                folder_name=path + folder_year + '/' + folder + '/'
                ftp.cwd(folder_name)
                print("Now in folder ..." + folder_name)
            
                # Get All Files in this folder
                files=[]
                files = ftp.nlst()
                
                # Print out the files
                for file in files:
                    print("Downloading..." + file)
                    ftp.retrbinary("RETR " + file ,open("D:/OIB/AR/" + folder_year + '/' + folder + "/" + file, 'wb').write)
                    #ftp.retrbinary("RETR " + file ,open("download/to/your/directory/" + file, 'wb').write)
    
    ftp.close()
    end = datetime.now()
    diff = end - start
    print('All files downloaded for ' + str(diff.seconds) + 's')
############################ Download old AR data #############################


########################## Download 2010-2014 AR data #########################
#Code from: https://gist.github.com/nasrulhazim/cfd5f01e3b261b09d54f721cc1a7c50d
if (year_to_download=='2010_2014'):
    
    from ftplib import FTP
    from datetime import datetime
    
    print('Initialisation ...')
    #Set data we want to download
    download_images='FALSE'
    download_mat='TRUE'
    
    #Load the data we have to download
    f = open('C:/Users/jullienn/Documents/working_environement/iceslabs_MacFerrin/data/2010_2014_data_download.txt','r')
    lines = [line.strip() for line in f.readlines() if len(line.strip()) > 0]
    f.close()
    
    #Extract the folder year month day nb
    df_yearmonthdaynb=pd.DataFrame({'complete_date':[None]*len(lines),
                                    'YMDnb':[None]*len(lines),
                                    'traces_nb':[None]*len(lines)})
    
    #Create the dictionnary
    yearmonthdaynb_dict = {k: {} for k in list(lines)}
    
    i=0
    all_data_to_download=[]
    for indiv_trace in lines:
        #Create the pd dataframe
        df_yearmonthdaynb['complete_date'][i]=indiv_trace
        df_yearmonthdaynb['YMDnb'][i]=indiv_trace[0:11]
        df_yearmonthdaynb['traces_nb'][i]=indiv_trace[12:19]
        i=i+1
        
        #Create nested dictionaries
        yearmonthdaynb_dict[indiv_trace]= {k: {} for k in list(['complete_date',
                                                               'YMDnb',
                                                               'traces_nb',
                                                               'vector_traces_nb'])}
        #Fill in the dictionnary
        yearmonthdaynb_dict[indiv_trace]['complete_date']=indiv_trace
        yearmonthdaynb_dict[indiv_trace]['YMDnb']=indiv_trace[0:11]
        yearmonthdaynb_dict[indiv_trace]['traces_nb']=indiv_trace[12:19]
        
        #List all the individual trace nb
        traces_nb=indiv_trace[12:19]
        begin_trace_nb=int(traces_nb[0:3])
        end_trace_nb=int(traces_nb[4:7])
        
        #Create a vector from begin_trace_nb to end_trace_nb
        vect_traces_nb=np.arange(begin_trace_nb,end_trace_nb+1,1)
        
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
            all_data_to_download=np.append(all_data_to_download,indiv_trace[0:11]+'_'+trace_nb_indiv_str)
        
        #Store this vector into the yearmonthdaynb_dict dictionnary
        yearmonthdaynb_dict[indiv_trace]['vector_traces_nb']=vect_str
    
    print('Logging on the cresis website ...')
    
    start = datetime.now()
    ftp = FTP('data.cresis.ku.edu')
    ftp.login()
    
    #path='/data/accum/old_format/2003/'
    path='/data/accum/'
    ftp.cwd(path)
    
    # Get folders_years name
    folders_years = ftp.nlst()
    
    #folders I want to download
    list_download=['2010_Greenland_P3','2011_Greenland_P3','2012_Greenland_P3',
                   '2013_Greenland_P3','2014_Greenland_P3']
    
    print('Starting to download ...')
    for folder_year in folders_years:
        
        if (folder_year in list(list_download)):
            print('    Downloading ' + folder_year)
            
            #Go to folder year
            folder_year_name=[]
            
            if (download_images=='TRUE'):
                folder_year_name=path + folder_year + '/images/'
                print('    Downloading the images ...')
                
            if(download_mat=='TRUE'):
                folder_year_name=path + folder_year + '/CSARP_qlook/'
                print('    Downloading the .mat files ...')
            
            #Go to folder CSARP_standard
            ftp.cwd(folder_year_name)
            
            # For this particular year, get folders name
            folders=[]
            folders = ftp.nlst()
            #Go to the folder
            
            #Loop over the folders, and download the data used by MacFerrin et al., 2019 in this folder
            for folder in folders:
                if (not(folder in list(df_yearmonthdaynb['YMDnb']))):
                    continue
                
                print('        '+folder+' holds data to download')
                folder_name=[]
                folder_name=folder_year_name + folder + '/'
                ftp.cwd(folder_name)
                #print("Now in folder " + folder_name)
            
                # Get all Files in this folder
                files=[]
                files = ftp.nlst()
                
                # Print out and download the files
                if (download_images=='TRUE'):
                    #Define the path to save for image saving
                    path_to_save='C:/Users/jullienn/Documents/working_environement/iceslabs_MacFerrin/data/' + folder_year + '/' + 'images/'
                    
                    for file in files:
                        if (not(file[0:15] in list(all_data_to_download))):
                            #Not a data to download
                            continue
                        else:
                            if (os.path.isfile(path_to_save + file)):
                                #If the file have already been downloaded, continue
                                print('                '+file+' have already been downloaded. Continue ...')
                                continue
                            #Download images
                            print("                Downloading... " + file)
                            ftp.retrbinary("RETR " + file ,open(path_to_save + file, 'wb').write)
                            
                            #Display some information about progress
                            perc_display=np.where(all_data_to_download==file[0:15])[0][0]/len(all_data_to_download)*100
                            print(str(perc_display)+' % done')
                
                if (download_mat=='TRUE'):
                    #Define the path to save for mat file saving
                    path_to_save='C:/Users/jullienn/Documents/working_environement/iceslabs_MacFerrin/data/' + folder_year + '/' + 'CSARP_qlook/' + folder + '/'
                    
                    #Create the directory to store the data if does not exist yet
                    #this is from: https://thispointer.com/how-to-create-a-directory
                    #-in-python/#:~:text=Python%27s%20OS%20module%20provides%20an%20
                    #another%20function%20to%20create%20a%20directories%20i.e.&text=.
                    #makedirs(path)-,os.,mkdir%20%2Dp%20command%20in%20linux.
                    if (not os.path.exists(path_to_save)):
                        os.mkdir(path_to_save)
                        print('            Created the directory'+path_to_save)
                    
                    if (folder_year=='2012_Greenland_P3'):
                        #Download Data_img data
                        for file in files:
                            if (not(file[12:27] in list(all_data_to_download))):
                                #Not a data to download
                                continue
                            else:
                                
                                if (os.path.isfile(path_to_save + file)):
                                    #If the file have already been downloaded, continue
                                    print('                '+file+' have already been downloaded. Continue ...')
                                    continue
                                #Download .mat files
                                print("                Downloading... " + file)
                                ftp.retrbinary("RETR " + file ,open(path_to_save + file, 'wb').write)
                                
                                #Display some information about progress
                                perc_display=np.where(all_data_to_download==file[12:27])[0][0]/len(all_data_to_download)*100
                                print(str(perc_display)+' % done')
                    
                    else:
                        for file in files:
                            if (not(file[5:20] in list(all_data_to_download))):
                                #Not a data to download
                                continue
                            else:
                                
                                if (os.path.isfile(path_to_save + file)):
                                    #If the file have already been downloaded, continue
                                    print('                '+file+' have already been downloaded. Continue ...')
                                    continue
                                #Download .mat files
                                print("                Downloading... " + file)
                                ftp.retrbinary("RETR " + file ,open(path_to_save + file, 'wb').write)
                                
                                #Display some information about progress
                                perc_display=np.where(all_data_to_download==file[5:20])[0][0]/len(all_data_to_download)*100
                                print(str(perc_display)+' % done')
        
        else:
            print('Not 2010-11-12-13-14, continue')
            continue
            
    ftp.close()
    end = datetime.now()
    diff = end - start
    print('All files downloaded for ' + str(diff.seconds) + 's')

########################## Download 2010-2014 AR data #########################

############################# Download 2017 AR data #############################
#Code from: https://gist.github.com/nasrulhazim/cfd5f01e3b261b09d54f721cc1a7c50d
if (year_to_download=='2017'):
    #pdb.set_trace()
    
    #Define the folder to not download
    to_not_download=['20170329_02','20170330_01','20170330_02','20170330_03','20170330_05',
                 '20170403_04','20170403_06','20170406_01','20170406_02','20170407_01',
                 '20170407_02','20170407_04','20170426_01','20170426_02','20170428_01',
                 '20170506_02','20170510_01']
    
    from ftplib import FTP
    from datetime import datetime
    
    #Set data we want to download
    download_images='FALSE'
    download_mat='TRUE'
    
    start = datetime.now()
    ftp = FTP('data.cresis.ku.edu')
    ftp.login()
    
    #path='/data/accum/old_format/2003/'
    path='/data/accum/'
    ftp.cwd(path)
    
    # Get folders_years name
    folders_years = ftp.nlst()
    
    #Set the path where to save data
    path_save='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/'
    
    for folder_year in folders_years:
        if (folder_year == '2017_Greenland_P3'):
            #pdb.set_trace()
            
            print('Downloading 2017 data')
            
            #Go to folder year
            folder_year_name=[]
            
            if (download_images=='TRUE'):
                folder_year_name=path + folder_year + '/images/'
                path_save=path_save + folder_year + '/images/'
            
            if(download_mat=='TRUE'):
                folder_year_name=path + folder_year + '/CSARP_standard/'
                path_save=path_save + folder_year + '/CSARP_qlook/'
            
            #Go to folder to download
            ftp.cwd(folder_year_name)
            
            # For this particular year, get folders name
            folders=[]
            folders = ftp.nlst()
            
            #Loop over the folders, and download all the data in this folder
            for folder in folders:
                
                if folder in list(to_not_download):
                    #Not Greenland, do not bother download it
                    continue
                
                if (folder == '20170322_04'):
                    print('Folder 20170322_04 does not exist for .mat files, continue')
                    continue
                                
                folder_name=[]
                folder_name=folder_year_name + folder + '/'
                ftp.cwd(folder_name)
                print("Now in folder ..." + folder_name)
                
                #Update the folder name to store .mat file data
                if(download_mat=='TRUE'):
                    path_to_save=path_save + folder + '/'
                    
                    #If the folder does not exist, create it
                    if not(os.path.isdir(path_to_save)):
                        os.mkdir(path_to_save)
                
                if (download_images=='TRUE'):
                    path_to_save=path_save
                
                # Get all Files in the folder in the cresis data repository
                files=[]
                files = ftp.nlst()
                
                # Print out and download the files
                if (download_images=='TRUE'):
                    for file in files:
                        if (os.path.isfile(path_to_save + file)):
                            #If the file have already been downloaded, continue
                            print(file+' have already been downloaded. Continue ...')
                            continue
                        
                        print("Downloading..." + file)
                        ftp.retrbinary("RETR " + file ,open(path_to_save + file, 'wb').write)
                    
                if (download_mat=='TRUE'):
                    for file in files:
                        if (file[0:9]=='Data_2017'):
                            if (os.path.isfile(path_to_save + file)):
                                #If the file have already been downloaded, continue
                                print(file+' have already been downloaded. Continue ...')
                                continue
                            #Grab only the files starting by 'Data_2017...'
                            print("Downloading..." + file)
                            ftp.retrbinary("RETR " + file ,open(path_to_save + file, 'wb').write)
                        else:
                            print('This is a file Data_img ...')
                            #This is data starting by 'Data_img...', we do not want that
                            continue
                
        else:
            print('This is not 2017, continue')
            
    ftp.close()
    end = datetime.now()
    diff = end - start
    print('All files downloaded for ' + str(diff.seconds) + 's')

############################# Download 2017 AR data #############################


############################# Download 2018 AR data #############################
#Code from: https://gist.github.com/nasrulhazim/cfd5f01e3b261b09d54f721cc1a7c50d
if (year_to_download=='2018'):
    #pdb.set_trace()
    
    #Define the folder to not download
    to_not_download=['20180315_01','20180315_02','20180315_03','20180315_04']
    
    from ftplib import FTP
    from datetime import datetime
    
    #Set data we want to download
    download_images='TRUE'
    download_mat='FALSE'
    
    start = datetime.now()
    ftp = FTP('data.cresis.ku.edu')
    ftp.login()
    
    #path='/data/accum/old_format/2003/'
    path='/data/accum/'
    ftp.cwd(path)
    
    # Get folders_years name
    folders_years = ftp.nlst()
    
    #Set the path where to save data
    path_save='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/'
    
    for folder_year in folders_years:
        if (folder_year == '2018_Greenland_P3'):
            #pdb.set_trace()
            
            print('Downloading 2018 data')
            
            #Go to folder year
            folder_year_name=[]
            
            if (download_images=='TRUE'):
                folder_year_name=path + folder_year + '/images/'
                path_save=path_save + folder_year + '/images/'
            
            if(download_mat=='TRUE'):
                folder_year_name=path + folder_year + '/CSARP_qlook/'
                path_save=path_save + folder_year + '/CSARP_qlook/'
            
            #Go to folder to download
            ftp.cwd(folder_year_name)
            
            # For this particular year, get folders name
            folders=[]
            folders = ftp.nlst()
            
            #Loop over the folders, and download all the data in this folder
            for folder in folders:
                
                if folder in list(to_not_download):
                    #Not Greenland, do not bother download it
                    continue
                
                folder_name=[]
                folder_name=folder_year_name + folder + '/'
                ftp.cwd(folder_name)
                print("Now in folder ..." + folder_name)
                
                #Update the folder name to store .mat file data
                if(download_mat=='TRUE'):
                    path_to_save=path_save + folder + '/'
                    
                    #If the folder does not exist, create it
                    if not(os.path.isdir(path_to_save)):
                        os.mkdir(path_to_save)
                
                if (download_images=='TRUE'):
                    path_to_save=path_save
                
                # Get all Files in the folder in the cresis data repository
                files=[]
                files = ftp.nlst()
                
                # Print out and download the files
                if (download_images=='TRUE'):
                    for file in files:
                        if (os.path.isfile(path_to_save + file)):
                            #If the file have already been downloaded, continue
                            print(file+' have already been downloaded. Continue ...')
                            continue
                        
                        print("Downloading..." + file)
                        ftp.retrbinary("RETR " + file ,open(path_to_save + file, 'wb').write)
                    
                if (download_mat=='TRUE'):
                    for file in files:
                        if (file[0:9]=='Data_2018'):
                            if (os.path.isfile(path_to_save + file)):
                                #If the file have already been downloaded, continue
                                print(file+' have already been downloaded. Continue ...')
                                continue
                            #Grab only the files starting by 'Data_2017...'
                            print("Downloading..." + file)
                            ftp.retrbinary("RETR " + file ,open(path_to_save + file, 'wb').write)
                        else:
                            print('This is a file Data_img ...')
                            #This is data starting by 'Data_img...', we do not want that
                            continue
                
        else:
            print('This is not 2018, continue')
            
    ftp.close()
    end = datetime.now()
    diff = end - start
    print('All files downloaded for ' + str(diff.seconds) + 's')

############################# Download 2018 AR data #############################
