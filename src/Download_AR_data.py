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

year_to_download='2010_2014'

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
    
    #Set data we want to download
    download_images='TRUE'
    download_mat='FALSE'
    
    #Load the data we have to download
    f = open('C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/2010_2014_data_download.txt','r')
    lines = [line.strip() for line in f.readlines() if len(line.strip()) > 0]
    f.close()
    
    #Extract the folder year month day nb
    x_empty=[None]*len(lines)
    
    df_yearmonthdaynb=pd.DataFrame({'complete_date':x_empty,
                                    'YMDnb':x_empty,
                                    'traces_nb':x_empty})
    
    for i in range (0,len(lines)):
        df_yearmonthdaynb['complete_date'][i]=lines[i]
        df_yearmonthdaynb['YMDnb'][i]=lines[i][0:11]
        df_yearmonthdaynb['traces_nb'][i]=lines[i][12:19]
    
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
    
    for folder_year in folders_years:
        pdb.set_trace()
        if (folder_year in list(list_download)):
            print('Downloading ' + folder_year)
            
            #Go to folder year
            folder_year_name=[]
            
            if (download_images=='TRUE'):
                folder_year_name=path + folder_year + '/images/'
            
            if(download_mat=='TRUE'):
                folder_year_name=path + folder_year + '/CSARP_standard/'
            
            #Go to folder CSARP_standard
            ftp.cwd(folder_year_name)
            
            # For this particular year, get folders name
            folders=[]
            folders = ftp.nlst()
            #Go to the folder
            
            pdb.set_trace()
            #Loop over the folders, and download the data used by MacFerrin et al., 2019 in this folder
            for folder in folders:
                if (not(folder in list(df_yearmonthdaynb['YMDnb']))):
                    print(folder+' does not hold data to download, continue')
                    continue
                
                print(folder+' holds data to download')
                folder_name=[]
                folder_name=folder_year_name + folder + '/'
                ftp.cwd(folder_name)
                print("Now in folder " + folder_name)
            
                # Get all Files in this folder
                files=[]
                files = ftp.nlst()
                
                pdb.set_trace()
                
                # Print out and download the files
                if (download_images=='TRUE'):
                    for file in files:
                        if (os.path.isfile('C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/' + folder_year + '/' + folder + '/' + file)):
                            #If the file have already been downloaded, continue
                            print(file+' have already been downloaded. Continue ...')
                            continue
                        #Grab only the files starting by 'Data_2017...'
                        print("Downloading..." + file)
                        ftp.retrbinary("RETR " + file ,open('C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data' + folder_year + folder + "/" + file, 'wb').write)
                    
                if (download_mat=='TRUE'):
                    #!!!!! FOR 2012, download the Data_img, and rename them as Data_
                    if (folder_year=='2012_Greenland_P3'):
                        print('This is 2012')
                        #Download Data_img data
                    else:
                        for file in files:
                            if (file[0:9]=='Data_2017'):
                                if (os.path.isfile('C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data' + folder_year + folder + "/" + file)):
                                    #If the file have already been downloaded, continue
                                    print(file+' have already been downloaded. Continue ...')
                                    continue
                                #Grab only the files starting by 'Data_2017...'
                                print("Downloading..." + file)
                                ftp.retrbinary("RETR " + file ,open('C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data' + folder_year + folder + "/" + file, 'wb').write)
                            else:
                                print('This is a file Data_img ...')
                                #This is data starting by 'Data_img...', we do not want that
                                continue
                
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
    
    for folder_year in folders_years:
        pdb.set_trace()
        if (folder_year == '2017_Greenland_P3'):
            print('Downloading 2017 data')
            
            #Go to folder year
            folder_year_name=[]
            
            if (download_images=='TRUE'):
                folder_year_name=path + folder_year + '/images/'
            
            if(download_mat=='TRUE'):
                folder_year_name=path + folder_year + '/CSARP_standard/'
            
            #Go to folder CSARP_standard
            ftp.cwd(folder_year_name)
            
            #Set the working directory to load the 2017 data file
            os.chdir('C:\\Users\\jullienn\\Documents\\working_environment\\iceslabs_MacFerrin')
            #Check that we are in the right working directy
            print(os.getcwd())
            # For this particular year, get folders name where we have data in SW Greenland
            #Actually not working, so create the folder array manually
            #folders = pd.read_csv('C:\\Users\\jullienn\\Documents\\working_environment\\iceslabs_MacFerrin\\download_2017_SW.csv')
            folders=['20170421_01',
                     '20170422_01',
                     '20170424_01',
                     '20170429_01',
                     '20170501_02',
                     '20170502_01',
                     '20170505_02',
                     '20170506_01',
                     '20170508_02',
                     '20170510_02',
                     '20170511_01']
            #pdb.set_trace()
            #Loop over the folders, and download all the data in this folder
            for folder in folders:
                folder_name=[]
                folder_name=folder_year_name + folder + '/'
                ftp.cwd(folder_name)
                print("Now in folder ..." + folder_name)
            
                # Get all Files in this folder
                files=[]
                files = ftp.nlst()
                
                # Print out and download the files
                if (download_images=='TRUE'):
                    for file in files:
                        if (os.path.isfile("D:/OIB/AR/" + folder_year + '/' + folder + "/" + file)):
                            #If the file have already been downloaded, continue
                            print(file+' have already been downloaded. Continue ...')
                            continue
                        #Grab only the files starting by 'Data_2017...'
                        print("Downloading..." + file)
                        ftp.retrbinary("RETR " + file ,open("D:/OIB/AR/" + folder_year + '/' + folder + "/" + file, 'wb').write)
                    
                if (download_mat=='TRUE'):
                    for file in files:
                        if (file[0:9]=='Data_2017'):
                            if (os.path.isfile("D:/OIB/AR/" + folder_year + '/' + folder + "/" + file)):
                                #If the file have already been downloaded, continue
                                print(file+' have already been downloaded. Continue ...')
                                continue
                            #Grab only the files starting by 'Data_2017...'
                            print("Downloading..." + file)
                            ftp.retrbinary("RETR " + file ,open("D:/OIB/AR/" + folder_year + '/' + folder + "/" + file, 'wb').write)
                        else:
                            print('This is a file Data_img ...')
                            #This is data starting by 'Data_img...', we do not want that
                            continue
                
        else:
            print('This is not 2017')
            
    ftp.close()
    end = datetime.now()
    diff = end - start
    print('All files downloaded for ' + str(diff.seconds) + 's')

############################# Download 2017 AR data #############################
