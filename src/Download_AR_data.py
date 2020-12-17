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
############################# Download old AR data #############################
##Code from: https://gist.github.com/nasrulhazim/cfd5f01e3b261b09d54f721cc1a7c50d
#
#from ftplib import FTP
#from datetime import datetime
#
#
#start = datetime.now()
#ftp = FTP('data.cresis.ku.edu')
#ftp.login()
#
##path='/data/accum/old_format/2003/'
#path='/data/accum/old_format/'
#ftp.cwd(path)
#
## Get folders_years name
#folders_years = ftp.nlst()
#
#for folder_year in folders_years:
#    
#    if (folder_year == '2006'):
#        print('Break because year 2006')
#        break
#    elif (folder_year == '2009'):
#        print('Break because year 2009')
#        break
#    elif (folder_year == '2011'):
#        print('Break because year 2011')
#        break
#    else:
#        print('This is either the year 2002 or 2003')
#        folder_year_name=[]
#        folder_year_name=path + folder_year + '/'
#        ftp.cwd(folder_year_name)
#        
#        # For this particular year, get folders name
#        folders=[]
#        folders = ftp.nlst()
#        #Go to the folder
#        
#        for folder in folders:
#            folder_name=[]
#            folder_name=path + folder_year + '/' + folder + '/'
#            ftp.cwd(folder_name)
#            print("Now in folder ..." + folder_name)
#        
#            # Get All Files in this folder
#            files=[]
#            files = ftp.nlst()
#            
#            # Print out the files
#            for file in files:
#                print("Downloading..." + file)
#                ftp.retrbinary("RETR " + file ,open("D:/OIB/AR/" + folder_year + '/' + folder + "/" + file, 'wb').write)
#                #ftp.retrbinary("RETR " + file ,open("download/to/your/directory/" + file, 'wb').write)
#
#ftp.close()
#end = datetime.now()
#diff = end - start
#print('All files downloaded for ' + str(diff.seconds) + 's')
############################# Download old AR data #############################


############################# Download 2017 AR data #############################
#Code from: https://gist.github.com/nasrulhazim/cfd5f01e3b261b09d54f721cc1a7c50d

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
        pdb.set_trace()
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
