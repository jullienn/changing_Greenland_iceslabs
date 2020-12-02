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

        folder_year_name=[]
        folder_year_name=path + folder_year + '/'
        ftp.cwd(folder_year_name)
        
        # For this particular year, get folders name where we have data in SW Greenland
        open ...
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
    else:
        print('This is not 2017')
        
ftp.close()
end = datetime.now()
diff = end - start
print('All files downloaded for ' + str(diff.seconds) + 's')

############################# Download 2017 AR data #############################
