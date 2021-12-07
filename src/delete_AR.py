# -*- coding: utf-8 -*-
"""
Created on Fri May 28 18:02:33 2021

@author: JullienN

This code delete AR data for 2017 and 2018 according to data to kept that are
listed in the file data_2017_2018.xlsx
"""

import os
from os import listdir
import os.path
from os.path import isfile, join
import pdb
import numpy as np

#Here choose mat files or jpg files
to_delete_is='mat' #can either be 'jpg' of 'mat'

#Open excel file with refined 2017 data
data2017_to_be_kept='data_2018_toberun.txt'
f = open('C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/Exclusion_folder/'+data2017_to_be_kept,'r')
lines = [line.strip() for line in f.readlines() if len(line.strip()) > 0]
f.close()

#Define all the individual names of files from this list
indiv_filenames=[]
for file_sequence in lines:
    file_start=file_sequence[12:15]
    file_end=file_sequence[16:19]
    
    #number_of_files=int(file_end)-int(file_start)+1
    
    for i in range(int(file_start),int(file_end)+1):
        if (i<10):
            constructed_filename=file_sequence[0:12]+'00'+str(i)
        elif ((i>9) and (i<100)):
            constructed_filename=file_sequence[0:12]+'0'+str(i)
        else:
            constructed_filename=file_sequence[0:12]+str(i)

        indiv_filenames=np.append(indiv_filenames,constructed_filename)

path_data='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data'
os.chdir(path_data)

# Read the folders of data
folders_year = [ f.name for f in os.scandir(path_data) if f.is_dir() ]

#Loop over the folder
for folder in folders_year:
    if (not(folder in list(['2018_Greenland_P3']))):#,'2018_Greenland_P3']))):
        print('Not 2017, continue')
        continue
    else:
        print('This is year ',folder)
        if (to_delete_is=='mat'):
            #Define path
            path_to_scan=path_data+'/'+folder+'/CSARP_qlook/'
            
            #Read the daily folders
            folders_day = [ f.name for f in os.scandir(path_to_scan) if f.is_dir() ]
            
            for folder_day in folders_day:
                #Update path
                path_indiv_day=path_to_scan+folder_day
                #List all the files
                onlyfiles = [f for f in listdir(path_indiv_day) if isfile(join(path_indiv_day, f))]
                
                for indiv_file in onlyfiles:
                    if (not(indiv_file[5:20] in list(indiv_filenames))):
                        #Delete the file
                        #This is from https://stackoverflow.com/questions/6996603/how-to-delete-a-file-or-folder
                        if os.path.isfile(path_indiv_day+'/'+indiv_file):
                            os.remove(path_indiv_day+'/'+indiv_file)
                        else:
                            #Print an error
                            print("Error: %s file not found" % indiv_file)
                    else:
                        #We keep the file
                        continue
        
        elif (to_delete_is=='jpg'):
            #Define paths
            path_to_scan=path_data+'/'+folder+'/images/'
            
            #List all the files
            onlyfiles = [f for f in listdir(path_to_scan) if isfile(join(path_to_scan, f))]
            for indiv_file in onlyfiles:
                if (not(indiv_file[0:15] in list(indiv_filenames))):
                    #Delete the file
                    #This is from https://stackoverflow.com/questions/6996603/how-to-delete-a-file-or-folder
                    if os.path.isfile(path_to_scan+indiv_file):
                        os.remove(path_to_scan+indiv_file)
                    else:
                        #Print an error
                        print("Error: %s file not found" % indiv_file)
                else:
                    #We keep this file
                    continue
        else:
            print('Change format of file!')

