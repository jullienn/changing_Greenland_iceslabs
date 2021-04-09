# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 18:01:56 2021

@author: JullienN
"""

#Adapted from https://stackoverflow.com/questions/37467561/renaming-multiple-files-in-a-directory-using-python

import os
from os import listdir
from os.path import isfile, join
import pdb

pdb.set_trace()

#Define the working environment
path= 'C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data'
os.chdir(path) # relative path: scripts dir is under Lab

# Read the years of data
folder_years = [ f.name for f in os.scandir(path) if f.is_dir() ]

for folder_year in folder_years:
    print(folder_year)
    if (folder_year in list(['2012_Greenland_P3'])):
        
        #Go into the yearly folders 
        folder_year_name=path+'/'+folder_year+'/CSARP_qlook'
        os.chdir(folder_year_name)

        # Read the days of this specific year
        folder_days = [ f.name for f in os.scandir(folder_year_name) if f.is_dir() ]
        
        for folder_day in folder_days:
            print('Now in year',folder_year,'day',folder_day)
            
            #Go into the daily folders 
            folder_day_name=folder_year_name+'/'+folder_day
            os.chdir(folder_day_name)
            
            # Read the files of this specific day
            onlyfiles = [f for f in listdir(folder_day_name) if isfile(join(folder_day_name, f))]
            #pdb.set_trace()
            
            for indiv_file in onlyfiles:
                print('Treating file',indiv_file)
                #pdb.set_trace()
                
                #Rename the files
                #os.rename(os.path.join(path, file), os.path.join(path,file[0:5]+file[12:31]))
    else:
        print('Not desired folder, continue')
        continue
    

print('End of processing')