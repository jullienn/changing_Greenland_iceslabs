# -*- coding: utf-8 -*-
"""
Created on Mon Mar 29 09:16:52 2021

@author: JullienN
"""

#Import librairies
import pdb
import pandas as pd
import numpy as np
import scipy.io

#Define the path where the data are stored
path_data='C:/Users/jullienn/Documents/working_environement/iceslabs_MacFerrin/data/'
#1. Open exclusion file
path_exclusion_folder='C:/Users/jullienn/Documents/working_environement/iceslabs_MacFerrin/data/Exclusion_folder/txt/'

f = open(path_exclusion_folder+'LAKES_AND_OTHER_EXCLUSIONS.txt','r')
lakes_and_other_exclusions = [line.strip() for line in f.readlines() if len(line.strip()) > 0]
f.close()

#2. Retreive only the dates from this exclusion file
#Loop over the lakes_and_other_exclusions file to keep only the dates
all_data=[]
lakes_and_other_exclusions_20102014=[]

for indiv_file in list(lakes_and_other_exclusions):
    if (indiv_file[0:4]=='2017'):
        continue
    all_data=np.append(all_data,indiv_file[0:19])
    lakes_and_other_exclusions_20102014=np.append(lakes_and_other_exclusions_20102014,indiv_file)

#3. Loop over the dates of the exclusion file
loc=0
for indiv_trace in list(all_data):
    print('Track being read: '+indiv_trace)
    # ---------------- This is from 'Download_AR_data.py' ------------------- #

    #List all the individual trace nb of that date
    traces_nb=indiv_trace[12:19]
    begin_trace_nb=int(traces_nb[0:3])
    end_trace_nb=int(traces_nb[4:7])
    
    #Create a vector from begin_trace_nb to end_trace_nb
    vect_traces_nb=np.arange(begin_trace_nb,end_trace_nb+1,1)
    
    #Create the right suffixe for immediate file opening
    vect_str=[]
    lat_append=[]
    lon_append=[]
    
    for trace_nb_indiv in list(vect_traces_nb):
        if (trace_nb_indiv<10):
            trace_nb_indiv_str='00'+str(trace_nb_indiv)
        elif ((trace_nb_indiv>=10) and (trace_nb_indiv<100)):
            trace_nb_indiv_str='0'+str(trace_nb_indiv)
        else:
            trace_nb_indiv_str=str(trace_nb_indiv)
        
        # ---------------- This is from 'Download_AR_data.py' ------------------- #
        #Define folders where to load the data
        folder_data=path_data+indiv_trace[0:4]+'_Greenland_P3/CSARP_qlook/'+indiv_trace[0:11]
        #Load files of that trace
        fdata=scipy.io.loadmat(folder_data+'/Data_'+indiv_trace[0:11]+'_'+trace_nb_indiv_str+'.mat')
        #Load lat and lon
        lat=fdata['Latitude']
        lon=fdata['Longitude']
        
        #Append the data to each other for pixel identification
        lat_append=np.append(lat_append,lat)
        lon_append=np.append(lon_append,lon)
    
    #Load the pixels exclusions
    line_of_interest=lakes_and_other_exclusions_20102014[loc]
    
    if (len(line_of_interest)==19):
        print('No exclusion for '+indiv_trace)
    else:
        #Several scenarios:
        #1. -222 : means that all the pixel from the begining until 222 are excluded
        #2. 234-287 : means that all the pixels between 234 and 287 are excluded
        #3. 567- : meand that all the pixels from 567 until the end of the trace are excluded.
        
        #Know whether the pixel is included or not in the exclusion
        
        pdb.set_trace()
        print('lala')
        
        index_exclusion=[]
        
        for i in range(0,len(line_of_interest.split())-1):
            #-1 because the first one is the date
            #if len(line_of_interest.partition(" "))>= 2, then there is a least one exclusion to consider
            #i+1
            
            exclusion_to_consider=line_of_interest.split()[i+1]
            #pdb.set_trace()
            if (exclusion_to_consider[0]=="-"):
                #Scenario 1
                begin_pixel=0
                end_pixel=int(exclusion_to_consider.partition("-")[2])
            elif (exclusion_to_consider[-1]=="-"):
                #Scenario 3
                begin_pixel=int(exclusion_to_consider.partition("-")[0])
                end_pixel=len(lon_append)
            elif (exclusion_to_consider.partition("-")[1]=='-'):
                #Scenario 2
                begin_pixel=int(exclusion_to_consider.partition("-")[0])
                end_pixel=int(exclusion_to_consider.partition("-")[2])
            else:
                print('Problem in identifying the exclusion')
            pdb.set_trace()
            print('lala')
            #I check, OK the above works!!
            
            #Create a vector from begin_pixel to end_pixel
            vect_begin_end=np.arange(begin_pixel,end_pixel+1,1)
            
            #Append all the indexes of exclusion
            index_exclusion=np.append(index_exclusion,vect_begin_end).astype(int)
        
        pdb.set_trace()
        #Extract the corresponding lon/lat exclusion
        lat_exclusion=lat[0,index_exclusion]
        lon_exclusion=lon[0,index_exclusion]
        
        #Save the lat/lon with the corresponding date in a dictionnary
            
    
    
    
    
    #Update loc index
    loc=loc+1
    
        
#3. In this loop, open individual files and store retreive the associated lat/lon

#4. Save the associated lat/lon (save the vector of all the lat and lon in between)