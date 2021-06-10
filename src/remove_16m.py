# -*- coding: utf-8 -*-
"""
Created on Thu Jun 10 10:18:09 2021

@author: JullienN
"""

from IPython.display import Image
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from PIL import Image
from os import listdir
from os.path import isfile, join
import os
import pdb

##############################################################################
################### Define function for ice lenses logging ###################
##############################################################################
#This function if adapted from https://stackoverflow.com/questions/37363755/python-mouse-click-coordinates-as-simply-as-possible
def onclick(event):
    #This functions print and save the x and y coordinates in pixels!
    print(event.xdata, event.ydata)
    #Fill in the file to log on the information
    #filename_flog='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/2002_2003_radar_slice/flog_icelenses_alldates.txt'
    #f_log = open(filename_flog, "a")
    #f_log.write(str(round(event.xdata,2))+','+str(round(event.ydata,2))+'\n')
    #f_log.close() #Close the quality assessment file when we’re done!
##############################################################################
################### Define function for ice lenses logging ###################
##############################################################################


#Define the working environment
path= 'C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/exported/refine_exclusion/'

f = open('C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/Exclusion_folder/data_2017_toberun.txt','r')
dates_surf_2018 = [line.strip() for line in f.readlines() if len(line.strip()) > 0]
f.close()

for indiv_file in list(dates_surf_2018):
    print(indiv_file)
    #Open only the _orig_CUTOFF_-0.45_THRESHOLD_000 files of interest
    filename=indiv_file+'_orig_CUTOFF_-0.45_THRESHOLD_000.png'
    
    #Open image
    file_to_show = Image.open(path+filename).convert("L")
    arr_file_to_show = np.asarray(file_to_show)
    
    #Convert arr_file_to_show into a boolean (from 0 and 255 to 1 and 0).
    #When pixel=255, it means white. When pixel=0, it means black. I want to
    #count the number of times I have blacks pixels (ice), thus I need to
    #transform 255 into 0 and 0 in 1. This is what is done below.
    
    arr_file_to_show_boolean=arr_file_to_show==0
    arr_file_to_show_boolean=arr_file_to_show_boolean*1
    
    #Retreive the 16m depth index
    depth20m=arr_file_to_show_boolean.shape[0]
    depth16m=arr_file_to_show_boolean.shape[0]*16/20
    
    exclusion=[]
    exclusion_vect=np.zeros(arr_file_to_show_boolean.shape[1])
    exclusion_vect[:]=np.nan
        
    #Check if ice content >16m on each vertial column. If yes, create an exclusion pixel
    for i in range(0,arr_file_to_show_boolean.shape[1]):
        if (np.sum(arr_file_to_show_boolean[:,i])>=depth16m):
            #more than 16m of ice in this vertical column, retreive exclusion pixel
            exclusion=np.append(exclusion,i)
            exclusion_vect[i]=i
    
    #Create a visualisation exclusion vector
    #Prepare plot
    fig, (ax1) = plt.subplots(1, 1)#, gridspec_kw={'width_ratios': [1, 3]})
    fig.suptitle(indiv_file)
    ax1.imshow(arr_file_to_show, cmap='gray', vmin=0, vmax=255)
    fig.canvas.mpl_connect('button_press_event', onclick)
    
    #Show the exclusion
    ax1.plot(exclusion_vect,np.ones(arr_file_to_show_boolean.shape[1]))   
    
    plt.show()
    pdb.set_trace()
    
