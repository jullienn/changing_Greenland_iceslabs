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

#Define the working environment
path= 'C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/exported/refine_exclusion/'

f = open('C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/Exclusion_folder/data_2018_toberun.txt','r')
dates_surf_2018 = [line.strip() for line in f.readlines() if len(line.strip()) > 0]
f.close()

#create a file text storing the indiv_file name and the exclusions
path_exclusionfile_store='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/Exclusion_folder/exclusion_16m_2017_2018/'
f_exclusions = open(path_exclusionfile_store+'exclusion_16m_2018.txt', "w")

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
    
    #Show the exclusion
    ax1.plot(exclusion_vect,np.ones(arr_file_to_show_boolean.shape[1]),color='red')
                
    #Save the figure
    fig_name=[]
    plt.savefig(path_exclusionfile_store+'exclusion_'+indiv_file,dpi=1000)
    plt.close()
    
    #Compute the exclusions limits
    
    #itialize the start and end pixel
    start_pix=np.nan
    end_pix=np.nan
    #initalize the suite of exclusion pixels
    exclusion_suite=np.nan
    
    for i in range(0,exclusion_vect.size):
        if (np.isnan(exclusion_vect[i])):
            #pixel of interest is nan. Several possibilities
            #1. Exclusion suite have been defined, we are at the end, i.e.
            #start pixel is not nan but end pixel is
            if (not(np.isnan(start_pix)) and np.isnan(end_pix)):
                #pdb.set_trace()
                #identify end pixel
                end_pix=i-1
                #if the start pix and end pix are identical, do not store and continue
                if (start_pix==end_pix):
                    #reset start and end pixels to nan
                    start_pix=np.nan
                    end_pix=np.nan
                    continue
                
                #store the exclusion suite just generated
                exclusion_suite=np.append(exclusion_suite,[str(start_pix)+'-'+str(end_pix)])
                #reset start and end pixels to nan
                start_pix=np.nan
                end_pix=np.nan
            else:
                #2. No exclusion suite have been identifies yet
                continue
        else:
            #pixel of interest not a nan, suite to save
            if not(np.isnan(start_pix)):
                #A suite is currently being identified
                if (i==(exclusion_vect.size-1)):
                    #it means the last pixel is also and exclusion but the process have missed it
                    #identify end pixel
                    end_pix=i
                    
                    #store the last exclusion suite
                    exclusion_suite=np.append(exclusion_suite,[str(start_pix)+'-'+str(end_pix)])
                    
                    #reset start and end pixels to nan
                    start_pix=np.nan
                    end_pix=np.nan
                else:
                    #we aleary are in a suite and this is not the end, continue
                    continue
            else:
                #Suite not being identified, this is a new one => define the start pixel
                start_pix=i
                #pdb.set_trace()
    
    if (len(str(exclusion_suite))==3):
        #No exclusion for this date
        f_exclusions.write(indiv_file+'\n')
    else:
        #Exclusion for this date
        #Remove the nan from the exclusion_suite
        list_exclusion_suite=list(exclusion_suite)
        list_exclusion_suite=list_exclusion_suite[1:] #remove the first one which is a nan
        
        f_exclusions.write(indiv_file+' ')
        #Store the exclusions for each date
        for indiv_excl in list_exclusion_suite:
            f_exclusions.write(indiv_excl+' ')
        
        f_exclusions.write('\n')

#Close the exclusion_16m file
f_exclusions.close()

