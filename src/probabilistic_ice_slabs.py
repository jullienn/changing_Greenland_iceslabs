# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 08:26:05 2021

@author: jullienn
"""

def _export_to_8bit_array(array):
    #This is from MacFerrin et al., 2019, IceBridgeGPR_Manager_v2.py
    '''In order to export a function to a PNG image, use this funciton to
    export to an 8 bit unsigned integer array of scaled values.'''

    output_array = np.zeros(array.shape, dtype=np.uint8)
    excluded_mask = np.isnan(array)

    range_min = 0
    range_max = 2**8 - 1
    # Get the data minimum and maximum while cutting off 0.5% of outliers
    nonzero_values = array[~excluded_mask]
    data_cutoff_min = np.percentile(nonzero_values,  0.5)
    data_cutoff_max = np.percentile(nonzero_values, 99.5)

    export_array_rescaled = (array - data_cutoff_min) / (data_cutoff_max - data_cutoff_min) * range_max
    # Round to integer values
    export_array_rescaled_int = np.rint(export_array_rescaled)
    # Saturate at top & bottom
    export_array_rescaled_int[export_array_rescaled_int < range_min] = range_min
    export_array_rescaled_int[export_array_rescaled_int > range_max] = range_max
    # Set all numpy.nan values to zero
    export_array_rescaled_int[excluded_mask] = range_min
    # plug into the integer array (conversion from larger to smaller integers)
    output_array[:,:] = export_array_rescaled_int[:,:]
    
    return output_array

#1. Open the data
#2. Loop over dates and load data
#3. Loop over quantiles and calculate probability
#4. Create probabilitic slabs plots and excels files

import pandas as pd
import numpy as np
import pdb
import pickle
import png

import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

import time
import os.path
import glob

#I. Define path, open datetracks and define desired quantiles
#Define path where to pick roll corrected data
'''
path_quantiles_data='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2010_2018/custom_threshold_method/pickles/'
'''
path_quantiles_data='/flash/jullienn/data/threshold_processing_output/pickles/'

#Identify all the datetraces to process
'''
path_datetrack='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/'
'''
path_datetrack='/flash/jullienn/data/threshold_processing/'
datetrack_toread = np.asarray(pd.read_csv(path_datetrack+'datetrack_20102018.txt', header=None))

#Define the desired quantiles over which we will loop
desired_quantiles=np.round(np.arange(0.63,0.82,0.01),2)

#intialize counter to 0
count_time=0

#II. Loop over these datetracks, and perform probability calculation:
for indiv_trace in datetrack_toread:
    
    #pdb.set_trace()
    #If pickle files have already been created, do not process and continue
    filename_to_check='/flash/jullienn/data/threshold_processing_output/probability_iceslabs/pickles/'+indiv_trace[0]+'*'
    
    if (len(glob.glob(filename_to_check))>0):
        print(indiv_trace[0],': files already existent, move on to the next date')
        continue
    
    #To access advance
    start = time.time()
    print(indiv_trace[0])
    print(count_time/len(datetrack_toread)*100,'%')
    
    #Must open the first quantile (0.63) to know the size 
    filename_quantile_open=indiv_trace[0]+'_SG1_cutoffisquantile_'+str(0.63)+'_threshold_350.pickle'
                         
    #Open the corresponding quantile 0.63 file
    f_quantile = open(path_quantiles_data+filename_quantile_open, "rb")
    indiv_quantile063_slice = pickle.load(f_quantile)
    f_quantile.close()
    
    #Set the probabilistic_slice to zeros
    probabilistic_slice=np.zeros((indiv_quantile063_slice.shape[0],indiv_quantile063_slice.shape[1]))
    
    #Loop over the quantiles, load data and perform probability calculation
    for indiv_quantile in desired_quantiles:
        #print(str('%.2f' % indiv_quantile))
        
        #Define filename of quantiles of interest
        filename_quantile_open=indiv_trace[0]+'_SG1_cutoffisquantile_'+str('%.2f' % indiv_quantile)+'_threshold_350.pickle'
                    
        #Open the corresponding quantile file
        f_quantile = open(path_quantiles_data+filename_quantile_open, "rb")
        indiv_quantile_slice=pickle.load(f_quantile)
        f_quantile.close()
        
        #Add up the numbers
        probabilistic_slice=probabilistic_slice+indiv_quantile_slice
    
    #pdb.set_trace()
    
    #Divide the probabilistic_slice by the number of quantiles to have a probability map
    probabilistic_slice=probabilistic_slice/len(desired_quantiles)
    
    #Save the image
    '''
    fig_name='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2010_2018/custom_threshold_method/images/'+indiv_trace[0]+'_probability_iceslabs_presence.png'
    '''
    fig_name='/flash/jullienn/data/threshold_processing_output/probability_iceslabs/images/'+indiv_trace[0]+'_probability_iceslabs_presence.png'
    
    '''
    #Traditional was of plotting, depreciated here
    fig, (ax1) = plt.subplots(1, 1)
    
    #Plot custom threshold ice slabs identification
    cb=ax1.imshow(probabilistic_slice,cmap=plt.get_cmap('Blues'))#,norm=divnorm)
    ax1.title.set_text(indiv_trace[0]+' - ice slabs presence probability (quantile 0.63-0.81)')
    
    plt.show()
    
    #Save the figure
    plt.savefig(fig_name,dpi=2000)
    plt.close(fig)
    '''
    
    '''
    #If one wants to check that _export_to_8bit_array does not more than
    #rescaling the image between 0 and 255
    
    probabilistic_slice_png_toplot=_export_to_8bit_array(probabilistic_slice)
    
    fig, (ax1,ax2) = plt.subplots(2, 1)
    ax1.imshow(probabilistic_slice,cmap=plt.get_cmap('Blues'))
    ax2.imshow(probabilistic_slice_png_toplot,cmap=plt.get_cmap('Blues'))
    plt.show()
    '''
    
    #Prepare matrix for png plot. (1-probabilistic_slice) because 1 is white
    #out of the function _export_to_8bit_array, and I want black
    probabilistic_slice_png=_export_to_8bit_array((1-probabilistic_slice))
    
    #Save the image
    png_to_save=png.from_array(probabilistic_slice_png, mode='L')
    png_to_save.save(fig_name)
    
    #Save the pickle file
    '''
    filename_tosave='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2010_2018/custom_threshold_method/pickles/prob/'+indiv_trace[0]+'_probability_iceslabs_presence.pickle'
    '''
    filename_tosave='/flash/jullienn/data/threshold_processing_output/probability_iceslabs/pickles/'+indiv_trace[0]+'_probability_iceslabs_presence.pickle'
    
    outfile= open(filename_tosave, "wb" )
    pickle.dump(probabilistic_slice,outfile)
    outfile.close()

print('End of probabilistic processing')