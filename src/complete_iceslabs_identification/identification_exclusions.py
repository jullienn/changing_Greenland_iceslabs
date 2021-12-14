# -*- coding: utf-8 -*-
"""
Created on Fri Nov 26 20:03:29 2021

@author: jullienn
"""
#This code is inspired from 'display_raw_2010_2014.py

#_gaussian function taken from IceBridgeGPR_Manager_v2.py
# Define a quick guassian function to scale the cutoff mask above
def _gaussian(x,mu,sigma):
    return np.exp(-np.power((x-mu)/sigma, 2.)/2.)

#This function have been taken from 'IceBridgeGPR_Manager_v2.py
def kernel_function(traces_input,suggested_pixel):
    #pdb.set_trace()
    
    #Compute log10
    traces = traces_input
    traces = np.log10(traces)
    
    # We do not have the original indicies to use as a starter so we use our suggestion for surface picking start
    
    # 3) Perform surface pick crawling threshold behavior mask (assume a step-change analysis [goes from weak->strong at surface], and continuity of surface in original file.)
    # Create a step-change mask to optimze where the returns transition from "dark" to "bright"
    MASK_RADIUS = 50
    vertical_span_mask = np.empty([MASK_RADIUS*2,], dtype=np.float)
    vertical_span_mask[:MASK_RADIUS] = -1.0
    vertical_span_mask[MASK_RADIUS:] = +3.0
    
    vertical_span_mask = vertical_span_mask * _gaussian(np.arange(vertical_span_mask.shape[0]),mu=(MASK_RADIUS-5),sigma=(float(MASK_RADIUS)/3.0))
    
    # Expand the shape to handle array broadcasting below
    vertical_span_mask.shape = vertical_span_mask.shape[0], 1
    
    # This is the vertical window size of the extent of the search.  Should be bigger than any jump from one surface pixel to the next.
    MASK_SEARCH_RADIUS = 150
    
    improved_indices = np.zeros(traces.shape[1], dtype='int64')
    #pdb.set_trace()
    #traces.shape[1] indeed correspond to the horizontal distance
    
    # Start at the left with the hand-picked "suggested surface pick" in the ICEBRIDGE_SURFACE_PICK_SUGGESTIONS_FILE as starting point
    
    last_best_index = suggested_pixel
     
    #pdb.set_trace()
    # A template graph to use, just have to add in the center vertical index at each point and go from there.
    search_indices_template = np.sum(np.indices((vertical_span_mask.shape[0], 2*MASK_SEARCH_RADIUS)),axis=0) - MASK_SEARCH_RADIUS - MASK_RADIUS
    for i in range(traces.shape[1]):
        # Create an array of indices spanning the top-to-bottom of the MASK_SEARCH_RADIUS, and fanning out MASK_RADIUS above and below that point.
        search_indices = search_indices_template + last_best_index
        # Handle overflow indices if below zero or above max (shouldn't generally happen)... just assign to the top or bottom pixel
        search_indices[search_indices < 0] = 0
        search_indices[search_indices >= traces.shape[0]] = traces.shape[0]-1
        
        bestfit_sum = np.sum(traces[:,i][search_indices] * vertical_span_mask, axis=0)
        
        assert bestfit_sum.shape[0] == 2*MASK_SEARCH_RADIUS
        
        # Get the best fit (with the highest value from the transformation fit)
        last_best_index = search_indices[MASK_RADIUS,np.argmax(bestfit_sum)]
        improved_indices[i] = last_best_index
        
    #If there are pixels with particularly strong echo that are being erroneously
    #picked up as the surface, erase most the little "jump" artifacts in
    #the surface picker.
    improved_indices = _get_rid_of_false_surface_jumps(improved_indices)
    
    #I do not use any mask so I think I shouldn't need to use that:
    ###### Must re-expand the surface indices to account for masked values (filled w/ nan)
    ##### improved_indices_expanded = self._refill_array(improved_indices, surface_maskname)
    
    #pdb.set_trace()
    return improved_indices

def _radar_slice_indices_above_and_below(meters_cutoff_above, meters_cutoff_below,depths):
    #pdb.set_trace()

    delta_distance = np.mean(depths[1:] - depths[:-1])
    idx_above = int(np.round(float(meters_cutoff_above) / delta_distance))
    # Add one to the index below to include that last pixel when array-slicing
    idx_below = int(np.round(float(meters_cutoff_below) / delta_distance)) + 1

    return idx_above, idx_below

def _return_radar_slice_given_surface(traces,
                                      depths,
                                      surface_indices,
                                      meters_cutoff_above,
                                      meters_cutoff_below):
    '''From this radar track, return a "slice" of the image above and below the surface by
    (meters_cutoff_above, meters_cutoff_below), respectively.

    Return value:
    A ((idx_below+idx_above), numtraces]-sized array of trace sample values.
    '''
    #pdb.set_trace()
    idx_above, idx_below = _radar_slice_indices_above_and_below(meters_cutoff_above, meters_cutoff_below,depths)

    output_traces = np.empty((idx_above + idx_below, traces.shape[1]), dtype=traces.dtype)
    bottom_indices = np.zeros(shape=(1,traces.shape[1]))
    
    for i,s in enumerate(surface_indices):
        try:
            output_traces[:,i] = traces[(s-idx_above):(s+idx_below), i]
            bottom_indices[0,i]=(s+idx_below)
        except ValueError:
            # If the surf_i is too close to one end of the array or the other, it extends beyond the edge of the array and breaks.
            if s < idx_above:
                start, end = None, idx_above+idx_below
            elif s > (traces.shape[0] - idx_below):
                start, end = traces.shape[0] - (idx_above + idx_below), None
            else:
                # SHouldn't get here.
                print(i, s, traces.shape)
                assert False
            output_traces[:,i] = traces[start:end, i]
            bottom_indices[0,i]=end
    return output_traces, bottom_indices


def _get_rid_of_false_surface_jumps(surface_indices):
    '''Some of the 2011 files especially, have strong echos that are errantly being picked up as the surface.  Find these big "jumps", and get rid of them.  Use the suggested surface instead.'''
    improved_surface = surface_indices.copy()
    
    jumps = improved_surface[1:] - improved_surface[:-1]
    # Substitute any large jumps with brightest pixel in a window of original surface.  Do this until large jumps either go away or have all been corrected to original surface.
    for i in range(len(jumps)):
        
        # Slope windowsize = number of pixels we use to average the previous slope.
        slope_windowsize = 10
        if i < slope_windowsize:
            continue
        mean_slope = np.mean(np.array(jumps[i-slope_windowsize:i], dtype=np.float))

        # Find the difference of this slope from the last five stops
        difference_from_mean_slope = jumps[i] - mean_slope
        # Ignore if it's jumped less than 3 from the mean recent slope, or less than 50% greater than the mean slope at this time.
        if (difference_from_mean_slope < 5) or (difference_from_mean_slope < (1.5*mean_slope)):
            continue

        # tune settings
        jump_lookahead = 20 # Number of pixels to look ahead to see if we can find a matching down-jump
        if i+jump_lookahead > len(jumps):
            jump_lookahead = len(jumps) - i

        # This is how close the surface on the "other side" of the jump must be to the original slope to be considered for it.
        jump_magnitude_threshold = 1.10

        # See if we can find a point in the near future that would approximate the current slope.
        slopes_ahead = np.cumsum(jumps[i:i+jump_lookahead]) / np.arange(1,jump_lookahead+1)
        opposite_match = np.argmax(slopes_ahead <= (mean_slope * jump_magnitude_threshold))
        
        if opposite_match > 0:
            # We found a match, onward!
            opposite_match_index = i + opposite_match
            for j in range(i+1,opposite_match_index+1):
                improved_surface[j] = np.round(improved_surface[i] + float(improved_surface[opposite_match_index+1] - improved_surface[i])*(j-i)/(opposite_match_index+1-i))    
            # now recompute jumps
            jumps = improved_surface[1:] - improved_surface[:-1]
            continue

        # IF THE ABOVE DIDN'T WORK, TRY THE 'JUMP' TECHNIQUE, SEEING WHETHER AN ANOMALOUS 'JUMP' IS COUNTERBALANCED BY AN
        # OPPOSITE AND (APPROXIMATELY) EQUAL JUMP IN THE OPPOSITE DIRECTION.
        # Don't worry about any trends less than 12 pixels.  Hills do that.
        jump = jumps[i]
        if abs(jump) < 5:
            continue

        # tune settings
        jump_lookahead = 50 # Number of pixels to look ahead to see if we can find a matching down-jump
        jump_magnitude_threshold = 0.50 # What fraction of the original jump the new jump has to be (in the opposite direction) to qualify.

        # see if we can find a jump in the near-future that crosses this threshold in the other direction.  If so, we've found our counter-part
        if jump < 0:
            opposite_jump_index = np.argmax((jumps[i:i+jump_lookahead]) > (-jump*jump_magnitude_threshold))
        elif jump > 0:
            opposite_jump_index = np.argmax((jumps[i:i+jump_lookahead]) < (-jump*jump_magnitude_threshold))

        if opposite_jump_index > 0:
            opposite_jump_index += i
        else: # If we didn't find a partner opposite offset, skip and move along.
            continue

        # Linearly interpolate, get to the closest pixel
        try:
            for j in range(i+1,opposite_jump_index+1):
                improved_surface[j] = np.round(improved_surface[i] + float(improved_surface[opposite_jump_index+1] - improved_surface[i])*(j-i)/(opposite_jump_index+1-i))
        except IndexError:
            print("i", i, "j", j, "opposite_jump_index", opposite_jump_index, improved_surface.shape, jumps.shape)
            # Break the program here.
            100/0

        # now recompute jumps
        jumps = improved_surface[1:] - improved_surface[:-1]
        continue
    return improved_surface

#Function taken from IceBridgeGPR_Manager_v2.py

def _export_to_8bit_array(array):
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

#Import packages

import scipy.io
import rasterio
from rasterio.plot import show
import matplotlib.pyplot as plt
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
from pysheds.grid import Grid
import osgeo.ogr as ogr
import osgeo.osr as osr
from pyproj import Transformer
import matplotlib.gridspec as gridspec
import png
import glob

obvious_identification='FALSE'
identification_after_depth_correction='FALSE'
identification_dry_firn_exclusions='FALSE'
generate_exclusion_files='TRUE'

if (generate_exclusion_files=='TRUE'):
        
    #Generate exclusion files according to logboog of exclusions
    path_excel='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/intial_selection_20172018/'
    excel_exclusions=pd.read_csv(path_excel+'logbook_2017_2018_data_processing.csv',sep=';',skiprows=1)
    
    
    #Extract datetrack
    datetrack_to_export=excel_exclusions['datetrack_tobeprocessed']
    #Get rid of dates that have been deleted
    datetrack_to_export=datetrack_to_export[~datetrack_to_export.isnull()]
    #Save the exclusion file
    filename_flog='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/intial_selection_20172018/exclusions/datetrack_20172018.txt'
    f_log = open(filename_flog, "a")
    for i in range(0,len(datetrack_to_export)):
        f_log.write(str(datetrack_to_export.iloc[i])+'\n')
    f_log.close()
    
    
    #Extract obvisous exclusions dataframe
    obvious_exclusions_to_export=excel_exclusions[['datetrack_tobeprocessed','Obvious_exclusions']]
    #Get rid of dates that have been deleted
    obvious_exclusions_to_export=obvious_exclusions_to_export[~obvious_exclusions_to_export['datetrack_tobeprocessed'].isnull()]
    #Replace NaNs by empty
    obvious_exclusions_to_export = obvious_exclusions_to_export.replace(np.nan, '', regex=True)
    #Save the exclusion file
    filename_flog='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/intial_selection_20172018/exclusions/obvious_exclusions_20172018.txt'
    f_log = open(filename_flog, "a")
    for i in range(0,len(obvious_exclusions_to_export['datetrack_tobeprocessed'])):
        f_log.write(str(obvious_exclusions_to_export['datetrack_tobeprocessed'].iloc[i])+' '+str(obvious_exclusions_to_export['Obvious_exclusions'].iloc[i])+'\n')
    f_log.close()
    
    
    #Extract Deletion_ablation_zone
    ablation_zone_exclusions_to_export=excel_exclusions[['datetrack_tobeprocessed','Deletion_ablation_zone']]
    #Get rid of dates that have been deleted
    ablation_zone_exclusions_to_export=ablation_zone_exclusions_to_export[~ablation_zone_exclusions_to_export['datetrack_tobeprocessed'].isnull()]
    #Replace NaNs by empty
    ablation_zone_exclusions_to_export = ablation_zone_exclusions_to_export.replace(np.nan, '', regex=True)
    #Save the exclusion file
    filename_flog='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/intial_selection_20172018/exclusions/ablation_zone_exclusions_20172018.txt'
    f_log = open(filename_flog, "a")
    for i in range(0,len(ablation_zone_exclusions_to_export['datetrack_tobeprocessed'])):
        f_log.write(str(ablation_zone_exclusions_to_export['datetrack_tobeprocessed'].iloc[i])+' '+str(ablation_zone_exclusions_to_export['Deletion_ablation_zone'].iloc[i])+'\n')
    f_log.close()
    
    
    #Extract fail_roll_correction
    fail_roll_correction_to_export=excel_exclusions[['datetrack_tobeprocessed','Exclusions_1st_run']]
    #Get rid of dates that have been deleted
    fail_roll_correction_to_export=fail_roll_correction_to_export[~fail_roll_correction_to_export['datetrack_tobeprocessed'].isnull()]
    #Replace NaNs by empty
    fail_roll_correction_to_export = fail_roll_correction_to_export.replace(np.nan, '', regex=True)
    #Save the exclusion file
    filename_flog='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/intial_selection_20172018/exclusions/fail_roll_correction_exclusions_20172018.txt'
    f_log = open(filename_flog, "a")
    for i in range(0,len(fail_roll_correction_to_export['datetrack_tobeprocessed'])):
        f_log.write(str(fail_roll_correction_to_export['datetrack_tobeprocessed'].iloc[i])+' '+str(fail_roll_correction_to_export['Exclusions_1st_run'].iloc[i])+'\n')
    f_log.close()
    
    
    #Extract dry_firn_and_other_exclusions
    fail_dry_firn_to_export=excel_exclusions[['datetrack_tobeprocessed','dry_firn_and_other_exclusions']]
    #Get rid of dates that have been deleted
    fail_dry_firn_to_export=fail_dry_firn_to_export[~fail_dry_firn_to_export['datetrack_tobeprocessed'].isnull()]
    #Replace NaNs by empty
    fail_dry_firn_to_export = fail_dry_firn_to_export.replace(np.nan, '', regex=True)
    #Save the exclusion file
    filename_flog='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/intial_selection_20172018/exclusions/dry_firn_and_other_exclusions_20172018.txt'
    f_log = open(filename_flog, "a")
    for i in range(0,len(fail_dry_firn_to_export['datetrack_tobeprocessed'])):
        f_log.write(str(fail_dry_firn_to_export['datetrack_tobeprocessed'].iloc[i])+' '+str(fail_dry_firn_to_export['dry_firn_and_other_exclusions'].iloc[i])+'\n')
    f_log.close()

pdb.set_trace()
#Compute the speed (Modified Robin speed):
# self.C / (1.0 + (coefficient*density_kg_m3/1000.0))
v= 299792458 / (1.0 + (0.734*0.873/1000.0))

#Open, read and close the file of suggested surface picks
f = open('C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/intial_selection_20172018/exclusions/datetrack_20172018.txt','r')
data_20172018 = [line.strip() for line in f.readlines() if len(line.strip()) > 0]
f.close()

#Define path where data are stored
path_data='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/'


if (obvious_identification=='TRUE'):
    
    count=0
    suggested_surface_pixel=[]
    
    #Loop over the dates of the 2017-2018 selection
    for indiv_trace in list(data_20172018):
        
        if (not(indiv_trace)=='20180430_01_139_144'):
            continue
        else:
            pdb.set_trace()
            
        #Set radar_echo_dimensions to empty
        radar_echo_dimensions=[]
    
        print(count/len(list(data_20172018))*100,' %')
        
        #Define the suite of indiv file to open
        nb_indiv_file_to_produce=int(indiv_trace[16:20])-int(indiv_trace[12:15])
        
        #Define path data to open
        path_data_open=path_data+indiv_trace[0:4]+'_Greenland_P3/CSARP_qlook/'+indiv_trace[0:11]+'/'
        
        for j in range(0,nb_indiv_file_to_produce+1):
            indiv_file_nb=int(indiv_trace[12:15])+j
            if (indiv_file_nb<10):
                fname_toload='Data_'+indiv_trace[0:11]+'_00'+str(indiv_file_nb)+'.mat'
            elif ((indiv_file_nb>=10) and (indiv_file_nb<100)):
                fname_toload='Data_'+indiv_trace[0:11]+'_0'+str(indiv_file_nb)+'.mat'
            else:
                fname_toload='Data_'+indiv_trace[0:11]+'_'+str(indiv_file_nb)+'.mat'
    
            #Open the corresponding data
            with h5py.File(path_data_open+fname_toload, 'r') as f:
                #Select radar echogram
                radar_echo=f['Data'][:].transpose() #2017 data should be transposed
                #Save horizontal dimension of radar slice
                radar_echo_dimensions=np.append(radar_echo_dimensions,radar_echo.shape[1])
                
            #Append data to each other
            if (j==0):
                #Initialize the appended radar echogram
                radar_echo_suite=radar_echo
                #Retrieve the start of surface identification
                with h5py.File(path_data_open+fname_toload, 'r') as f:
                    #Select radar echogram
                    surface_start=f['Surface'][:]
                    time_variable=f['Time'][:].transpose()
            else:
                radar_echo_suite=np.concatenate((radar_echo_suite,radar_echo),axis=1)
            #time=8373
         
        #Pick the surface
        #We can use the surface from f['Surface'][:], where the resulting is in Time
        #dimension. The time is not perfectly matching, so use where
        ind_starting_pixel=np.argmax(time_variable>surface_start[0][0])
        
        #Save surf pick for SURFACE_STARTING_PICKS_Suggestions.txt file
        suggested_surface_pixel=np.append(suggested_surface_pixel,ind_starting_pixel)
        
        #Identify the surface indices
        surface_indices=kernel_function(radar_echo_suite, ind_starting_pixel)
        
        #Compute the depths
        #self.SAMPLE_DEPTHS = self.radar_speed_m_s * self.SAMPLE_TIMES / 2.0
        depths = v * time_variable / 2.0
        
        #I.d. Select the radar slice
        #Get our slice (30 meters as currently set)
        radar_slice, bottom_indices = _return_radar_slice_given_surface(radar_echo_suite,
                                                                        depths,
                                                                        surface_indices,
                                                                        meters_cutoff_above=0,
                                                                        meters_cutoff_below=30)
        #Convert radar slice into log10
        radar_slice=np.log10(radar_slice)
        
        #Where inf in radar slice, replace by nan
        radar_slice[np.isinf(radar_slice)]=np.nan
        
        #To export slice
        slice_to_export=_export_to_8bit_array(radar_slice)
            
        #If radar_echo_dimensions larger than 1, introduce marker to differentiate between
        #the individual files
        if (len(radar_echo_dimensions)>1):
            #Get rid of the last index in radar_echo_dimensions
            radar_echo_dimensions=radar_echo_dimensions[:-1]
            #Mark the limits of the individual files by black vertical lines
            for index_to_mark in np.cumsum(radar_echo_dimensions):
                slice_to_export[:,int(index_to_mark)]=np.ones(slice_to_export.shape[0])*0
        
        '''
        #Plot the figure
        fig, (ax1) = plt.subplots()#, gridspec_kw={'width_ratios': [1, 3]})
        ax1.set_title(indiv_trace)
        ax1.imshow(radar_slice,cmap='gray')
        ax1.vlines(np.cumsum(radar_echo_dimensions), 0, radar_slice.shape[0])
        ax1.set_aspect(4)
        
        figManager = plt.get_current_fig_manager()
        figManager.window.showMaximized()
        plt.show()
        '''
        pdb.set_trace()
    
        #Save the image
        path_save_png='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/intial_selection_20172018/figures_check_iceslabs_presence/'
    
        png_to_save=png.from_array(slice_to_export, mode='L')
        png_to_save.save(path_save_png+indiv_trace+'_raw_slice.png')
        
        count=count+1

if (identification_after_depth_correction == 'TRUE'):
    count=0
    count_display=0
    
    #Define path of depth corrected
    path_depth_corrected=path_data+'exported/Depth_Corrected_Picklefiles/'
    path_boolean=path_data+'exported/Boolean Array Picklefiles/'
    
    #Loop over the dates of the 2017-2018 selection
    for indiv_trace in list(data_20172018):
        '''
        if (indiv_trace[0:4]=='2017'):
            print('2017, continue')
            continue
        '''
        if (count_display<140):
            count_display=count_display+1
            continue
        
        #Let's work with depth corrected
        print(count/len(list(data_20172018))*100,' %')
        
        #Define filename depth corrected
        filename_depth_corrected=indiv_trace+'_DEPTH_CORRECTED.pickle'
        
        #Open depth corrected pickles files
        f_depth_corrected = open(path_depth_corrected+filename_depth_corrected, "rb")
        depth_corrected_file = pickle.load(f_depth_corrected)
        f_depth_corrected.close()
        
        
        #Define the boolean filename
        filename_boolean=indiv_trace+'_SG1_CUTOFF_-0.45_THRESHOLD_000.pickle'
        #Open boolean pickles files
        f_boolean = open(path_boolean+filename_boolean, "rb")
        boolean_file = pickle.load(f_boolean)
        f_boolean.close()
        
        #Select the first 30m of the slice:
        
        #Define path data to open time variable
        path_data_open=path_data+indiv_trace[0:4]+'_Greenland_P3/CSARP_qlook/'+indiv_trace[0:11]+'/'
        #Open time variable
        with h5py.File(path_data_open+'Data_'+indiv_trace[0:15]+'.mat', 'r') as f:
                    #Select time variable
                    time_variable=f['Time'][:].transpose()        
        #calculate depth
        depths = v * time_variable / 2.0
        
        #Reset depths to 0
        depths=depths-depths[0]
        
        #Identify index where time > 30 m
        ind_lower_30m=np.where(depths<30)[0]
        depth_corrected_30m=depth_corrected_file[ind_lower_30m,:]
        
        #Plot roll corrected pickle files
        fig, (ax1,ax2) = plt.subplots(2,1)#, gridspec_kw={'width_ratios': [1, 3]})
        ax1.set_title(indiv_trace+' - first 30m')
        ax1.imshow(depth_corrected_30m,cmap='gray')
        ax1.set_aspect(4)
        
        ax2.imshow(boolean_file,cmap='gray_r')
        ax2.set_aspect(4)
        
        figManager = plt.get_current_fig_manager()
        figManager.window.showMaximized()
        plt.show()
        pdb.set_trace()
        
if (identification_dry_firn_exclusions == 'TRUE'):
    count=0
    count_display=0
    
    #Define path of depth corrected
    path_depth_corrected=path_data+'exported/Depth_Corrected_Picklefiles/'
    path_probability='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/intial_selection_20172018/figures_iceslabs_probability/pickles/'
    path_probability_old='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/intial_selection_20172018/figures_iceslabs_probability/old_version/pickles/'

    #Loop over the dates of the 2017-2018 selection
    for indiv_trace in list(data_20172018):
        
        '''
        if (indiv_trace[0:4]=='2017'):
            print('2017, continue')
            continue
        '''
        if (count_display<90):
            count_display=count_display+1
            continue
        
        #Let's work with depth corrected
        print(count/len(list(data_20172018))*100,' %')
        
        #Define filename depth corrected
        filename_depth_corrected=indiv_trace+'_DEPTH_CORRECTED.pickle'
        
        #Open depth corrected pickles files
        f_depth_corrected = open(path_depth_corrected+filename_depth_corrected, "rb")
        depth_corrected_file = pickle.load(f_depth_corrected)
        f_depth_corrected.close()
        
        
        #Define the probability filename
        filename_probability=indiv_trace+'_probability_iceslabs_presence.pickle'
        #Open probability pickles files
        f_probability = open(path_probability+filename_probability, "rb")
        probability_file = pickle.load(f_probability)
        f_probability.close()

        #Define the old probability filename
        filename_probability_old=indiv_trace+'_probability_iceslabs_presence.pickle'
        
        #Check if old probability file exist            
        if (len(glob.glob(path_probability_old+filename_probability_old))>0):
            #Open old probability pickles files
            f_probability_old = open(path_probability_old+filename_probability_old, "rb")
            probability_old_file = pickle.load(f_probability_old)
            f_probability_old.close()
        
        #Select the first 30m of the slice:
        
        #Define path data to open time variable
        path_data_open=path_data+indiv_trace[0:4]+'_Greenland_P3/CSARP_qlook/'+indiv_trace[0:11]+'/'
        #Open time variable
        with h5py.File(path_data_open+'Data_'+indiv_trace[0:15]+'.mat', 'r') as f:
                    #Select time variable
                    time_variable=f['Time'][:].transpose()        
        #calculate depth
        depths = v * time_variable / 2.0
        
        #Reset depths to 0
        depths=depths-depths[0]
        
        #Identify index where time > 30 m
        ind_lower_30m=np.where(depths<30)[0]
        depth_corrected_30m=depth_corrected_file[ind_lower_30m,:]
        
        #Plot roll corrected pickle files
        fig, (ax1,ax2,ax3) = plt.subplots(3,1)#, gridspec_kw={'width_ratios': [1, 3]})
        ax1.set_title(indiv_trace+' - depth corrected 30m')
        ax1.imshow(depth_corrected_30m,cmap='gray')
        ax1.set_aspect(4)
        
        ax2.set_title(indiv_trace+' - probability')
        ax2.imshow(probability_file,cmap='gray_r')
        ax2.set_aspect(4)
        
        if (len(glob.glob(path_probability_old+filename_probability_old))>0):
            ax3.set_title(indiv_trace+' - OLD probability')
            ax3.imshow(probability_old_file,cmap='gray_r')
            ax3.set_aspect(4)
        
        figManager = plt.get_current_fig_manager()
        figManager.window.showMaximized()
        plt.show()
        pdb.set_trace()
        


        count=count+1













 