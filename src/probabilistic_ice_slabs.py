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
    #pdb.set_trace()
    data_cutoff_min = 0#np.percentile(nonzero_values,  0.5)
    data_cutoff_max = 1#np.percentile(nonzero_values, 99.5)

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

path_quantiles_data='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2010_2018/ii_out_from_iceslabs_processing_jullien.py/pickles/'
'''
path_quantiles_data='/flash/jullienn/data/threshold_processing_output/pickles/'
'''
#Identify all the datetraces to process

path_datetrack='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/'
'''
path_datetrack='/flash/jullienn/data/threshold_processing/'
'''
datetrack_toread = np.asarray(pd.read_csv(path_datetrack+'datetrack_20102018.txt', header=None))

#Define the desired quantiles over which we will loop
desired_quantiles=np.round(np.arange(0.63,0.82,0.01),2)

#intialize counter to 0
count_time=0

#II. Loop over these datetracks, and perform probability calculation:
for indiv_trace in datetrack_toread:
    
    '''
    #pdb.set_trace()
    #If pickle files have already been created, do not process and continue
    filename_to_check='/flash/jullienn/data/threshold_processing_output/probability_iceslabs/pickles/'+indiv_trace[0]+'*'
    
    if (len(glob.glob(filename_to_check))>0):
        print(indiv_trace[0],': files already existent, move on to the next date')
        continue
    '''
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
    
        
    #Divide the probabilistic_slice by the number of quantiles to have a probability map
    probabilistic_slice=probabilistic_slice/len(desired_quantiles)
    
    #Save the image
    
    fig_name='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2010_2018/iii_out_from_probabilistic_iceslabs.py/images/'+indiv_trace[0]+'_probability_iceslabs_presence.png'
    '''
    fig_name='/flash/jullienn/data/threshold_processing_output/probability_iceslabs/images/'+indiv_trace[0]+'_probability_iceslabs_presence.png'
    '''
    
    '''
    #Traditional was of plotting, depreciated here
    fig, (ax1) = plt.subplots(1, 1)
    
    #Plot custom threshold ice slabs identification
    cb=ax1.imshow(final_probability_slice,cmap=plt.get_cmap('Blues'))#,norm=divnorm)
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
    ax2.imshow(final_probability_slice,cmap=plt.get_cmap('Blues'))
    plt.show()
    '''
    
    if (indiv_trace not in list(['20110416_01_053_055'])):
        continue
    else:
        pdb.set_trace()
    
    #Prepare matrix for png plot. (1-probabilistic_slice) because 1 is white
    #out of the function _export_to_8bit_array, and I want black
    probabilistic_slice_png=_export_to_8bit_array((1-probabilistic_slice))
    
    #Save the image
    #pdb.set_trace()
    png_to_save=png.from_array(probabilistic_slice_png, mode='L')
    png_to_save.save(fig_name)
    
    #Save the pickle file
    
    filename_tosave='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2010_2018/iii_out_from_probabilistic_iceslabs.py/pickles/'+indiv_trace[0]+'_probability_iceslabs_presence.pickle'
    '''
    filename_tosave='/flash/jullienn/data/threshold_processing_output/probability_iceslabs/pickles/'+indiv_trace[0]+'_probability_iceslabs_presence.pickle'
    '''
    outfile= open(filename_tosave, "wb" )
    pickle.dump(probabilistic_slice,outfile)
    outfile.close()

print('End of probabilistic processing')






##############################################################################
###              Generate en excel file of ice slabs thickness             ###
##############################################################################
# This is from IceBridgeGPR_Manager_v2.py



def export_ice_layer_lat_lon_distance_thicknesses(self):
    '''Export to a CSV, ice layer lat,lon,& thicknesses.  Omit all zero values.'''

    tracks = self.compile_icebridge_tracks_with_ice_lenses()
    #Once I am out of this, I just store the suite of commands to execute in tracks but no actual data sotred in it

    fout = open(ICEBRDIGE_ICE_LAYER_OUTPUT_CSV_FILE, 'w')
    #pdb.set_trace()
    header = "Track_name,Tracenumber,lat,lon,alongtrack_distance_m,20m_ice_content_m\n"
    fout.write(header)
            
    for track in tracks:
        
        print(track.NAME, end=' ')
        '''
        if (not(track.NAME == '20170511_01_010_025')):
            continue
        else:
            pdb.set_trace()
        '''
        lats, lons, distances, ice_contents = track.return_ice_layers_lat_lon_distance_thickness(masked=False)
        # The one "really long" track has artifacts in the center that aren't real ice layers.  Filter these out.
        if track.NAME == "20120412_01_095_095":
            ice_contents[0:9000] = 0.0

        assert len(lats) == len(lons) == len(ice_contents)
        tracenums = numpy.arange(len(lats), dtype=numpy.int)

        tracecount = 0
        for lat, lon, tracenum, distance, ice_content in zip(lats, lons, tracenums, distances, ice_contents):
            # Record ONLY traces that have > 1 m ice content in them.  We're not interested in thinner stuff here.
            # If it has > 16 m ice content (80%), we also omit it, just to keep pure ice out of it.
            if 1.0 <= ice_content <= 16.0:
                line = "{0},{1},{2},{3},{4},{5}\n".format(track.NAME, tracenum, lat, lon, distance, ice_content)
                fout.write(line)
                tracecount += 1
        print(tracecount, "of", len(lats), "traces.")
        print()

    fout.close()
    #pdb.set_trace()
    print("Exported", os.path.split(ICEBRDIGE_ICE_LAYER_OUTPUT_CSV_FILE)[-1])


ICEBRIDGE_ICELENS_QUICKLOOK_FOLDER = r'C:\Users\jullienn\Documents\working_environment\iceslabs_MacFerrin\data\2018_Greenland_P3\images'

def compile_icebridge_tracks_with_ice_lenses(self, quicklook_directory = ICEBRIDGE_ICELENS_QUICKLOOK_FOLDER):

    '''Similar to "::compile_icebridge_tracks()", this will complile a dictionary of IceBridge_Track objects.
    However, instead of using all tracks and all files within the tracks, this will peruse a directory containing only files that
    have been flagged to be part of IceBridge files that may (or do) contain ice lenses in them.

    The directory contains only the quicklook PDF files that visualize the tracks.
    '''
    track_count_N = 0

    # 1) Open the directory, get a list of all the files in there
    # 2) Peruse files, separate out "map" files from "echo" files, use just the echo files.
    all_files = [f for f in os.listdir(quicklook_directory) if f.find("1echo.jpg") != -1]

    # 3) Compile lists of all sequential adjoining files (part of the same sequence) that will create tracks.
    track_list_dict = {}
    current_list = []
    for f in all_files:
        # If we have a new list, OR if the last file listed in the current list has the
        # same FILE_ID, and the file number is one greater than the last one, append it onto the list.
        if len(current_list) == 0 or \
            (f[0:11] == current_list[-1][0:11] and (int(f[12:15]) - int(current_list[-1][12:15])) == 1):

            current_list.append(f)
        else:
            # Convert the track_id substring into a TRACK_ID integer, append it.
            track_id = current_list[0][0:11] + current_list[0][11:15] + current_list[-1][11:15]
            # Stick it into the dictionary
            track_list_dict[track_id] = current_list

            # Create a new "current list" with the new file in it.
            current_list = [f]
            track_count_N += 1

    # Make sure the very last one gets on there.
    track_list_dict[current_list[0][0:11] + current_list[0][11:15] + current_list[-1][11:15]] = current_list
    track_count_N += 1
    # We have just created a list where the track id are stores

    if self.verbose:
        print(track_count_N, "tracks from", len(all_files), "files.\n")


    # 4) Create a (subset) track from each file.
    self.track_names = list(track_list_dict.keys())
    self.track_names.sort()

    tracks = [IceBridgeGPR_Track_v2(self.h5file, name, verbose=self.verbose) for name in self.track_names]

    # 5) Save the dictionary.
    self.tracks = tracks

    return tracks
    


def return_ice_layers_lat_lon_distance_thickness(self, masked=False):

    '''Once we have boolean ice layers calculated, return the latitude,
    longitude, elevation, and ice thickness for each trace.  If masked=True,
    return them masked out.  If masked=False, don't bother masking them.'''

    lats, lons = self.return_coordinates_lat_lon()
    # So far I have read the data and stored the lat and lon coordinates
    boolean_traces = self.get_processed_traces(datatype="boolean_ice_layers") ====>>> This is the SG+ cutoff -0.45, continuits 350 file 

    depths = self.get_sample_depths(trace_array = boolean_traces)
    depth_delta_m = numpy.mean(depths[1:] - depths[:-1])
    distances = numpy.cumsum(self.compute_distances())
    distances = numpy.append([0.0], distances)

    # Number of pixels times the thickness of each pixel
    ice_content_m = numpy.sum(boolean_traces, axis=0) * depth_delta_m

    if masked:
        mask = self.get_boolean_ice_mask()
        lats = lats[mask]
        lons = lons[mask]
        distances = distances[mask]
        ice_content_m = ice_content_m[mask]

    return lats, lons, distances, ice_content_m

ICEBRDIGE_ICE_LAYER_OUTPUT_CSV_FILE=ICEBRDIGE_ICE_LAYER_OUTPUT_CSV_FILE = os.path.join(ICEBRIDGE_EXPORT_FOLDER,"txt\Ice_Layer_Output_Thicknesses_2018.csv")


##############################################################################
###              Generate en excel file of ice slabs thickness             ###
##############################################################################
