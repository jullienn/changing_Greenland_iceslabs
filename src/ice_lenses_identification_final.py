# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 12:02:59 2021

#Description of the code:
    
Identification of ice lenses in 2002-2003 data

@author: JullienN
"""
#Import packages

import scipy.io
import rasterio
from rasterio.plot import show
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

from pysheds.grid import Grid

import osgeo.ogr as ogr
import osgeo.osr as osr

from pyproj import Transformer

import matplotlib.gridspec as gridspec
import matplotlib.animation as animation
##############################################################################
############################## Define variables ##############################
##############################################################################
#Choose the dates

#SW Greenland
date_a='jun04_02proc_4'
#date_a='may12_03_36'

#NW Greenland
#date_a='may13_03_29'
#date_a='may14_03_51'
#date_a='may18_02_0'

#SW Greenland
#date_a='may12_03_1'
#date_a='jun04_02proc_53'

#NW Greenland
#date_a='may30_02_51'
#date_a='may13_03_29'

#date_a='may30_02_51'
#date_a='may14_03_51'

#NW Greenland
#date_a='may30_02_51'
#date_a='may15_03_37'

#date_a='jun04_02proc_52'
#date_a='may12_03_1'

#date_a='may24_02_25'
#date_a='may15_03_37'

#date_a='may24_02_25'
#date_a='may15_03_37'

#Define variables
dt = 2.034489716724874e-09 #Timestep for 2002/2003 traces
t0 = 0; # Unknown so set to zero

#Compute the speed (Modified Robin speed):
# self.C / (1.0 + (coefficient*density_kg_m3/1000.0))
v= 299792458 / (1.0 + (0.734*0.873/1000.0))

#Create the dataframe that define the dates who have experienced a surface
#picking improvement
df_dates_surf_pick=pd.DataFrame({'dates_surf_pick_impr':pd.Series(['may24_02_23','may24_02_24','may24_02_25',
                                                                   'may30_02_2','may30_02_4','may30_02_5','may30_02_6',
                                                                   'may30_02_7','may30_02_13','may30_02_14','may30_02_15',
                                                                   'may30_02_50','may30_02_51'])})

plot_slice_and_improved_slice='TRUE'
display_only_potential_ice_slabs='TRUE'
technique='perc_2p5_97p5'
#perc_2p5_97p5

#N defines the number of different colors I want to use for the elevation plot
N=10
##############################################################################
############################## Define variables ##############################
##############################################################################


##############################################################################
############# Define kernel function for surface identification ##############
##############################################################################
#_gaussian function taken from IceBridgeGPR_Manager_v2.py
# Define a quick guassian function to scale the cutoff mask above
def _gaussian(x,mu,sigma):
    return np.exp(-np.power((x-mu)/sigma, 2.)/2.)

#This function have been taken from 'IceBridgeGPR_Manager_v2.py
def kernel_function(traces_input,suggested_pixel):
    
    traces = traces_input
    #Do not take the log10 of traces because 'data have been detrented in the log domain' according to John Paden's email, so I guess they are already log10!
    #traces = np.log10(traces)
    
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
    MASK_SEARCH_RADIUS = 40
    
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
##############################################################################
############# Define kernel function for surface identification ##############
##############################################################################

##############################################################################
################## Define functions for radar slice picking ##################
##############################################################################
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
##############################################################################
################## Define functions for radar slice picking ##################
##############################################################################

##############################################################################
################### Define function for radargram display ####################
##############################################################################
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
    data_cutoff_min = np.percentile(nonzero_values,  5)
    data_cutoff_max = np.percentile(nonzero_values, 95)

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
##############################################################################
################### Define function for radargram display ####################
##############################################################################

##############################################################################
############# Define function for depth correction of the traces #############
##############################################################################

def perform_depth_correction(traces_all,depths_all,surface_indices,trace_name,export, max_depth_m = 100):
    
    #We perform the depth correction over the first 100m below the surface, so
    #select the 100m slice before
    
    #Get our slice (100 meters for depth correction)
    traces, bottom_indices_100m = _return_radar_slice_given_surface(traces_all,
                                                                    depths_all,
                                                                    surface_indices,
                                                                    meters_cutoff_above=0,
                                                                    meters_cutoff_below=100)
    
    #Create an array for the depths ranging from 0 to 100m
    depths=depths_all[np.arange(0,(np.where(np.round(depths_all)==100)[-1][-1]+1))]
    
    # Use array broadcasting here.
    #pdb.set_trace()
    depths_expanded = np.zeros(traces.shape, dtype=depths.dtype)
    # Use array broadcasting to copy the depths into all the trace values
    depths.shape = depths.shape[0],1
    depths_expanded[:] = depths
    depths.shape = depths.shape[0]

    assert traces.shape == depths_expanded.shape
    #pdb.set_trace()
    # 1) Get the exponential curve fit
    def exfunc(y,A,B,C):
        return A * np.exp(B * y) + C

    popt, pcov = scipy.optimize.curve_fit(exfunc, depths_expanded.flatten(), traces.flatten(),
                                          bounds=((-np.inf, -np.inf, -np.inf),
                                                  ( np.inf,0,0)),
                                          max_nfev=1000000)

    A,B,C = popt
    print(popt)

    # Correct the traces and normalize them.
    # Original function is Z = A * e^(By) + C
    # Inverse function to normalize AND get rid of heteroscedasticitiy is 0 = ((Z - C)/A * e^(-By) - 1.0) * e^(By)
    traces_norm = ((traces - C) / A * np.exp(-B * depths_expanded) - 1.0) * np.exp(B * depths_expanded)
    # Then divide by the standard deviation of the traces to have them normalized for variance
    # All traces  for all tracks will have a MEAN of zero and a STDDEV of 1
    traces_norm = traces_norm / (np.std(traces_norm))
    
    return traces_norm
##############################################################################
############# Define function for depth correction of the traces #############
##############################################################################

##############################################################################
############### Define function for discrete colorbar display ###############
##############################################################################
def discrete_cmap(N, base_cmap=None):
    """Create an N-bin discrete colormap from the specified input map"""
    #This piece of code is from: https://gist.github.com/jakevdp/91077b0cae40f8f8244a
    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:

    base = pyplot.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)
##############################################################################
############### Define function for discrete colorbar display ###############
##############################################################################

#Open, read and close the file of suggested surface picks
f = open('C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/Exclusion_folder/txt/SURFACE_STARTING_PICKS_Suggestions_2002_2003.txt','r')
lines = [line.strip() for line in f.readlines() if len(line.strip()) > 0]
f.close()

#Open, read and close the potential ice slabs rescale file
f = open('C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/2002_2003_slice_and_loc/potential_iceslabs_rescale.txt','r')
potential_iceslabs = [line.strip() for line in f.readlines() if len(line.strip()) > 0]
f.close()

#Open the DEM
grid = Grid.from_raster("C:/Users/jullienn/Documents/working_environment/greenland_topo_data/elevations/greenland_dem_mosaic_100m_v3.0.tif",data_name='dem')
#Minnor slicing on borders to enhance colorbars
elevDem=grid.dem[:-1,:-1]              
#Scale the colormap
divnorm = mcolors.DivergingNorm(vmin=0, vcenter=1250, vmax=2500)

#Define the working environment
path= 'C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data'
os.chdir(path) # relative path: scripts dir is under Lab

# Read the years of data
folder_years = [ f.name for f in os.scandir(path) if f.is_dir() ]

for folder_year in folder_years:
    for_radar_average=[]
    
    if (folder_year in list(['2002','2003'])):
        print('Treating the year',folder_year)
        
        #Go into the yearly folders 
        folder_year_name=path+'/'+folder_year
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
                
                if (folder_day=='jun04'):
                    if (not(str(indiv_file.replace(".mat",""))==date_a)):
                        print('Not the date we are interested in, continue')
                        continue
                else:
                    if (not(str(indiv_file.replace("_aggregated",""))==date_a)):
                        print('Not the date we are interested in, continue')
                        continue
                print('Treating file',indiv_file)
                
                if (folder_day=='jun04'):
                    
                    fdata= scipy.io.loadmat(folder_day_name+'/'+indiv_file)
                    #Select radar echogram and corresponding lat/lon
                    radar_echo=fdata['data']
                    lat=fdata['latitude']
                    lon=fdata['longitude']
                    pdb.set_trace()

                else:
                    #Open the file and read it
                    f_agg = open(folder_day_name+'/'+indiv_file, "rb")
                    data = pickle.load(f_agg)
                    f_agg.close()
                                    
                    #Select radar echogram and corresponding lat/lon
                    radar_echo=data['radar_echogram']
                    
                    latlontime=data['latlontime']
                    lat=latlontime['lat_gps']
                    lon=latlontime['lon_gps']
                    
                    pdb.set_trace()
                
                #Select the first 30m of radar echogram
                #1. Compute the vertical resolution
                #a. Time computation according to John Paden's email.
                Nt = radar_echo.shape[0]
                Time = t0 + dt*np.arange(1,Nt+1)
                #b. Calculate the depth:
                #self.SAMPLE_DEPTHS = self.radar_speed_m_s * self.SAMPLE_TIMES / 2.0
                depths = v * Time / 2.0
            
                #If plot_slice_and_improved_slice is set to 'TRUE', then plot the
                #radar echogram slice and the imprived radar slice of that date and save it
                if (plot_slice_and_improved_slice=='TRUE'):
                    
                    #II.a. Plot radar track
                    #II.a.1. Reproject the track from WGS 84 to EPSG 3413
                    #Some index have lat and lon equal to 0 because of jumps in data aggregation.
                    #Replace these 0 by NaNs
                    
                    if (not(folder_day=='jun04')):
                        lat.replace(0, np.nan, inplace=True)
                        lon.replace(0, np.nan, inplace=True)
                    
                    #Transform the longitudes. The longitudes are ~46 whereas they should be ~-46! So add a '-' in front of lon
                    lon=-lon
                    
                    #Transform the coordinated from WGS84 to EPSG:3413
                    #Example from: https://pyproj4.github.io/pyproj/stable/examples.html
                    transformer = Transformer.from_crs("EPSG:4326", "EPSG:3413", always_xy=True)
                    points=transformer.transform(np.array(lon),np.array(lat))
                    
                    lon_3413=points[0]
                    lat_3413=points[1]
                    
                    if (str(indiv_file.replace("_aggregated",""))=='may12_03_36'):
                        lon_3413=np.flipud(lon_3413)
                        lat_3413=np.flipud(lat_3413)
                        
                    if (str(indiv_file.replace("_aggregated",""))=='may14_03_51'):
                        lon_3413=np.flipud(lon_3413)
                        lat_3413=np.flipud(lat_3413)                        

                    if (str(indiv_file.replace("_aggregated",""))=='may13_03_29'):
                        lon_3413=np.flipud(lon_3413)
                        lat_3413=np.flipud(lat_3413)
                    
                    if (str(indiv_file.replace(".mat",""))=='jun04_02proc_53'):
                        lon_3413=np.fliplr(lon_3413)
                        lat_3413=np.fliplr(lat_3413)
                        
                    if (str(indiv_file.replace("_aggregated",""))=='may30_02_51'):
                        lon_3413=np.flipud(lon_3413)
                        lat_3413=np.flipud(lat_3413)
                        
                    if (str(indiv_file.replace("_aggregated",""))=='may15_03_37'):
                        lon_3413=np.flipud(lon_3413)
                        lat_3413=np.flipud(lat_3413)                        
                    
                    if (str(indiv_file.replace(".mat",""))=='jun04_02proc_52'):
                        lon_3413=np.fliplr(lon_3413)
                        lat_3413=np.fliplr(lat_3413) 
                        
                    if (str(indiv_file.replace("_aggregated",""))=='may24_02_25'):
                        lon_3413=np.flipud(lon_3413)
                        lat_3413=np.flipud(lat_3413)
                    
                    
                    #II.a.2 Create the subplot
                    pyplot.figure(figsize=(48,40))
                    #Change label font
                    pyplot.rcParams.update({'font.size': 5})
                    #fig, (ax1, ax2) = pyplot.subplots(1, 2)
                    fig, (ax1, ax2) = pyplot.subplots(2, 1)#, gridspec_kw={'width_ratios': [1, 3]})

                    fig.suptitle(indiv_file.replace("_aggregated",""))
                    
                    #Subplot N°1:
                    #II.a.3. Plot dem
                    cb1=ax1.imshow(elevDem, extent=grid.extent,cmap=discrete_cmap(N,'cubehelix_r'),norm=divnorm)
                    cbar1=fig.colorbar(cb1, ax=[ax1], location='left')
                    cbar1.set_label('Elevation [m]', fontsize=5)
                    ax1.grid()
                    ax1.set_title('Radar echogram localisation',fontsize=5)
                    
                    #II.a.4. Plot the tracks
                    ax1.scatter(lon_3413, lat_3413,s=0.1)
                    
                    if (folder_day=='jun04'):
                        ax1.scatter(lon_3413[0,0],lat_3413[0,0],c='m',s=1) #Plot the start in green
                    else:
                        ax1.scatter(lon_3413[0],lat_3413[0],c='m',s=1)
                    
                    ax1.grid()
                    
                    if (folder_day=='jun04'):
                        ax1.set_xlim(lon_3413[0,0]-100000, lon_3413[0,0]+100000)
                        ax1.set_ylim(lat_3413[0,0]-100000, lat_3413[0,0]+100000)
                    else:
                        ax1.set_xlim(lon_3413[0]-100000, lon_3413[0]+100000)
                        ax1.set_ylim(lat_3413[0]-100000, lat_3413[0]+100000)
                    
                    #II.b. Plot the radar slice (first 30m of radar echogram)
                    # Load the suggested pixel for the specific date
                    for date_pix in lines:
                        if (folder_day=='jun04'):
                            if (date_pix.partition(" ")[0]==str(indiv_file.replace(".mat",""))):
                                suggested_pixel=int(date_pix.partition(" ")[2])
                                #If it has found its suggested pixel, leave the loop
                                continue   
                        else:
                            if (date_pix.partition(" ")[0]==str(indiv_file.replace("_aggregated",""))):
                                suggested_pixel=int(date_pix.partition(" ")[2])
                                #If it has found its suggested pixel, leave the loop
                                continue
                    
                    #I.b. Get the surface indices
                    
                    if (indiv_file.replace("_aggregated","") in list(df_dates_surf_pick['dates_surf_pick_impr'])):
                        #I.b.1. If already semi automatically generated, read the file
                        print(indiv_file+' have a semi-automatic improved file: use it!')
                        
                        #Construct the fiename of the wanted file
                        filename_improved_indices=[]
                        path_improved_indices='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/2002_2003_radar_slice/surf_'
                        filename_improved_indices=path_improved_indices+indiv_file+'.txt'
                        
                        #Open, read and close the file of surface picks
                        fsurf = open(filename_improved_indices,'r')
                        lines_fsurf = [line.strip() for line in fsurf.readlines() if len(line.strip()) > 0]
                        fsurf.close()
                        
                        #Store the surface indices into the right variable as int64
                        surface_indices=np.asarray(lines_fsurf,dtype=np.int64)
                        
                    else:
                        #I.b.2. If not already semi automatically generated, call
                        #the kernel_function to pick the surface
                        surface_indices=kernel_function(radar_echo, suggested_pixel)
                        
                    ##I.c. Perform depth correction
                    #depth_corrected_traces=perform_depth_correction(radar_echo, depths, surface_indices, indiv_file.replace("_aggregated",""), 'FALSE')
                    depth_corrected_traces=radar_echo
                    
                    #I.d. Select the radar slice
                    #Define the uppermost and lowermost limits
                    meters_cutoff_above=0
                    meters_cutoff_below=30
                    
                    #Redefine the 'surface_indices' variable: the surface have just been picked
                    #for the depth correction, so now we want to pick from the top down to
                    #30m depth, thus the 'surface_indices' must be [0,0,...,0]!!
                    #surface_indices=np.zeros(surface_indices.shape[0],dtype=np.int64)
                    
                    #Get our slice (30 meters as currently set)
                    radar_slice, bottom_indices = _return_radar_slice_given_surface(depth_corrected_traces,
                                                                    depths,
                                                                    surface_indices,
                                                                    meters_cutoff_above=meters_cutoff_above,
                                                                    meters_cutoff_below=meters_cutoff_below)

                        
                    #The range have already been computed, plot the data:
                    if (technique=='perc_25_75'):
                        if (folder_year=='2002'):
                            perc_lower_end=-0.08318485583215623
                            perc_upper_end=0.09414986209628376
                        elif (folder_year=='2003'):
                            perc_lower_end=-0.08488332785270308
                            perc_upper_end=0.09654050592743407
                    elif (technique=='perc_5_95'):
                        if (folder_year=='2002'):
                            perc_lower_end=-0.2870889087496134
                            perc_upper_end=0.3138722799744009
                        elif (folder_year=='2003'):
                            perc_lower_end=-0.31927843730229416
                            perc_upper_end=0.3682849426401127
                    elif (technique=='perc_05_995'):
                        if (folder_year=='2002'):
                            perc_lower_end=-2.1488917418616134
                            perc_upper_end=2.650167679823621
                        elif (folder_year=='2003'):
                            perc_lower_end=-1.661495950494564
                            perc_upper_end=1.9431298622848088
                    elif (technique=='perc_2p5_97p5'):
                        if (folder_year=='2002'):
                            perc_lower_end=-0.5709792307554173
                            perc_upper_end=0.7082634842114803
                        elif (folder_year=='2003'):
                            perc_lower_end=-0.6061610403154447
                            perc_upper_end=0.7572821079440079      
                            
                    
                    if (str(indiv_file.replace("_aggregated",""))=='may12_03_36'):
                        radar_slice=np.fliplr(radar_slice)
                    
                    if (str(indiv_file.replace("_aggregated",""))=='may14_03_51'):
                        radar_slice=np.fliplr(radar_slice)                    
                    
                    if (str(indiv_file.replace("_aggregated",""))=='may13_03_29'):
                        radar_slice=np.fliplr(radar_slice)  
                        
                    if (str(indiv_file.replace(".mat",""))=='jun04_02proc_53'):
                        radar_slice=np.fliplr(radar_slice)

                    if (str(indiv_file.replace("_aggregated",""))=='may30_02_51'):
                        radar_slice=np.fliplr(radar_slice)
                        
                    if (str(indiv_file.replace("_aggregated",""))=='may15_03_37'):
                        radar_slice=np.fliplr(radar_slice)
                        
                    if (str(indiv_file.replace(".mat",""))=='jun04_02proc_52'):
                        radar_slice=np.fliplr(radar_slice)
                        
                    if (str(indiv_file.replace("_aggregated",""))=='may24_02_25'):
                        radar_slice=np.fliplr(radar_slice)                        
                    
                    #Generate the pick for vertical distance display
                    ticks_yplot=np.arange(0,radar_slice.shape[0],20)
                    
                    ### Plot the data by changing the range
                    #Plot the radar slice
                    cb2=ax2.pcolor(radar_slice,cmap=pyplot.get_cmap('gray'))#,norm=divnorm)
                    ax2.invert_yaxis() #Invert the y axis = avoid using flipud.
                    ax2.set_aspect('equal') # X scale matches Y scale
                    #In order to display the depth, I used the example 1 of the
                    #following site: https://www.geeksforgeeks.org/matplotlib-axes-axes-set_yticklabels-in-python/
                    ax2.set_yticks(ticks_yplot) 
                    ax2.set_yticklabels(np.round(depths[ticks_yplot]))
                    ax2.set_title('Radar echogram slice, rescaling: '+technique+' percentiles',fontsize=5)
                    ax2.set_ylabel('Depth [m]')
                    ax2.set_xlabel('Horizontal distance')
                        
                    #Colorbar custom
                    cb2.set_clim(perc_lower_end,perc_upper_end)
                    cbar2=fig.colorbar(cb2, ax=[ax2], location='right')
                    cbar2.set_label('Signal strength')
                    pyplot.show()                    
                    
                    ##Create figure or radar signal evolution
                    #start_vertical_cut=370
                    #end_vertical_cut=420
                    
                    start_vertical_cut=384
                    end_vertical_cut=388
                    
                    lonely_cut=391
                    #II.a.2 Create the subplot
                    pyplot.figure(figsize=(48,40))
                    #Change label font
                    pyplot.rcParams.update({'font.size': 5})
                    fig, (ax1, ax2) = pyplot.subplots(1, 2)#, gridspec_kw={'width_ratios': [1, 3]})

                    fig.suptitle(indiv_file.replace("_aggregated",""))
                    
                    #Plot the radar slice
                    cb1=ax1.pcolor(radar_slice[:,start_vertical_cut-30:end_vertical_cut+30],cmap=pyplot.get_cmap('gray'))#,norm=divnorm)
                    
                    ax1.plot(np.ones(radar_slice[:,start_vertical_cut].size)*(start_vertical_cut-(start_vertical_cut-30)),np.arange(0,radar_slice[:,start_vertical_cut].size),color='red')
                    ax1.plot(np.ones(radar_slice[:,end_vertical_cut].size)*(end_vertical_cut-(start_vertical_cut-30)),np.arange(0,radar_slice[:,end_vertical_cut].size),color='red')

                    #ax1.plot(np.ones(radar_slice[:,lonely_cut].size)*(lonely_cut-start_vertical_cut),np.arange(0,radar_slice[:,lonely_cut].size),color='red')
                    ax1.invert_yaxis() #Invert the y axis = avoid using flipud.
                    ax1.set_aspect('equal') # X scale matches Y scale
                    #In order to display the depth, I used the example 1 of the
                    #following site: https://www.geeksforgeeks.org/matplotlib-axes-axes-set_yticklabels-in-python/
                    ax1.set_yticks(ticks_yplot) 
                    ax1.set_yticklabels(np.round(depths[ticks_yplot]))
                    ax1.set_title('Radar echogram slice - vertical cut, rescaling: '+technique+' percentiles',fontsize=5)
                    ax1.set_ylabel('Depth [m]')
                    ax1.set_xlabel('Horizontal distance')
                        
                    #Colorbar custom
                    cb1.set_clim(perc_lower_end,perc_upper_end)
                    cbar1=fig.colorbar(cb1, ax=[ax1], location='right')
                    cbar1.set_label('Signal strength')
                    
                    #Plot the radar strength evolution
                    #ax2.plot(radar_slice[:,start_vertical_cut+16],np.arange(0,radar_slice[:,vertical_cut].size))#,norm=divnorm)
                    #ax2.plot(radar_slice[:,lonely_cut],np.arange(0,radar_slice[:,lonely_cut].size),color='red')
                    
                    #If several vertical lines to plot and plot average
                    for i in range (start_vertical_cut,end_vertical_cut+1):
                        print(i)
                        ax2.plot(radar_slice[:,i],np.arange(0,radar_slice[:,i].size),color='0.7')
                    ax2.plot(np.mean(radar_slice[:,start_vertical_cut:end_vertical_cut],1),np.arange(0,radar_slice[:,i].size),color='red')
                    
                    #ax2.set_ylim(0,30)
                    
                    #max_plot=np.maximum(np.absolute(np.min(radar_slice[:,lonely_cut])),np.absolute(np.max(radar_slice[:,lonely_cut])))
                    #max_plot=max_plot+0.1*max_plot
                    #ax2.set_xlim(-max_plot,max_plot)
                    ax2.invert_yaxis() #Invert the y axis = avoid using flipud.
                    ax2.set_yticks(ticks_yplot) 
                    ax2.set_yticklabels(np.round(depths[ticks_yplot]))                    
                    ax2.set_title('Radar signal strength',fontsize=5)
                    ax2.set_ylabel('Depth [m]')
                    ax2.set_xlabel('Radar signal strength')
                    ax2.grid()
                    
                    pdb.set_trace()
                    
                    #Create an animated plot
                    
                    fig, ax = plt.subplots()
                    
                    
                    for i in range (start_vertical_cut,end_vertical_cut+1):
                        print(i)
                        ax2.plot(radar_slice[:,i],np.arange(0,radar_slice[:,i].size),color='0.7')
                                        
                    ax2.plot(np.mean(radar_slice[:,start_vertical_cut:end_vertical_cut],1),np.arange(0,radar_slice[:,i].size),color='red')
                    
                    x = np.arange(0, 2*np.pi, 0.01)
                    line, = ax.plot(x, np.sin(x))
                    
                    
                    def animate(i):
                        line.set_ydata(np.sin(x + i / 50))  # update the data.
                        return line,
                    
                    
                    ani = animation.FuncAnimation(
                    fig, animate, interval=20, blit=True, save_count=50)

                    
                    
print('End of processing')