# -*- coding: utf-8 -*-
"""
Created on Sat Dec  5 11:34:11 2020

@author: Jullien Nicolas
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

##############################################################################
############################## Define variables ##############################
##############################################################################
dt = 2.034489716724874e-09 #Timestep for 2002/2003 traces
t0 = 0; # Unknown so set to zero

#Compute the speed (Modified Robin speed):
# self.C / (1.0 + (coefficient*density_kg_m3/1000.0))
v= 299792458 / (1.0 + (0.734*0.873/1000.0))


#In the following sections, we do perform depth correction of the traces
plot_radar_echogram_slice='FALSE'
plot_slice_and_loc='FALSE'

#In the following sections, we DO NOT perform depth correction of the traces
surf_pick_selection='FALSE'
raw_radar_echograms='FALSE'
plot_radar_loc='FALSE'
plot_original_slice_and_cutted_slice='FALSE'
plot_slice_and_loc_rescaled='TRUE'

#Which resclaing procedure we want. Choose among:
#perc_25_75, perc_5_95, perc_2p5_97p5, perc_05_995
technique='perc_5_95'

#N defines the number of different colors I want to use for the elevation plot
N=10

#Create the dataframe that define the dates who have experienced a surface
#picking improvement
df_dates_surf_pick=pd.DataFrame({'dates_surf_pick_impr':pd.Series(['may24_02_23','may24_02_24','may24_02_25',
                                                                   'may30_02_2','may30_02_4','may30_02_5','may30_02_6',
                                                                   'may30_02_7','may30_02_13','may30_02_14','may30_02_15',
                                                                   'may30_02_50','may30_02_51'])})
 
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

def _return_radar_slice_given_surface_below_surf(traces,
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
    
    output_traces = np.empty((-idx_above + idx_below, traces.shape[1]), dtype=traces.dtype)
    bottom_indices = np.zeros(shape=(1,traces.shape[1]))
    
    for i,s in enumerate(surface_indices):
        try:
            output_traces[:,i] = traces[(s+idx_above):(s+idx_below), i]
            bottom_indices[0,i]=(s+idx_below)
        except ValueError:
            print('Should not come here')
            # If the surf_i is too close to one end of the array or the other, it extends beyond the edge of the array and breaks.
            if s < idx_above:
                start, end = None, idx_above+idx_below
            elif s > (traces.shape[0] + idx_below):
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

    if (export=='TRUE'):
        ###################################################
        ## Depth-correction and normalization PLOT
        ###################################################
        # We don't need to plot all the traces, just a subset (100,000 will do)
        if traces.size > 100000:
            # Subset to only plot 100000 (?) of the points
            gap = int(traces.size / 100000)
            traces_subset = traces.flatten()[::gap]
            # Contract the variability of the points to have ~ the same variability as the original points, for display only
            norm_subset = (traces_norm.flatten()/4.0)[::gap]
            depths_subset = depths_expanded.flatten()[::gap]
        else:
            traces_subset = traces.flatten()
            norm_subset = (traces_norm/4.0).flatten()
            depths_subset = depths_expanded.flatten()

        curve_fit_y = exfunc(depths, *popt)
        # 2) Make a plot, save it.
        fig = pyplot.figure(figsize=(5,3))
        # Plot the corrected points below, in pink/red
        pyplot.plot(depths_subset, norm_subset, "o", ms=1, color="salmon", fillstyle="full", mec="salmon")
        pyplot.axhline(y=0,color="red",ls="--",label="corrected")

        # Plot the original points atop, in blue
        pyplot.plot(depths_subset, traces_subset, "o", ms=1, color="lightblue", fillstyle="full", mec="lightblue")
        pyplot.plot(depths, curve_fit_y, color="blue",label="uncorrected")

        ax = fig.axes[0]

        equation_string = "$\Omega(y) = {0:0.3f} ".format(A) + "\cdot e^{" + "{0:0.5f}\cdot y".format(B) + "}" + "{0:0.3f}$".format(C)
        pyplot.text(0.04,0.10,equation_string,
                 horizontalalignment="left",
                 verticalalignment="center",
                 transform=ax.transAxes)

        # Plot legend
        handles, labels = ax.get_legend_handles_labels()
        # Even thought we plotted the corrected first, put the uncorrected first in the legend
        handles = handles[::-1]
        labels = labels[::-1]
        ax.legend(handles, labels, loc="upper right", fontsize="x-small", markerscale=0.70)

        # Title and axis labels
        pyplot.title(trace_name)
        pyplot.xlabel("Depth $y$ (m)")
        pyplot.ylabel("GPR $\Omega$ (dB)")

        # Begin: Added on September 16, 2020 to fit MacFerrins' figures
        #ax.set_xlim(0,100)
        #ax.set_ylim(-8,2)
        # End: Added on September 16, 2020 to fit MacFerrins' figures

        pyplot.tight_layout()
        figname = os.path.join('C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/depth_correction', trace_name + "_DEPTH_CURVE_PLOT.png")
        pyplot.savefig(figname, dpi=600)
        print("Exported", os.path.split(figname)[1])
        pyplot.cla()
        pyplot.close()

        #######################################
        ### Export picklefile
        #######################################
        #traces_norm_inflated = self._refill_array(traces_norm, mask)
        #
        #f = open(self.FNAME_depth_corrected_picklefile, 'wb')
        #pickle.dump(traces_norm_inflated, f)
        #f.close()
        #print("Exported", os.path.split(self.FNAME_depth_corrected_picklefile)[-1])
        #
        ## Save to object
        #self.TRACES_depth_corrected = traces_norm_inflated

        #######################################
        ### Export corrected image
        #######################################
        #cutoff_30m = depths[(depths <= 30.0)].size
        #traces_export = traces_norm_inflated[:cutoff_30m, :]
        #self.export_image(traces_export,"_XDEPTHCORRECT_AFTER")

    # 3) Return depth-correction parameters
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

#Open the DEM
grid = Grid.from_raster("C:/Users/jullienn/Documents/working_environment/greenland_topo_data/elevations/greenland_dem_mosaic_100m_v3.0.tif",data_name='dem')
#Minnor slicing on borders to enhance colorbars
elevDem=grid.dem[:-1,:-1]              
#Scale the colormap
divnorm = mcolors.DivergingNorm(vmin=0, vcenter=1250, vmax=2500)

#Open, read and close the file of suggested surface picks
f = open('C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/Exclusion_folder/txt/SURFACE_STARTING_PICKS_Suggestions_2002_2003.txt','r')
lines = [line.strip() for line in f.readlines() if len(line.strip()) > 0]
f.close()

#Open, read and close the potential ice slabs rescale file
f = open('C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/2002_2003_slice_and_loc/potential_iceslabs_rescale.txt','r')
potential_iceslabs = [line.strip() for line in f.readlines() if len(line.strip()) > 0]
f.close()

#Define the working environment
path= 'C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data'
os.chdir(path) # relative path: scripts dir is under Lab

# Read the years of data
folder_years = [ f.name for f in os.scandir(path) if f.is_dir() ]

for folder_year in folder_years:
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
                print('Treating file',indiv_file)

                #pdb.set_trace()
                #If indiv_file is the quality file, continue
                if (indiv_file[0:7]==('quality')):
                    #pdb.set_trace()
                    continue
                
                if (folder_day=='jun04'):
                    
                    fdata= scipy.io.loadmat(folder_day_name+'/'+indiv_file)
                    #Select radar echogram and corresponding lat/lon
                    radar_echo=fdata['data']
                    lat=fdata['latitude']
                    lon=fdata['longitude']
                    #pdb.set_trace()

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
                
                #Select the first 30m of radar echogram
                #1. Compute the vertical resolution
                #a. Time computation according to John Paden's email.
                Nt = radar_echo.shape[0]
                Time = t0 + dt*np.arange(1,Nt+1)
                #b. Calculate the depth:
                #self.SAMPLE_DEPTHS = self.radar_speed_m_s * self.SAMPLE_TIMES / 2.0
                depths = v * Time / 2.0
                
                #If raw_radar_echograms is set to 'TRUE', then plot the raw
                #radar echogram of that date and save it
                if (raw_radar_echograms=='TRUE'):
                    #If file have already been created, continue
                    filename_to_check='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/2002_2003_radar_raw_echogram/'+indiv_file+'.png'
                    if (os.path.isfile(filename_to_check)):
                        print('Figure already existent, move on to the next date')
                        continue
                    
                    #Generate the pick for vertical distance display
                    ticks_yplot=np.arange(0,radar_echo.shape[0],200)
                    
                    #Plot the raw radar echogram
                    pyplot.figure(figsize=(48,40))
                    
                    #Change label font
                    pyplot.rcParams.update({'font.size': 40})
                    #pdb.set_trace()
                    #pyplot.figure()
                    color_map=pyplot.pcolor(radar_echo,cmap=pyplot.get_cmap('gray'))#,norm=divnorm)
                    pyplot.gca().invert_yaxis() #Invert the y axis = avoid using flipud.
                    pyplot.yticks(ticks=ticks_yplot,labels=(np.round(depths[ticks_yplot])))
                    pyplot.ylabel('Depth [m]')
                    pyplot.xlabel('Horizontal distance')
                    pyplot.title('Raw radar echogram: '+indiv_file.replace("_aggregated",""))
                    cbar=pyplot.colorbar()
                    cbar.set_label('Signal strength')
                    
                    #Create the figure name
                    fig_name=[]
                    fig_name='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/2002_2003_radar_raw_echogram/'+indiv_file+'.png'
                    
                    #Save the figure
                    pyplot.savefig(fig_name)
                    pyplot.clf()
                    
                    #pdb.set_trace()
                    
                    continue
                
                #If plot_radar_echogram_slice is set to 'TRUE', then plot the slice
                #radar echogram of that date and save it
                if (plot_radar_echogram_slice=='TRUE'):
                    #If file have already been created, continue
                    filename_to_check='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/2002_2003_radar_slice/'+indiv_file+'.png'
                    if (os.path.isfile(filename_to_check)):
                        print('Figure already existent, move on to the next date')
                        continue
                    
                    #I. Process and plot radar echogram
                    #I.a. Load the surface suggestion pick (there is no 'Surface'
                    # variable in 2002/2003 dataset such as 2010/2014 datset).

                    # Load the suggested pixel for the specific date
                    #pdb.set_trace()
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
                    #pdb.set_trace()
                    
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
                    
                    #I.c. Perform depth correction
                    depth_corrected_traces=perform_depth_correction(radar_echo, depths, surface_indices, indiv_file.replace("_aggregated",""), 'FALSE')

                    #I.d. Select the radar slice
                    #Define the uppermost and lowermost limits
                    meters_cutoff_above=0
                    meters_cutoff_below=30
                    
                    #Redefine the 'surface_indices' variable: the surface have just been picked
                    #for the depth correction, so now we want to pick from the top down to
                    #30m depth, thus the 'surface_indices' must be [0,0,...,0]!!
                    surface_indices=np.zeros(surface_indices.shape[0],dtype=np.int64)
                    
                    #Get our slice (30 meters as currently set)
                    radar_slice, bottom_indices = _return_radar_slice_given_surface(depth_corrected_traces,
                                                                    depths,
                                                                    surface_indices,
                                                                    meters_cutoff_above=meters_cutoff_above,
                                                                    meters_cutoff_below=meters_cutoff_below)
                    
                    # I have taken and adatped the functions '_return_radar_slice_given_surface' and
                    # '_radar_slice_indices_above_and_below' from 'IceBridgeGPR_Manager_v2.py'
                    # and it seems to correctly selecting the slice! I did not manually check
                    # by looking in the variables it the job was done correctly but I have
                    # checked several variables such as idx_above, idx_below, output traces
                    # and it seems okay to me!
    
                    ##############################################################
                    ############### Begin explanations on pcolor #################
                    
                    #Explainations on how pcolor works. I convinced myself doing a small example that I plotted.
                    #Further explanations: https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.pcolor.html#:~:text=pcolor()%20displays%20all%20columns,Similarly%20for%20the%20rows.
                    #Example: I have a matrix Z that I want to plot
                    
                    #Z=[10,3,24,70,                            40,5,48,22
                    #    2,6,87,21,      ----pcolor(Z)--->     2,6,87,21
                    #    40,5,48,22]                           10,3,24,70
                    
                    #I must use np.flipud(), or invert the y axis to display the data from top to bottom.
                    
                    ################ End explanations on pcolor ##################
                    ##############################################################
                    
                    #Generate the pick for vertical distance display
                    ticks_yplot=np.arange(0,radar_slice.shape[0],20)
                    #pdb.set_trace()
                    #I.d. Plot the radar slice (first 30m of radar echogram)
                    #pdb.set_trace()
                    #Plot the data            
                    
                    pyplot.figure(figsize=(40,10))
                    
                    #Change label font
                    pyplot.rcParams.update({'font.size': 40})
                    
                    color_map=pyplot.pcolor(radar_slice,cmap=pyplot.get_cmap('gray'))#,norm=divnorm)
                    pyplot.yticks(ticks=ticks_yplot,labels=(np.round(depths[ticks_yplot])))
                    pyplot.gca().invert_yaxis() #Imvert the y axis = avoid using flipud.
                    pyplot.gca().set_aspect('equal') # X scale matches Y scale
                    pyplot.ylabel('Depth [m]')
                    pyplot.xlabel('Horizontal distance')
                    #pyplot.yticks(ticks=ticks_yplot,labels=labels_yplot)
                    #pyplot.xticks(fontsize=20)
                    #pyplot.yticks(fontsize=20)
                    #pyplot.ylim(0, 200)
                    pyplot.title('Radar echogram slice: '+indiv_file.replace("_aggregated",""))

                    cbar=pyplot.colorbar()
                    cbar.set_label('Signal strength')
                    #pyplot.show()
                    
                    #Create the figure name
                    fig_name=[]
                    fig_name='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/2002_2003_radar_slice/'+indiv_file+'.png'
                    
                    #Save the figure
                    pyplot.savefig(fig_name,dpi=500)
                    pyplot.clf()
                    
                    continue
                
                #If plot_radar_loc is set to 'TRUE', then plot the location of
                #radar echogram of that date and save it
                if (plot_radar_loc=='TRUE'):
                    #pdb.set_trace()
                    #If file have already been created, continue
                    filename_to_check='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/2002_2003_data_localisation/'+indiv_file+'.png'
                    if (os.path.isfile(filename_to_check)):
                        print('Figure already existent, move on to the next date')
                        continue
                    
                    #II. Plot radar echogram localisation
                    #II.a. Plot dem
                    pyplot.figure(figsize=(48,40))
                    
                    #Change label font
                    pyplot.rcParams.update({'font.size': 40})
                    
                    pyplot.imshow(elevDem, extent=grid.extent,cmap='hot_r',norm=divnorm)
                    pyplot.colorbar(label='Elevation [m]')
                    pyplot.grid()
                    pyplot.title('Radar echogram location: '+indiv_file.replace("_aggregated",""))
   
                    #II.b. Plot radar track
                    #II.b.1. Reproject the track from WGS 84 to EPSG 3413
                    #Some index have lat and lon equal to 0 because of jumps in data aggregation.
                    #Replace these 0 by NaNs
                    
                    #pdb.set_trace()
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
                    
                    #II.b.2. Plot the tracks
                    pyplot.scatter(lon_3413, lat_3413)
                    #pyplot.scatter(lon_3413[0],lat_3413[0],c='m') #Plot the start in green
                    pyplot.grid()
                    #pyplot.show()
                    
                    #Create the figure name
                    fig_name=[]
                    fig_name='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/2002_2003_data_localisation/'+indiv_file+'.png'
                    
                    #Save the figure
                    pyplot.savefig(fig_name)
                    pyplot.clf()
                    
                    continue
                
                #If plot_slice_and_loc is set to 'TRUE', then plot the location of
                #radar echogram AND the radar slice of that date and save it
                if (plot_slice_and_loc=='TRUE'):
                    ##If file have already been created, continue
                    #filename_to_check='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/2002_2003_slice_and_loc/'+indiv_file+'.png'
                    #if (os.path.isfile(filename_to_check)):
                    #    print('Figure already existent, move on to the next date')
                    #    continue
                                        
                    #Subplot N°1:
                    #I. Process and plot radar echogram
                    #I.a. Load the surface suggestion pick (there is no 'Surface'
                    # variable in 2002/2003 dataset such as 2010/2014 datset).
                    
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
                        
                    #I.c. Perform depth correction
                    depth_corrected_traces=perform_depth_correction(radar_echo, depths, surface_indices, indiv_file.replace("_aggregated",""), 'FALSE')

                    #I.d. Select the radar slice
                    #Define the uppermost and lowermost limits
                    meters_cutoff_above=0
                    meters_cutoff_below=30
                    
                    #Redefine the 'surface_indices' variable: the surface have just been picked
                    #for the depth correction, so now we want to pick from the top down to
                    #30m depth, thus the 'surface_indices' must be [0,0,...,0]!!
                    surface_indices=np.zeros(surface_indices.shape[0],dtype=np.int64)
                    
                    #Get our slice (30 meters as currently set)
                    radar_slice, bottom_indices = _return_radar_slice_given_surface(depth_corrected_traces,
                                                                    depths,
                                                                    surface_indices,
                                                                    meters_cutoff_above=meters_cutoff_above,
                                                                    meters_cutoff_below=meters_cutoff_below)
                    # I have taken and adatped the functions '_return_radar_slice_given_surface' and
                    # '_radar_slice_indices_above_and_below' from 'IceBridgeGPR_Manager_v2.py'
                    # and it seems to correctly selecting the slice! I did not manually check
                    # by looking in the variables it the job was done correctly but I have
                    # checked several variables such as idx_above, idx_below, output traces
                    # and it seems okay to me!
                    
                    #Rescale the radar slice as MacFerrin et al. 2019
                    #1.If required to go through _export_to_8bit_array as a vector
                    #radar_slice_rescaled=_export_to_8bit_array(np.ndarray.flatten(radar_slice))
                    #radar_slice_rescaled_mat = radar_slice_rescaled.reshape(radar_slice.shape[0],radar_slice.shape[1])
                    
                    #2.If not required to go through _export_to_8bit_array as a vector
                    radar_slice_rescaled_mat=_export_to_8bit_array(radar_slice)
    
                    ##############################################################
                    ############### Begin explanations on pcolor #################
                    
                    #Explainations on how pcolor works. I convinced myself doing a small example that I plotted.
                    #Further explanations: https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.pcolor.html#:~:text=pcolor()%20displays%20all%20columns,Similarly%20for%20the%20rows.
                    #Example: I have a matrix Z that I want to plot
                    
                    #Z=[10,3,24,70,                            40,5,48,22
                    #    2,6,87,21,      ----pcolor(Z)--->     2,6,87,21
                    #    40,5,48,22]                           10,3,24,70
                    
                    #I must use np.flipud(), or invert the y axis to display the data from top to bottom.
                    
                    ################ End explanations on pcolor ##################
                    ##############################################################
                    
                    #II. Plot radar echogram localisation
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
                        ax1.scatter(lon_3413[0,0],lat_3413[0,0],c='m',s=0.1) #Plot the start in green
                    else:
                        ax1.scatter(lon_3413[0],lat_3413[0],c='m',s=0.1)
                    
                    ax1.grid()
                    
                    if (folder_day=='jun04'):
                        ax1.set_xlim(lon_3413[0,0]-500000, lon_3413[0,0]+500000)
                        ax1.set_ylim(lat_3413[0,0]-500000, lat_3413[0,0]+500000)
                    else:
                        ax1.set_xlim(lon_3413[0]-500000, lon_3413[0]+500000)
                        ax1.set_ylim(lat_3413[0]-500000, lat_3413[0]+500000)
                    
                    #II.b. Plot the radar slice (first 30m of radar echogram)
                    
                    #Subplot N°2:
                    #Create the y vector for plotting
                    ticks_yplot=np.arange(0,radar_slice.shape[0],20)
                    
                    #Plot the radar slice
                    cb=ax2.pcolor(radar_slice_rescaled_mat,cmap=pyplot.get_cmap('gray'))#,norm=divnorm)
                    ax2.invert_yaxis() #Invert the y axis = avoid using flipud.
                    ax2.set_aspect('equal') # X scale matches Y scale
                    #In order to display the depth, I used the example 1 of the
                    #following site: https://www.geeksforgeeks.org/matplotlib-axes-axes-set_yticklabels-in-python/
                    ax2.set_yticks(ticks_yplot) 
                    ax2.set_yticklabels(np.round(depths[ticks_yplot]))
                    ax2.set_title('Radar echogram slice, rescaled from 0 to 256 - 5-95%',fontsize=5)
                    ax2.set_ylabel('Depth [m]')
                    ax2.set_xlabel('Horizontal distance')
                    #cbar=fig.colorbar(cb)
                    #cbar.set_label('Signal strength', fontsize=5)
                    #fig.tight_layout()
                    pyplot.show()
                    
                    #Create the figure name
                    #fig_name=[]
                    #fig_name='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/2002_2003_slice_and_loc/'+indiv_file+'_5_95_after_depth_corr.png'
                    
                    ##Save the figure
                    #pyplot.savefig(fig_name,dpi=500)
                    #pyplot.clf()
                    ##Plot the data
                    pdb.set_trace()
                    
                    continue
                
                #If surf_pick_selection is set to 'TRUE', then plot the raw radar
                # echogram with the surface overlayed AND the radar slice of
                #that date and save it
                if (surf_pick_selection=='TRUE'):
                    #If file have already been created, continue
                    filename_to_check='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/2002_2003_raw_and_slice/'+indiv_file+'.png'
                    if (os.path.isfile(filename_to_check)):
                        print('Figure already existent, move on to the next date')
                        continue
                    
                    #I. Process and radar echogram
                    #I.a. Load the surface suggestion pick (there is no 'Surface'
                    # variable in 2002/2003 dataset such as 2010/2014 datset).

                    # Load the suggested pixel for the specific date
                    #pdb.set_trace()
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
                
                    #I.c. Select the radar slice
                    #Define the uppermost and lowermost limits
                    meters_cutoff_above=0
                    meters_cutoff_below=30
                    
                    #Get our slice (30 meters as currently set)
                    radar_slice, bottom_indices = _return_radar_slice_given_surface(radar_echo,
                                                                    depths,
                                                                    surface_indices,
                                                                    meters_cutoff_above=meters_cutoff_above,
                                                                    meters_cutoff_below=meters_cutoff_below)
                    # I have taken and adatped the functions '_return_radar_slice_given_surface' and
                    # '_radar_slice_indices_above_and_below' from 'IceBridgeGPR_Manager_v2.py'
                    # and it seems to correctly selecting the slice! I did not manually check
                    # by looking in the variables it the job was done correctly but I have
                    # checked several variables such as idx_above, idx_below, output traces
                    # and it seems okay to me!
                    
                    #Create the y vector for plotting
                    ticks_yplot=np.arange(0,radar_echo.shape[0],200)
                    
                    #II. Plot radar echogram localisation
                    #II.a. Create the subplot
                    #pdb.set_trace()
                    pyplot.figure(figsize=(48,40))
                    pyplot.rcParams.update({'font.size': 5})
                    fig, (ax1, ax2) = pyplot.subplots(2, 1)#, gridspec_kw={'width_ratios': [1, 3]})
                    fig.suptitle(indiv_file.replace("_aggregated",""))
                    
                    #II.b. Plot the raw radar slice with the surface highlighted
                    #II.b.1. Plot the raw radar echogram
                    cb_1=ax1.pcolor(radar_echo,cmap=pyplot.get_cmap('gray'))#,norm=divnorm)
                    
                    #II.b.2. Display the picked surface over it
                    #pdb.set_trace()
                    ax1.plot(np.arange(0,radar_echo.shape[1]),surface_indices,linestyle='--',color='red',linewidth=0.1)
                    #II.b.3. Display the bootom picked surface over it
                    ax1.plot(np.arange(0,radar_echo.shape[1]),np.asarray(bottom_indices).flatten(),linestyle='--',color='magenta',linewidth=0.1)
                    #II.b.4 Set plot properties
                    ax1.invert_yaxis() #Invert the y axis = avoid using flipud.
                    #ax1.set_aspect('equal') # X scale matches Y scale
                    ax1.set_yticks(ticks_yplot) 
                    ax1.set_yticklabels(np.round(depths[ticks_yplot]))
                    ax1.set_title('Raw radar echogram',fontsize=5)
                    ax1.set_ylabel('Depth [m]')
                    ax1.set_xlabel('Horizontal distance')                    
                    cbar1=fig.colorbar(cb_1, ax=[ax1], location='right')
                    cbar1.set_label('Signal strength', fontsize=5)
                    #pdb.set_trace()
                    #II.c. Plot the radar slice (first 30m of radar echogram)
                    
                    #Create the y vector for plotting
                    ticks_yplot=np.arange(0,radar_slice.shape[0],20)
                    
                    #Plot the radar slice
                    cb=ax2.pcolor(radar_slice,cmap=pyplot.get_cmap('gray'))#,norm=divnorm)
                    ax2.invert_yaxis() #Invert the y axis = avoid using flipud.
                    ax2.set_aspect('equal') # X scale matches Y scale
                    #In order to display the depth, I used the example 1 of the
                    #following site: https://www.geeksforgeeks.org/matplotlib-axes-axes-set_yticklabels-in-python/
                    ax2.set_yticks(ticks_yplot) 
                    ax2.set_yticklabels(np.round(depths[ticks_yplot]))
                    ax2.set_title('Radar echogram slice',fontsize=5)
                    ax2.set_ylabel('Depth [m]')
                    ax2.set_xlabel('Horizontal distance')
                    cbar=fig.colorbar(cb)
                    cbar.set_label('Signal strength', fontsize=5)

                    #Create the figure name
                    fig_name=[]
                    fig_name='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/2002_2003_raw_and_slice/'+indiv_file+'.png'
                    
                    #Save the figure
                    pyplot.savefig(fig_name,dpi=500)
                    pyplot.clf()
                    #Plot the data

                #If plot_original_slice_and_cutted_slice is set to 'TRUE', then
                #plot the radar slice of that date and as well as the radar slice
                #that is free from the very bright surface, save it
                if (plot_original_slice_and_cutted_slice=='TRUE'):
                    #If file have already been created, continue
                    filename_to_check='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/2002_2003_slice_and_cutted_slice/'+indiv_file+'.png'
                    if (os.path.isfile(filename_to_check)):
                        print('Figure already existent, move on to the next date')
                        continue
                                        
                    #Subplot N°1:
                    #I. Process and plot radar echogram
                    #I.a. Load the surface suggestion pick (there is no 'Surface'
                    # variable in 2002/2003 dataset such as 2010/2014 datset).
                    
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
                    
                    #I.c. Select the original radar slice
                    #Define the uppermost and lowermost limits
                    meters_cutoff_above=0
                    meters_cutoff_below=30
                    
                    #Get our slice (30 meters as currently set)
                    original_radar_slice, bottom_indices = _return_radar_slice_given_surface(radar_echo,
                                                                    depths,
                                                                    surface_indices,
                                                                    meters_cutoff_above=meters_cutoff_above,
                                                                    meters_cutoff_below=meters_cutoff_below)
                    # I have taken and adatped the functions '_return_radar_slice_given_surface' and
                    # '_radar_slice_indices_above_and_below' from 'IceBridgeGPR_Manager_v2.py'
                    # and it seems to correctly selecting the slice! I did not manually check
                    # by looking in the variables it the job was done correctly but I have
                    # checked several variables such as idx_above, idx_below, output traces
                    # and it seems okay to me!
                    
                    #Rescale the radar slice as MacFerrin et al. 2019
                    #1.If required to go through _export_to_8bit_array as a vector
                    #radar_slice_rescaled=_export_to_8bit_array(np.ndarray.flatten(radar_slice))
                    #radar_slice_rescaled_mat = radar_slice_rescaled.reshape(radar_slice.shape[0],radar_slice.shape[1])
                    
                    #2.If not required to go through _export_to_8bit_array as a vector
                    original_radar_slice_rescaled_mat=_export_to_8bit_array(original_radar_slice)

                    #I.d. Select the cutted radar slice
                    #Define the uppermost and lowermost limits
                    meters_cutoff_above=4
                    meters_cutoff_below=30
                    
                    #Get our slice (30 meters as currently set)
                    cutted_radar_slice, bottom_indices = _return_radar_slice_given_surface_below_surf(radar_echo,
                                                                    depths,
                                                                    surface_indices,
                                                                    meters_cutoff_above=meters_cutoff_above,
                                                                    meters_cutoff_below=meters_cutoff_below)
                    # I have taken and adatped the functions '_return_radar_slice_given_surface' and
                    # '_radar_slice_indices_above_and_below' from 'IceBridgeGPR_Manager_v2.py'
                    # and it seems to correctly selecting the slice! I did not manually check
                    # by looking in the variables it the job was done correctly but I have
                    # checked several variables such as idx_above, idx_below, output traces
                    # and it seems okay to me!
                    
                    #Rescale the radar slice as MacFerrin et al. 2019
                    #1.If required to go through _export_to_8bit_array as a vector
                    #radar_slice_rescaled=_export_to_8bit_array(np.ndarray.flatten(radar_slice))
                    #radar_slice_rescaled_mat = radar_slice_rescaled.reshape(radar_slice.shape[0],radar_slice.shape[1])
                    
                    #2.If not required to go through _export_to_8bit_array as a vector
                    cutted_radar_slice_rescaled_mat=_export_to_8bit_array(cutted_radar_slice)
    
                    #II. Plot radar echograms
                    #II.a Create the subplot
                    pyplot.figure(figsize=(48,40))
                    #Change label font
                    pyplot.rcParams.update({'font.size': 5})
                    #fig, (ax1, ax2) = pyplot.subplots(1, 2)
                    fig, (ax1, ax2) = pyplot.subplots(2, 1)#, gridspec_kw={'width_ratios': [1, 3]})

                    fig.suptitle(indiv_file.replace("_aggregated",""))
                    
                    #Subplot N°1:
                    #II.b. Plot the original radar slice (first 30m of radar echogram)
                    #Create the y vector for plotting
                    ticks_yplot=np.arange(0,original_radar_slice.shape[0],20)
                    
                    #Plot the radar slice
                    cb=ax1.pcolor(original_radar_slice_rescaled_mat,cmap=pyplot.get_cmap('gray'))#,norm=divnorm)
                    ax1.invert_yaxis() #Invert the y axis = avoid using flipud.
                    ax1.set_aspect('equal') # X scale matches Y scale
                    #In order to display the depth, I used the example 1 of the
                    #following site: https://www.geeksforgeeks.org/matplotlib-axes-axes-set_yticklabels-in-python/
                    ax1.set_yticks(ticks_yplot) 
                    ax1.set_yticklabels(np.round(depths[ticks_yplot]))
                    ax1.set_title('Radar echogram slice, rescaled from 0 to 256, from the surface',fontsize=5)
                    ax1.set_ylabel('Depth [m]')
                    ax1.set_xlabel('Horizontal distance')
                    
                    #Subplot N°2:
                    #II.c. Plot the cutted radar slice (from number defined by 
                    #meters_cutoff_above to 30m deep)
                    #Create the y vector for plotting
                    floored_depths=np.floor(depths)
                    depths_bool=(floored_depths==meters_cutoff_above)
                    depths_trues=[i for i, x in enumerate(depths_bool) if x]
                    ticks_yplot=np.arange(depths_trues[0],cutted_radar_slice_rescaled_mat.shape[0],20)
                    
                    #Plot the radar slice
                    cb=ax2.pcolor(cutted_radar_slice_rescaled_mat,cmap=pyplot.get_cmap('gray'))#,norm=divnorm)
                    ax2.invert_yaxis() #Invert the y axis = avoid using flipud.
                    ax2.set_aspect('equal') # X scale matches Y scale
                    #In order to display the depth, I used the example 1 of the
                    #following site: https://www.geeksforgeeks.org/matplotlib-axes-axes-set_yticklabels-in-python/
                    ax2.set_yticks(ticks_yplot) 
                    ax2.set_yticklabels(np.round(depths[ticks_yplot]))
                    ax2.set_title('Radar echogram slice, rescaled from 0 to 256, from 4m below the surface',fontsize=5)
                    ax2.set_ylabel('Depth [m]')
                    ax2.set_xlabel('Horizontal distance')
                    #cbar=fig.colorbar(cb)
                    #cbar.set_label('Signal strength', fontsize=5)
                    #fig.tight_layout()
                    #pyplot.show()
                    
                    #Create the figure name
                    fig_name=[]
                    fig_name='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/2002_2003_slice_and_cutted_slice/'+indiv_file+'.png'
                    
                    #Save the figure
                    pyplot.savefig(fig_name,dpi=500)
                    pyplot.clf()
                    #Plot the data
                    #pdb.set_trace()
                    
                    continue

                #If plot_slice_and_loc_rescaled is set to 'TRUE', then plot the location of
                #radar echogram AND the rescaled radar slice of that date and save it
                if (plot_slice_and_loc_rescaled=='TRUE'):
                    #If file have already been created, continue
                    filename_to_check='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/2002_2003_slice_and_loc_rescaled/'+indiv_file.replace("_aggregated","")+technique+'.png'
                    if (os.path.isfile(filename_to_check)):
                        print('Figure already existent, move on to the next date')
                        continue
                    
                    #If not a potential ice slabs, continue
                    if (folder_day=='jun04'):
                        if (not(indiv_file.replace(".mat","") in list(potential_iceslabs))):
                            print(indiv_file.replace(".mat",""),'not a potential iceslabs: continue')
                            continue
                    else:
                        if (not(indiv_file.replace("_aggregated","") in list(potential_iceslabs))):
                            print(indiv_file.replace("_aggregated",""),'not a potential iceslabs: continue')
                            continue
                                        
                    #Subplot N°1:
                    #I. Process and plot radar echogram
                    #I.a. Load the surface suggestion pick (there is no 'Surface'
                    # variable in 2002/2003 dataset such as 2010/2014 datset).
                    
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
                        
                    #I.c. Perform depth correction
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
                    # I have taken and adatped the functions '_return_radar_slice_given_surface' and
                    # '_radar_slice_indices_above_and_below' from 'IceBridgeGPR_Manager_v2.py'
                    # and it seems to correctly selecting the slice! I did not manually check
                    # by looking in the variables it the job was done correctly but I have
                    # checked several variables such as idx_above, idx_below, output traces
                    # and it seems okay to me!
                    
                    #DO NOT rescale the radar slice as MacFerrin et al. 2019
                    #radar_slice_rescaled_mat=_export_to_8bit_array(radar_slice)
    
                    ##############################################################
                    ############### Begin explanations on pcolor #################
                    
                    #Explainations on how pcolor works. I convinced myself doing a small example that I plotted.
                    #Further explanations: https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.pcolor.html#:~:text=pcolor()%20displays%20all%20columns,Similarly%20for%20the%20rows.
                    #Example: I have a matrix Z that I want to plot
                    
                    #Z=[10,3,24,70,                            40,5,48,22
                    #    2,6,87,21,      ----pcolor(Z)--->     2,6,87,21
                    #    40,5,48,22]                           10,3,24,70
                    
                    #I must use np.flipud(), or invert the y axis to display the data from top to bottom.
                    
                    ################ End explanations on pcolor ##################
                    ##############################################################
                    
                    #II. Plot radar echogram localisation
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
                        ax1.scatter(lon_3413[0,0],lat_3413[0,0],c='m',s=0.1) #Plot the start in green
                    else:
                        ax1.scatter(lon_3413[0],lat_3413[0],c='m',s=0.1)
                    
                    ax1.grid()
                    
                    if (folder_day=='jun04'):
                        ax1.set_xlim(lon_3413[0,0]-500000, lon_3413[0,0]+500000)
                        ax1.set_ylim(lat_3413[0,0]-500000, lat_3413[0,0]+500000)
                    else:
                        ax1.set_xlim(lon_3413[0]-500000, lon_3413[0]+500000)
                        ax1.set_ylim(lat_3413[0]-500000, lat_3413[0]+500000)
                    
                    #II.b. Plot the radar slice (first 30m of radar echogram)
                    
                    #Subplot N°2:
                    #Create the y vector for plotting
                    ticks_yplot=np.arange(0,radar_slice.shape[0],20)
                    
                    #Select the lower and upper end for rescaling.
                    #These data are from 'yearly_averaged_signal.py' function
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
                    
                    #Plot the radar slice
                    cb2=ax2.pcolor(radar_slice,cmap=pyplot.get_cmap('gray'))#,norm=divnorm)
                    ax2.invert_yaxis() #Invert the y axis = avoid using flipud.
                    ax2.set_aspect('equal') # X scale matches Y scale
                    #In order to display the depth, I used the example 1 of the
                    #following site: https://www.geeksforgeeks.org/matplotlib-axes-axes-set_yticklabels-in-python/
                    ax2.set_yticks(ticks_yplot) 
                    ax2.set_yticklabels(np.round(depths[ticks_yplot]))
                    ax2.set_title('Radar echogram slice, not depth correct, rescaling: '+technique+' percentiles',fontsize=5)
                    ax2.set_ylabel('Depth [m]')
                    ax2.set_xlabel('Horizontal distance')
                    
                    cb2.set_clim(perc_lower_end,perc_upper_end)
                    cbar2=fig.colorbar(cb2, ax=[ax2], location='right')
                    cbar2.set_label('Signal strength')
                    #fig.tight_layout()
                    #pyplot.show()
                    
                    ##Create the figure name
                    fig_name=[]
                    fig_name='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/2002_2003_slice_and_loc_rescaled/'+indiv_file.replace("_aggregated","")+technique+'.png'
                    
                    ##Save the figure
                    pyplot.savefig(fig_name,dpi=500)
                    pyplot.clf()
                    #Plot the data
                    #pdb.set_trace()
                    
                    continue                    
    else:
        print('Folder',folder_year,', continue ...')
        continue

#Open files




#Visualise the data

            #Store everything into one dictionnary (matrix and vectors of data)
            #if (folder_day == 'jun04'):
            #    dic_file_being_read=[]
            #    dic_file_being_read = { "trace_id" : indiv_file.replace(".mat",""),
            #         "radar_echogram" : file_being_read['data'],
            #         "latlontime" : df_final}
            #else:
            #    dic_file_being_read=[]
            #    dic_file_being_read = { "trace_id" : indiv_file.replace(".mat",""),
            #         "radar_echogram" : file_being_read['filtfin'],
            #         "latlontime" : df_final}