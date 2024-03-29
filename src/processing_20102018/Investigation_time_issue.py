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


plot_radar_echogram_slice='TRUE'

dt = 2.034489716724874e-09 #Timestep for 2002/2003 traces
t0 = 0; # Unknown so set to zero

#Compute the speed (Modified Robin speed):
# self.C / (1.0 + (coefficient*density_kg_m3/1000.0))
v= 299792458 / (1.0 + (0.734*0.873/1000.0))


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
    #pdb.set_trace()
    
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
##############################################################################
################## Define functions for radar slice picking ##################
##############################################################################

##############################################################################
############# Define function for depth correction of the traces #############
##############################################################################
#Function taken from IceBridgeGPR_Manager_v2.py

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

##############################################################################
################### Define function for radargram display ####################
##############################################################################

##############################################################################
################### Define function for ice lenses logging ###################
##############################################################################
#This function if adapted from https://stackoverflow.com/questions/37363755/python-mouse-click-coordinates-as-simply-as-possible
def onclick(event):
    #This functions print and save the x and y coordinates in pixels!
    print(event.xdata, event.ydata)
    #Fill in the file to log on the information
    filename_flog='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/2002_2003_radar_slice/flog_icelenses_alldates.txt'
    f_log = open(filename_flog, "a")
    f_log.write(str(round(event.xdata,2))+','+str(round(event.ydata,2))+'\n')
    f_log.close() #Close the quality assessment file when we’re done!
##############################################################################
################### Define function for ice lenses logging ###################
##############################################################################

#Open, read and close the file of suggested surface picks
f = open('C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/Exclusion_folder/txt/SURFACE_STARTING_PICKS_Suggestions.txt','r')
lines = [line.strip() for line in f.readlines() if len(line.strip()) > 0]
f.close()

#Define the working environment
path= 'C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data'
os.chdir(path) # relative path: scripts dir is under Lab

# Read the years of data
folder_years = [ f.name for f in os.scandir(path) if f.is_dir() ]

for folder_year in folder_years:
    if (folder_year in list(['2002','2003'])):
        continue
    
    if (folder_year == '2014_Greenland_P3'):
        print('Treating the year',folder_year)

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
                ##Load data
                #fdata_p0= scipy.io.loadmat(folder_day_name+'/'+indiv_file)
                #fdata_p1= scipy.io.loadmat(folder_day_name+'/Data_20100507_01_009.mat')
                #fdata_p2= scipy.io.loadmat(folder_day_name+'/Data_20100507_01_010.mat')
                #
                #pdb.set_trace()
                #
                ##Select radar echogram
                #radar_echo_p0=fdata_p0['Data']
                #radar_echo_p1=fdata_p1['Data']
                #radar_echo_p2=fdata_p2['Data']
                #
                #radar_echo_p0_p1=np.append(radar_echo_p0,radar_echo_p1,axis=1)
                #radar_echo=np.append(radar_echo_p0_p1,radar_echo_p2,axis=1)
                #
                ##Tranform in log10
                #radar_echo=np.log10(radar_echo)

                if (folder_year=='2010_Greenland_P3'):
                    
                    plot_name1='20100508_01_114'
                    plot_name2='20100508_01_114'
                    plot_name3='20100508_01_115'
                    
                if (folder_year=='2011_Greenland_P3'):
                    
                    plot_name1='20110419_01_008'
                    plot_name2='20110419_01_009'
                    plot_name3='20110419_01_010'
                    
                if (folder_year=='2012_Greenland_P3'):
                    
                    plot_name1='20120418_01_129'
                    plot_name2='20120418_01_130'
                    plot_name3='20120418_01_131'
                                        
                if (folder_year=='2013_Greenland_P3'):
                    
                    plot_name1='20130405_01_165'
                    plot_name2='20130405_01_166'
                    plot_name3='20130405_01_167'
                    
                #if (folder_year=='2014_Greenland_P3'):
                #    
                #    plot_name1='20140424_01_002'
                #    plot_name2='20140424_01_003'
                #    plot_name3='20140424_01_004'
                #pdb.set_trace()
                
                #Load data
                fdata1= scipy.io.loadmat(folder_day_name+'/Data_'+plot_name1+'.mat')
                radar_echo1=fdata1['Data']
                time_echo1=fdata1['Time']
                    
                fdata2= scipy.io.loadmat(folder_day_name+'/Data_'+plot_name2+'.mat')
                radar_echo2=fdata2['Data']
                time_echo2=fdata2['Time']
                    
                fdata3= scipy.io.loadmat(folder_day_name+'/Data_'+plot_name3+'.mat')
                radar_echo3=fdata3['Data']
                time_echo3=fdata3['Time']
                
                #Plot the data
                
                #Create the subplot
                pyplot.figure(figsize=(48,40))
                pyplot.rcParams.update({'font.size': 5})
                fig, (ax1, ax2) = pyplot.subplots(1, 2)#, gridspec_kw={'width_ratios': [1, 3]})
    
                fig.suptitle(str(plot_name1))
    
                #Plot the radar slice
                cb1=ax1.pcolor(np.log10(radar_echo1),cmap=pyplot.get_cmap('gray'))#,norm=divnorm)
                ax1.invert_yaxis() #Invert the y axis = avoid using flipud.
                ax1.set_aspect('equal') # X scale matches Y scale
                ax1.set_title('log10(radar echo)')
                ax1.set_ylabel('Depth [m]')
                ax1.set_xlabel('Horizontal distance')
                cbar1=fig.colorbar(cb1, ax=[ax1], location='left')
                cbar1.set_label('Signal strength')
    
                ax2.plot(time_echo1)
                ax2.grid()
                ax2.set_title('Time')
                ax2.set_xlabel('1:length(time)')
                ax2.set_ylabel('Time [s]')
                
                fig_name=[]
                fig_name='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/investigation_2012_2013_bug/'+folder_year[0:4]+'_example.png'
                #Save the figure
                pyplot.savefig(fig_name,dpi=500)
                pyplot.clf()
                    
                #Create the subplot
                pyplot.figure(figsize=(48,40))
                pyplot.rcParams.update({'font.size': 5})
                fig, (ax1, ax2, ax3) = pyplot.subplots(1, 3)#, gridspec_kw={'width_ratios': [1, 3]})
    
                fig.suptitle('Time variable')
                    
                #Subplot N°1:
                ax1.plot(time_echo1)
                ax1.grid()
                ax1.set_title(str(plot_name1))
                ax1.set_xlabel('1:length(time)')
                ax1.set_ylabel('Time [s]')
                    
                #Subplot N°2:
                ax2.plot(time_echo2)
                ax2.grid()
                ax2.set_title(str(plot_name2))
                ax2.set_xlabel('1:length(time)')
                ax2.set_ylabel('Time [s]')
    
                #Subplot N°1:
                ax3.plot(time_echo3)
                ax3.grid()
                ax3.set_title(str(plot_name3))
                ax3.set_xlabel('1:length(time)')
                ax3.set_ylabel('Time [s]')
                
                fig_name=[]
                fig_name='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/investigation_2012_2013_bug/'+folder_year[0:4]+'_time_example.png'
                #Save the figure
                pyplot.savefig(fig_name,dpi=500)
                pyplot.clf()
                    
                pdb.set_trace()
                break
                
                #Select the first 30m of radar echogram
                #1. Compute the vertical resolution
                #a. Time computation according to John Paden's email.
                Nt = radar_echo.shape[0]
                Time = t0 + dt*np.arange(1,Nt+1)
                #b. Calculate the depth:
                #self.SAMPLE_DEPTHS = self.radar_speed_m_s * self.SAMPLE_TIMES / 2.0
                depths = v * Time / 2.0
                
                #If plot_radar_echogram_slice is set to 'TRUE', then plot the slice
                #radar echogram of that date and save it
                if (plot_radar_echogram_slice=='TRUE'):
                    
                    if (indiv_file == 'Data_20100507_01_008.mat'):
                        suggested_pixel= 1850
                    
                    surface_indices=kernel_function(radar_echo, suggested_pixel)
                    
                    #I.d. Select the radar slice
                    #Define the uppermost and lowermost limits
                    meters_cutoff_above=0
                    meters_cutoff_below=30

                    #Get our slice (30 meters as currently set)
                    radar_slice, bottom_indices = _return_radar_slice_given_surface(radar_echo,
                                                                    depths,
                                                                    surface_indices,
                                                                    meters_cutoff_above=meters_cutoff_above,
                                                                    meters_cutoff_below=meters_cutoff_below)
                    

                    #2.If not required to go through _export_to_8bit_array as a vector
                    #radar_slice=_export_to_8bit_array(radar_slice)

                    #Generate the pick for vertical distance display
                    ticks_yplot=np.arange(0,radar_slice.shape[0],20)
                    
                    #I.d. Plot the radar slice (first 30m of radar echogram)
                    #pdb.set_trace()
                    #Plot the data            
                    
                    fig=pyplot.figure(figsize=(40,10))
                    
                    #Change label font
                    pyplot.rcParams.update({'font.size': 20})
                    
                    color_map=pyplot.pcolor(radar_slice,cmap=pyplot.get_cmap('gray'))#,norm=divnorm)
                    pyplot.yticks(ticks=ticks_yplot,labels=(np.round(depths[ticks_yplot])))
                    pyplot.gca().invert_yaxis() #Imvert the y axis = avoid using flipud.
                    pyplot.gca().set_aspect('equal') # X scale matches Y scale
                    pyplot.ylabel('Depth [m]')
                    pyplot.xlabel('Horizontal distance')
                    #pyplot.clim(lowerb_plot,upperb_plot)
                    #pyplot.yticks(ticks=ticks_yplot,labels=labels_yplot)
                    #pyplot.xticks(fontsize=20)
                    #pyplot.yticks(fontsize=20)
                    #pyplot.ylim(0, 200)
                    pyplot.title('Radar echogram slice: Data_20100507_01_008_010')

                    #cbar=pyplot.colorbar()
                    #cbar.set_label('Signal strength')

                    pyplot.show()
                    pdb.set_trace()
                    
                    continue
                    
    else:
        print('Folder',folder_year,', continue ...')
        continue