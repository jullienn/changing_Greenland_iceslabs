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

surf_pick_selection='TRUE'

#N defines the number of different colors I want to use for the elevation plot
N=10
##############################################################################
############################## Define variables ##############################
##############################################################################

def onclick(event):
    #This functions print and save the x and y coordinates in pixels!
    print(event.xdata, event.ydata)

##############################################################################
############# Define kernel function for surface identification ##############
##############################################################################
#_gaussian function taken from IceBridgeGPR_Manager_v2.py
# Define a quick guassian function to scale the cutoff mask above
def _gaussian(x,mu,sigma):
    return np.exp(-np.power((x-mu)/sigma, 2.)/2.)

#This function have been taken from 'IceBridgeGPR_Manager_v2.py
def modified_kernel_function(traces_input,suggested_pixel,drift_start):
    
    traces = traces_input
    #Do not take the log10 of traces because 'data have been detrented in the log domain' according to John Paden's email, so I guess they are already log10!
    #traces = np.log10(traces)
        
    ##########################################################################
    #           Plot base figure for evolution of surface picking            #
    ##########################################################################
    
    #Create the y vector for plotting
    ticks_yplot=np.arange(0,traces.shape[0],200)
    
    #II. Plot radar echogram localisation
    #II.a. Create the subplot
    #pdb.set_trace()
    fig = pyplot.figure(figsize=(6,5))
    pyplot.rcParams.update({'font.size': 5})
    fig.suptitle('Raw radar echogram: '+indiv_file.replace("_aggregated",""),fontsize=5)
    
    #II.b. Plot the raw radar slice with the surface highlighted
    #II.b.1. Plot the raw radar echogram
    cb_1=pyplot.pcolor(traces,cmap=pyplot.get_cmap('gray'))#,norm=divnorm)
    pyplot.yticks(ticks_yplot, np.round(depths[ticks_yplot]))
    #pyplot.set_title('Raw radar echogram',fontsize=5)
    pyplot.ylabel('Depth [m]')
    pyplot.xlabel('Horizontal distance [ ]')                    
    cbar1=fig.colorbar(cb_1)
    cbar1.set_label('Signal strength [dB]', fontsize=5)
    
    pyplot.axis([0,traces.shape[1],1250,3500])
    
    pyplot.show()
    
    ##########################################################################
    #           Plot base figure for evolution of surface picking            #
    ##########################################################################
             
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
        
        #####################################################################
        #               Plot the evolution of surface picking               #
        #####################################################################
        #II.b.2. Display the picked surface over it
        pyplot.plot(np.arange(0,traces.shape[1]),improved_indices,linestyle='--',color='red',linewidth=0.1)        
        
        fig.canvas.mpl_connect('button_press_event', onclick)
        pyplot.show()
        #pdb.set_trace()
        #####################################################################
        #               Plot the evolution of surface picking               #
        #####################################################################
        #Pick up the pixel where it start to drift
        if (i>=drift_start[0]):
            pdb.set_trace()
        
        #Change the last index of imprived_indices manually if I see a drift:
        #   This line is to execute directly in the console.
        #improved_indices[i]=
        #last_best_index
    
    pdb.set_trace()
    #If there are pixels with particularly strong echo that are being erroneously
    #picked up as the surface, erase most the little "jump" artifacts in
    #the surface picker.
    improved_indices = _get_rid_of_false_surface_jumps(improved_indices)
    
    #I do not use any mask so I think I shouldn't need to use that:
    ###### Must re-expand the surface indices to account for masked values (filled w/ nan)
    ##### improved_indices_expanded = self._refill_array(improved_indices, surface_maskname)
    
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

#Define the working environment
path= 'C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data'
os.chdir(path) # relative path: scripts dir is under Lab

#Create the dataframe that define the dates that need the surface picking
#to be semi-automatic as well as the pixel start to drift
df_dates_surf_pick=pd.DataFrame({'dates_surf_pick_impr':pd.Series(['may24_02_23','may24_02_24','may24_02_25','may30_02_1',
                                                                   'may30_02_2','may30_02_4','may30_02_5','may30_02_6',
                                                                   'may30_02_7','may30_02_13','may30_02_14','may30_02_15',
                                                                   'may30_02_50','may30_02_51']),
                                 'start_drift_vect':pd.Series([410,150,800,8,
                                                               470,550,30,35,
                                                               40,20,150,30,
                                                               180,450])})

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
                
                #If indiv_file is the quality file, continue
                if (indiv_file[0:7]==('quality')):
                    #pdb.set_trace()
                    continue
                
                if (not(indiv_file.replace("_aggregated","") in list(df_dates_surf_pick['dates_surf_pick_impr']))):
                    print('Do need to improve surface picking of '+indiv_file)
                    continue
                
                if (indiv_file == 'may30_02_1_aggregated'):
                    print('Will do this one later on '+indiv_file)
                    continue         
                
                filename_to_check='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/2002_2003_radar_slice/'+indiv_file.replace("_aggregated","")+'_improved_surfpick.png'
                if (os.path.isfile(filename_to_check)):
                    print('Surface have already been improved, continue')
                    continue

                print('Start surface picking improvement procedure of '+indiv_file)
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
                
                #If surf_pick_selection is set to 'TRUE', then plot the raw radar
                # echogram with the surface overlayed AND the radar slice of
                #that date and save it
                if (surf_pick_selection=='TRUE'):
                    ##If file have already been created, continue
                    #filename_to_check='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/2002_2003_raw_and_slice/'+indiv_file+'.png'
                    #if (os.path.isfile(filename_to_check)):
                    #    print('Figure already existent, move on to the next date')
                    #    continue
                    
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
                    
                    #Grab the pixel where it start to drift
                    start_drift=[]
                    start_drift=list(df_dates_surf_pick[df_dates_surf_pick['dates_surf_pick_impr']==indiv_file.replace("_aggregated","")]['start_drift_vect'])
                    pdb.set_trace()
                                        
                    #I.b. Call the modified_kernel_function to pick the surface
                    surface_indices=modified_kernel_function(radar_echo, suggested_pixel,start_drift)
                    
                    #I.c. Save the surface_indices for this specific date
                    
                    #Log the date we are dealing with in the ice lenses location file
                    filename_surf='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/2002_2003_radar_slice/surf_'+indiv_file.replace(".mat","")+'.txt'
                    
                    f_surf = open(filename_surf, "a")
                    for i in range(0,len(surface_indices)):
                        f_surf.write(str(surface_indices[i])+'\n')
                    f_surf.close() #Close the file
                    
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
                    # I have taken and adatped the functions '_return_radar_slice_given_surface' and
                    # '_radar_slice_indices_above_and_below' from 'IceBridgeGPR_Manager_v2.py'
                    # and it seems to correctly selecting the slice! I did not manually check
                    # by looking in the variables it the job was done correctly but I have
                    # checked several variables such as idx_above, idx_below, output traces
                    # and it seems okay to me!
                    pdb.set_trace()
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
                    
                    pyplot.show()
                    
                    ##Create the figure name
                    #fig_name=[]
                    #fig_name='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/2002_2003_raw_and_slice/'+indiv_file+'.png'
                    
                    ##Save the figure
                    #pyplot.savefig(fig_name,dpi=500)
                    #pyplot.clf()
                    
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