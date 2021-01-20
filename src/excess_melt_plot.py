# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 14:29:11 2021

@author: JullienN

Read and plot the excess melt calculation made by Andrew Tedstone
"""

#Import librairies

#Import packages

import scipy.io
import matplotlib.pyplot as plt
import rasterio as rio
import geopandas as gpd
from os import listdir
import os
from os.path import isfile, join
import numpy as np
import xarray as xr
import pickle
import pdb
import pandas as pd
import pyproj

##############################################################################
############################## Define variables ##############################
##############################################################################
dt = 2.034489716724874e-09 #Timestep for 2002/2003 traces
t0 = 0; # Unknown so set to zero

#Compute the speed (Modified Robin speed):
# self.C / (1.0 + (coefficient*density_kg_m3/1000.0))
v= 299792458 / (1.0 + (0.734*0.873/1000.0))

#Define path of the working directory
global_path='C:/Users/jullienn/Documents/working_environment/'

#Define the desired year to plot
desired_year='2001'
generate_raw_excess_melt='FALSE'
generate_excess_melt_traces='FALSE'
generate_excess_melt_traces_with_slices='TRUE'
##############################################################################
############################## Define variables ##############################
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

    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)
##############################################################################
############### Define function for discrete colorbar display ###############
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

##########################################################################
###                      Load Greenland DEM and contours               ###
##########################################################################

dem_filename = 'C:/Users/jullienn/Documents/working_environment/greenland_topo_data/elevations/greenland_dem_mosaic_100m_v3.0.tif'
contours_filename = 'C:/Users/jullienn/Documents/working_environment/greenland_topo_data/elevations/greenland_contours_100m_v3.0.shp'

with rio.open(dem_filename) as fh:
    dem = fh.read()
    print(fh)
    dem_bounds = fh.bounds
    dem_crs = fh.crs
contours = gpd.read_file(contours_filename)
# Make sure that the coordinate reference system matches the image
contours = contours.to_crs(dem_crs)
         
##########################################################################
###                      Load Greenland DEM and contours               ###
##########################################################################


##########################################################################
###                          Load excess melt data   	               ###
##########################################################################
#Define path
path='C:/Users/jullienn/Documents/working_environment/excess_melt/'

#Load the data
data_path = path+'excess_melt_mbyear.nc'
DS = xr.open_dataset(data_path)

#Extract coordinates
lat_M_e=DS.x.data
lon_M_e=DS.y.data

#Load melt data
melt_data= DS.M_e
##########################################################################
###                          Load excess melt data   	               ###
##########################################################################

if (generate_excess_melt_traces=='TRUE'):
    #Create the figures of old dataset traces on top of annual excess melt figures
    
    #Define the working environment for loading the data
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
                    
                    #If file have already been created, continue
                    filename_to_check='C:/Users/jullienn/Documents/working_environment/excess_melt/figures_excess_melt/'+desired_year+'/year_'+desired_year+'_'+indiv_file.replace("_aggregated","")+'.png'
                    if (os.path.isfile(filename_to_check)):
                        print('Figure already existent, move on to the next date')
                        continue
    
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
                        
                        #Transform lat and lon variables into series to be able
                        #to create the pandas dataframes:
                        lat=pd.Series(lat.flatten())
                        lon=pd.Series(lon.flatten())
    
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
                    
                    ######################################################
                    ###     Convert lat/lon into MAR's projection      ###
                    #This piece of code is adapted from Andrew
                    
                    #1. Create a dataframe of traces with lat/lon in EPSG:4326
                    df_latlon=pd.DataFrame({'lon':pd.Series(lon),
                                            'lat':pd.Series(lat)})
                    
                    #2. Create the geodataframe
                    df_latlon_gdf = gpd.GeoDataFrame(df_latlon, geometry=gpd.points_from_xy(df_latlon.lon, df_latlon.lat), crs=4326)
                    
                    #3. Convert the data from EPSG:4326 to MAR's projection
                    mar_proj4 = '+proj=stere +lat_0=70.5 +lon_0=-40 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'
                    df_latlon_gdf = df_latlon_gdf.to_crs(mar_proj4)
                    
                    ###     Convert lat/lon into MAR's projection      ###
                    ######################################################    
                    plt.rcParams.update({'font.size': 40})
                    plt.figure(figsize=(48,40))
                    
                    ax = plt.subplot()
                    plot_melt=DS.M_e.sel(time=desired_year).plot(ax=ax)
                    
                    plot_melt.set_clim(0,1000)
                    plot_melt.set_cmap(discrete_cmap(10, 'hot_r'))
                    #plt.colorbar(label='Excess melt [mm w.e./year]')
                    
                    df_latlon_gdf.plot(ax=ax, marker='o', markersize=1, zorder=45, color='blue')
                    
                    #Display begining of tracks
                    begining_traces=df_latlon_gdf.loc[0:10]
                    begining_traces.plot(ax=ax, marker='o', markersize=1, zorder=45, color='magenta')
                    
                    #Display title
                    plt.title(indiv_file.replace("_aggregated","")+' with excess melt year: '+desired_year)
                    
                    ax.set_xlim(df_latlon_gdf['geometry'].iloc[0].x-500000,df_latlon_gdf['geometry'].iloc[0].x+500000)
                    ax.set_ylim(df_latlon_gdf['geometry'].iloc[0].y-300000,df_latlon_gdf['geometry'].iloc[0].y+300000)
                    
                    #Create the figure name
                    fig_name=[]
                    fig_name='C:/Users/jullienn/Documents/working_environment/excess_melt/figures_excess_melt/'+desired_year+'/year_'+desired_year+'_'+indiv_file.replace("_aggregated","")+'.png'
                    
                    #Save the figure
                    plt.savefig(fig_name)
                    plt.clf()

if (generate_excess_melt_traces_with_slices=='TRUE'):
    #1. Deal only with potential_ice_slabs locations
    #2. Plot on the top the excess melt map with traces location
    #3. Plot on the bottom the slice
    
    #Open, read and close the potential ice slabs file
    f = open('C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/2002_2003_slice_and_loc/potential_iceslabs.txt','r')
    potential_iceslabs = [line.strip() for line in f.readlines() if len(line.strip()) > 0]
    f.close()
    
    #Open, read and close the file of suggested surface picks
    f = open('C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/Exclusion_folder/txt/SURFACE_STARTING_PICKS_Suggestions_2002_2003.txt','r')
    lines = [line.strip() for line in f.readlines() if len(line.strip()) > 0]
    f.close()
    
    #Create the dataframe that define the dates who have experienced a surface
    #picking improvement
    df_dates_surf_pick=pd.DataFrame({'dates_surf_pick_impr':pd.Series(['may24_02_23','may24_02_24','may24_02_25',
                                                                   'may30_02_2','may30_02_4','may30_02_5','may30_02_6',
                                                                   'may30_02_7','may30_02_13','may30_02_14','may30_02_15',
                                                                   'may30_02_50','may30_02_51'])})
    
    #Define the working environment for loading the data
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
                    
                    #If indiv_file is the quality file, continue
                    if (indiv_file[0:7]==('quality')):
                        #pdb.set_trace()
                        continue
                    
                    #If the indiv_file is not in the potential ice slabs file,
                    #move on to the next date
                    if (folder_day=='jun04'):
                        if (not(indiv_file.replace(".mat","") in list(potential_iceslabs))):
                            print(indiv_file+' is not a potential ice slab, move on.')
                            continue
                    else:
                        if (not(indiv_file.replace("_aggregated","") in list(potential_iceslabs))):
                            print(indiv_file+' is not a potential ice slab, move on.')
                            continue
                        
                    #Print which file is processed
                    print('Treating file',indiv_file)
                    
                    ##If file have already been created, continue
                    #filename_to_check='C:/Users/jullienn/Documents/working_environment/excess_melt/figures_excess_melt/'+desired_year+'/year_'+desired_year+'_'+indiv_file.replace("_aggregated","")+'.png'
                    #if (os.path.isfile(filename_to_check)):
                    #    print('Figure already existent, move on to the next date')
                    #    continue

                    #Load the data                    
                    if (folder_day=='jun04'):
                        fdata= scipy.io.loadmat(folder_day_name+'/'+indiv_file)
                        #Select radar echogram and corresponding lat/lon
                        radar_echo=fdata['data']
                        lat=fdata['latitude']
                        lon=fdata['longitude']
                        
                        #Transform lat and lon variables into series to be able
                        #to create the pandas dataframes:
                        lat=pd.Series(lat.flatten())
                        lon=pd.Series(lon.flatten())
    
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

                    #II. Plot radar echogram localisation over excess melt map
                    #Some index have lat and lon equal to 0 because of jumps in data aggregation.
                    #Replace these 0 by NaNs
                    
                    if (not(folder_day=='jun04')):
                        lat.replace(0, np.nan, inplace=True)
                        lon.replace(0, np.nan, inplace=True)
                    
                    #Transform the longitudes. The longitudes are ~46 whereas they should be ~-46! So add a '-' in front of lon
                    lon=-lon
                    
                    ##########################################################
                    ###             SUBPLOT 1: EXCESS MELT plot            ###
                    ##########################################################
                    
                    ######################################################
                    ###     Convert lat/lon into MAR's projection      ###
                    #This piece of code is adapted from Andrew
                    
                    #1. Create a dataframe of traces with lat/lon in EPSG:4326
                    df_latlon=pd.DataFrame({'lon':pd.Series(lon),
                                            'lat':pd.Series(lat)})
                    
                    #2. Create the geodataframe
                    df_latlon_gdf = gpd.GeoDataFrame(df_latlon, geometry=gpd.points_from_xy(df_latlon.lon, df_latlon.lat), crs=4326)
                    
                    #3. Convert the data from EPSG:4326 to MAR's projection
                    mar_proj4 = '+proj=stere +lat_0=70.5 +lon_0=-40 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'
                    df_latlon_gdf = df_latlon_gdf.to_crs(mar_proj4)
                    
                    ###     Convert lat/lon into MAR's projection      ###
                    ###################################################### 
                    
                    #Plot excess melt map and traces localisation
                    plt.rcParams.update({'font.size': 10})
                    plt.figure(figsize=(48,40))
                    fig, (ax1, ax2) = plt.subplots(2, 1)
                    
                    plot_melt=DS.M_e.sel(time=desired_year).plot(ax=ax1)
                    
                    plot_melt.set_clim(0,1000)
                    plot_melt.set_cmap(discrete_cmap(10, 'hot_r'))
                    #plt.colorbar(label='Excess melt [mm w.e./year]')
                    
                    df_latlon_gdf.plot(ax=ax1, marker='o', markersize=1, zorder=45, color='blue')
                    
                    #Display begining of tracks
                    begining_traces=df_latlon_gdf.loc[0:10]
                    begining_traces.plot(ax=ax1, marker='o', markersize=1, zorder=45, color='magenta')
                    
                    ax1.set_xlim(df_latlon_gdf['geometry'].iloc[0].x-500000,df_latlon_gdf['geometry'].iloc[0].x+500000)
                    ax1.set_ylim(df_latlon_gdf['geometry'].iloc[0].y-300000,df_latlon_gdf['geometry'].iloc[0].y+300000)
                    
                    #Display title
                    plt.title(indiv_file.replace("_aggregated","")+' with excess melt year: '+desired_year)
                    
                    ##########################################################
                    ###             SUBPLOT 1: EXCESS MELT plot            ###
                    ##########################################################
                    
                    ##########################################################
                    ###             SUBPLOT 2: Radar slice plot            ###
                    ##########################################################
                    
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
                    
                    #Select the first 30m of radar echogram
                    #1. Compute the vertical resolution
                    #a. Time computation according to John Paden's email.
                    Nt = radar_echo.shape[0]
                    Time = t0 + dt*np.arange(1,Nt+1)
                    #b. Calculate the depth:
                    #self.SAMPLE_DEPTHS = self.radar_speed_m_s * self.SAMPLE_TIMES / 2.0
                    depths = v * Time / 2.0
                    
                    #I.c. Perform depth correction
                    #depth_corrected_traces=perform_depth_correction(radar_echo, depths, surface_indices, indiv_file.replace("_aggregated",""), 'FALSE')
                    
                    #I.d. Select the radar slice
                    #Define the uppermost and lowermost limits
                    meters_cutoff_above=0
                    meters_cutoff_below=30
                    
                    #Redefine the 'surface_indices' variable: the surface have just been picked
                    #for the depth correction, so now we want to pick from the top down to
                    #30m depth, thus the 'surface_indices' must be [0,0,...,0]!!
                    #surface_indices=np.zeros(surface_indices.shape[0],dtype=np.int64)
                    
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
                    
                    #Rescale the radar slice as MacFerrin et al. 2019
                    #If not required to go through _export_to_8bit_array as a vector
                    radar_slice_rescaled_mat=_export_to_8bit_array(radar_slice)
                    
                    #II.b. Plot the radar slice (first 30m of radar echogram)
                    
                    #Subplot N°2:
                    #Create the y vector for plotting
                    ticks_yplot=np.arange(0,radar_slice.shape[0],20)
                    
                    #Plot the radar slice
                    cb=ax2.pcolor(radar_slice_rescaled_mat,cmap=plt.get_cmap('gray'))#,norm=divnorm)
                    ax2.invert_yaxis() #Invert the y axis = avoid using flipud.
                    ax2.set_aspect('equal') # X scale matches Y scale
                    #In order to display the depth, I used the example 1 of the
                    #following site: https://www.geeksforgeeks.org/matplotlib-axes-axes-set_yticklabels-in-python/
                    ax2.set_yticks(ticks_yplot) 
                    ax2.set_yticklabels(np.round(depths[ticks_yplot]))
                    ax2.set_title('Radar echogram slice, rescaled from 0 to 256 - 5-95% - No depth correction',fontsize=5)
                    ax2.set_ylabel('Depth [m]')
                    ax2.set_xlabel('Horizontal distance')
                    #cbar=fig.colorbar(cb)
                    #cbar.set_label('Signal strength', fontsize=5)
                    #fig.tight_layout()
                    plt.show()        
                    
                    ##########################################################
                    ###             SUBPLOT 2: Radar slice plot            ###
                    ##########################################################                     
                    
                    ##Create the figure name
                    #fig_name=[]
                    #fig_name='C:/Users/jullienn/Documents/working_environment/excess_melt/figures_excess_melt/'+desired_year+'/year_'+desired_year+'_'+indiv_file.replace("_aggregated","")+'.png'
                    #
                    ##Save the figure
                    #plt.savefig(fig_name)
                    #plt.clf()
                    
                    pdb.set_trace()

if (generate_raw_excess_melt=='TRUE'):
    #Generate and save the raw annual excess melt figures
    for year in list(np.arange(1990,2020)):
        
        #Define the year
        wanted_year=str(year)
        
        #Select the data associated with the wanted year
        melt_year = melt_data.sel(time=wanted_year)
        melt_year_np = melt_year.values
        
        melt_year_plot=np.asarray(melt_year_np)
        melt_year_plot=melt_year_plot[0,0,:,:,0]
        
        #Plot dem and contours elevation
        plt.rcParams.update({'font.size': 20})
        plt.figure(figsize=(48,40))
        ax = plt.subplot(111)
        dem_extent = (dem_bounds[0], dem_bounds[2], dem_bounds[1], dem_bounds[3])
        plt.imshow(np.squeeze(np.flipud(melt_year_plot)), extent=dem_extent,cmap=discrete_cmap(5,'hot_r'))
        
        plt.colorbar(label='Excess melt [mm w.e./year]')
        plt.clim(0,1000)
        ax.grid()
        #contours.plot(ax=ax, edgecolor='black')
        plt.title('Excess melt plot, year: '+wanted_year)
    
        #Create the figure name
        fig_name=[]
        fig_name='C:/Users/jullienn/Documents/working_environment/excess_melt/figures_excess_melt/excess_melt_'+wanted_year+'.png'
                            
        #Save the figure
        plt.savefig(fig_name)
        plt.clf()








