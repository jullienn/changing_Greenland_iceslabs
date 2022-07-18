# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 14:25:01 2021

#Description of the code:
    
This code compared radar slices of two dates

@author: JullienN
"""
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


def compute_distances(eastings,northings):
    #This function is from plot_2002_2003.py, which was originally taken from MacFerrin et al., 2019
    '''Compute the distance (in m here, not km as written originally) of the traces in the file.'''
    # C = sqrt(A^2  + B^2)
    distances = np.power(np.power((eastings[1:] - eastings[:-1]),2) + np.power((northings[1:] - northings[:-1]),2), 0.5)

    #Calculate the cumsum of the distances
    cumsum_distances=np.nancumsum(distances)
    #Seeting the first value of the cumsum to be zero as it is the origin
    return_cumsum_distances=np.zeros(eastings.shape[0])
    return_cumsum_distances[1:eastings.shape[0]]=cumsum_distances

    return return_cumsum_distances

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

##############################################################################
############################## Define variables ##############################
##############################################################################

#Open, read and close the file of suggested surface picks
f = open('C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/Exclusion_folder/txt/SURFACE_STARTING_PICKS_Suggestions_2002_2003.txt','r')
lines = [line.strip() for line in f.readlines() if len(line.strip()) > 0]
f.close()

############################ Load GrIS elevation ##############################
#This is from fig3_paper_iceslabs.py
#Open and display satelite image behind map
from pyproj import CRS
import rioxarray as rxr
import cartopy.crs as ccrs
import geopandas as gpd
#This section of displaying sat data was coding using tips from
#https://www.earthdatascience.org/courses/use-data-open-source-python/intro-raster-data-python/raster-data-processing/reproject-raster/
#https://towardsdatascience.com/visualizing-satellite-data-using-matplotlib-and-cartopy-8274acb07b84

path_DEM='C:/Users/jullienn/switchdrive/Private/research/backup_Aglaja/working_environment/greenland_topo_data/elevations/'
#Load DEM for display
GrIS_DEM_display = rxr.open_rasterio(path_DEM+'greenland_dem_mosaic_100m_v3.0.tif',
                              masked=True).squeeze()
#Load DEM for elevation pick up
GrIS_DEM = rasterio.open(path_DEM+'greenland_dem_mosaic_100m_v3.0.tif')

###################### From Tedstone et al., 2022 #####################
#from plot_map_decadal_change.py
# Define the CartoPy CRS object.
crs = ccrs.NorthPolarStereo(central_longitude=-45., true_scale_latitude=70.)
# This can be converted into a `proj4` string/dict compatible with GeoPandas
crs_proj4 = crs.proj4_init
###################### From Tedstone et al., 2022 #####################

#Define extents
#ease_extent = [-500000, 100000, -3000000, -2000000]
#ease_extent = [west limit, east limit., south limit, north limit]
extent_DEM = [np.asarray(GrIS_DEM_display.x[0]), np.asarray(GrIS_DEM_display.x[-1]), np.asarray(GrIS_DEM_display.y[-1]), np.asarray(GrIS_DEM_display.y[0])]

#Load DEM clipped over the SW
GrIS_DEM_display_SW = rxr.open_rasterio(path_DEM+'SW_zoom/greenland_dem_mosaic_100m_v3.0_SW.tif',
                              masked=True).squeeze()
extent_DEM_SW = [np.asarray(GrIS_DEM_display_SW.x[0]), np.asarray(GrIS_DEM_display_SW.x[-1]), np.asarray(GrIS_DEM_display_SW.y[-1]), np.asarray(GrIS_DEM_display_SW.y[0])]

############################ Load GrIS elevation ##############################

### -------------------------- Load shapefiles --------------------------- ###
#Load Rignot et al., 2016 Greenland drainage bassins
path_rignotetal2016_GrIS_drainage_bassins='C:/Users/jullienn/switchdrive/Private/research/backup_Aglaja/working_environment/greenland_topo_data/GRE_Basins_IMBIE2_v1.3/'
GrIS_drainage_bassins=gpd.read_file(path_rignotetal2016_GrIS_drainage_bassins+'GRE_Basins_IMBIE2_v1.3_EPSG_3413.shp',rows=slice(51,57,1)) #the regions are the last rows of the shapefile

#Extract indiv regions and create related indiv shapefiles
NO_rignotetal=GrIS_drainage_bassins[GrIS_drainage_bassins.SUBREGION1=='NO']
NE_rignotetal=GrIS_drainage_bassins[GrIS_drainage_bassins.SUBREGION1=='NE']
SE_rignotetal=GrIS_drainage_bassins[GrIS_drainage_bassins.SUBREGION1=='SE']
SW_rignotetal=GrIS_drainage_bassins[GrIS_drainage_bassins.SUBREGION1=='SW']
CW_rignotetal=GrIS_drainage_bassins[GrIS_drainage_bassins.SUBREGION1=='CW']
NW_rignotetal=GrIS_drainage_bassins[GrIS_drainage_bassins.SUBREGION1=='NW']
### -------------------------- Load shapefiles --------------------------- ###

### ---------------------- Load 2010-18 ice slabs  ------------------------ ###
path_df_with_elevation='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2002_2018/final_excel/' 
#Load 2010-2018 high estimate
f_20102018_high = open(path_df_with_elevation+'high_estimate/df_20102018_with_elevation_high_estimate_rignotetalregions', "rb")
df_2010_2018_high = pickle.load(f_20102018_high)
f_20102018_high.close()
### ---------------------- Load 2010-18 ice slabs  ------------------------ ###

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
                    fig = plt.figure(figsize=(10,4))
                    gs = gridspec.GridSpec(10, 4)
                    gs.update(wspace=0.001)
                    #gs.update(wspace=0.001)
                    #projection set up from https://stackoverflow.com/questions/33942233/how-do-i-change-matplotlibs-subplot-projection-of-an-existing-axis
                    ax1 = plt.subplot(gs[0:10, 0:2],projection=crs)
                    ax2 = plt.subplot(gs[0:10, 2:4])
                    #Change label font
                    plt.rcParams.update({'font.size': 5})
                    
                    ease_extent = [-500000, 100000, -3000000, -2000000]
                    ax1.set_extent(extent_DEM, crs=crs)
                    #Subplot N°1:
                    #II.a.3. Plot dem
                    ax1.imshow(GrIS_DEM_display[:,:], extent=extent_DEM, transform=crs,cmap='gray', origin='upper') #NIR
                    ax1.gridlines(color='gray', linestyle='--')
                    ax1.coastlines()
                    
                    ax1.set_xlim(ease_extent[0],ease_extent[1])
                    ax1.set_xlim(ease_extent[2],ease_extent[3])

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
                    cb2=ax2.pcolor(radar_slice,cmap=plt.get_cmap('gray'))#,norm=divnorm)
                    ax2.invert_yaxis() #Invert the y axis = avoid using flipud.
                    ax2.set_aspect('equal') # X scale matches Y scale
                    #In order to display the depth, I used the example 1 of the
                    #following site: https://www.geeksforgeeks.org/matplotlib-axes-axes-set_yticklabels-in-python/
                    ax2.set_yticks(ticks_yplot) 
                    ax2.set_yticklabels(np.round(depths[ticks_yplot]))
                    ax2.set_ylabel('Depth [m]')
                    ax2.set_xlabel('Horizontal distance')
                        
                    #Colorbar custom
                    cb2.set_clim(perc_lower_end,perc_upper_end)
                    cbar2=fig.colorbar(cb2, ax=[ax2], location='right')
                    cbar2.set_label('Signal strength')
                    plt.show()
                                        
                    #Compute distance
                    distances=compute_distances(lon_3413[0],lat_3413[0])/1000
                    
                    # ---------------------------- Extract elevation ------------------------ #
                    #Create the vector for elevations storing
                    vect_for_elevation=[]
                    #This is from extract_elevation.py
                    for i in range(0,len(lon_3413[0])):
                        #This is from https://gis.stackexchange.com/questions/190423/getting-pixel-values-at-single-point-using-rasterio
                        for val in GrIS_DEM.sample([(lon_3413[0][i], lat_3413[0][i])]):
                            #Calculate the corresponding elevation
                            vect_for_elevation=np.append(vect_for_elevation,val)
                    # ---------------------------- Extract elevation ------------------------ #
                    
                    #II.a.2 Create the subplot
                    fig = plt.figure(figsize=(10,12))
                    gs = gridspec.GridSpec(10, 12)
                    plt.rcParams.update({'font.size': 15})
                    
                    #projection set up from https://stackoverflow.com/questions/33942233/how-do-i-change-matplotlibs-subplot-projection-of-an-existing-axis
                    ax1 = plt.subplot(gs[0:3, 0:11])
                    ax2 = plt.subplot(gs[3:10, 0:2])
                    ax3 = plt.subplot(gs[3:10, 2:4])
                    ax4 = plt.subplot(gs[3:10, 4:6])
                    ax5 = plt.subplot(gs[3:10, 6:8])
                    ax6 = plt.subplot(gs[3:10, 8:10])
                    #projection set up from https://stackoverflow.com/questions/33942233/how-do-i-change-matplotlibs-subplot-projection-of-an-existing-axis
                    #ax7_general = plt.subplot(gs[3:6, 10:12],projection=crs)
                    ax7 = plt.subplot(gs[3:10, 10:12],projection=crs)
                    gs.update(wspace=0.1)
                    gs.update(hspace=4.5)

                    #Plot the radar slice
                    cb1=ax1.pcolor(np.arange(0,radar_slice.shape[1]),depths[0:99],radar_slice,cmap=plt.get_cmap('gray'))
                                        
                    #Display the distance
                    xticks_plot=np.array([0,200,400,600,800,999])                    
                    ax1.set_xticks(xticks_plot)
                    
                    #Display elevation on the top xticklabels
                    #This is from https://stackoverflow.com/questions/19884335/matplotlib-top-bottom-ticks-different "Zaus' reply"
                    ax1_elev = ax1.secondary_xaxis('top')
                    ax1_elev.set_xticks(xticks_plot)
                    
                    #Display xticks labels
                    ax1.set_xticklabels(np.round(distances[xticks_plot]).astype(int))
                    ax1_elev.set_xticklabels((np.round(vect_for_elevation[xticks_plot]).astype(int)))

                    #Insert colorbar in plot from https://matplotlib.org/stable/gallery/axes_grid1/demo_colorbar_with_inset_locator.html
                    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
                    axins1 = inset_axes(ax1,
                    width="2%",  # width = 50% of parent_bbox width
                    height="100%",  # height : 5%
                    loc='center right',
                    borderpad=-3) #For borderpad, fropm https://stackoverflow.com/questions/13310594/positioning-the-colorbar

                    ax1.invert_yaxis() #Invert the y axis = avoid using flipud.
                    #ax1.set_aspect('equal') # X scale matches Y scale
                    ax1.set_ylabel('Depth [m]')
                    ax1.set_xlabel('Distance [km]')
                    ax1_elev.set_xlabel('Elevation [m]')
                    ax1.text(0.01, 0.78,'a',ha='center', va='center', transform=ax1.transAxes,fontsize=15,zorder=10,weight='bold',color='black')#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot

                    #Colorbar custom
                    cb1.set_clim(perc_lower_end,perc_upper_end)
                    cbar1=fig.colorbar(cb1, cax=axins1, orientation = 'vertical')
                    cbar1.set_label('Signal strength [dB]')
                    
                    ##Create figure or radar signal evolution
                    #start_vertical_cut=370
                    #end_vertical_cut=420
                    
                    color_plot='#1a9850'
                    start_vertical_cut=277
                    end_vertical_cut=start_vertical_cut+10
                    ax1.axvline(start_vertical_cut,color=color_plot)
                    ax1.axvline(end_vertical_cut,color=color_plot)
                    #If several vertical lines to plot and plot average
                    for i in range (start_vertical_cut,end_vertical_cut+1):
                        ax2.plot(radar_slice[:,i],np.arange(0,radar_slice[:,i].size),color='0.7')
                    ax2.plot(np.mean(radar_slice[:,start_vertical_cut:end_vertical_cut],1),np.arange(0,radar_slice[:,i].size),color=color_plot)
                    ax2.invert_yaxis() #Invert the y axis
                    ax2.grid()
                    ax2.set_xlim(-6,6)
                    ax2.set_ylim(radar_slice.shape[0],-5)
                    ax2.text(0.04, 0.97,'b',ha='center', va='center', transform=ax2.transAxes,fontsize=15,zorder=10,weight='bold',color='black')#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot
                    #Display sectors of interest
                    ax7.scatter(lon_3413[0][int((start_vertical_cut+end_vertical_cut)/2)],lat_3413[0][int((start_vertical_cut+end_vertical_cut)/2)],s=200,marker='.',linewidths=0,color=color_plot,zorder=10)

                    color_plot='#2166ac'
                    start_vertical_cut=384
                    end_vertical_cut=start_vertical_cut+10
                    ax1.axvline(start_vertical_cut,color=color_plot)
                    ax1.axvline(end_vertical_cut,color=color_plot)
                    #If several vertical lines to plot and plot average
                    for i in range (start_vertical_cut,end_vertical_cut+1):
                        ax3.plot(radar_slice[:,i],np.arange(0,radar_slice[:,i].size),color='0.7')
                    ax3.plot(np.mean(radar_slice[:,start_vertical_cut:end_vertical_cut],1),np.arange(0,radar_slice[:,i].size),color=color_plot)
                    ax3.invert_yaxis() #Invert the y axis
                    ax3.grid()
                    ax3.set_xlim(-6,6)
                    ax3.set_ylim(radar_slice.shape[0],-5)
                    ax3.set_yticklabels([])
                    ax3.text(0.04, 0.97,'c',ha='center', va='center', transform=ax3.transAxes,fontsize=15,zorder=10,weight='bold',color='black')#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot
                    #Display sectors of interest
                    ax7.scatter(lon_3413[0][int((start_vertical_cut+end_vertical_cut)/2)],lat_3413[0][int((start_vertical_cut+end_vertical_cut)/2)],s=200,marker='.',linewidths=0,color=color_plot,zorder=10)

                    color_plot='#feb24c'
                    start_vertical_cut=465
                    end_vertical_cut=start_vertical_cut+10
                    ax1.axvline(start_vertical_cut,color=color_plot)
                    ax1.axvline(end_vertical_cut,color=color_plot)
                    #If several vertical lines to plot and plot average
                    for i in range (start_vertical_cut,end_vertical_cut+1):
                        ax4.plot(radar_slice[:,i],np.arange(0,radar_slice[:,i].size),color='0.7')
                    ax4.plot(np.mean(radar_slice[:,start_vertical_cut:end_vertical_cut],1),np.arange(0,radar_slice[:,i].size),color=color_plot)
                    ax4.invert_yaxis() #Invert the y axis
                    ax4.grid()
                    ax4.set_xlim(-6,6)
                    ax4.set_ylim(radar_slice.shape[0],-5)
                    ax4.set_yticklabels([])
                    ax4.text(0.04, 0.97,'d',ha='center', va='center', transform=ax4.transAxes,fontsize=15,zorder=10,weight='bold',color='black')#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot
                    #Display sectors of interest
                    ax7.scatter(lon_3413[0][int((start_vertical_cut+end_vertical_cut)/2)],lat_3413[0][int((start_vertical_cut+end_vertical_cut)/2)],s=200,marker='.',linewidths=0,color=color_plot,zorder=10)

                    color_plot='#c51b8a'
                    start_vertical_cut=490
                    end_vertical_cut=start_vertical_cut+10
                    ax1.axvline(start_vertical_cut,color=color_plot)
                    ax1.axvline(end_vertical_cut,color=color_plot)
                    #If several vertical lines to plot and plot average
                    for i in range (start_vertical_cut,end_vertical_cut+1):
                        ax5.plot(radar_slice[:,i],np.arange(0,radar_slice[:,i].size),color='0.7')
                    ax5.plot(np.mean(radar_slice[:,start_vertical_cut:end_vertical_cut],1),np.arange(0,radar_slice[:,i].size),color=color_plot)
                    ax5.invert_yaxis() #Invert the y axis
                    ax5.grid()
                    ax5.set_xlim(-6,6)
                    ax5.set_ylim(radar_slice.shape[0],-5)
                    ax5.set_yticklabels([])
                    ax5.text(0.04, 0.97,'e',ha='center', va='center', transform=ax5.transAxes,fontsize=15,zorder=10,weight='bold',color='black')#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot
                    #Display sectors of interest
                    ax7.scatter(lon_3413[0][int((start_vertical_cut+end_vertical_cut)/2)],lat_3413[0][int((start_vertical_cut+end_vertical_cut)/2)],s=200,marker='.',linewidths=0,color=color_plot,zorder=10)

                    color_plot='#54278f'
                    start_vertical_cut=800
                    end_vertical_cut=start_vertical_cut+10
                    ax1.axvline(start_vertical_cut,color=color_plot)
                    ax1.axvline(end_vertical_cut,color=color_plot)
                    #If several vertical lines to plot and plot average
                    for i in range (start_vertical_cut,end_vertical_cut+1):
                        ax6.plot(radar_slice[:,i],np.arange(0,radar_slice[:,i].size),color='0.7')
                    ax6.plot(np.mean(radar_slice[:,start_vertical_cut:end_vertical_cut],1),np.arange(0,radar_slice[:,i].size),color=color_plot)
                    ax6.invert_yaxis() #Invert the y axis
                    ax6.grid()
                    ax6.set_xlim(-6,6)
                    ax6.set_ylim(radar_slice.shape[0],-5)
                    ax6.set_yticklabels([])
                    ax6.text(0.04, 0.97,'f',ha='center', va='center', transform=ax6.transAxes,fontsize=15,zorder=10,weight='bold',color='black')#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot
                    #Display sectors of interest
                    ax7.scatter(lon_3413[0][int((start_vertical_cut+end_vertical_cut)/2)],lat_3413[0][int((start_vertical_cut+end_vertical_cut)/2)],s=200,marker='.',linewidths=0,color=color_plot,zorder=10)

                    #Finalize
                    ax2.set_yticks(ticks_yplot) 
                    ax2.set_yticklabels(np.round(depths[ticks_yplot]))                    
                    ax2.set_ylabel('Depth [m]')
                    ax4.set_xlabel('Radar signal strength [dB]')
                                        
                    ### ------------------ Display map ------------------- ###
                    '''
                    #Display regions
                    SW_rignotetal.plot(ax=ax7,color='white', edgecolor='black',linewidth=0.5) 
                    CW_rignotetal.plot(ax=ax7,color='white', edgecolor='black',linewidth=0.5)
                    '''
                    #Display contours
                    cont=ax7.contour(GrIS_DEM_display_SW[:,:], levels=np.arange(1500,2750,250), extent=extent_DEM_SW, transform=crs, origin='upper', colors=['#8c510a'],linewidth=0.25)
                    
                    #Display 2010-2018 high end ice slabs
                    ax7.scatter(df_2010_2018_high['lon_3413'],df_2010_2018_high['lat_3413'],s=7,marker='.',linewidths=0,color='#4575b4')

                    #Display the 2002-2003 trace of interest
                    ax7.scatter(lon_3413,lat_3413,s=7,marker='.',linewidths=0,color='black')

                    #Show KAN_U
                    ax7.scatter(-89205.404,-2522571.489,s=15,c='#b2182b',zorder=10)
                    
                    #Add legend
                    from matplotlib.patches import Patch
                    from matplotlib.lines import Line2D

                    #Custom legend myself,  line2D from https://stackoverflow.com/questions/39500265/how-to-manually-create-a-legend, marker from https://stackoverflow.com/questions/47391702/how-to-make-a-colored-markers-legend-from-scratch
                    legend_elements = [Line2D([0], [0], label='Elevation contours', color='#8c510a'),
                                       Line2D([0], [0], label='Radargram 02/06/2002', color='k'),
                                       Line2D([0], [0], marker='o', linestyle='none', label='2010-2018 ice slabs', color='#4575b4'),
                                       Line2D([0], [0], marker='o', linestyle='none', label='KAN_U', color='#b2182b')]

                    ax7.legend(handles=legend_elements,loc='lower center',fontsize=10)
                    
                    '''
                    lgnd = ax7.legend(loc='lower center',fontsize=10)
                    # Plot legend. This is from https://stackoverflow.com/questions/24706125/setting-a-fixed-size-for-points-in-legend
                    lgnd.legendHandles[0]._sizes = [100]
                    lgnd.legendHandles[1]._sizes = [100]
                    lgnd.legendHandles[2]._sizes = [20]
                    #ax7.gridlines(color='gray', linestyle='--')
                    #ax7.coastlines()
                    '''
                    #Set x lims
                    ax7.set_xlim(-187234.32225556672, -62652.03297527041)
                    ax7.set_ylim(-2757515.136754996, -2485033.0174200167)
                    
                    #Display lat/lon lines in map
                    gl=ax7.gridlines(draw_labels=True, xlocs=[-47, -48, -49], ylocs=[65, 66, 67], x_inline=False, y_inline=False,linewidth=0.5)
                    #Customize lat and lon labels
                    gl.ylabels_left = False
                    gl.xlabels_top = False
                    
                    #Add elevation contour values
                    ax7.text(0.24, 0.72,'1500', ha='center', va='center', rotation=70,transform=ax7.transAxes,fontsize=10,color='#8c510a')#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot
                    ax7.text(0.48, 0.74,'1750', ha='center', va='center', rotation=70,transform=ax7.transAxes,fontsize=10,color='#8c510a')#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot
                    ax7.text(0.85, 0.74,'2000', ha='center', va='center', rotation=70,transform=ax7.transAxes,fontsize=10,color='#8c510a')#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot
                    ax7.text(0.89, 0.35,'2250', ha='center', va='center', rotation=60,transform=ax7.transAxes,fontsize=10,color='#8c510a')#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot
                    
                    ax7.text(0.04, 0.97,'g',ha='center', va='center', transform=ax7.transAxes,fontsize=15,zorder=10,weight='bold',color='black')#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot
                    ### ------------------ Display map ------------------- ###

                    plt.show()
                    pdb.set_trace()
                    
                    #Save the figure
                    plt.savefig('C:/Users/jullienn/switchdrive/Private/research/RT1/figures/S2/v3/figS2.png',dpi=500)
                    
                    
print('End of processing')