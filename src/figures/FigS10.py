# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 10:04:37 2022

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

import pickle
import scipy.io
import numpy as np
import pdb
import h5py
import matplotlib.pyplot as plt
import pandas as pd
import geopandas as gpd
import matplotlib.gridspec as gridspec
from pyproj import Transformer
import seaborn as sns
sns.set_theme(style="whitegrid")
import matplotlib.image as mpimg
from pyproj import Transformer

#Define palette for time periods, this is from fig2_paper_icelsabs.py
#This is from https://www.python-graph-gallery.com/33-control-colors-of-boxplot-seaborn
#my_pal = {'2010': "#fdd49e", '2011': "#fc8d59", '2012': "#fc8d59", '2013':"#d7301f",'2014':"#d7301f",'2017':"#7f0000",'2018':"#7f0000"}
my_pal = {'2010': "#9ecae1", '2011': "#6baed6", '2012': "#6baed6", '2013':"#3182bd", '2014':"#3182bd", '2017':"#d73027", '2018':"#d73027"}

### -------------------------- Load shapefiles --------------------------- ###
#Load Rignot et al., 2016 Greenland drainage bassins
path_rignotetal2016_GrIS_drainage_bassins='C:/Users/jullienn/switchdrive/Private/research/backup_Aglaja/working_environment/greenland_topo_data/GRE_Basins_IMBIE2_v1.3/'
GrIS_drainage_bassins=gpd.read_file(path_rignotetal2016_GrIS_drainage_bassins+'GRE_Basins_IMBIE2_v1.3.shp',rows=slice(51,57,1)) #the regions are the last rows of the shapefile

#Extract indiv regions and create related indiv shapefiles
SW_rignotetal=GrIS_drainage_bassins[GrIS_drainage_bassins.SUBREGION1=='SW']
CW_rignotetal=GrIS_drainage_bassins[GrIS_drainage_bassins.SUBREGION1=='CW']
### -------------------------- Load shapefiles --------------------------- ###

#Define transformer for coordinates transform from "EPSG:4326" to "EPSG:3413"
transformer = Transformer.from_crs("EPSG:4326", "EPSG:3413", always_xy=True)

panel_a={2010:'empty',
      2011:['Data_20110516_01_009.mat','Data_20110516_01_010.mat'],
      2012:'empty',
      2013:['Data_20130402_01_008.mat'],
      2014:'empty',
      2017:['Data_20170412_01_075.mat','Data_20170412_01_076.mat'],
      2018:'empty'}

panel_b={2010:['Data_20100517_02_001.mat','Data_20100517_02_002.mat'],
      2011:['Data_20110502_01_171.mat'],
      2012:['Data_20120516_01_002.mat'],
      2013:['Data_20130419_01_004.mat','Data_20130419_01_005.mat'],
      2014:['Data_20140507_03_007.mat','Data_20140507_03_008.mat'], #test with 20140514_02_087_089 and 20140515_02_173_175 also
      2017:['Data_20170417_01_171.mat','Data_20170417_01_172.mat','Data_20170417_01_173.mat','Data_20170417_01_174.mat'],
      2018:'empty'}

'''
panelb:
2010  is reversed
2011  is reversed
2014  is reversed
2017  is reversed
'''

panel_c={2010:['Data_20100507_01_008.mat','Data_20100507_01_009.mat','Data_20100507_01_010.mat'],
      2011:['Data_20110426_01_009.mat','Data_20110426_01_010.mat','Data_20110426_01_011.mat'],
      2012:'empty',
      2013:'empty',
      2014:['Data_20140421_01_009.mat','Data_20140421_01_010.mat','Data_20140421_01_011.mat','Data_20140421_01_012.mat','Data_20140421_01_013.mat'],
      2017:['Data_20170424_01_008.mat','Data_20170424_01_009.mat','Data_20170424_01_010.mat','Data_20170424_01_011.mat','Data_20170424_01_012.mat','Data_20170424_01_013.mat','Data_20170424_01_014.mat'],
      2018:'empty'}

panel_d={2010:['Data_20100508_01_114.mat','Data_20100508_01_115.mat'],
      2011:['Data_20110419_01_008.mat','Data_20110419_01_009.mat','Data_20110419_01_010.mat'],
      2012:['Data_20120418_01_129.mat','Data_20120418_01_130.mat','Data_20120418_01_131.mat'],
      2013:['Data_20130405_01_165.mat','Data_20130405_01_166.mat','Data_20130405_01_167.mat'],
      2014:['Data_20140424_01_002.mat','Data_20140424_01_003.mat','Data_20140424_01_004.mat'],
      2017:['Data_20170422_01_168.mat','Data_20170422_01_169.mat','Data_20170422_01_170.mat','Data_20170422_01_171.mat'],
      2018:['Data_20180427_01_170.mat','Data_20180427_01_171.mat','Data_20180427_01_172.mat']}

'''
panel d
2010  is reversed
2012  is reversed
2013  is reversed
2017  is reversed
2018  is reversed
'''
#This one is collocated with FS1, 2, 3.
panel_e={2010:'empty',
      2011:'empty',
      2012:['Data_20120423_01_137.mat','Data_20120423_01_138.mat'],
      2013:['Data_20130409_01_010.mat','Data_20130409_01_011.mat','Data_20130409_01_012.mat'],
      2014:'empty',
      2017:'empty',
      2018:['Data_20180421_01_004.mat','Data_20180421_01_005.mat','Data_20180421_01_006.mat','Data_20180421_01_007.mat']}

'''panel_e
2012  is reversed
'''

panel_f={2010:['Data_20100513_01_001.mat','Data_20100513_01_002.mat'],
      2011:['Data_20110411_01_116.mat','Data_20110411_01_117.mat','Data_20110411_01_118.mat'],
      2012:['Data_20120428_01_125.mat','Data_20120428_01_126.mat'],
      2013:'empty',
      2014:['Data_20140408_11_024.mat','Data_20140408_11_025.mat','Data_20140408_11_026.mat'],
      2017:['Data_20170508_02_165.mat','Data_20170508_02_166.mat','Data_20170508_02_167.mat','Data_20170508_02_168.mat','Data_20170508_02_169.mat','Data_20170508_02_170.mat','Data_20170508_02_171.mat'],
      2018:'empty'}

'''
panel_f
2011  is reversed
2012  is reversed
2014  is reversed
2017  is reversed
'''
#Define the panel to study
investigation_year=panel_f

plt.rcParams.update({'font.size': 20})

fig = plt.figure()
gs = gridspec.GridSpec(30, 101)
gs.update(wspace=0.1)
gs.update(wspace=0.5)


if (investigation_year==panel_a):
    ax_sector = plt.subplot(gs[0:1, 0:100])
    ax2 = plt.subplot(gs[1:5, 0:100])
    ax4 = plt.subplot(gs[5:9, 0:100])
    ax6 = plt.subplot(gs[9:13, 0:100])
    axc = plt.subplot(gs[1:13, 100:101])
elif (investigation_year==panel_b):
    ax_sector = plt.subplot(gs[0:1, 0:100])
    ax1 = plt.subplot(gs[1:5, 0:100])
    ax2 = plt.subplot(gs[5:9, 0:100])
    ax3 = plt.subplot(gs[9:13, 0:100])
    ax4 = plt.subplot(gs[13:17, 0:100])
    ax5 = plt.subplot(gs[17:21, 0:100])
    ax6 = plt.subplot(gs[21:25, 0:100])
    axc = plt.subplot(gs[1:25, 100:101])
elif (investigation_year==panel_c):
    ax_sector = plt.subplot(gs[0:1, 0:100])
    ax1 = plt.subplot(gs[1:5, 0:100])
    ax2 = plt.subplot(gs[5:9, 0:100])
    ax5 = plt.subplot(gs[9:13, 0:100])
    ax6 = plt.subplot(gs[13:17, 0:100])
    axc = plt.subplot(gs[1:17, 100:101])
elif (investigation_year==panel_d):
    ax_sector = plt.subplot(gs[0:1, 0:100])
    ax1 = plt.subplot(gs[1:5, 0:100])
    ax2 = plt.subplot(gs[5:9, 0:100])
    ax3 = plt.subplot(gs[9:13, 0:100])
    ax4 = plt.subplot(gs[13:17, 0:100])
    ax5 = plt.subplot(gs[17:21, 0:100])
    ax6 = plt.subplot(gs[21:25, 0:100])
    ax7 = plt.subplot(gs[25:29, 0:100])
    axc = plt.subplot(gs[1:29, 100:101])
elif (investigation_year==panel_e):
    ax_sector = plt.subplot(gs[0:1, 0:100])
    ax3 = plt.subplot(gs[1:5, 0:100])
    ax4 = plt.subplot(gs[5:9, 0:100])
    ax7 = plt.subplot(gs[9:13, 0:100])
    axc = plt.subplot(gs[1:13, 100:101])
elif (investigation_year==panel_f):
    ax_sector = plt.subplot(gs[0:1, 0:100])
    ax1 = plt.subplot(gs[1:5, 0:100])
    ax2 = plt.subplot(gs[5:9, 0:100])
    ax3 = plt.subplot(gs[9:13, 0:100])
    ax5 = plt.subplot(gs[13:17, 0:100])
    ax6 = plt.subplot(gs[17:21, 0:100])
    axc = plt.subplot(gs[1:21, 100:101])
else:
    print('Wrong transect name input')

#Define paths
path_data_Jullien='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2002_2018/i_out_from_IceBridgeGPR_Manager_v2.py/pickles_and_images/'
path_data='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/'

#Compute the speed (Modified Robin speed):
# self.C / (1.0 + (coefficient*density_kg_m3/1000.0))
v= 299792458

dataframe={}

#lon_for_min_max
lon_for_min_max=[]

for single_year in investigation_year.keys():
    #print(single_year)
        
    #If no data, continue
    if (investigation_year[single_year]=='empty'):
        #print('No data for year '+str(single_year)+', continue')
        continue
    
    ###1. Load the images:
    start_date_track=investigation_year[single_year][0]
    end_date_track=investigation_year[single_year][-1]
    date_track=start_date_track[5:20]+'_'+end_date_track[17:20]
    
    #pat_to_use=path_data_MacFerrin+date_track+'_1'
    pat_to_use=path_data_Jullien+date_track+'_X'
    
    filename_image=pat_to_use+'DEPTHCORRECT_AFTER.png'
    #Load image
    img = mpimg.imread(filename_image)
    
    filename_depth_corrected=path_data_Jullien+'Depth_Corrected_Picklefiles/'+date_track+'_Depth_CORRECTED.pickle'
    #Open files
    f_depth_corrected = open(filename_depth_corrected, "rb")
    depth_corr = pickle.load(f_depth_corrected)
    f_depth_corrected.close()

    ###2. Load the latitude and longitude
    lat_appended=[]
    lon_appended=[]
    
    for indiv_file_load in investigation_year[single_year]:
        
        #Create the path
        path_raw_data=path_data+str(single_year)+'_Greenland_P3/CSARP_qlook/'+indiv_file_load[5:16]+'/'

        #Load data
        if (single_year>=2014):
            
            fdata_filename = h5py.File(path_raw_data+indiv_file_load)
            lat_filename=fdata_filename['Latitude'][:,:]
            lon_filename=fdata_filename['Longitude'][:,:]
            time_filename=fdata_filename['Time'][:,:]
            
        else:
            fdata_filename = scipy.io.loadmat(path_raw_data+indiv_file_load)
            lat_filename = fdata_filename['Latitude']
            lon_filename = fdata_filename['Longitude']
            time_filename = fdata_filename['Time']
            
        #Append data
        lat_appended=np.append(lat_appended,lat_filename)
        lon_appended=np.append(lon_appended,lon_filename)
        
    #Check whether the data are acquired ascending or descending elevation wise.
    #I choose the ascending format. For the data that are descending, reverse them
    #To have ascending data, the longitude should increase
    #(-48 = low elevations, -46 = higher elevations)
    
    if (np.sum(np.diff(lon_appended))<0):
        #It is going toward lower elevations, thus flip left-right
        #(or up-down) all the data!
        print(single_year,' is reversed')
        lat_appended=np.flipud(lat_appended)
        lon_appended=np.flipud(lon_appended)
        img=np.fliplr(img)
        depth_corr=np.fliplr(depth_corr)
    
    #Calculate the depth from the time
    #########################################################################
    # From plot_2002_2003.py - BEGIN
    #########################################################################
    depth_check = v * time_filename / 2.0
    
    #If 2014, transpose the vector
    if (str(single_year)>='2014'):
        depth_check=np.transpose(depth_check)
    
    #Reset times to zero! This is from IceBridgeGPR_Manager_v2.py
    if (depth_check[10]<0):
        #depth_check[10] so that I am sure that the whole vector is negative and
        #not the first as can be for some date were the proccessing is working
        depth_check=depth_check+abs(depth_check[0])
        depth = depth_check
    else:
        depth = depth_check
    
    if (str(single_year) in list(['2011','2012','2014','2017','2018'])):
        if (depth_check[10]>1):
            #depth_check[10] so that I am sure that the whole vector is largely positive and
            #not the first as can be for some date were the proccessing is working
            depth_check=depth_check-abs(depth_check[0])
            depth = depth_check
    
    #Transform the coordinates from WGS84 to EPSG:3413
    #Example from: https://pyproj4.github.io/pyproj/stable/examples.html
    points=transformer.transform(np.array(lon_appended),np.array(lat_appended))
    lon_3413=points[0]
    lat_3413=points[1]
    
    #Calculate distances
    distances_with_start_transect=compute_distances(lon_3413,lat_3413)
    
    #Store reunited lat/lon, slice output and mask in a dictionnary:
    dataframe[str(single_year)]={'lat_appended':lat_appended,
                                 'lon_appended':lon_appended,
                                 'lat_3413':lat_3413,
                                 'lon_3413':lon_3413,
                                 'distances':distances_with_start_transect,
                                 'depth':depth,
                                 'img':img,
                                 'depth_corr':depth_corr}
    
    #Append for lon min/max computation
    lon_for_min_max=np.append(lon_for_min_max,lon_appended)

#Define min and max
x_start=np.min(lon_for_min_max)
x_end=np.max(lon_for_min_max)

for single_year in investigation_year.keys():
    
    if (investigation_year[single_year]=='empty'):
        continue
    print(single_year)
    
    if (investigation_year==panel_a):
        start_transect=-64.6886
        end_transect=-63.88337773831948
        vmin_plot=-4.5
        vmax_plot=4.5
        #limits of ice slabs development sectors
        start_well_developed=-64.67344699111679
        end_well_developed=-64.5116804275458
        end_in_development=-64.07045269370637
        end_in_initiation=-64.07045269370637
        
    elif (investigation_year==panel_b):
        start_transect=-66.8557
        end_transect=-65.48369681346244
        vmin_plot=-4.5
        vmax_plot=4.5
        #limits of ice slabs development sectors
        start_well_developed=-66.57146840643074
        end_well_developed=-66.57146840643074
        end_in_development=-66.00486187318322
        end_in_initiation=-65.17987346178424
        
    elif (investigation_year==panel_c):
        start_transect=-47.566978989211904
        end_transect=-46.088871815847654
        vmin_plot=-4.5
        vmax_plot=4.5
        #limits of ice slabs development sectors
        start_well_developed=-47.85901106902736
        end_well_developed=-47.16406520774126
        end_in_development=-47.11308507589178
        end_in_initiation=-46.632074831791776
        
    elif (investigation_year==panel_d):
        start_transect=-47.70785561652585
        end_transect=-46.41555609606877
        vmin_plot=-4.5
        vmax_plot=4.5
        #limits of ice slabs development sectors
        start_well_developed=-47.82252402970681
        end_well_developed=-47.33395800638842
        end_in_development=-46.94903911890173
        end_in_initiation=-46.75989755131823
        
    elif (investigation_year==panel_e):
        start_transect=-47.42328169558814
        end_transect=-46.56546212605787    
        vmin_plot=-4.5
        vmax_plot=4.5
        #limits of ice slabs development sectors
        start_well_developed=-47.42328169558814
        end_well_developed=-47.11742506657421
        end_in_development=-47.00359606697867
        end_in_initiation= -46.87457243557124 #end ice slabs in 2013: -46.63497339965015
        
    elif (investigation_year==panel_f):
        start_transect=-48.21060856534727
        end_transect=-46.88764316176339
        vmin_plot=-4.5
        vmax_plot=4.5
        #limits of ice slabs development sectors
        start_well_developed=-48.21058671285614
        end_well_developed=-47.72093027678061
        end_in_development=-47.50105332110282
        end_in_initiation=-47.12849888596365
        
    else:
        print('Wrong transect name input')

    if (single_year==2010):
        ax_plot=ax1
        color_toplot="#4dac26"
    elif (single_year==2011):
        ax_plot=ax2
        color_toplot="#0571b0"
    elif (single_year==2012):
        ax_plot=ax3
        color_toplot="#e31a1c"
    elif (single_year==2013):
        ax_plot=ax4
        color_toplot="#e31a1c"
    elif (single_year==2014):
        ax_plot=ax5
        color_toplot="#2171b5"
    elif (single_year==2017):
        ax_plot=ax6
        color_toplot="#2171b5"
    elif (single_year==2018):
        ax_plot=ax7
        color_toplot="#e31a1c"
        ax_plot.set_xlabel('Longitude [°]')
        #Activate ticks xlabel
        ax_plot.xaxis.tick_bottom()
    else:
        print('year not know')
    
    #Load data
    X=dataframe[str(single_year)]['lon_appended']
    Y=np.arange(0,100,100/dataframe[str(single_year)]['depth_corr'].shape[0])
    C=dataframe[str(single_year)]['depth_corr']
    
    #plot data
    cb=ax_plot.pcolor(X, Y, C,cmap=plt.get_cmap('gray'),zorder=-1,vmin=vmin_plot, vmax=vmax_plot)
    ax_plot.invert_yaxis() #Invert the y axis = avoid using flipud.
    
    #Activate ticks ylabel
    ax_plot.yaxis.tick_left()
    
    #Set lims
    ax_plot.set_ylim(20,0)
    
    #Set yticklabels
    ax_plot.set_yticks([0,10,20])
    ax_plot.set_yticklabels(['0','10',''],fontsize=25)
        
    #Set transect limits
    ax_plot.set_xlim(start_transect,end_transect)
    
    #Get rid of xticklabels
    ax_plot.set_xticklabels([])
    
    #Display year
    ax_plot.text(0.97, 0.775,str(single_year),ha='center', va='center', transform=ax_plot.transAxes,weight='bold',fontsize=25,color=my_pal[str(single_year)])#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot
    
    ''' Old way of stages of ice slabs development sector display
    #Display vertical lines showing the limits of the ice slabs development sectors
    #Start of well developed sector
    if (investigation_year==panel_a):
        ax_plot.axvline(x=start_well_developed,color='green',linestyle='--',linewidth=2) 
    else:
        ax_plot.axvline(x=start_transect+0.002,color='green',linestyle='--',linewidth=2) 
    
    if (investigation_year==panel_b):
        ax_plot.axvline(x=end_transect-0.002,color='red',linestyle='--',linewidth=2)
    else:
        ax_plot.axvline(x=end_in_initiation,color='red',linestyle='--',linewidth=2)
    '''
    if (investigation_year==panel_e):
        #Display start and end of Fig. 4
        ax_plot.axvline(x=-47.196803171288636,color='green',linewidth=2)
        ax_plot.axvline(x=-46.74361490492269,color='green',linewidth=2)

    #End of well developed sector
    ax_plot.axvline(x=end_well_developed,color='white',linestyle='--',linewidth=2)    
    #End of in development sector
    
    ax_plot.axvline(x=end_in_development,color='white',linestyle='--',linewidth=2)
    #End of in initiation sector
    
if (investigation_year==panel_a):
    ax4.set_ylabel('Depth [m]',fontsize=25)
    ax6.set_yticklabels(['0','10','20'],fontsize=25)
    ticks_through=ax6.get_xticks()
    year_ticks=2017
    ax_tick_plot=ax6
    ax_top=ax2
    
elif (investigation_year==panel_b):
    ax4.set_ylabel('Depth [m]',fontsize=25)
    ax6.set_yticklabels(['0','10','20'],fontsize=25)
    ticks_through=ax6.get_xticks()
    year_ticks=2017
    ax_tick_plot=ax6
    ax_top=ax1

elif (investigation_year==panel_c):
    ax5.set_ylabel('Depth [m]',fontsize=25)
    ax6.set_yticklabels(['0','10','20'],fontsize=25)
    ticks_through=ax6.get_xticks()
    year_ticks=2017
    ax_tick_plot=ax6
    ax_top=ax1

elif (investigation_year==panel_d):
    ax4.set_ylabel('Depth [m]',fontsize=25)
    ax7.set_yticklabels(['0','10','20'],fontsize=25)
    ticks_through=ax7.get_xticks()
    year_ticks=2018
    ax_tick_plot=ax7
    ax_top=ax1
    
elif (investigation_year==panel_e):
    ax4.set_ylabel('Depth [m]',fontsize=25)
    ax7.set_yticklabels(['0','10','20'],fontsize=25)
    ticks_through=ax7.get_xticks()
    year_ticks=2018
    ax_tick_plot=ax7
    ax_top=ax3
    
elif (investigation_year==panel_f):
    ax3.set_ylabel('Depth [m]',fontsize=25)
    ax6.set_yticklabels(['0','10','20'],fontsize=25)
    ticks_through=ax6.get_xticks()
    year_ticks=2017
    ax_tick_plot=ax6
    ax_top=ax1

else:
    print('Wrong transect name input')

#Display colorbar. This is from FigS1.py
cbar_depth=fig.colorbar(cb, cax=axc, aspect=5)#aspect is from https://stackoverflow.com/questions/33443334/how-to-decrease-colorbar-width-in-matplotlib

cbar_depth.set_ticks(cbar_depth.get_ticks()[1:-1])
cbar_depth.set_ticklabels(cbar_depth.get_ticks(),fontsize=25)
cbar_depth.set_label('Radar signal strength [dB]',fontsize=25)

plot_dist=[]
for indiv_tick in ticks_through:
    lon_diff=[]
    lon_diff=np.abs(dataframe[str(year_ticks)]['lon_appended']-indiv_tick)
    index_min=np.argmin(lon_diff)
    if (lon_diff[index_min]>0.2):
        plot_dist=np.append(plot_dist,999)
    else:
        plot_dist=np.append(plot_dist,dataframe[str(year_ticks)]['distances'][index_min]/1000-dataframe[str(year_ticks)]['distances'][np.argmin(np.abs(dataframe[str(year_ticks)]['lon_appended']-start_transect))]/1000)


ax_tick_plot.xaxis.set_ticks_position('bottom') 
ax_tick_plot.set_xticklabels(np.round(plot_dist).astype(int),fontsize=25)
ax_tick_plot.set_xlabel('Distance [km]',fontsize=25)

#Display stages of ice slabs development sectors
ax_sector.axvspan(start_well_developed, end_well_developed, facecolor='#000000')
ax_sector.axvspan(end_well_developed, end_in_development, facecolor='#7f7f7f')
ax_sector.axvspan(end_in_development, end_in_initiation, facecolor='#b8b8b8')
ax_sector.set_xlim(start_transect,end_transect)
ax_sector.set_ylim(0.97,1.05)
ax_sector.axis('off')


ax_top.set_title('Transect F',fontsize=25,pad=30)

figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()

plt.show()

pdb.set_trace()

#Estimate distance from start of transect fron any longitude: dataframe[str(2018)]['distances'][np.argmin(np.abs(dataframe[str(2018)]['lon_appended']+np.abs(lon_of_interest)))]-dataframe[str(2018)]['distances'][np.argmin(np.abs(dataframe[str(2018)]['lon_appended']-start_transect))]

#Save the figure
plt.savefig('C:/Users/jullienn/switchdrive/Private/research/RT1/figures/S7/v6/figS7_panelf.png',dpi=300,bbox_inches='tight')
#bbox_inches is from https://stackoverflow.com/questions/32428193/saving-matplotlib-graphs-to-image-as-full-screen)
