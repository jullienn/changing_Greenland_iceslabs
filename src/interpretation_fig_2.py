# -*- coding: utf-8 -*-
"""
Created on Thu Jan 20 10:41:00 2022

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


### -------------------------- Load shapefiles --------------------------- ###
#Load Rignot et al., 2016 Greenland drainage bassins
path_rignotetal2016_GrIS_drainage_bassins='C:/Users/jullienn/switchdrive/Private/research/backup_Aglaja/working_environment/greenland_topo_data/GRE_Basins_IMBIE2_v1.3/'
GrIS_drainage_bassins=gpd.read_file(path_rignotetal2016_GrIS_drainage_bassins+'GRE_Basins_IMBIE2_v1.3.shp',rows=slice(51,57,1)) #the regions are the last rows of the shapefile

#Extract indiv regions and create related indiv shapefiles
SW_rignotetal=GrIS_drainage_bassins[GrIS_drainage_bassins.SUBREGION1=='SW']
CW_rignotetal=GrIS_drainage_bassins[GrIS_drainage_bassins.SUBREGION1=='CW']
### -------------------------- Load shapefiles --------------------------- ###
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

panel_b_new=

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
investigation_year=panel_b_new

fig = plt.figure()
gs = gridspec.GridSpec(30, 10)
gs.update(wspace=0.1)
#gs.update(wspace=0.001)
ax1 = plt.subplot(gs[0:4, 0:10])
ax2 = plt.subplot(gs[4:8, 0:10])
ax3 = plt.subplot(gs[8:12, 0:10])
ax4 = plt.subplot(gs[12:16, 0:10])
ax5 = plt.subplot(gs[16:20, 0:10])
ax6 = plt.subplot(gs[20:24, 0:10])
ax7 = plt.subplot(gs[24:28, 0:10])

#Fig_prob
fig_prob = plt.figure()
gs = gridspec.GridSpec(30, 10)
gs.update(wspace=0.1)
#gs.update(wspace=0.001)
ax1_prob = plt.subplot(gs[0:4, 0:10])
ax2_prob = plt.subplot(gs[4:8, 0:10])
ax3_prob = plt.subplot(gs[8:12, 0:10])
ax4_prob = plt.subplot(gs[12:16, 0:10])
ax5_prob = plt.subplot(gs[16:20, 0:10])
ax6_prob = plt.subplot(gs[20:24, 0:10])
ax7_prob = plt.subplot(gs[24:28, 0:10])

#Fig_raw
fig_raw = plt.figure()
gs = gridspec.GridSpec(30, 10)
gs.update(wspace=0.1)
#gs.update(wspace=0.001)
ax1_raw = plt.subplot(gs[0:4, 0:10])
ax2_raw = plt.subplot(gs[4:8, 0:10])
ax3_raw = plt.subplot(gs[8:12, 0:10])
ax4_raw = plt.subplot(gs[12:16, 0:10])
ax5_raw = plt.subplot(gs[16:20, 0:10])
ax6_raw = plt.subplot(gs[20:24, 0:10])
ax7_raw = plt.subplot(gs[24:28, 0:10])

#Define paths
path_data_MacFerrin='C:/Users/jullienn/switchdrive/Private/research/RT1/MacFerrin_FigShare/IceBridge_Accumulation_Radar/IceBridge_Accumulation_Radar/Images_PostProcessed_Greyscale/'
path_data_Jullien='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2010_2018/i_out_from_IceBridgeGPR_Manager_v2.py/pickles_and_images/'
path_data_Jullien_probability='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2010_2018/iii_out_from_probabilistic_iceslabs.py/images/'
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
    
    #Load probability
    path_prob=path_data_Jullien_probability+date_track+'_probability_iceslabs_presence.png'
    #Load probability image
    img_prob = mpimg.imread(path_prob)
    
    #Load raw data after surface picking
    filename_surfpick_data=path_data_Jullien+date_track+'_0m_30m_BESTFIT_V1.png'
    #Open the test file
    surfpick_data = mpimg.imread(filename_surfpick_data)

    ###2. Load the latitude and longitude
    lat_appended=[]
    lon_appended=[]
    
    for indiv_file_load in investigation_year[single_year]:
        
        #Create the path
        path_raw_data=path_data+str(single_year)+'_Greenland_P3/CSARP_qlook/'+indiv_file_load[5:16]+'/'

        # If 20180421_01_007, different path
        if (indiv_file_load=='Data_20180421_01_007.mat'):
            path_raw_data='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/2018_Greenland_P3/'
        
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
        img_prob=np.fliplr(img_prob)
        surfpick_data=np.fliplr(surfpick_data)
    
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
        
    #Store reunited lat/lon, slice output and mask in a dictionnary:
    dataframe[str(single_year)]={'lat_appended':lat_appended,
                                 'lon_appended':lon_appended,
                                 'depth':depth,
                                 'img':img,
                                 'img_prob':img_prob,
                                 'surfpick_data':surfpick_data}
    
    #Append for lon min/max computation
    lon_for_min_max=np.append(lon_for_min_max,lon_appended)

#Define min and max
x_start=np.min(lon_for_min_max)
x_end=np.max(lon_for_min_max)

for single_year in investigation_year.keys():
    
    if (investigation_year[single_year]=='empty'):
        continue
    print(single_year)
    
    if (single_year==2010):
        ax_plot=ax1
        ax_prob=ax1_prob
        ax_raw=ax1_raw
        color_toplot="#4dac26"
        ax_plot.set_xticklabels([])
    elif (single_year==2011):
        ax_plot=ax2
        ax_prob=ax2_prob
        ax_raw=ax2_raw
        color_toplot="#0571b0"
        ax_plot.set_xticklabels([])
    elif (single_year==2012):
        ax_plot=ax3
        ax_prob=ax3_prob
        ax_raw=ax3_raw
        color_toplot="#e31a1c"
        ax_plot.set_xticklabels([])
    elif (single_year==2013):
        ax_plot=ax4
        ax_prob=ax4_prob
        ax_raw=ax4_raw
        color_toplot="#e31a1c"
        ax_plot.set_xticklabels([])
    elif (single_year==2014):
        ax_plot=ax5
        ax_prob=ax5_prob
        ax_raw=ax5_raw
        color_toplot="#2171b5"
        ax_plot.set_xticklabels([])
    elif (single_year==2017):
        ax_plot=ax6
        ax_prob=ax6_prob
        ax_raw=ax6_raw
        color_toplot="#2171b5"
        ax_plot.set_xticklabels([])
    elif (single_year==2018):
        ax_plot=ax7
        ax_prob=ax7_prob
        ax_raw=ax7_raw
        color_toplot="#e31a1c"
        ax_plot.set_xlabel('Longitude [°]')
        #Activate ticks xlabel
        ax_plot.xaxis.tick_bottom()
    else:
        print('year not know')
    
    #Load data
    X=dataframe[str(single_year)]['lon_appended']
    
    Y=np.arange(0,30,30/dataframe[str(single_year)]['img'].shape[0])
    Y_prob=np.arange(0,20,20/dataframe[str(single_year)]['img_prob'].shape[0])
    Y_raw=np.arange(0,30,30/dataframe[str(single_year)]['surfpick_data'].shape[0])
    
    C=dataframe[str(single_year)]['img']
    C_prob=dataframe[str(single_year)]['img_prob']
    C_raw=dataframe[str(single_year)]['surfpick_data']
  
    #plot data
    cb=ax_plot.pcolor(X, Y, C,cmap=plt.get_cmap('gray'),zorder=-1)#,norm=divnorm)
    ax_plot.invert_yaxis() #Invert the y axis = avoid using flipud.  
    
    ax_prob.pcolor(X, Y_prob, C_prob,cmap=plt.get_cmap('gray'),zorder=-1)#,norm=divnorm)
    ax_prob.invert_yaxis() #Invert the y axis = avoid using flipud.
    
    ax_raw.pcolor(X, Y_raw, C_raw,cmap=plt.get_cmap('gray'),zorder=-1)#,norm=divnorm)
    ax_raw.invert_yaxis() #Invert the y axis = avoid using flipud.
        
    #Set xlim
    ax_plot.set_xlim(x_start,x_end)
    ax_prob.set_xlim(x_start,x_end)
    ax_raw.set_xlim(x_start,x_end)
    
    #Set ylim
    ax_plot.set_ylim(20,0)
    ax_raw.set_ylim(20,0)
    
    #Activate ticks ylabel
    ax_plot.yaxis.tick_left()
    ax_prob.yaxis.tick_left()

ax1.set_ylabel('Depth [m]')
ax1_prob.set_ylabel('Depth [m]')

figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()

plt.show()

pdb.set_trace()
##############################################################################
###                              Show dry firn                             ###
##############################################################################

path_higher=path_data+'fig_2_investigation/'

list_added=list(['Data_20110516_01_011.mat','Data_20130402_01_009.mat','Data_20170412_01_077.mat',
                'Data_20100517_02_003.mat','Data_20110502_01_170.mat','Data_20120516_01_003.mat','Data_20130419_01_006.mat','Data_20140507_03_006.mat','Data_20170417_01_170.mat',
                'Data_20100507_01_011.mat','Data_20110426_01_012.mat','Data_20140421_01_014.mat','Data_20170424_01_015.mat',
                'Data_20100508_01_113.mat','Data_20110419_01_011.mat','Data_20120418_01_128.mat','Data_20130405_01_164.mat','Data_20140424_01_005.mat','Data_20170422_01_167.mat','Data_20180427_01_169.mat',
                'Data_20120423_01_136.mat','Data_20130409_01_013.mat','Data_20180421_01_008.mat',
                'Data_20100513_01_003.mat','Data_20110411_01_115.mat','Data_20120428_01_124.mat','Data_20140408_11_023.mat','Data_20170508_02_164.mat'])

panel_a={2010:'empty',
      2011:['Data_20110516_01_009.mat','Data_20110516_01_010.mat','Data_20110516_01_011.mat'],
      2012:'empty',
      2013:['Data_20130402_01_008.mat','Data_20130402_01_009.mat'],
      2014:'empty',
      2017:['Data_20170412_01_076.mat','Data_20170412_01_077.mat'],
      2018:'empty'}

panel_b={2010:['Data_20100517_02_001.mat','Data_20100517_02_002.mat','Data_20100517_02_003.mat'],
      2011:['Data_20110502_01_170.mat','Data_20110502_01_171.mat'],
      2012:['Data_20120516_01_002.mat','Data_20120516_01_003.mat'],
      2013:['Data_20130419_01_004.mat','Data_20130419_01_005.mat','Data_20130419_01_006.mat'],
      2014:['Data_20140507_03_006.mat','Data_20140507_03_007.mat','Data_20140507_03_008.mat'], #test with 20140514_02_087_089 and 20140515_02_173_175 also
      2017:['Data_20170417_01_170.mat','Data_20170417_01_171.mat','Data_20170417_01_172.mat','Data_20170417_01_173.mat','Data_20170417_01_174.mat'],
      2018:'empty'}

panel_c={2010:['Data_20100507_01_008.mat','Data_20100507_01_009.mat','Data_20100507_01_010.mat'],
      2011:['Data_20110426_01_009.mat','Data_20110426_01_010.mat','Data_20110426_01_011.mat','Data_20110426_01_012.mat'],
      2012:'empty',
      2013:'empty',
      2014:['Data_20140421_01_009.mat','Data_20140421_01_010.mat','Data_20140421_01_011.mat','Data_20140421_01_012.mat','Data_20140421_01_013.mat','Data_20140421_01_014.mat'],
      2017:['Data_20170424_01_008.mat','Data_20170424_01_009.mat','Data_20170424_01_010.mat','Data_20170424_01_011.mat','Data_20170424_01_012.mat','Data_20170424_01_013.mat','Data_20170424_01_014.mat','Data_20170424_01_015.mat'],
      2018:'empty'}

panel_d={2010:['Data_20100508_01_113.mat','Data_20100508_01_114.mat','Data_20100508_01_115.mat'],
      2011:['Data_20110419_01_008.mat','Data_20110419_01_009.mat','Data_20110419_01_010.mat','Data_20110419_01_011.mat'],
      2012:['Data_20120418_01_128.mat','Data_20120418_01_129.mat','Data_20120418_01_130.mat','Data_20120418_01_131.mat'],
      2013:['Data_20130405_01_164.mat','Data_20130405_01_165.mat','Data_20130405_01_166.mat','Data_20130405_01_167.mat'],
      2014:['Data_20140424_01_002.mat','Data_20140424_01_003.mat','Data_20140424_01_004.mat','Data_20140424_01_005.mat'],
      2017:['Data_20170422_01_167.mat','Data_20170422_01_168.mat','Data_20170422_01_169.mat','Data_20170422_01_170.mat','Data_20170422_01_171.mat'],
      2018:['Data_20180427_01_169.mat','Data_20180427_01_170.mat','Data_20180427_01_171.mat','Data_20180427_01_172.mat']}

#This one is collocated with FS1, 2, 3.
panel_e={2010:'empty',
      2011:'empty',
      2012:['Data_20120423_01_136.mat','Data_20120423_01_137.mat','Data_20120423_01_138.mat'],
      2013:['Data_20130409_01_010.mat','Data_20130409_01_011.mat','Data_20130409_01_012.mat','Data_20130409_01_013.mat'],
      2014:'empty',
      2017:'empty',
      2018:['Data_20180421_01_004.mat','Data_20180421_01_005.mat','Data_20180421_01_006.mat','Data_20180421_01_007.mat','Data_20180421_01_008.mat']}

panel_f={2010:['Data_20100513_01_001.mat','Data_20100513_01_002.mat','Data_20100513_01_003.mat'],
      2011:['Data_20110411_01_115.mat','Data_20110411_01_116.mat','Data_20110411_01_117.mat','Data_20110411_01_118.mat'],
      2012:['Data_20120428_01_124.mat','Data_20120428_01_125.mat','Data_20120428_01_126.mat'],
      2013:'empty',
      2014:['Data_20140408_11_023.mat','Data_20140408_11_024.mat','Data_20140408_11_025.mat','Data_20140408_11_026.mat'],
      2017:['Data_20170508_02_164.mat','Data_20170508_02_165.mat','Data_20170508_02_166.mat','Data_20170508_02_167.mat','Data_20170508_02_168.mat','Data_20170508_02_169.mat','Data_20170508_02_170.mat','Data_20170508_02_171.mat'],
      2018:'empty'}

panel_a_surfpick={2011:[1500],
                  2013:[1668],
                  2017:[600]}
panel_b_surfpick={2010:[2300],
                  2011:[3360],
                  2012:[1570],
                  2013:[1800],
                  2014:[1400],
                  2017:[2500]}
panel_c_surfpick={2010:[1996],
                  2011:[1850],
                  2014:[910],
                  2017:[3480]}
panel_d_surfpick={2010:[1823],
                  2011:[1877],
                  2012:[1550],
                  2013:[1620],
                  2014:[900],
                  2017:[1900],
                  2018:[870]}
panel_e_surfpick={2012:[1540],
                  2013:[1620],
                  2018:[360]}
panel_f_surfpick={2010:[3660],
                  2011:[1240],
                  2012:[1560],
                  2014:[900],
                  2017:[485]}

#Fig_elongated
fig_l = plt.figure()
gs = gridspec.GridSpec(30, 10)
gs.update(wspace=0.1)
#gs.update(wspace=0.001)
ax1_l = plt.subplot(gs[0:4, 0:10])
ax2_l = plt.subplot(gs[4:8, 0:10])
ax3_l = plt.subplot(gs[8:12, 0:10])
ax4_l = plt.subplot(gs[12:16, 0:10])
ax5_l = plt.subplot(gs[16:20, 0:10])
ax6_l = plt.subplot(gs[20:24, 0:10])
ax7_l = plt.subplot(gs[24:28, 0:10])

investigation_year=panel_b
investigation_year_surfpick=panel_b_surfpick
lon_for_min_max=[]

#Load appended radar for additional check about ice slabs continuity
for single_year in investigation_year.keys():

    if (investigation_year[single_year]=='empty'):
        continue
    print(single_year)
    
    radar_appended=[]
    lon_appended=[]
    lon_bars=[]
    
    for indiv_file_load in investigation_year[single_year]:
        print(indiv_file_load)
        indiv_radar=[]
        indiv_lon=[]
        
        if (indiv_file_load in list(list_added)):
            #Create the path
            path_raw_data=path_higher+'panel_b/'
        else:
            #Create the path
            path_raw_data=path_data+str(single_year)+'_Greenland_P3/CSARP_qlook/'+indiv_file_load[5:16]+'/'
                
        #Load data
        if (single_year>=2014):
            fdata_filename = h5py.File(path_raw_data+indiv_file_load)
            indiv_radar=np.transpose(fdata_filename['Data'][:,:])
            indiv_lon=fdata_filename['Longitude'][:,:]
        else:
            fdata_filename = scipy.io.loadmat(path_raw_data+indiv_file_load)
            indiv_radar = fdata_filename['Data']
            indiv_lon=fdata_filename['Longitude']
        if (len(radar_appended)==0):
            radar_appended=indiv_radar
            lon_appended=indiv_lon
        else:
            radar_appended=np.append(radar_appended,indiv_radar,axis=1)
            lon_appended=np.append(lon_appended,indiv_lon)
        
        #Extract lon[-1]
        lon_bars=np.append(lon_bars,indiv_lon[0][-1])
    
    #Store data
    dataframe[str(single_year)]['radar_appended']=radar_appended
    '''
    C=radar_appended
    Y=dataframe[str(single_year)]['depth']
    #Identify surface pick
    #Plot data
    figraw, (ax1raw) = plt.subplots(1, 1)#, gridspec_kw={'width_ratios': [1, 3]})
    figraw.suptitle(str(single_year))
    
    start_display=0
    end_display=2000
    
    cb=ax1raw.pcolor(lon_appended, Y[start_display:end_display,], np.log10(C[start_display:end_display,]),cmap=plt.get_cmap('gray'),zorder=-1)#,norm=divnorm)
    plt.show()
    '''
    #Surface pickup
    suggested_pixel=investigation_year_surfpick[single_year][0]    
    surface_indices=kernel_function(radar_appended, suggested_pixel)

    #Get our slice (30 meters as currently set)
    radar_slice, bottom_indices = _return_radar_slice_given_surface(radar_appended,
                                                    dataframe[str(single_year)]['depth'],
                                                    surface_indices,
                                                    meters_cutoff_above=0,
                                                    meters_cutoff_below=20)
        
    if (single_year==2010):
        ax_plot=ax1_l
        ax_plot.set_xticklabels([])
    elif (single_year==2011):
        ax_plot=ax2_l
        ax_plot.set_xticklabels([])
    elif (single_year==2012):
        ax_plot=ax3_l
        ax_plot.set_xticklabels([])
    elif (single_year==2013):
        ax_plot=ax4_l
        ax_plot.set_xticklabels([])
    elif (single_year==2014):
        ax_plot=ax5_l
        ax_plot.set_xticklabels([])
    elif (single_year==2017):
        ax_plot=ax6_l
        ax_plot.set_xticklabels([])
    elif (single_year==2018):
        ax_plot=ax7_l
        ax_plot.set_xlabel('Longitude [°]')
        #Activate ticks xlabel
        ax_plot.xaxis.tick_bottom()
    else:
        print('year not know')
        
    cb=ax_plot.pcolor(lon_appended, np.arange(0,radar_slice.shape[0],1), np.log10(radar_slice),cmap=plt.get_cmap('gray'),zorder=-1)#,norm=divnorm)    
    for xc in lon_bars:
        ax_plot.axvline(x=xc)
    ax_plot.invert_yaxis()
    #Append for lon min/max computation
    lon_for_min_max=np.append(lon_for_min_max,lon_appended)

'''
lon_for_min_max=[]
lon_for_min_max=[-66.86,-66.63]
'''
#Set xmin and xmax
ax1_l.set_xlim(np.min(lon_for_min_max),np.max(lon_for_min_max))
ax2_l.set_xlim(np.min(lon_for_min_max),np.max(lon_for_min_max))
ax3_l.set_xlim(np.min(lon_for_min_max),np.max(lon_for_min_max))
ax4_l.set_xlim(np.min(lon_for_min_max),np.max(lon_for_min_max))
ax5_l.set_xlim(np.min(lon_for_min_max),np.max(lon_for_min_max))
ax6_l.set_xlim(np.min(lon_for_min_max),np.max(lon_for_min_max))
ax7_l.set_xlim(np.min(lon_for_min_max),np.max(lon_for_min_max))

ax1.set_xlim(np.min(lon_for_min_max),np.max(lon_for_min_max))
ax2.set_xlim(np.min(lon_for_min_max),np.max(lon_for_min_max))
ax3.set_xlim(np.min(lon_for_min_max),np.max(lon_for_min_max))
ax4.set_xlim(np.min(lon_for_min_max),np.max(lon_for_min_max))
ax5.set_xlim(np.min(lon_for_min_max),np.max(lon_for_min_max))
ax6.set_xlim(np.min(lon_for_min_max),np.max(lon_for_min_max))
ax7.set_xlim(np.min(lon_for_min_max),np.max(lon_for_min_max))

ax1_prob.set_xlim(np.min(lon_for_min_max),np.max(lon_for_min_max))
ax2_prob.set_xlim(np.min(lon_for_min_max),np.max(lon_for_min_max))
ax3_prob.set_xlim(np.min(lon_for_min_max),np.max(lon_for_min_max))
ax4_prob.set_xlim(np.min(lon_for_min_max),np.max(lon_for_min_max))
ax5_prob.set_xlim(np.min(lon_for_min_max),np.max(lon_for_min_max))
ax6_prob.set_xlim(np.min(lon_for_min_max),np.max(lon_for_min_max))
ax7_prob.set_xlim(np.min(lon_for_min_max),np.max(lon_for_min_max))

ax1_raw.set_xlim(np.min(lon_for_min_max),np.max(lon_for_min_max))
ax2_raw.set_xlim(np.min(lon_for_min_max),np.max(lon_for_min_max))
ax3_raw.set_xlim(np.min(lon_for_min_max),np.max(lon_for_min_max))
ax4_raw.set_xlim(np.min(lon_for_min_max),np.max(lon_for_min_max))
ax5_raw.set_xlim(np.min(lon_for_min_max),np.max(lon_for_min_max))
ax6_raw.set_xlim(np.min(lon_for_min_max),np.max(lon_for_min_max))
ax7_raw.set_xlim(np.min(lon_for_min_max),np.max(lon_for_min_max))

plt.show()


    