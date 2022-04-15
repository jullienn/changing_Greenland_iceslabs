# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 16:59:43 2021

@author: jullienn
"""

def exfunc(y,A,B,C):
    return A * np.exp(B * y) + C

def apply_normalisation(roll_corrected_array,mask, depth):
    #Where mask is False
    index_false=np.where(mask==False)[0]
    #Where mask is True
    index_true=np.where(mask==True)[0]
    
    #Prepare the mask for appliance
    mask.shape=mask.shape[0],1
    mask_matrix=np.repeat(mask,roll_corrected_array.shape[0],axis=1)
    mask_matrix=np.transpose(mask_matrix)
    #Apply the mask
    roll_corrected_array_masked=roll_corrected_array[mask_matrix]
    roll_corrected_array_masked=roll_corrected_array_masked.reshape((roll_corrected_array.shape[0],(roll_corrected_array.shape[1])-len(index_false)))
    
    #Create the depth expanded array
    depths_expanded=np.repeat(depth,roll_corrected_array_masked.shape[1])
    
    popt, pcov = scipy.optimize.curve_fit(exfunc, depths_expanded, roll_corrected_array_masked.flatten(),
                                                  bounds=((-np.inf, -np.inf, -np.inf),
                                                          ( np.inf,          0,          0)),
                                                  max_nfev=1000000)
    #Retrieve constants
    A,B,C = popt

    #Apply normalisation
    traces=roll_corrected_array_masked

    #pdb.set_trace()
    # Correct the traces and normalize them.
    # Original function is Z = A * e^(By) + C
    # Inverse function to normalize AND get rid of heteroscedasticitiy is 0 = ((Z - C)/A * e^(-By) - 1.0) * e^(By)
    traces_norm = ((traces - C) / A * np.exp(-B * depth) - 1.0) * np.exp(B * depth)
    # Then divide by the standard deviation of the traces to have them normalized for variance
    # All traces  for all tracks will have a MEAN of zero and a STDDEV of 1
    '''
    traces_norm = traces_norm / (np.std(traces_norm))
    '''
    #reconstruct the array with the NaNs
    traces_norm_full=np.zeros((roll_corrected_array.shape[0],roll_corrected_array.shape[1]))
    traces_norm_full[:,index_false]=np.nan
    traces_norm_full[:,index_true]=traces_norm

    return traces_norm_full

def identify_ice_lenses(traces,dry_firn_normalisation,depth,mask,datetrack,quantile_investigation,desired_quantiles):
    #traces_20m should be the 20m traces
    
    
    # We identified the minimum signal-cutoff and continuity-threshold values for each algorithm.  They
    # gave very close results.  Produce one image from each Algorithm and we will evaluate which one worked best on all the datasets.
    # These sets of cutoffs were produced from examination done in validate_reference_track_w_in_situ_data(),
    # and plot_validation_data_and_find_minima()
    
    '''
    #Original
    ALGORITHMS = ("orig","SG1","SG1")
    CUTOFFS = (-0.45, -0.45, -0.45)
    THRESHOLDS = (0, 0, 350)
    '''
    
    '''
    #Custom identification with distribution:
    # ---- -0.08927652699581005 is quantile(iceslabs,0.5)
    # ---- -0.04 is the hollow within the 2 distributions
    # ---- -0.00439656779575809 is quantile(iceslabs,0.75)
    ALGORITHMS = ("SG1","SG1","SG1")
    CUTOFFS = (-0.08927652699581005,-0.04,-0.00439656779575809)
    THRESHOLDS = (350, 350, 350)
    '''
    
    #pdb.set_trace()
    #Investigate custom threshold sensitivity
    ALGORITHMS = tuple(np.repeat("SG1",len(quantile_investigation)))
    CUTOFFS = tuple(quantile_investigation)
    THRESHOLDS = tuple(np.repeat(350,len(quantile_investigation)))
    
    #Initalize count to 0 for cutoff names
    count=0
    names_cutoff=desired_quantiles
    
    
    for algorithm, cutoff, continuity_threshold in zip(ALGORITHMS, CUTOFFS, THRESHOLDS):
        #pdb.set_trace()
        #Retrieve cutoff name
        cutoff_q=names_cutoff[count]
        
        # Apply the cutoff.
        boolean_traces = (traces <= cutoff)
        
        # Get rid of truth values in the top 2 pixels.  These are most often artifacts and should not be included
        boolean_traces[:2,:] = False

        # Apply the filter.
        if algorithm == "orig":
            pass
        elif algorithm == "SG1":
            boolean_traces = boolean_shrink_and_grow(boolean_traces, N=1)
        elif algorithm == "SG2":
            boolean_traces = boolean_shrink_and_grow(boolean_traces, N=2)
        elif algorithm == "S1":
            boolean_traces = boolean_grow_by_1(boolean_shrink_by_1(boolean_traces, N=1), N=1)
        elif algorithm == "S2":
            boolean_traces = boolean_grow_by_1(boolean_shrink_by_1(boolean_traces, N=2), N=2)

        # Perform the continuity thresholding.
        group_id_array, group_size_dict = caluculate_icelens_connectedness(boolean_traces)
        ice_lenses_above_cutoff_size = np.zeros(boolean_traces.shape, dtype=np.bool)

        # Set each group of pixels in that category to True, only for groups larger than the cutoff
        for group_ID in [ID for (ID,size) in list(group_size_dict.items()) if size >= continuity_threshold]:
            ice_lenses_above_cutoff_size[group_id_array == group_ID] = True
        
        #Inspect how looks like ice_lenses_above_cutoff_size
        '''
        if export:
            traces_refilled = self._refill_array(ice_lenses_above_cutoff_size, mask)

            fname_ext = "_{0}_CUTOFF_{1:0.2f}_THRESHOLD_{2:0>3d}.png".format(algorithm, cutoff, continuity_threshold)
            self.export_boolean_image(~traces_refilled, fname_ext)
        '''
        
        boolean_slabs=ice_lenses_above_cutoff_size
        
        #Reconstruct full slices with exclusions
        boolean_full_slabs=reconstruct_with_NaNs(dry_firn_normalisation,depth,mask,boolean_slabs)
        #Ok that works!
        
        #pdb.set_trace()
        #Save as pickle file     
        filename_tosave='C:/Users/jullienn/switchdrive/Private/research/RT1/masking_iceslabs/quantiles_threshold_application/'+datetrack+'_'+algorithm+'_cutoff_'+str(np.round(cutoff_q,2))+'_threshold_'+str(continuity_threshold)+'.pickle'
        outfile= open(filename_tosave, "wb" )
        pickle.dump(boolean_full_slabs,outfile)
        outfile.close()
        
        
        #Update count
        count=count+1
        
        
        
        '''
        #Save figures here
        
        fig1, (ax1) = plt.subplots(1, 1)
        fig1.suptitle(single_year)
        
        cb1=ax1.pcolor(np.arange(0,boolean_full_slabs.shape[1]),np.arange(0,boolean_full_slabs.shape[0]),boolean_full_slabs,cmap=plt.get_cmap('gray'))#,norm=divnorm)
        #fig1.colorbar(cb1,ax=ax1)
        ax1.invert_yaxis() #Invert the y axis = avoid using flipud.
        ax1.title.set_text('boolean full slabs')
        ax1.set_aspect(3)

        #Create the figure name
        fig_name=[]
        fig_name='C:/Users/jullienn/switchdrive/Private/research/RT1/dry_firn_normalisation/'+single_year+'_'+algorithm+'_cutoff_'+str(cutoff)+'_threshold_'+str(continuity_threshold)+'.png'
        
        #"_{0}_CUTOFF_{1:0.2f}_THRESHOLD_{2:0>3d}.png".format(algorithm, cutoff, continuity_threshold)
        
        #Save the figure
        plt.savefig(fig_name,dpi=2000)
        plt.clf()
        '''
    return 

'''       
    def export_boolean_image(array, image_label=""):
        #print('-------------------- ENTERING export_boolean_image --------------------')
'''        '''Create a black-and-white boolean image of this track.''' '''
        outfilename = self.NAME + ("_" if (len(image_label)>0 and image_label[0] != "_") else "") + image_label + (".png" if ((len(image_label) < 4) or (len(image_label) > 4 and image_label[-4:] != '.png')) else "")
        outfilepath = os.path.join(ICEBRIDGE_EXPORT_FOLDER, outfilename)

        png_file = png.from_array(array, mode="L;1")
        png_file.save(outfilepath)

        if self.VERBOSE:
            print("Exported", outfilename)
        #print('-------------------- OUT export_boolean_image --------------------')

        return        
'''        


def boolean_shrink_by_1(orig, N=1):
    #print('-------------------- ENTERING _boolean_shrink_by_1 --------------------')

    # The & operator will "shrink" the True values of the array.
    # If that pixel or any adjacent pixel (L,R,U,D) is not true, it will make that pixel not true.
    for _ in range(N):
        new = orig.copy()
        new[ :  , :-1] = new[ :  , :-1] & orig[ :  ,1:  ] # SUBSET_LEFT
        new[ :  ,1:  ] = new[ :  ,1:  ] & orig[ :  , :-1] # SUBSET_RIGHT
        new[ :-1, :  ] = new[ :-1, :  ] & orig[1:  , :  ] # SUBSET_UP
        new[1:  , :  ] = new[1:  , :  ] & orig[ :-1, :  ] # SUBSET_DOWN
        orig = new
    #print('-------------------- OUT _boolean_shrink_by_1 --------------------')

    return new

def boolean_grow_by_1(orig, N=1):
    #print('-------------------- ENTERING _boolean_grow_by_1 --------------------')

    # The | operator will "grow" the True values of the array.
    # If that pixel or any adjacent pixel (L,R,U,D) is true, it will make that pixel True.
    for _ in range(N):
        new = orig.copy()
        new[ :  , :-1] = new[ :  , :-1] | orig[ :  ,1:  ] # SUBSET_LEFT
        new[ :  ,1:  ] = new[ :  ,1:  ] | orig[ :  , :-1] # SUBSET_RIGHT
        new[ :-1, :  ] = new[ :-1, :  ] | orig[1:  , :  ] # SUBSET_UP
        new[1:  , :  ] = new[1:  , :  ] | orig[ :-1, :  ] # SUBSET_DOWN
        orig = new
    #print('-------------------- OUT _boolean_grow_by_1 --------------------')

    return new

def boolean_shrink_and_grow(boolean_array, N=1):
    #print('-------------------- ENTERING _boolean_shrink_and_grow --------------------')

    '''Take a boolean T/F array (assuming True == ice lenses) and use a
    "shrink and grow" method to get rid of noise.  Shrink, then grow the pixels
    by N steps to get rid of small errant strands and noise.'''
    array = boolean_array
    # SHRINK N TIMES
    for _ in range(N):
        array = boolean_shrink_by_1(array)
    for _ in range(N*2):
        array = boolean_grow_by_1(array)
    for _ in range(N):
        array = boolean_shrink_by_1(array)
    #print('-------------------- OUT _boolean_shrink_and_grow --------------------')

    return array

def caluculate_icelens_connectedness(boolean_image):

    '''A parent function that iterates over an image until all pixels have found which "group" they belong to.
    Return an int array of group_ID numbers (zero are empty pixels), and a dictionary of (ID:size) pairs.'''
    group_id_array = np.zeros(boolean_image.shape, dtype=np.int)
    visited_mask_empty = np.zeros(boolean_image.shape, dtype=np.bool)
    # Visited mask cumulative -- a boolean array of all the pixels we've visited.  Starts out empty, should match boolean_image in the end
    visited_mask_cumulative = visited_mask_empty.copy()
    # Keeps track of how many pixels are in each group.
    group_size_dict = {}

    # Keep iterating while there are still pixels we haven't visited.
    still_to_visit = boolean_image
    current_group_id = 1
    while np.count_nonzero(still_to_visit) > 0:
        # Grab the next non-zero pixel location in still_to_visit
        next_location = np.unravel_index(np.argmax(still_to_visit), boolean_image.shape)
        # Make a copy of the empty mask for this recursive trip.
        this_visited_mask = visited_mask_empty.copy()
        # Recurse!
        icelens_connectedness_iterator_subfunction(boolean_image, this_visited_mask, next_location)
        # Mark all pixels in the big map with this group's id
        group_id_array[this_visited_mask] = current_group_id
        # Save the size of this group to the dictionary
        this_group_size = np.count_nonzero(this_visited_mask)
        group_size_dict[current_group_id] = this_group_size
        # Add these pixels to the cumulative "visited" pixels
        visited_mask_cumulative = visited_mask_cumulative | this_visited_mask
        # Subtract them from the pixels still to visit.
        still_to_visit = (boolean_image & (~visited_mask_cumulative))
        # Add one to the current running group ID
        current_group_id += 1

    return group_id_array, group_size_dict
    
    
def icelens_connectedness_iterator_subfunction(boolean_image, visited_mask, pixel_coords):

    '''An iterative function for finding all connected pixels in a region of the image.'''
    # If THIS pixel is not an ice layer OR this pixel has already been visited, return zero.
    pixel_coords_list = [pixel_coords]
    visit_directions = [0] # 0=right, 1=down, 2=left, 3=up

    fright = lambda y,x:(y,x+1)
    fdown  = lambda y,x:(y+1,x)
    fleft  = lambda y,x:(y,x-1)
    fup    = lambda y,x:(y-1,x)
    direction_dict = {0:fright, 1:fdown, 2:fleft, 3:fup}

    in_bounds = lambda y,x: (x>=0 and y>=0 and x<boolean_image.shape[1] and y<boolean_image.shape[0])

    while len(pixel_coords_list) > 0:
        direction = visit_directions[-1]
        y,x = pixel_coords_list[-1]
        visited_mask[y,x] = True

        if (0 <= direction <= 3):
            next_loc = direction_dict[direction](y,x)
            # We'll look the next direction at THIS pixel
            visit_directions[-1] += 1
            # If it's in bounds, still in the ice lens and not yet visited, go there next.
            if in_bounds(*next_loc) and boolean_image[next_loc] and (not visited_mask[next_loc]):
                pixel_coords_list.append(next_loc)
                visit_directions.append(0)
        elif (direction == 4):
            # Done!  Pop back to the last pixel.
            pixel_coords_list.pop()
            visit_directions.pop()
        else: # Shouldn't get here.
            assert False

    return visited_mask

def select_20m_slice_without_NaNs(dry_firn_normalisation,depth,mask):
    #I should input slices without the NaNs in it 
    #Where mask is False
    index_false=np.where(mask==False)[0]
    #Where mask is True
    index_true=np.where(mask==True)[0]
    
    #Prepare the mask for appliance
    mask.shape=mask.shape[0],1
    mask_matrix=np.repeat(mask,dry_firn_normalisation.shape[0],axis=1)
    mask_matrix=np.transpose(mask_matrix)
    #Apply the mask
    dry_firn_normalisation_masked=dry_firn_normalisation[mask_matrix]
    dry_firn_normalisation_masked=dry_firn_normalisation_masked.reshape((dry_firn_normalisation.shape[0],(dry_firn_normalisation.shape[1])-len(index_false)))
    
    #Keep only the first 20m of slices
    boolean_20m=depth <= 20
    ind_20m=np.where(boolean_20m==True)[0]
    
    traces=dry_firn_normalisation_masked
    traces_20m=traces[ind_20m,:]
    
    return traces_20m


def reconstruct_with_NaNs(dry_firn_normalisation,depth,mask,boolean_iceslabs):
    #pdb.set_trace()
    #Where mask is False
    index_false=np.where(mask==False)[0]
    #Where mask is True
    index_true=np.where(mask==True)[0]
    
    #reconstruct the array with the NaNs
    traces_20m_full=np.zeros((boolean_iceslabs.shape[0],dry_firn_normalisation.shape[1]))
    traces_20m_full[:,index_false]=False#np.nan
    traces_20m_full[:,index_true]=boolean_iceslabs
    
    return traces_20m_full

def extract_surface_return(slice_roll_corrected):
    
    #slice_roll_corrected is dataframe['year_of_interest']['roll_corrected']
    
    # --- Remove the average
    #Let's say we take the 1 top pixels
    surface_return=slice_roll_corrected[0,]
    #substract the average of surface_return to the whole radar slice
    roll_corrected_after_surf_removal=slice_roll_corrected-np.nanmean(surface_return)

    '''
    # --- Remove the top at each column
    roll_corrected_after_surf_removal=np.empty((slice_roll_corrected.shape[0],slice_roll_corrected.shape[1]))
    
    #Set half the size of the moving window
    size_mov_window=5
    
    for i in range(0,slice_roll_corrected.shape[1]):
        #pdb.set_trace()
        
        if (i==0):
            mov_window=slice_roll_corrected[0,0]
        elif (i==(slice_roll_corrected.shape[1]-1)):
            mov_window=slice_roll_corrected[0,-1]
        elif (i<(size_mov_window)):
            mov_window=slice_roll_corrected[0,0:2*i+1]
        elif (i>=(slice_roll_corrected.shape[1]-size_mov_window)):
            mov_window=slice_roll_corrected[0,(i-(slice_roll_corrected.shape[1]-i-1)):slice_roll_corrected.shape[1]]
        else:
            mov_window=slice_roll_corrected[0,(i-size_mov_window):(i+size_mov_window+1)]
        
        #Remove the moving window to the data
        roll_corrected_after_surf_removal[:,i]=slice_roll_corrected[:,i]-np.nanmean(mov_window)
    '''   
    return roll_corrected_after_surf_removal

import pickle
import matplotlib.pyplot as plt
import numpy as np
import scipy.io
import h5py
import scipy.optimize
import pdb
from PIL import Image
from sklearn.metrics.cluster import contingency_matrix
import pandas as pd

create_pickle='TRUE'
display_pickle='FALSE'
display_plots_quick_check='FALSE'

investigation_quantile='TRUE'
gaussian_calibration='TRUE'
#For quantile investigation plotting
show_reference_trace='TRUE'
show_case_study='TRUE'

#1. Open roll corrected of the specific year
'''
investigation_year={2010:'empty',
                    2011:'empty',
                    2012:['Data_20120423_01_137.mat','Data_20120423_01_138.mat'],
                    2013:['Data_20130409_01_010.mat','Data_20130409_01_011.mat','Data_20130409_01_012.mat'],
                    2014:'empty',
                    2017:'empty',
                    2018:['Data_20180421_01_004.mat','Data_20180421_01_005.mat','Data_20180421_01_006.mat','Data_20180421_01_007.mat']}

'''
'''
investigation_year={2010:['Data_20100513_01_001.mat','Data_20100513_01_002.mat'],
                    2011:['Data_20110411_01_116.mat','Data_20110411_01_117.mat','Data_20110411_01_118.mat'],
                    2012:['Data_20120428_01_125.mat','Data_20120428_01_126.mat'],
                    2013:'empty',
                    2014:['Data_20140408_11_024.mat','Data_20140408_11_025.mat','Data_20140408_11_026.mat'],
                    2017:['Data_20170508_02_165.mat','Data_20170508_02_166.mat','Data_20170508_02_167.mat','Data_20170508_02_168.mat','Data_20170508_02_169.mat','Data_20170508_02_170.mat','Data_20170508_02_171.mat'],
                    2018:'empty'}
'''
'''
#7years case study
investigation_year={2010:['Data_20100508_01_114.mat','Data_20100508_01_115.mat'],
                    2011:['Data_20110419_01_008.mat','Data_20110419_01_009.mat','Data_20110419_01_010.mat'],
                    2012:['Data_20120418_01_129.mat','Data_20120418_01_130.mat','Data_20120418_01_131.mat'],
                    2013:['Data_20130405_01_165.mat','Data_20130405_01_166.mat','Data_20130405_01_167.mat'],
                    2014:['Data_20140424_01_002.mat','Data_20140424_01_003.mat','Data_20140424_01_004.mat'],
                    2017:['Data_20170422_01_168.mat','Data_20170422_01_169.mat','Data_20170422_01_170.mat','Data_20170422_01_171.mat'],
                    2018:['Data_20180427_01_170.mat','Data_20180427_01_171.mat','Data_20180427_01_172.mat']}
'''
'''
#Calibration track in MacFerrin et al, 2019
investigation_year={2010:'empty',
                    2011:'empty',
                    2012:['Data_20120423_01_137.mat','Data_20120423_01_138.mat'],
                    2013:['Data_20130409_01_010.mat','Data_20130409_01_011.mat','Data_20130409_01_012.mat'],
                    2014:'empty',
                    2017:'empty',
                    2018:['Data_20180421_01_004.mat','Data_20180421_01_005.mat','Data_20180421_01_006.mat','Data_20180421_01_007.mat']}
#2014 and 2017 almost colocated
'''

#Calibration track in MacFerrin et al, 2019
investigation_year={2010:'empty',
                    2011:'empty',
                    2012:'empty',
                    2013:['Data_20130409_01_010.mat','Data_20130409_01_011.mat','Data_20130409_01_012.mat'],
                    2014:'empty',
                    2017:'empty',
                    2018:'empty'}
#2014 and 2017 almost colocated


#Investigation failing ice slabs likelihood
#list_trace=list(['20110416_01_053_055':['Data_20110416_01_053.mat','Data_20110416_01_054.mat','Data_20110416_01_055.mat']
                #'20120421_01_052_052':['Data_20120421_01_052.mat']
                #'20130423_01_125_125':['Data_20130423_01_125.mat']
                #'20130423_01_127_127':['Data_20130423_01_127.mat']
                #'20130426_01_089_089':['Data_20130426_01_089.mat']
                #'20140419_01_016_017':['Data_20140419_01_016.mat','Data_20140419_01_017.mat']
                #'20140419_01_028_028':['Data_20140419_01_028.mat']
                #'20140419_03_075_075':['Data_20140419_03_075.mat']
                #'20140516_02_031_034':['Data_20140516_02_031.mat','Data_20140516_02_032.mat','Data_20140516_02_033.mat','Data_20140516_02_034.mat']
                #'20180419_02_032_033':['Data_20180419_02_032.mat','Data_20180419_02_033.mat']
                #'20180419_02_035_036':['Data_20180419_02_035.mat','Data_20180419_02_036.mat']
                #'20180425_01_166_169':['Data_20180425_01_166.mat','Data_20180425_01_167.mat','Data_20180425_01_168.mat','Data_20180425_01_169.mat']
                #'20180427_01_170_172':['Data_20180427_01_170.mat','Data_20180427_01_171.mat','Data_20180427_01_172.mat']
                #'20180429_01_008_014'])['Data_20180429_01_008.mat','Data_20180429_01_009.mat','Data_20180429_01_010.mat','Data_20180429_01_011.mat','Data_20180429_01_012.mat','Data_20180429_01_013.mat','Data_20180429_01_014.mat']
#One date in 2018 where processing is good for comparison: '20180427_01_004_006':['Data_20180427_01_004.mat','Data_20180427_01_005.mat','Data_20180427_01_006.mat']
'''
investigation_year={2010:'empty',
                    2011:'empty',
                    2012:'empty',
                    2013:'empty',
                    2014:'empty',
                    2017:'empty',
                    2018:['Data_20180427_01_004.mat','Data_20180427_01_005.mat','Data_20180427_01_006.mat']}
'''

if (create_pickle == 'TRUE'):
    ##############################################################################
    ###                             Load data                                  ###
    ##############################################################################
    #Compute the speed (Modified Robin speed):
    # self.C / (1.0 + (coefficient*density_kg_m3/1000.0))
    v= 299792458 / (1.0 + (0.734*0.873/1000.0))
    
    path_data='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/'
    path_roll_corrected='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2010_2018/i_out_from_IceBridgeGPR_Manager_v2.py/pickles_and_images/Roll_Corrected_Picklefiles/'
    path_mask='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2010_2018/i_out_from_IceBridgeGPR_Manager_v2.py/pickles_and_images/Boolean_Array_Picklefiles/'
    dataframe={}

    #Define the desired quantiles
    desired_quantiles=np.arange(0.63,0.82,0.01)
    filename_quantiles='C:/Users/jullienn/switchdrive/Private/research/RT1/masking_iceslabs/quantiles_threshold_application/quantile_file_'+str(np.round(desired_quantiles[0],2))+'_'+str(np.round(desired_quantiles[-1],2))+'.txt'

    for single_year in investigation_year.keys():
        print(single_year)
        
        #If no data, continue
        if (investigation_year[single_year]=='empty'):
            print('No data for year '+str(single_year)+', continue')
            continue
        
        ###1. Load the mask and depth_corrected files:
        start_date_track=investigation_year[single_year][0]
        end_date_track=investigation_year[single_year][-1]
        date_track=start_date_track[5:20]+'_'+end_date_track[17:20]
        
        #Define filename roll, depth corrected data and mask
        filename_roll_corrected=date_track+'_ROLL_CORRECTED.pickle'
        filename_mask=date_track+'_mask.pickle'
                     
        #Open the roll corrected file
        f_roll_corrected = open(path_roll_corrected+filename_roll_corrected, "rb")
        roll_corrected = pickle.load(f_roll_corrected)
        f_roll_corrected.close()
        
        #Open the mask file
        f_mask = open(path_mask+filename_mask, "rb")
        mask = pickle.load(f_mask)
        f_mask.close()
        
        ###2. Load the latitude and longitude
        lat_appended=[]
        lon_appended=[]
            
        for indiv_file_load in investigation_year[single_year]:
            print(indiv_file_load)
            
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
            
            lat_appended=np.flipud(lat_appended)
            lon_appended=np.flipud(lon_appended)
            roll_corrected=np.fliplr(roll_corrected)
            mask=np.flipud(mask)
            
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
        
        #########################################################################
        # From plot_2002_2003.py - END
        #########################################################################
        
        #Store reunited lat/lon, slice output and mask in a dictionnary:
        dataframe[str(single_year)]={'lat_appended':lat_appended,
                                     'lon_appended':lon_appended,
                                     'depth':depth,
                                     'roll_corrected':roll_corrected,
                                     'roll_corrected_after_surf_removal':np.nan,
                                     'depth_corrected_after_surf_removal_without_norm':np.nan,
                                     'mask':mask,
                                     'datetrack':date_track}
        
        ##############################################################################
        ###                          Load and organise data                        ###
        ##############################################################################
    pdb.set_trace()
    #3. Extract surface return and perform depth correction without normalisation
    for single_year in investigation_year.keys():
        print('--- Perform depth correction ---')
        #If no data, continue
        if (investigation_year[single_year]=='empty'):
            print('No data for year '+str(single_year)+', continue')
            continue
        
        print(single_year)
        #Extract surface return
        
        dataframe[str(single_year)]['roll_corrected_after_surf_removal']=extract_surface_return(dataframe[str(single_year)]['roll_corrected'])
        #Perform depth correction + normalisation
        #if 2010 or 2011, depth is from 0:428
        if (str(single_year) in list(['2010','2011'])):
            dataframe[str(single_year)]['depth_corrected_after_surf_removal_without_norm']=apply_normalisation(dataframe[str(single_year)]['roll_corrected_after_surf_removal'],dataframe[str(single_year)]['mask'],dataframe[str(single_year)]['depth'][0:428])
        else:
            dataframe[str(single_year)]['depth_corrected_after_surf_removal_without_norm']=apply_normalisation(dataframe[str(single_year)]['roll_corrected_after_surf_removal'],dataframe[str(single_year)]['mask'],dataframe[str(single_year)]['depth'][0:201])
        
        #Save as the depth corrected trace as pickle file     
        filename_tosave='C:/Users/jullienn/switchdrive/Private/research/RT1/masking_iceslabs/quantiles_threshold_application/'+dataframe[str(single_year)]['datetrack']+'_Depth_Corrected_surf_removal.pickle'
        outfile= open(filename_tosave, "wb" )
        pickle.dump(dataframe[str(single_year)]['depth_corrected_after_surf_removal_without_norm'],outfile)
        outfile.close()
        
        #Save depth corrected, depth and lat/lon
        filename_tosave='C:/Users/jullienn/switchdrive/Private/research/RT1/masking_iceslabs/quantiles_threshold_application/'+dataframe[str(single_year)]['datetrack']+'_dataframeforS2.pickle'
        outfile= open(filename_tosave, "wb" )
        pickle.dump(dataframe[str(single_year)],outfile)
        outfile.close()
        
        print('Exporting '+dataframe[str(single_year)]['datetrack']+' depth corrected pickle file')
        
    if (display_plots_quick_check=='TRUE'):
        
        #Display results
        fig1, (ax1,ax2,ax3,ax4,ax5,ax6,ax7) = plt.subplots(7, 1)
        fig1.suptitle('Roll corrected')
    
        cb1=ax1.imshow(dataframe['2010']['roll_corrected'],cmap=plt.get_cmap('gray'))#,norm=divnorm)
        fig1.colorbar(cb1,ax=ax1)
        ax1.set_ylim(129,0)
        cb1.set_clim(-12.0,0)
    
        cb2=ax2.imshow(dataframe['2011']['roll_corrected'],cmap=plt.get_cmap('gray'))#,norm=divnorm)
        fig1.colorbar(cb2,ax=ax2)
        ax2.set_ylim(129,0)
        cb2.set_clim(-12.0,0)
        
        cb3=ax3.imshow(dataframe['2012']['roll_corrected'],cmap=plt.get_cmap('gray'))#,norm=divnorm)
        fig1.colorbar(cb3,ax=ax3)
        ax3.set_ylim(61,0)
        cb3.set_clim(-12.0,0)
       
        cb4=ax4.imshow(dataframe['2013']['roll_corrected'],cmap=plt.get_cmap('gray'))#,norm=divnorm)
        fig1.colorbar(cb4,ax=ax4)
        ax4.set_ylim(61,0)
        cb4.set_clim(-12.0,0)
        
        cb5=ax5.imshow(dataframe['2014']['roll_corrected'],cmap=plt.get_cmap('gray'))#,norm=divnorm)
        fig1.colorbar(cb5,ax=ax5)
        ax5.set_ylim(61,0)
        cb5.set_clim(-12.0,0)
        
        cb6=ax6.imshow(dataframe['2017']['roll_corrected'],cmap=plt.get_cmap('gray'))#,norm=divnorm)
        fig1.colorbar(cb6,ax=ax6)
        ax6.set_ylim(61,0)
        cb6.set_clim(-12.0,0)
        
        cb7=ax7.imshow(dataframe['2018']['roll_corrected'],cmap=plt.get_cmap('gray'))#,norm=divnorm)
        fig1.colorbar(cb7,ax=ax7)
        ax7.set_ylim(61,0)
        cb7.set_clim(-12.0,0)
        
        plt.show()
        
        fig2, (ax1,ax2,ax3,ax4,ax5,ax6,ax7) = plt.subplots(7, 1)
        fig2.suptitle('Roll corrected after removal of surface return')
    
        cb1=ax1.imshow(dataframe['2010']['roll_corrected_after_surf_removal'],cmap=plt.get_cmap('gray'))#,norm=divnorm)
        fig2.colorbar(cb1,ax=ax1)
        ax1.set_ylim(129,0)
        cb1.set_clim(-2.0,2)
      
        cb2=ax2.imshow(dataframe['2011']['roll_corrected_after_surf_removal'],cmap=plt.get_cmap('gray'))#,norm=divnorm)
        fig2.colorbar(cb2,ax=ax2)
        ax2.set_ylim(129,0)
        cb2.set_clim(-2.0,2)
        
        cb3=ax3.imshow(dataframe['2012']['roll_corrected_after_surf_removal'],cmap=plt.get_cmap('gray'))#,norm=divnorm)
        fig2.colorbar(cb3,ax=ax3)
        ax3.set_ylim(61,0)
        cb3.set_clim(-2.0,2)
        
        cb4=ax4.imshow(dataframe['2013']['roll_corrected_after_surf_removal'],cmap=plt.get_cmap('gray'))#,norm=divnorm)
        fig2.colorbar(cb4,ax=ax4)
        ax4.set_ylim(61,0)
        cb4.set_clim(-2.0,2)
       
        cb5=ax5.imshow(dataframe['2014']['roll_corrected_after_surf_removal'],cmap=plt.get_cmap('gray'))#,norm=divnorm)
        fig2.colorbar(cb5,ax=ax5)
        ax5.set_ylim(61,0)
        cb5.set_clim(-2.0,2)
        
        cb6=ax6.imshow(dataframe['2017']['roll_corrected_after_surf_removal'],cmap=plt.get_cmap('gray'))#,norm=divnorm)
        fig2.colorbar(cb6,ax=ax6)
        ax6.set_ylim(61,0)
        cb6.set_clim(-2.0,2)
        
        cb7=ax7.imshow(dataframe['2018']['roll_corrected_after_surf_removal'],cmap=plt.get_cmap('gray'))#,norm=divnorm)
        fig2.colorbar(cb7,ax=ax7)
        ax7.set_ylim(61,0)
        cb7.set_clim(-2.0,2)
        
        plt.show()
        
        pdb.set_trace()
        
        #Display roll corrected, removal of surface return, depth corrected without normalisation
        fig1, (ax3,ax4,ax7) = plt.subplots(3, 1)
        fig1.suptitle('Depth corrected without normalisation, surface return removal')
        '''
        cb1=ax1.imshow(dataframe['2010']['depth_corrected_after_surf_removal_without_norm'],cmap=plt.get_cmap('gray'))#,norm=divnorm)
        fig1.colorbar(cb1,ax=ax1)
        ax1.set_ylim(129,0)
        #cb1.set_clim(-2.0,0)
    
        cb2=ax2.imshow(dataframe['2011']['depth_corrected_after_surf_removal_without_norm'],cmap=plt.get_cmap('gray'))#,norm=divnorm)
        fig1.colorbar(cb2,ax=ax2)
        ax2.set_ylim(129,0)
        #cb2.set_clim(-12.0,0)
        '''
        cb3=ax3.imshow(dataframe['2012']['depth_corrected_after_surf_removal_without_norm'],cmap=plt.get_cmap('gray'))#,norm=divnorm)
        fig1.colorbar(cb3,ax=ax3)
        ax3.set_ylim(61,0)
        #cb3.set_clim(-12.0,0)
       
        cb4=ax4.imshow(dataframe['2013']['depth_corrected_after_surf_removal_without_norm'],cmap=plt.get_cmap('gray'))#,norm=divnorm)
        fig1.colorbar(cb4,ax=ax4)
        ax4.set_ylim(61,0)
        #cb4.set_clim(-12.0,0)
        '''
        cb5=ax5.imshow(dataframe['2014']['depth_corrected_after_surf_removal_without_norm'],cmap=plt.get_cmap('gray'))#,norm=divnorm)
        fig1.colorbar(cb5,ax=ax5)
        ax5.set_ylim(61,0)
        #cb5.set_clim(-12.0,0)
        
        cb6=ax6.imshow(dataframe['2017']['depth_corrected_after_surf_removal_without_norm'],cmap=plt.get_cmap('gray'))#,norm=divnorm)
        fig1.colorbar(cb6,ax=ax6)
        ax6.set_ylim(61,0)
        #cb6.set_clim(-12.0,0)
        '''
        cb7=ax7.imshow(dataframe['2018']['depth_corrected_after_surf_removal_without_norm'],cmap=plt.get_cmap('gray'))#,norm=divnorm)
        fig1.colorbar(cb7,ax=ax7)
        ax7.set_ylim(61,0)
        #cb7.set_clim(-12.0,0)
        plt.show()
    
    if (gaussian_calibration=='TRUE'):
        ### -------- Gaussian of iceslabs VS non iceslabs
        
        #Open the conservative mask
        path_mask='C:/Users/jullienn/switchdrive/Private/research/RT1/masking_iceslabs/'
        boolean_mask = Image.open(path_mask+'binary_conservative_mask_20130409_01_010_012_XDEPTHCORRECT_AFTER.png').convert("L")
        arr_boolean_mask = np.asarray(boolean_mask)
        
        #Keep only the first 20m of the 30m boolean mask
        boolean_20m=dataframe['2013']['depth'] <= 20
        ind_20m=np.where(boolean_20m==True)[0]
        arr_boolean_mask_20m=arr_boolean_mask[ind_20m,:]
        
        #Display the mask on 2013 data
        fig1, (ax1) = plt.subplots(1, 1)
        fig1.suptitle('Depth corrected with 2013 conservative ask on top of it - 30m slice')
        cb1=ax1.imshow(dataframe['2013']['depth_corrected_after_surf_removal_without_norm'],cmap=plt.get_cmap('gray'))#,norm=divnorm)
        fig1.colorbar(cb1,ax=ax1)
        ax1.set_ylim(61,0)
        ax1.imshow(arr_boolean_mask,cmap=plt.get_cmap('autumn'), alpha=0.1)#,norm=divnorm)
        plt.show()
                
        #Apply mask to depth corrected traces
        #-- a. extract the 20m of depth corrected traces
        depth_corr_20m=dataframe['2013']['depth_corrected_after_surf_removal_without_norm'][ind_20m,:]
        #-- b. convert the image ask into a boolean
        boolean_to_apply=~(arr_boolean_mask_20m==255)
        #-- c. apply it to ice slabs
        iceslabs=depth_corr_20m[boolean_to_apply]
        #-- d. apply it to the remaining, i.e. dry firn
        dry_firn=depth_corr_20m[~boolean_to_apply]
        
        #Display the distribution of the signals
        fig1, (ax1) = plt.subplots(1, 1)
        #fig1.suptitle('Signal strenght distribution - Calibration trace in MacFerrin et al., 2019')
        ax1.hist(iceslabs,bins=500,density=True,label='Ice')
        ax1.hist(dry_firn,bins=500,density=True,alpha=0.2,label='Dry firn')
        ax1.legend()
        ax1.set_xlabel('Radar signal strength [dB]')
        ax1.set_ylabel('Probability density [ ]')
        
        #Define quantiles for investigation of accuracy
        quantile_investigation=np.quantile(iceslabs,desired_quantiles)
        
        #Add low and high quantile as dashed lines
        ax1.axvline(x=quantile_investigation[0],linestyle='--',color='k')
        ax1.axvline(x=quantile_investigation[-1],linestyle='--',color='k')

        plt.show()
        
        #Save for supplementary figure plotting
        #Save as pickle file     
        filename_tosave='C:/Users/jullienn/switchdrive/Private/research/RT1/masking_iceslabs/quantiles_threshold_application/referece_iceslabs_distrib.pickle'
        outfile= open(filename_tosave, "wb" )
        pickle.dump(iceslabs,outfile)
        outfile.close()
        
        filename_tosave='C:/Users/jullienn/switchdrive/Private/research/RT1/masking_iceslabs/quantiles_threshold_application/referece_dry_firn_distrib.pickle'
        outfile= open(filename_tosave, "wb" )
        pickle.dump(dry_firn,outfile)
        outfile.close()
                
        '''
        f_quantiles = open(filename_quantiles, "w")
        f_quantiles.write(str(np.round(desired_quantiles,2))+'\n')
        f_quantiles.write(str(quantile_investigation))
        f_quantiles.close() #Close the quantile file when we’re done!
        '''
        '''
        plt.savefig('C:/Users/jullienn/switchdrive/Private/research/RT1/figures/S2/figS2.png',dpi=300)
        '''
    else:
        #pdb.set_trace()
        #Then open the quantile file created previously        
        quantile_file = pd.read_csv(filename_quantiles, sep=" ", header=None)
        quantile_investigation=np.asarray(quantile_file.iloc[1])
    
    '''
    Not usefull anymore but kept in case needed
    boolean_mask_nb=np.empty((depth_corr_20m.shape[0],depth_corr_20m.shape[1]))
    boolean_mask_nb[:] = np.nan
    
    where_slabs=(arr_boolean_mask_20m==0)
    boolean_mask_nb[where_slabs]=depth_corr_20m[where_slabs]
    
    #Get rid of NaNs for plotting the boxplot
    boolean_mask_nb_flatten=np.ndarray.flatten(boolean_mask_nb)
    nan_index=np.isnan(boolean_mask_nb_flatten)
    boolean_mask_nb_flatten_without_nans=boolean_mask_nb_flatten[~nan_index]
    '''
    #pdb.set_trace()
    #5. apply smoothing and thresholding
    for single_year in investigation_year.keys():
        print('--- Creating the pickle files ---')
        #If no data, continue
        if (investigation_year[single_year]=='empty'):
            print('No data for year '+str(single_year)+', continue')
            continue
        
        print(single_year)
        #pdb.set_trace()
        
        iceslabs_to_extract=dataframe[str(single_year)]['depth_corrected_after_surf_removal_without_norm']
        depth=dataframe[str(single_year)]['depth']
        mask=dataframe[str(single_year)]['mask']
        datetrack=dataframe[str(single_year)]['datetrack']
        
        #Extrac the 20m slices and get rid of exclusions
        traces_20m=select_20m_slice_without_NaNs(iceslabs_to_extract,depth,mask)
        #Identify ice slabs
        boolean_slabs=identify_ice_lenses(traces_20m,iceslabs_to_extract,depth,mask,datetrack,quantile_investigation,desired_quantiles)
        #Ok all of the above works!!!! Great :). Just need to save the figures of display them and that's it
        
    print('end')

pdb.set_trace()

if (investigation_quantile=='TRUE'):
    
    #pdb.set_trace()
    
    #Compute the speed (Modified Robin speed):
    # self.C / (1.0 + (coefficient*density_kg_m3/1000.0))
    v= 299792458 / (1.0 + (0.734*0.873/1000.0))
    
    #Define paths
    path_data='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/'
    path_boolean_remove_surf='C:/Users/jullienn/switchdrive/Private/research/RT1/masking_iceslabs/quantiles_threshold_application/'
    path_depth_corrected='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2010_2018/pickles_and_images/Depth_Corrected_Picklefiles/'
    
    if (show_reference_trace=='TRUE'):
        
        #Define the year
        single_year=2013
        
        #Construct the date for data loading
        start_date_track=investigation_year[single_year][0]
        end_date_track=investigation_year[single_year][-1]
        date_track=start_date_track[5:20]+'_'+end_date_track[17:20]
        
        #Define filename roll, depth corrected data and mask
        filename_depth_corrected=date_track+'_DEPTH_CORRECTED_surf_removal.pickle'
                     
        #Open the roll corrected file
        f_depth_corrected = open(path_boolean_remove_surf+filename_depth_corrected, "rb")
        depth_corrected = pickle.load(f_depth_corrected)
        f_depth_corrected.close()
        
        #Define the quantiles to open
        quantiles_open=np.round(np.arange(0.63,0.82,0.01),2)
        
        #pdb.set_trace()
        #Set dataframe
        dataframe={}
        dataframe={k: {} for k in list(quantiles_open)}
        
        #Open the quantile files
        for indiv_quantile in dataframe.keys():
            #Define filename of boolean files 
            filename_boolean_quantile=date_track+'_SG1_cutoff_'+str(indiv_quantile)+'_threshold_350.pickle'
            
            #Load boolean files
            f_boolean_quantile = open(path_boolean_remove_surf+filename_boolean_quantile, "rb")
            dataframe[indiv_quantile] = pickle.load(f_boolean_quantile)
            f_boolean_quantile.close()
        
        #########################################################################
        # From plot_2002_2003.py - END
        #########################################################################
        
        ###2. Load the latitude and longitude
        lat_appended=[]
        lon_appended=[]
            
        for indiv_file_load in investigation_year[single_year]:
            print(indiv_file_load)
            
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
            
            lat_appended=np.flipud(lat_appended)
            lon_appended=np.flipud(lon_appended)
            
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
            
        pdb.set_trace()
        #Store reunited lat/lon, slice output and mask in a dictionnary:
        dataframe['lat_appended']=lat_appended
        dataframe['lon_appended']=lon_appended
        dataframe['depth']=depth
        dataframe['datetrack']=date_track
        dataframe['depth_corrected']=depth_corrected
        
        #Open the conservative mask
        path_mask='C:/Users/jullienn/switchdrive/Private/research/RT1/masking_iceslabs/'
        boolean_mask = Image.open(path_mask+'binary_conservative_mask_20130409_01_010_012_XDEPTHCORRECT_AFTER.png').convert("L")
        arr_boolean_mask = np.asarray(boolean_mask)
        
        #Keep only the first 20m of the 30m boolean mask
        boolean_20m=dataframe['depth'] <= 20
        ind_20m=np.where(boolean_20m==True)[0]
        arr_boolean_mask_20m=arr_boolean_mask[ind_20m,:]
        # Convert the image mask into a boolean
        arr_boolean_mask_20m=~arr_boolean_mask_20m
        arr_boolean_mask_20m[arr_boolean_mask_20m==255]=1
        
        dataframe['mask_truth']=arr_boolean_mask_20m
        ##############################################################################
        ###                          Load and organise data                        ###
        ##############################################################################
        
        #pdb.set_trace()
        
        #Extract overall accuracy and plot quatile VS accuracy
    
        #Contigency table to perform here
        
        #Define the Overall accuracy vector
        OA=np.zeros(len(quantiles_open))
        PA_dryfirn=np.zeros(len(quantiles_open))
        UA_dryfirn=np.zeros(len(quantiles_open))
        PA_iceslabs=np.zeros(len(quantiles_open))
        UA_iceslabs=np.zeros(len(quantiles_open))
        
        #Use sklearn.metrics.cluster.contingency_matrix, cf. https://scikit-learn.org/stable/modules/generated/sklearn.metrics.cluster.contingency_matrix.html
        #loop over the different quantiles
        for i in range(0,len(quantiles_open)):
            print(quantiles_open[i])
            
            cont_matrix=contingency_matrix(dataframe['mask_truth'], dataframe[quantiles_open[i]])
            
            if (cont_matrix.shape[1]==1):
                #Only one class, store nan
                #0 is dry firn, 1 is slabs
                OA[i]=np.nan
                #error of omissions = along the columns
                PA_dryfirn[i]=np.nan
                UA_dryfirn[i]=np.nan
                
                PA_iceslabs[i]=np.nan
                UA_iceslabs[i]=np.nan
            else:
                #0 is dry firn, 1 is slabs
                OA[i]=(cont_matrix[0,0]+cont_matrix[1,1])/np.sum(cont_matrix)
                #error of omissions = along the columns
                PA_dryfirn[i]=cont_matrix[0,0]/(np.sum(cont_matrix[0,:]))
                UA_dryfirn[i]=cont_matrix[0,0]/(np.sum(cont_matrix[:,0]))
                
                PA_iceslabs[i]=cont_matrix[1,1]/(np.sum(cont_matrix[1,:]))
                UA_iceslabs[i]=cont_matrix[1,1]/(np.sum(cont_matrix[:,1]))
            
        
        fig, (ax1,ax2,ax3) = plt.subplots(1, 3)
        fig.suptitle('Accuracy VS quantile')
        
        #Plot OA
        ax1.plot(quantiles_open,OA)
        ax1.grid()
        ax1.set_title("Overall accuracy")
        ax1.set_xlabel('Quantile [ ]')
        ax1.set_ylabel('Accuracy [%]')
        
        #Plot PA
        #Producer's Accuracy is the map accuracy from the point of view of the map
        #maker (the producer). This is how often are real features on the ground
        #correctly shown on the classified map or the probability that a certain
        #land cover of an area on the ground is classified as such
        #(from http://gsp.humboldt.edu/olm_2019/courses/GSP_216_Online/lesson6-2/metrics.html)
        # ---> This is this quantity we are interested in !!
        ax2.plot(quantiles_open,PA_dryfirn,label='dry firn')
        ax2.plot(quantiles_open,PA_iceslabs,label='ice slabs')
        ax2.set_title("Producer's accuracy")
        ax2.legend()
        ax2.set_xlabel('Quantile [ ]')
        ax2.set_ylabel('Accuracy [%]')
        ax2.grid()
        
        #Plot UA
        #the User's accuracy essentially tells use how often the class on the map
        #will actually be present on the ground. This is referred to as reliability
        #(from http://gsp.humboldt.edu/olm_2019/courses/GSP_216_Online/lesson6-2/metrics.html)
        
        ax3.plot(quantiles_open,UA_dryfirn,label='dry firn')
        ax3.plot(quantiles_open,UA_iceslabs,label='iceslabs')
        ax3.set_title("User's accuracy")
        ax3.legend()
        ax3.set_xlabel('Quantile [ ]')
        ax3.set_ylabel('Accuracy [%]')
        ax3.grid()
        
        plt.show()
        pdb.set_trace()
        
        #Display the resulting slabs identification
        #path_savefig='C:/Users/jullienn/switchdrive/Private/research/RT1/masking_iceslabs/quantile_investigation/'
        
        for i in range(0,len(quantiles_open)):
            
            #Define fig name
            fig_name=path_savefig+'referencetrace_quant_'+str(quantiles_open[i])+'.png'
            
            #Prepare the plot
            fig, (ax1) = plt.subplots(1, 1)
            fig.suptitle('Custom threshold: quantile'+str(quantiles_open[i])+' of ice slabs distribution, SG1, 350 continuity from 2013 trace in MF2019')
            figManager = plt.get_current_fig_manager()
            figManager.window.showMaximized()
            
            #Replace where dry firn by nan so that overlay plot can be possible
            quantile_to_plot=dataframe[quantiles_open[i]]
            quantile_to_plot[quantile_to_plot==0]=np.nan
            
            #Plot custom threshold ice slabs identification
            ax1.imshow(dataframe['depth_corrected'],cmap=plt.get_cmap('gray'))#,norm=divnorm)
            ax1.imshow(quantile_to_plot,cmap=plt.get_cmap('gray'))#,norm=divnorm)
            ax1.imshow(dataframe['mask_truth'],cmap=plt.get_cmap('OrRd'), alpha=0.2)
            ax1.title.set_text(dataframe['datetrack']+' - quantile: '+str(quantiles_open[i]))
            ax1.set_xlim(0,2500)
            ax1.set_ylim(41,0)
            ax1.set_aspect(4)
            
            plt.setp(ax1.get_xticklabels(), visible=False)
            ax1.set_yticks(np.linspace(0,41,3))
            ax1.set_yticklabels(list(np.linspace(0,20,3)))
            ax1.set_ylabel('Depth [m]')
            
            pdb.set_trace()
            
            #Save the figure
            plt.savefig(fig_name,dpi=2000)
            plt.clf()
        
        #Plot the depth correctes traces
        #pdb.set_trace()
        
        #Define fig name
        fig_name=path_savefig+'referencetrace_depth_corr'+'.png'
        
        #Prepare the plot
        fig, (ax1) = plt.subplots(1, 1)
        fig.suptitle('Depth corrected traces, 2013 trace in MF2019')
        figManager = plt.get_current_fig_manager()
        figManager.window.showMaximized()
            
        ax1.imshow(dataframe['depth_corrected'],cmap=plt.get_cmap('gray'))#,norm=divnorm)
        ax1.imshow(dataframe['mask_truth'],cmap=plt.get_cmap('OrRd'), alpha=0.2)
        ax1.title.set_text(dataframe['datetrack'] +' - depth corrected')
        ax1.set_xlim(0,2500)
        ax1.set_ylim(41,0)
        ax1.set_aspect(4)
        plt.setp(ax1.get_xticklabels(), visible=False)
        ax1.set_yticks(np.linspace(0,41,3))
        ax1.set_yticklabels(list(np.linspace(0,20,3)))
        ax1.set_ylabel('Depth [m]')
        
        #Save the figure
        plt.savefig(fig_name,dpi=2000)
        plt.clf()
    
    if (show_case_study=='TRUE'):
        #pdb.set_trace()
        
        #Set dataframe
        dataframe={}
        dataframe={k: {} for k in list(investigation_year.keys())}
            
        for single_year in investigation_year.keys():
            
            #If no data, continue
            if (investigation_year[single_year]=='empty'):
                print('No data for year '+str(single_year)+', continue')
                continue
        
            print(single_year)
            #Define the quantiles to open
            quantiles_open=np.round(np.arange(0.63,0.83,0.01),2)
            
            #Set dataframe
            dataframe[single_year]={k: {} for k in list(quantiles_open)}
            
            #Construct the date for data loading
            start_date_track=investigation_year[single_year][0]
            end_date_track=investigation_year[single_year][-1]
            date_track=start_date_track[5:20]+'_'+end_date_track[17:20]
            
            #Define filename roll, depth corrected data and mask
            filename_depth_corrected=date_track+'_DEPTH_CORRECTED.pickle'
                         
            #Open the roll corrected file
            f_depth_corrected = open(path_depth_corrected+filename_depth_corrected, "rb")
            depth_corrected = pickle.load(f_depth_corrected)
            f_depth_corrected.close()
            
            #pdb.set_trace()
            
            #Open the quantile files
            for indiv_quantile in quantiles_open:
                #Define filename of boolean files 
                filename_boolean_quantile=date_track+'_SG1_cutoff_'+str(indiv_quantile)+'_threshold_350.pickle'
                
                #Load boolean files
                f_boolean_quantile = open(path_boolean_remove_surf+filename_boolean_quantile, "rb")
                dataframe[single_year][indiv_quantile] = pickle.load(f_boolean_quantile)
                f_boolean_quantile.close()
            
            #########################################################################
            # From plot_2002_2003.py - END
            #########################################################################
            #pdb.set_trace()
            ###2. Load the latitude and longitude
            lat_appended=[]
            lon_appended=[]
                
            for indiv_file_load in investigation_year[single_year]:
                print(indiv_file_load)
                
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
                
                lat_appended=np.flipud(lat_appended)
                lon_appended=np.flipud(lon_appended)
                depth_corrected=np.fliplr(depth_corrected)
                
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
                
            #pdb.set_trace()
            #Store reunited lat/lon, slice output and mask in a dictionnary:
            dataframe[single_year]['lat_appended']=lat_appended
            dataframe[single_year]['lon_appended']=lon_appended
            dataframe[single_year]['depth']=depth
            #Store depth corrected traces into dataframe
            dataframe[single_year]['datetrack']=date_track
            dataframe[single_year]['depth_corrected']=depth_corrected
            
            ##############################################################################
            ###                          Load and organise data                        ###
            ##############################################################################
            
        #pdb.set_trace()
        path_savefig='C:/Users/jullienn/switchdrive/Private/research/RT1/masking_iceslabs/quantile_investigation/'
        
        #Plot data
        #Loop on the quantiles
        for i in range(0,len(quantiles_open)):
            #pdb.set_trace()
            print('Quantile: ',quantiles_open[i])
            #Prepare the plot
            fig, (ax1,ax2,ax3,ax4,ax5,ax6,ax7) = plt.subplots(7, 1)
            fig.suptitle('Custom threshold: quantile '+str(quantiles_open[i])+' of ice slabs distribution, SG1, 350 continuity from 2013 trace in MF2019')
            figManager = plt.get_current_fig_manager()
            figManager.window.showMaximized()

            #Define fig name
            fig_name=path_savefig+'ref_year_overlap_quant_'+str(quantiles_open[i])+'.png'
            
            #Loop on the years
            for single_year in investigation_year.keys():
                
                #If no data, continue
                if (investigation_year[single_year]=='empty'):
                    print('No data for year '+str(single_year)+', continue')
                    continue
                
                print(single_year)
                
                if (single_year==2010):
                    ax_plotting=ax1
                    ylim_down=86
                elif (single_year==2011):
                    ax_plotting=ax2
                    ylim_down=86
                elif (single_year==2012):
                    ax_plotting=ax3
                    ylim_down=41
                elif (single_year==2013):
                    ax_plotting=ax4
                    ylim_down=41
                elif (single_year==2014):
                    ax_plotting=ax5
                    ylim_down=41
                elif (single_year==2017):
                    ax_plotting=ax6
                    ylim_down=41
                elif (single_year==2018):
                    ax_plotting=ax7
                    ylim_down=41
                else:
                    print('Year not existing')
                
                #pdb.set_trace()
                
                #Replace where dry firn by nan so that overlay plot can be possible
                quantile_to_plot=dataframe[single_year][quantiles_open[i]]
                quantile_to_plot[quantile_to_plot==0]=np.nan
                
                #Find x and y lim for imshow plotting
                xlimits_start=np.where(dataframe[single_year]['lon_appended'] >= -47.8)[0][0]
                xlimits_end=np.where(dataframe[single_year]['lon_appended'] <= -46.8)[0][-1]
                #xlimits=(dataframe[single_year]['lon_appended'] >= -47.8)&(dataframe[single_year]['lon_appended'] <= -46.8)
                
                #Plot custom threshold ice slabs identification
                ax_plotting.imshow(dataframe[single_year]['depth_corrected'],cmap=plt.get_cmap('gray'))#,norm=divnorm)
                ax_plotting.imshow(quantile_to_plot,cmap=plt.get_cmap('autumn'),alpha=0.1)#,norm=divnorm)
                #ax_plotting.title.set_text(dataframe[single_year]['datetrack'])
                ax_plotting.set_xlim(xlimits_start,xlimits_end)
                ax_plotting.set_ylim(ylim_down,0)
                ax_plotting.set_aspect(np.abs(xlimits_start-xlimits_end)/ylim_down/17)
                #pdb.set_trace()
                plt.setp(ax_plotting.get_xticklabels(), visible=False)
                ax_plotting.set_yticks(np.linspace(0,ylim_down,3))
                ax_plotting.set_yticklabels(list(np.linspace(0,20,3)))
                ax_plotting.set_ylabel('Depth [m]')
                
                #pdb.set_trace()
                #change figsize to make it bigger while saving

                
            #plt.show()
            pdb.set_trace()
            #Save the figure
            plt.savefig(fig_name,dpi=2000)
            plt.clf()
            
            
pdb.set_trace()
    
if (display_pickle=='TRUE'):
    
    #Compute the speed (Modified Robin speed):
    # self.C / (1.0 + (coefficient*density_kg_m3/1000.0))
    v= 299792458 / (1.0 + (0.734*0.873/1000.0))
    
    #Define paths
    path_data='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/'
    path_depth_corrected='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2010_2018/pickles_and_images/Depth_Corrected_Picklefiles/'
    path_boolean_MacFerrin='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2010_2018/pickles_and_images/Boolean_Array_Picklefiles/'
    path_boolean_remove_surf='C:/Users/jullienn/switchdrive/Private/research/RT1/remove_surface_return/'

    dataframe={}
    
    for single_year in investigation_year.keys():
        print(single_year)
        
        #If no data, continue
        if (investigation_year[single_year]=='empty'):
            print('No data for year '+str(single_year)+', continue')
            continue
        
        ###1. Load the mask and depth_corrected files:
        start_date_track=investigation_year[single_year][0]
        end_date_track=investigation_year[single_year][-1]
        date_track=start_date_track[5:20]+'_'+end_date_track[17:20]
        
        #Define filename depth corrected and boolean files from Macferrin
        filename_depth_corrected=date_track+'_DEPTH_CORRECTED.pickle'
        filename_boolean_MF=date_track+'_SG1_CUTOFF_-0.45_THRESHOLD_350.pickle'
        
        #Open the depth corrected file from Macferrin
        f_depth_corrected = open(path_depth_corrected+filename_depth_corrected, "rb")
        depth_corrected = pickle.load(f_depth_corrected)
        f_depth_corrected.close()
        
        #Load boolean files from Macferrin
        f_boolean_MF = open(path_boolean_MacFerrin+filename_boolean_MF, "rb")
        boolean_MF = pickle.load(f_boolean_MF)
        f_boolean_MF.close()


        #Define filename of boolean files 
        filename_boolean_c=date_track+'_SG1_cutoff_0.65_threshold_350.pickle'
        
        #Load boolean files
        f_boolean_c = open(path_boolean_remove_surf+filename_boolean_c, "rb")
        boolean_c_DF = pickle.load(f_boolean_c)
        f_boolean_c.close()
        
        #########################################################################
        # From plot_2002_2003.py - END
        #########################################################################
        
        ###2. Load the latitude and longitude
        lat_appended=[]
        lon_appended=[]
            
        for indiv_file_load in investigation_year[single_year]:
            print(indiv_file_load)
            
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
            
            lat_appended=np.flipud(lat_appended)
            lon_appended=np.flipud(lon_appended)
            depth_corrected=np.fliplr(depth_corrected)
            boolean_MF=np.fliplr(boolean_MF)
            #Dry firn boolean do not need to be flipped because they are from depthcorrected which where already flipped over
            
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
                                     'depth_corrected':depth_corrected,
                                     'SG1_CUTOFF_-0.45_THRESHOLD_350':boolean_MF,
                                     'remove_surf_proc':boolean_c_DF,
                                     'datetrack':date_track}
        
        
        #Calculate the min longitude
        min_lon_temp=np.ndarray.min(lon_appended)
        if(min_lon_temp<min_lon):
            min_lon=min_lon_temp
            #print('Min is:'+str(min_lon))
        
        #Calculate the max longitude
        max_lon_temp=np.ndarray.max(lon_appended)
        if(max_lon_temp>max_lon):
            max_lon=max_lon_temp
        ##############################################################################
        ###                          Load and organise data                        ###
        ##############################################################################
    
    '''
    pdb.set_trace()
    
    fig1, (ax1,ax2,ax3,ax4,ax5) = plt.subplots(5, 1)
    fig1.suptitle('2018 - Custom threshold from calibration trace in MacFerrin et al., 2019')
    
    cb1=ax1.imshow(dataframe['2018']['depth_corrected'],cmap=plt.get_cmap('gray'))#,norm=divnorm)
    #fig1.colorbar(cb1,ax=ax1)
    ax1.set_ylim(61,0)
    ax1.set_xlim(1000,4520) 
    ax1.set_title('Depth corrected from MacFeerin et al., 2019')
    
    ax2.imshow(dataframe['2018']['remove_surf_proc1'],cmap=plt.get_cmap('gray_r'))#,norm=divnorm)
    ax2.set_ylim(61,0)
    ax2.set_xlim(1000,4520) 
    ax2.set_title('Custom threshold: quantile 0.5 of ice slabs distribution, SG1, 350 continuity from 2013 trace in MF2019')
    
    ax3.imshow(dataframe['2018']['remove_surf_proc2'],cmap=plt.get_cmap('gray_r'))#,norm=divnorm)
    ax3.set_ylim(61,0)
    ax3.set_xlim(1000,4520) 
    ax3.set_title('Custom threshold: hollow bewteen ice slabs  and dry firn distribution, SG1, 350 continuity from 2013 trace in MF2019')
    
    ax4.imshow(dataframe['2018']['remove_surf_proc3'],cmap=plt.get_cmap('gray_r'))#,norm=divnorm)
    ax4.set_ylim(61,0)
    ax4.set_xlim(1000,4520) 
    ax4.set_title('Custom threshold: quantile 0.75 of ice slabs distribution, SG1, 350 continuity from 2013 trace in MF2019')
    
    ax5.imshow(dataframe['2018']['SG1_CUTOFF_-0.45_THRESHOLD_350'],cmap=plt.get_cmap('gray_r'))#,norm=divnorm)
    ax5.set_ylim(61,0)
    ax5.set_xlim(1000,4520)   
    ax5.set_title('MacFerrin et al., 2019: SG1_CUTOFF_-0.45_THRESHOLD_350')
    
    #ax1.set_ylim(129,0) for 2010, 2011
    #ax1.set_ylim(61,0) for remaining years
    
    plt.show()
    
    pdb.set_trace()
   '''
    
    
    #Prepare figures
    figde, (ax1de,ax2de,ax3de,ax4de,ax5de,ax6de,ax7de) = plt.subplots(7, 1)
    figde.suptitle('Depth corrected from MacFerrin et al., 2019')
    
    figd, (ax1d,ax2d,ax3d,ax4d,ax5d,ax6d,ax7d) = plt.subplots(7, 1)
    figd.suptitle('Custom threshold: quantile 0.65 of ice slabs distribution, SG1, 350 continuity from 2013 trace in MF2019')
    
    figm, (ax1m,ax2m,ax3m,ax4m,ax5m,ax6m,ax7m) = plt.subplots(7, 1)
    figm.suptitle('MacFerrin et al., 2019: SG1_CUTOFF_-0.45_THRESHOLD_350')
    
    for single_year in investigation_year.keys():
        
        #If no data, continue
        if (investigation_year[single_year]=='empty'):
            print('No data for year '+str(single_year)+', continue')
            continue
        
        print(single_year)
        
        #Associate axs for plotting
        if (single_year==2010):
            ax_plotting_dry=ax1d
            ax_plotting_MF=ax1m
            ax_plotting_depth=ax1de
        elif (single_year==2011):
            ax_plotting_dry=ax2d
            ax_plotting_MF=ax2m
            ax_plotting_depth=ax2de
        elif (single_year==2012):
            ax_plotting_dry=ax3d
            ax_plotting_MF=ax3m
            ax_plotting_depth=ax3de
        elif (single_year==2013):
            ax_plotting_dry=ax4d
            ax_plotting_MF=ax4m
            ax_plotting_depth=ax4de
        elif (single_year==2014):
            ax_plotting_dry=ax5d
            ax_plotting_MF=ax5m
            ax_plotting_depth=ax5de
        elif (single_year==2017):
            ax_plotting_dry=ax6d
            ax_plotting_MF=ax6m
            ax_plotting_depth=ax6de
        elif (single_year==2018):
            ax_plotting_dry=ax7d
            ax_plotting_MF=ax7m
            ax_plotting_depth=ax7de
        else:
            print('Year not existing')
        
        '''
        ax_plotting_dry.imshow(dataframe[str(single_year)]['remove_surf_proc3'],cmap=plt.get_cmap('gray_r'))
        
        ax_plotting_MF.imshow(dataframe[str(single_year)]['SG1_CUTOFF_-0.45_THRESHOLD_350'],cmap=plt.get_cmap('gray_r'))      
        
        ax_plotting_depth.imshow(dataframe[str(single_year)]['depth_corrected'],cmap=plt.get_cmap('gray'))
        if (str(single_year) in list(['2010','2011'])):
            ax_plotting_depth.set_ylim(86,0)
        else:
            ax_plotting_depth.set_ylim(41,0)
        '''
        
        #Plot custom threshold ice slabs identification
        X=dataframe[str(single_year)]['lon_appended']
        Y=np.arange(0,20,20/dataframe[str(single_year)]['SG1_CUTOFF_-0.45_THRESHOLD_350'].shape[0])
        
        ax_plotting_dry.pcolor(X,Y,dataframe[str(single_year)]['remove_surf_proc'],cmap=plt.get_cmap('gray_r'))#,norm=divnorm)
        ax_plotting_dry.title.set_text(dataframe[str(single_year)]['datetrack'])
        ax_plotting_dry.set_xlim(-47.8,-46.8)
        ax_plotting_dry.set_ylim(20,0)
        ax_plotting_dry.set_aspect(0.002)
        
        
        #Plot MacFerrin's boolean
        ax_plotting_MF.pcolor(X,Y,dataframe[str(single_year)]['SG1_CUTOFF_-0.45_THRESHOLD_350'],cmap=plt.get_cmap('gray_r'))#,norm=divnorm)
        ax_plotting_MF.title.set_text(dataframe[str(single_year)]['datetrack'])
        ax_plotting_MF.set_xlim(-47.8,-46.8)
        ax_plotting_MF.set_ylim(20,0)
        ax_plotting_MF.set_aspect(0.002)
        
        #Plot Depth corrected traces from MacFerrin et al., 2019
        Y=np.arange(0,100,100/dataframe[str(single_year)]['depth_corrected'].shape[0])
        
        ax_plotting_depth.pcolor(X,Y,dataframe[str(single_year)]['depth_corrected'],cmap=plt.get_cmap('gray'))#,norm=divnorm)
        ax_plotting_depth.title.set_text(dataframe[str(single_year)]['datetrack'])
        ax_plotting_depth.set_xlim(-47.8,-46.8)
        ax_plotting_depth.set_ylim(20,0)
        ax_plotting_depth.set_aspect(0.002)
        
        '''
        #Display my remove_surf boolean file
        X=dataframe[str(single_year)]['lon_appended']
        Y=np.arange(0,20,20/dataframe[str(single_year)]['remove_surf_SG1_CUTOFF_-0.45_THRESHOLD_350'].shape[0])
        
        ax_plotting_dry.pcolor(X,Y,dataframe[str(single_year)]['remove_surf_SG1_CUTOFF_-0.45_THRESHOLD_350'],cmap=plt.get_cmap('gray_r'))#,norm=divnorm)
        ax_plotting_dry.invert_yaxis() #Invert the y axis = avoid using flipud.
        ax_plotting_dry.title.set_text(dataframe[str(single_year)]['datetrack'])
        #ax_plotting_dry.set_xlim(min_lon,max_lon)
        
        #Display MacFerrin's boolean file
        X=dataframe[str(single_year)]['lon_appended']
        Y=np.arange(0,20,20/dataframe[str(single_year)]['SG1_CUTOFF_-0.45_THRESHOLD_350'].shape[0])
        
        ax_plotting_MF.pcolor(X,Y,~dataframe[str(single_year)]['SG1_CUTOFF_-0.45_THRESHOLD_350'],cmap=plt.get_cmap('gray'))#,norm=divnorm)
        ax_plotting_MF.invert_yaxis() #Invert the y axis = avoid using flipud.
        ax_plotting_MF.title.set_text(dataframe[str(single_year)]['datetrack'])
        #ax_plotting_MF.set_xlim(min_lon,max_lon)

        #Display MacFerrin's depth corrected traces
        X=dataframe[str(single_year)]['lon_appended']
        Y=np.arange(0,100,100/dataframe[str(single_year)]['depth_corrected'].shape[0])
        
        ax_plotting_depth.pcolor(X,Y,dataframe[str(single_year)]['depth_corrected'],cmap=plt.get_cmap('gray'))#,norm=divnorm)
        ax_plotting_depth.invert_yaxis() #Invert the y axis = avoid using flipud.
        ax_plotting_depth.title.set_text(dataframe[str(single_year)]['datetrack'])
        #ax_plotting_MF.set_xlim(min_lon,max_lon)
        ax_plotting_depth.set_ylim(20,0)
        '''
        
        
    plt.show()










    