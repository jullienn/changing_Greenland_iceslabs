# -*- coding: utf-8 -*-
"""
Created on Wed Nov  3 15:21:23 2021

@author: jullienn
"""


def reconstruct_with_NaNs(full_slice,mask,slice_without_NaNs):
    #Where mask is False
    index_false=np.where(mask==False)[0]
    #Where mask is True
    index_true=np.where(mask==True)[0]
    
    #reconstruct the array with the NaNs
    traces_20m_full=np.zeros((slice_without_NaNs.shape[0],full_slice.shape[1]))
    traces_20m_full[:,index_false]=np.nan
    traces_20m_full[:,index_true]=slice_without_NaNs
    
    return traces_20m_full


def select_20m_slice_without_NaNs(slice_input,depth,mask):
    #I should input slices without the NaNs in it 
    #Where mask is False
    index_false=np.where(mask==False)[0]
    #Where mask is True
    index_true=np.where(mask==True)[0]
    
    #Prepare the mask for appliance
    mask.shape=mask.shape[0],1
    mask_matrix=np.repeat(mask,slice_input.shape[0],axis=1)
    mask_matrix=np.transpose(mask_matrix)
    #Apply the mask
    slice_input_masked=slice_input[mask_matrix]
    slice_input_masked=slice_input_masked.reshape((slice_input.shape[0],(slice_input.shape[1])-len(index_false)))
    
    #Keep only the first 20m of slices
    boolean_20m=depth <= 20
    ind_20m=np.where(boolean_20m==True)[0]
    
    traces=slice_input_masked
    traces_20m=traces[ind_20m,:]
    
    return traces_20m

def identify_ice_lenses(traces,slices_depth_corrected_after_surf_removal_without_norm,depth,mask,datetrack,quantile_investigation,desired_quantiles):
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
        
        print('      Creating the pickle files of quantile',cutoff_q)
        
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
        ice_lenses_above_cutoff_size = np.zeros(boolean_traces.shape, dtype=bool)

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
        boolean_full_slabs=reconstruct_with_NaNs(slices_depth_corrected_after_surf_removal_without_norm,mask,boolean_slabs)
        #Ok that works!
        
        #format cutoff_q for name saving
        cutoff_q_save=str(np.round(cutoff_q,2))
        if (len(cutoff_q_save)==3):
            cutoff_q_save=cutoff_q_save+'0'
        
        #pdb.set_trace()
        #Save as pickle file
        
        filename_tosave='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2010_2018/ii_out_from_iceslabs_processing_jullien.py/pickles/'+datetrack+'_'+algorithm+'_cutoffisquantile_'+cutoff_q_save+'_threshold_'+str(continuity_threshold)+'.pickle'
        '''
        filename_tosave='/flash/jullienn/data/threshold_processing_output/pickles/'+datetrack+'_'+algorithm+'_cutoffisquantile_'+cutoff_q_save+'_threshold_'+str(continuity_threshold)+'.pickle'
        '''
        outfile= open(filename_tosave, "wb" )
        pickle.dump(boolean_full_slabs,outfile)
        outfile.close()
        
        #Save figures here
        #Define fig name
        
        fig_name='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2010_2018/ii_out_from_iceslabs_processing_jullien.py/images/'+datetrack+'_'+algorithm+'_cutoffisquantile_'+cutoff_q_save+'_threshold_'+str(continuity_threshold)+'.png'
        '''
        fig_name='/flash/jullienn/data/threshold_processing_output/images/'+datetrack+'_'+algorithm+'_cutoffisquantile_'+cutoff_q_save+'_threshold_'+str(continuity_threshold)+'.png'
        '''
        #Prepare the plot
        fig, (ax1) = plt.subplots(1, 1)
        '''
        figManager = plt.get_current_fig_manager()
        figManager.window.showMaximized()
        '''
        
        #Replace where dry firn by nan so that overlay plot can be possible
        quantile_to_plot=boolean_full_slabs
        quantile_to_plot[boolean_full_slabs==0]=np.nan
        
        #Plot custom threshold ice slabs identification
        ax1.imshow(slices_depth_corrected_after_surf_removal_without_norm,cmap=plt.get_cmap('gray'))#,norm=divnorm)
        ax1.imshow(quantile_to_plot,cmap=plt.get_cmap('autumn'),alpha=0.1)#,norm=divnorm)
        ax1.title.set_text(datetrack+' - quantile: '+cutoff_q_save)
        ax1.set_ylim(boolean_full_slabs.shape[0],0)
        #ax1.set_aspect(boolean_full_slabs.shape[1]/boolean_full_slabs.shape[0]/17)
        
        plt.setp(ax1.get_xticklabels(), visible=False)
        ax1.set_yticks(np.linspace(0,boolean_full_slabs.shape[0],3))
        ax1.set_yticklabels(list(np.linspace(0,20,3)))
        ax1.set_ylabel('Depth [m]')
                
        #Save the figure
        plt.savefig(fig_name,dpi=2000)
        plt.close(fig)
        
        #Update count
        count=count+1
        
    return

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
    group_id_array = np.zeros(boolean_image.shape, dtype=int)
    visited_mask_empty = np.zeros(boolean_image.shape, dtype=bool)
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


import pickle
import matplotlib.pyplot as plt
import numpy as np
import pdb
import pandas as pd
import scipy.io
import h5py
import sklearn.preprocessing

export_rescaled_pickes='TRUE'

v= 299792458 / (1.0 + (0.734*0.873/1000.0))

#List of traces where iceslabs likelihood identification have failed
list_trace=list(['20110416_01_053_055','20120421_01_052_052','20130423_01_125_125',
                  '20130423_01_127_127','20130426_01_089_089','20140419_01_016_017',
                  '20140419_01_028_028','20140419_03_075_075','20140516_02_031_034',
                  '20180419_02_032_033','20180419_02_035_036','20180425_01_166_169',
                  '20180427_01_170_172','20180429_01_008_014','20180427_01_004_006'])#The processing on the last one works well, this is for comparison

#Define paths
path_depth_corrected='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2010_2018/ii_out_from_iceslabs_processing_jullien.py/pickles/'
path_datetrack='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/'
path_data='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/'
path_mask='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2010_2018/i_out_from_IceBridgeGPR_Manager_v2.py/pickles_and_images/Boolean_Array_Picklefiles/'

#Read datetracks
datetrack_toread = np.asarray(pd.read_csv(path_datetrack+'datetrack_20102018.txt', header=None))

#Open the quantile file over which we will loop
quantiles_file=np.arange(0.63,0.82,0.01)
filename_quantiles='C:/Users/jullienn/switchdrive/Private/research/RT1/masking_iceslabs/quantiles_threshold_application/quantile_file_'+str(np.round(quantiles_file[0],2))+'_'+str(np.round(quantiles_file[-1],2))+'.txt'
'''
filename_quantiles='/flash/jullienn/data/threshold_processing/quantile_file_'+str(np.round(quantiles_file[0],2))+'_'+str(np.round(quantiles_file[-1],2))+'.txt'
'''
quantile_file = np.asarray(pd.read_csv(filename_quantiles, sep=" ", header=None))

appended_radar_slices=[]

#Loop over the dates
for indiv_date in datetrack_toread:
    
    print(indiv_date[0])
    
    #Define filename
    filename_depth_corrected=indiv_date[0]+'_Depth_Corrected_surf_removal.pickle'
    
    #Open the depth corrected file of the corresponding date
    f_depth_corrected = open(path_depth_corrected+filename_depth_corrected, "rb")
    slice_depth_corrected = pickle.load(f_depth_corrected)
    f_depth_corrected.close()
    
    '''
    #Set ylim
    if indiv_date[0][0:4] in list(['2010','2011']):
        custom_ylim=129
    else:
        custom_ylim=61
    #ylim = 129 if 2010, 2011, else it is 61
    
    #Extract slice to plot
    slice_depth_corrected_toplot=slice_depth_corrected[0:custom_ylim,:]

    #Reshape into a vector
    slice_depth_corrected_array=np.asarray(slice_depth_corrected_toplot).reshape(-1)

    #Extract the value of the 5th and 95th quantile
    appended_radar_slices=np.append(appended_radar_slices,slice_depth_corrected_array)
    '''
    #Reshape into a vector
    slice_depth_corrected_array=np.asarray(slice_depth_corrected).reshape(-1)
    
    #Extract the value of the 5th and 95th quantile
    appended_radar_slices=np.append(appended_radar_slices,slice_depth_corrected_array)
    
#open all the 20m depth corrected pickles files and extract the distributution of signal return
#then rescale the dates which are problematic with the e.g. 5-95 percentiles of the ditribution
#Thsi requires to create the deèth corrected files

for indiv_file in list_trace:
    #pdb.set_trace()
    
    #Define filename
    filename_depth_corrected=indiv_file+'_Depth_Corrected_surf_removal.pickle'
    #Open the depth corrected file of the corresponding date
    f_depth_corrected = open(path_depth_corrected+filename_depth_corrected, "rb")
    slice_depth_corrected = pickle.load(f_depth_corrected)
    f_depth_corrected.close()
    
    
    filename_mask=indiv_file+'_mask.pickle'
    #Open the mask file of the corresponding date
    f_mask = open(path_mask+filename_mask, "rb")
    mask = pickle.load(f_mask)
    f_mask.close()  
    
    ###2. Load the time
    #Create list of dates
    start_trace=int(indiv_file[12:15])
    end_trace=int(indiv_file[16:19])
    
    single_year=int(indiv_file[0:4])
    
    lat_appended=[]
    lon_appended=[]
        
    for nb_trace in np.arange(start_trace,end_trace+1,1):
        
        #Reconstruct the name of the file
        if (nb_trace<10):
            indiv_file_load='Data_'+indiv_file[0:12]+'00'+str(nb_trace)+'.mat'
        elif ((nb_trace>=10)&(nb_trace<100)):
            indiv_file_load='Data_'+indiv_file[0:12]+'0'+str(nb_trace)+'.mat'
        else:
            indiv_file_load='Data_'+indiv_file[0:12]+str(nb_trace)+'.mat'
        
        print('  ',indiv_file_load)
        #pdb.set_trace()
        #Create the path
        path_raw_data=path_data+str(single_year)+'_Greenland_P3/CSARP_qlook/'+indiv_file_load[5:16]+'/'
        
        #Load data
        if (single_year>=2014):
            
            fdata_filename = h5py.File(path_raw_data+indiv_file_load)
            time_filename=fdata_filename['Time'][:,:]
            
        else:
            fdata_filename = scipy.io.loadmat(path_raw_data+indiv_file_load)
            time_filename = fdata_filename['Time']
    
    #3. Calculate the depth from the time
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
        
    #Identify the index where depth <=20m
    boolean_20m=depth <= 20
    ind_20m=np.where(boolean_20m==True)[0]

    #Extract slice to plot
    slice_depth_corrected_20m=slice_depth_corrected[ind_20m,:]
    
    #Extract NaN from slice
    slice_depth_corrected_20m_without_nans=select_20m_slice_without_NaNs(slice_depth_corrected_20m,depth,mask)

    #Apply rescaling
    rescaled_slice=sklearn.preprocessing.minmax_scale(slice_depth_corrected_20m_without_nans, feature_range=(np.nanpercentile(appended_radar_slices,5),np.nanpercentile(appended_radar_slices,95)))
    
    #Reconsruct full slice with NaNs
    full_rescaled_slice=reconstruct_with_NaNs(slice_depth_corrected_20m,mask,rescaled_slice)
    
    #Display the results
    
    #Plot to briefly check
    fig, (ax1,ax2) = plt.subplots(2, 1)
    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()
    fig.suptitle(indiv_file)
    
    cb=ax1.imshow(slice_depth_corrected_20m,cmap=plt.get_cmap('gray'))
    fig.colorbar(cb,ax=ax1)
    ax1.set_ylim(ind_20m[-1],0)
    ax1.set_title('Depth corrected trace')
       
    cb2=ax2.imshow(full_rescaled_slice,cmap=plt.get_cmap('gray'))
    fig.colorbar(cb2,ax=ax2)
    ax2.set_ylim(ind_20m[-1],0)
    ax2.set_title('Rescaled slice 5-95th percentiles')
    
    plt.show()
    
    if (export_rescaled_pickes=='TRUE'):
        #Save the depth corrected rescaled slices as pickle files
        filename_tosave=path_depth_corrected+indiv_file[0:19]+'_Depth_Corrected_surf_removal_rescaled.pickle'
        outfile= open(filename_tosave, "wb" )
        pickle.dump(full_rescaled_slice,outfile)
        outfile.close()
        print('   Exporting '+indiv_file[0:19]+' depth corrected pickle file')
    
    ##########################################################################
    ###                From iceslabs_processing_jullien.py                 ###
    ##########################################################################
    #7.Perform ice slabs identification (thresholding and smoothing)
    print('   Perform iceslabs identification')
    
    #Extract the 20m slices and get rid of exclusions
    traces_20m=select_20m_slice_without_NaNs(full_rescaled_slice,depth,mask)
    #Identify ice slabs
        
    boolean_slabs=identify_ice_lenses(traces_20m,full_rescaled_slice,depth,mask,indiv_file[0:19],quantile_file[1,:],quantile_file[0,:])
    ##########################################################################
    ###                From iceslabs_processing_jullien.py                 ###
    ##########################################################################  
    
    
    
    

#Plot the distribution of 2010-2018 radar signal strenght
fig, (ax1) = plt.subplots(1, 1)
figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()

ax1.set_title('2010-2018 depth corrected radar signal strength distribution')
ax1.hist(appended_radar_slices,bins=100)
plt.show()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    