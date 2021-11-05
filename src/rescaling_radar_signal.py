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


import pickle
import matplotlib.pyplot as plt
import numpy as np
import pdb
import pandas as pd
import scipy.io
import h5py
import sklearn.preprocessing

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
    rescaled_slice=sklearn.preprocessing.minmax_scale(slice_depth_corrected_20m_without_nans, feature_range=(np.nanpercentile(appended_radar_slices,0.05),np.nanpercentile(appended_radar_slices,0.95)))
    
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
    
    #Save the depth corrected rescaled slices as pickle files
    filename_tosave=path_depth_corrected+indiv_file[0:19]+'_Depth_Corrected_surf_removal_rescaled.pickle'
    outfile= open(filename_tosave, "wb" )
    pickle.dump(full_rescaled_slice,outfile)
    outfile.close()
    print('   Exporting '+indiv_file[0:19]+' depth corrected pickle file')

#Plot the distribution of 2010-2018 radar signal strenght
fig, (ax1) = plt.subplots(1, 1)
figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()

ax1.set_title('2010-2018 depth corrected radar signal strength distribution')
ax1.hist(appended_radar_slices,bins=100)
plt.show()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    