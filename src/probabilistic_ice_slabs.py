# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 08:26:05 2021

@author: jullienn
"""

def _export_to_8bit_array(array):
    #This is from MacFerrin et al., 2019, IceBridgeGPR_Manager_v2.py
    '''In order to export a function to a PNG image, use this funciton to
    export to an 8 bit unsigned integer array of scaled values.'''

    output_array = np.zeros(array.shape, dtype=np.uint8)
    excluded_mask = np.isnan(array)

    range_min = 0
    range_max = 2**8 - 1
    # Get the data minimum and maximum while cutting off 0.5% of outliers
    nonzero_values = array[~excluded_mask]
    #pdb.set_trace()
    data_cutoff_min = 0#np.percentile(nonzero_values,  0.5)
    data_cutoff_max = 1#np.percentile(nonzero_values, 99.5)

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

########################## From plot_2002_2003.py #############################
def compute_distances(eastings,northings):
    #This part of code is from MacFerrin et al., 2019
    '''Compute the distance (in m here, not km as written originally) of the traces in the file.'''
    # C = sqrt(A^2  + B^2)
    distances = np.power(np.power((eastings[1:] - eastings[:-1]),2) + np.power((northings[1:] - northings[:-1]),2), 0.5)
    #Calculate the cumsum of the distances
    cumsum_distances=np.nancumsum(distances)
    #Seeting the first value of the cumsum to be zero as it is the origin
    return_cumsum_distances=np.zeros(eastings.shape[0])
    return_cumsum_distances[1:eastings.shape[0]]=cumsum_distances
    
    return return_cumsum_distances
########################## From plot_2002_2003.py #############################

#1. Open the data
#2. Loop over dates and load data
#3. Loop over quantiles and calculate probability
#4. Create probabilitic slabs plots and excels files

import pandas as pd
import numpy as np
import pdb
import pickle
import png

import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

import time
import os.path
import glob

generate_probability_iceslabs_files='TRUE'
apply_dry_firn_exclusions='TRUE'
generate_excel_file='FALSE'

#Identify all the datetraces to process
'''
path_datetrack='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/Exclusion_folder/'
'''
path_datetrack='/flash/jullienn/data/threshold_processing/'

datetrack_toread = np.asarray(pd.read_csv(path_datetrack+'datetrack_20102018.txt', header=None))

#Read dry firn exclusions file
DF_exclusions = pd.read_csv(path_datetrack+'dry_firn_removal.txt',sep="\n", header=None)
#This is from https://stackoverflow.com/questions/55129640/read-csv-into-a-dataframe-with-varying-row-lengths-using-pandas
DF_exclusions = DF_exclusions[0].str.split(' ', expand=True)
#Create colnames
len_exclusions = [str(i) for i in np.arange(0,DF_exclusions.shape[1]-1)] #from https://www.geeksforgeeks.org/python-converting-all-strings-in-list-to-integers/
colnames=['indiv_dates']+len_exclusions #from https://stackoverflow.com/questions/1720421/how-do-i-concatenate-two-lists-in-python
DF_exclusions.columns=colnames

if (generate_probability_iceslabs_files=='TRUE'):
    #I. Define path, open datetracks and define desired quantiles
    #Define path where to pick roll corrected data
    '''
    path_quantiles_data='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/to_delete/'
    '''
    path_quantiles_data='/flash/jullienn/data/threshold_processing_output/pickles/'
    
    #Define the desired quantiles over which we will loop
    desired_quantiles=np.round(np.arange(0.65,0.79,0.01),2)
    
    #intialize counter to 0
    count_time=0
        
    #II. Loop over these datetracks, and perform probability calculation:
    for indiv_trace in datetrack_toread:
        
        #We want to process only 2017
        if (not(indiv_trace[0][0:4] in list(['2017']))):
            print(indiv_trace[0],' not 2017, continue')
            continue
        
        #If pickle files have already been created, do not process and continue
        filename_to_check='/flash/jullienn/data/threshold_processing_output/probability_iceslabs/before_DF_appliance/pickles/'+indiv_trace[0]+'*'
        #filename_to_check='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2010_2018/iii_out_from_probabilistic_iceslabs.py/pickles/'+indiv_trace[0]+'_probability_iceslabs_presence.pickle'
        #filename_to_check='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/to_delete/'+indiv_trace[0]+'_probability_iceslabs_presence.pickle'
        
        if (len(glob.glob(filename_to_check))>0):
            print(indiv_trace[0],': file already existent, move on to the next date')
            continue
        
        #To access advance
        start = time.time()
        print(indiv_trace[0])
        print(count_time/len(datetrack_toread)*100,'%')
        
        #Must open the first quantile (0.65) to know the size 
        filename_quantile_open=indiv_trace[0]+'_SG1_cutoffisquantile_'+str(0.65)+'_threshold_350.pickle'
                             
        #Open the corresponding quantile 0.65 file
        f_quantile = open(path_quantiles_data+filename_quantile_open, "rb")
        indiv_quantile065_slice = pickle.load(f_quantile)
        f_quantile.close()
        
        #Set the probabilistic_slice to zeros
        probabilistic_slice=np.zeros((indiv_quantile065_slice.shape[0],indiv_quantile065_slice.shape[1]))
                
        #Loop over the quantiles, load data and perform probability calculation
        for indiv_quantile in desired_quantiles:
            #print(str('%.2f' % indiv_quantile))
            #Define filename of quantiles of interest
            filename_quantile_open=indiv_trace[0]+'_SG1_cutoffisquantile_'+str('%.2f' % indiv_quantile)+'_threshold_350.pickle'
            
            #Open the corresponding quantile file
            f_quantile = open(path_quantiles_data+filename_quantile_open, "rb")
            indiv_quantile_slice=pickle.load(f_quantile)
            f_quantile.close()
            
            #Add up the numbers
            probabilistic_slice=probabilistic_slice+indiv_quantile_slice
        
        #Divide the probabilistic_slice by the number of quantiles to have a probability map
        probabilistic_slice=probabilistic_slice/len(desired_quantiles)
        
        #Save the image
        '''
        fig_name=path_quantiles_data+'/prob/'+indiv_trace[0]+'_probability_iceslabs_presence.png'
        '''
        fig_name='/flash/jullienn/data/threshold_processing_output/probability_iceslabs/before_DF_appliance/images/'+indiv_trace[0]+'_probability_iceslabs_presence.png'
        
        #Prepare matrix for png plot. (1-probabilistic_slice) because 1 is white
        #out of the function _export_to_8bit_array, and I want black
        probabilistic_slice_png=_export_to_8bit_array((1-probabilistic_slice))
        
        #Save the image
        png_to_save=png.from_array(probabilistic_slice_png, mode='L')
        png_to_save.save(fig_name)
        
        #Save the pickle file
        '''
        filename_tosave=path_quantiles_data+'/prob/'+indiv_trace[0]+'_probability_iceslabs_presence.pickle'
        '''
        filename_tosave='/flash/jullienn/data/threshold_processing_output/probability_iceslabs/before_DF_appliance/pickles/'+indiv_trace[0]+'_probability_iceslabs_presence.pickle'
        
        outfile= open(filename_tosave, "wb" )
        pickle.dump(probabilistic_slice,outfile)
        outfile.close()
                
        ##################### Apply dry firn exclusions #######################
        #Apply dry firn exclusions here and save results as images and pickles
        if (apply_dry_firn_exclusions=='TRUE'):
            print('starting dry firn removal')
            
            #1. Extract exclusions
            #Extract total length of trace for dry firn mask generation
            total_len_trace=probabilistic_slice.shape[1]
            
            #Extract indiv_date and corresponding exclusions
            indiv_date=np.asarray(DF_exclusions[DF_exclusions['indiv_dates']==indiv_trace[0]])[0,0]
            exclusions=np.asarray(DF_exclusions[DF_exclusions['indiv_dates']==indiv_trace[0]])[0,1:]
            
            #Get rid of None in exclusions
            exclusions=exclusions[~(exclusions==None)]
            
            print(indiv_date)
            print(exclusions)
            print(total_len_trace)
            
            #2. Gather all exclusions
            #Loop through the exclusions
            for indiv_exclusion in exclusions:
                #Set indiv_mask to empty
                indiv_mask=[]
                
                if (indiv_exclusion.partition("-")[0]==''):
                    #-N case
                    indiv_mask='0-'+indiv_exclusion.partition("-")[-1]
                elif (indiv_exclusion.partition("-")[1]=='-'):
                    if (indiv_exclusion.partition("-")[-1]==''):
                        #N- case
                        indiv_mask=indiv_exclusion.partition("-")[0]+'-'+str(total_len_trace)
                    else:
                        #N-N case
                        indiv_mask=indiv_exclusion
                else:
                    print('Error!!')
                    break
                #3. Apply exclusions to the slice
                print(indiv_mask)
                start_exclusion=int(indiv_mask.partition("-")[0])
                end_exclusion=int(indiv_mask.partition("-")[-1])
                
                probabilistic_slice[:,start_exclusion:end_exclusion]=0
                probabilistic_slice_png[:,start_exclusion:end_exclusion]=255
            
            #4. Save final result
            
            #Save the image
            '''
            fig_name=path_quantiles_data+'/prob/'+indiv_trace[0]+'_probability_iceslabs_presence_after_DF.png'
            '''
            fig_name='/flash/jullienn/data/threshold_processing_output/probability_iceslabs/after_DF_appliance/images/'+indiv_trace[0]+'_probability_iceslabs_presence_after_DF.png'
            
            #Save the image
            png_to_save=png.from_array(probabilistic_slice_png, mode='L')
            png_to_save.save(fig_name)
            
            #Save the pickle file
            '''
            filename_tosave=path_quantiles_data+'/prob/'+indiv_trace[0]+'_probability_iceslabs_presence_after_DF.pickle'
            '''
            filename_tosave='/flash/jullienn/data/threshold_processing_output/probability_iceslabs/after_DF_appliance/pickles/'+indiv_trace[0]+'_probability_iceslabs_presence_after_DF.pickle'
            
            outfile= open(filename_tosave, "wb" )
            pickle.dump(probabilistic_slice,outfile)
            outfile.close()
        
        ##################### Apply dry firn exclusions #######################

    print('End of probabilistic processing')

if (generate_excel_file=='TRUE'):
    print('Generate excel file')

    ##############################################################################
    ###              Generate en excel file of ice slabs thickness             ###
    ##############################################################################
    # This is inspired from IceBridgeGPR_Manager_v2.py from MacFerrin et al., 2019
    
    import scipy.io
    import h5py
    from pyproj import Transformer
    
    #1. Open the excel file to populate
    #2. Loop over all the dates and perform filling of excel file at each iteration
    #3. Close the excel file
    
    #Define speed
    v= 299792458 / (1.0 + (0.734*0.873/1000.0))
    
    #Define path where data are stored
    '''
    path_probability_iceslabs='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2010_2018/iii_out_from_probabilistic_iceslabs.py/pickles/'
    path_data='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/'
    '''
    path_probability_iceslabs='/flash/jullienn/data/threshold_processing_output/probability_iceslabs/after_DF_appliance/pickles/'
    path_data='/flash/jullienn/data/threshold_processing/'
    
    #Define filename
    '''
    filename_excel_output='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2010_2018/iii_out_from_probabilistic_iceslabs.py/Ice_Layer_Output_Thicknesses_Likelihood_2010_2018_jullienetal2021.csv'
    '''
    filename_excel_output='/flash/jullienn/data/threshold_processing_output/probability_iceslabs/after_DF_appliance/Ice_Layer_Output_Thicknesses_2010_2018_jullienetal2021_low_estimate.csv'
    
    #Open filename (same procedure as MacFerrin et al., 2019)
    fout = open(filename_excel_output, 'w')
    header = "Track_name,Tracenumber,lat,lon,alongtrack_distance_m,20m_ice_content_m,likelihood\n"
    fout.write(header)
    
    #Loop over the traces
    for indiv_trace in datetrack_toread:
        
        #Open probability ice slabs pickle file
        filename_probability_open=indiv_trace[0]+'_probability_iceslabs_presence_after_DF.pickle'
        
        f_probability = open(path_probability_iceslabs+filename_probability_open, "rb")
        indiv_probability_slice=pickle.load(f_probability)
        f_probability.close()
                
        #Open depths and lat/lon
        #Create list of dates
        start_trace=int(indiv_trace[0][12:15])
        end_trace=int(indiv_trace[0][16:19])
        
        single_year=int(indiv_trace[0][0:4])
        
        lat_appended=[]
        lon_appended=[]
            
        for nb_trace in np.arange(start_trace,end_trace+1,1):
            
            #Reconstruct the name of the file
            if (nb_trace<10):
                indiv_file_load='Data_'+indiv_trace[0][0:12]+'00'+str(nb_trace)+'.mat'
            elif ((nb_trace>=10)&(nb_trace<100)):
                indiv_file_load='Data_'+indiv_trace[0][0:12]+'0'+str(nb_trace)+'.mat'
            else:
                indiv_file_load='Data_'+indiv_trace[0][0:12]+str(nb_trace)+'.mat'
            
            #pdb.set_trace()
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
        
                         
        #Transform the coordinated from WGS84 to EPSG:3413 for distance calculation
        #Example from: https://pyproj4.github.io/pyproj/stable/examples.html
        transformer = Transformer.from_crs("EPSG:4326", "EPSG:3413", always_xy=True)
        points=transformer.transform(np.array(lon_appended),np.array(lat_appended))
        
        #Reset the lat_3413 and lon_3413 to empty vectors.
        lon_3413=[]
        lat_3413=[]
        
        lon_3413=points[0]
        lat_3413=points[1]
        
        #Compute distances
        distances=compute_distances(lon_3413,lat_3413)
        
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
            
        #Compute depth_delta_m
        depth_delta_m = np.mean(depth[1:] - depth[:-1])
                
        #Let's transform the probabilistic ice slabs into an ice content
        #We must derive a low end and high end of ice slabs likelihood
        #for low end: slabs identified in 15 quantiles out of 15 => likelihood = 15/15=1
        #for high end: slabs identified in 1 quantile out of 15 => likelihood = 1/15 = 0.06667
        index_prob=indiv_probability_slice>=0.0665
        
        #Create slice full of nans
        slice_for_calculation=np.zeros((indiv_probability_slice.shape[0],indiv_probability_slice.shape[1]))
        #fill in slice_for_calculation by ones where likelihood >= chosen high/low end
        slice_for_calculation[index_prob]=1
        # Number of pixels times the thickness of each pixel
        ice_content_m = np.sum(slice_for_calculation, axis=0) * depth_delta_m
        
        #Derive likelihood averaged on the column
        slice_for_likelihood=np.zeros((indiv_probability_slice.shape[0],indiv_probability_slice.shape[1]))
        slice_for_likelihood[:]=np.nan
        slice_for_likelihood[index_prob]=indiv_probability_slice[index_prob]
        columnal_likelihoods=np.nanmean(slice_for_likelihood,axis=0)
        
        #Use the same names as MacFerrin et al., 2019
        lats=lat_appended
        lons=lon_appended
        ice_contents = ice_content_m
        
        # The one "really long" track has artifacts in the center that aren't real ice layers.  Filter these out.
        if indiv_trace[0] == "20120412_01_095_095":
            ice_contents[0:9000] = 0.0
            columnal_likelihoods[0:9000] = np.nan
            
        assert len(lats) == len(lons) == len(ice_contents) == len(columnal_likelihoods)
        tracenums = np.arange(len(lats), dtype=np.int64)
    
        tracecount = 0
        for lat, lon, tracenum, distance, ice_content, columnal_likelihood in zip(lats, lons, tracenums, distances, ice_contents, columnal_likelihoods):
            # Record ONLY traces that have > 1 m ice content in them.  We're not interested in thinner stuff here.
            # If it has > 16 m ice content (80%), we also omit it, just to keep pure ice out of it.
            if 1.0 <= ice_content <= 16.0:
                line = "{0},{1},{2},{3},{4},{5},{6}\n".format(indiv_trace[0], tracenum, lat, lon, distance, ice_content, columnal_likelihood)
                fout.write(line)
                tracecount += 1
        print(tracecount, "of", len(lats), "traces.")
        print()
    
    fout.close()
    print('End of probabilistic processing')

    ##############################################################################
    ###              Generate en excel file of ice slabs thickness             ###
    ##############################################################################
