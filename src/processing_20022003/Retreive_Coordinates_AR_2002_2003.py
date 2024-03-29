# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 15:35:11 2020

@author: Nicolas Jullien
"""
import scipy.io
#import rasterio
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

##############################################################################
############################# Data manipulation ##############################
##############################################################################
import os
#pdb.set_trace()

#Define the path for working
path='D://OIB//AR'
os.chdir(path) # relative path: scripts dir is under Lab
os.getcwd()

# 1. Read the years to process
folder_years = [ f.name for f in os.scandir(path) if f.is_dir() ]

#Create the issue dataframe
df_issue=pd.DataFrame({'name_file_issue':pd.Series(['may11_03_40.mat','may11_03_9.mat','may12_03_13.mat','may12_03_17.mat',
                                                    'may12_03_38.mat','may12_03_39.mat','may12_03_44.mat','may13_03_30.mat',
                                                    'may14_03_11.mat','may14_03_15.mat','may14_03_52.mat','may15_03_3.mat',
                                                    'may15_03_37.mat','may18_02_30.mat','may09_03_24.mat','may09_03_42.mat',
                                                    'may11_03_15.mat','may11_03_18.mat']),
                        'index_issue':pd.Series([940,740,80,990,
                                                 850,40,310,330,
                                                 10,40,50,330,
                                                 750,230,320,330,
                                                 10,560])})

for folder_year in folder_years:
    folder_year_name=path+'//'+folder_year
    
    #Read the days for this specific year
    folder_days=[]
    folder_days = [ f.name for f in os.scandir(folder_year_name) if f.is_dir() ]
    #pdb.set_trace()
    
    for folder_day in folder_days:
        folder_day_name=folder_year_name+'//'+folder_day
        print('Folder name:')
        print(folder_day_name)
        #pdb.set_trace()
        
        #Go into the daily folders 
        os.chdir(folder_day_name)
        
        #Define the path
        mypath = folder_day_name+'//'
        
        #Save the filenames present in folder of interest (here May 13 2003)
        onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
        
        #If files have already been created, do not process and continue
        filename_to_check='D://OIB//2002_2003_export//'+folder_year+'//'+folder_day+'//'+onlyfiles[0].replace(".mat","")+"_aggregated"
        if (os.path.isfile(filename_to_check)):
            print('Aggregated files already existent, move on to the next date')
            continue
        
        if (folder_year=='2017_Greenland_P3'):
            print('2017, no need to process, break. Enf of 2002/2003 dataset aggregation.')
            break
        
        if (folder_day=='jun04'):
            print('There is no need to aggregate the data for june 04 2002 because they are already in the form of lat/lon/data. However, the vertical dimension is 1401 VS 4095 for others, so there might be an issue with this date in terms of vertical resolution!')
            continue
        
        ############################# MASTER FILE #############################
        # Load the master file of that date
        path_gps='D://OIB//2002_2003_associated_files//gps_coord//'+folder_year+'//'+folder_day+'//'+folder_day
        
        if (folder_year=='2003'):
            filename_gps=[]
            filename_gps=path_gps+'_03_gpslatlontime.mat'
            master_file= scipy.io.loadmat(filename_gps)
        elif (folder_year=='2002'):
            if(folder_day=='may18'):
                filename_gps=[]
                filename_gps=path_gps+'_02_gps_latlontime.mat'
                master_file= scipy.io.loadmat(filename_gps)
            else:
                filename_gps=[]
                filename_gps=path_gps+'_02_latlontime.mat'
                master_file= scipy.io.loadmat(filename_gps)
        
        #Create the dataframe of the master file
        df_master_file=pd.DataFrame({'seconds':pd.Series(np.ndarray.flatten(np.transpose(master_file['seconds']))),
                                     'useconds':pd.Series(np.ndarray.flatten(np.transpose(master_file['useconds']))),
                                     'lat_gps':pd.Series(np.ndarray.flatten(np.transpose(master_file['lat_gps']))),
                                     'lon_gps':pd.Series(np.ndarray.flatten(np.transpose(master_file['lon_gps']))),
                                     'time_gps':pd.Series(np.ndarray.flatten(np.transpose(master_file['time_gps']))),
                                     'seconds_gps':pd.Series(np.ndarray.flatten(np.transpose(master_file['seconds']))),
                                     'index_gps':pd.Series(np.arange(0,master_file['time_gps'].size,1)),
                                     'timearr_dec':pd.Series(np.zeros(master_file['time_gps'].size))})

        #Set the seconds column to be the index of the masterfile dataframe
        df_master_file=df_master_file.set_index('seconds')
        ############################# MASTER FILE #############################
        
        ############################## INDIVIDUAL FILES ##############################
        #Go into the daily folders 
        os.chdir(folder_day_name)
        
        #Create the quality assessment file
        filename_fquality='D://OIB//2002_2003_export//'+folder_year+'//'+folder_day+'//'+'quality_'+folder_day+'_'+folder_year+'.txt'
        f_quality = open(filename_fquality, "w")
        
        # Herebelow is a summary of what I save in this quality file
        #1. date of the file assessed
        #2. correspondance in df_final between begin 'timearr' and begin 'timearr_dec' -> 1 is matching, 0 is no match
        #3. correspondance in df_final between end 'timearr' and end 'timearr_dec' -> 1 is matching, 0 is no match
        #4. correspondance between begin 'timearr_dec' in df_final and begin 'timearr' of df_file_being_read -> 1 is matching, 0 is no match
        #5. correspondance between end 'timearr_dec' in df_final and end 'timearr' of df_file_being_read for which we have echogram data (defined by the size of filftin)-> 1 is matching, 0 is no match
        #6. length of df_final and file_being_read['filtfin'] should be the same -> 1 is yes, 0 is no
        #7. the first 3 rows of floor(df_final['timearr_dec']) and df_final['seconds_gps'] must be identical. If they are, the value 3 should be stored
        #8. the last 3 rows of floor(df_final['timearr_dec']) and df_final['seconds_gps'] must be identical. If they are, the value 3 should be stored
        #9. size of file_being_read[filtfin].shape[1]
        #10. size of df_final.shape[0]
        #pdb.set_trace()
        
        #Create the column names
        f_quality.write('date,B_match_dftimearr_dftimearrdec,E_match_dftimearr_dftimearrdec,B_match_dftimearrdec_filetimearr,E_match_dftimearrdec_filetimearr,length_match_df_file,B0to2_df_timearrdec_df_secondsgps,Em1tom3_df_timearrdec_df_secondsgps,file_being_read[filtfin].shape[1],df_final.shape[0]\n')
        #Loop over any file in the folder date and do the operations of joining in the loop
        for indiv_file in onlyfiles:
            print('Now treating the file:')
            print(join(mypath,indiv_file))
            #pdb.set_trace()
            
            if(indiv_file=='may14_03_18.mat'):
                #The echogram of this file is empty!
                print('The file '+indiv_file+' is empty! Continue ...')
                continue
            
            filename_to_check='D://OIB//2002_2003_export//'+folder_year+'//'+folder_day+'//'+indiv_file.replace(".mat","")+"_aggregated"

            if (os.path.isfile(filename_to_check)):
                #If aggregated file already exist, go to the next date.
                print(indiv_file.replace(".mat","")+"_aggregated"+' file exist, move on to the next date')
                continue

            #Load the file
            file_being_read=[]
            file_being_read=scipy.io.loadmat(join(mypath,indiv_file))
            
            #We have to round the timearr to the lower integer
            int_file_being_read = np.floor(file_being_read['timearr'])
            #Transform decimals into integer
            int_file_being_read=int_file_being_read.astype(int)
            
            #Create the dataframe
            df_file_being_read=pd.DataFrame({'timearr':pd.Series(np.ndarray.flatten(np.transpose(int_file_being_read))),
                                               'index_vector':pd.Series(np.arange(0,int_file_being_read.size,1)),
                                               'timearr_trace':pd.Series(np.ndarray.flatten(np.transpose(file_being_read['timearr']))),
                                               'timearr_floor':pd.Series(np.ndarray.flatten(np.transpose(int_file_being_read)))})
         
            #Set the timearr column to be the index of the dataframe
            df_file_being_read=df_file_being_read.set_index('timearr')
            
            #Rename the column 'timearr_trace' to 'timearr'
            df_file_being_read.columns = ['index_vector', 'timearr','timearr_floor']
            #pdb.set_trace()
            
            ###################################################################
            ### Investigate on the files having issues
            ### 1. Problem type 1:
                
            ### may11_03_40,1,1,1,1,1,3,0: size(filtfin)=(4095,940), size(timearr)=(1,1000), reset of timearr at index 940 => only data in the first 940 index!
            ### may11_03_9, 1,1,1,1,1,3,0: size(filtfin)=(4095,740), size(timearr)=(1,1000), reset of timearr at index 740 => only data in the first 740 index!
            ### may12_03_13,1,1,1,1,1,3,0: size(filtfin)=(4095,80), size(timearr)=(1,1000), reset of timearr at index 80 => only data in the first 80 index!
            ### may12_03_17,1,1,1,1,1,3,0: size(filtfin)=(4095,990), size(timearr)=(1,1000), reset of timearr at index 990 => only data in the first 990 index!
            ### may12_03_38,1,1,1,1,1,3,0: size(filtfin)=(4095,850), size(timearr)=(1,1000), reset of timearr at index 850 => only data in the first 850 index!
            ### may12_03_39,1,1,1,1,1,3,0: size(filtfin)=(4095,40), size(timearr)=(1,1000), reset of timearr at index 40 => only data in the first 40 index! Attention, another reset at index=850.
            ### may12_03_44,1,1,1,1,1,3,0: size(filtfin)=(4095,310), size(timearr)=(1,1000), reset of timearr at index 310 => only data in the first 310 index!
            ### may13_03_30,1,1,1,1,1,3,0: size(filtfin)=(4095,330), size(timearr)=(1,1000), reset of timearr at index 330 => only data in the first 330 index!
            ### may14_03_11,1,1,1,1,1,3,0: size(filtfin)=(4095,10), size(timearr)=(1,1000), reset of timearr at index 10 => only data in the first 10 index!
            ### may14_03_15,1,1,1,1,1,3,0: size(filtfin)=(4095,40), size(timearr)=(1,1000), reset of timearr at index 40 => only data in the first 40 index!
            ### may14_03_52,1,1,1,1,1,3,0: size(filtfin)=(4095,50), size(timearr)=(1,1000), reset of timearr at index 50 => only data in the first 50 index!
            ### may15_03_3, 1,1,1,1,1,3,0: size(filtfin)=(4095,330), size(timearr)=(1,1000), reset of timearr at index 330 => only data in the first 330 index!
            ### may15_03_37,1,1,1,1,1,3,0: size(filtfin)=(4095,750), size(timearr)=(1,1000), reset of timearr at index 750 => only data in the first 750 index!
            ### may18_02_30 : size(filtfin)=(4095,230), size(timearr)=(1,1000), reset of timearr at index 230 => only data in the first 230 index!
            ### may09_03_24 : size(filtfin)=(4095,320), size(timearr)=(1,1000), reset of timearr at index 320 => only data in the first 320 index!
            ### may09_03_42 : size(filtfin)=(4095,330), size(timearr)=(1,1000), reset of timearr at index 330 => only data in the first 330 index!
            ### may11_03_15 : size(filtfin)=(4095,10), size(timearr)=(1,1000), reset of timearr at index 10 => only data in the first 10 index!
            ### may11_03_18 : size(filtfin)=(4095,560), size(timearr)=(1,1000), reset of timearr at index 560 => only data in the first 560 index!
            ### -> Action needed for these files!
            
            ### 2. Problem type 2:
            
            ### may11_03_28,1,1,1,1,1,3,2: the last 'seconds' of df_file_being_read is absent in df_master_file (=jump in df_master_file <=>last row is filled with 0). This is why the last quality index=2. This date is safe and reliable!
            ### may13_03_3, 1,1,1,1,1,3,2: the last 'seconds' of df_file_being_read is absent in df_master_file (=jump in df_master_file <=>last row is filled with 0). This is why the last quality index=2. This date is safe and reliable!
            ### -> No action needed for these files, they are clean
            #pdb.set_trace()
            
            if indiv_file in list(df_issue['name_file_issue']):
                print('We are dealing with the issue file: '+indiv_file)
                #pdb.set_trace()
                
                #Identify the last index where we have data                
                row_of_interest=df_issue[df_issue['name_file_issue']==indiv_file]
                index_of_interest=np.array(row_of_interest['index_issue'])
                
                #1. Identify the value of beginign and end for selection in master file
                #Select the first and last seconds of the timearr
                value_begin=df_file_being_read['timearr_floor'].iloc[0]
                value_end=df_file_being_read['timearr_floor'].iloc[index_of_interest[0]-1] 
            
                #Sort out whether the first 'seconds' value is single or doubled
                dupli=df_file_being_read['timearr_floor']
                duplicateRowsDF = dupli.duplicated()

            else:
                #Do the classical processing steps
                #1. Identify the value of beginign and end for selection in master file
                #Select the first and last seconds of the timearr
                value_begin=df_file_being_read['timearr_floor'].iloc[0]
                value_end=df_file_being_read['timearr_floor'].iloc[-1]

                #Sort out whether the first 'seconds' value is single or doubled
                dupli=df_file_being_read['timearr_floor']
                duplicateRowsDF = dupli.duplicated()
                                
            #We need index_begin, but not index_end. However, I keep the calculation for index_end in case.
            if (duplicateRowsDF.iloc[1]):
                #if (duplicateRowsDF.iloc[1] is True, the first value is double
                index_begin=df_master_file[df_master_file['seconds_gps']==value_begin]['index_gps'].iloc[0]
            elif (~(duplicateRowsDF.iloc[1])):
                # if duplicateRowsDF.iloc[1] is False, the first value is single
                index_begin=df_master_file[df_master_file['seconds_gps']==value_begin]['index_gps'].iloc[1]
                
            ################## I am keeping that just in case #################
            #if (duplicateRowsDF.iloc[-1]):
            #    #if (duplicateRowsDF.iloc[-1] is True, the last value is double
            #    index_end=df_master_file[df_master_file['seconds_gps']==value_end]['index_gps'].iloc[-1]+1
            #elif (~(duplicateRowsDF.iloc[-1])):
            #    # if duplicateRowsDF.iloc[-1] is False, the last value is single
            #    index_end=df_master_file[df_master_file['seconds_gps']==value_end]['index_gps'].iloc[-2]+1
            ################### I am keeping that just in case #################
       
            #2. Associate the 'seconds_gps' to its decimal value from 'timearr'
            #I have to do this inside a for loop and check at any iteration because sometimes I have gaps in seconds!!
            
            #I need to know the length of the datafile having duplicates in it   
            join_doubled_duplicates=[]
            #join along the index
            join_doubled_duplicates=df_file_being_read.join(df_master_file, how='left',lsuffix='_time_arr', rsuffix='_gps_latlontime')
            
            join_duplicates=[]
            join_duplicates=join_doubled_duplicates.drop_duplicates(subset=['index_gps'])
                
            loc_df=index_begin
            i_timearr=0
            i_ite=0
        
            #1.Create the empty iterative dataset
            iterative_dataset=pd.DataFrame({'seconds':pd.Series(np.zeros(1020)),
                                            'useconds':pd.Series(np.zeros(1020)),
                                            'lat_gps':pd.Series(np.zeros(1020)),
                                            'lon_gps':pd.Series(np.zeros(1020)),
                                            'time_gps':pd.Series(np.zeros(1020)),
                                            'seconds_gps':pd.Series(np.zeros(1020)),
                                            'index_gps':pd.Series(np.zeros(1020)),
                                            'timearr_dec':pd.Series(np.zeros(1020)),
                                            'jump':pd.Series(np.zeros(1020))})
            #pdb.set_trace()
            
            #for i in range(0,len(join_duplicates),1):
            while (i_timearr<(file_being_read['filtfin'].shape[1])):
                
                if ((df_master_file['seconds_gps'].iloc[loc_df])==(df_file_being_read['timearr_floor'].iloc[i_timearr])):
                    df_master_file['timearr_dec'].iloc[loc_df]=df_file_being_read['timearr'].iloc[i_timearr]
                    
                    #Create the iterative dataset
                    iterative_dataset.iloc[i_ite]=df_master_file.iloc[loc_df]
                    iterative_dataset['timearr_dec'].iloc[i_ite]=df_file_being_read['timearr'].iloc[i_timearr]
                    iterative_dataset['jump'].iloc[i_ite]=0 #If jump, then = 1
        
                    i_ite=i_ite+1
                    i_timearr=i_timearr+1
                    loc_df=loc_df+1
                    
                elif ((df_master_file['seconds_gps'].iloc[loc_df])!=(df_file_being_read['timearr_floor'].iloc[i_timearr])):
                    #pdb.set_trace()
                    #Possibility 1: the jump is in df_master_file
                    if ((df_master_file['seconds_gps'].iloc[loc_df])>(df_file_being_read['timearr_floor'].iloc[i_timearr])):
                        #Create a empty lat/lon point which will be matched
                        #1.Create the iterative dataset
                        iterative_dataset['seconds'].iloc[i_ite]=0
                        iterative_dataset['useconds'].iloc[i_ite]=0
                        iterative_dataset['lat_gps'].iloc[i_ite]=0
                        iterative_dataset['lon_gps'].iloc[i_ite]=0
                        iterative_dataset['time_gps'].iloc[i_ite]=0
                        iterative_dataset['seconds_gps'].iloc[i_ite]=0
                        iterative_dataset['index_gps'].iloc[i_ite]=0
                        iterative_dataset['timearr_dec'].iloc[i_ite]=df_file_being_read['timearr'].iloc[i_timearr]
                        iterative_dataset['jump'].iloc[i_ite]=0 #If jump, then = 1
           
                        i_ite=i_ite+1
                        i_timearr=i_timearr+1
                        #We do not iterate the index in the df_masterfile
                        
                        #Store the information in df_file_being_read[i_time_arr] on this new line
                        #print('WE HAVE BEEN HERE!')
                        #pdb.set_trace()
                        #We have jumped one time arr, so not update of i_timearr!
                        #i_timearr=i_timearr-1
                        #df_master_file['timearr_dec'].iloc[loc_df]=df_file_being_read['timearr'].iloc[i_timearr]
                        #loc_df=loc_df+1
                        #i_timearr=i_timearr+1
                    
                    #Possibility 2: the jump is in df_file_being_read
                    elif ((df_master_file['seconds_gps'].iloc[loc_df])<(df_file_being_read['timearr_floor'].iloc[i_timearr])):
                        #Here must be calculated the mean with the previous index on the previous index
                        # Add a column where I specify this is a jump here!!
        
                        i_timearr=i_timearr-1
                        iterative_dataset.iloc[i_ite]=df_master_file.iloc[loc_df]
                        iterative_dataset['timearr_dec'].iloc[i_ite]=df_file_being_read['timearr'].iloc[i_timearr]
                        iterative_dataset['jump'].iloc[i_ite]=1 #If jump, then = 1
                        
                        i_ite=i_ite+1
                        i_timearr=i_timearr+1
                        loc_df=loc_df+1
                        #We do not iterate the index in the df_masterfile
                #print("loc_df: ", loc_df)
                #print("i_timearr: ", i_timearr)
       
            #3. Make the correspondance between timearr and seconds and join datasets
            #Set the timearr_floor column to be the index of the dataframe to avoid ambiguity in the pd.merge
            df_file_being_read=df_file_being_read.set_index('timearr_floor')
            
            merged_inner=[]
            merged_inner = pd.merge(left=df_file_being_read, right=iterative_dataset, left_on='timearr', right_on='timearr_dec')
            
            #Now that I have added the jump, I can simple remove the lines having jump=1 and then I am good!
            #Find the location where I have the jump=0 (i.e. not 1!) and obtain the dataset!
            df_final=merged_inner.loc[merged_inner['jump']==0]
            #This is done and it is working!!
            #pdb.set_trace()
            
            #Create the different index for quality assessment
            B_match_dftimearr_dftimearrdec=(df_final['timearr'].iloc[0] == df_final['timearr_dec'].iloc[0]).astype(int)
            E_match_dftimearr_dftimearrdec=(df_final['timearr'].iloc[-1] == df_final['timearr_dec'].iloc[-1]).astype(int)
            B_match_dftimearrdec_filetimearr=(df_final['timearr_dec'].iloc[0] == df_file_being_read['timearr'].iloc[0]).astype(int)
            
            if indiv_file in list(df_issue['name_file_issue']):
                E_match_dftimearrdec_filetimearr=(df_final['timearr_dec'].iloc[-1] == df_file_being_read['timearr'].iloc[index_of_interest[0]-1]).astype(int)
            else:
                E_match_dftimearrdec_filetimearr=(df_final['timearr_dec'].iloc[-1] == df_file_being_read['timearr'].iloc[-1]).astype(int)

            length_match_df_file=int(len(df_final) == (file_being_read['filtfin'].shape[1]))
        
            df_timearrdec_df_secondsgps_0=(np.floor(df_final['timearr_dec'].iloc[0]) == df_final['seconds_gps'].iloc[0]).astype(int)
            df_timearrdec_df_secondsgps_1=(np.floor(df_final['timearr_dec'].iloc[1]) == df_final['seconds_gps'].iloc[1]).astype(int)
            df_timearrdec_df_secondsgps_2=(np.floor(df_final['timearr_dec'].iloc[2]) == df_final['seconds_gps'].iloc[2]).astype(int)
            B0to2_df_timearrdec_df_secondsgps=df_timearrdec_df_secondsgps_0+df_timearrdec_df_secondsgps_1+df_timearrdec_df_secondsgps_2
            
            df_timearrdec_df_secondsgps_m1=(np.floor(df_final['timearr_dec'].iloc[-1]) == df_final['seconds_gps'].iloc[-1]).astype(int)
            df_timearrdec_df_secondsgps_m2=(np.floor(df_final['timearr_dec'].iloc[-2]) == df_final['seconds_gps'].iloc[-2]).astype(int)
            df_timearrdec_df_secondsgps_m3=(np.floor(df_final['timearr_dec'].iloc[-3]) == df_final['seconds_gps'].iloc[-3]).astype(int)
            Em1tom3_df_timearrdec_df_secondsgps=df_timearrdec_df_secondsgps_m1+df_timearrdec_df_secondsgps_m2+df_timearrdec_df_secondsgps_m3
            
            #Writting in the quality assessment file
            f_quality.write(str(indiv_file.replace(".mat",""))+','+str(B_match_dftimearr_dftimearrdec)+','+str(E_match_dftimearr_dftimearrdec)+','+str(B_match_dftimearrdec_filetimearr)+','+str(E_match_dftimearrdec_filetimearr)+','+str(length_match_df_file)+','+str(B0to2_df_timearrdec_df_secondsgps)+','+str(Em1tom3_df_timearrdec_df_secondsgps)+','+str(file_being_read['filtfin'].shape[1])+','+str(df_final.shape[0])+'\n')
        
            #Select only the variables of interest for the data storage
            df_final=df_final.drop(['seconds','timearr_dec','jump'],axis=1)
            
            #Store everything into one dictionnary (matrix and vectors of data)
            if (folder_day == 'jun04'):
                dic_file_being_read=[]
                dic_file_being_read = { "trace_id" : indiv_file.replace(".mat",""),
                     "radar_echogram" : file_being_read['data'],
                     "latlontime" : df_final}
            else:
                dic_file_being_read=[]
                dic_file_being_read = { "trace_id" : indiv_file.replace(".mat",""),
                     "radar_echogram" : file_being_read['filtfin'],
                     "latlontime" : df_final}
            
            #pdb.set_trace()
            
            #Save the dictionary into a picke file
            filename_tosave='D://OIB//2002_2003_export//'+folder_year+'//'+folder_day+'//'+indiv_file.replace(".mat","")+"_aggregated"
            outfile= open(filename_tosave, "wb" )
            pickle.dump(dic_file_being_read,outfile)
            outfile.close()
            
        f_quality.close() #Close the quality assessment file when we’re done!

##############################################################################
############################# Data manipulation ##############################
##############################################################################