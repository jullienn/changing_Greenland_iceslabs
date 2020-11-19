# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 15:35:11 2020

@author: Nicolas Jullien
"""
import scipy.io
import rasterio
from matplotlib import pyplot
import numpy as np
import h5py
import matplotlib.colors as mcolors
import pandas as pd
from os import listdir
from os.path import isfile, join
import pdb


############################# Download old AR data #############################
##Code from: https://gist.github.com/nasrulhazim/cfd5f01e3b261b09d54f721cc1a7c50d
#
#from ftplib import FTP
#from datetime import datetime
#
#
#start = datetime.now()
#ftp = FTP('data.cresis.ku.edu')
#ftp.login()
#
##path='/data/accum/old_format/2003/'
#path='/data/accum/old_format/'
#ftp.cwd(path)
#
## Get folders_years name
#folders_years = ftp.nlst()
#
#for folder_year in folders_years:
#    
#    if (folder_year == '2006'):
#        print('Break because year 2006')
#        break
#    elif (folder_year == '2009'):
#        print('Break because year 2009')
#        break
#    elif (folder_year == '2011'):
#        print('Break because year 2011')
#        break
#    else:
#        print('This is either the year 2002 or 2003')
#        folder_year_name=[]
#        folder_year_name=path + folder_year + '/'
#        ftp.cwd(folder_year_name)
#        
#        # For this particular year, get folders name
#        folders=[]
#        folders = ftp.nlst()
#        #Go to the folder
#        
#        for folder in folders:
#            folder_name=[]
#            folder_name=path + folder_year + '/' + folder + '/'
#            ftp.cwd(folder_name)
#            print("Now in folder ..." + folder_name)
#        
#            # Get All Files in this folder
#            files=[]
#            files = ftp.nlst()
#            
#            # Print out the files
#            for file in files:
#                print("Downloading..." + file)
#                ftp.retrbinary("RETR " + file ,open("D:/OIB/AR/" + folder_year + '/' + folder + "/" + file, 'wb').write)
#                #ftp.retrbinary("RETR " + file ,open("download/to/your/directory/" + file, 'wb').write)
#
#ftp.close()
#end = datetime.now()
#diff = end - start
#print('All files downloaded for ' + str(diff.seconds) + 's')
############################# Download old AR data #############################



##############################################################################
############################# Data manipulation ##############################
##############################################################################
# 1. Let's start with 2003 dataset

################################ MASTER FILE ################################
# Start by loading the master file of that date
master_file= scipy.io.loadmat('C://Users//Nicolas Jullien//Documents//PhD//iceslabs_processing//iceslabs_MacFerrin//data//2003//may13_03_gpslatlontime.mat')
master_file['seconds']
master_file['useconds']
master_file['lat_gps']
master_file['lon_gps']
master_file['time_gps']
    
#Create the dataframe  of the master file       
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

################################ MASTER FILE ################################

############################## INDIVIDUAL FILES ##############################

#Define the path
mypath = 'C://Users//Nicolas Jullien//Documents//PhD//iceslabs_processing//iceslabs_MacFerrin//data//2003//may13//'

#Save the filenames present in folder of interest (here May 13 2003)
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
pdb.set_trace()
#Loop over any file in the folder date and do the operations of joining in the loop
for indiv_file in onlyfiles:
    print(join(mypath,indiv_file))
    
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
    pdb.set_trace()

    #1. select the begining of the df_master_file of concern
    #Select the first and last seconds of the timearr
    value_begin=df_file_being_read['timearr_floor'].iloc[0]
    value_end=df_file_being_read['timearr_floor'].iloc[-1]
    
    #Sort out whether the first 'seconds' value is single or doubled
    dupli=df_file_being_read['timearr_floor']
    duplicateRowsDF = dupli.duplicated()

    if (duplicateRowsDF.iloc[1]):
        #if (duplicateRowsDF.iloc[1] is True, the first value is double
        index_begin=df_master_file[df_master_file['seconds_gps']==value_begin]['index_gps'].iloc[0]
    elif (~(duplicateRowsDF.iloc[1])):
        # if duplicateRowsDF.iloc[1] is False, the first value is single
        index_begin=df_master_file[df_master_file['seconds_gps']==value_begin]['index_gps'].iloc[1]
        
    if (duplicateRowsDF.iloc[-1]):
        #if (duplicateRowsDF.iloc[-1] is True, the last value is double
        index_end=df_master_file[df_master_file['seconds_gps']==value_end]['index_gps'].iloc[-1]+1
    elif (~(duplicateRowsDF.iloc[-1])):
        # if duplicateRowsDF.iloc[-1] is False, the last value is single
        index_end=df_master_file[df_master_file['seconds_gps']==value_end]['index_gps'].iloc[-2]+1
    
    #Select the corresponding slice of df_master_file:
    df_slice=df_master_file.iloc[index_begin:index_end]
    
    #Check if we have duplicates of missing values in this slice df compared to the file being read
    merged_check=[]
    merged_check = pd.merge(left=df_file_being_read, right=df_slice, left_on='timearr_floor', right_on='seconds_gps')
    
    check_join_duplicates=[]
    check_join_duplicates=merged_check.drop_duplicates(subset=['time_gps'],keep='first')
    
    pdb.set_trace()
    
    #2. Associate the 'seconds_gps' to its decimal value from 'timearr'
    #I have to do this inside a for loop and check at any iteration because sometimes I have gaps in seconds!!
    #count=0
    #seconds_gps_stored=df_master_file['seconds_gps'].iloc[index_begin]
    
    
    
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
    iterative_dataset=pd.DataFrame({'seconds':pd.Series(np.zeros(1010)),
                                    'useconds':pd.Series(np.zeros(1010)),
                                    'lat_gps':pd.Series(np.zeros(1010)),
                                    'lon_gps':pd.Series(np.zeros(1010)),
                                    'time_gps':pd.Series(np.zeros(1010)),
                                    'seconds_gps':pd.Series(np.zeros(1010)),
                                    'index_gps':pd.Series(np.zeros(1010)),
                                    'timearr_dec':pd.Series(np.zeros(1010)),
                                    'jump':pd.Series(np.zeros(1010))})
    pdb.set_trace()
    
    #for i in range(0,len(join_duplicates),1):
    while (i_timearr<len(df_file_being_read)):

        #if (i_timearr>=len(df_file_being_read)):
        #    print('break out')
        #    break
        
        #len(df_file_being_read)
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
            pdb.set_trace()
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
            if ((df_master_file['seconds_gps'].iloc[loc_df])<(df_file_being_read['timearr_floor'].iloc[i_timearr])):
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
        
        print("i_timearr: ", i_timearr)
        #print("loc_df: ", loc_df)
                
    #3. Make the correspondance between timearr and seconds and join datasets
    #Set the timearr_floor column to be the index of the dataframe to avoid ambiguity in the pd.merge
    df_file_being_read=df_file_being_read.set_index('timearr_floor')
    
    merged_inner=[]
    merged_inner = pd.merge(left=df_file_being_read, right=iterative_dataset, left_on='timearr', right_on='timearr_dec')
    
    #Now that I have added the jump, I can simple remove the lines having jump=1 and then I am good!
    #Find the location where I have the jump=0 (i.e. not 1!) and obtain the dataset!
    df_final=merged_inner.loc[merged_inner['jump']==0]
    #This is done and it is working!!
    pdb.set_trace()
  
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
    #Add a quality accessment: first and last timearr in df_final must fullfil
    #the first and last from timearr that are present in the df_file_begin_read.
    #Also, make sure the length between df_final and df_file_being_read are identical
    #Also, check that the length of merged_inner=2*1000-(nb of jumps)
    #Should I check also the number of times I do not have data?
    #Save the results in a separate text file.
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
    
    #Select only the variables of interest for the data storage
    df_for_storage=
    #Store everything into one dictionnary (matrix and vectors of data)
    dic_file_being_read=[]
    dic_file_being_read = { "trace_id" : indiv_file.replace(".mat",""),
         "radar_echogram" : file_being_read['filtfin'],
         "latlontime" : result_join_without_duplicates }
    pdb.set_trace()

        

#        #initialize the count
#        count=0
        
        
#        df_master_file['seconds_gps'][df_master_file['seconds_gps']==value_begin]
        
#        count=count+1
    
    
#        df_master_file[df_master_file['seconds_gps']==value_begin].index.item()
#        int(df_master_file[df_master_file['seconds_gps']==value_begin]['index_gps'].index[0])
    
#    int(df[df['A']==5].index[0])
##############################################################################
############################# Data manipulation ##############################
##############################################################################    
    
    
#pyplot.figure()
#pyplot.plot((np.arange(0,len(np.array(result_join_without_duplicates['index_vector'])),1)),np.array(result_join_without_duplicates['index_vector']))
#pyplot.show()
    
    
#duplicateRowsDF = result_join[result_join.duplicated()]
#print("Duplicate Rows except first occurrence based on all columns are :")
#print(duplicateRowsDF)

#duplicateRowsDF = result_join_without_duplicates[result_join_without_duplicates.duplicated(['index_vector'])]
#print("Duplicate Rows based on a single column are:", duplicateRowsDF, sep='\n')

#double_duplicateRowsDF = duplicateRowsDF[duplicateRowsDF.duplicated(['index_vector'])]
#print("Duplicate Rows based on a single column are:", duplicateRowsDF, sep='\n')

#pyplot.Annotation(
#    np.array(result_join_without_duplicates['index_vector']), xy, kwargs)