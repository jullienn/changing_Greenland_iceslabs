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
        FV = 'double'
        index_begin=df_master_file[df_master_file['seconds_gps']==value_begin]['index_gps'].iloc[0]
        #index_end=df_master_file[df_master_file['seconds_gps']==value_end]['index_gps'].iloc[-2]
    elif (~(duplicateRowsDF.iloc[1])):
        # if duplicateRowsDF.iloc[1] is False, the first value is single
        FV = 'single'
        index_begin=df_master_file[df_master_file['seconds_gps']==value_begin]['index_gps'].iloc[1]
        #index_end=df_master_file[df_master_file['seconds_gps']==value_end]['index_gps'].iloc[-1]
    
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
    
    pdb.set_trace()
    
    for i in range(0,len(join_duplicates),1):
        
        if (i_timearr>=len(df_file_being_read)):
            print('break out')
            break
        
        #len(df_file_being_read)
        if ((df_master_file['seconds_gps'].iloc[loc_df])==(df_file_being_read['timearr_floor'].iloc[i_timearr])):
            df_master_file['timearr_dec'].iloc[loc_df]=df_file_being_read['timearr'].iloc[i_timearr]
            i_timearr=i_timearr+1
            loc_df=loc_df+1
            
        elif ((df_master_file['seconds_gps'].iloc[loc_df])!=(df_file_being_read['timearr_floor'].iloc[i_timearr])):
            #We have jumped one time arr, so not update of i_timearr!
            
            i_timearr=i_timearr-1
            #df_master_file['timearr_dec'].iloc[loc_df]=df_file_being_read['timearr'].iloc[i_timearr]
            loc_df=loc_df+1
            i_timearr=i_timearr+1
        
        print("i_timearr: ", i_timearr)
        #print("loc_df: ", loc_df)



    ##### THERE IS A PROBLEM WITH THE DATE 'may13_03_11.mat' !!! SOLVE THE ISSUE
    #For the jump in seconds, do compute the mean of lat, lon, time_gps between the 2 to end up with one data point
    pdb.set_trace()
    
    #3. Make the correspondance between timearr and seconds and join datasets
    #Set the timearr_floor column to be the index of the dataframe to avoid ambiguity in the pd.merge
    df_file_being_read=df_file_being_read.set_index('timearr_floor')
    
    merged_inner=[]
    merged_inner = pd.merge(left=df_file_being_read, right=df_master_file, left_on='timearr', right_on='timearr_dec')
    
    
    
    
    result_join=[]
    
    #join along the index
    result_join=df_file_being_read.join(df_master_file, how='left',lsuffix='_time_arr', rsuffix='_gps_latlontime')
    
    # Be careful about the remove duplicates with respect to the time_gps. This
    # is the only solution I found to indeed remove duplicates but maybe check
    # on other datasets to be sure that it does the right thing!
       
    result_join_without_duplicates=[]
    result_join_without_duplicates=result_join.drop_duplicates(subset=['index_gps'])
    
    #Store everything into one dictionnary (matrix and vectors of data)
    dic_file_being_read=[]
    dic_file_being_read = { "trace_id" : indiv_file.replace(".mat",""),
         "radar_echogram" : file_being_read['filtfin'],
         "latlontime" : result_join_without_duplicates }
    pdb.set_trace()

        
#        ith_seconds_gps=seconds_gps_stored
#        if (count<1):
#            #This is the first time we have this seconds_gps
#            df_master_file['timearr_dec']=df_file_being_read['timearr_floor'].iloc[i]
#
#            count=count+1
#        else:
#            #This is the second time we have this seconds_gps
#            
#            count=0
#            
#    
#    master_df_of_interest['seconds_gps_decimal'].iloc[i]=
#
#    #Create the sliced gps dataset
#    master_df_of_interest=df_master_file.iloc[index_begin:index_end]

        
        
        


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