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


############################ Download old AR data #############################
#Code from: https://gist.github.com/nasrulhazim/cfd5f01e3b261b09d54f721cc1a7c50d

from ftplib import FTP
from datetime import datetime


start = datetime.now()
ftp = FTP('data.cresis.ku.edu')
ftp.login()

#path='/data/accum/old_format/2003/'
path='/data/accum/old_format/'
ftp.cwd(path)

# Get folders_years name
folders_years = ftp.nlst()

for folder_year in folders_years:
    
    if (folder_year == '2006'):
        print('Break because year 2006')
        break
    elif (folder_year == '2009'):
        print('Break because year 2009')
        break
    elif (folder_year == '2011'):
        print('Break because year 2011')
        break
    else:
        print('This is either the year 2002 or 2003')
        folder_year_name=[]
        folder_year_name=path + folder_year + '/'
        ftp.cwd(folder_year_name)
        
        # For this particular year, get folders name
        folders=[]
        folders = ftp.nlst()
        #Go to the folder
        
        for folder in folders:
            folder_name=[]
            folder_name=path + folder_year + '/' + folder + '/'
            ftp.cwd(folder_name)
            print("Now in folder ..." + folder_name)
        
            # Get All Files in this folder
            files=[]
            files = ftp.nlst()
            
            # Print out the files
            for file in files:
                print("Downloading..." + file)
                ftp.retrbinary("RETR " + file ,open("D:/OIB/AR/" + folder_year + '/' + folder + "/" + file, 'wb').write)
                #ftp.retrbinary("RETR " + file ,open("download/to/your/directory/" + file, 'wb').write)

ftp.close()
end = datetime.now()
diff = end - start
print('All files downloaded for ' + str(diff.seconds) + 's')
############################ Download old AR data #############################



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
'seconds_gps':pd.Series(np.ndarray.flatten(np.transpose(master_file['seconds'])))})
    
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
                                       'timearr_trace':pd.Series(np.ndarray.flatten(np.transpose(file_being_read['timearr'])))})
 
    #Set the timearr column to be the index of the dataframe
    df_file_being_read=df_file_being_read.set_index('timearr')
    
    #Rename the column 'timearr_trace' to 'timearr'
    df_file_being_read.columns = ['index_vector', 'timearr']
    pdb.set_trace()
    #Make the correspondance between timearr and seconds and join datasets
    result_join=[]
    
    #join along the index
    result_join=df_file_being_read.join(df_master_file, how='left',lsuffix='_time_arr', rsuffix='_gps_latlontime')
    
    # Be careful about the remove duplicates with respect to the time_gps. This
    # is the only solution I found to indeed remove duplicates but maybe check
    # on other datasets to be sure that it does the right thing!
    result_join_without_duplicates=[]
    result_join_without_duplicates=result_join.drop_duplicates(subset=['time_gps'])
    
    #Store everything into one dictionnary (matrix and vectors of data)
    dic_file_being_read=[]
    dic_file_being_read = { "trace_id" : indiv_file.replace(".mat",""),
         "radar_echogram" : file_being_read['filtfin'],
         "latlontime" : result_join_without_duplicates }
    
    pdb.set_trace()
    
    
    #there is an issue with the matching between the time! Some final match have size of 1001, other 1002, other 1000 (which should be the case for all)
    
    
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