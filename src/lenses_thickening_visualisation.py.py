# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 16:07:33 2021

@author: JullienN
"""
import pickle
import scipy.io
import numpy as np
import pdb
import h5py
from matplotlib import pyplot

##############################################################################
###                     Define what we want to show                        ###
##############################################################################
plot_boolean= 'plot_boolean_SG1_cut045_th350' # can be 'plot_boolean_orig_cut045_th000',
              #'plot_boolean_SG1_cut045_th000', 'plot_boolean_SG1_cut045_th350'     
plot_years_overlay='FALSE'
plot_depth_corrected_single='FALSE'
plot_depth_corrected_subplot='TRUE'
plot_boolean_subplot='TRUE'

#Define the years and data to investigate:
investigation_year={2010:['Data_20100508_01_114.mat','Data_20100508_01_115.mat'],
                    2011:['Data_20110419_01_008.mat','Data_20110419_01_009.mat','Data_20110419_01_010.mat'],
                    2012:['Data_20120418_01_129.mat','Data_20120418_01_130.mat','Data_20120418_01_131.mat'],
                    2013:['Data_20130405_01_165.mat','Data_20130405_01_166.mat','Data_20130405_01_167.mat'],
                    2014:['Data_20140424_01_002.mat','Data_20140424_01_003.mat','Data_20140424_01_004.mat'],
                    2017:['Data_20170422_01_168.mat','Data_20170422_01_169.mat','Data_20170422_01_170.mat','Data_20170422_01_171.mat'],
                    2018:['Data_20180427_01_170.mat','Data_20180427_01_171.mat','Data_20180427_01_172.mat']}

#investigation_year={2010:['Data_20100513_01_001.mat','Data_20100513_01_002.mat'],
#                    2011:['Data_20110411_01_116.mat','Data_20110411_01_117.mat','Data_20110411_01_118.mat'],
#                    2012:['Data_20120428_01_125.mat','Data_20120428_01_126.mat'],
#                    2013:'empty',
#                    2014:['Data_20140408_11_024.mat','Data_20140408_11_025.mat','Data_20140408_11_026.mat'],
#                    2017:['Data_20170508_02_168.mat','Data_20170508_02_169.mat','Data_20170508_02_170.mat','Data_20170508_02_171.mat']}


#investigation_year={2010:['Data_20100507_01_008.mat','Data_20100507_01_009.mat','Data_20100507_01_010.mat'],
#                    2011:['Data_20110426_01_009.mat','Data_20110426_01_010.mat','Data_20110426_01_011.mat'],
#                    2012:'empty',
#                    2013:'empty',
#                    2014:['Data_20140421_01_009.mat','Data_20140421_01_010.mat','Data_20140421_01_011.mat','Data_20140421_01_012.mat','Data_20140421_01_013.mat'],
#                    2017:['Data_20170424_01_008.mat','Data_20170424_01_009.mat','Data_20170424_01_010.mat','Data_20170424_01_011.mat','Data_20170424_01_012.mat']}

##############################################################################
###                     Define what we want to show                        ###
##############################################################################

##############################################################################
###                               Define paths                             ###
##############################################################################
#Define the general path as a function of the year
path_data='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/'
#Define the path for masks
path_mask=path_data+'exported/temp_for_overlap/Boolean Array Picklefiles/'
#Define path for depth corrected
path_depth_corrected=path_data+'exported/temp_for_overlap/Depth_Corrected_Picklefiles/'
#Define the path for boolean
path_boolean=path_data+'exported/temp_for_overlap/Boolean Array Picklefiles/'
##############################################################################
###                               Define paths                             ###
##############################################################################


##############################################################################
###                          Load and organise data                        ###
##############################################################################
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
    
    #pdb.set_trace()
    
    #Define filename depth corrected data
    filename_depth_corrected=date_track+'_DEPTH_CORRECTED.pickle'
              
    #Define boolean filename
    if (plot_boolean=='plot_boolean_orig_cut045_th000'):
        filename_boolean=date_track+'_orig_CUTOFF_-0.45_THRESHOLD_000.pickle'
        
    if (plot_boolean=='plot_boolean_SG1_cut045_th000'):
        filename_boolean=date_track+'_SG1_CUTOFF_-0.45_THRESHOLD_000.pickle'
        
    if (plot_boolean=='plot_boolean_SG1_cut045_th350'):
        filename_boolean=date_track+'_SG1_CUTOFF_-0.45_THRESHOLD_350.pickle'
    
    #Open the depth corrected file
    f_depth_corrected = open(path_depth_corrected+filename_depth_corrected, "rb")
    radar = pickle.load(f_depth_corrected)
    f_depth_corrected.close()
    
    #Open boolean file
    f_boolean = open(path_boolean+filename_boolean, "rb")
    boolean_file = pickle.load(f_boolean)
    f_boolean.close()
    
    #Open mask file
    f_mask = open(path_mask+date_track+'_mask.pickle', "rb")
    data_mask = pickle.load(f_mask)
    f_mask.close()
    
    #Create the title for the figures
    file_for_title=filename_depth_corrected
    file_for_title=file_for_title.partition("_")[2]
    file_for_title=file_for_title.partition("_")[2]
    file_for_title=file_for_title.partition("_")[2]
    file_for_title=file_for_title.partition("_")[2]
    file_for_title=file_for_title.partition(".pickle")[0]
    
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
            
        else:
            fdata_filename = scipy.io.loadmat(path_raw_data+indiv_file_load)
            lat_filename = fdata_filename['Latitude']
            lon_filename = fdata_filename['Longitude']
        
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
        radar=np.fliplr(radar)
        data_mask=np.flipud(data_mask)
    
    #Store reunited lat/lon, slice output and mask in a dictionnary:
    dataframe[str(single_year)]={'lat_appended':lat_appended,
                                 'lon_appended':lon_appended,
                                 'radar':radar,
                                 'boolean':boolean_file,
                                 'mask':data_mask}
    
##############################################################################
###                          Load and organise data                        ###
##############################################################################

#Plot the results:
#1. Create a plot with the minimum and maximum extent of the traces for all year
#2. Overlay (x% transparent the boolean traces on top of each other)    
   
##############################################################################
###                 Identifiy the lower and upper bound                    ###
##############################################################################
 
#Find the min and max longitudes present in dataframe:
min_lon=0
max_lon=-180

for single_year in investigation_year.keys():
    
    #If no data, continue
    if (investigation_year[single_year]=='empty'):
        print('No data for year '+str(single_year)+', continue')
        continue
    
    start_date_track=investigation_year[single_year][0]
    end_date_track=investigation_year[single_year][-1]
    date_track=start_date_track[5:20]+'_'+end_date_track[17:20]
    
    #Calculate the min longitude
    min_lon_temp=np.ndarray.min(dataframe[str(single_year)]['lon_appended'])
    if(min_lon_temp<min_lon):
        min_lon=min_lon_temp
        #print('Min is:'+str(min_lon))
    
    #Calculate the max longitude
    max_lon_temp=np.ndarray.max(dataframe[str(single_year)]['lon_appended'])
    if(max_lon_temp>max_lon):
        max_lon=max_lon_temp
        #print('Max is:'+str(max_lon))

##############################################################################
###                 Identifiy the lower and upper bound                    ###
##############################################################################

##############################################################################
###                                 Plot data                              ###
##############################################################################
 
pdb.set_trace()

if (plot_years_overlay=='TRUE'):
    #Create an empty radar slice to plot data over it
    empty_slice=np.empty((dataframe[str(single_year)]['radar'].shape[0],5000))
    #Plot the data:
    
    pyplot.figure(figsize=(48,40))
    #Change label font
    pyplot.rcParams.update({'font.size': 40})
    color_map=pyplot.pcolor(empty_slice,cmap=pyplot.get_cmap('gray'))#,norm=divnorm)
    pyplot.gca().invert_yaxis() #Invert the y axis = avoid using flipud.
    #pyplot.yticks(ticks=ticks_yplot,labels=(np.round(depths[ticks_yplot])))
    pyplot.ylabel('Depth [m]')
    pyplot.xlabel('Horizontal distance - lon [°]')
    #pyplot.title('Raw radar echogram: '+indiv_file.replace("_aggregated",""))
    cbar=pyplot.colorbar()
    cbar.set_label('Signal strength')


#Prepare the plot for all years display
if (plot_depth_corrected_subplot=='TRUE'):
    pyplot.rcParams['axes.linewidth'] = 0.1 #set the value globally
    pyplot.rcParams['xtick.major.width']=0.1
    pyplot.rcParams['ytick.major.width']=0.1
    
    pyplot.figure(figsize=(40,20))
    pyplot.rcParams.update({'font.size': 2})
    fig1, (ax1s,ax2s,ax3s,ax4s,ax5s,ax6s,ax7s) = pyplot.subplots(7, 1)
    
if (plot_boolean_subplot=='TRUE'):
    pyplot.rcParams['axes.linewidth'] = 0.1 #set the value globally
    pyplot.rcParams['xtick.major.width']=0.1
    pyplot.rcParams['ytick.major.width']=0.1
    
    pyplot.figure(figsize=(40,20))
    pyplot.rcParams.update({'font.size': 2})
    fig2, (ax1b,ax2b,ax3b,ax4b,ax5b,ax6b,ax7b) = pyplot.subplots(7, 1)
        
for single_year in investigation_year.keys():
    
    print(str(single_year))
    
    #If no data, continue
    if (investigation_year[single_year]=='empty'):
        print('No data for year '+str(single_year)+', continue')
        continue
    
    start_date_track=investigation_year[single_year][0]
    end_date_track=investigation_year[single_year][-1]
    date_track=start_date_track[5:20]+'_'+end_date_track[17:20]
    
    if (plot_depth_corrected_single=='TRUE'):
        pyplot.rcParams['axes.linewidth'] = 0.1 #set the value globally
        pyplot.rcParams['xtick.major.width']=0.1
        pyplot.rcParams['ytick.major.width']=0.1
        
        pyplot.figure(figsize=(40,20))
        pyplot.rcParams.update({'font.size': 2})
        fig, (ax1) = pyplot.subplots(1, 1)
        
        #fig.suptitle(str(plot_name1))
        X=dataframe[str(single_year)]['lon_appended']
        Y=np.arange(0,100,100/dataframe[str(single_year)]['radar'].shape[0])
        C=dataframe[str(single_year)]['radar'].astype(float)
                
        cb1=ax1.pcolor(X, Y, C,cmap=pyplot.get_cmap('gray'))#,norm=divnorm)
        ax1.invert_yaxis() #Invert the y axis = avoid using flipud.
        ax1.set_aspect(0.0025) # X scale matches Y scale
        ax1.set_title(date_track+' Depth corrected')
        ax1.set_ylabel('Depth [m]')
        ax1.set_xlabel('Longitude [°]')
        #pdb.set_trace()
        #ax1.set_xlim(-47.9,-46.8)
        ax1.set_xlim(min_lon,max_lon)
        ax1.set_ylim(30,0)
        
        cbar1=fig.colorbar(cb1, ax=[ax1], location='right',shrink=0.12,aspect=10,pad=0.01)
        cbar1.set_label('Signal strength')

        ##Create the figure name
        #fig_name=[]
        #fig_name='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/2010_2014_thickening/'+date_track+'_depthcorrected.png'
        
        pyplot.show()
        
        ##Save the figure
        #pyplot.savefig(fig_name,dpi=2000)
        #pyplot.clf()

    if (plot_depth_corrected_subplot=='TRUE'):
        if (single_year==2010):
            ax_plotting=ax1s
        elif (single_year==2011):
            ax_plotting=ax2s
        elif (single_year==2012):
            ax_plotting=ax3s
        elif (single_year==2013):
            ax_plotting=ax4s
        elif (single_year==2014):
            ax_plotting=ax5s
        elif (single_year==2017):
            ax_plotting=ax6s
        elif (single_year==2018):
            ax_plotting=ax7s
        else:
            print('Year not existing')
            
        #fig.suptitle(str(plot_name1))
        X=dataframe[str(single_year)]['lon_appended']
        Y=np.arange(0,100,100/dataframe[str(single_year)]['radar'].shape[0])
        C=dataframe[str(single_year)]['radar'].astype(float)
                
        cb=ax_plotting.pcolor(X, Y, C,cmap=pyplot.get_cmap('gray'))#,norm=divnorm)
        ax_plotting.invert_yaxis() #Invert the y axis = avoid using flipud.
        ax_plotting.set_aspect(0.0025) # X scale matches Y scale
        ax_plotting.set_title(str(single_year)+' Depth corrected')
        ax_plotting.set_ylabel('Depth [m]')
        ax_plotting.set_xlabel('Longitude [°]')
        ax_plotting.set_xlim(min_lon,max_lon)
        ax_plotting.set_ylim(30,0)
        
        cbar=fig1.colorbar(cb, ax=[ax_plotting], location='right',shrink=0.12,aspect=10,pad=0.01)
        cbar.set_label('Signal strength')

        ##Create the figure name
        #fig_name=[]
        #fig_name='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/2010_2014_thickening/'+date_track+'_depth_corrected.png'
        
        pyplot.show()
        
        ##Save the figure
        #pyplot.savefig(fig_name,dpi=2000)
        #pyplot.clf()

    if (plot_boolean_subplot=='TRUE'):
        if (single_year==2010):
            ax_plotting=ax1b
        elif (single_year==2011):
            ax_plotting=ax2b
        elif (single_year==2012):
            ax_plotting=ax3b
        elif (single_year==2013):
            ax_plotting=ax4b
        elif (single_year==2014):
            ax_plotting=ax5b
        elif (single_year==2017):
            ax_plotting=ax6b
        elif (single_year==2018):
            ax_plotting=ax7b
        else:
            print('Year not existing')
            
        #fig.suptitle(str(plot_name1))
        X=dataframe[str(single_year)]['lon_appended']
        Y=np.arange(0,100,100/dataframe[str(single_year)]['boolean'].shape[0])
        C=dataframe[str(single_year)]['boolean'].astype(float)
        
        cb=ax_plotting.pcolor(X, Y, C,cmap=pyplot.get_cmap('gray'))#,norm=divnorm)
        ax_plotting.invert_yaxis() #Invert the y axis = avoid using flipud.
        ax_plotting.set_aspect(0.0025) # X scale matches Y scale
        ax_plotting.set_title(str(single_year)+' '+plot_boolean)
        ax_plotting.set_ylabel('Depth [m]')
        ax_plotting.set_xlabel('Longitude [°]')
        ax_plotting.set_xlim(min_lon,max_lon)
        ax_plotting.set_ylim(30,0)
        
        cbar=fig1.colorbar(cb, ax=[ax_plotting], location='right',shrink=0.12,aspect=10,pad=0.01)
        cbar.set_label('Signal strength')
        
        ##Create the figure name
        #fig_name=[]
        #fig_name='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/2010_2014_thickening/'+date_track+'_'+plot_boolean+'.png'
        
        pyplot.show()
        
        ##Save the figure
        #pyplot.savefig(fig_name,dpi=2000)
        #pyplot.clf()        
    else:
        
        pyplot.rcParams['axes.linewidth'] = 0.1 #set the value globally
        pyplot.rcParams['xtick.major.width']=0.1
        pyplot.rcParams['ytick.major.width']=0.1
        
        pyplot.figure(figsize=(40,20))
        pyplot.rcParams.update({'font.size': 3})
        fig, (ax1) = pyplot.subplots(1, 1)
    
        #fig.suptitle(str(plot_name1))
        X=dataframe[str(single_year)]['lon_appended']
        Y=np.arange(0,20,20/dataframe[str(single_year)]['radar'].shape[0])
        C=dataframe[str(single_year)]['radar'].astype(float)
        C[C==0]=np.nan
        
        cb1=ax1.pcolor(X, Y, C,cmap=pyplot.get_cmap('gray'))#,alpha=0.1,edgecolor='none')#,norm=divnorm)
        ax1.invert_yaxis() #Invert the y axis = avoid using flipud.
        ax1.set_aspect(0.0025) # X scale matches Y scale
        ax1.set_title(date_track+' '+plot_boolean)
        ax1.set_ylabel('Depth [m]')
        ax1.set_xlabel('Longitude [°]')        
        #ax1.set_xlim(-47.9,-46.8)
        ax1.set_xlim(min_lon,max_lon)
        #ax1.set_ylim(30,0)
        
        cbar1=fig.colorbar(cb1, ax=[ax1], location='right',shrink=0.12,aspect=10,pad=0.01)
        cbar1.set_label('Signal strength')
            
        pyplot.show()
        
        ##Create the figure name
        #fig_name=[]
        #fig_name='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/2010_2014_thickening/'+date_track+'_'+file_for_title+'.png'
    
        ##Save the figure
        #pyplot.savefig(fig_name,dpi=2000)
        #pyplot.clf()
        
        pdb.set_trace()

##############################################################################
###                                 Plot data                              ###
##############################################################################


pdb.set_trace()

#Replace the FALSE by NaN, then plot different colors? or semi transparent?

single_year=2011
X=dataframe[str(single_year)]['lon_appended']
Y=np.arange(0,20,20/dataframe[str(single_year)]['radar'].shape[0])
C=dataframe[str(single_year)]['radar'].astype(float)
C[C==0]=np.nan
color_map=pyplot.pcolor(X, Y, C,cmap=pyplot.get_cmap('gray'),alpha=0.1)#,norm=divnorm)
pdb.set_trace()

single_year=2012
X=dataframe[str(single_year)]['lon_appended']
Y=np.arange(0,20,20/dataframe[str(single_year)]['radar'].shape[0])
C=dataframe[str(single_year)]['radar'].astype(float)
C[C==0]=np.nan
color_map=pyplot.pcolor(X, Y, C,cmap=pyplot.get_cmap('summer'),alpha=0.1)#,norm=divnorm)
pdb.set_trace()

single_year=2013
X=dataframe[str(single_year)]['lon_appended']
Y=np.arange(0,20,20/dataframe[str(single_year)]['radar'].shape[0])
C=dataframe[str(single_year)]['radar'].astype(float)
C[C==0]=np.nan
color_map=pyplot.pcolor(X, Y, C,cmap=pyplot.get_cmap('autumn'),alpha=0.1)#,norm=divnorm)
pdb.set_trace()

single_year=2014
X=dataframe[str(single_year)]['lon_appended']
Y=np.arange(0,20,20/dataframe[str(single_year)]['radar'].shape[0])
C=dataframe[str(single_year)]['radar'].astype(float)
C[C==0]=np.nan
color_map=pyplot.pcolor(X, Y, C,cmap=pyplot.get_cmap('winter'),alpha=0.1)#,norm=divnorm)

pyplot.show()
pdb.set_trace()

##Create the figure name
#fig_name=[]
#fig_name='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/2002_2003_radar_raw_echogram/'+indiv_file+'.png'

##Save the figure
#pyplot.savefig(fig_name)
#pyplot.clf()
                    
#pdb.set_trace()
    
    
    
#Work with the booelan and loop over the year and add every year
#on top of each opther, then plot the intensity
    
    
    
        
#lat1.shape[1]+lat2.shape[1]+lat3.shape[1]= data.shape[1] => great news!

#To do:
#    1. append data and lat/lon1 OK
#    2. load data from 2010 to 2014 OK
#    3. select as a function of lat/lon box OK
#    4. plot the data OK
#    5. compare! OK
#    6. Decide whether I should work with depth corrected, or other post-precessed files that are available in 'Boolean Array Picklefiles' folder OK