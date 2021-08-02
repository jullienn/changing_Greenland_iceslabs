# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 16:07:33 2021

@author: JullienN
"""

##############################################################################
############### Define function for discrete colorbar display ###############
##############################################################################
def discrete_cmap(N, base_cmap=None):
    """Create an N-bin discrete colormap from the specified input map"""
    #This piece of code is from: https://gist.github.com/jakevdp/91077b0cae40f8f8244a
    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:
    
    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)
##############################################################################
############### Define function for discrete colorbar display ###############
##############################################################################


def compute_distances(lon,lat):
    
    #Transform the coordinated from WGS84 to EPSG:3413
    #Example from: https://pyproj4.github.io/pyproj/stable/examples.html
    transformer = Transformer.from_crs("EPSG:4326", "EPSG:3413", always_xy=True)
    points=transformer.transform(np.array(lon),np.array(lat))
    
    #Reset the lat_3413 and lon_3413 to empty vectors.
    lon_3413=[]
    lat_3413=[]
    
    lon_3413=points[0]
    lat_3413=points[1]


    #This part of code is from MacFerrin et al., 2019
    '''Compute the distance (in m here, not km as written originally) of the traces in the file.'''
    # C = sqrt(A^2  + B^2)
    distances = np.power(np.power((lon_3413[1:] - lon_3413[:-1]),2) + np.power((lat_3413[1:] - lat_3413[:-1]),2), 0.5)
    #Calculate the cumsum of the distances
    cumsum_distances=np.nancumsum(distances)
    #Seeting the first value of the cumsum to be zero as it is the origin
    return_cumsum_distances=np.zeros(lon_3413.shape[0])
    return_cumsum_distances[1:lon_3413.shape[0]]=cumsum_distances
    
    return return_cumsum_distances

import pickle
import scipy.io
import numpy as np
import pdb
import h5py
from matplotlib import pyplot
from PIL import Image
from pyproj import Transformer
import xarray as xr

##############################################################################
###                     Define what we want to show                        ###
##############################################################################
plot_boolean= 'plot_boolean_SG1_cut045_th350' # can be 'plot_boolean_orig_cut045_th000',
              #'plot_boolean_SG1_cut045_th000', 'plot_boolean_SG1_cut045_th350'     
plot_years_overlay='FALSE'
plot_depth_corrected_single='FALSE'
plot_depth_corrected_subplot='TRUE'
plot_boolean_subplot='TRUE'
plot_images_subplot='FALSE'
yearly_comparison_indiv='FALSE'
yearly_comparison_ref='TRUE'
cumulative_comparison='FALSE'

#Define the years and data to investigate:
 

#Less good candidate for ice slabs filling!!
investigation_year={2010:['Data_20100508_01_114.mat','Data_20100508_01_115.mat'],
                    2011:['Data_20110419_01_008.mat','Data_20110419_01_009.mat','Data_20110419_01_010.mat'],
                    2012:['Data_20120418_01_129.mat','Data_20120418_01_130.mat','Data_20120418_01_131.mat'],
                    2013:['Data_20130405_01_165.mat','Data_20130405_01_166.mat','Data_20130405_01_167.mat'],
                    2014:['Data_20140424_01_002.mat','Data_20140424_01_003.mat','Data_20140424_01_004.mat'],
                    2017:['Data_20170422_01_168.mat','Data_20170422_01_169.mat','Data_20170422_01_170.mat','Data_20170422_01_171.mat'],
                    2018:['Data_20180427_01_170.mat','Data_20180427_01_171.mat','Data_20180427_01_172.mat']}
'''
#Very good candidate for ice slabs filling!!
investigation_year={2010:['Data_20100513_01_001.mat','Data_20100513_01_002.mat'],
                    2011:['Data_20110411_01_116.mat','Data_20110411_01_117.mat','Data_20110411_01_118.mat'],
                    2012:['Data_20120428_01_125.mat','Data_20120428_01_126.mat'],
                    2013:'empty',
                    2014:['Data_20140408_11_024.mat','Data_20140408_11_025.mat','Data_20140408_11_026.mat'],
                    2017:['Data_20170508_02_165.mat','Data_20170508_02_166.mat','Data_20170508_02_167.mat','Data_20170508_02_168.mat','Data_20170508_02_169.mat','Data_20170508_02_170.mat','Data_20170508_02_171.mat']}

'''
'''
#investigation_year={2010:['Data_20100507_01_008.mat','Data_20100507_01_009.mat','Data_20100507_01_010.mat'],
#                    2011:['Data_20110426_01_009.mat','Data_20110426_01_010.mat','Data_20110426_01_011.mat'],
#                    2012:'empty',
#                    2013:'empty',
#                    2014:['Data_20140421_01_009.mat','Data_20140421_01_010.mat','Data_20140421_01_011.mat','Data_20140421_01_012.mat','Data_20140421_01_013.mat'],
#                    2017:['Data_20170424_01_008.mat','Data_20170424_01_009.mat','Data_20170424_01_010.mat','Data_20170424_01_011.mat','Data_20170424_01_012.mat']}

#This one is collocated with FS1, 2, 3.
#investigation_year={2010:'empty',
#                    2011:'empty',
#                    2012:['Data_20120423_01_137.mat','Data_20120423_01_138.mat'],
#                    2013:['Data_20130409_01_010.mat','Data_20130409_01_011.mat','Data_20130409_01_012.mat'],
#                    2014:'empty',
#                    2017:'empty',
#                    2018:['Data_20180421_01_004.mat','Data_20180421_01_005.mat','Data_20180421_01_006.mat','Data_20180421_01_007.mat']}

#investigation_year={2010:'empty',
#                    2011:'empty',
#                    2012:['Data_20120423_01_006.mat','Data_20120423_01_007.mat'],
#                    2013:'empty',
#                    2014:'empty',
#                    2017:['Data_20170505_02_008.mat','Data_20170505_02_009.mat','Data_20170505_02_010.mat'],
#                    2018:'empty'}

#investigation_year={2010:'empty',
#                    2011:'empty',
#                    2012:['Data_20120423_01_137.mat','Data_20120423_01_138.mat'],
#                    2013:['Data_20130409_01_010.mat','Data_20130409_01_011.mat','Data_20130409_01_012.mat'],
#                    2014:'empty',
#                    2017:'empty',
#                    2018:['Data_20180421_01_004.mat','Data_20180421_01_005.mat','Data_20180421_01_006.mat','Data_20180421_01_007.mat']}

#investigation_year={2010:['Data_20100512_04_073.mat','Data_20100512_04_074.mat'],
#                    2011:'empty',
#                    2012:'empty',
#                    2013:'empty',
#                    2014:'empty',
#                    2017:'empty',
#                    2018:['Data_20180425_01_166.mat','Data_20180425_01_167.mat','Data_20180425_01_168.mat','Data_20180425_01_169.mat']}

#investigation_year={2010:'empty',
#                    2011:'empty',
#                    2012:'empty',
#                    2013:['Data_20130405_01_011.mat','Data_20130405_01_012.mat','Data_20130405_01_013.mat'],
#                    2014:['Data_20140424_03_046.mat','Data_20140424_03_047.mat','Data_20140424_03_048.mat'],
#                    2017:['Data_20170422_01_007.mat','Data_20170422_01_008.mat','Data_20170422_01_009.mat'],
#                    2018:'empty'}

#investigation_year={2010:'empty',
#                    2011:['Data_20110407_01_009.mat','Data_20110407_01_010.mat','Data_20110407_01_011.mat','Data_20110407_01_012.mat','Data_20110407_01_013.mat'],
#                    2012:'empty',
#                    2013:'empty',
#                    2014:'empty',
#                    2017:['Data_20170510_02_151.mat','Data_20170510_02_152.mat','Data_20170510_02_153.mat','Data_20170510_02_154.mat','Data_20170510_02_155.mat','Data_20170510_02_156.mat'],
#                    2018:'empty'}

#investigation_year={2010:['Data_20100514_02_035.mat','Data_20100514_02_036.mat','Data_20100514_02_037.mat','Data_20100514_02_038.mat','Data_20100514_02_039.mat',],
#                    2011:['Data_20110406_01_144.mat','Data_20110406_01_145.mat','Data_20110406_01_146.mat',],
#                    2012:['Data_20120413_01_006.mat','Data_20120413_01_007.mat','Data_20120413_01_008.mat','Data_20120413_01_009.mat','Data_20120413_01_010.mat','Data_20120413_01_011.mat','Data_20120413_01_012.mat'],
#                    2013:['Data_20130404_01_139.mat','Data_20130404_01_140.mat','Data_20130404_01_141.mat','Data_20130404_01_142.mat','Data_20130404_01_143.mat','Data_20130404_01_144.mat','Data_20130404_01_145.mat'],
#                    2014:['Data_20140409_10_057.mat','Data_20140409_10_058.mat','Data_20140409_10_059.mat','Data_20140409_10_060.mat','Data_20140409_10_061.mat','Data_20140409_10_062.mat','Data_20140409_10_063.mat','Data_20140409_10_064.mat','Data_20140409_10_065.mat','Data_20140409_10_066.mat'],
#                    2017:['Data_20170429_01_148.mat','Data_20170429_01_149.mat','Data_20170429_01_150.mat','Data_20170429_01_151.mat','Data_20170429_01_152.mat','Data_20170429_01_153.mat','Data_20170429_01_154.mat'],
#                    2018:'empty'}

investigation_year={2010:['Data_20100514_02_001.mat','Data_20100514_02_002.mat','Data_20100514_02_003.mat'],
                    2011:['Data_20110406_01_108.mat','Data_20110406_01_109.mat','Data_20110406_01_110.mat','Data_20110406_01_111.mat','Data_20110406_01_112.mat'],
                    2012:['Data_20120429_01_016.mat','Data_20120429_01_017.mat','Data_20120429_01_018.mat','Data_20120429_01_019.mat','Data_20120429_01_020.mat'],
                    2013:['Data_20130404_01_103.mat','Data_20130404_01_104.mat','Data_20130404_01_105.mat','Data_20130404_01_106.mat','Data_20130404_01_107.mat'],
                    2014:['Data_20140409_10_022.mat','Data_20140409_10_023.mat','Data_20140409_10_024.mat','Data_20140409_10_025.mat'],
                    2017:['Data_20170429_01_112.mat','Data_20170429_01_113.mat','Data_20170429_01_114.mat','Data_20170429_01_115.mat','Data_20170429_01_116.mat','Data_20170429_01_117.mat'],
                    2018:['Data_20180430_01_103.mat','Data_20180430_01_104.mat','Data_20180430_01_105.mat','Data_20180430_01_106.mat','Data_20180430_01_107.mat','Data_20180430_01_108.mat',]}

investigation_year={2010:'empty',
                    2011:['Data_20110422_02_001.mat','Data_20110422_02_002.mat','Data_20110422_02_003.mat','Data_20110422_02_004.mat','Data_20110422_02_005.mat','Data_20110422_02_006.mat'],
                    2012:'empty',
                    2013:['Data_20130406_01_078.mat','Data_20130406_01_079.mat','Data_20130406_01_080.mat','Data_20130406_01_081.mat','Data_20130406_01_082.mat','Data_20130406_01_083.mat','Data_20130406_01_084.mat','Data_20130406_01_085.mat','Data_20130406_01_086.mat'],
                    2014:'empty',
                    2017:['Data_20170501_02_077.mat','Data_20170501_02_078.mat','Data_20170501_02_079.mat','Data_20170501_02_080.mat','Data_20170501_02_081.mat','Data_20170501_02_082.mat','Data_20170501_02_083.mat','Data_20170501_02_084.mat','Data_20170501_02_085.mat','Data_20170501_02_086.mat','Data_20170501_02_087.mat'],
                    2018:'empty'}


investigation_year={2010:['Data_20100514_02_035.mat','Data_20100514_02_036.mat','Data_20100514_02_037.mat','Data_20100514_02_038.mat','Data_20100514_02_039.mat'],
                    2011:'empty',
                    2012:'empty',
                    2013:['Data_20130404_01_139.mat','Data_20130404_01_140.mat','Data_20130404_01_141.mat','Data_20130404_01_142.mat','Data_20130404_01_143.mat','Data_20130404_01_144.mat','Data_20130404_01_145.mat'],
                    2014:['Data_20140409_10_057.mat','Data_20140409_10_058.mat','Data_20140409_10_059.mat','Data_20140409_10_060.mat','Data_20140409_10_061.mat','Data_20140409_10_062.mat','Data_20140409_10_063.mat','Data_20140409_10_064.mat','Data_20140409_10_065.mat','Data_20140409_10_066.mat'],
                    2017:['Data_20170429_01_148.mat','Data_20170429_01_149.mat','Data_20170429_01_150.mat','Data_20170429_01_151.mat','Data_20170429_01_152.mat','Data_20170429_01_153.mat','Data_20170429_01_154.mat'],
                    2018:['Data_20180430_01_139.mat','Data_20180430_01_140.mat','Data_20180430_01_141.mat','Data_20180430_01_142.mat','Data_20180430_01_143.mat','Data_20180430_01_144.mat','Data_20180430_01_145.mat']}



investigation_year={2010:'empty',
                    2011:['Data_20110422_02_070.mat','Data_20110422_02_071.mat','Data_20110422_02_072.mat','Data_20110422_02_073.mat','Data_20110422_02_074.mat','Data_20110422_02_075.mat','Data_20110422_02_076.mat'],
                    2012:['Data_20120429_01_016.mat','Data_20120429_01_017.mat','Data_20120429_01_018.mat','Data_20120429_01_019.mat','Data_20120429_01_020.mat'],
                    2013:['Data_20130406_01_146.mat','Data_20130406_01_147.mat','Data_20130406_01_148.mat','Data_20130406_01_149.mat','Data_20130406_01_150.mat','Data_20130406_01_151.mat','Data_20130406_01_152.mat'],
                    2014:'empty',
                    2017:['Data_20170501_02_006.mat','Data_20170501_02_007.mat','Data_20170501_02_008.mat','Data_20170501_02_009.mat','Data_20170501_02_010.mat','Data_20170501_02_011.mat','Data_20170501_02_012.mat','Data_20170501_02_013.mat','Data_20170501_02_014.mat','Data_20170501_02_015.mat','Data_20170501_02_016.mat'],
                    2018:'empty'}




investigation_year={2010:['Data_20100514_02_009.mat','Data_20100514_02_010.mat'],
                    2011:'empty',
                    2012:'empty',
                    2013:['Data_20130404_01_113.mat','Data_20130404_01_114.mat','Data_20130404_01_115.mat','Data_20130404_01_116.mat','Data_20130404_01_117.mat','Data_20130404_01_118.mat','Data_20130404_01_119.mat'],
                    2014:['Data_20140409_10_033.mat'],
                    2017:['Data_20170429_01_122.mat','Data_20170429_01_123.mat','Data_20170429_01_124.mat','Data_20170429_01_125.mat','Data_20170429_01_126.mat','Data_20170429_01_127.mat','Data_20170429_01_128.mat','Data_20170429_01_129.mat'],
                    2018:'empty'}

'''





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
#Define the path for boolean images
path_boolean_images=path_data+'exported/temp_for_overlap/'

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
              
    #Define boolean filename and boolean image filename
    if (plot_boolean=='plot_boolean_orig_cut045_th000'):
        filename_boolean=date_track+'_orig_CUTOFF_-0.45_THRESHOLD_000.pickle'
        filename_boolean_image=date_track+'_orig_CUTOFF_-0.45_THRESHOLD_000.png'
        
    if (plot_boolean=='plot_boolean_SG1_cut045_th000'):
        filename_boolean=date_track+'_SG1_CUTOFF_-0.45_THRESHOLD_000.pickle'
        filename_boolean_image=date_track+'_SG1_CUTOFF_-0.45_THRESHOLD_000.png'
        
    if (plot_boolean=='plot_boolean_SG1_cut045_th350'):
        filename_boolean=date_track+'_SG1_CUTOFF_-0.45_THRESHOLD_350.pickle'
        filename_boolean_image=date_track+'_SG1_CUTOFF_-0.45_THRESHOLD_350.png'
    
    #Open the depth corrected file
    f_depth_corrected = open(path_depth_corrected+filename_depth_corrected, "rb")
    radar = pickle.load(f_depth_corrected)
    f_depth_corrected.close()
    
    #Open boolean file
    f_boolean = open(path_boolean+filename_boolean, "rb")
    boolean_file = pickle.load(f_boolean)
    f_boolean.close()
    
    #Open boolean image
    boolean_image = Image.open(path_boolean_images+filename_boolean_image).convert("L")
    arr_boolean_image = np.asarray(boolean_image)
    
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
        boolean_file=np.fliplr(boolean_file)
        arr_boolean_image=np.fliplr(arr_boolean_image)
        data_mask=np.flipud(data_mask)
    
    #Store reunited lat/lon, slice output and mask in a dictionnary:
    dataframe[str(single_year)]={'lat_appended':lat_appended,
                                 'lon_appended':lon_appended,
                                 'radar':radar,
                                 'boolean':boolean_file,
                                 'boolean_image':arr_boolean_image,
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
    pyplot.rcParams.update({'font.size': 5})
    fig1, (ax1s,ax2s,ax3s,ax4s,ax5s,ax6s,ax7s) = pyplot.subplots(7, 1)
    
if (plot_boolean_subplot=='TRUE'):
    pyplot.rcParams['axes.linewidth'] = 0.1 #set the value globally
    pyplot.rcParams['xtick.major.width']=0.1
    pyplot.rcParams['ytick.major.width']=0.1
    
    pyplot.figure(figsize=(40,20))
    pyplot.rcParams.update({'font.size': 5})
    fig2, (ax1b,ax2b,ax3b,ax4b,ax5b,ax6b,ax7b) = pyplot.subplots(7, 1)
    
if (plot_images_subplot=='TRUE'):
    pyplot.rcParams['axes.linewidth'] = 0.1 #set the value globally
    pyplot.rcParams['xtick.major.width']=0.1
    pyplot.rcParams['ytick.major.width']=0.1
    
    pyplot.figure(figsize=(40,20))
    pyplot.rcParams.update({'font.size': 5})
    fig3, (ax1i,ax2i,ax3i,ax4i,ax5i,ax6i,ax7i) = pyplot.subplots(7, 1)
        
for single_year in investigation_year.keys():
    
    print(str(single_year))
    
    #If no data, continue
    if (investigation_year[single_year]=='empty'):
        print('No data for year '+str(single_year)+', continue')
        continue
    
    start_date_track=investigation_year[single_year][0]
    end_date_track=investigation_year[single_year][-1]
    date_track=start_date_track[5:20]+'_'+end_date_track[17:20]
    
    #Calculate distances (in m)
    distances=compute_distances(dataframe[str(single_year)]['lon_appended'],dataframe[str(single_year)]['lat_appended'])
    
    #Convert distances from m to km
    distances=distances/1000
    
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
        #ax1.set_aspect(0.0025) # X scale matches Y scale
        ax1.set_title(date_track+' Depth corrected')
        ax1.set_ylabel('Depth [m]')
        ax1.set_xlabel('Longitude [°]')
        #pdb.set_trace()
        ax1.set_xlim(min_lon,max_lon)
        ax1.set_ylim(dataframe[str(single_year)]['boolean'].shape[0],0)
        
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
        
        #If no data, continue
        if (investigation_year[single_year]=='empty'):
            print('No data for year '+str(single_year)+', continue')
            continue
    
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
        ax_plotting.set_aspect(0.001) # X scale matches Y scale
        ax_plotting.set_title(str(single_year)+' Depth corrected')
        ax_plotting.set_ylabel('Depth [m]')
        ax_plotting.set_xlabel('Longitude [°]]')
        ax_plotting.set_xlim(-48,-46.6)
        ax_plotting.set_ylim(30,0)
        
        #ax_plotting.set_xticks(distances)
        #ax_plotting.set_xticklabels(distances)
        
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
        
        #If no data, continue
        if (investigation_year[single_year]=='empty'):
            print('No data for year '+str(single_year)+', continue')
            continue
        
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
        Y=np.arange(0,20,20/dataframe[str(single_year)]['boolean'].shape[0])
        C=dataframe[str(single_year)]['boolean'].astype(float)
        
        cb=ax_plotting.pcolor(X, Y, C,cmap=pyplot.get_cmap('gray_r'))#,norm=divnorm)
        ax_plotting.invert_yaxis() #Invert the y axis = avoid using flipud.
        ax_plotting.set_aspect(0.001) # X scale matches Y scale
        ax_plotting.set_title(str(single_year)+' '+plot_boolean)
        ax_plotting.set_ylabel('Depth [m]')
        ax_plotting.set_xlabel('Longitude [°]')
        ax_plotting.set_xlim(min_lon,max_lon)
        ax_plotting.set_ylim(20,0)
        
        cbar=fig1.colorbar(cb, ax=[ax_plotting], location='right',shrink=0.12,aspect=10,pad=0.01)
        cbar.set_label('Signal strength')
        
        ##Create the figure name
        #fig_name=[]
        #fig_name='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/2010_2014_thickening/'+date_track+'_'+plot_boolean+'.png'
        
        pyplot.show()
        
        ##Save the figure
        #pyplot.savefig(fig_name,dpi=2000)
        #pyplot.clf()        

    if (plot_images_subplot=='TRUE'):
        
        #If no data, continue
        if (investigation_year[single_year]=='empty'):
            print('No data for year '+str(single_year)+', continue')
            continue
        
        if (single_year==2010):
            ax_plotting=ax1i
        elif (single_year==2011):
            ax_plotting=ax2i
        elif (single_year==2012):
            ax_plotting=ax3i
        elif (single_year==2013):
            ax_plotting=ax4i
        elif (single_year==2014):
            ax_plotting=ax5i
        elif (single_year==2017):
            ax_plotting=ax6i
        elif (single_year==2018):
            ax_plotting=ax7i
        else:
            print('Year not existing')
        
        
        #Prepare the plot
        #ax_plotting.set_xticks(ticks=np.linspace(min_lon,max_lon,dataframe[str(single_year)]['boolean'].shape[1]))

        X=dataframe[str(single_year)]['lon_appended']
        Y=np.arange(0,20,20/dataframe[str(single_year)]['boolean'].shape[0])
        C=dataframe[str(single_year)]['boolean_image'].astype(float)
        
        cb=ax_plotting.pcolor(X, Y, C,cmap=pyplot.get_cmap('gray'))#,norm=divnorm)
        ax_plotting.invert_yaxis() #Invert the y axis = avoid using flipud.
        ax_plotting.set_aspect(0.001) # X scale matches Y scale
        ax_plotting.set_title(str(single_year)+' '+plot_boolean)
        ax_plotting.set_ylabel('Depth [m]')
        ax_plotting.set_xlabel('Longitude [°]')
        ax_plotting.set_xlim(min_lon,max_lon)
        ax_plotting.set_ylim(20,0)
        
        cbar=fig1.colorbar(cb, ax=[ax_plotting], location='right',shrink=0.12,aspect=10,pad=0.01)
        cbar.set_label('Signal strength')
        
        ##Create the figure name
        #fig_name=[]
        #fig_name='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/2010_2014_thickening/'+date_track+'_'+plot_boolean+'.png'
        
        pyplot.show()
        ##Save the figure
        #pyplot.savefig(fig_name,dpi=2000)
        #pyplot.clf()    

'''
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
'''
##############################################################################
###                                 Plot data                              ###
##############################################################################

import matplotlib.pyplot as plt
import rasterio as rio
import geopandas as gpd
##########################################################################
###                          Load excess melt data   	               ###
##########################################################################
#This is from the code excess_melt.py
#Define path
path='C:/Users/jullienn/Documents/working_environment/excess_melt/'

#Load the data
data_path = path+'excess_melt_mbyear.nc'
DS = xr.open_dataset(data_path)

#Extract coordinates
lat_M_e=DS.x.data
lon_M_e=DS.y.data

#Load melt data
melt_data= DS.M_e

##########################################################################
###                          Load excess melt data   	               ###
##########################################################################

##########################################################################
###                      Load Greenland DEM and contours               ###
##########################################################################
#This is from the code excess_melt.py
dem_filename = 'C:/Users/jullienn/switchdrive/Private/research/backup_Aglaja/working_environment/greenland_topo_data/elevations/greenland_dem_mosaic_100m_v3.0.tif'
contours_filename = 'C:/Users/jullienn/switchdrive/Private/research/backup_Aglaja/working_environment/greenland_topo_data/elevations/greenland_contours_100m_v3.0.shp'

with rio.open(dem_filename) as fh:
    dem = fh.read()
    print(fh)
    dem_bounds = fh.bounds
    dem_crs = fh.crs
contours = gpd.read_file(contours_filename)
# Make sure that the coordinate reference system matches the image
contours = contours.to_crs(dem_crs)

##########################################################################
###                      Load Greenland DEM and contours               ###
##########################################################################


##########################################################################
###                          Plot excess melt data   	               ###
##########################################################################
#This is from the code excess_melt.py

#Plot the excess melt with ice bridge track on top.
#Note: I will have to work with excess melt differences: typically if I show 2012: I will have to do 2011-2010.
generate_raw_excess_melt='TRUE'

if (generate_raw_excess_melt=='TRUE'):
    #Generate and save the raw annual excess melt figures
    for year in list(np.arange(1990,2020)):
        
        #Define the year
        wanted_year=str(year)
        
        #Select the data associated with the wanted year
        melt_year = melt_data.sel(time=wanted_year)
        melt_year_np = melt_year.values
        
        melt_year_plot=np.asarray(melt_year_np)
        melt_year_plot=melt_year_plot[0,0,:,:,0]
        
        #Plot dem and contours elevation
        plt.rcParams.update({'font.size': 20})
        plt.figure(figsize=(48,40))
        ax = plt.subplot(111)
        dem_extent = (dem_bounds[0], dem_bounds[2], dem_bounds[1], dem_bounds[3])
        plt.imshow(np.squeeze(np.flipud(melt_year_plot)), extent=dem_extent,cmap=discrete_cmap(5,'hot_r'))
        
        plt.colorbar(label='Excess melt [mm w.e./year]')
        plt.clim(0,1500)
        ax.grid()
        #contours.plot(ax=ax, edgecolor='black')
        plt.title('Excess melt plot, year: '+wanted_year)
        plt.show()

##########################################################################
###                          Plot excess melt data   	              ###
##########################################################################


##########################################################################
###                Plot difference excess melt data   	              ###
##########################################################################
cum_excess=np.zeros((446,240))
    
for year in list(dataframe.keys()):
    #Transform the coordinates from WGS84 to EPSG:3413
    #Example from: https://pyproj4.github.io/pyproj/stable/examples.html
    transformer = Transformer.from_crs("EPSG:4326", "EPSG:3413", always_xy=True)
    points=transformer.transform(np.array(dataframe[str(year)]['lon_appended']),np.array(dataframe[str(year)]['lat_appended']))
    
    #Reset the lat_3413 and lon_3413 to empty vectors.
    lon_3413=[]
    lat_3413=[]
    
    lon_3413=points[0]
    lat_3413=points[1]
    
    #Select the data associated with the reference year
    melt_year_ref = melt_data.sel(time='2009')
    melt_year_np_ref = melt_year_ref.values
    
    melt_year_plot_ref=np.asarray(melt_year_np_ref)
    melt_year_plot_ref=melt_year_plot_ref[0,0,:,:,0]
    
    if (yearly_comparison_ref=='TRUE'):
            
        if (str(year) in list(['2010','2011','2012','2013','2014','2018'])):
            #Select the data associated with the wanted year -1
            melt_year = melt_data.sel(time=str(int(year)-1))
            melt_year_np = melt_year.values
            
            melt_year_plot=np.asarray(melt_year_np)
            melt_year_plot=melt_year_plot[0,0,:,:,0]
            
            #Calculate the difference between the melt year of interest compared to the reference year
            diff_melt_year_plot = melt_year_plot-melt_year_plot_ref
            
            #Create title for plot
            title_to_plot='Year: '+year+', Excess melt: '+str(int(year)-1)+'-2009'
            
        elif (year == str(2017)):
            
            for indiv_year in range(2014,2017):
                print(indiv_year)
                #Select the data associated with the wanted year
                melt_year = melt_data.sel(time=str(indiv_year))
                melt_year_np = melt_year.values
                
                melt_year_plot=np.asarray(melt_year_np)
                melt_year_plot=melt_year_plot[0,0,:,:,0]
                
                #Calculate the difference between the melt year of interest compared to the reference year
                temp_diff=melt_year_plot-melt_year_plot_ref
                
                #As we are only interested in positive difference of excess melt,
                #set to zero where temp_diff<0
                
                x_loc,y_loc=np.where(temp_diff<0)
                temp_diff[x_loc,y_loc]=0
                
                #Calculate the cumulative difference between the year of interest and the reference years
                cum_excess=cum_excess+temp_diff
        
            #Store in the right variable for plotting
            diff_melt_year_plot=cum_excess/3 #Divide it by 3 because sum of 3 years
            
            #Create title for plot
            title_to_plot='Year: '+year+', Excess melt: (abs(2014-2009)+abs(2015-2009)+abs(2016-2009))/3'
        else:
            print('Year input is not known')
            break
    
    if (yearly_comparison_indiv=='TRUE'):
        if (year==str(2010)):
            #Select the data associated with the wanted year -1
            melt_year = melt_data.sel(time=str(int(year)-1))
            melt_year_np = melt_year.values
            
            melt_year_plot=np.asarray(melt_year_np)
            melt_year_plot=melt_year_plot[0,0,:,:,0]
            
            #Calculate the difference between the two years
            diff_melt_year_plot = melt_year_plot-melt_year_plot_ref
            
            #Create title for plot
            title_to_plot='Year: '+year+', Excess melt: 2009-2009'
            
        elif (str(year) in list(['2011','2012','2013','2014','2018'])):
            #Select the data associated with the wanted year -1
            melt_year = melt_data.sel(time=str(int(year)-1))
            melt_year_np = melt_year.values
            
            melt_year_plot=np.asarray(melt_year_np)
            melt_year_plot=melt_year_plot[0,0,:,:,0]
            
            #Select the data associated with the wanted year -2
            melt_year_minus = melt_data.sel(time=str(int(year)-2))
            melt_year_minus_np = melt_year_minus.values
            
            melt_year_plot_minus=np.asarray(melt_year_minus_np)
            melt_year_plot_minus=melt_year_plot_minus[0,0,:,:,0]
            
            #Calculate the difference between the two years
            diff_melt_year_plot = melt_year_plot-melt_year_plot_minus
            
            #Create title for plot
            title_to_plot='Year: '+year+', Excess melt: '+str(int(year)-1)+'-'+str(int(year)-2)
            
        elif (year == str(2017)):
            
            for indiv_year in range(2014,2017):
                print(indiv_year)
                #Select the data associated with the wanted year
                melt_year = melt_data.sel(time=str(indiv_year))
                melt_year_np = melt_year.values
                
                melt_year_plot=np.asarray(melt_year_np)
                melt_year_plot=melt_year_plot[0,0,:,:,0]
                
                #Select the data associated with the previous year
                melt_year_minus = melt_data.sel(time=str((indiv_year-1)))
                melt_year_minus_np = melt_year_minus.values
                
                melt_year_plot_minus=np.asarray(melt_year_minus_np)
                melt_year_plot_minus=melt_year_plot_minus[0,0,:,:,0]
                
                #Calculate the difference between the melt year of interest compared to the reference year
                temp_diff=melt_year_plot-melt_year_plot_minus
                
                #As we are only interested in positive difference of excess melt,
                #set to zero where temp_diff<0
                
                x_loc,y_loc=np.where(temp_diff<0)
                temp_diff[x_loc,y_loc]=0
                
                #Calculate the cumulative difference between the year of interest and the reference years
                cum_excess=cum_excess+temp_diff
        
            #Store in the right variable for plotting
            pdb.set_trace()
            diff_melt_year_plot=cum_excess
            
            #Create title for plot
            title_to_plot='Year: '+year+', Excess melt: abs(2014-2013)+abs(2015-2014)+abs(2016-2015)'
        else:
            print('Year input is not known')
            break
            
    if (cumulative_comparison=='TRUE'):
        for indiv_year in range(2009,int(year)):
            #Select the data associated with the wanted year
            melt_year = melt_data.sel(time=str(indiv_year))
            melt_year_np = melt_year.values
            
            melt_year_plot=np.asarray(melt_year_np)
            melt_year_plot=melt_year_plot[0,0,:,:,0]
            
            #Calculate the cumulative difference between the year of interest and the reference years
            cum_excess=cum_excess+melt_year_plot-melt_year_plot_ref
        
        #Store in the right variable for plotting
        diff_melt_year_plot=cum_excess
    
    #Plot dem bounds and excess melt
    plt.rcParams.update({'font.size': 10})
    plt.figure(figsize=(48,40))
    ax = plt.subplot(111)
    dem_extent = (dem_bounds[0], dem_bounds[2], dem_bounds[1], dem_bounds[3])
    plt.imshow(np.squeeze(np.flipud(diff_melt_year_plot)), extent=dem_extent,cmap=discrete_cmap(17,'RdBu_r'))
    
    #Plot ice bridge trace over excess melt
    plt.plot(lon_3413,lat_3413,marker='o', markersize=1, zorder=45, color='blue')
    #Zoom over the trace
    plt.xlim(lon_3413[int(np.round(lat_3413.size/2))]-500000,lon_3413[int(np.round(lat_3413.size/2))]+500000)
    plt.ylim(lat_3413[int(np.round(lat_3413.size/2))]-300000,lat_3413[int(np.round(lat_3413.size/2))]+300000)
                    
    plt.colorbar(label='Cumulative excess melt difference[mm w.e./year]')
    plt.clim(-1500,1500)
    ax.grid()
    #contours.plot(ax=ax, edgecolor='black')
    #plt.title('Trace year: '+str(year)+'. Cumulative: 2009-'+str(int(year)-1)+' -x.2009')
    plt.title(title_to_plot)
    plt.show()

pdb.set_trace()

##########################################################################
###                Plot difference excess melt data   	              ###
##########################################################################

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
    
    

    
    
        
#lat1.shape[1]+lat2.shape[1]+lat3.shape[1]= data.shape[1] => great news!

#To do:
#    1. append data and lat/lon1 OK
#    2. load data from 2010 to 2014 OK
#    3. select as a function of lat/lon box OK
#    4. plot the data OK
#    5. compare! OK
#    6. Decide whether I should work with depth corrected, or other post-precessed files that are available in 'Boolean Array Picklefiles' folder OK