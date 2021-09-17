# -*- coding: utf-8 -*-
"""
Created on Tue Sep 14 11:45:45 2021

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

#Import librairies
import pandas as pd
import datetime
from scipy import spatial
import pdb
import numpy as np
from pyproj import Transformer
import matplotlib.pyplot as plt
from pysheds.grid import Grid
import matplotlib.colors as mcolors

#Load the spatial aggregated data. All the points within a radius of 100m are averaged
path='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2010_2018/excel_files/'
df_2010_2018_spatially_aggregated = pd.read_csv(path+'jullien_etal_20102018_spatial_aggregation_grid_1000.csv',delimiter=';',decimal=',')

#Load all 2010-2018 data without spatial aggregation
df_2010_2018 = pd.read_csv(path+'jullienetal_20102018.csv',delimiter=';',decimal=',')
#Transform the coordinated from WGS84 to EPSG:3413
transformer = Transformer.from_crs("EPSG:4326", "EPSG:3413", always_xy=True)
points=transformer.transform(np.asarray(df_2010_2018["lon"]),np.asarray(df_2010_2018["lat"]))

#Store lat/lon in 3413
df_2010_2018['lon_3413']=points[0]
df_2010_2018['lat_3413']=points[1]

#Visualize the spatial aggregation process
########################## Load GrIS elevation ##########################
#Open the DEM
grid = Grid.from_raster("C:/Users/jullienn/switchdrive/Private/research/backup_Aglaja/working_environment/greenland_topo_data/elevations/greenland_dem_mosaic_100m_v3.0.tif",data_name='dem')
#Minnor slicing on borders to enhance colorbars
elevDem=grid.dem[:-1,:-1]              
#Scale the colormap
divnorm = mcolors.DivergingNorm(vmin=0, vcenter=1250, vmax=2500)
########################## Load GrIS elevation ##########################

fig, (ax1) = plt.subplots(1, 1)#, gridspec_kw={'width_ratios': [1, 3]})
fig.suptitle('Spatial aggregation illustration')
#Display DEM
cb1=ax1.imshow(elevDem, extent=grid.extent,cmap=discrete_cmap(10,'cubehelix_r'),alpha=0.5,norm=divnorm)
cbar1=fig.colorbar(cb1, ax=[ax1], location='left')
cbar1.set_label('Elevation [m]')
#Display 2010-2018 data
plt.scatter(df_2010_2018[df_2010_2018.Track_name.str[:4]=='2010']['lon_3413'],df_2010_2018[df_2010_2018.Track_name.str[:4]=='2010']['lat_3413'],s=0.1,color='#525252',label='2010-2014')
plt.scatter(df_2010_2018[df_2010_2018.Track_name.str[:4]=='2011']['lon_3413'],df_2010_2018[df_2010_2018.Track_name.str[:4]=='2011']['lat_3413'],s=0.1,color='#525252')
plt.scatter(df_2010_2018[df_2010_2018.Track_name.str[:4]=='2012']['lon_3413'],df_2010_2018[df_2010_2018.Track_name.str[:4]=='2012']['lat_3413'],s=0.1,color='#525252')
plt.scatter(df_2010_2018[df_2010_2018.Track_name.str[:4]=='2013']['lon_3413'],df_2010_2018[df_2010_2018.Track_name.str[:4]=='2013']['lat_3413'],s=0.1,color='#525252')
plt.scatter(df_2010_2018[df_2010_2018.Track_name.str[:4]=='2014']['lon_3413'],df_2010_2018[df_2010_2018.Track_name.str[:4]=='2014']['lat_3413'],s=0.1,color='#525252')
plt.scatter(df_2010_2018[df_2010_2018.Track_name.str[:4]=='2017']['lon_3413'],df_2010_2018[df_2010_2018.Track_name.str[:4]=='2017']['lat_3413'],s=0.1,color='#d9d9d9',label='2017-2018')
plt.scatter(df_2010_2018[df_2010_2018.Track_name.str[:4]=='2018']['lon_3413'],df_2010_2018[df_2010_2018.Track_name.str[:4]=='2018']['lat_3413'],s=0.1,color='#d9d9d9')
#Display the aggregation data
plt.scatter(df_2010_2018_spatially_aggregated['avg_lon_3413'],df_2010_2018_spatially_aggregated['avg_lat_3413'],color='green',label='2010-2018 spatially aggregated')
plt.legend()
plt.show()


#Loop over the keys and create 1 dataframe per year. Where no data for this particular year, store a nan

#Create empty arrays
nan_array=np.zeros((1,7))
nan_array[:]=np.nan

array_2010=np.zeros((1,7))
array_2010[:]=np.nan

array_2011=np.zeros((1,7))
array_2011[:]=np.nan

array_2012=np.zeros((1,7))
array_2012[:]=np.nan

array_2013=np.zeros((1,7))
array_2013[:]=np.nan

array_2014=np.zeros((1,7))
array_2014[:]=np.nan

array_2017=np.zeros((1,7))
array_2017[:]=np.nan

array_2018=np.zeros((1,7))
array_2018[:]=np.nan

#Loop over the keys
for indiv_key in np.unique(df_2010_2018_spatially_aggregated.key):
    
    #Select data where key
    df_key=df_2010_2018_spatially_aggregated[df_2010_2018_spatially_aggregated.key==indiv_key]
    
    # ----- 2010
    if (not(2010 in np.asarray(df_key.year))):
        #No data for this year, store NaNs
        array_2010=np.append(array_2010,nan_array,axis=0)
    else:
        #There are data for this year, store them
        array_2010=np.append(array_2010,np.asarray(df_key[df_key.year==2010]),axis=0)
    # ----- 2010
    
    # ----- 2011
    if (not(2011 in np.asarray(df_key.year))):
        #No data for this year, store NaNs
        array_2011=np.append(array_2011,nan_array,axis=0)
    else:
        #There are data for this year, store them
        array_2011=np.append(array_2011,np.asarray(df_key[df_key.year==2011]),axis=0)
    # ----- 2011
    
    # ----- 2012
    if (not(2012 in np.asarray(df_key.year))):
        #No data for this year, store NaNs
        array_2012=np.append(array_2012,nan_array,axis=0)
    else:
        #There are data for this year, store them
        array_2012=np.append(array_2012,np.asarray(df_key[df_key.year==2012]),axis=0)
    # ----- 2012
    
    # ----- 2013
    if (not(2013 in np.asarray(df_key.year))):
        #No data for this year, store NaNs
        array_2013=np.append(array_2013,nan_array,axis=0)
    else:
        #There are data for this year, store them
        array_2013=np.append(array_2013,np.asarray(df_key[df_key.year==2013]),axis=0)
    # ----- 2013
    
    # ----- 2014
    if (not(2014 in np.asarray(df_key.year))):
        #No data for this year, store NaNs
        array_2014=np.append(array_2014,nan_array,axis=0)
    else:
        #There are data for this year, store them
        array_2014=np.append(array_2014,np.asarray(df_key[df_key.year==2014]),axis=0)
    # ----- 2014
    
    # ----- 2017
    if (not(2017 in np.asarray(df_key.year))):
        #No data for this year, store NaNs
        array_2017=np.append(array_2017,nan_array,axis=0)
    else:
        #There are data for this year, store them
        array_2017=np.append(array_2017,np.asarray(df_key[df_key.year==2017]),axis=0)
    # ----- 2017
    
    # ----- 2018
    if (not(2018 in np.asarray(df_key.year))):
        #No data for this year, store NaNs
        array_2018=np.append(array_2018,nan_array,axis=0)
    else:
        #There are data for this year, store them
        array_2018=np.append(array_2018,np.asarray(df_key[df_key.year==2018]),axis=0)
    # ----- 2018
    
    print(indiv_key/len(np.unique(df_2010_2018_spatially_aggregated.key))*100,' %')

#Delete the first line of all the array_year because NaNs
array_2010=np.delete(array_2010,0,0)
array_2011=np.delete(array_2011,0,0)
array_2012=np.delete(array_2012,0,0)
array_2013=np.delete(array_2013,0,0)
array_2014=np.delete(array_2014,0,0)
array_2017=np.delete(array_2017,0,0)
array_2018=np.delete(array_2018,0,0)

#Store array as dataframes
df_spatially_aggregated_2010=pd.DataFrame(data=array_2010,
                                          columns=['index', 'avg_20m_icecontent', 'std_20m_icecontent', 'avg_lat_3413','avg_lon_3413', 'year', 'key'])
df_spatially_aggregated_2011=pd.DataFrame(data=array_2011,
                                          columns=['index', 'avg_20m_icecontent', 'std_20m_icecontent', 'avg_lat_3413','avg_lon_3413', 'year', 'key'])
df_spatially_aggregated_2012=pd.DataFrame(data=array_2012,
                                          columns=['index', 'avg_20m_icecontent', 'std_20m_icecontent', 'avg_lat_3413','avg_lon_3413', 'year', 'key'])
df_spatially_aggregated_2013=pd.DataFrame(data=array_2013,
                                          columns=['index', 'avg_20m_icecontent', 'std_20m_icecontent', 'avg_lat_3413','avg_lon_3413', 'year', 'key'])
df_spatially_aggregated_2014=pd.DataFrame(data=array_2014,
                                          columns=['index', 'avg_20m_icecontent', 'std_20m_icecontent', 'avg_lat_3413','avg_lon_3413', 'year', 'key'])
df_spatially_aggregated_2017=pd.DataFrame(data=array_2017,
                                          columns=['index', 'avg_20m_icecontent', 'std_20m_icecontent', 'avg_lat_3413','avg_lon_3413', 'year', 'key'])
df_spatially_aggregated_2018=pd.DataFrame(data=array_2018,
                                          columns=['index', 'avg_20m_icecontent', 'std_20m_icecontent', 'avg_lat_3413','avg_lon_3413', 'year', 'key'])

#Plot differences
diff_to_plot=df_spatially_aggregated_2011-df_spatially_aggregated_2010

fig, (ax1) = plt.subplots(1, 1)#, gridspec_kw={'width_ratios': [1, 3]})
fig.suptitle('Spatial aggregation difference 2011-2010')
#Display DEM
cb1=ax1.imshow(elevDem, extent=grid.extent,cmap=discrete_cmap(10,'cubehelix_r'),alpha=0.5,norm=divnorm)
cbar1=fig.colorbar(cb1, ax=[ax1], location='left')
cbar1.set_label('Elevation [m]')
# Make the norm for difference plotting
divnorm_diff = mcolors.DivergingNorm(vmin=-5, vcenter=0, vmax=5)
#Display 2010 data
plt.scatter(df_2010_2018[df_2010_2018.Track_name.str[:4]=='2010']['lon_3413'],df_2010_2018[df_2010_2018.Track_name.str[:4]=='2010']['lat_3413'],s=0.1,color='#525252',label='2010')
#Display 2011 data
plt.scatter(df_2010_2018[df_2010_2018.Track_name.str[:4]=='2011']['lon_3413'],df_2010_2018[df_2010_2018.Track_name.str[:4]=='2011']['lat_3413'],s=0.1,color='#d9d9d9',label='2011')
#Display the difference between 2011 and 2010 if aggregated data
sc= ax1.scatter(df_spatially_aggregated_2011['avg_lon_3413'],df_spatially_aggregated_2011['avg_lat_3413'],c=diff_to_plot['avg_20m_icecontent'],cmap='seismic_r',norm=divnorm_diff)
cbar=fig.colorbar(sc)
cbar.set_label('Difference in iceslabs thickness', fontsize=15)
plt.legend()

#Allows to open plot in full size directly
figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()

plt.show()


diff_to_plot=df_spatially_aggregated_2012-df_spatially_aggregated_2011

fig, (ax1) = plt.subplots(1, 1)#, gridspec_kw={'width_ratios': [1, 3]})
fig.suptitle('Spatial aggregation difference 2012-2011')
#Display DEM
cb1=ax1.imshow(elevDem, extent=grid.extent,cmap=discrete_cmap(10,'cubehelix_r'),alpha=0.5,norm=divnorm)
cbar1=fig.colorbar(cb1, ax=[ax1], location='left')
cbar1.set_label('Elevation [m]')
# Make the norm for difference plotting
divnorm_diff = mcolors.DivergingNorm(vmin=-5, vcenter=0, vmax=5)
#Display 2011 data
plt.scatter(df_2010_2018[df_2010_2018.Track_name.str[:4]=='2011']['lon_3413'],df_2010_2018[df_2010_2018.Track_name.str[:4]=='2011']['lat_3413'],s=0.1,color='#525252',label='2011')
#Display 2012 data
plt.scatter(df_2010_2018[df_2010_2018.Track_name.str[:4]=='2012']['lon_3413'],df_2010_2018[df_2010_2018.Track_name.str[:4]=='2012']['lat_3413'],s=0.1,color='#d9d9d9',label='2012')
#Display the difference between 2011 and 2010 if aggregated data
sc= ax1.scatter(df_spatially_aggregated_2012['avg_lon_3413'],df_spatially_aggregated_2012['avg_lat_3413'],c=diff_to_plot['avg_20m_icecontent'],cmap='seismic_r',norm=divnorm_diff)
cbar=fig.colorbar(sc)
cbar.set_label('Difference in iceslabs thickness', fontsize=15)
plt.legend()

#Allows to open plot in full size directly
figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()

plt.show()



diff_to_plot=df_spatially_aggregated_2013-df_spatially_aggregated_2012

fig, (ax1) = plt.subplots(1, 1)#, gridspec_kw={'width_ratios': [1, 3]})
fig.suptitle('Spatial aggregation difference 2013-2012')
#Display DEM
cb1=ax1.imshow(elevDem, extent=grid.extent,cmap=discrete_cmap(10,'cubehelix_r'),alpha=0.5,norm=divnorm)
cbar1=fig.colorbar(cb1, ax=[ax1], location='left')
cbar1.set_label('Elevation [m]')
# Make the norm for difference plotting
divnorm_diff = mcolors.DivergingNorm(vmin=-5, vcenter=0, vmax=5)
#Display 2012 data
plt.scatter(df_2010_2018[df_2010_2018.Track_name.str[:4]=='2012']['lon_3413'],df_2010_2018[df_2010_2018.Track_name.str[:4]=='2012']['lat_3413'],s=0.1,color='#525252',label='2012')
#Display 2013 data
plt.scatter(df_2010_2018[df_2010_2018.Track_name.str[:4]=='2013']['lon_3413'],df_2010_2018[df_2010_2018.Track_name.str[:4]=='2013']['lat_3413'],s=0.1,color='#d9d9d9',label='2013')
#Display the difference between 2011 and 2010 if aggregated data
sc= ax1.scatter(df_spatially_aggregated_2013['avg_lon_3413'],df_spatially_aggregated_2013['avg_lat_3413'],c=diff_to_plot['avg_20m_icecontent'],cmap='seismic_r',norm=divnorm_diff)
cbar=fig.colorbar(sc)
cbar.set_label('Difference in iceslabs thickness', fontsize=15)
plt.legend()

#Allows to open plot in full size directly
figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()

plt.show()



diff_to_plot=df_spatially_aggregated_2014-df_spatially_aggregated_2013

fig, (ax1) = plt.subplots(1, 1)#, gridspec_kw={'width_ratios': [1, 3]})
fig.suptitle('Spatial aggregation difference 2014-2013')
#Display DEM
cb1=ax1.imshow(elevDem, extent=grid.extent,cmap=discrete_cmap(10,'cubehelix_r'),alpha=0.5,norm=divnorm)
cbar1=fig.colorbar(cb1, ax=[ax1], location='left')
cbar1.set_label('Elevation [m]')
# Make the norm for difference plotting
divnorm_diff = mcolors.DivergingNorm(vmin=-5, vcenter=0, vmax=5)
#Display 2013 data
plt.scatter(df_2010_2018[df_2010_2018.Track_name.str[:4]=='2013']['lon_3413'],df_2010_2018[df_2010_2018.Track_name.str[:4]=='2013']['lat_3413'],s=0.1,color='#525252',label='2013')
#Display 2014 data
plt.scatter(df_2010_2018[df_2010_2018.Track_name.str[:4]=='2014']['lon_3413'],df_2010_2018[df_2010_2018.Track_name.str[:4]=='2014']['lat_3413'],s=0.1,color='#d9d9d9',label='2014')
#Display the difference between 2011 and 2010 if aggregated data
sc= ax1.scatter(df_spatially_aggregated_2014['avg_lon_3413'],df_spatially_aggregated_2014['avg_lat_3413'],c=diff_to_plot['avg_20m_icecontent'],cmap='seismic_r',norm=divnorm_diff)
cbar=fig.colorbar(sc)
cbar.set_label('Difference in iceslabs thickness', fontsize=15)
plt.legend()

#Allows to open plot in full size directly
figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()

plt.show()



diff_to_plot=df_spatially_aggregated_2017-df_spatially_aggregated_2014

fig, (ax1) = plt.subplots(1, 1)#, gridspec_kw={'width_ratios': [1, 3]})
fig.suptitle('Spatial aggregation difference 2017-2014')
#Display DEM
cb1=ax1.imshow(elevDem, extent=grid.extent,cmap=discrete_cmap(10,'cubehelix_r'),alpha=0.5,norm=divnorm)
cbar1=fig.colorbar(cb1, ax=[ax1], location='left')
cbar1.set_label('Elevation [m]')
# Make the norm for difference plotting
divnorm_diff = mcolors.DivergingNorm(vmin=-5, vcenter=0, vmax=5)
#Display 2014 data
plt.scatter(df_2010_2018[df_2010_2018.Track_name.str[:4]=='2014']['lon_3413'],df_2010_2018[df_2010_2018.Track_name.str[:4]=='2014']['lat_3413'],s=0.1,color='#525252',label='2014')
#Display 2017 data
plt.scatter(df_2010_2018[df_2010_2018.Track_name.str[:4]=='2017']['lon_3413'],df_2010_2018[df_2010_2018.Track_name.str[:4]=='2017']['lat_3413'],s=0.1,color='#d9d9d9',label='2017')
#Display the difference between 2011 and 2010 if aggregated data
sc= ax1.scatter(df_spatially_aggregated_2017['avg_lon_3413'],df_spatially_aggregated_2017['avg_lat_3413'],c=diff_to_plot['avg_20m_icecontent'],cmap='seismic_r',norm=divnorm_diff)
cbar=fig.colorbar(sc)
cbar.set_label('Difference in iceslabs thickness', fontsize=15)
plt.legend()

#Allows to open plot in full size directly
figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()

plt.show()



diff_to_plot=df_spatially_aggregated_2018-df_spatially_aggregated_2017

fig, (ax1) = plt.subplots(1, 1)#, gridspec_kw={'width_ratios': [1, 3]})
fig.suptitle('Spatial aggregation difference 2018-2017')
#Display DEM
cb1=ax1.imshow(elevDem, extent=grid.extent,cmap=discrete_cmap(10,'cubehelix_r'),alpha=0.5,norm=divnorm)
cbar1=fig.colorbar(cb1, ax=[ax1], location='left')
cbar1.set_label('Elevation [m]')
# Make the norm for difference plotting
divnorm_diff = mcolors.DivergingNorm(vmin=-5, vcenter=0, vmax=5)
#Display 2017 data
plt.scatter(df_2010_2018[df_2010_2018.Track_name.str[:4]=='2017']['lon_3413'],df_2010_2018[df_2010_2018.Track_name.str[:4]=='2017']['lat_3413'],s=0.1,color='#525252',label='2017')
#Display 2018 data
plt.scatter(df_2010_2018[df_2010_2018.Track_name.str[:4]=='2018']['lon_3413'],df_2010_2018[df_2010_2018.Track_name.str[:4]=='2018']['lat_3413'],s=0.1,color='#d9d9d9',label='2018')
#Display the difference between 2011 and 2010 if aggregated data
sc= ax1.scatter(df_spatially_aggregated_2018['avg_lon_3413'],df_spatially_aggregated_2018['avg_lat_3413'],c=diff_to_plot['avg_20m_icecontent'],cmap='seismic_r',norm=divnorm_diff)
cbar=fig.colorbar(sc)
cbar.set_label('Difference in iceslabs thickness', fontsize=15)
plt.legend()

#Allows to open plot in full size directly
figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()

plt.show()





#2017-2010
diff_to_plot=df_spatially_aggregated_2017-df_spatially_aggregated_2010

fig, (ax1) = plt.subplots(1, 1)#, gridspec_kw={'width_ratios': [1, 3]})
fig.suptitle('Spatial aggregation difference 2017-2010')
#Display DEM
cb1=ax1.imshow(elevDem, extent=grid.extent,cmap=discrete_cmap(10,'cubehelix_r'),alpha=0.5,norm=divnorm)
cbar1=fig.colorbar(cb1, ax=[ax1], location='left')
cbar1.set_label('Elevation [m]')
# Make the norm for difference plotting
divnorm_diff = mcolors.DivergingNorm(vmin=-5, vcenter=0, vmax=5)
#Display 2017 data
plt.scatter(df_2010_2018[df_2010_2018.Track_name.str[:4]=='2017']['lon_3413'],df_2010_2018[df_2010_2018.Track_name.str[:4]=='2017']['lat_3413'],s=0.1,color='#525252',label='2017')
#Display 2010 data
plt.scatter(df_2010_2018[df_2010_2018.Track_name.str[:4]=='2010']['lon_3413'],df_2010_2018[df_2010_2018.Track_name.str[:4]=='2010']['lat_3413'],s=0.1,color='#d9d9d9',label='2010')
#Display the difference between 2017 and 2010 if aggregated data
sc= ax1.scatter(df_spatially_aggregated_2017['avg_lon_3413'],df_spatially_aggregated_2017['avg_lat_3413'],c=diff_to_plot['avg_20m_icecontent'],cmap='seismic_r',norm=divnorm_diff)
cbar=fig.colorbar(sc)
cbar.set_label('Difference in iceslabs thickness', fontsize=15)
plt.legend()

#Allows to open plot in full size directly
figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()

plt.show()


#2018-2010
#this is the 7 years overlaping case study!!
diff_to_plot=df_spatially_aggregated_2018-df_spatially_aggregated_2010

fig, (ax1) = plt.subplots(1, 1)#, gridspec_kw={'width_ratios': [1, 3]})
fig.suptitle('Spatial aggregation difference 2018-2010')
#Display DEM
cb1=ax1.imshow(elevDem, extent=grid.extent,cmap=discrete_cmap(10,'cubehelix_r'),alpha=0.5,norm=divnorm)
cbar1=fig.colorbar(cb1, ax=[ax1], location='left')
cbar1.set_label('Elevation [m]')
# Make the norm for difference plotting
divnorm_diff = mcolors.DivergingNorm(vmin=-5, vcenter=0, vmax=5)
#Display 2018 data
plt.scatter(df_2010_2018[df_2010_2018.Track_name.str[:4]=='2018']['lon_3413'],df_2010_2018[df_2010_2018.Track_name.str[:4]=='2018']['lat_3413'],s=0.1,color='#525252',label='2018')
#Display 2010 data
plt.scatter(df_2010_2018[df_2010_2018.Track_name.str[:4]=='2010']['lon_3413'],df_2010_2018[df_2010_2018.Track_name.str[:4]=='2010']['lat_3413'],s=0.1,color='#d9d9d9',label='2010')
#Display the difference between 2017 and 2010 if aggregated data
sc= ax1.scatter(df_spatially_aggregated_2018['avg_lon_3413'],df_spatially_aggregated_2018['avg_lat_3413'],c=diff_to_plot['avg_20m_icecontent'],cmap='seismic_r',norm=divnorm_diff)
cbar=fig.colorbar(sc)
cbar.set_label('Difference in iceslabs thickness', fontsize=15)
plt.legend()

ax1.set_xlim(-150000,-80000)
ax1.set_ylim(-2460000,-2410000)

#Allows to open plot in full size directly
figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()

plt.show()

'''
##############################################################################
###     This part of code is from OIB1014_iceslabs_script.py and adapted   ###
##############################################################################

import datetime as dt
import statsmodels.api as sm
import datetime
from scipy import stats
from mpl_toolkits.mplot3d import Axes3D

# --------------- Calculate the slope of each aggregated point --------------- #
# Create empty vectors
slope_vect=np.zeros(len(df_2010_2018_spatially_aggregated))
slope_vect[:]=np.nan

rvalue_vect=np.zeros(len(df_2010_2018_spatially_aggregated))
rvalue_vect[:]=np.nan

pvalue_vect=np.zeros(len(df_2010_2018_spatially_aggregated))
pvalue_vect[:]=np.nan

stderr_vect=np.zeros(len(df_2010_2018_spatially_aggregated))
stderr_vect[:]=np.nan

lat_vect=np.zeros(len(df_2010_2018_spatially_aggregated))
lat_vect[:]=np.nan

lon_vect=np.zeros(len(df_2010_2018_spatially_aggregated))
lon_vect[:]=np.nan

key_vect=np.zeros(len(df_2010_2018_spatially_aggregated))
key_vect[:]=np.nan

nb_vect=np.zeros(len(df_2010_2018_spatially_aggregated))
nb_vect[:]=np.nan

#Calculate the slope:
#count for datastorage
count=0
#Loop over the keys
for k in range(0,len(pd.unique(df_2010_2018_spatially_aggregated['key']))):
    #Define the key bounds
    key_bound=str('key=='+str(k))
    
    #Filter the data accoring to key bin
    df_temp=df_2010_2018_spatially_aggregated.query(key_bound)
    
    #Check if there are data:
    check=df_temp['year'].to_numpy()
    if (len(np.unique(check))>1): #there are more than one year of data, growth
    #rate calcultation possible!
        # -------------- Calculate slope using all the data -------------- #
        #Convert x and y to array for the slope calculation
        ice_thickness_array_temp=df_temp['avg_20m_icecontent'].to_numpy()
        year_array_temp=df_temp['year'].to_numpy()

        #Calculate the growth rate - any data point is considered in the
        #growth rate calculation
        slope, intercept, r_value, p_value, std_err = stats.linregress(year_array_temp,ice_thickness_array_temp)

        ##plot the data and fitted slope
        # plt.plot(year_array_temp,ice_thickness_array_temp,'o',label='original data')
        # plt.plot(year_array_temp,intercept+slope*year_array_temp,'r',label='fitted line')
        # plt.legend()
        # plt.show()
        # -------------- Calculate slope using all the data -------------- #

        #Store the results from the linear regression fonction into df_slope
        slope_vect[count]=slope
        rvalue_vect[count]=r_value
        pvalue_vect[count]=p_value
        stderr_vect[count]=std_err
        nb_vect[count]=len(np.unique(check))

        lat_vect[count]=df_temp.avg_lat_3413.iloc[0]
        lon_vect[count]=df_temp.avg_lon_3413.iloc[0]
        key_vect[count]=df_temp.key.iloc[0]

        count=count+1
    
    print(k/len(pd.unique(df_2010_2018_spatially_aggregated['key']))*100,' %')

#Create a dataframe where the slope calculation are stored
column_names = ["lat_3413","lon_3413","slope","key","nb_year","r_value","p_value","std_err"]
df_slope = pd.DataFrame(columns = column_names)
#Store the resulting slope, p_value, r_value, std_err into a dataframe
df_slope.lat_3413=lat_vect
df_slope.lon_3413=lon_vect
df_slope.key=key_vect
df_slope.nb_year=nb_vect
df_slope.slope=slope_vect
df_slope.r_value=rvalue_vect
df_slope.p_value=pvalue_vect
df_slope.std_err=stderr_vect

'''
'''
# --------------------- Build one colormap combining two -------------------- #
#This code is from:
# https://stackoverflow.com/questions/31051488/combining-two-matplotlib-colormaps
# sample the colormaps that you want to use. Use 128 from each so we get 256
# colors in total
colors1 = plt.cm.winter(np.linspace(0., 1, 128))
colors2 = plt.cm.autumn_r(np.linspace(0, 1, 128))

# combine them and build a new colormap
colors = np.vstack((colors1, colors2))
mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)
# --------------------- Build one colormap combining two -------------------- #
'''
'''

# ------ Display slope
fig, (ax1) = plt.subplots(1, 1)#, gridspec_kw={'width_ratios': [1, 3]})
fig.suptitle('Spatial aggregation slope')
#Display DEM
cb1=ax1.imshow(elevDem, extent=grid.extent,cmap=discrete_cmap(10,'cubehelix_r'),alpha=0.5,norm=divnorm)
cbar1=fig.colorbar(cb1, ax=[ax1], location='left')
cbar1.set_label('Elevation [m]')
# Make the norm for slope plotting
divnorm_slope = mcolors.DivergingNorm(vmin=-3, vcenter=0, vmax=3)
#Display 2017-2018 data
plt.scatter(df_2010_2018[df_2010_2018.Track_name.str[:4]=='2017']['lon_3413'],df_2010_2018[df_2010_2018.Track_name.str[:4]=='2017']['lat_3413'],s=0.1,color='#d9d9d9',label='2017-2018')
plt.scatter(df_2010_2018[df_2010_2018.Track_name.str[:4]=='2018']['lon_3413'],df_2010_2018[df_2010_2018.Track_name.str[:4]=='2018']['lat_3413'],s=0.1,color='#d9d9d9')
#Display the aggregation data
sc= ax1.scatter(df_slope['lon_3413'],df_slope['lat_3413'],c=df_slope['slope'],s=df_slope.nb_year*20,cmap='seismic_r',norm=divnorm_slope,label='2010-2018 slope')
cbar=fig.colorbar(sc)
cbar.set_label('Growth rate [m/year]', fontsize=15)
plt.legend()
ax1.set_xlim(-150000,-50000)
ax1.set_ylim(-2540000,-2325000)
plt.show()

#Keep only where more than 2 years of data
df_slope_sup2years=df_slope[df_slope.nb_year>2]

# ------ Display slope for overlaping years >2 
fig, (ax1) = plt.subplots(1, 1)#, gridspec_kw={'width_ratios': [1, 3]})
fig.suptitle('Spatial aggregation slope, overlapping > 2 years')
#Display DEM
cb1=ax1.imshow(elevDem, extent=grid.extent,cmap=discrete_cmap(10,'cubehelix_r'),alpha=0.5,norm=divnorm)
cbar1=fig.colorbar(cb1, ax=[ax1], location='left')
cbar1.set_label('Elevation [m]')
#Display 2017-2018 data
plt.scatter(df_2010_2018[df_2010_2018.Track_name.str[:4]=='2017']['lon_3413'],df_2010_2018[df_2010_2018.Track_name.str[:4]=='2017']['lat_3413'],s=0.1,color='#d9d9d9',label='2017-2018')
plt.scatter(df_2010_2018[df_2010_2018.Track_name.str[:4]=='2018']['lon_3413'],df_2010_2018[df_2010_2018.Track_name.str[:4]=='2018']['lat_3413'],s=0.1,color='#d9d9d9')
#Display the aggregation data
sc= ax1.scatter(df_slope_sup2years['lon_3413'],df_slope_sup2years['lat_3413'],c=df_slope_sup2years['slope'],cmap='seismic_r',norm=divnorm_slope,label='2010-2018 slope, overlapping>2years')
cbar=fig.colorbar(sc)
cbar.set_label('Growth rate [m/year]', fontsize=15)
plt.legend()
ax1.set_xlim(-150000,-50000)
ax1.set_ylim(-2540000,-2325000)
plt.show()

#Keep only statistical significant values
df_slope_sign=df_slope[df_slope.p_value<0.05]

# ------ Display statisitcally significant slope
fig, (ax1) = plt.subplots(1, 1)#, gridspec_kw={'width_ratios': [1, 3]})
fig.suptitle('Spatial aggregation slope - statistically significant')
#Display DEM
cb1=ax1.imshow(elevDem, extent=grid.extent,cmap=discrete_cmap(10,'cubehelix_r'),alpha=0.5,norm=divnorm)
cbar1=fig.colorbar(cb1, ax=[ax1], location='left')
cbar1.set_label('Elevation [m]')
# Make the norm for slope plotting
divnorm_slope = mcolors.DivergingNorm(vmin=-3, vcenter=0, vmax=3)
#Display 2017-2018 data
plt.scatter(df_2010_2018[df_2010_2018.Track_name.str[:4]=='2017']['lon_3413'],df_2010_2018[df_2010_2018.Track_name.str[:4]=='2017']['lat_3413'],s=0.1,color='#d9d9d9',label='2017-2018')
plt.scatter(df_2010_2018[df_2010_2018.Track_name.str[:4]=='2018']['lon_3413'],df_2010_2018[df_2010_2018.Track_name.str[:4]=='2018']['lat_3413'],s=0.1,color='#d9d9d9')
#Display the aggregation data
sc= ax1.scatter(df_slope_sign['lon_3413'],df_slope_sign['lat_3413'],c=df_slope_sign['slope'],s=df_slope_sign.nb_year*20,cmap='seismic_r',norm=divnorm_slope,label='2010-2018 slope, statistically significant')
cbar=fig.colorbar(sc)
cbar.set_label('Growth rate [m/year]', fontsize=15)
plt.legend()
ax1.set_xlim(-150000,-50000)
ax1.set_ylim(-2540000,-2325000)
plt.show()

#Keep only where more than 2 years of data
df_slope_sign_sup2years=df_slope_sign[df_slope_sign.nb_year>2]

# ------ Display slope for overlaping years >2 
fig, (ax1) = plt.subplots(1, 1)#, gridspec_kw={'width_ratios': [1, 3]})
fig.suptitle('Spatial aggregation slope, statisitcally significant, overlapping > 2 years')
#Display DEM
cb1=ax1.imshow(elevDem, extent=grid.extent,cmap=discrete_cmap(10,'cubehelix_r'),alpha=0.5,norm=divnorm)
cbar1=fig.colorbar(cb1, ax=[ax1], location='left')
cbar1.set_label('Elevation [m]')
#Display 2017-2018 data
plt.scatter(df_2010_2018[df_2010_2018.Track_name.str[:4]=='2017']['lon_3413'],df_2010_2018[df_2010_2018.Track_name.str[:4]=='2017']['lat_3413'],s=0.1,color='#d9d9d9',label='2017-2018')
plt.scatter(df_2010_2018[df_2010_2018.Track_name.str[:4]=='2018']['lon_3413'],df_2010_2018[df_2010_2018.Track_name.str[:4]=='2018']['lat_3413'],s=0.1,color='#d9d9d9')
#Display the aggregation data
sc= ax1.scatter(df_slope_sign_sup2years['lon_3413'],df_slope_sign_sup2years['lat_3413'],c=df_slope_sign_sup2years['slope'],cmap='seismic_r',norm=divnorm_slope,label='2010-2018 slope, statistically significant, overlapping>2years')
cbar=fig.colorbar(sc)
cbar.set_label('Growth rate [m/year]', fontsize=15)
plt.legend()
ax1.set_xlim(-150000,-50000)
ax1.set_ylim(-2540000,-2325000)
plt.show()

##############################################################################
###     This part of code is from OIB1014_iceslabs_script.py and adapted   ###
##############################################################################

#Display where overlapping = 7 years.
df_slope_7years=df_slope[df_slope.nb_year==7]

# Make the norm for slope plotting
divnorm_slope = mcolors.DivergingNorm(vmin=-2, vcenter=0, vmax=2)
# ------ Display slope for overlaping years >2 
fig, (ax1) = plt.subplots(1, 1)#, gridspec_kw={'width_ratios': [1, 3]})
fig.suptitle('Spatial aggregation slope, overlapping = 7 years')
#Display DEM
cb1=ax1.imshow(elevDem, extent=grid.extent,cmap=discrete_cmap(10,'cubehelix_r'),alpha=0.5,norm=divnorm)
cbar1=fig.colorbar(cb1, ax=[ax1], location='left')
cbar1.set_label('Elevation [m]')
#Display 2017-2018 data
plt.scatter(df_2010_2018[df_2010_2018.Track_name.str[:4]=='2017']['lon_3413'],df_2010_2018[df_2010_2018.Track_name.str[:4]=='2017']['lat_3413'],s=0.1,color='#d9d9d9',label='2017-2018')
plt.scatter(df_2010_2018[df_2010_2018.Track_name.str[:4]=='2018']['lon_3413'],df_2010_2018[df_2010_2018.Track_name.str[:4]=='2018']['lat_3413'],s=0.1,color='#d9d9d9')
#Display the 7 years aggregation data
sc= ax1.scatter(df_slope_7years['lon_3413'],df_slope_7years['lat_3413'],c=df_slope_7years['slope'],cmap='seismic_r',norm=divnorm_slope,label='2010-2018 slope, overlapping = 7years')
cbar=fig.colorbar(sc)
cbar.set_label('Growth rate [m/year]', fontsize=15)
plt.legend()
ax1.set_xlim(-115000,-80000)
ax1.set_ylim(-2470000,-2440000)
plt.show()

'''








