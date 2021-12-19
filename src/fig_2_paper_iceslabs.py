# -*- coding: utf-8 -*-
"""
Created on Sun Dec 19 12:14:06 2021

@author: jullienn
"""


#######################################################################
###       Thickening analysis using spatially aggregated files      ###
#######################################################################   

###     This is from iceslabs_20102018_thickening_analysis.py       ###

#Import librairies
import datetime
from scipy import spatial

#Load the spatial aggregated data. All the points within a radius of 100m are averaged
path='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2010_2018/excel_spatial_aggreation_and_other/final_excel/prob00/'
df_2010_2018_spatially_aggregated = pd.read_csv(path+'jullien_etal_20102018_spatial_aggregation_grid_1000_prob00.csv',delimiter=';',decimal=',')

#Load all 2010-2018 data without spatial aggregation
df_2010_2018_csv = pd.read_csv(path+'Ice_Layer_Output_Thicknesses_2010_2018_jullienetal2021_prob00.csv',delimiter=',',decimal='.')
#Transform the coordinated from WGS84 to EPSG:3413
transformer = Transformer.from_crs("EPSG:4326", "EPSG:3413", always_xy=True)
points=transformer.transform(np.asarray(df_2010_2018_csv["lon"]),np.asarray(df_2010_2018_csv["lat"]))

#Store lat/lon in 3413
df_2010_2018_csv['lon_3413']=points[0]
df_2010_2018_csv['lat_3413']=points[1]

#Visualize the spatial aggregation process
########################## Load GrIS elevation ##########################
#Open the DEM
grid = Grid.from_raster("C:/Users/jullienn/switchdrive/Private/research/backup_Aglaja/working_environment/greenland_topo_data/elevations/greenland_dem_mosaic_100m_v3.0.tif",data_name='dem')
#Minnor slicing on borders to enhance colorbars
elevDem=grid.dem[:-1,:-1]              
#Scale the colormap
#divnorm = mcolors.DivergingNorm(vmin=0, vcenter=1250, vmax=2500)
########################## Load GrIS elevation ##########################

'''
fig, (ax1) = plt.subplots(1, 1)#, gridspec_kw={'width_ratios': [1, 3]})
fig.suptitle('Spatial aggregation illustration')
#Display DEM
cb1=ax1.imshow(elevDem, extent=grid.extent,cmap=discrete_cmap(10,'cubehelix_r'),alpha=0.5,norm=divnorm)
cbar1=fig.colorbar(cb1, ax=[ax1], location='left')
cbar1.set_label('Elevation [m]')
#Display 2010-2018 data
plt.scatter(df_2010_2018_csv[df_2010_2018_csv.Track_name.str[:4]=='2010']['lon_3413'],df_2010_2018_csv[df_2010_2018_csv.Track_name.str[:4]=='2010']['lat_3413'],s=0.1,color='#525252',label='2010-2014')
plt.scatter(df_2010_2018_csv[df_2010_2018_csv.Track_name.str[:4]=='2011']['lon_3413'],df_2010_2018_csv[df_2010_2018_csv.Track_name.str[:4]=='2011']['lat_3413'],s=0.1,color='#525252')
plt.scatter(df_2010_2018_csv[df_2010_2018_csv.Track_name.str[:4]=='2012']['lon_3413'],df_2010_2018_csv[df_2010_2018_csv.Track_name.str[:4]=='2012']['lat_3413'],s=0.1,color='#525252')
plt.scatter(df_2010_2018_csv[df_2010_2018_csv.Track_name.str[:4]=='2013']['lon_3413'],df_2010_2018_csv[df_2010_2018_csv.Track_name.str[:4]=='2013']['lat_3413'],s=0.1,color='#525252')
plt.scatter(df_2010_2018_csv[df_2010_2018_csv.Track_name.str[:4]=='2014']['lon_3413'],df_2010_2018_csv[df_2010_2018_csv.Track_name.str[:4]=='2014']['lat_3413'],s=0.1,color='#525252')
plt.scatter(df_2010_2018_csv[df_2010_2018_csv.Track_name.str[:4]=='2017']['lon_3413'],df_2010_2018_csv[df_2010_2018_csv.Track_name.str[:4]=='2017']['lat_3413'],s=0.1,color='#d9d9d9',label='2017-2018')
plt.scatter(df_2010_2018_csv[df_2010_2018_csv.Track_name.str[:4]=='2018']['lon_3413'],df_2010_2018_csv[df_2010_2018_csv.Track_name.str[:4]=='2018']['lat_3413'],s=0.1,color='#d9d9d9')
#Display the aggregation data
plt.scatter(df_2010_2018_spatially_aggregated['avg_lon_3413'],df_2010_2018_spatially_aggregated['avg_lat_3413'],color='green',label='2010-2018 spatially aggregated')
plt.legend()
plt.show()
'''

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



list_high_end=list(['2002-2003','2010','2011-2012','2013-2014','2017-2018'])
plot_thickness_high_end(df_2010_2018,df_spatially_aggregated_2017,df_spatially_aggregated_2010,elevDem,grid,slice_lon_summary,lat_slices,list_high_end)

'''
list_high_end=list(['2002-2003','2011-2012','2017-2018'])
plot_thickness_high_end(df_2010_2018,df_spatially_aggregated_2017,df_spatially_aggregated_2011,elevDem,grid,slice_lon_summary,lat_slices,list_high_end)

list_high_end=list(['2002-2003','2011-2012','2017-2018'])
plot_thickness_high_end(df_2010_2018,df_spatially_aggregated_2017,df_spatially_aggregated_2012,elevDem,grid,slice_lon_summary,lat_slices,list_high_end)

list_high_end=list(['2002-2003','2013-2014','2017-2018'])
plot_thickness_high_end(df_2010_2018,df_spatially_aggregated_2017,df_spatially_aggregated_2013,elevDem,grid,slice_lon_summary,lat_slices,list_high_end)

list_high_end=list(['2002-2003','2013-2014','2017-2018'])
plot_thickness_high_end(df_2010_2018,df_spatially_aggregated_2017,df_spatially_aggregated_2014,elevDem,grid,slice_lon_summary,lat_slices,list_high_end)

list_high_end=list(['2002-2003','2013-2014','2011-2012'])
plot_thickness_high_end(df_2010_2018,df_spatially_aggregated_2013,df_spatially_aggregated_2012,elevDem,grid,slice_lon_summary,lat_slices,list_high_end)

list_high_end=list(['2002-2003','2011-2012','2010'])
plot_thickness_high_end(df_2010_2018,df_spatially_aggregated_2011,df_spatially_aggregated_2010,elevDem,grid,slice_lon_summary,lat_slices,list_high_end)
'''
###     This is from iceslabs_20102018_thickening_analysis.py       ###

#######################################################################
###       Thickening analysis using spatially aggregated files      ###
#######################################################################