# -*- coding: utf-8 -*-
"""
Created on Thu Sep 16 09:20:00 2021

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

#Load all 2010-2018 data
#path='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2010_2018/excel_files/'
path='/home/jullienn/data/iceslabs/'
df_2010_2018 = pd.read_csv(path+'Ice_Layer_Output_Thicknesses_2010_2018_jullienetal2021.csv',delimiter=',',decimal='.')

#Select individual years and save corresponding files
df_2010=df_2010_2018[df_2010_2018.Track_name.str[:4]=='2010']
df_2011=df_2010_2018[df_2010_2018.Track_name.str[:4]=='2011']
df_2012=df_2010_2018[df_2010_2018.Track_name.str[:4]=='2012']
df_2013=df_2010_2018[df_2010_2018.Track_name.str[:4]=='2013']
df_2014=df_2010_2018[df_2010_2018.Track_name.str[:4]=='2014']
df_2017=df_2010_2018[df_2010_2018.Track_name.str[:4]=='2017']
df_2018=df_2010_2018[df_2010_2018.Track_name.str[:4]=='2018']

df_2010.to_csv(path_or_buf=path+'Ice_Layer_Output_Thicknesses_2010_jullienetal2021.csv',header=True)
df_2011.to_csv(path_or_buf=path+'Ice_Layer_Output_Thicknesses_2011_jullienetal2021.csv',header=True)
df_2012.to_csv(path_or_buf=path+'Ice_Layer_Output_Thicknesses_2012_jullienetal2021.csv',header=True)
df_2013.to_csv(path_or_buf=path+'Ice_Layer_Output_Thicknesses_2013_jullienetal2021.csv',header=True)
df_2014.to_csv(path_or_buf=path+'Ice_Layer_Output_Thicknesses_2014_jullienetal2021.csv',header=True)
df_2017.to_csv(path_or_buf=path+'Ice_Layer_Output_Thicknesses_2017_jullienetal2021.csv',header=True)
df_2018.to_csv(path_or_buf=path+'Ice_Layer_Output_Thicknesses_2018_jullienetal2021.csv',header=True)


#Create a year column
year=df_2010_2018.Track_name.str[:4]
df_2010_2018['year']=year.astype(str).astype(int)

#Transform the coordinated from WGS84 to EPSG:3413
transformer = Transformer.from_crs("EPSG:4326", "EPSG:3413", always_xy=True)
points=transformer.transform(np.asarray(df_2010_2018["lon"]),np.asarray(df_2010_2018["lat"]))

#Store lat/lon in 3413
df_2010_2018['lon_3413']=points[0]
df_2010_2018['lat_3413']=points[1]

#Visualize the spatial aggregation process
########################## Load GrIS elevation ##########################
#Open the DEM
#path_raster='C:/Users/jullienn/switchdrive/Private/research/backup_Aglaja/working_environment/greenland_topo_data/elevations/'
path_raster='/home/jullienn/data/iceslabs/'
grid = Grid.from_raster(path_raster+"greenland_dem_mosaic_100m_v3.0.tif",data_name='dem')
#Minnor slicing on borders to enhance colorbars
elevDem=grid.dem[:-1,:-1]              
#Scale the colormap
divnorm = mcolors.DivergingNorm(vmin=0, vcenter=1250, vmax=2500)
########################## Load GrIS elevation ##########################
'''
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
plt.legend()
plt.show()
'''

#Create the grid
a=1000 #length of one side of the cell
b=1000 #lenght of the other side of the cell

lat_vector=np.arange(np.min(df_2010_2018['lat_3413']),np.max(df_2010_2018['lat_3413']),a)
lon_vector=np.arange(np.min(df_2010_2018['lon_3413']),np.max(df_2010_2018['lon_3413']),b)

grid_for_aggregation=np.meshgrid(lon_vector,lat_vector)
points=np.c_[grid_for_aggregation[0].ravel(), grid_for_aggregation[1].ravel()]

#Radius to search is: r = sqrt(a^2+b^2)/2. We round up to the next integer to be conservative
radius_search=np.ceil(np.sqrt(a*a+b*b)/2)
######################################################################
###         Begin from spatial_aggregation_2010_2018.py            ###
#Define coordinates tuples
XY_iceslabs = np.array((df_2010_2018['lon_3413'],df_2010_2018['lat_3413'])).swapaxes(0,1)

#Define the look up tree
tree = spatial.cKDTree(XY_iceslabs)

print('Build the lookup tree')
#This is from : https://medium.com/nam-r/10-essential-operations-for-spatial-data-in-python-4603d933bdda, point 7
#Find all the point that lies within the circle of radius r centered around the point of reference
neigh_list=[]
for count, value in enumerate(points):
    neigh_list.append(tree.query_ball_point(value,r=radius_search))
    #print(count/len(points)*100,'%')

#Create the vectors to store the aggregated data
avg_20m_icecontent=np.zeros(len(neigh_list))
avg_20m_icecontent[:]=np.nan

std_20m_icecontent=np.zeros(len(neigh_list))
std_20m_icecontent[:]=np.nan

avg_lat_3413=np.zeros(len(neigh_list))
avg_lat_3413[:]=np.nan

avg_lon_3413=np.zeros(len(neigh_list))
avg_lon_3413[:]=np.nan

year_data=np.zeros(len(neigh_list))
year_data[:]=np.nan

keys=np.zeros(len(neigh_list))
keys[:]=np.nan

#Set key to zero
key = 0

#initialize count to 0
count=0

#Initialise ite for lat/lon of grid retrieval
ite=0

#Initialize the verification vector
verification=np.nan

print('Make the correspondance')

for index_list in neigh_list:

    #retreive all data which index belong to index_list
    df_temp=df_2010_2018.iloc[index_list]
    
    #print(ite/len(neigh_list)*100,' %')
    
    if (len(df_temp) == 0):
        #no data at this grid point, continue
        #Update ite
        ite=ite+1
        continue
    else:
        #data in this grid point, average data
        #Loop over the years for aggregation
        for year in list(np.unique(df_temp['year'])):
            #average the 20m ice content for each year
            ice_avg=np.average(df_temp[df_temp['year']==year]['20m_ice_content_m'])
            
            #calculate std of 20m ice content for each year
            std_avg=np.std(df_temp[df_temp['year']==year]['20m_ice_content_m'])
            
            #take the lat/lon of the gris point
            lat_avg=points[ite][1]
            lon_avg=points[ite][0]
            
            #Store the results
            avg_20m_icecontent[count]=ice_avg
            std_20m_icecontent[count]=std_avg
            avg_lat_3413[count]=lat_avg
            avg_lon_3413[count]=lon_avg
            year_data[count]=year
            keys[count]=key
            
            #Update count
            count=count+1
            
            #List for verification
            verification=np.append(verification,index_list)
        
        #Update ite
        ite=ite+1
        
        #Update key
        key=key+1


#Delete all the NaNs in verification vector
verification_without_nans=verification[~np.isnan(verification)]
#Compare the length of verification_without_nans with length of df_2010_2018
data_not_aggregated=(len(df_2010_2018)-len(np.unique(verification_without_nans)))/len(df_2010_2018)
#Display in command window
print(data_not_aggregated,' % of the data have not been aggreated')

#Remove nans
avg_20m_icecontent=avg_20m_icecontent[~np.isnan(avg_20m_icecontent)]
std_20m_icecontent=std_20m_icecontent[~np.isnan(std_20m_icecontent)]
avg_lat_3413=avg_lat_3413[~np.isnan(avg_lat_3413)]
avg_lon_3413=avg_lon_3413[~np.isnan(avg_lon_3413)]
year_data=year_data[~np.isnan(year_data)]
keys=keys[~np.isnan(keys)]

#Create a pd dataframe
df_spatial_aggregation=pd.DataFrame({'avg_20m_icecontent':avg_20m_icecontent,
                                     'std_20m_icecontent':std_20m_icecontent,
                                     'avg_lat_3413':avg_lat_3413,
                                     'avg_lon_3413':avg_lon_3413,
                                     'year':year_data,
                                     'key':keys})


#Save df as excel file
filename_to_save='jullien_etal_20102018_spatial_aggregation_grid_'+str(a)+'.xlsx'
df_spatial_aggregation.to_excel(path+filename_to_save)

print(' --- End of processing ---')

###         End from spatial_aggregation_2010_2018.py            ###
######################################################################



