# -*- coding: utf-8 -*-
"""
Created on Mon Sep 13 10:26:33 2021

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

#Custom radius
custom_radius=1000
#With 1000, 0.025% of the data do not belong to any aggregation point
#With 500, 0.02% of the data do not belong to any aggregation point
#With 100, 0.016% of the data do not belong to any aggregation point

#Load the excel file
path='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2010_2018/excel_files/'
df_2010_2018 = pd.read_csv(path+'jullienetal_20102018.csv',delimiter=';',decimal=',')

#Create a year column
year=df_2010_2018.Track_name.str[:4]
df_2010_2018['year']=year.astype(str).astype(int)

#Transform the coordinated from WGS84 to EPSG:3413
#Example from: https://pyproj4.github.io/pyproj/stable/examples.html
transformer = Transformer.from_crs("EPSG:4326", "EPSG:3413", always_xy=True)
points=transformer.transform(np.asarray(df_2010_2018["lon"]),np.asarray(df_2010_2018["lat"]))

#Store lat/lon in 3413
df_2010_2018['lon_3413']=points[0]
df_2010_2018['lat_3413']=points[1]

#Define coordinates tuples
XY_iceslabs = np.array((df_2010_2018['lon_3413'],df_2010_2018['lat_3413'])).swapaxes(0,1)

#Define the look up tree
tree = spatial.cKDTree(XY_iceslabs)

#This is from : https://medium.com/nam-r/10-essential-operations-for-spatial-data-in-python-4603d933bdda, point 7
#Find all the point that lies within the circle of radius r centered around the point of reference
neigh_list=[]
for count, value in enumerate(XY_iceslabs):
    neigh_list.append(tree.query_ball_point(value,r=custom_radius))

#Create a vector of nan of length of neigh_list
keys=np.empty(len(neigh_list))
keys[:]=np.nan

#Initialise the list for storing
neigh_list_unique=[]

#Initialize the verification vector
verification=np.nan

#For step identificastion
step=np.empty(len(neigh_list))
step[:]=np.nan

i=0
for index in neigh_list:
    step[i]=len(index)
    i=i+1

#Keep only unique keys
for i in range(0,len(neigh_list),int(np.quantile(step,0.1))):
    index = neigh_list[i]
    #Copy the list of index and transform into numpy array
    index_to_copy=np.asarray(index).astype(float)
    
    #Loop over the indices of this list of index
    for indiv_index in index:
        #pdb.set_trace()
        if (len(keys[keys==indiv_index])==1):
            #Key have already been encountered, delete it from
            index_to_copy[index_to_copy==indiv_index]=np.nan
        else:
            #Key have not been encountered yet, store it into the keys vector
            keys[indiv_index]=indiv_index
     
    #pdb.set_trace()
    #Transform index_to_copy into a list and store it
    neigh_list_unique.append(list(index_to_copy))
    
    #List for verification
    verification=np.append(verification,index_to_copy)
    
    print(i/len(neigh_list)*100,' %')

#Delete all the NaNs in verification vector
verification_without_nans=verification[~np.isnan(verification)]
#Compare the length of verification_without_nans with length of df_2010_2018
data_not_aggregated=(len(df_2010_2018)-len(verification_without_nans))/len(df_2010_2018)
#Display in command window
print(data_not_aggregated,' % of the data have not been aggreated')

'''
#Create the dictionnary to store the data
dict_aggregation={k: {} for k in np.arange(0,len(neigh_list_unique))}
'''

#Create the vectors to store the data
avg_20m_icecontent=np.zeros(len(neigh_list))
avg_20m_icecontent[:]=np.nan

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

#Loop over neigh_list_unique and average
for i in range(0,len(neigh_list_unique)):
    
    #Remove NaNs if any
    index_for_aggregation=np.asarray(neigh_list_unique[i])
    index_for_aggregation=index_for_aggregation[~np.isnan(index_for_aggregation)]
    
    '''
    #Continue building dictionnary
    dict_aggregation[i]={k: {} for k in list(df_2010_2018['year'].unique())}
    '''
    
    #Retreive all data belonging to this aggregation point
    df_temp=df_2010_2018.iloc[index_for_aggregation]
    
    #Loop over the years for aggregation
    for year in list(np.unique(df_temp['year'])):
        #average the 20m ice content for each year
        ice_avg=np.average(df_temp[df_temp['year']==year]['20m_ice_content_m'])
        #average the lat/lon of THE WHOLE (no matter the year)
        lat_avg=np.average(df_temp['lat_3413'])
        lon_avg=np.average(df_temp['lon_3413'])
        
        #Store the results
        avg_20m_icecontent[count]=ice_avg
        avg_lat_3413[count]=lat_avg
        avg_lon_3413[count]=lon_avg
        year_data[count]=year
        keys[count]=key
        
        #Update count
        count=count+1
        
        '''
        dict_aggregation[i][year]=np.array([ice_avg,
                                            lat_avg,
                                            lon_avg,
                                            key])
        '''
        
    print(i/len(neigh_list_unique)*100,' %')
    
    #Update key
    key=key+1
        
#Remove nans
avg_20m_icecontent=avg_20m_icecontent[~np.isnan(avg_20m_icecontent)]
avg_lat_3413=avg_lat_3413[~np.isnan(avg_lat_3413)]
avg_lon_3413=avg_lon_3413[~np.isnan(avg_lon_3413)]
year_data=year_data[~np.isnan(year_data)]
keys=keys[~np.isnan(keys)]

#Create a pd dataframe
df_spatial_aggregation=pd.DataFrame({'avg_20m_icecontent':avg_20m_icecontent,
                                     'avg_lat_3413':avg_lat_3413,
                                     'avg_lon_3413':avg_lon_3413,
                                     'year':year_data,
                                     'key':keys})

'''
#Save df as excel file
filename_to_save='jullien_etal_20102018_spatial_aggregation'+str(custom_radius)+'.xlsx'
df_spatial_aggregation.to_excel('C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2010_2018/excel_files/'+filename_to_save)
'''
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

#Select the index to show
index_to_show=1

#Load raw data
index_toshow=np.asarray(neigh_list_unique[index_to_show])
index_toshow=index_toshow[~np.isnan(index_toshow)]
df_toshow=df_2010_2018.iloc[index_toshow]

#Load aggregated data
df_agg=df_spatial_aggregation[df_spatial_aggregation['key']==index_to_show]

#Display the whole 2010-2018 data
plt.scatter(df_2010_2018['lon_3413'],df_2010_2018['lat_3413'],color='blue')

#Display the raw data
plt.scatter(df_toshow['lon_3413'],df_toshow['lat_3413'],color='red')

#Display the aggregation data
plt.scatter(df_agg['avg_lon_3413'],df_agg['avg_lat_3413'],color='green')
plt.show()


#Display the spatial aggregated dataset
fig, (ax1) = plt.subplots(1, 1)#, gridspec_kw={'width_ratios': [1, 3]})
fig.suptitle('Spatial aggregation dataset')

#Display DEM
cb1=ax1.imshow(elevDem, extent=grid.extent,cmap=discrete_cmap(10,'cubehelix_r'),alpha=0.5,norm=divnorm)
cbar1=fig.colorbar(cb1, ax=[ax1], location='left')
cbar1.set_label('Elevation [m]')

#Display the whole 2010-2018 data
plt.scatter(df_2010_2018['lon_3413'],df_2010_2018['lat_3413'],color='blue')

#Display the 2010-2018 spatially aggregated dataset
plt.scatter(df_spatial_aggregation['avg_lon_3413'],df_spatial_aggregation['avg_lat_3413'],color='green')

plt.show()

'''
#Create a grid
X_grid=np.arange(np.min(XY_iceslabs[:,0]),np.max(XY_iceslabs[:,0])+899300,100)
Y_grid=np.arange(np.min(XY_iceslabs[:,1]),np.max(XY_iceslabs[:,1]),100)

XY_grid=np.array((X_grid,Y_grid)).swapaxes(0,1)

#Create a lookup tree associated with this grid
tree_grid = spatial.cKDTree(XY_grid)





#execute lookup tree
distances, indices = tree_grid.query(XY_iceslabs, k=[2], eps=10, p=2, distance_upper_bound=100)

#Remove the last index
indices[indices==len(indices)]=0

#Flatten the vector
indices=np.ndarray.flatten(indices)

#Plot all the data
plt.plot(XY_iceslabs[:,0], XY_iceslabs[:,1], '.')

#Test only with a subset
indices_subset=indices[0:50000]

#Loop over subset to check results
for index in indices_subset:
    
    nearby_points = XY_iceslabs[index,:]
    print(index)
    print(nearby_points)
    plt.plot(nearby_points[0],nearby_points[1], '.')
    plt.show()
    pdb.set_trace()
    
plt.show()





pdb.set_trace()

#search_radius: the radius of points to return, shall broadcast to the length of x.
search_radius=100

#Epsilon: Approximate search. Branches of the tree are not explored if their
#nearest points are further than r / (1 + eps), and branches are added in bulk
#if their furthest points are nearer than r * (1 + eps).

results = tree.query_ball_point(XY_iceslabs, r=100, p=2, eps=0)
pdb.set_trace()




plt.plot(XY_iceslabs[:,0], XY_iceslabs[:,1], '.')
for results in tree.query_ball_point(XY_iceslabs[0:10000:,], search_radius, p=2):
    nearby_points = XY_iceslabs[results]
    plt.plot(nearby_points[:,0], nearby_points[:,1], 'o')





# Check each point against tree
# Might be possible to remove the loop here
distances, indices = tree.query(np.array([data_iceslabs["lon"], data_iceslabs["lat"]]).swapaxes(0,1),
                                k=1,
                                distance_upper_bound=float(config.get('filtering','search_radius_m')),
                                eps=float(config.get('filtering', 'search_radius_epsilon')))

rlims['dist'] = distances
rlims['i'] = indices

rlims_keep = rlims[np.isfinite(rlims.dist)]



from scipy import spatial
x, y = np.mgrid[0:4, 0:4]
points = np.c_[x.ravel(), y.ravel()]
tree = spatial.cKDTree(points)
tree.query_ball_point([2, 0], 1)

import matplotlib.pyplot as plt
points = np.asarray(points)
plt.plot(points[:,0], points[:,1], '.')
for results in tree.query_ball_point(([2, 0], [3, 3]), 0.5):
    nearby_points = points[results]
    plt.plot(nearby_points[:,0], nearby_points[:,1], 'o')
plt.margins(0.1, 0.1)
plt.show()
'''