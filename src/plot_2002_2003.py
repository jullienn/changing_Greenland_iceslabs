# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 17:18:07 2021

@author: JullienN
"""
#Import packages
import rasterio
from rasterio.plot import show
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pandas as pd
from os import listdir
from os.path import isfile, join
import pickle
from pysheds.grid import Grid
import pdb
import numpy as np
from pyproj import Transformer

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


pdb.set_trace()
########################## Load GrIS elevation ##########################
#Open the DEM
grid = Grid.from_raster("C:/Users/jullienn/Documents/working_environment/greenland_topo_data/elevations/greenland_dem_mosaic_100m_v3.0.tif",data_name='dem')
#Minnor slicing on borders to enhance colorbars
elevDem=grid.dem[:-1,:-1]              
#Scale the colormap
divnorm = mcolors.DivergingNorm(vmin=0, vcenter=1250, vmax=2500)
########################## Load GrIS elevation ##########################

################# Load 2002-2003 flightlines coordinates ################
path_data='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/icelens_identification'

#Open the file and read it
f_flightlines = open(path_data+'/metadata_coord_2002_2003', "rb")
all_2002_3_flightlines = pickle.load(f_flightlines)
f_flightlines.close()

lat_all=[]
lon_all=[]

for year in list(all_2002_3_flightlines.keys()):
    for days in list(all_2002_3_flightlines[year].keys()):
        for indiv_file in list(all_2002_3_flightlines[year][days].keys()):
            if (indiv_file[0:7]=='quality'):
                continue
            else:
                lat_all=np.append(lat_all,all_2002_3_flightlines[year][days][indiv_file][0])
                lon_all=np.append(lon_all,all_2002_3_flightlines[year][days][indiv_file][1])

################# Load 2002-2003 flightlines coordinates ################

################### Load 2002-2003 ice lenses location ##################
#Open the file and read it
f_icelens_flightlines = open(path_data+'/metadata_coord_icelens_2002_2003', "rb")
icelens_2002_3_flightlines = pickle.load(f_icelens_flightlines)
f_icelens_flightlines.close()

lat_icelens=[]
lon_icelens=[]

for year in list(icelens_2002_3_flightlines.keys()):
    for days in list(icelens_2002_3_flightlines[year].keys()):
        for indiv_file in list(icelens_2002_3_flightlines[year][days].keys()):
            print(indiv_file)
            if (indiv_file[0:7]=='quality'):
                print('Quality file, continue')
                continue
            elif (not(bool(icelens_2002_3_flightlines[year][days][indiv_file]))):
                print('No ice lens, continue')
                continue
            else:
                lat_icelens=np.append(lat_icelens,icelens_2002_3_flightlines[year][days][indiv_file][0])
                lon_icelens=np.append(lon_icelens,icelens_2002_3_flightlines[year][days][indiv_file][1])
################### Load 2002-2003 ice lenses location ##################

################### Load 2010-2014 ice slabs location ##################
#Load the data
filename_MacFerrin= 'C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/icelens_identification/indiv_traces_icelenses/MacFerrin_etal2019_iceslabs.xlsx'
#Read MacFerrin data thanks to https://stackoverflow.com/questions/65254535/xlrd-biffh-xlrderror-excel-xlsx-file-not-supported
df_MacFerrin = pd.read_excel(filename_MacFerrin,engine='openpyxl')

#Transform the coordinated from WGS84 to EPSG:3413
#Example from: https://pyproj4.github.io/pyproj/stable/examples.html
transformer = Transformer.from_crs("EPSG:4326", "EPSG:3413", always_xy=True)
points=transformer.transform(np.array(df_MacFerrin.lon),np.array(df_MacFerrin.lat))

lon_3413_MacFerrin=points[0]
lat_3413_MacFerrin=points[1]

################### Load 2010-2014 ice slabs location ##################

################################### Plot ##################################
plt.rcParams.update({'font.size': 5})
fig, (ax1) = plt.subplots(1, 1)#, gridspec_kw={'width_ratios': [1, 3]})
fig.suptitle('Insert title')
cb1=ax1.imshow(elevDem, extent=grid.extent,cmap=discrete_cmap(10,'cubehelix_r'),norm=divnorm)
cbar1=fig.colorbar(cb1, ax=[ax1], location='left')
cbar1.set_label('Elevation [m]', fontsize=5)
#ax1.grid()
ax1.set_title('Ice lenses and slabs mapping',fontsize=5)

#Plot all the 2002-2003 flightlines
ax1.scatter(lon_all, lat_all,s=0.1,color='lightgrey')

#Plot all the 2002-2003 icelenses
ax1.scatter(lon_3413_MacFerrin, lat_3413_MacFerrin,s=0.1,color='slateblue')

#Plot all the 2002-2003 icelenses
ax1.scatter(lon_icelens, lat_icelens,s=0.1,color='blue')
################################### Plot ##################################
