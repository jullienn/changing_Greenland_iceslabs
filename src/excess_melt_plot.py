# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 14:29:11 2021

@author: JullienN

Read and plot the excess melt calculation made by Andrew Tedstone
"""

#Import librairies

#Import packages

import matplotlib.pyplot as plt
import rasterio as rio
import geopandas as gpd
import os
import numpy as np
import xarray as xr

#Define path of the working directory
global_path='C:/Users/jullienn/Documents/working_environment/'

#Define the desired year to plot
wanted_year='2006'

##########################################################################
###                      Load Greenland DEM and contours               ###
##########################################################################

dem_filename = 'C:/Users/jullienn/Documents/working_environment/greenland_topo_data/elevations/greenland_dem_mosaic_100m_v3.0.tif'
contours_filename = 'C:/Users/jullienn/Documents/working_environment/greenland_topo_data/elevations/greenland_contours_100m_v3.0.shp'

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
###                          Load excess melt data   	               ###
##########################################################################

#Define path
path='C:/Users/jullienn/Documents/working_environment/excess_melt/'

#Load the data
data_path = path+'excess_melt_mbyear.nc'
DS = xr.open_dataset(data_path)

#Extract coordinates
lat=DS.x.data
lon=DS.y.data

#Load melt data
melt_data= DS.M_e

#Select a single year
melt_year = melt_data.sel(time=wanted_year)
melt_year_np = melt_year.values

melt_year_plot=np.asarray(melt_year_np)
melt_year_plot=melt_year_plot[0,0,:,:,0]

##########################################################################
###                          Load excess melt data   	               ###
##########################################################################

##########################################################################
###                          Plot excess melt data   	               ###
##########################################################################

##Where excess melt is above 500 mm w.e./year, equal to 1, else equal to 0
#melt_year_plot_binary=melt_year_plot[melt_year_plot<500]=0
#melt_year_plot_binary=melt_year_plot[~(melt_year_plot<500)]=1

#Plot dem and contours elevation
plt.figure()
ax = plt.subplot(111)
dem_extent = (dem_bounds[0], dem_bounds[2], dem_bounds[1], dem_bounds[3])
plt.imshow(np.squeeze(np.flipud(melt_year_plot)), extent=dem_extent,cmap='plasma_r')
plt.colorbar(label='Excess melt [mm w.e./year]')
plt.clim(0,1000)
contours.plot(ax=ax, edgecolor='black')
plt.title('Excess melt plot, year: '+wanted_year)

##Create the figure name
#fig_name=[]
#fig_name='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/2002_2003_data_localisation/'+indiv_file+'.png'
                    
##Save the figure
#pyplot.savefig(fig_name)
#pyplot.clf()



























##########################################################################
###                          Plot excess melt data   	               ###
##########################################################################





#da = DS.T_z.sel(z=500,method='nearest').dropna(dim='time') - 273.15  # or .ffill(dim='time')

