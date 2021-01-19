# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 14:29:11 2021

@author: JullienN

Read and plot the excess melt calculation made by Andrew Tedstone
"""

#Import librairies

#Import packages

import xarray as xr
import rasterio
from rasterio.plot import show
from matplotlib import pyplot
import numpy as np
from pysheds.grid import Grid

from osgeo import gdal
import elevation 

#Define path of the working directory
global_path='C:/Users/jullienn/Documents/working_environment/'

#Define the desired year to plot
wanted_year='2006'

##########################################################################
###                   Load Greenland satellite image	and DEM            ###
##########################################################################

#Load satellite image
path_image=global_path+'greenland_topo_data/satellite_image/Greenland_natural_90m.tif'

#im = rasterio.open(path_image)
#show(im)
#show(im.read(), transform=im.transform)

#Load the DEM
grid = Grid.from_raster("C:/Users/jullienn/Documents/working_environment/greenland_topo_data/elevations/greenland_dem_mosaic_100m_v3.0.tif",data_name='dem')
#Minnor slicing on borders to enhance colorbars
elevDem=grid.dem[:-1,:-1]              

##########################################################################
###                      Load Greenland satellite image	               ###
##########################################################################

##########################################################################
###                          Compute DEM contour	                       ###
##########################################################################
#Following the procedure from: https://www.earthdatascience.org/tutorials/visualize-digital-elevation-model-contours-matplotlib/

#1. Read the data
filename="C:/Users/jullienn/Documents/working_environment/greenland_topo_data/elevations/greenland_dem_mosaic_100m_v3.0.tif"
gdal_data = gdal.Open(filename)
gdal_band = gdal_data.GetRasterBand(1)
nodataval = gdal_band.GetNoDataValue()

#2. Some adjustments
# convert to a numpy array
data_array = gdal_data.ReadAsArray().astype(np.float)

# replace missing values if necessary
if np.any(data_array == nodataval):
    data_array[data_array == nodataval] = np.nan

#3. Plot quickly the data
#Plot out data with Matplotlib's 'contour'
fig = pyplot.figure(figsize = (12, 8))
ax = fig.add_subplot(111)
pyplot.contour(data_array, cmap = "viridis", 
            levels = list(range(0, 3300, 1000)))
pyplot.title("Elevation Contours of GrIS")
cbar = pyplot.colorbar()
pyplot.gca().set_aspect('equal', adjustable='box')
pyplot.show()


##########################################################################
###                          Compute DEM contour	                       ###
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

#Plot excess melt
pyplot.figure(figsize=(48,40))   
#Change label font
pyplot.rcParams.update({'font.size': 20})
pyplot.grid()
pyplot.imshow(np.flipud(melt_year_plot),extent=grid.extent,cmap='plasma_r')
pyplot.clim(0,1000)
pyplot.colorbar(label='Excess melt [mm w.e./year]')
pyplot.title('Excess melt plot, year: '+wanted_year)
pyplot.contour(np.flipud(data_array),extent=grid.extent,cmap = "viridis", 
            levels = list(range(0, 3300, 1000)))
                    
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
