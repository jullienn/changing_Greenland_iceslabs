# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 13:07:58 2022

@author: jullienn

Mainly inspired from Fig3and4.py
The choice of projection and the matching with satellite image and coordinates of 
radar track was made using the following references:
https://scitools.org.uk/cartopy/docs/latest/tutorials/understanding_transform.html
https://scitools.org.uk/cartopy/docs/v0.13/crs/projections.html
https://scitools.org.uk/cartopy/docs/v0.14/crs/index.html
https://epsg.io/32622
"""
from pyproj import CRS
import rioxarray as rxr
import cartopy.crs as ccrs
import numpy as np
import pdb
import matplotlib.pyplot as plt
from pyproj import Transformer
import pandas as pd

#Case study in Fig. 3
dictionnary_case_study={2010:['Data_20100515_01_007.mat','Data_20100515_01_008.mat','Data_20100515_01_009.mat'],
                        2011:['Data_20110408_01_087.mat','Data_20110408_01_088.mat','Data_20110408_01_089.mat',
                              'Data_20110408_01_090.mat','Data_20110408_01_091.mat','Data_20110408_01_092.mat',
                              'Data_20110408_01_093.mat','Data_20110408_01_094.mat','Data_20110408_01_095.mat',
                              'Data_20110408_01_096.mat','Data_20110408_01_097.mat','Data_20110408_01_098.mat',
                              'Data_20110408_01_099.mat','Data_20110408_01_100.mat','Data_20110408_01_101.mat',
                              'Data_20110408_01_102.mat','Data_20110408_01_103.mat'],
                        2012:['Data_20120423_01_137.mat','Data_20120423_01_138.mat'],
                        2013:['Data_20130409_01_010.mat','Data_20130409_01_011.mat','Data_20130409_01_012.mat'],
                        2014:['Data_20140416_05_035.mat','Data_20140416_05_036.mat','Data_20140416_05_037.mat'],
                        2017:['Data_20170502_01_171.mat','Data_20170502_01_172.mat','Data_20170502_01_173.mat'],
                        2018:['Data_20180421_01_004.mat','Data_20180421_01_005.mat','Data_20180421_01_006.mat','Data_20180421_01_007.mat']}

#Load all 2010-2018 data without spatial aggregation
path='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2002_2018/final_excel/high_estimate/clipped/'
df_2010_2018_csv = pd.read_csv(path+'Ice_Layer_Output_Thicknesses_2010_2018_jullienetal2023_high_estimate_cleaned.csv',delimiter=',',decimal='.')

#Define transformer for coordinates transform from "EPSG:4326" to "EPSG:32622" (=UTM 22N)
transformer = Transformer.from_crs("EPSG:4326", "EPSG:32622", always_xy=True)
#Transform the coordinated from "EPSG:4326" to "EPSG:32622"
points=transformer.transform(np.asarray(df_2010_2018_csv["lon"]),np.asarray(df_2010_2018_csv["lat"]))
#Store lat/lon in 32622
df_2010_2018_csv['lon_32622']=points[0]
df_2010_2018_csv['lat_32622']=points[1]

#Define empty dictionnary to store data related only to the case study
colnames=list(df_2010_2018_csv.keys())
colnames.append('UndersidesAcc')
df_casestudy=pd.DataFrame(columns=colnames)

#Loop over the years and flag where in the sector of ice accretion on the undersides of ice slabs
for year in dictionnary_case_study.keys():
    if (str(year) in list(['2010','2011','2014','2017'])):
        continue
    print(year)
    #Select data for the trace
    df_casestudy_temp=df_2010_2018_csv[df_2010_2018_csv['Track_name']==dictionnary_case_study[year][0][5:20]+'_'+dictionnary_case_study[year][-1][17:20]]
    #Do not keep where lon<-47.4233 (2012 is the constraining year) and lon> -46.2981 (2018 is the constraining year)
    df_casestudy_temp=df_casestudy_temp[np.logical_and(df_casestudy_temp['lon']>=-47.4233,df_casestudy_temp['lon']<=-46.2981)]
    #Where lon<-47.07 and lon> -47.11, ice accretion on the undersides of ice slabs is observed: flag UndersidesAcc with a 1
    df_casestudy_temp['UndersidesAcc']=np.asarray(np.logical_and(df_casestudy_temp['lon']>=-47.11,df_casestudy_temp['lon']<=-47.07)).astype(int)
    #Store the dataframe with UndersidesAcc in df_casestudy
    df_casestudy=df_casestudy.append(df_casestudy_temp)

# Define CRS of th satelite image as being the base for plotting
###################### From Tedstone et al., 2022 #####################
#from plot_map_decadal_change.py
# Define the CartoPy CRS object.
crs = ccrs.UTM(zone=22)
# This can be converted into a `proj4` string/dict compatible with GeoPandas
crs_proj4 = crs.proj4_init
###################### From Tedstone et al., 2022 #####################

#Load satelite NIR band
#This section of displaying sat data was coding using tips from
#https://www.earthdatascience.org/courses/use-data-open-source-python/intro-raster-data-python/raster-data-processing/reproject-raster/
#https://towardsdatascience.com/visualizing-satellite-data-using-matplotlib-and-cartopy-8274acb07b84

path_satellite='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/satellite_image/'
'''
#Load 2021 satelite data
sat_image = rxr.open_rasterio(path_satellite+'T22WFV_20210823T145759_B2348.tif',
                              masked=True).squeeze()#Spatial red also in UTM 22N
'''
path_2016='T22WFV_20160824T145915/S2A_MSIL1C_20160824T145912_N0204_R125_T22WFV_20160824T145915.SAFE/GRANULE/L1C_T22WFV_A006128_20160824T145915/IMG_DATA/'
#Load 2016 satelite data
sat_image = rxr.open_rasterio(path_satellite+path_2016+'T22WFV_20160824T145912_B08.jp2',
                              masked=True).squeeze() #No need to reproject satelite image
#Define extents based on the satelite image
extent_image = [np.asarray(sat_image.x[0]), np.asarray(sat_image.x[-1]), np.asarray(sat_image.y[-1]), np.asarray(sat_image.y[0])]#[west limit, east limit., south limit, north limit]

#Define fontsize
plt.rcParams.update({'font.size': 20})

#Prepare figure
plt.figure(figsize=(14,10))
#Define the projection of the axis
ax = plt.axes(projection=crs)
#Display satelite image
ax.imshow(sat_image, extent=extent_image, transform=crs,cmap='Blues_r', origin='upper',zorder=1) #NIR
'''
ax.imshow(sat_image[3,:,:], extent=extent_image, transform=crs,cmap='Blues_r', origin='lower',zorder=1) #NIR
'''
#Display radar tracks with the right coordinates "EPSG:32622"
ax.scatter(df_casestudy[df_casestudy['UndersidesAcc']==0]['lon_32622'],df_casestudy[df_casestudy['UndersidesAcc']==0]['lat_32622'],color='black',transform=crs,label='2012, 2013, 2018 radar track')
#Display where ice accretion on the undersides of ice slabs
ax.scatter(df_casestudy[df_casestudy['UndersidesAcc']==1]['lon_32622'],df_casestudy[df_casestudy['UndersidesAcc']==1]['lat_32622'],color='green',transform=crs,label='Ice accretion on the underside of ice slab')
#Display KAN_U in "EPSG:32622"
KAN_U_coord=transformer.transform([-47.0253],[67.0003])
ax.scatter(KAN_U_coord[0],KAN_U_coord[1],c='#b2182b',label='KAN_U',zorder=10,transform=crs)

#Set axis limits
ax.set_xlim(663836, 675955)
ax.set_ylim(7429308, 7438218)
#ax.set_extent(extent_image, crs=crs) 

###################### From Tedstone et al., 2022 #####################
#from plot_map_decadal_change.py
# x0, x1, y0, y1
gl=ax.gridlines(draw_labels=True, xlocs=[-47, -47.1, -47.2], ylocs=[66.95, 67], x_inline=False, y_inline=False,linewidth=0.5)
#Customize lat labels
gl.ylabels_right = False
gl.xlabels_bottom = False
ax.axis('off')
ax.legend(loc='lower right')
###################### From Tedstone et al., 2022 #####################

figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()

plt.show()

pdb.set_trace()
#Save the figure
plt.savefig('C:/Users/jullienn/switchdrive/Private/research/RT1/figures/S9/v1/figS9.png',dpi=500)



