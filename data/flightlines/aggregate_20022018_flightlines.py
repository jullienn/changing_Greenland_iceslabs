# -*- coding: utf-8 -*-
"""
Created on Mon Dec 13 14:10:54 2021

@author: JullienN
"""

#Import packages
import pickle
import pandas as pd
import numpy as np
import pdb
from pyproj import Transformer
import geopandas as gpd
from shapely.geometry import Point, Polygon

#Load Rignot et al., 2016 Greenland drainage bassins
path_rignotetal2016_GrIS='C:/Users/jullienn/switchdrive/Private/research/backup_Aglaja/working_environment/greenland_topo_data/GRE_IceSheet_IMBIE2/GRE_IceSheet_IMBIE2/'
GrIS_rignotetal2016=gpd.read_file(path_rignotetal2016_GrIS+'GRE_IceSheet_IMBIE2_v1_EPSG3413.shp',rows=slice(1,2,1)) #the regions are the last rows of the shapefile
GrIS_mask=GrIS_rignotetal2016[GrIS_rignotetal2016.SUBREGION1=='ICE_SHEET']

################# Load 2002-2003 flightlines coordinates ################
path_20022003_flightlines='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/icelens_identification/'
#path_20022003_flightlines='/flash/jullienn/flightlines/data/'

#Open the file and read it
f_flightlines = open(path_20022003_flightlines+'metadata_coord_2002_2003', "rb")
all_2002_3_flightlines = pickle.load(f_flightlines)
f_flightlines.close()

#Define empty datasets
lat_all=[]
lon_all=[]

for year in list(all_2002_3_flightlines.keys()):
    for days in list(all_2002_3_flightlines[year].keys()):
        for indiv_file in list(all_2002_3_flightlines[year][days].keys()):
            if (indiv_file[0:7]=='quality'):
                continue
            else:
                print(indiv_file)
                lat_all=np.append(lat_all,all_2002_3_flightlines[year][days][indiv_file][0])
                lon_all=np.append(lon_all,all_2002_3_flightlines[year][days][indiv_file][1])

#Create a dataframe with 2002-2003 flightlines
flightlines_20022003=pd.DataFrame(lat_all,columns=['lat_3413'])
flightlines_20022003['lon_3413']=lon_all
flightlines_20022003['coords'] = list(zip(flightlines_20022003['lon_3413'],flightlines_20022003['lat_3413']))
flightlines_20022003['coords'] = flightlines_20022003['coords'].apply(Point)
points = gpd.GeoDataFrame(flightlines_20022003, geometry='coords', crs="EPSG:3413")
pointInPolys = gpd.tools.sjoin(points, GrIS_mask, op="within", how='left') #This is from https://www.matecdev.com/posts/point-in-polygon.html
flightlines_20022003_GrIS = points[pointInPolys.SUBREGION1=='ICE_SHEET']
#Add the year
flightlines_20022003_GrIS['str_year']=np.asarray(['2002-2003']*len(flightlines_20022003_GrIS))
flightlines_20022003_GrIS['year']=np.asarray([20022003]*len(flightlines_20022003_GrIS))
################# Load 2002-2003 flightlines coordinates ################

################# Load 2010-2018 flightlines coordinates ################
path_flightlines='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/flightlines/'

flightlines_2010=pd.read_csv(path_flightlines+'2010_Greenland_P3.csv',decimal='.',sep=',')
flightlines_2011=pd.read_csv(path_flightlines+'2011_Greenland_P3.csv',decimal='.',sep=',')
flightlines_2012=pd.read_csv(path_flightlines+'2012_Greenland_P3.csv',decimal='.',sep=',')
flightlines_2013=pd.read_csv(path_flightlines+'2013_Greenland_P3.csv',decimal='.',sep=',')
flightlines_2014=pd.read_csv(path_flightlines+'2014_Greenland_P3.csv',decimal='.',sep=',')
flightlines_2017=pd.read_csv(path_flightlines+'2017_Greenland_P3.csv',decimal='.',sep=',')
flightlines_2018=pd.read_csv(path_flightlines+'2018_Greenland_P3.csv',decimal='.',sep=',')
################# Load 2010-2018 flightlines coordinates ################

################ Transform 2010-2018 coordinates into EPSG:3413 ###############
#Transform the coordinates from WGS84 to EPSG:3413
transformer = Transformer.from_crs("EPSG:4326", "EPSG:3413", always_xy=True)

#2010
points=transformer.transform(np.asarray(flightlines_2010["LON"]),np.asarray(flightlines_2010["LAT"]))
flightlines_2010['lon_3413']=points[0]
flightlines_2010['lat_3413']=points[1]
flightlines_2010['coords'] = list(zip(flightlines_2010['lon_3413'],flightlines_2010['lat_3413']))
flightlines_2010['coords'] = flightlines_2010['coords'].apply(Point)
points = gpd.GeoDataFrame(flightlines_2010, geometry='coords', crs="EPSG:3413")
pointInPolys = gpd.tools.sjoin(points, GrIS_mask, op="within", how='left') #This is from https://www.matecdev.com/posts/point-in-polygon.html
flightlines_2010_GrIS = points[pointInPolys.SUBREGION1=='ICE_SHEET']
#Add the year
flightlines_2010_GrIS['str_year']=np.asarray(['2010']*len(flightlines_2010_GrIS))
flightlines_2010_GrIS['year']=np.asarray([2010]*len(flightlines_2010_GrIS))

#2011
points=transformer.transform(np.asarray(flightlines_2011["LON"]),np.asarray(flightlines_2011["LAT"]))
flightlines_2011['lon_3413']=points[0]
flightlines_2011['lat_3413']=points[1]
flightlines_2011['coords'] = list(zip(flightlines_2011['lon_3413'],flightlines_2011['lat_3413']))
flightlines_2011['coords'] = flightlines_2011['coords'].apply(Point)
points = gpd.GeoDataFrame(flightlines_2011, geometry='coords', crs="EPSG:3413")
pointInPolys = gpd.tools.sjoin(points, GrIS_mask, op="within", how='left') #This is from https://www.matecdev.com/posts/point-in-polygon.html
flightlines_2011_GrIS = points[pointInPolys.SUBREGION1=='ICE_SHEET']
#Add the year
flightlines_2011_GrIS['str_year']=np.asarray(['2011-2012']*len(flightlines_2011_GrIS))
flightlines_2011_GrIS['year']=np.asarray([2011]*len(flightlines_2010_GrIS))

#2012
points=transformer.transform(np.asarray(flightlines_2012["LON"]),np.asarray(flightlines_2012["LAT"]))
flightlines_2012['lon_3413']=points[0]
flightlines_2012['lat_3413']=points[1]
flightlines_2012['coords'] = list(zip(flightlines_2012['lon_3413'],flightlines_2012['lat_3413']))
flightlines_2012['coords'] = flightlines_2012['coords'].apply(Point)
points = gpd.GeoDataFrame(flightlines_2012, geometry='coords', crs="EPSG:3413")
pointInPolys = gpd.tools.sjoin(points, GrIS_mask, op="within", how='left') #This is from https://www.matecdev.com/posts/point-in-polygon.html
flightlines_2012_GrIS = points[pointInPolys.SUBREGION1=='ICE_SHEET']
#Add the year
flightlines_2012_GrIS['str_year']=np.asarray(['2011-2012']*len(flightlines_2012_GrIS))
flightlines_2012_GrIS['year']=np.asarray([2012]*len(flightlines_2012_GrIS))

#2013
points=transformer.transform(np.asarray(flightlines_2013["LON"]),np.asarray(flightlines_2013["LAT"]))
flightlines_2013['lon_3413']=points[0]
flightlines_2013['lat_3413']=points[1]
flightlines_2013['coords'] = list(zip(flightlines_2013['lon_3413'],flightlines_2013['lat_3413']))
flightlines_2013['coords'] = flightlines_2013['coords'].apply(Point)
points = gpd.GeoDataFrame(flightlines_2013, geometry='coords', crs="EPSG:3413")
pointInPolys = gpd.tools.sjoin(points, GrIS_mask, op="within", how='left') #This is from https://www.matecdev.com/posts/point-in-polygon.html
flightlines_2013_GrIS = points[pointInPolys.SUBREGION1=='ICE_SHEET']
#Add the year
flightlines_2013_GrIS['str_year']=np.asarray(['2013-2014']*len(flightlines_2013_GrIS))
flightlines_2013_GrIS['year']=np.asarray([2013]*len(flightlines_2012_GrIS))

#2014
points=transformer.transform(np.asarray(flightlines_2014["LON"]),np.asarray(flightlines_2014["LAT"]))
flightlines_2014['lon_3413']=points[0]
flightlines_2014['lat_3413']=points[1]
flightlines_2014['coords'] = list(zip(flightlines_2014['lon_3413'],flightlines_2014['lat_3413']))
flightlines_2014['coords'] = flightlines_2014['coords'].apply(Point)
points = gpd.GeoDataFrame(flightlines_2014, geometry='coords', crs="EPSG:3413")
pointInPolys = gpd.tools.sjoin(points, GrIS_mask, op="within", how='left') #This is from https://www.matecdev.com/posts/point-in-polygon.html
flightlines_2014_GrIS = points[pointInPolys.SUBREGION1=='ICE_SHEET']
#Add the year
flightlines_2014_GrIS['str_year']=np.asarray(['2013-2014']*len(flightlines_2014_GrIS))
flightlines_2014_GrIS['year']=np.asarray([2014]*len(flightlines_2012_GrIS))

#2017
points=transformer.transform(np.asarray(flightlines_2017["LON"]),np.asarray(flightlines_2017["LAT"]))
flightlines_2017['lon_3413']=points[0]
flightlines_2017['lat_3413']=points[1]
flightlines_2017['coords'] = list(zip(flightlines_2017['lon_3413'],flightlines_2017['lat_3413']))
flightlines_2017['coords'] = flightlines_2017['coords'].apply(Point)
points = gpd.GeoDataFrame(flightlines_2017, geometry='coords', crs="EPSG:3413")
pointInPolys = gpd.tools.sjoin(points, GrIS_mask, op="within", how='left') #This is from https://www.matecdev.com/posts/point-in-polygon.html
flightlines_2017_GrIS = points[pointInPolys.SUBREGION1=='ICE_SHEET']
#Add the year
flightlines_2017_GrIS['str_year']=np.asarray(['2017-2018']*len(flightlines_2017_GrIS))
flightlines_2017_GrIS['year']=np.asarray([2017]*len(flightlines_2017_GrIS))

#2018
points=transformer.transform(np.asarray(flightlines_2018["LON"]),np.asarray(flightlines_2018["LAT"]))
flightlines_2018['lon_3413']=points[0]
flightlines_2018['lat_3413']=points[1]
flightlines_2018['coords'] = list(zip(flightlines_2018['lon_3413'],flightlines_2018['lat_3413']))
flightlines_2018['coords'] = flightlines_2018['coords'].apply(Point)
points = gpd.GeoDataFrame(flightlines_2018, geometry='coords', crs="EPSG:3413")
pointInPolys = gpd.tools.sjoin(points, GrIS_mask, op="within", how='left') #This is from https://www.matecdev.com/posts/point-in-polygon.html
flightlines_2018_GrIS = points[pointInPolys.SUBREGION1=='ICE_SHEET']
#Add the year
flightlines_2018_GrIS['str_year']=np.asarray(['2017-2018']*len(flightlines_2018_GrIS))
flightlines_2018_GrIS['year']=np.asarray([2018]*len(flightlines_2018_GrIS))

################ Transform 2010-2018 coordinates into EPSG:3413 ###############

#Append all the flightlines together
flightlines_20102018_GrIS=flightlines_2010_GrIS
flightlines_20102018_GrIS=flightlines_20102018_GrIS.append(flightlines_2011_GrIS)
flightlines_20102018_GrIS=flightlines_20102018_GrIS.append(flightlines_2012_GrIS)
flightlines_20102018_GrIS=flightlines_20102018_GrIS.append(flightlines_2013_GrIS)
flightlines_20102018_GrIS=flightlines_20102018_GrIS.append(flightlines_2014_GrIS)
flightlines_20102018_GrIS=flightlines_20102018_GrIS.append(flightlines_2017_GrIS)
flightlines_20102018_GrIS=flightlines_20102018_GrIS.append(flightlines_2018_GrIS)

################# Aggregate 2002-2003 and 2010-2018 together ################
flightlines_20022018_GrIS=flightlines_20022003_GrIS.append(flightlines_20102018_GrIS)
################# Aggregate 2002-2003 and 2010-2018 together ################

#Save the generated file!
path_save='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/flightlines/'
flightlines_20022018_GrIS.to_csv(path_save+'flightlines_20022018_GrIS.csv')



