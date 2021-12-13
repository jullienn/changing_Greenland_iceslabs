# -*- coding: utf-8 -*-
"""
Created on Mon Dec 13 14:10:54 2021

@author: JullienN
"""


#Load 2002-2003 flightlines
################# Load 2002-2003 flightlines coordinates ################
path_20022003_flightlines='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/icelens_identification'

#Open the file and read it
f_flightlines = open(path_20022003_flightlines+'/metadata_coord_2002_2003', "rb")
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




################# Load 2002-2003 flightlines coordinates ################
check_belonging_GrIS = 'TRUE'

if (check_belonging_GrIS == 'TRUE'):
    #Load 2010-2018 flightlines
    path_flightlines='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/flightlines/'
    
    #flightlines_2010=pd.read_csv(path_flightlines+'2010_Greenland_P3.csv',decimal='.',sep=',')
    flightlines_2011=pd.read_csv(path_flightlines+'2011_Greenland_P3.csv',decimal='.',sep=',')
    flightlines_2012=pd.read_csv(path_flightlines+'2012_Greenland_P3.csv',decimal='.',sep=',')
    flightlines_2013=pd.read_csv(path_flightlines+'2013_Greenland_P3.csv',decimal='.',sep=',')
    flightlines_2014=pd.read_csv(path_flightlines+'2014_Greenland_P3.csv',decimal='.',sep=',')
    flightlines_2017=pd.read_csv(path_flightlines+'2017_Greenland_P3.csv',decimal='.',sep=',')
    flightlines_2018=pd.read_csv(path_flightlines+'2018_Greenland_P3.csv',decimal='.',sep=',')
    
    #Append all the flightlines together
    flightlines_20102018=flightlines_2011
    #all_flightlines=all_flightlines.append(flightlines_2011)
    flightlines_20102018=flightlines_20102018.append(flightlines_2012)
    flightlines_20102018=flightlines_20102018.append(flightlines_2013)
    flightlines_20102018=flightlines_20102018.append(flightlines_2014)
    flightlines_20102018=flightlines_20102018.append(flightlines_2017)
    flightlines_20102018=flightlines_20102018.append(flightlines_2018)
    
    #Transform the coordinated from WGS84 to EPSG:3413
    transformer = Transformer.from_crs("EPSG:4326", "EPSG:3413", always_xy=True)
    points=transformer.transform(np.asarray(flightlines_20102018["LON"]),np.asarray(flightlines_20102018["LAT"]))
    
    #Store lat/lon in 3413
    flightlines_20102018['lon_3413']=points[0]
    flightlines_20102018['lat_3413']=points[1]
    
    '''
    #Clip data to Rignot wt al., 2018 GrIS mask. If do not belong to, to not consider
    #Load Rignot et al., 2016 Greenland drainage bassins
    path_rignotetal2016_GrIS='C:/Users/jullienn/switchdrive/Private/research/backup_Aglaja/working_environment/greenland_topo_data/GRE_IceSheet_IMBIE2/GRE_IceSheet_IMBIE2/'
    GrIS_rignotetal2016=gpd.read_file(path_rignotetal2016_GrIS+'GRE_IceSheet_IMBIE2_v1.shp',rows=slice(1,2,1)) #the regions are the last rows of the shapefile
    GrIS_mask=GrIS_rignotetal2016[GrIS_rignotetal2016.SUBREGION1=='ICE_SHEET']

    #Create empty dataframe
    flightlines_20102018_GrIS=pd.DataFrame()
        
    count=0
    for i in range(0,len(flightlines_20102018)):
        print(count/len(flightlines_20102018)*100,' %')
    
        #select the point i
        single_point=Point(flightlines_20102018.iloc[i]['LON'],flightlines_20102018.iloc[i]['LAT'])
        
        #Do the identification between the point i and the regional shapefiles
        #From: https://automating-gis-processes.github.io/CSC18/lessons/L4/point-in-polygon.html
        check_GrIS=np.asarray(GrIS_mask.contains(single_point)).astype(int)
        
        #The point belongs to the GrIS, keep it
        if (np.sum(check_GrIS)>0):
            flightlines_20102018_GrIS[i]=flightlines_20102018[i]
        
        count=count+1
        
    #Save the generated file! Do that with 2002-2003 also?
    path_save='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/flightlines/'
    flightlines_20102018_GrIS.to_csv(path_save+'flightlines_20102018_GrIS.csv')
    
    
    pdb.set_trace()
    '''
else:
    #Load the data
    print('Data loaded')