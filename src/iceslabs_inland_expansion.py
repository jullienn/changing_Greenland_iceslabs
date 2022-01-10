# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 16:07:38 2021

@author: jullienn
"""
from scipy.spatial import Delaunay
from shapely.ops import cascaded_union, polygonize
from shapely.geometry import  MultiLineString

def alpha_shape(points, alpha):
    #This function is from https://gist.github.com/dwyerk/10561690
    """
    Compute the alpha shape (concave hull) of a set
    of points.
    @param points: Iterable container of points.
    @param alpha: alpha value to influence the
        gooeyness of the border. Smaller numbers
        don't fall inward as much as larger numbers.
        Too large, and you lose everything!
    """
    if len(points) < 4:
        # When you have a triangle, there is no sense
        # in computing an alpha shape.
        return geometry.MultiPoint(list(points)).convex_hull

    coords = np.array([p.coords[:][0] for p in points.coords])
    tri = Delaunay(coords)
    triangles = coords[tri.vertices]
    a = ((triangles[:,0,0] - triangles[:,1,0]) ** 2 + (triangles[:,0,1] - triangles[:,1,1]) ** 2) ** 0.5
    b = ((triangles[:,1,0] - triangles[:,2,0]) ** 2 + (triangles[:,1,1] - triangles[:,2,1]) ** 2) ** 0.5
    c = ((triangles[:,2,0] - triangles[:,0,0]) ** 2 + (triangles[:,2,1] - triangles[:,0,1]) ** 2) ** 0.5
    s = ( a + b + c ) / 2.0
    areas = (s*(s-a)*(s-b)*(s-c)) ** 0.5
    circums = a * b * c / (4.0 * areas)
    filtered = triangles[circums < (1.0 / alpha)]
    edge1 = filtered[:,(0,1)]
    edge2 = filtered[:,(1,2)]
    edge3 = filtered[:,(2,0)]
    edge_points = np.unique(np.concatenate((edge1,edge2,edge3)), axis = 0).tolist()
    m = MultiLineString(edge_points)
    triangles = list(polygonize(m))
    
    return cascaded_union(triangles), edge_points

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

def calcul_elevation(lon,lat,data_dem,yOrigin,pixelHeight,pixelWidth,index_lon_zero):
    
    if (np.isnan(lon) or np.isnan(lat)):
        #elev_all=np.append(elev_all,np.nan)
        elevation=np.nan
    else:
        #The origin is top left corner!!
        #y will always be negative
        row = int((yOrigin - lat ) / pixelHeight)
        if (lon<0):
            # if x negative
            col = index_lon_zero-int((-lon-0) / pixelWidth)
        elif (lon>0):
            # if x positive
            col = index_lon_zero+int((lon-0) / pixelWidth)
        #Read the elevation
        elevation=data_dem[row][col]
    
    return elevation

    
def concave_hull_computation(df_in_use,dictionnaries_convexhullmasks,ax1c,do_plot,input_file):
    from descartes.patch import PolygonPatch

    #Prepare for convex hull intersection
    df_in_use['coords'] = list(zip(df_in_use['lon_3413'],df_in_use['lat_3413']))
    df_in_use['coords'] = df_in_use['coords'].apply(Point)
    
    #Set summary_area
    summary_area={k: {} for k in list(['2011-2012','2017-2018'])}
    
    #Loop over time period
    for time_period in list(['2017-2018','2011-2012']):
        print(time_period)
        #Set color for plotting
        if (time_period == '2017-2018'):
            if (input_file=='low_end'):
                col_year='#fee0d2'
                set_alpha=0.2
            elif(input_file=='high_end'):
                col_year='#de2d26'
                set_alpha=0.5
            else:
                print('Input file not known, break')
            #Select data of the corresponding time period
            df_time_period=df_in_use[df_in_use['str_year']=='2011-2012'].append(df_in_use[df_in_use['str_year']=='2017-2018'])
            
        elif(time_period == '2011-2012'):
            if (input_file=='low_end'):
                col_year='#deebf7'
                set_alpha=0.2
            elif(input_file=='high_end'):
                col_year='#3182bd'
                set_alpha=0.5
            else:
                print('Input file not known, break')
            #Select data of the corresponding time period
            df_time_period=df_in_use[df_in_use['str_year']==time_period]
        else:
            print('Time period not known')
            break
        
        #Set summary_area
        summary_area[time_period]={k: {} for k in list(['NE','NO','NW','CW','SW'])}
                
        #Loop over each region and do the hull for each region of the IS
        for region in list(np.unique(df_time_period['key_shp'])):
            print('   ',region)
            #Select the corresponding region
            df_time_period_region=df_time_period[df_time_period['key_shp']==region]
            #Select point coordinates
            points = gpd.GeoDataFrame(df_time_period_region, geometry='coords', crs="EPSG:3413")

            if (region in list(['SE','Out'])):
                #do not compute, continue
                continue
            #reset area region to 0
            area_region=0
            
            # Perform spatial join to match points and polygons
            for convex_hull_mask in dictionnaries_convexhullmasks[region].keys():
                print('      ',convex_hull_mask)
                pointInPolys = gpd.tools.sjoin(points, dictionnaries_convexhullmasks[region][convex_hull_mask], op="within", how='left') #This is from https://www.matecdev.com/posts/point-in-polygon.html
                #Keep only matched point
                if (region in list(['CW','SW'])):
                    pnt_matched = points[pointInPolys.SUBREGION1==region]
                else:
                    pnt_matched = points[pointInPolys.id==1]
                
                if (len(pnt_matched)>1):
                    #pdb.set_trace()
                    #this function is from https://gist.github.com/dwyerk/10561690
                    concave_hull, edge_points= alpha_shape(pnt_matched, 0.00002) #0.00005 is a bit too aggresive, 0.00001 is a bit generous     
                    
                    alpha_play=0.00002
                    while (len(concave_hull.bounds)==0):
                        print('Empty alpha! Iterate until fit is made')
                        #Could not match a poylgon with specified alpha. Decrease alpha until poylgon can be matched
                        #update alpha_play
                        alpha_play=alpha_play/2
                        print(alpha_play)
                        concave_hull, edge_points= alpha_shape(pnt_matched, alpha_play) #0.00005 is a bit too aggresive, 0.00001 is a bit generous
                        
                    if (do_plot=='TRUE'):
                        patch1 = PolygonPatch(concave_hull, zorder=2, alpha=set_alpha,color=col_year)
                        ax1c.add_patch(patch1)
                        #ax1c.scatter(pnt_matched.lon_3413,pnt_matched.lat_3413,zorder=3)
                        plt.show()
                    
                    #Update area_region
                    area_region=area_region+concave_hull.area #IS THAT CORRECT???
            
            #Store total area per region and per time period
            summary_area[time_period][region]=area_region  
            
    return summary_area


def plot_fig1(df_all,flightlines_20022018,df_2010_2018_low,df_2010_2018_high,df_firn_aquifer_all,df_thickness_likelihood_20102018):   
    plot_fig_S1='FALSE'
    plot_panela='FALSE'
    plot_panelb='FALSE'
    plot_panelc='TRUE'
    
    if (plot_fig_S1 == 'TRUE'):
        # -------------------------------- FIG S1 --------------------------------
        #prepare the figure
        fig, (ax1) = plt.subplots(1, 1)#, gridspec_kw={'width_ratios': [1, 3]})
        fig.set_size_inches(40,20)
        
        #Display GrIS drainage bassins
        NO_rignotetal.plot(ax=ax1,color='white', edgecolor='black')
        NE_rignotetal.plot(ax=ax1,color='white', edgecolor='black') 
        SE_rignotetal.plot(ax=ax1,color='white', edgecolor='black') 
        SW_rignotetal.plot(ax=ax1,color='white', edgecolor='black') 
        CW_rignotetal.plot(ax=ax1,color='white', edgecolor='black') 
        NW_rignotetal.plot(ax=ax1,color='white', edgecolor='black') 
        
        #Display 2002-2018 flightlines
        plt.scatter(flightlines_20022018['lon_3413'],flightlines_20022018['lat_3413'],s=0.001,color='#d9d9d9',label='Flightlines')#,label='2002-2003')
        
        #Display firn aquifers
        plt.scatter(df_firn_aquifer_all['lon_3413'],df_firn_aquifer_all['lat_3413'],s=1,color='#238b45',label='Firn aquifers')
        
        #Display 2010-2018 iceslabs thickness
        plt.scatter(df_all[df_all.Track_name.str[:4]=='2010']['lon_3413'],df_all[df_all.Track_name.str[:4]=='2010']['lat_3413'],c=df_all[df_all.Track_name.str[:4]=='2010']['20m_ice_content_m'],s=1,cmap=plt.get_cmap('Blues'),label='2010-2018 ice slabs')
        plt.scatter(df_all[df_all.Track_name.str[:4]=='2011']['lon_3413'],df_all[df_all.Track_name.str[:4]=='2011']['lat_3413'],c=df_all[df_all.Track_name.str[:4]=='2011']['20m_ice_content_m'],s=1,cmap=plt.get_cmap('Blues'))
        plt.scatter(df_all[df_all.Track_name.str[:4]=='2012']['lon_3413'],df_all[df_all.Track_name.str[:4]=='2012']['lat_3413'],c=df_all[df_all.Track_name.str[:4]=='2012']['20m_ice_content_m'],s=1,cmap=plt.get_cmap('Blues'))
        plt.scatter(df_all[df_all.Track_name.str[:4]=='2013']['lon_3413'],df_all[df_all.Track_name.str[:4]=='2013']['lat_3413'],c=df_all[df_all.Track_name.str[:4]=='2013']['20m_ice_content_m'],s=1,cmap=plt.get_cmap('Blues'))
        plt.scatter(df_all[df_all.Track_name.str[:4]=='2014']['lon_3413'],df_all[df_all.Track_name.str[:4]=='2014']['lat_3413'],c=df_all[df_all.Track_name.str[:4]=='2014']['20m_ice_content_m'],s=1,cmap=plt.get_cmap('Blues'))
        plt.scatter(df_all[df_all.Track_name.str[:4]=='2017']['lon_3413'],df_all[df_all.Track_name.str[:4]=='2017']['lat_3413'],c=df_all[df_all.Track_name.str[:4]=='2017']['20m_ice_content_m'],s=1,cmap=plt.get_cmap('Blues'))
        lik_blues=plt.scatter(df_all[df_all.Track_name.str[:4]=='2018']['lon_3413'],df_all[df_all.Track_name.str[:4]=='2018']['lat_3413'],c=df_all[df_all.Track_name.str[:4]=='2018']['20m_ice_content_m'],s=1,cmap=plt.get_cmap('Blues'))
        
        cbar_blues = plt.colorbar(lik_blues)
        cbar_blues.set_label('Columnal ice content [m]')
        
        #Display 2002-2003 iceslabs
        plt.scatter(df_all[df_all.str_year=='2002-2003']['lon_3413'],df_all[df_all.str_year=='2002-2003']['lat_3413'],s=0.1,color='#8c6bb1',label='2002-2003 ice slabs')
        
        #Display region name on panel a 
        ax1.text(NO_rignotetal.centroid.x,NO_rignotetal.centroid.y+20000,np.asarray(NO_rignotetal.SUBREGION1)[0])
        ax1.text(NE_rignotetal.centroid.x,NE_rignotetal.centroid.y+20000,np.asarray(NE_rignotetal.SUBREGION1)[0])
        ax1.text(SE_rignotetal.centroid.x,SE_rignotetal.centroid.y+20000,np.asarray(SE_rignotetal.SUBREGION1)[0])
        ax1.text(SW_rignotetal.centroid.x,SW_rignotetal.centroid.y+20000,np.asarray(SW_rignotetal.SUBREGION1)[0])
        ax1.text(CW_rignotetal.centroid.x,CW_rignotetal.centroid.y+20000,np.asarray(CW_rignotetal.SUBREGION1)[0])
        ax1.text(NW_rignotetal.centroid.x,NW_rignotetal.centroid.y+20000,np.asarray(NW_rignotetal.SUBREGION1)[0])
        
        ax1.set_xlabel('Easting [m]')
        ax1.set_ylabel('Northing [m]')
    
        # Plot legend. This is from https://stackoverflow.com/questions/24706125/setting-a-fixed-size-for-points-in-legend
        lgnd = plt.legend(loc="lower right", scatterpoints=1, fontsize=10)
        lgnd.legendHandles[0]._sizes = [30]
        lgnd.legendHandles[1]._sizes = [30]
        lgnd.legendHandles[2]._sizes = [30]
        lgnd.legendHandles[3]._sizes = [30]
        
        ax1.set_xlim(-645000,855000)
        ax1.set_ylim(-3330000,-785000)
        
        figManager = plt.get_current_fig_manager()
        figManager.window.showMaximized()
        
        #Save the figure
        plt.savefig('C:/Users/jullienn/switchdrive/Private/research/RT1/figures/supp/v1/figS1.png',dpi=1000)
        # -------------------------------- FIG S1 --------------------------------
    
    #pdb.set_trace()
    
    if (plot_panela=='TRUE'):
        # -------------------------------- PANEL A --------------------------------
        panel_a_save='TRUE'

        #prepare the figure
        fig, (ax1) = plt.subplots(1, 1)#, gridspec_kw={'width_ratios': [1, 3]})
        fig.set_size_inches(40,20)
        
        #Display GrIS drainage bassins
        NO_rignotetal.plot(ax=ax1,color='white', edgecolor='black')
        NE_rignotetal.plot(ax=ax1,color='white', edgecolor='black') 
        SE_rignotetal.plot(ax=ax1,color='white', edgecolor='black') 
        SW_rignotetal.plot(ax=ax1,color='white', edgecolor='black') 
        CW_rignotetal.plot(ax=ax1,color='white', edgecolor='black') 
        NW_rignotetal.plot(ax=ax1,color='white', edgecolor='black') 
        
        #Display 2002-2018 flightlines
        plt.scatter(flightlines_20022018['lon_3413'],flightlines_20022018['lat_3413'],s=0.001,color='#d9d9d9',label='Flightlines')#,label='2002-2003')
        
        #Display firn aquifers
        plt.scatter(df_firn_aquifer_all['lon_3413'],df_firn_aquifer_all['lat_3413'],s=1,color='#238b45',label='Firn aquifers')
    
        #Display 2010-2018 iceslabs
        plt.scatter(df_all[df_all.Track_name.str[:4]=='2010']['lon_3413'],df_all[df_all.Track_name.str[:4]=='2010']['lat_3413'],s=1,color='#3690c0',label='2010-2014 ice slabs')
        plt.scatter(df_all[df_all.Track_name.str[:4]=='2011']['lon_3413'],df_all[df_all.Track_name.str[:4]=='2011']['lat_3413'],s=1,color='#3690c0')
        plt.scatter(df_all[df_all.Track_name.str[:4]=='2012']['lon_3413'],df_all[df_all.Track_name.str[:4]=='2012']['lat_3413'],s=1,color='#3690c0')
        plt.scatter(df_all[df_all.Track_name.str[:4]=='2013']['lon_3413'],df_all[df_all.Track_name.str[:4]=='2013']['lat_3413'],s=1,color='#3690c0')
        plt.scatter(df_all[df_all.Track_name.str[:4]=='2014']['lon_3413'],df_all[df_all.Track_name.str[:4]=='2014']['lat_3413'],s=1,color='#3690c0')
        plt.scatter(df_all[df_all.Track_name.str[:4]=='2017']['lon_3413'],df_all[df_all.Track_name.str[:4]=='2017']['lat_3413'],s=1,color='#a6bddb',label='2017-2018 ice slabs')
        plt.scatter(df_all[df_all.Track_name.str[:4]=='2018']['lon_3413'],df_all[df_all.Track_name.str[:4]=='2018']['lat_3413'],s=1,color='#a6bddb')
        
        #Display 2002-2003 iceslabs
        plt.scatter(df_all[df_all.str_year=='2002-2003']['lon_3413'],df_all[df_all.str_year=='2002-2003']['lat_3413'],s=1,color='#0570b0',label='2002-2003 ice slabs')
        
        #Display region name on panel a 
        ax1.text(NO_rignotetal.centroid.x,NO_rignotetal.centroid.y+20000,np.asarray(NO_rignotetal.SUBREGION1)[0])
        ax1.text(NE_rignotetal.centroid.x,NE_rignotetal.centroid.y+20000,np.asarray(NE_rignotetal.SUBREGION1)[0])
        ax1.text(SE_rignotetal.centroid.x,SE_rignotetal.centroid.y+20000,np.asarray(SE_rignotetal.SUBREGION1)[0])
        ax1.text(SW_rignotetal.centroid.x,SW_rignotetal.centroid.y+20000,np.asarray(SW_rignotetal.SUBREGION1)[0])
        ax1.text(CW_rignotetal.centroid.x,CW_rignotetal.centroid.y+20000,np.asarray(CW_rignotetal.SUBREGION1)[0])
        ax1.text(NW_rignotetal.centroid.x,NW_rignotetal.centroid.y+20000,np.asarray(NW_rignotetal.SUBREGION1)[0])
                
        ax1.set_xlabel('Easting [m]')
        ax1.set_ylabel('Northing [m]')
    
        # Plot legend. This is from https://stackoverflow.com/questions/24706125/setting-a-fixed-size-for-points-in-legend
        lgnd = plt.legend(loc="lower right", scatterpoints=1, fontsize=10)
        lgnd.legendHandles[0]._sizes = [30]
        lgnd.legendHandles[1]._sizes = [30]
        lgnd.legendHandles[2]._sizes = [30]
        lgnd.legendHandles[3]._sizes = [30]
        lgnd.legendHandles[4]._sizes = [30]
        
        ax1.set_xlim(-645000,855000)
        ax1.set_ylim(-3330000,-785000)
        
        figManager = plt.get_current_fig_manager()
        figManager.window.showMaximized()
        
        if (panel_a_save == 'TRUE'):
            #Save the figure
            plt.savefig('C:/Users/jullienn/switchdrive/Private/research/RT1/figures/fig1/v2/fig1_panel_a.png',dpi=1000)
     
        plt.close(fig)
        # -------------------------------- PANEL A --------------------------------

    #pdb.set_trace()
    
    if (plot_panelb=='TRUE'):
        panel_b_save='TRUE'
        
        # -------------------------------- PANEL B --------------------------------    
        #Define panel names
        labels = ['NE', 'NO', 'NW', 'CW', 'SW']
    
        #Stack max_elev_mean data for barplot
        dplot_20022003=[dict_summary['NE']['2002-2003']['max_elev_median'],dict_summary['NO']['2002-2003']['max_elev_median'],
                        dict_summary['NW']['2002-2003']['max_elev_median'],dict_summary['CW']['2002-2003']['max_elev_median'],
                        dict_summary['SW']['2002-2003']['max_elev_median']]
        
        dplot_2010=[dict_summary['NE']['2010']['max_elev_median'],dict_summary['NO']['2010']['max_elev_median'],
                    dict_summary['NW']['2010']['max_elev_median'],dict_summary['CW']['2010']['max_elev_median'],
                    dict_summary['SW']['2010']['max_elev_median']]
        
        dplot_20112012=[dict_summary['NE']['2011-2012']['max_elev_median'],dict_summary['NO']['2011-2012']['max_elev_median'],
                        dict_summary['NW']['2011-2012']['max_elev_median'],dict_summary['CW']['2011-2012']['max_elev_median'],
                        dict_summary['SW']['2011-2012']['max_elev_median']]
        
        dplot_20132014=[dict_summary['NE']['2013-2014']['max_elev_median'],dict_summary['NO']['2013-2014']['max_elev_median'],
                        dict_summary['NW']['2013-2014']['max_elev_median'],dict_summary['CW']['2013-2014']['max_elev_median'],
                        dict_summary['SW']['2013-2014']['max_elev_median']]
        
        dplot_20172018=[dict_summary['NE']['2017-2018']['max_elev_median'],dict_summary['NO']['2017-2018']['max_elev_median'],
                        dict_summary['NW']['2017-2018']['max_elev_median'],dict_summary['CW']['2017-2018']['max_elev_median'],
                        dict_summary['SW']['2017-2018']['max_elev_median']]
        
        #Stack max_elev_mean data for barplot
        dplotstd_20022003=[dict_summary['NE']['2002-2003']['max_elev_std'],dict_summary['NO']['2002-2003']['max_elev_std'],
                        dict_summary['NW']['2002-2003']['max_elev_std'],dict_summary['CW']['2002-2003']['max_elev_std'],
                        dict_summary['SW']['2002-2003']['max_elev_std']]
        
        dplotstd_2010=[dict_summary['NE']['2010']['max_elev_std'],dict_summary['NO']['2010']['max_elev_std'],
                    dict_summary['NW']['2010']['max_elev_std'],dict_summary['CW']['2010']['max_elev_std'],
                    dict_summary['SW']['2010']['max_elev_std']]
        
        dplotstd_20112012=[dict_summary['NE']['2011-2012']['max_elev_std'],dict_summary['NO']['2011-2012']['max_elev_std'],
                        dict_summary['NW']['2011-2012']['max_elev_std'],dict_summary['CW']['2011-2012']['max_elev_std'],
                        dict_summary['SW']['2011-2012']['max_elev_std']]
        
        dplotstd_20132014=[dict_summary['NE']['2013-2014']['max_elev_std'],dict_summary['NO']['2013-2014']['max_elev_std'],
                        dict_summary['NW']['2013-2014']['max_elev_std'],dict_summary['CW']['2013-2014']['max_elev_std'],
                        dict_summary['SW']['2013-2014']['max_elev_std']]
        
        dplotstd_20172018=[dict_summary['NE']['2017-2018']['max_elev_std'],dict_summary['NO']['2017-2018']['max_elev_std'],
                        dict_summary['NW']['2017-2018']['max_elev_std'],dict_summary['CW']['2017-2018']['max_elev_std'],
                        dict_summary['SW']['2017-2018']['max_elev_std']]
        
        
        #Stack data for maximum elevation difference calculation
        max_elev_diff_NE=[dict_summary['NE']['2002-2003']['max_elev_median'],dict_summary['NE']['2010']['max_elev_median'],
                          dict_summary['NE']['2011-2012']['max_elev_median'],dict_summary['NE']['2013-2014']['max_elev_median'],
                          dict_summary['NE']['2017-2018']['max_elev_median']]
        
        max_elev_diff_NO=[dict_summary['NO']['2002-2003']['max_elev_median'],dict_summary['NO']['2010']['max_elev_median'],
                          dict_summary['NO']['2011-2012']['max_elev_median'],dict_summary['NO']['2013-2014']['max_elev_median'],
                          dict_summary['NO']['2017-2018']['max_elev_median']]
        
        max_elev_diff_NW=[dict_summary['NW']['2002-2003']['max_elev_median'],dict_summary['NW']['2010']['max_elev_median'],
                          dict_summary['NW']['2011-2012']['max_elev_median'],dict_summary['NW']['2013-2014']['max_elev_median'],
                          dict_summary['NW']['2017-2018']['max_elev_median']]
        
        max_elev_diff_CW=[dict_summary['CW']['2002-2003']['max_elev_median'],dict_summary['CW']['2010']['max_elev_median'],
                          dict_summary['CW']['2011-2012']['max_elev_median'],dict_summary['CW']['2013-2014']['max_elev_median'],
                          dict_summary['CW']['2017-2018']['max_elev_median']]
        
        max_elev_diff_SW=[dict_summary['SW']['2002-2003']['max_elev_median'],dict_summary['SW']['2010']['max_elev_median'],
                          dict_summary['SW']['2011-2012']['max_elev_median'],dict_summary['SW']['2013-2014']['max_elev_median'],
                          dict_summary['SW']['2017-2018']['max_elev_median']]
    
        #Barplot inspired from https://stackoverflow.com/questions/10369681/how-to-plot-bar-graphs-with-same-x-coordinates-side-by-side-dodged
        #Arguments for barplot
        width = 0.1# the width of the bars: can also be len(x) sequence
        N=5 #Number of regions
        ind= np.arange(N) #Position of regions
            
        fig, ax = plt.subplots()
        ax.bar(ind, dplot_20022003, width, label='2002-2003',color='#c6dbef', yerr= dplotstd_20022003) #yerr=men_std
        ax.bar(ind+1*width, dplot_2010, width, label='2010',color='#9ecae1', yerr= dplotstd_2010)
        ax.bar(ind+2*width, dplot_20112012, width, label='2011-2012',color='#6baed6', yerr= dplotstd_20112012)
        ax.bar(ind+3*width, dplot_20132014, width, label='2013-2014',color='#3182bd', yerr= dplotstd_20132014)
        ax.bar(ind+4*width, dplot_20172018, width, label='2017-2018',color='#08519c', yerr= dplotstd_20172018)
        ax.set_xticks(ind + 2*width)
        ax.set_xticklabels(labels)
        ax.set_ylim(1000,2050)
        
        ax.text(ind[0],np.nanmax(max_elev_diff_NE)+50,str(int(np.round(np.nanmax(max_elev_diff_NE)-np.nanmin(max_elev_diff_NE))))+' m')
        ax.text(ind[1],np.nanmax(max_elev_diff_NO)+50,str(int(np.round(np.nanmax(max_elev_diff_NO)-np.nanmin(max_elev_diff_NO))))+' m')
        ax.text(ind[2],np.nanmax(max_elev_diff_NW)+50,str(int(np.round(np.nanmax(max_elev_diff_NW)-np.nanmin(max_elev_diff_NW))))+' m')
        ax.text(ind[3],np.nanmax(max_elev_diff_CW)+50,str(int(np.round(np.nanmax(max_elev_diff_CW)-np.nanmin(max_elev_diff_CW))))+' m')
        ax.text(ind[4],np.nanmax(max_elev_diff_SW)+50,str(int(np.round(np.nanmax(max_elev_diff_SW)-np.nanmin(max_elev_diff_SW))))+' m')
        
        ax.set_ylabel('Elevation [m]')
        ax.set_title('Median of ice slabs maximum elevation per slice')
        
        #Custom legend myself
        from matplotlib.patches import Patch
        from matplotlib.lines import Line2D
        
        legend_elements = [Patch(facecolor='#c6dbef', alpha=0.5,label='2002-2003'),
                           Patch(facecolor='#9ecae1', alpha=0.5,label='2010'),
                           Patch(facecolor='#6baed6', alpha=0.5,label='2011-2012'),
                           Patch(facecolor='#3182bd', alpha=0.5,label='2013-2014'),
                           Patch(facecolor='#08519c', alpha=0.5,label='2017-2018'),
                           Line2D([0], [0], color='k', lw=2, label='Standard deviation around the mean')]
        
        ax.legend(handles=legend_elements,loc='upper center')
        plt.legend()
        plt.show()
        
        if (panel_b_save == 'TRUE'):
            #Save the figure
            plt.savefig('C:/Users/jullienn/switchdrive/Private/research/RT1/figures/fig1/v2/fig1_panel_b.png',dpi=1000)
        
        plt.close(fig)
        # -------------------------------- PANEL B --------------------------------    
    
    #pdb.set_trace()

    if (plot_panelc=='TRUE'):
        
        hull_computation='TRUE'
        likelihood_display='TRUE'
        
        panel_c_save='TRUE'

        # -------------------------------- PANELS C -------------------------------        
        if (hull_computation=='TRUE'):
            #Panel C
            #Load convex hull mask over which convex hull must be computed
            path_convexhull_masks='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2010_2018/shapefiles/'
            
            dictionnaries_convexhullmasks = {k: {} for k in list(['NE','NO','NW','CW','SW'])}
            dictionnaries_convexhullmasks['NE']={k: {} for k in list(['NE_CH_1','NE_CH_2','NE_CH_3','NE_CH_4'])}
            dictionnaries_convexhullmasks['NO']={k: {} for k in list(['NO_CH_1','NO_CH_2','NO_CH_3','NO_CH_4','NO_CH_5','NO_CH_6','NO_CH_7'])}
            dictionnaries_convexhullmasks['NW']={k: {} for k in list(['NW_CH_1','NW_CH_2','NW_CH_3','NW_CH_4','NW_CH_5','NW_CH_6'])}
            dictionnaries_convexhullmasks['CW']={k: {} for k in list(['CW_CH_1'])}
            dictionnaries_convexhullmasks['SW']={k: {} for k in list(['SW_CH_1'])}
            
            for indiv_region in dictionnaries_convexhullmasks.keys():
                for indiv_file in dictionnaries_convexhullmasks[indiv_region].keys():
                    print('convex_hull_'+indiv_file[0:2]+indiv_file[5:8]+'.shp')
                    if (indiv_region == 'CW'):
                        dictionnaries_convexhullmasks[indiv_region][indiv_file]=CW_rignotetal
                    elif (indiv_region == 'SW'):
                        dictionnaries_convexhullmasks[indiv_region][indiv_file]=SW_rignotetal
                    else:
                        dictionnaries_convexhullmasks[indiv_region][indiv_file]=gpd.read_file(path_convexhull_masks+'convex_hull_'+indiv_file[0:2]+indiv_file[5:8]+'.shp')
            
            #Generate shapefile from iceslabs data. This si from https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.ConvexHull.html
            from scipy.spatial import ConvexHull
            from shapely import geometry
            from shapely.ops import unary_union
        
            #prepare the figure
            figc, (ax1c) = plt.subplots(1, 1)#, gridspec_kw={'width_ratios': [1, 3]})
            figc.suptitle('GrIS ice slabs extent - high estimate')
            
            #Display GrIS drainage bassins
            NO_rignotetal.plot(ax=ax1c,color='white', edgecolor='black')
            NE_rignotetal.plot(ax=ax1c,color='white', edgecolor='black') 
            SE_rignotetal.plot(ax=ax1c,color='white', edgecolor='black') 
            SW_rignotetal.plot(ax=ax1c,color='white', edgecolor='black') 
            CW_rignotetal.plot(ax=ax1c,color='white', edgecolor='black') 
            NW_rignotetal.plot(ax=ax1c,color='white', edgecolor='black')
            
            '''
            #Load ELA shapefile
            path_ELA='C:/Users/jullienn/switchdrive/Private/research/backup_Aglaja/working_environment/greenland_topo_data/ELA/'
            ELA=gpd.read_file(path_ELA+'RACMO.3p2_ERA5_3h_FGRN055.1km.1980-2001.3413.accum.sieve10.vect_EPSG3413_right.shp')
            above_ELA=ELA[ELA['fid']==2279]
            below_ELA=ELA[ELA['DN']==0]
            
            #Plot below ELA
            below_ELA.plot(ax=ax1c,color='red', edgecolor='black')
            
            plt.scatter(df_2010_2018_high[df_2010_2018_high.Track_name.str[:4]=='2010']['lon_3413'],df_2010_2018_high[df_2010_2018_high.Track_name.str[:4]=='2010']['lat_3413'],s=1,color='#3690c0',label='2010-2014 ice slabs')
            plt.scatter(df_2010_2018_high[df_2010_2018_high.Track_name.str[:4]=='2011']['lon_3413'],df_2010_2018_high[df_2010_2018_high.Track_name.str[:4]=='2011']['lat_3413'],s=1,color='#3690c0')
            plt.scatter(df_2010_2018_high[df_2010_2018_high.Track_name.str[:4]=='2012']['lon_3413'],df_2010_2018_high[df_2010_2018_high.Track_name.str[:4]=='2012']['lat_3413'],s=1,color='#3690c0')
            plt.scatter(df_2010_2018_high[df_2010_2018_high.Track_name.str[:4]=='2013']['lon_3413'],df_2010_2018_high[df_2010_2018_high.Track_name.str[:4]=='2013']['lat_3413'],s=1,color='#3690c0')
            plt.scatter(df_2010_2018_high[df_2010_2018_high.Track_name.str[:4]=='2014']['lon_3413'],df_2010_2018_high[df_2010_2018_high.Track_name.str[:4]=='2014']['lat_3413'],s=1,color='#3690c0')
            plt.scatter(df_2010_2018_high[df_2010_2018_high.Track_name.str[:4]=='2017']['lon_3413'],df_2010_2018_high[df_2010_2018_high.Track_name.str[:4]=='2017']['lat_3413'],s=1,color='#a6bddb',label='2017-2018 ice slabs')
            plt.scatter(df_2010_2018_high[df_2010_2018_high.Track_name.str[:4]=='2018']['lon_3413'],df_2010_2018_high[df_2010_2018_high.Track_name.str[:4]=='2018']['lat_3413'],s=1,color='#a6bddb')
            
            plt.legend()
            plt.show()
            
            pdb.set_trace()
            
            ####################### From concave hull computation #######################
            #Prepare for convex hull intersection
            df_2010_2018_high['coords'] = list(zip(df_2010_2018_high['lon_3413'],df_2010_2018_high['lat_3413']))
            df_2010_2018_high['coords'] = df_2010_2018_high['coords'].apply(Point)
            
            points = gpd.GeoDataFrame(df_2010_2018_high, geometry='coords', crs="EPSG:3413")


            pointInPolys = gpd.tools.sjoin(points, below_ELA, op="within", how='left') #This is from https://www.matecdev.com/posts/point-in-polygon.html
            pnt_matched = points[pointInPolys.SUBREGION1==region]

            '''
            '''
            if (region in list(['SE','Out'])):
                #do not compute, continue
                continue
            #reset area region to 0
            area_region=0
            
            # Perform spatial join to match points and polygons
            for convex_hull_mask in dictionnaries_convexhullmasks[region].keys():
                print('      ',convex_hull_mask)
                pointInPolys = gpd.tools.sjoin(points, dictionnaries_convexhullmasks[region][convex_hull_mask], op="within", how='left') #This is from https://www.matecdev.com/posts/point-in-polygon.html
                #Keep only matched point
                if (region in list(['CW','SW'])):
                    pnt_matched = points[pointInPolys.SUBREGION1==region]
                else:
                    pnt_matched = points[pointInPolys.id==1]
            '''
            ####################### From concave hull computation #######################
            
            #Calculate concave hull and extract low and high end areas
            do_plot='FALSE'
            high_end_summary=concave_hull_computation(df_2010_2018_high,dictionnaries_convexhullmasks,ax1c,do_plot,'high_end')
            do_plot='FALSE'
            low_end_summary=concave_hull_computation(df_2010_2018_low,dictionnaries_convexhullmasks,ax1c,do_plot,'low_end')
                        
            #Display area change on the figure
            for region in list(['NE','NO','NW','CW','SW']):
                if (region =='NE'):
                    polygon_for_text=NE_rignotetal
                elif(region =='NO'):
                    polygon_for_text=NO_rignotetal
                elif(region =='NW'):
                    polygon_for_text=NW_rignotetal
                elif(region =='CW'):
                    polygon_for_text=CW_rignotetal
                elif(region =='SW'):
                    polygon_for_text=SW_rignotetal
                else:
                    print('Region not known')
                
                low_end_change=(int((low_end_summary['2017-2018'][region]-low_end_summary['2011-2012'][region])/low_end_summary['2011-2012'][region]*100))
                high_end_change=(int((high_end_summary['2017-2018'][region]-high_end_summary['2011-2012'][region])/high_end_summary['2011-2012'][region]*100))
                
                #Display region name on panel c
                ax1c.text(polygon_for_text.centroid.x+25000,polygon_for_text.centroid.y+20000,region)
                
                #Compute and display relative change
                ax1c.text(polygon_for_text.centroid.x,polygon_for_text.centroid.y,'[+'+str(low_end_change)+' : +'+str(high_end_change)+'] %')
                
            ax1c.set_xlabel('Easting [m]')
            ax1c.set_ylabel('Northing [m]')
            '''
            #Custom legend myself
            legend_elements = [Patch(facecolor='#3182bd', alpha=0.5,label='2011-2012'),
                               Patch(facecolor='#de2d26', alpha=0.5,label='2017-2018')]
            ax1c.legend(handles=legend_elements,loc='best')
            plt.legend()
            '''
            #Plot data
            #ax1c.scatter(df_all.lon_3413,df_all.lat_3413,s=0.1,zorder=3)
            
            '''
            if (plot_save == 'TRUE'):
                #Save the figure
                pdb.set_trace()
                #plt.savefig('C:/Users/jullienn/switchdrive/Private/research/RT1/figures/fig1_panels_c.png',dpi=2000)
                #plt.close(fig)
            '''
        
        if (likelihood_display=='TRUE'):
            #prepare the figure
            figc, (ax1c) = plt.subplots(1, 1)
            figc.set_size_inches(20,10)
            
            #Display GrIS drainage bassins
            NO_rignotetal.plot(ax=ax1c,color='white', edgecolor='black')
            NE_rignotetal.plot(ax=ax1c,color='white', edgecolor='black') 
            SE_rignotetal.plot(ax=ax1c,color='white', edgecolor='black') 
            SW_rignotetal.plot(ax=ax1c,color='white', edgecolor='black') 
            CW_rignotetal.plot(ax=ax1c,color='white', edgecolor='black') 
            NW_rignotetal.plot(ax=ax1c,color='white', edgecolor='black')
            
            #Flightlines
            # --- 2011-2012
            plt.scatter(flightlines_20022018[flightlines_20022018.str_year==2011]['lon_3413'],
                        flightlines_20022018[flightlines_20022018.str_year==2011]['lat_3413'],
                        s=1,c='#d9d9d9')
            plt.scatter(flightlines_20022018[flightlines_20022018.str_year==2012]['lon_3413'],
                        flightlines_20022018[flightlines_20022018.str_year==2012]['lat_3413'],
                        s=1,c='#d9d9d9',label='flightlines 2011-2012')
            
            # --- 2017-2018
            plt.scatter(flightlines_20022018[flightlines_20022018.str_year==2017]['lon_3413'],
                        flightlines_20022018[flightlines_20022018.str_year==2017]['lat_3413'],
                        s=1,c='#969696')
            plt.scatter(flightlines_20022018[flightlines_20022018.str_year==2018]['lon_3413'],
                        flightlines_20022018[flightlines_20022018.str_year==2018]['lat_3413'],
                        s=1,c='#969696',label='flightlines 2017-2018')
            
            #Likelihood
            # --- 2011-2012
            plt.scatter(df_thickness_likelihood_20102018[df_thickness_likelihood_20102018.Track_name.str[:4]=='2011']['lon_3413'],
                        df_thickness_likelihood_20102018[df_thickness_likelihood_20102018.Track_name.str[:4]=='2011']['lat_3413'],
                        c=df_thickness_likelihood_20102018[df_thickness_likelihood_20102018.Track_name.str[:4]=='2011']['likelihood'],
                        s=2,cmap=plt.get_cmap('Blues'))
            
            lik_blues=plt.scatter(df_thickness_likelihood_20102018[df_thickness_likelihood_20102018.Track_name.str[:4]=='2012']['lon_3413'],
                        df_thickness_likelihood_20102018[df_thickness_likelihood_20102018.Track_name.str[:4]=='2012']['lat_3413'],
                        c=df_thickness_likelihood_20102018[df_thickness_likelihood_20102018.Track_name.str[:4]=='2012']['likelihood'],
                        s=2,cmap=plt.get_cmap('Blues'),label='Likelihood 2011-2012')
            
            # --- 2017-2018            
            plt.scatter(df_thickness_likelihood_20102018[df_thickness_likelihood_20102018.Track_name.str[:4]=='2017']['lon_3413'],
                        df_thickness_likelihood_20102018[df_thickness_likelihood_20102018.Track_name.str[:4]=='2017']['lat_3413'],
                        c=df_thickness_likelihood_20102018[df_thickness_likelihood_20102018.Track_name.str[:4]=='2017']['likelihood'],
                        s=1.5,cmap=plt.get_cmap('Reds'))
            lik_reds=plt.scatter(df_thickness_likelihood_20102018[df_thickness_likelihood_20102018.Track_name.str[:4]=='2018']['lon_3413'],
                        df_thickness_likelihood_20102018[df_thickness_likelihood_20102018.Track_name.str[:4]=='2018']['lat_3413'],
                        c=df_thickness_likelihood_20102018[df_thickness_likelihood_20102018.Track_name.str[:4]=='2018']['likelihood'],
                        s=1.5,cmap=plt.get_cmap('Reds'),label='Likelihood 2017-2018')
            
            '''
            # Plot legend. This is from https://stackoverflow.com/questions/24706125/setting-a-fixed-size-for-points-in-legend
            lgnd = plt.legend(loc="best", scatterpoints=1, fontsize=10)
            lgnd.legendHandles[0]._sizes = [30]
            lgnd.legendHandles[1]._sizes = [30]
            lgnd.legendHandles[2]._sizes = [30]
            lgnd.legendHandles[3]._sizes = [30]
            '''         
            pdb.set_trace()
            
            '''
            import matplotlib.patches as patches
            from matplotlib.patches import Patch
            
            #Custom legend myself
            legend_elements = [Patch(facecolor='#d9d9d9',label='flightlines 2011-2012'),
                               Patch(facecolor='#969696',label='flightlines 2017-2018'),
                               Patch(facecolor='#2171b5',label='Likelihood 2011-2012'),
                               Patch(facecolor='#cb181d',label='Likelihood 2017-2018')]

            ax1c.legend(handles=legend_elements)
            '''
            plt.legend()
    
            ax1c.set_xlabel('Easting [m]')
            ax1c.set_ylabel('Northing [m]')
            
            cbar_reds = plt.colorbar(lik_reds,location = 'right')
            cbar_reds.set_label('Columnal average likelihood - 2017-2018')
            cbar_blues = plt.colorbar(lik_blues,location = 'left')
            cbar_blues.set_label('Columnal average likelihood - 2011-2012')
            
            #Save figure according to different regions
            for region in list(['NE','NO','NW','CW','SW']):
                if (region =='NE'):
                    polygon_for_zoom=NE_rignotetal
                    ax1c.set_xlim(283000,670000)
                    ax1c.set_ylim(-1880000,-940000)
                    ax1c.text(350000,-1600000,'NE',fontsize=25)
                    
                    ################ DISPLAY AREA CHANGE #####################
                    low_end_change=(int((low_end_summary['2017-2018'][region]-low_end_summary['2011-2012'][region])/low_end_summary['2011-2012'][region]*100))
                    high_end_change=(int((high_end_summary['2017-2018'][region]-high_end_summary['2011-2012'][region])/high_end_summary['2011-2012'][region]*100))
                    
                    #Compute and display relative change
                    ax1c.text(350000,-1650000,'[+'+str(low_end_change)+' : +'+str(high_end_change)+'] %')
                    ################ DISPLAY AREA CHANGE #####################
                    
                elif(region =='NO'):
                    polygon_for_zoom=NO_rignotetal
                    ax1c.set_xlim(-605000,302000)
                    ax1c.set_ylim(-1215000,-785000)
                    ax1c.text(-90000,-1075000,'NO',fontsize=25)
                    
                    ################ DISPLAY AREA CHANGE #####################
                    low_end_change=(int((low_end_summary['2017-2018'][region]-low_end_summary['2011-2012'][region])/low_end_summary['2011-2012'][region]*100))
                    high_end_change=(int((high_end_summary['2017-2018'][region]-high_end_summary['2011-2012'][region])/high_end_summary['2011-2012'][region]*100))
                    
                    #Compute and display relative change
                    ax1c.text(-90000,-1025000,'[+'+str(low_end_change)+' : +'+str(high_end_change)+'] %')
                    ################ DISPLAY AREA CHANGE #####################
                    
                elif(region =='NW'):
                    polygon_for_zoom=NW_rignotetal
                    ax1c.set_xlim(-610000,-189000)
                    ax1c.set_ylim(-1140000,-1985000)
                    ax1c.text(-300000,-1410000,'NW',fontsize=25)
                    
                    ################ DISPLAY AREA CHANGE #####################
                    low_end_change=(int((low_end_summary['2017-2018'][region]-low_end_summary['2011-2012'][region])/low_end_summary['2011-2012'][region]*100))
                    high_end_change=(int((high_end_summary['2017-2018'][region]-high_end_summary['2011-2012'][region])/high_end_summary['2011-2012'][region]*100))
                    
                    #Compute and display relative change
                    ax1c.text(-300000,-1460000,'[+'+str(low_end_change)+' : +'+str(high_end_change)+'] %')
                    ################ DISPLAY AREA CHANGE #####################
                    
                elif(region =='CW'):
                    polygon_for_zoom=CW_rignotetal
                    ax1c.set_xlim(-259000,-60500)
                    ax1c.set_ylim(-2385000,-1935000)
                    ax1c.text(-130000,-2061500,'CW',fontsize=25)
                    
                    ################ DISPLAY AREA CHANGE #####################
                    low_end_change=(int((low_end_summary['2017-2018'][region]-low_end_summary['2011-2012'][region])/low_end_summary['2011-2012'][region]*100))
                    high_end_change=(int((high_end_summary['2017-2018'][region]-high_end_summary['2011-2012'][region])/high_end_summary['2011-2012'][region]*100))
                    
                    #Compute and display relative change
                    ax1c.text(-130000,-2101500,'[+'+str(low_end_change)+' : +'+str(high_end_change)+'] %')
                    ################ DISPLAY AREA CHANGE #####################
                    
                elif(region =='SW'):
                    polygon_for_zoom=SW_rignotetal
                    ax1c.set_xlim(-265000,-55600)
                    ax1c.set_ylim(-2899000,-2370000)
                    ax1c.text(-200000,-2775000,'SW',fontsize=25)
                    
                    ################ DISPLAY AREA CHANGE #####################
                    low_end_change=(int((low_end_summary['2017-2018'][region]-low_end_summary['2011-2012'][region])/low_end_summary['2011-2012'][region]*100))
                    high_end_change=(int((high_end_summary['2017-2018'][region]-high_end_summary['2011-2012'][region])/high_end_summary['2011-2012'][region]*100))
                    
                    #Compute and display relative change
                    ax1c.text(-200000,-2825000,'[+'+str(low_end_change)+' : +'+str(high_end_change)+'] %')
                    ################ DISPLAY AREA CHANGE #####################
                    
                else:
                    print('Region not known')
                            
                if (panel_c_save == 'TRUE'):
                    #Save the figure
                    plt.savefig('C:/Users/jullienn/switchdrive/Private/research/RT1/figures/fig1/v2/fig1_panel_c_'+region+'.png',dpi=500)
                    
        # -------------------------------- PANELS C -------------------------------        
    
    pdb.set_trace()

#Import packages
#import rasterio
#from rasterio.plot import show
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pandas as pd
from os import listdir
from os.path import isfile, join
import pickle
#from pysheds.grid import Grid
import pdb
import numpy as np
from pyproj import Transformer
import matplotlib.gridspec as gridspec
import scipy.io
from osgeo import gdal
import geopandas as gpd  # Requires the pyshp package

from matplotlib.colors import ListedColormap, BoundaryNorm
from shapely.geometry import Point, Polygon
from matplotlib.patches import Patch

create_elevation_dictionaries='FALSE'
#pdb.set_trace()
'''
########################## Load GrIS elevation ##########################
#Open the DEM
grid = Grid.from_raster("C:/Users/jullienn/switchdrive/Private/research/backup_Aglaja/working_environment/greenland_topo_data/elevations/greenland_dem_mosaic_100m_v3.0.tif",data_name='dem')
#Minnor slicing on borders to enhance colorbars
elevDem=grid.dem[:-1,:-1]              
#Scale the colormap
#divnorm = mcolors.DivergingNorm(vmin=0, vcenter=1250, vmax=2500)
########################## Load GrIS elevation ##########################
'''

############################ Load DEM information ############################
#Extract elevation from DEM to associated with coordinates. This piece of code
#is from https://gis.stackexchange.com/questions/221292/retrieve-pixel-value-with-geographic-coordinate-as-input-with-gdal
driver = gdal.GetDriverByName('GTiff')
filename_raster = "C:/Users/jullienn/switchdrive/Private/research/backup_Aglaja/working_environment/greenland_topo_data/elevations/greenland_dem_mosaic_100m_v3.0.tif" #path to raster

dataset_dem = gdal.Open(filename_raster)
band = dataset_dem.GetRasterBand(1)

cols = dataset_dem.RasterXSize
rows = dataset_dem.RasterYSize

transform_elev = dataset_dem.GetGeoTransform()

xOrigin = transform_elev[0]
yOrigin = transform_elev[3]
pixelWidth = transform_elev[1]
pixelHeight = -transform_elev[5]

data_dem = band.ReadAsArray(0, 0, cols, rows)

#Define ther zero for longitude:
#Where lon==0 is in may09_03_15:
    #lon[886]=23.53372773084396 and lon[887]=-40.08804568537925
    #lat[886]=-3120053.856912824, lat[887]=-3120048.666364133
avg_lon_zero=(23.53372773084396+-40.08804568537925)/2
index_lon_zero=int((avg_lon_zero-xOrigin) / pixelWidth)
############################ Load DEM information ############################

### -------------------------- Load shapefiles --------------------------- ###
#from https://gis.stackexchange.com/questions/113799/how-to-read-a-shapefile-in-python
path_IceBridgeArea_Shape='C:/Users/jullienn/switchdrive/Private/research/backup_Aglaja/working_environment/greenland_topo_data/IceBridge Area Shapefiles/IceBridge Area Shapefiles/'
IceBridgeArea_Shape=gpd.read_file(path_IceBridgeArea_Shape+'IceBridgeArea_Shape.shp')

path_regional_masks='C:/Users/jullienn/switchdrive/Private/research/backup_Aglaja/working_environment/greenland_topo_data/masks_for_2002_2003_calculations'

NW_icecap_greenland_mask=gpd.read_file(path_regional_masks+'/NW_icecap_greenland_mask_3413.shp')
NW_north_greenland_mask=gpd.read_file(path_regional_masks+'/NW_north_greenland_mask_3413.shp')
NW_west_greenland_mask=gpd.read_file(path_regional_masks+'/NW_west_greenland_mask_3413.shp')
SW_lower_greenland_mask=gpd.read_file(path_regional_masks+'/SW_lower_greenland_mask_3413.shp')
SW_middle_greenland_mask=gpd.read_file(path_regional_masks+'/SW_middle_greenland_mask_3413.shp')
SW_upper_greenland_mask=gpd.read_file(path_regional_masks+'/SW_upper_greenland_mask_3413.shp')

#Load Rignot et al., 2016 Greenland drainage bassins
path_rignotetal2016_GrIS_drainage_bassins='C:/Users/jullienn/switchdrive/Private/research/backup_Aglaja/working_environment/greenland_topo_data/GRE_Basins_IMBIE2_v1.3/'
GrIS_drainage_bassins=gpd.read_file(path_rignotetal2016_GrIS_drainage_bassins+'GRE_Basins_IMBIE2_v1.3_EPSG_3413.shp',rows=slice(51,57,1)) #the regions are the last rows of the shapefile

#Extract indiv regions and create related indiv shapefiles
NO_rignotetal=GrIS_drainage_bassins[GrIS_drainage_bassins.SUBREGION1=='NO']
NE_rignotetal=GrIS_drainage_bassins[GrIS_drainage_bassins.SUBREGION1=='NE']
SE_rignotetal=GrIS_drainage_bassins[GrIS_drainage_bassins.SUBREGION1=='SE']
SW_rignotetal=GrIS_drainage_bassins[GrIS_drainage_bassins.SUBREGION1=='SW']
CW_rignotetal=GrIS_drainage_bassins[GrIS_drainage_bassins.SUBREGION1=='CW']
NW_rignotetal=GrIS_drainage_bassins[GrIS_drainage_bassins.SUBREGION1=='NW']
### -------------------------- Load shapefiles --------------------------- ###

if (create_elevation_dictionaries == 'TRUE'):
    
    ################# Load 2002-2003 flightlines coordinates ################
    path_data='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/icelens_identification'
    
    #Open the file and read it
    f_flightlines = open(path_data+'/metadata_coord_2002_2003', "rb")
    all_2002_3_flightlines = pickle.load(f_flightlines)
    f_flightlines.close()
    ################# Load 2002-2003 flightlines coordinates ################
    
    lat_all=[]
    lon_all=[]
    #elev_all=[]
    
    elevation_dictionnary = {k: {} for k in list(['2002','2003'])}
    
    for year in list(all_2002_3_flightlines.keys()):
        
        elevation_dictionnary[year]={k: {} for k in list(all_2002_3_flightlines[year].keys())}
        
        for days in list(all_2002_3_flightlines[year].keys()):
            
            elevation_dictionnary[year][days]={k: {} for k in list(all_2002_3_flightlines[year][days].keys())}
            
            for indiv_file in list(all_2002_3_flightlines[year][days].keys()):
                if (indiv_file[0:7]=='quality'):
                    continue
                else:
                    print(indiv_file)
                    lat_all=np.append(lat_all,all_2002_3_flightlines[year][days][indiv_file][0])
                    lon_all=np.append(lon_all,all_2002_3_flightlines[year][days][indiv_file][1])
                    #Extract the elevation:
                    lat_elev=[]
                    lon_elev=[]
                    if (days=='jun04'):
                        lat_elev=np.transpose(all_2002_3_flightlines[year][days][indiv_file][0])
                        lon_elev=np.transpose(all_2002_3_flightlines[year][days][indiv_file][1])
                    else:
                        lat_elev=all_2002_3_flightlines[year][days][indiv_file][0]
                        lon_elev=all_2002_3_flightlines[year][days][indiv_file][1]
                    
                    latlon_tuple=[]
                    latlon_tuple=list(zip(lon_elev,lat_elev))
                    
                    elev_indiv_file=[]
                    for indiv_coord in latlon_tuple:
                        if (np.isnan(indiv_coord[0]) or np.isnan(indiv_coord[1])):
                            #elev_all=np.append(elev_all,np.nan)
                            elev_indiv_file=np.append(elev_indiv_file,np.nan)
                        else:
                            #The origin is top left corner!!
                            #y will always be negative
                            row = int((yOrigin - indiv_coord[1] ) / pixelHeight)
                            if (indiv_coord[0]<0):
                                # if x negative
                                col = index_lon_zero-int((-indiv_coord[0]-0) / pixelWidth)
                            elif (indiv_coord[0]>0):
                                # if x positive
                                col = index_lon_zero+int((indiv_coord[0]-0) / pixelWidth)
                            #Read the elevation
                            #elev_all=np.append(elev_all,data_dem[row][col])
                            elev_indiv_file=np.append(elev_indiv_file,data_dem[row][col])
                    
                    #Store data into the dictionnary
                    elevation_dictionnary[year][days][indiv_file]=elev_indiv_file
                    
    ################# Load 2002-2003 flightlines coordinates ################
    
    ################### Load 2002-2003 ice lenses location ##################
    #Open the file and read it
    f_icelens_flightlines = open(path_data+'/metadata_coord_icelens_2002_2003_26022020', "rb")
    icelens_2002_3_flightlines = pickle.load(f_icelens_flightlines)
    f_icelens_flightlines.close()
    
    lat_icelens=[]
    lon_icelens=[]
    colorcode_icelens=[]
    Track_name=[]
    
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
                    colorcode_icelens=np.append(colorcode_icelens,icelens_2002_3_flightlines[year][days][indiv_file][2])
                    #Create an empty vector of strings
                    Track_name=np.append(Track_name,[indiv_file for x in range(0,len(icelens_2002_3_flightlines[year][days][indiv_file][0]))])
    
    #Create a dataframe out of it
    df_2002_2003=pd.DataFrame(lat_icelens, columns =['lat_3413'])
    df_2002_2003['lon_3413']=lon_icelens
    df_2002_2003['colorcode_icelens']=colorcode_icelens
    df_2002_2003['Track_name']=Track_name
    ################### Load 2002-2003 ice lenses location ##################
    
    ################### Load 2010-2018 ice slabs location ##################
    
    #Load the data
    filename_20102018='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2010_2018/excel_spatial_aggreation_and_other/final_excel/low_estimate/Ice_Layer_Output_Thicknesses_2010_2018_jullienetal2021_low_estimate.csv'
    df_20102018 = pd.read_csv(filename_20102018, sep=",", decimal='.')
    
    #Transform the coordinated from WGS84 to EPSG:3413
    #Example from: https://pyproj4.github.io/pyproj/stable/examples.html
    transformer = Transformer.from_crs("EPSG:4326", "EPSG:3413", always_xy=True)
    points=transformer.transform(np.array(df_20102018.lon),np.array(df_20102018.lat))
    
    lon_3413_20102018=points[0]
    lat_3413_20102018=points[1]
    ################### Load 2010-2018 ice slabs location ##################
    
    #######################################################################
    ###          Inland expansion of iceslabs from 2002 to 2018         ###
    #######################################################################
    
    #Prepare plot
    fig, (ax1) = plt.subplots(1, 1)#, gridspec_kw={'width_ratios': [1, 3]})
    fig.suptitle('Iceslabs area overview')
    
    #Display DEM
    cb1=ax1.imshow(elevDem, extent=grid.extent,cmap=discrete_cmap(10,'cubehelix_r'),alpha=0.5)#,norm=divnorm)
    cbar1=fig.colorbar(cb1, ax=[ax1], location='left')
    cbar1.set_label('Elevation [m]')
    
    #Display the shapefile
    IceBridgeArea_Shape.plot(ax=ax1)
    
    #Plot all the 2002-2003 icelenses according to their condifence color
    #1. Red
    ax1.scatter(df_2002_2003[df_2002_2003['colorcode_icelens']==-1]['lon_3413'], df_2002_2003[df_2002_2003['colorcode_icelens']==-1]['lat_3413'],s=1,facecolors='#c9662c', edgecolors='none')
    #2. Orange
    ax1.scatter(df_2002_2003[df_2002_2003['colorcode_icelens']==0]['lon_3413'], df_2002_2003[df_2002_2003['colorcode_icelens']==0]['lat_3413'],s=1,facecolors='#fed976', edgecolors='none')
    #3. Green
    ax1.scatter(df_2002_2003[df_2002_2003['colorcode_icelens']==1]['lon_3413'], df_2002_2003[df_2002_2003['colorcode_icelens']==1]['lat_3413'],s=1,facecolors='#238b45', edgecolors='none')
    ##Purple
    #ax1.scatter(lon_icelens[colorcode_icelens==2], lat_icelens[colorcode_icelens==2],s=1,facecolors='purple', edgecolors='none')
    
    #Correct zoom
    ax1.set_xlim(-650000,900000)
    ax1.set_ylim(-3360000,-650000)
    
    plt.show()

    # compare min and max of lat/lon of the track with respect to shapefile
    
    #Store lat/lon 3413
    df_20102018['lat_3413']=lat_3413_20102018
    df_20102018['lon_3413']=lon_3413_20102018
    
    #Initialise the elevation and shapefile belonging column
    df_20102018['key_shp']=np.nan
    df_20102018['elevation']=np.nan
    df_20102018['year']=np.nan
    
    #I. Load regional shapefile I have created on QGIS
    #Done before the if statement
    
    #II. Do the intersection between the mask and 2010-2018 data and keep only the matching one
    
    #This part of code is from 'refine_location_2017_2018.py'
    #Loop over all data point to check whether it belongs to one of the four shapefile
    for i in range(0,lon_3413_20102018.size):
        #select the point i
        single_point=Point(lon_3413_20102018[i],lat_3413_20102018[i])
        
        #Do the identification between the point i and the regional shapefiles
        #From: https://automating-gis-processes.github.io/CSC18/lessons/L4/point-in-polygon.html
        check_NO_rignotetal=np.asarray(NO_rignotetal.contains(single_point)).astype(int)
        check_NE_rignotetal=np.asarray(NE_rignotetal.contains(single_point)).astype(int)
        check_SE_rignotetal=np.asarray(SE_rignotetal.contains(single_point)).astype(int)
        check_SW_rignotetal=np.asarray(SW_rignotetal.contains(single_point)).astype(int)
        check_CW_rignotetal=np.asarray(CW_rignotetal.contains(single_point)).astype(int)
        check_NW_rignotetal=np.asarray(NW_rignotetal.contains(single_point)).astype(int)
    
        #Associated the point of interest to its regional shapefile in data_iceslabs
        if (np.sum(check_NO_rignotetal)>0):
            df_20102018['key_shp'][i]='NO'
        elif (np.sum(check_NE_rignotetal)>0):
            df_20102018['key_shp'][i]='NE'
        elif (np.sum(check_SE_rignotetal)>0):
            df_20102018['key_shp'][i]='SE'
        elif (np.sum(check_SW_rignotetal)>0):
            df_20102018['key_shp'][i]='SW'
        elif (np.sum(check_CW_rignotetal)>0):
            df_20102018['key_shp'][i]='CW'
        elif (np.sum(check_NW_rignotetal)>0):
            df_20102018['key_shp'][i]='NW'
        else:
            df_20102018['key_shp'][i]='Out'
        
        #Calculate the corresponding elevation
        df_20102018['elevation'][i]=calcul_elevation(df_20102018['lon_3413'][i],df_20102018['lat_3413'][i],data_dem,yOrigin,pixelHeight,pixelWidth,index_lon_zero)
        #Add the year
        df_20102018['year'][i]=int(df_20102018['Track_name'][i][0:4])
    
        #Monitor the process
        print(i/lon_3413_20102018.size*100,'%')
    
    #Save the dictionary into a picke file
    filename_tosave='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2010_2018/excel_spatial_aggreation_and_other/df_20102018_with_elevation_low_estimate_rignotetalregions'
    outfile= open(filename_tosave, "wb" )
    pickle.dump(df_20102018,outfile)
    outfile.close()
    
    '''
    #Display the keys
    fig, (ax1) = plt.subplots(1, 1)#, gridspec_kw={'width_ratios': [1, 3]})
    fig.suptitle('Iceslabs keys')
    
    #Display DEM
    cb1=ax1.imshow(elevDem, extent=grid.extent,cmap=discrete_cmap(10,'cubehelix_r'),alpha=0.5,norm=divnorm)
    cbar1=fig.colorbar(cb1, ax=[ax1], location='left')
    cbar1.set_label('Elevation [m]')
    
    #Display the shapefile
    check_NO_rignotetal.plot(ax=ax1)
    check_NE_rignotetal.plot(ax=ax1)
    check_SE_rignotetal.plot(ax=ax1)
    check_SW_rignotetal.plot(ax=ax1)
    check_CW_rignotetal.plot(ax=ax1)
    check_NW_rignotetal.plot(ax=ax1)
    
    #Display the data as a function of their belonging keys
    ax1.scatter(df_20102018['lon_3413'][df_20102018['key_shp']=='NO'],df_20102018['lat_3413'][df_20102018['key_shp']=='NO'],facecolors='orange')
    ax1.scatter(df_20102018['lon_3413'][df_20102018['key_shp']=='NE'],df_20102018['lat_3413'][df_20102018['key_shp']=='NE'],facecolors='blue')
    ax1.scatter(df_20102018['lon_3413'][df_20102018['key_shp']=='SE'],df_20102018['lat_3413'][df_20102018['key_shp']=='SE'],facecolors='purple')
    ax1.scatter(df_20102018['lon_3413'][df_20102018['key_shp']=='SW'],df_20102018['lat_3413'][df_20102018['key_shp']=='SW'],facecolors='red')
    ax1.scatter(df_20102018['lon_3413'][df_20102018['key_shp']=='CW'],df_20102018['lat_3413'][df_20102018['key_shp']=='CW'],facecolors='green')
    ax1.scatter(df_20102018['lon_3413'][df_20102018['key_shp']=='NW'],df_20102018['lat_3413'][df_20102018['key_shp']=='NW'],facecolors='k')
    '''
    #III. Do the intersection between the mask and 2002-2003 data
    
    #Initialise the shapefile belonging column
    df_2002_2003['key_shp']=np.nan
    df_2002_2003['elevation']=np.nan
    df_2002_2003['year']=np.nan
    
    #This part of code is from 'refine_location_2017_2018.py'
    #Loop over all data point to check whether it belongs to one of the four shapefile
    for i in range(0,len(df_2002_2003)):
        #select the point i
        single_point=Point(df_2002_2003['lon_3413'].iloc[i],df_2002_2003['lat_3413'].iloc[i])
        
        #Do the identification between the point i and the regional shapefiles
        #From: https://automating-gis-processes.github.io/CSC18/lessons/L4/point-in-polygon.html
        check_NO_rignotetal=np.asarray(NO_rignotetal.contains(single_point)).astype(int)
        check_NE_rignotetal=np.asarray(NE_rignotetal.contains(single_point)).astype(int)
        check_SE_rignotetal=np.asarray(SE_rignotetal.contains(single_point)).astype(int)
        check_SW_rignotetal=np.asarray(SW_rignotetal.contains(single_point)).astype(int)
        check_CW_rignotetal=np.asarray(CW_rignotetal.contains(single_point)).astype(int)
        check_NW_rignotetal=np.asarray(NW_rignotetal.contains(single_point)).astype(int)
        
        #Associated the point of interest to its regional shapefile in data_iceslabs
        if (np.sum(check_NO_rignotetal)>0):
            df_2002_2003['key_shp'].iloc[i]='NO'
        elif (np.sum(check_NE_rignotetal)>0):
            df_2002_2003['key_shp'].iloc[i]='NE'
        elif (np.sum(check_SE_rignotetal)>0):
            df_2002_2003['key_shp'].iloc[i]='SE'
        elif (np.sum(check_SW_rignotetal)>0):
            df_2002_2003['key_shp'].iloc[i]='SW'
        elif (np.sum(check_CW_rignotetal)>0):
            df_2002_2003['key_shp'].iloc[i]='CW'
        elif (np.sum(check_NW_rignotetal)>0):
            df_2002_2003['key_shp'].iloc[i]='NW'
        else:
            df_2002_2003['key_shp'].iloc[i]='Out'
        
        #Calculate the corresponding elevation
        df_2002_2003['elevation'].iloc[i]=calcul_elevation(df_2002_2003['lon_3413'].iloc[i],df_2002_2003['lat_3413'].iloc[i],data_dem,yOrigin,pixelHeight,pixelWidth,index_lon_zero)
        
        #Add the year
        if (df_2002_2003['Track_name'][i][6:8] == '02'):
            year_to_write=2002
        elif (df_2002_2003['Track_name'][i][6:8] == '03'):
            year_to_write=2003
        else:
            print('Year not known, error')
            break
        
        df_2002_2003['year'][i]=year_to_write
        
        #Monitor the process
        print(i/len(df_2002_2003)*100,'%')
    
    #Only work with green slabs
    df_2002_2003_green=df_2002_2003[df_2002_2003['colorcode_icelens']==1]
    '''
    #Display the data as a function of their belonging keys
    ax1.scatter(df_2002_2003_green['lon_3413'][df_2002_2003_green['key_shp']=='NO'],df_2002_2003_green['lat_3413'][df_2002_2003_green['key_shp']=='NO'],facecolors='brown')
    ax1.scatter(df_2002_2003_green['lon_3413'][df_2002_2003_green['key_shp']=='NE'],df_2002_2003_green['lat_3413'][df_2002_2003_green['key_shp']=='NE'],facecolors='cyan')
    ax1.scatter(df_2002_2003_green['lon_3413'][df_2002_2003_green['key_shp']=='SE'],df_2002_2003_green['lat_3413'][df_2002_2003_green['key_shp']=='SE'],facecolors='pink')
    ax1.scatter(df_2002_2003_green['lon_3413'][df_2002_2003_green['key_shp']=='SW'],df_2002_2003_green['lat_3413'][df_2002_2003_green['key_shp']=='SW'],facecolors='yellow')
    ax1.scatter(df_2002_2003_green['lon_3413'][df_2002_2003_green['key_shp']=='CW'],df_2002_2003_green['lat_3413'][df_2002_2003_green['key_shp']=='CW'],facecolors='olive')
    ax1.scatter(df_2002_2003_green['lon_3413'][df_2002_2003_green['key_shp']=='NW'],df_2002_2003_green['lat_3413'][df_2002_2003_green['key_shp']=='NW'],facecolors='gray')
    
    plt.show()
    '''
    #Save the dictionary into a picke file
    filename_tosave='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2010_2018/excel_spatial_aggreation_and_other/df_2002_2003_with_elevation_rignotetalregions'
    outfile= open(filename_tosave, "wb" )
    pickle.dump(df_2002_2003,outfile)
    outfile.close()
    
else:
    #Dictionnaries have already been created, load them
    path_df_with_elevation='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2010_2018/excel_spatial_aggreation_and_other/' 
    
    #Load 2002-2003
    f_20022003 = open(path_df_with_elevation+'df_2002_2003_with_elevation_prob00_rignotetalregions', "rb")
    #f_20022003 = open(path_df_with_elevation+'df_2002_2003_with_elevation_prob00', "rb")
    df_2002_2003 = pickle.load(f_20022003)
    f_20022003.close()
    
    #Load 2010-2018
    f_20102018 = open(path_df_with_elevation+'df_20102018_with_elevation_prob00_rignotetalregions', "rb")
    #f_20102018 = open(path_df_with_elevation+'df_20102018_with_elevation_prob00', "rb")
    df_2010_2018 = pickle.load(f_20102018)
    f_20102018.close()
    
    #Load 2010-2018 high estimate
    f_20102018_high = open(path_df_with_elevation+'df_20102018_with_elevation_high_estimate_rignotetalregions', "rb")
    df_2010_2018_high = pickle.load(f_20102018_high)
    f_20102018_high.close()

    #Load 2010-2018 low estimate
    f_20102018_low = open(path_df_with_elevation+'df_20102018_with_elevation_low_estimate_rignotetalregions', "rb")
    df_2010_2018_low = pickle.load(f_20102018_low)
    f_20102018_low.close()

#Load Miege firn aquifer
path_firn_aquifer='C:/Users/jullienn/switchdrive/Private/research/backup_Aglaja/working_environment/greenland_topo_data/firn_aquifers_miege/'
df_firn_aquifer_2010 = pd.read_csv(path_firn_aquifer+'MiegeFirnAquiferDetections2010.csv',delimiter=',',decimal='.')
df_firn_aquifer_2011 = pd.read_csv(path_firn_aquifer+'MiegeFirnAquiferDetections2011.csv',delimiter=',',decimal='.')
df_firn_aquifer_2012 = pd.read_csv(path_firn_aquifer+'MiegeFirnAquiferDetections2012.csv',delimiter=',',decimal='.')
df_firn_aquifer_2013 = pd.read_csv(path_firn_aquifer+'MiegeFirnAquiferDetections2013.csv',delimiter=',',decimal='.')
df_firn_aquifer_2014 = pd.read_csv(path_firn_aquifer+'MiegeFirnAquiferDetections2014.csv',delimiter=',',decimal='.')

#Append all the miege aquifer files to each other
df_firn_aquifer_all=df_firn_aquifer_2010
df_firn_aquifer_all=df_firn_aquifer_all.append(df_firn_aquifer_2011)
df_firn_aquifer_all=df_firn_aquifer_all.append(df_firn_aquifer_2012)
df_firn_aquifer_all=df_firn_aquifer_all.append(df_firn_aquifer_2013)
df_firn_aquifer_all=df_firn_aquifer_all.append(df_firn_aquifer_2014)

#Transform miege coordinates from WGS84 to EPSG:3413
transformer = Transformer.from_crs("EPSG:4326", "EPSG:3413", always_xy=True)
points=transformer.transform(np.asarray(df_firn_aquifer_all["LONG"]),np.asarray(df_firn_aquifer_all["LAT"]))

#Store lat/lon in 3413
df_firn_aquifer_all['lon_3413']=points[0]
df_firn_aquifer_all['lat_3413']=points[1]


#Load columnal likelihood file likelihood
path_thickness_likelihood='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2010_2018/excel_spatial_aggreation_and_other/final_excel/high_estimate_and_columnal_likelihood/'
df_thickness_likelihood_20102018 = pd.read_csv(path_thickness_likelihood+'Ice_Layer_Output_Thicknesses_Likelihood_2010_2018_jullienetal2021.csv',delimiter=',',decimal='.')
#Transform miege coordinates from WGS84 to EPSG:3413
points=transformer.transform(np.asarray(df_thickness_likelihood_20102018["lon"]),np.asarray(df_thickness_likelihood_20102018["lat"]))

#Store lat/lon in 3413
df_thickness_likelihood_20102018['lon_3413']=points[0]
df_thickness_likelihood_20102018['lat_3413']=points[1]


#IV. From here on, work with the different periods separated by strong melting summers.
#    Work thus with 2002-2003 VS 2010 VS 2011-2012 VS 2013-2014 VS 2017-2018
#    Select the absolute low and absolute high of 2002-2003, 2010-2014 and 2017-2018

#Let's create ~10km latitudinal (resp. longitudinal) slices for SW Greenland (resp. NW Greenland)
#and calculate the low and high end in each slice for elevation difference:

#1. Create the latitudinal (resp. longitudinal) slices
############ Is this grid to change?? This is about 10km width but not even whether north or south!
lat_slices=np.linspace(-3000000,-1150000,int((np.abs(-3000000)-np.abs(-1150000))/10000))
lon_slices=np.linspace(-600000,650000,int((np.abs(650000)+np.abs(-600000))/10000))

#2. Select and store all the data belonging to the lon/lat slices in a dictionnary.
#### ------------------------- 2010-2018 -------------------------------- ####
#   Retreive and store min and max elevation of each slice in a dataframe
#   ----- Latitudinal slices

#Create a dictionnary where to store slices information
dict_lat_slice={}

#Create a dictionnary to store np arrays storing slices min and max elevation for each region
dict_lat_slices_summary={k: {} for k in list(df_2010_2018['key_shp'].unique())}

#loop over the regions, create the room for each time period in each region
for region in list(df_2010_2018['key_shp'].unique()):
    dict_lat_slices_summary[region]={k: {} for k in list(['2010','2011-2012','2013-2014','2017-2018'])}
    
    for time_period in dict_lat_slices_summary[region].keys():
        #Fill the dict_lat_slices_summary dictionnary with zeros
        dict_lat_slices_summary[region][time_period]=np.zeros((len(lat_slices),2))*np.nan

#Loop over each boundary of lat slices and store dataset related to slices
for i in range(1,len(lat_slices)):
    
    #Identify low and higher end of the slice
    low_bound=lat_slices[i-1]
    high_bound=lat_slices[i]
    
    #Select all the data belonging to this slice
    ind_slice=np.logical_and(np.array(df_2010_2018['lat_3413']>=low_bound),np.array(df_2010_2018['lat_3413']<high_bound))
    df_slice=df_2010_2018[ind_slice]
    
    #Store the associated df
    dict_lat_slice[str(int(lat_slices[i-1]))+' to '+str(int(lat_slices[i]))]=df_slice   
    
    #Loop over the regions present in df_slice
    for region in list(df_slice['key_shp'].unique()):
        #Select only the data belonging to this region
        df_region=df_slice[df_slice['key_shp']==region]
        
        #Loop over the different time periods (2010, 2011-2012, 2013-2014, 2017-2018)
        for time_period in list(['2010','2011-2012','2013-2014','2017-2018']):
            if (time_period == '2010'):
                df_region_period=df_region[df_region['year']==2010]
            elif (time_period == '2011-2012'):
                df_region_period=df_region[(df_region['year']>=2011) & (df_region['year']<=2012)]
            elif (time_period == '2013-2014'):
                df_region_period=df_region[(df_region['year']>=2013) & (df_region['year']<=2014)]
            elif (time_period == '2017-2018'):
                df_region_period=df_region[(df_region['year']>=2017) & (df_region['year']<=2018)]
            else:
                print('Time period not known, break')
                break
            #Identify min and max of each region and store them into a dataframe
            #Retreive the stored array
            array_region_indiv=dict_lat_slices_summary[region][time_period]
            #Store min and max of this regional slice
            array_region_indiv[i,0]=np.min(df_region_period['elevation'])
            array_region_indiv[i,1]=np.max(df_region_period['elevation'])
            #Store again data into dict_lat_slices_summary
            dict_lat_slices_summary[region][time_period]=array_region_indiv
            
#   ----- Longitudinal slices
#Create a dictionnary where to store slices information
dict_lon_slice={}

#Create a dictionnary to store np arrays storing slices min and max elevation for each region
dict_lon_slices_summary={k: {} for k in list(df_2010_2018['key_shp'].unique())}

#loop over the regions, create the room for each time period in each region
for region in list(df_2010_2018['key_shp'].unique()):
    dict_lon_slices_summary[region]={k: {} for k in list(['2010','2011-2012','2013-2014','2017-2018'])}
    
    for time_period in dict_lon_slices_summary[region].keys():
        #Fill the dict_lon_slices_summary dictionnary with zeros
        dict_lon_slices_summary[region][time_period]=np.zeros((len(lon_slices),2))*np.nan

#Loop over each boundary of lon slices and store dataset related to slices
for i in range(1,len(lon_slices)):
    
    #Identify low and higher end of the slice
    low_bound=lon_slices[i-1]
    high_bound=lon_slices[i]
    
    #Select all the data belonging to this slice
    ind_slice=np.logical_and(np.array(df_2010_2018['lon_3413']>=low_bound),np.array(df_2010_2018['lon_3413']<high_bound))
    df_slice=df_2010_2018[ind_slice]
    
    #Store the associated df
    dict_lon_slice[str(int(lon_slices[i-1]))+' to '+str(int(lon_slices[i]))]=df_slice   
    
    #Loop over the regions present in df_slice
    for region in list(df_slice['key_shp'].unique()):
        #Select only the data belonging to this region
        df_region=df_slice[df_slice['key_shp']==region]
        
        #Loop over the different time periods (2010, 2011-2012, 2013-2014, 2017-2018)
        for time_period in list(['2010','2011-2012','2013-2014','2017-2018']):
            if (time_period == '2010'):
                df_region_period=df_region[df_region['year']==2010]
            elif (time_period == '2011-2012'):
                df_region_period=df_region[(df_region['year']>=2011) & (df_region['year']<=2012)]
            elif (time_period == '2013-2014'):
                df_region_period=df_region[(df_region['year']>=2013) & (df_region['year']<=2014)]
            elif (time_period == '2017-2018'):
                df_region_period=df_region[(df_region['year']>=2017) & (df_region['year']<=2018)]
            else:
                print('Time period not known, break')
                break
            #Identify min and max of each region and store them into a dataframe
            #Retreive the stored array
            array_region_indiv=dict_lon_slices_summary[region][time_period]
            #Store min and max of this regional slice
            array_region_indiv[i,0]=np.min(df_region_period['elevation'])
            array_region_indiv[i,1]=np.max(df_region_period['elevation'])
            #Store again data into dict_lat_slices_summary
            dict_lon_slices_summary[region][time_period]=array_region_indiv

#### ------------------------- 2010-2018 -------------------------------- ####

#3. Associate each slice to its belonging region.
#   Not needed! Already present in dataframes!

#4. Calculate the average minimum and maximum of each region among the slices

#5. Flag the more or less perpendicularly crossing 2002-2003 flight lines and exclude the one not crossing
flag_low=['jun04_02proc_4.mat','jun04_02proc_36.mat','jun04_02proc_52.mat','jun04_02proc_53.mat',
      'may09_03_0_aggregated','may09_03_1_aggregated','may09_03_30_aggregated',
      'may09_03_37_aggregated','may11_03_8_aggregated','may11_03_12_aggregated',
      'may11_03_13_aggregated','may11_03_16_aggregated','may11_03_20_aggregated',
      'may11_03_21_aggregated','may11_03_38_aggregated','may11_03_39_aggregated',
      'may12_03_1_aggregated','may12_03_2_aggregated','may12_03_11_aggregated',
      'may12_03_15_aggregated','may12_03_36_aggregated','may13_03_30_aggregated',
      'may14_03_1_aggregated','may14_03_2_aggregated','may14_03_7_aggregated',
      'may14_03_8_aggregated','may14_03_20_aggregated','may14_03_21_aggregated',
      'may15_03_0_aggregated','may15_03_2_aggregated','may15_03_4_aggregated',
      'may15_03_9_aggregated','may18_02_27_aggregated']

flag_high=['jun04_02proc_4.mat','jun04_02proc_36.mat','jun04_02proc_52.mat','jun04_02proc_53.mat',
      'may09_03_0_aggregated','may09_03_1_aggregated','may09_03_30_aggregated',
      'may09_03_37_aggregated','may11_03_20_aggregated','may11_03_21_aggregated',
      'may11_03_37_aggregated','may11_03_38_aggregated','may12_03_1_aggregated',
      'may12_03_2_aggregated','may12_03_11_aggregated','may12_03_36_aggregated',
      'may13_03_30_aggregated','may14_03_1_aggregated','may14_03_2_aggregated',
      'may14_03_7_aggregated','may14_03_20_aggregated','may14_03_21_aggregated',
      'may15_03_2_aggregated','may15_03_4_aggregated','may15_03_9_aggregated',
      'may18_02_27_aggregated']

unique_flags=np.unique(np.append(flag_low,flag_high))

#6. Take the absolute min and max of all 2002-2003 ice slabs in a specific region
#A suite of 2002-2003 traces do not belong to different region, which ease coding

#Here are the traces. For consecutive ones, the ice slabs range elevation is distributed through consecutive traces
traces=[['jun04_02proc_4.mat'],
        ['jun04_02proc_36.mat'],
        ['jun04_02proc_52.mat','jun04_02proc_53.mat'],
        ['may09_03_0_aggregated','may09_03_1_aggregated'],
        ['may09_03_30_aggregated'],
        ['may09_03_37_aggregated'],
        ['may11_03_8_aggregated'],
        ['may11_03_12_aggregated','may11_03_13_aggregated'],
        ['may11_03_16_aggregated'],
        ['may11_03_20_aggregated','may11_03_21_aggregated'],
        ['may11_03_37_aggregated','may11_03_38_aggregated','may11_03_39_aggregated'],
        ['may12_03_1_aggregated','may12_03_2_aggregated'],
        ['may12_03_11_aggregated'],
        ['may12_03_15_aggregated'],
        ['may12_03_36_aggregated'],
        ['may13_03_30_aggregated'],
        ['may14_03_1_aggregated','may14_03_2_aggregated'],
        ['may14_03_7_aggregated','may14_03_8_aggregated'],
        ['may14_03_20_aggregated','may14_03_21_aggregated'],
        ['may15_03_0_aggregated'],
        ['may15_03_2_aggregated'],
        ['may15_03_4_aggregated'],
        ['may15_03_9_aggregated'],
        ['may18_02_27_aggregated']]

list_traces=[item for sublist in traces for item in sublist]

#Loop over the traces, check the flags and populate a low end and high end vector where applicable.
#If consecutive traces, consider the suite of traces!

#Create the dictionary
dict_summary_2002_2003={k: {} for k in list(df_2002_2003['key_shp'].unique())}

#Fill the dict_summary_2002_2003 dictionnary with a np.nan
for region in list(df_2002_2003['key_shp'].unique()):
    dict_summary_2002_2003[region]=np.zeros((len(traces),2))*np.nan

count=0
#Loop over the traces
for trace in traces:
    
    #Check whether we are dealing with single or consecutive traces
    if(len(trace)>1):
        #We are dealing with consecutive traces
        #Select the data related to the first trace
        data_trace=df_2002_2003[df_2002_2003['Track_name']==trace[0]]
        
        #loop over the traces and append data to each other, do not take the first one
        for indiv_trace in list(trace[1:]):
            #Select all the data related to this trace
            data_trace=data_trace.append(df_2002_2003[df_2002_2003['Track_name']==indiv_trace])
            
    else:
        #We are dealing with individual traces
        #Select all the data related to this trace
        data_trace=df_2002_2003[df_2002_2003['Track_name']==trace[0]]

    #Now my data_trace datasets are ready to be worked with
    #Keep only green ice slabs
    data_trace=data_trace[data_trace['colorcode_icelens']==1]
    
    if (len(data_trace)<1):
        #No green ice slabs, continue
        continue
    else:
        #Identify the region
        region=list(np.unique(data_trace['key_shp']))
        
        #Retreive the stored array
        array_region_indiv=dict_summary_2002_2003[region[0]]
        
        #Check the flag: shall we store data?
        if trace[0] in list(flag_low):
            #Store min in array_region_indiv
            array_region_indiv[count,0]=np.min(data_trace['elevation'])
        
        if trace[0] in list(flag_high):
            #Store max in array_region_indiv
            array_region_indiv[count,1]=np.max(data_trace['elevation'])
        
        #Update count
        count=count+1
        #Store again data into dict_lat_slices_summary
        dict_summary_2002_2003[region[0]]=array_region_indiv

#7. Do the elevation difference and eventually the corresponding distance calculation in each region
#Create a dictionnary where to store relevant information
dict_summary={k: {} for k in list(df_2010_2018['key_shp'].unique())}

#Loop over the regions
for region in list(df_2010_2018['key_shp'].unique()):
    
    #Continue building the dictionnary
    dict_summary[region]={k: {} for k in list(['2002-2003','2010','2011-2012','2013-2014','2017-2018'])}
    
    #Loop over the 5 time periods
    
    for time_period in list(['2002-2003','2010','2011-2012','2013-2014','2017-2018']):
        dict_summary[region][time_period]={k: {} for k in list(['max_elev_mean','max_elev_median','max_elev_std'])}
        
        #Take the average, median and std dev of high elevation where ice slabs have been
        #identified in this region, no matter the year in this specific time
        #period, and store relevant information
        
        if (time_period=='2002-2003'):
            #Retreive the corresponding matrix where data are stored
            dict_temp=dict_summary_2002_2003[region]
        else:
            #The dictionnary to select is different whether we are in north or south greenland
            if (region in list(['NO'])):
                dict_temp=dict_lon_slices_summary[region][time_period]
            else:
                dict_temp=dict_lat_slices_summary[region][time_period]
        
        #Calculate and store averages
        dict_summary[region][time_period]['max_elev_mean']=np.nanmean(dict_temp[:,1])
        dict_summary[region][time_period]['max_elev_median']=np.nanmedian(dict_temp[:,1])
        dict_summary[region][time_period]['max_elev_std']=np.nanstd(dict_temp[:,1])

'''
#Plot the inland expansion as a graph
#Display the keys
fig, axs = plt.subplots(2, 3)#, gridspec_kw={'width_ratios': [1, 3]})
fig.suptitle('Iceslabs inland progression')

axs = axs.ravel()

#count for subplot
i=0
#Loop over the region
for region in list(dict_summary.keys()):
    if (region == 'Out'):
        continue
    #Create empty vectors
    low_end=np.zeros(1)
    high_end=np.zeros(1)

    for time_period in list(dict_summary[region].keys()):
        low_end=np.append(low_end,dict_summary[region][time_period]['min_elev'])
        high_end=np.append(high_end,dict_summary[region][time_period]['max_elev'])
    
    #Remove zeros from low_end and high_end vectors
    low_end = low_end[~(low_end==0)]
    high_end = high_end[~(high_end==0)]
    
    #Plot the low end and high end of each region
    axs[i].plot(np.linspace(0,2,len(low_end)),low_end,label='Low end')
    axs[i].plot(np.linspace(0,2,len(high_end)),high_end,label='High end')
    
    #Set title
    axs[i].title.set_text(region)
    
    #Set x tick
    axs[i].set_xticks(np.linspace(0,2,len(high_end)))
    axs[i].set_xticklabels(list(dict_summary[region].keys()))
    
    axs[i].set_xlim(0,2)
    
    axs[i].grid()
    #Update count
    i=i+1
    
plt.legend()
plt.show()
'''
#######################################################################
###          Inland expansion of iceslabs from 2002 to 2018         ###
#######################################################################

#######################################################################
###  Violin plot - Inland expansion of iceslabs from 2002 to 2018   ###
#######################################################################
#Try the violin plot - Do not require any latitudinal and longitudinal averaging!

import seaborn as sns

sns.set_theme(style="whitegrid")

#Select 2002-2003 with green ice slabs only
df_2002_2003_green=df_2002_2003[df_2002_2003['colorcode_icelens']==1]

#Set the year for plotting
df_2002_2003_green['str_year']=["2002-2003" for x in range(len(df_2002_2003_green))]
df_2010_2018.loc[df_2010_2018['year']==2010,'str_year']=["2010" for x in range(len(df_2010_2018[df_2010_2018['year']==2010]))]
df_2010_2018.loc[df_2010_2018['year']==2011,'str_year']=["2011-2012" for x in range(len(df_2010_2018[df_2010_2018['year']==2011]))]
df_2010_2018.loc[df_2010_2018['year']==2012,'str_year']=["2011-2012" for x in range(len(df_2010_2018[df_2010_2018['year']==2012]))]
df_2010_2018.loc[df_2010_2018['year']==2013,'str_year']=["2013-2014" for x in range(len(df_2010_2018[df_2010_2018['year']==2013]))]
df_2010_2018.loc[df_2010_2018['year']==2014,'str_year']=["2013-2014" for x in range(len(df_2010_2018[df_2010_2018['year']==2014]))]
df_2010_2018.loc[df_2010_2018['year']==2017,'str_year']=["2017-2018" for x in range(len(df_2010_2018[df_2010_2018['year']==2017]))]
df_2010_2018.loc[df_2010_2018['year']==2018,'str_year']=["2017-2018" for x in range(len(df_2010_2018[df_2010_2018['year']==2018]))]

#Set the year for plotting in high estimate
df_2010_2018_high.loc[df_2010_2018_high['year']==2010,'str_year']=["2010" for x in range(len(df_2010_2018_high[df_2010_2018_high['year']==2010]))]
df_2010_2018_high.loc[df_2010_2018_high['year']==2011,'str_year']=["2011-2012" for x in range(len(df_2010_2018_high[df_2010_2018_high['year']==2011]))]
df_2010_2018_high.loc[df_2010_2018_high['year']==2012,'str_year']=["2011-2012" for x in range(len(df_2010_2018_high[df_2010_2018_high['year']==2012]))]
df_2010_2018_high.loc[df_2010_2018_high['year']==2013,'str_year']=["2013-2014" for x in range(len(df_2010_2018_high[df_2010_2018_high['year']==2013]))]
df_2010_2018_high.loc[df_2010_2018_high['year']==2014,'str_year']=["2013-2014" for x in range(len(df_2010_2018_high[df_2010_2018_high['year']==2014]))]
df_2010_2018_high.loc[df_2010_2018_high['year']==2017,'str_year']=["2017-2018" for x in range(len(df_2010_2018_high[df_2010_2018_high['year']==2017]))]
df_2010_2018_high.loc[df_2010_2018_high['year']==2018,'str_year']=["2017-2018" for x in range(len(df_2010_2018_high[df_2010_2018_high['year']==2018]))]

#Set the year for plotting in high estimates
df_2010_2018_low.loc[df_2010_2018_low['year']==2010,'str_year']=["2010" for x in range(len(df_2010_2018_low[df_2010_2018_low['year']==2010]))]
df_2010_2018_low.loc[df_2010_2018_low['year']==2011,'str_year']=["2011-2012" for x in range(len(df_2010_2018_low[df_2010_2018_low['year']==2011]))]
df_2010_2018_low.loc[df_2010_2018_low['year']==2012,'str_year']=["2011-2012" for x in range(len(df_2010_2018_low[df_2010_2018_low['year']==2012]))]
df_2010_2018_low.loc[df_2010_2018_low['year']==2013,'str_year']=["2013-2014" for x in range(len(df_2010_2018_low[df_2010_2018_low['year']==2013]))]
df_2010_2018_low.loc[df_2010_2018_low['year']==2014,'str_year']=["2013-2014" for x in range(len(df_2010_2018_low[df_2010_2018_low['year']==2014]))]
df_2010_2018_low.loc[df_2010_2018_low['year']==2017,'str_year']=["2017-2018" for x in range(len(df_2010_2018_low[df_2010_2018_low['year']==2017]))]
df_2010_2018_low.loc[df_2010_2018_low['year']==2018,'str_year']=["2017-2018" for x in range(len(df_2010_2018_low[df_2010_2018_low['year']==2018]))]

#Append all the dataframes together
df_all=df_2002_2003_green
df_all=df_all.append(df_2010_2018)

#Prepare plot
fig, axs = plt.subplots(2, 3)#, gridspec_kw={'width_ratios': [1, 3]})
fig.suptitle('Iceslabs inland progression')

axs = axs.ravel()
i=0

for region in df_all['key_shp'].unique():
    if (region == 'Out'):
        continue
    
    #Add 2010-2017 and 2017-2018!
    sns.violinplot(ax=axs[i], data=df_all[df_all['key_shp']==region], x="str_year", y="elevation",
               inner="quart", linewidth=1)
    sns.despine(left=True)

    #Set title
    axs[i].title.set_text(region)
    
    axs[i].grid()
    #Update count
    i=i+1
                        
    '''
    catplot is noce but did not managed to to subplots with it
    sns.catplot(data=df_2002_2003_green[df_2002_2003_green['key_shp']==region], kind="violin", x="year", y="elevation", hue="colorcode_icelens",ax = axs[i])
    '''
#######################################################################
###  Violin plot - Inland expansion of iceslabs from 2002 to 2018   ###
#######################################################################   

#######################################################################
###   Slice plot - Inland expansion of iceslabs from 2002 to 2018   ###
#######################################################################   
### ------------------------------ 2002-2003 ----------------------------- ###
#Create a dictionnary where to store slices information
dict_lat_slice_west={}

#Initialize the slice summary
slice_summary=np.zeros((len(lat_slices),5))*np.nan

#Initialize the slice summary lat
slice_lon_summary=np.zeros((len(lat_slices),5))*np.nan

#Loop over the traces
for trace in traces:
    
    #Check whether we are dealing with single or consecutive traces
    if(len(trace)>1):
        #We are dealing with consecutive traces
        #Select the data related to the first trace
        data_trace=df_2002_2003[df_2002_2003['Track_name']==trace[0]]
        
        #loop over the traces and append data to each other, do not take the first one
        for indiv_trace in list(trace[1:]):
            #Select all the data related to this trace
            data_trace=data_trace.append(df_2002_2003[df_2002_2003['Track_name']==indiv_trace])
            
    else:
        #We are dealing with individual traces
        #Select all the data related to this trace
        data_trace=df_2002_2003[df_2002_2003['Track_name']==trace[0]]

    #Now my data_trace datasets are ready to be worked with
    #Keep only green ice slabs
    data_trace=data_trace[data_trace['colorcode_icelens']==1]
    
    if (len(data_trace)<1):
        #No green ice slabs, continue
        continue
    else:
        #Check the flag: shall we store data? We are only interested in high end of ice slabs
        if trace[0] in list(flag_high):
            print('2002-2003 ice lens present for',trace[0])
            #Identify the max of that trace, and assign to the corresponding lat_slice
            max_to_store=np.max(data_trace['elevation'])
            
            #Index where maximum
            ind_max=np.where(data_trace['elevation']==np.max(data_trace['elevation']))
            
            #Identify corresponding lat and lon of maximum elevation
            lat_max=np.asarray(data_trace.iloc[ind_max]['lat_3413'])[0]
            lon_max=np.asarray(data_trace.iloc[ind_max]['lon_3413'])[0]
            
            #Check if a maximum have already been identified here. If yes, compare
            #the two. If latter > than former, store this new max. If not, continue
            if (np.isnan(slice_summary[i,0])):
                #No data for this slice yet, store the data
                #Identify to which slice it belongs to
                for i in range(1,len(lat_slices)):
                    if ((lat_max>=lat_slices[i-1]) and (lat_max<lat_slices[i])):
                        slice_summary[i,0]=max_to_store
                        slice_lon_summary[i,0]=lon_max
                    else:
                        continue
                        #store the coprresponding max in the corresponding slice
            else:
                print('Max already present')
                #Data for this slice alreadey present check
                if (max_to_store>slice_summary[i,0]):
                    print('Replace max by new max')
                    #Identify to which slice it belongs to
                    for i in range(1,len(lat_slices)):
                        print(lat_slices[i])
                        if ((lat_max>=lat_slices[i-1]) and (lat_max<lat_slices[i])):
                            slice_summary[i,0]=max_to_store
                            slice_lon_summary[i,0]=lon_max
                        else:
                            continue
                            #store the coprresponding max in the corresponding slice
                    
                

### ------------------------------ 2002-2003 ----------------------------- ###

### ------------------------------ 2010-2018 ----------------------------- ###
count_lat=0
#Loop over each boundary of lat slices and store dataset related to slices
for i in range(1,len(lat_slices)):
    
    #Identify low and higher end of the slice
    low_bound=lat_slices[i-1]
    high_bound=lat_slices[i]
    
    #Select all the data belonging to this lat slice
    ind_slice=np.logical_and(np.array(df_2010_2018['lat_3413']>=low_bound),np.array(df_2010_2018['lat_3413']<high_bound))
    df_slice=df_2010_2018[ind_slice]
    
    #Affine data by selecting only west greenland
    ind_slice=np.array(df_slice['lon_3413']<-50000)
    df_slice_latlon=df_slice[ind_slice]
    
    #Store the associated df
    dict_lat_slice_west[str(int(lat_slices[i-1]))+' to '+str(int(lat_slices[i]))]=df_slice_latlon
    
    #Loop over the different time periods (2010, 2011-2012, 2013-2014, 2017-2018)
    count_period=0
    
    for time_period in list(['2010','2011-2012','2013-2014','2017-2018']):
        if (time_period == '2010'):
            df_under_use=df_slice_latlon[df_slice_latlon['year']==2010]
            
            #If max in this slice of this time period is lower than max identified
            #in previous time period, store the max of previous time period
            if (np.max(df_under_use['elevation'])<=(np.nanmax(slice_summary[count_lat,:]))):
                #slice_summary[count_lat,1]=slice_summary[count_lat,0]
                #slice_lon_summary[count_lat,1]=slice_lon_summary[count_lat,0]
                
                slice_summary[count_lat,1]=np.nan
                slice_lon_summary[count_lat,1]=np.nan
            else:
                #store the new max elevation
                slice_summary[count_lat,1]=np.max(df_under_use['elevation'])
                
                if (len(df_under_use)>0):
                    #Data in this slice, can do the lon picking. Several point have the same elevation, take the eastern one (<=> the max)
                    slice_lon_summary[count_lat,1]=np.max(np.unique(df_under_use[df_under_use['elevation']==np.max(df_under_use['elevation'])]['lon_3413']))
                
        elif (time_period == '2011-2012'):
            df_under_use=df_slice_latlon[(df_slice_latlon['year']>=2011) & (df_slice_latlon['year']<=2012)]
            
            #If max in this slice of this time period is lower than max identified
            #in previous time period, store the max of previous time period
            if (np.max(df_under_use['elevation'])<=(np.nanmax(slice_summary[count_lat,:]))):
                #slice_summary[count_lat,2]=slice_summary[count_lat,1]
                #slice_lon_summary[count_lat,2]=slice_lon_summary[count_lat,1]
                
                slice_summary[count_lat,2]=np.nan
                slice_lon_summary[count_lat,2]=np.nan
            else:
                #store the new max elevation
                slice_summary[count_lat,2]=np.max(df_under_use['elevation'])
                
                if (len(df_under_use)>0):
                    #Data in this slice, can do the lon picking. Several point have the same elevation, take the eastern one (<=> the max)
                    slice_lon_summary[count_lat,2]=np.max(np.unique(df_under_use[df_under_use['elevation']==np.max(df_under_use['elevation'])]['lon_3413']))
                
        elif (time_period == '2013-2014'):
            df_under_use=df_slice_latlon[(df_slice_latlon['year']>=2013) & (df_slice_latlon['year']<=2014)]
            
            #If max in this slice of this time period is lower than max identified
            #in previous time period, store the max of previous time period
            if (np.max(df_under_use['elevation'])<=(np.nanmax(slice_summary[count_lat,:]))):
                #slice_summary[count_lat,3]=slice_summary[count_lat,2]
                #slice_lon_summary[count_lat,3]=slice_lon_summary[count_lat,2]
                slice_summary[count_lat,3]=np.nan
                slice_lon_summary[count_lat,3]=np.nan
            else:
                #store the new max elevation
                slice_summary[count_lat,3]=np.max(df_under_use['elevation'])
                
                if (len(df_under_use)>0):
                    #Data in this slice, can do the lon picking. Several point have the same elevation, take the eastern one (<=> the max)
                    slice_lon_summary[count_lat,3]=np.max(np.unique(df_under_use[df_under_use['elevation']==np.max(df_under_use['elevation'])]['lon_3413']))
                
        elif (time_period == '2017-2018'):
            df_under_use=df_slice_latlon[(df_slice_latlon['year']>=2017) & (df_slice_latlon['year']<=2018)]
            
            #If max in this slice of this time period is lower than max identified
            #in previous time period, store the max of previous time period
            if (np.max(df_under_use['elevation'])<=np.nanmax(slice_summary[count_lat,:])):
                #slice_summary[count_lat,4]=slice_summary[count_lat,3]
                #slice_lon_summary[count_lat,4]=slice_lon_summary[count_lat,3]
                slice_summary[count_lat,4]=np.nan
                slice_lon_summary[count_lat,4]=np.nan
            else:
                #store the new max elevation
                slice_summary[count_lat,4]=np.max(df_under_use['elevation'])
                
                if (len(df_under_use)>0):
                    #Data in this slice, can do the lon picking. Several point have the same elevation, take the eastern one (<=> the max)
                    slice_lon_summary[count_lat,4]=np.max(np.unique(df_under_use[df_under_use['elevation']==np.max(df_under_use['elevation'])]['lon_3413']))
               
        else:
            print('Time period not known, break')
            break
        
        #Update count
        count_period=count_period+1
    
    #Update count_lat
    count_lat=count_lat+1
### ------------------------------ 2010-2018 ----------------------------- ###

fig, (ax1,ax2) = plt.subplots(1,2)#, gridspec_kw={'width_ratios': [1, 3]})
fig.suptitle('Iceslabs inland progression')

ax1.set_title('Elevation')
ax1.plot(slice_summary[:,0],lat_slices,'.',label='2002-2003')
#ax1.plot(slice_summary[:,1],lat_slices,'.',label='2010')
ax1.plot(slice_summary[:,2],lat_slices,'.',label='2011-2012')
#ax1.plot(slice_summary[:,2],lat_slices,'.',label='2013-2014')
ax1.plot(slice_summary[:,4],lat_slices,'.',label='2017-2018')

ax2.set_title('Longitude')
ax2.plot(slice_lon_summary[:,0],lat_slices,'.',label='2002-2003')
#ax2.plot(slice_lon_summary[:,1],lat_slices,'.',label='2010')
ax2.plot(slice_lon_summary[:,2],lat_slices,'.',label='2011-2012')
#ax2.plot(slice_lon_summary[:,2],lat_slices,'.',label='2013-2014')
ax2.plot(slice_lon_summary[:,4],lat_slices,'.',label='2017-2018')


plt.legend()
plt.show()

#Prepare plot
fig, ax1 = plt.subplots()#, gridspec_kw={'width_ratios': [1, 3]})
fig.suptitle('Iceslabs difference 2017-2018 minus 2011-2012')

ax1.plot(slice_summary[:,4]-slice_summary[:,2],lat_slices,'.')

plt.show()

#######################################################################
###   Slice plot - Inland expansion of iceslabs from 2002 to 2018   ###
#######################################################################   

fig, ax1 = plt.subplots()
ax1.step(slice_summary[:,4]-slice_summary[:,0],lat_slices,label='2017-2018 minus 2002-2003',color='#238b45')
ax1.step(slice_summary[:,4]-slice_summary[:,1],lat_slices,label='2017-2018 minus 2010',color='#2b8cbe')
ax1.step(slice_summary[:,4]-slice_summary[:,2],lat_slices,label='2017-2018 minus 2011-2012',color='#c994c7')
ax1.step(slice_summary[:,4]-slice_summary[:,3],lat_slices,label='2017-2018 minus 2013-2014',color='#7a0177')
plt.show()

#Display Fig.1

path_flightlines='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/flightlines/'
flightlines_20022018=pd.read_csv(path_flightlines+'flightlines_20022018_GrIS.csv',decimal='.',sep=',')#,low_memory=False)
#flightlines_20022018_GrIS

#Transform the coordinates from WGS84 to EPSG:3413
transformer = Transformer.from_crs("EPSG:4326", "EPSG:3413", always_xy=True)
points=transformer.transform(np.asarray(flightlines_20022018["LON"]),np.asarray(flightlines_20022018["LAT"]))

#Store lat/lon in 3413
flightlines_20022018['lon_3413']=points[0]
flightlines_20022018['lat_3413']=points[1]

plot_fig1(df_all,flightlines_20022018,df_2010_2018_low,df_2010_2018_high,df_firn_aquifer_all,df_thickness_likelihood_20102018)

pdb.set_trace()


#Boxplot of max elevation per lat slice for each region for different time periods
#Loop over dict_lat_slices_summary

#######################################################################
###     Barplot of maximum elevation of ice slabs per time period   ###
#######################################################################
#Barplot per regions for each elevation slice
fig, (ax1) = plt.subplots()
ax1.bar(lat_slices,slice_lon_summary[:,4],width=10000,label='2017-2018')
ax1.bar(lat_slices,slice_lon_summary[:,3],width=10000,label='2013-2014')
ax1.bar(lat_slices,slice_lon_summary[:,2],width=10000,label='2012-2010')
ax1.bar(lat_slices,slice_lon_summary[:,1],width=10000,label='2010')
ax1.bar(lat_slices,slice_lon_summary[:,0],width=10000,label='2002-2003')
plt.legend()
plt.show()
#######################################################################
###     Barplot of maximum elevation of ice slabs per time period   ###
#######################################################################

fig, (ax1) = plt.subplots()
ax1.bar(lat_slices,slice_lon_summary[:,4],width=10000)
plt.show()
