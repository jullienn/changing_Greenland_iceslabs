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
                        ax1c.scatter(pnt_matched.lon_3413,pnt_matched.lat_3413,zorder=3,s=0.1)
                        plt.show()
                        pdb.set_trace()
                    
                    #Update area_region
                    area_region=area_region+concave_hull.area #I do not think this is correct. Look for 'area' on this webpage https://gist.github.com/dwyerk/10561690
            
            #Store total area per region and per time period
            summary_area[time_period][region]=area_region  
            
    return summary_area


def plot_pannels_supp(ax_plot,flightlines_20022018,df_firn_aquifer_all,df_all,time_period,label_panel):
    
    #Display GrIS drainage bassins
    NO_rignotetal.plot(ax=ax_plot,color='white', edgecolor='black',linewidth=0.5)
    NE_rignotetal.plot(ax=ax_plot,color='white', edgecolor='black',linewidth=0.5) 
    SE_rignotetal.plot(ax=ax_plot,color='white', edgecolor='black',linewidth=0.5) 
    SW_rignotetal.plot(ax=ax_plot,color='white', edgecolor='black',linewidth=0.5) 
    CW_rignotetal.plot(ax=ax_plot,color='white', edgecolor='black',linewidth=0.5) 
    NW_rignotetal.plot(ax=ax_plot,color='white', edgecolor='black',linewidth=0.5) 
    
    if (time_period=='2010'):
        #Issue, there are 2010 and '2010': take both
        #Display flightlines of this time period    
        ax_plot.scatter(flightlines_20022018[flightlines_20022018.str_year==2010]['lon_3413'],
                        flightlines_20022018[flightlines_20022018.str_year==2010]['lat_3413'],
                        s=1,marker='.',linewidths=0,c='#d9d9d9')
        
        #Display flightlines of this time period    
        ax_plot.scatter(flightlines_20022018[flightlines_20022018.str_year=='2010']['lon_3413'],
                        flightlines_20022018[flightlines_20022018.str_year=='2010']['lat_3413'],
                        s=1,marker='.',linewidths=0,c='#d9d9d9',label='Flightlines')
    else:
        #Display flightlines of this time period    
        ax_plot.scatter(flightlines_20022018[flightlines_20022018.str_year==time_period]['lon_3413'],
                        flightlines_20022018[flightlines_20022018.str_year==time_period]['lat_3413'],
                        s=1,marker='.',linewidths=0,c='#d9d9d9',label='Flightlines')
        
    '''
    #Display firn aquifers
    ax_plot.scatter(df_firn_aquifer_all['lon_3413'],df_firn_aquifer_all['lat_3413'],s=1,color='#238b45',label='Firn aquifers')
    '''
    
    if (time_period=='2002-2003'):
        #Display 2002-2003 iceslabs
        ax_plot.scatter(df_all[df_all.str_year=='2002-2003']['lon_3413'],
                        df_all[df_all.str_year=='2002-2003']['lat_3413'],
                        s=3,marker='.',color='#8c6bb1',linewidths=0,
                        label='2002-2003 ice slabs')
        color_legend_display='#8c6bb1'
    else:
        #Display iceslabs thickness of the corresponding time period
        lik_blues=ax_plot.scatter(df_all[df_all.str_year==time_period]['lon_3413'],
                                  df_all[df_all.str_year==time_period]['lat_3413'],
                                  c=df_all[df_all.str_year==time_period]['20m_ice_content_m'],
                                  s=3,marker='.',cmap=plt.get_cmap('Blues'),linewidths=0,
                                  label=time_period+' ice slabs')        
        color_legend_display='#4292c6'
    
    if (time_period=='2017-2018'):
        #Inspired from this https://matplotlib.org/stable/gallery/axes_grid1/demo_colorbar_with_inset_locator.html
        axins1 = inset_axes(ax_plot,
                            width="5%",  # width = 50% of parent_bbox width
                            height="100%",  # height : 5%
                            loc='lower left',
                            bbox_to_anchor=(1, 0., 1, 1),
                            bbox_transform=ax_plot.transAxes,
                            borderpad=0)
        
        cbar_blues=plt.colorbar(lik_blues, ax=ax_plot, cax=axins1, shrink=1,orientation='vertical')
        cbar_blues.set_label('Columnal ice content [m]')
    
    #Add label
    ax_plot.text(0, 1, label_panel,zorder=10, ha='center', va='center', transform=ax_plot.transAxes, weight='bold',fontsize=20)#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot

    '''
    #Legend display: this is from https://stackoverflow.com/questions/24706125/setting-a-fixed-size-for-points-in-legend
    # Create dummy Line2D objects for legend
    h1 = Line2D([0], [0], marker='o', markersize=np.sqrt(20), color='#d9d9d9', linestyle='None')
    h2 = Line2D([0], [0], marker='o', markersize=np.sqrt(20), color=color_legend_display, linestyle='None')
    
    # Plot legend.
    ax_plot.legend([h1, h2], ['Flightlines', time_period+' ice slabs'], loc="lower left", markerscale=2,
               scatterpoints=1, fontsize=10)
    '''
    
    ###################### From Tedstone et al., 2022 #####################
    #from plot_map_decadal_change.py
    # Define the CartoPy CRS object.
    crs = ccrs.NorthPolarStereo(central_longitude=-45., true_scale_latitude=70.)
    # This can be converted into a `proj4` string/dict compatible with GeoPandas
    crs_proj4 = crs.proj4_init
    # x0, x1, y0, y1
    ax_plot.set_extent([-692338, 916954, -3392187, -627732], crs=crs)
    gl=ax_plot.gridlines(draw_labels=True, xlocs=[-50, -35], ylocs=[65, 75], x_inline=False, y_inline=False, color='#969696',linewidth=0.5)
    
    #Customize lat labels
    gl.ylabels_right = False
            
    ax_plot.set_title(time_period,weight='bold',fontsize=20)
    ax_plot.axis('off')
    #scalebar.scale_bar(ax, (0, 0), 300, zorder=200)
    ###################### From Tedstone et al., 2022 #####################
    
    #Display region name on panel a 
    ax_plot.text(NO_rignotetal.centroid.x,NO_rignotetal.centroid.y+20000,np.asarray(NO_rignotetal.SUBREGION1)[0])
    ax_plot.text(NE_rignotetal.centroid.x,NE_rignotetal.centroid.y+20000,np.asarray(NE_rignotetal.SUBREGION1)[0])
    ax_plot.text(SE_rignotetal.centroid.x,SE_rignotetal.centroid.y+20000,np.asarray(SE_rignotetal.SUBREGION1)[0])
    ax_plot.text(SW_rignotetal.centroid.x,SW_rignotetal.centroid.y+20000,np.asarray(SW_rignotetal.SUBREGION1)[0])
    ax_plot.text(CW_rignotetal.centroid.x,CW_rignotetal.centroid.y+20000,np.asarray(CW_rignotetal.SUBREGION1)[0])
    ax_plot.text(NW_rignotetal.centroid.x,NW_rignotetal.centroid.y+20000,np.asarray(NW_rignotetal.SUBREGION1)[0])


def display_panels_c(ax1c,region_rignot,x0,x1,y0,y1,flightlines_20022018,df_thickness_likelihood_20102018,crs):
    
    ax1c.set_facecolor('white')

    #Display GrIS drainage bassins of the specific region
    region_rignot.plot(ax=ax1c,color='white', edgecolor='black',linewidth=0.5)
    
    #Flightlines
    # --- 2011-2012
    ax1c.scatter(flightlines_20022018[flightlines_20022018.str_year=='2011-2012']['lon_3413'],
                flightlines_20022018[flightlines_20022018.str_year=='2011-2012']['lat_3413'],
                s=0.1,marker='.',linewidths=0,c='#d9d9d9',label='flightlines 2011-2012')
    
    # --- 2017-2018
    ax1c.scatter(flightlines_20022018[flightlines_20022018.str_year=='2017-2018']['lon_3413'],
                flightlines_20022018[flightlines_20022018.str_year=='2017-2018']['lat_3413'],
                s=0.1,marker='.',linewidths=0,c='#969696',label='flightlines 2017-2018')
    
    #Likelihood
    # --- 2011-2012
    ax1c.scatter(df_thickness_likelihood_20102018[df_thickness_likelihood_20102018.Track_name.str[:4]=='2011']['lon_3413'],
                df_thickness_likelihood_20102018[df_thickness_likelihood_20102018.Track_name.str[:4]=='2011']['lat_3413'],
                c=df_thickness_likelihood_20102018[df_thickness_likelihood_20102018.Track_name.str[:4]=='2011']['likelihood'],
                s=15,marker='.',linewidths=0,cmap=plt.get_cmap('Blues'))
    
    lik_blues=ax1c.scatter(df_thickness_likelihood_20102018[df_thickness_likelihood_20102018.Track_name.str[:4]=='2012']['lon_3413'],
                df_thickness_likelihood_20102018[df_thickness_likelihood_20102018.Track_name.str[:4]=='2012']['lat_3413'],
                c=df_thickness_likelihood_20102018[df_thickness_likelihood_20102018.Track_name.str[:4]=='2012']['likelihood'],
                s=15,marker='.',linewidths=0,cmap=plt.get_cmap('Blues'))
    
    # --- 2017-2018            
    ax1c.scatter(df_thickness_likelihood_20102018[df_thickness_likelihood_20102018.Track_name.str[:4]=='2017']['lon_3413'],
                df_thickness_likelihood_20102018[df_thickness_likelihood_20102018.Track_name.str[:4]=='2017']['lat_3413'],
                c=df_thickness_likelihood_20102018[df_thickness_likelihood_20102018.Track_name.str[:4]=='2017']['likelihood'],
                s=3,marker='.',linewidths=0,cmap=plt.get_cmap('Reds'))
    lik_reds=ax1c.scatter(df_thickness_likelihood_20102018[df_thickness_likelihood_20102018.Track_name.str[:4]=='2018']['lon_3413'],
                df_thickness_likelihood_20102018[df_thickness_likelihood_20102018.Track_name.str[:4]=='2018']['lat_3413'],
                c=df_thickness_likelihood_20102018[df_thickness_likelihood_20102018.Track_name.str[:4]=='2018']['likelihood'],
                s=3,marker='.',linewidths=0,cmap=plt.get_cmap('Reds'))
    
    '''
    # Plot legend. This is from https://stackoverflow.com/questions/24706125/setting-a-fixed-size-for-points-in-legend
    lgnd = plt.legend(loc="best", scatterpoints=1)
    lgnd.legendHandles[0]._sizes = [30]
    lgnd.legendHandles[1]._sizes = [30]
    '''
    '''
    import matplotlib.patches as patches
    from matplotlib.patches import Patch
    
    #Custom legend myself
    legend_elements = [Patch(facecolor='#d9d9d9',label='flightlines 2011-2012'),
                       Patch(facecolor='#969696',label='flightlines 2017-2018'),
                       Patch(facecolor='#2171b5',label='Likelihood 2011-2012'),
                       Patch(facecolor='#cb181d',label='Likelihood 2017-2018')]

    ax1c.legend(handles=legend_elements)
    plt.legend()

    '''
    '''
    cbar_reds = plt.colorbar(lik_reds,location = 'right')
    cbar_reds.set_label('Columnal average likelihood - 2017-2018')
    cbar_blues = plt.colorbar(lik_blues,location = 'left')
    cbar_blues.set_label('Columnal average likelihood - 2011-2012')
    '''
    
    ###################### From Tedstone et al., 2022 #####################
    #from plot_map_decadal_change.py
    # x0, x1, y0, y1
    ax1c.set_extent([x0, x1, y0, y1], crs=crs)
    #ax1c.gridlines(draw_labels=True, xlocs=[-50, -35], ylocs=[65, 75], x_inline=False, y_inline=False,linewidth=0.5)
    #import scalebar
    #scalebar.scale_bar(ax1c, (0.65, 0.06), 300)
    ax1c.axis('off')
    ###################### From Tedstone et al., 2022 #####################
    return

def plot_fig1(df_all,flightlines_20022018,df_2010_2018_low,df_2010_2018_high,df_firn_aquifer_all,df_thickness_likelihood_20102018,dict_summary):   
    plot_fig_S1='FALSE'
    plot_panela='FALSE'
    plot_panelb='FALSE'
    plot_panelc='TRUE'
    
    if (plot_fig_S1 == 'TRUE'):
        # -------------------------------- FIG S1 --------------------------------
        ###################### From Tedstone et al., 2022 #####################
        #from plot_map_decadal_change.py
        # Define the CartoPy CRS object.
        crs = ccrs.NorthPolarStereo(central_longitude=-45., true_scale_latitude=70.)
        
        # This can be converted into a `proj4` string/dict compatible with GeoPandas
        crs_proj4 = crs.proj4_init
        ###################### From Tedstone et al., 2022 #####################
        
        fig = plt.figure(figsize=(14,50))
        gs = gridspec.GridSpec(7, 25)
        gs.update(wspace = 2.5)
        #gs.update(wspace=0.001)
        #projection set up from https://stackoverflow.com/questions/33942233/how-do-i-change-matplotlibs-subplot-projection-of-an-existing-axis
        ax1 = plt.subplot(gs[0:7, 0:5],projection=crs)
        ax2 = plt.subplot(gs[0:7, 5:10],projection=crs)
        ax3 = plt.subplot(gs[0:7, 10:15],projection=crs)
        ax4 = plt.subplot(gs[0:7, 15:20],projection=crs)
        ax5 = plt.subplot(gs[0:7, 20:25],projection=crs)
        
        ax1.set_facecolor('white')
        ax2.set_facecolor('white')
        ax3.set_facecolor('white')
        ax4.set_facecolor('white')
        ax5.set_facecolor('white')
                
        plot_pannels_supp(ax1,flightlines_20022018,df_firn_aquifer_all,df_all,'2002-2003','a')
        plot_pannels_supp(ax2,flightlines_20022018,df_firn_aquifer_all,df_all,'2010','b')
        plot_pannels_supp(ax3,flightlines_20022018,df_firn_aquifer_all,df_all,'2011-2012','c')
        plot_pannels_supp(ax4,flightlines_20022018,df_firn_aquifer_all,df_all,'2013-2014','d')
        plot_pannels_supp(ax5,flightlines_20022018,df_firn_aquifer_all,df_all,'2017-2018','e')
        
        figManager = plt.get_current_fig_manager()
        figManager.window.showMaximized()
        
        #Save the figure
        plt.savefig('C:/Users/jullienn/switchdrive/Private/research/RT1/figures/S1/v3/figS1.png',dpi=300)
        # -------------------------------- FIG S1 --------------------------------
    
    # --------------------------------- FIG 1 --------------------------------
    ###################### From Tedstone et al., 2022 #####################
    #from plot_map_decadal_change.py
    # Define the CartoPy CRS object.
    crs = ccrs.NorthPolarStereo(central_longitude=-45., true_scale_latitude=70.)
    # This can be converted into a `proj4` string/dict compatible with GeoPandas
    crs_proj4 = crs.proj4_init
    ###################### From Tedstone et al., 2022 #####################
        
    #Prepare Fig. 1
    fig = plt.figure(figsize=(14,50))
    gs = gridspec.GridSpec(20, 16)
    gs.update(wspace = 0.5)
    #gs.update(wspace=0.001)
    #projection set up from https://stackoverflow.com/questions/33942233/how-do-i-change-matplotlibs-subplot-projection-of-an-existing-axis
    axmap = plt.subplot(gs[0:20, 0:10],projection=crs)
    axNO = plt.subplot(gs[0:4, 10:13],projection=crs)
    axNE = plt.subplot(gs[0:4, 13:16],projection=crs)
    axNW = plt.subplot(gs[4:10, 10:12],projection=crs)
    axCW = plt.subplot(gs[4:10, 12:14],projection=crs)
    axSW = plt.subplot(gs[4:10, 14:16],projection=crs)
    axelev = plt.subplot(gs[10:20, 10:16])
    
    if (plot_panela=='TRUE'):
        # -------------------------------- PANEL A --------------------------------
        axmap.set_facecolor('white')

        #Display GrIS drainage bassins
        NO_rignotetal.plot(ax=axmap,color='white', edgecolor='black')
        NE_rignotetal.plot(ax=axmap,color='white', edgecolor='black') 
        SE_rignotetal.plot(ax=axmap,color='white', edgecolor='black') 
        SW_rignotetal.plot(ax=axmap,color='white', edgecolor='black') 
        CW_rignotetal.plot(ax=axmap,color='white', edgecolor='black') 
        NW_rignotetal.plot(ax=axmap,color='white', edgecolor='black')         

        #Display 2002-2018 flightlines
        axmap.scatter(flightlines_20022018['lon_3413'],flightlines_20022018['lat_3413'],s=0.1,marker='.',linewidths=0,color='#d9d9d9',label='Flightlines')#,label='2002-2003')
        
        #Display 2010-2014 iceslabs
        axmap.scatter(df_all[df_all.Track_name.str[:4]=='2010']['lon_3413'],df_all[df_all.Track_name.str[:4]=='2010']['lat_3413'],s=7,marker='.',linewidths=0,color='#4575b4',label='2010-2014 ice slabs')
        axmap.scatter(df_all[df_all.Track_name.str[:4]=='2011']['lon_3413'],df_all[df_all.Track_name.str[:4]=='2011']['lat_3413'],s=7,marker='.',linewidths=0,color='#4575b4')
        axmap.scatter(df_all[df_all.Track_name.str[:4]=='2012']['lon_3413'],df_all[df_all.Track_name.str[:4]=='2012']['lat_3413'],s=7,marker='.',linewidths=0,color='#4575b4')
        axmap.scatter(df_all[df_all.Track_name.str[:4]=='2013']['lon_3413'],df_all[df_all.Track_name.str[:4]=='2013']['lat_3413'],s=7,marker='.',linewidths=0,color='#4575b4')
        axmap.scatter(df_all[df_all.Track_name.str[:4]=='2014']['lon_3413'],df_all[df_all.Track_name.str[:4]=='2014']['lat_3413'],s=7,marker='.',linewidths=0,color='#4575b4')
        
        #Display 2017-2018 iceslabs
        axmap.scatter(df_all[df_all.Track_name.str[:4]=='2017']['lon_3413'],df_all[df_all.Track_name.str[:4]=='2017']['lat_3413'],s=3,marker='.',linewidths=0,color='#d73027',label='2017-2018 ice slabs')
        axmap.scatter(df_all[df_all.Track_name.str[:4]=='2018']['lon_3413'],df_all[df_all.Track_name.str[:4]=='2018']['lat_3413'],s=3,marker='.',linewidths=0,color='#d73027')
        
        #Display 2002-2003 iceslabs
        axmap.scatter(df_all[df_all.str_year=='2002-2003']['lon_3413'],df_all[df_all.str_year=='2002-2003']['lat_3413'],s=3,marker='.',linewidths=0,color='#2ECC71',label='2002-2003 ice slabs')
        
        #Display firn aquifers
        axmap.scatter(df_firn_aquifer_all['lon_3413'],df_firn_aquifer_all['lat_3413'],s=3,marker='.',linewidths=0,color='#807dba',label='Firn aquifers')
                
        #Display region name on panel a 
        axmap.text(NO_rignotetal.centroid.x,NO_rignotetal.centroid.y+20000,np.asarray(NO_rignotetal.SUBREGION1)[0])
        axmap.text(NE_rignotetal.centroid.x,NE_rignotetal.centroid.y+20000,np.asarray(NE_rignotetal.SUBREGION1)[0])
        axmap.text(SE_rignotetal.centroid.x,SE_rignotetal.centroid.y+20000,np.asarray(SE_rignotetal.SUBREGION1)[0])
        axmap.text(SW_rignotetal.centroid.x,SW_rignotetal.centroid.y+20000,np.asarray(SW_rignotetal.SUBREGION1)[0])
        axmap.text(CW_rignotetal.centroid.x,CW_rignotetal.centroid.y+20000,np.asarray(CW_rignotetal.SUBREGION1)[0])
        axmap.text(NW_rignotetal.centroid.x,NW_rignotetal.centroid.y+20000,np.asarray(NW_rignotetal.SUBREGION1)[0])
        
        # Plot legend. This is from https://stackoverflow.com/questions/24706125/setting-a-fixed-size-for-points-in-legend
        lgnd = axmap.legend(loc="lower right", scatterpoints=1)
        lgnd.legendHandles[0]._sizes = [30]
        lgnd.legendHandles[1]._sizes = [30]
        lgnd.legendHandles[2]._sizes = [30]
        lgnd.legendHandles[3]._sizes = [30]
        lgnd.legendHandles[4]._sizes = [30]
        
        pdb.set_trace()

        #Display coastlines
        axmap.coastlines(edgecolor='black',linewidth=0.75)
        
        ###################### From Tedstone et al., 2022 #####################
        #from plot_map_decadal_change.py
        # x0, x1, y0, y1
        axmap.set_extent([-692338, 916954, -3392187, -627732], crs=crs)
        axmap.gridlines(draw_labels=True, xlocs=[-50, -35], ylocs=[65, 75], x_inline=False, y_inline=False,color='#969696')
        #scalebar.scale_bar(ax, (0, 0), 300, zorder=200)
        axmap.axis('off')
        pdb.set_trace()
        ###################### From Tedstone et al., 2022 #####################        
        # -------------------------------- PANEL A --------------------------------

    if (plot_panelb=='TRUE'):
        
        # -------------------------------- PANEL B --------------------------------    
        #Define panel names
        labels = ['NO', 'NW', 'CW', 'SW']
    
        #Stack max_elev_mean data for barplot
        dplot_20022003=[dict_summary['NO']['2002-2003']['max_elev_median'],
                        dict_summary['NW']['2002-2003']['max_elev_median'],dict_summary['CW']['2002-2003']['max_elev_median'],
                        dict_summary['SW']['2002-2003']['max_elev_median']]
        
        dplot_2010=[dict_summary['NO']['2010']['max_elev_median'],
                    dict_summary['NW']['2010']['max_elev_median'],dict_summary['CW']['2010']['max_elev_median'],
                    dict_summary['SW']['2010']['max_elev_median']]
        
        dplot_20112012=[dict_summary['NO']['2011-2012']['max_elev_median'],
                        dict_summary['NW']['2011-2012']['max_elev_median'],dict_summary['CW']['2011-2012']['max_elev_median'],
                        dict_summary['SW']['2011-2012']['max_elev_median']]
        
        dplot_20132014=[dict_summary['NO']['2013-2014']['max_elev_median'],
                        dict_summary['NW']['2013-2014']['max_elev_median'],dict_summary['CW']['2013-2014']['max_elev_median'],
                        dict_summary['SW']['2013-2014']['max_elev_median']]
        
        dplot_20172018=[dict_summary['NO']['2017-2018']['max_elev_median'],
                        dict_summary['NW']['2017-2018']['max_elev_median'],dict_summary['CW']['2017-2018']['max_elev_median'],
                        dict_summary['SW']['2017-2018']['max_elev_median']]
        
        #Stack max_elev_mean data for barplot
        dplotstd_20022003=[dict_summary['NO']['2002-2003']['max_elev_std'],
                        dict_summary['NW']['2002-2003']['max_elev_std'],dict_summary['CW']['2002-2003']['max_elev_std'],
                        dict_summary['SW']['2002-2003']['max_elev_std']]
        
        dplotstd_2010=[dict_summary['NO']['2010']['max_elev_std'],
                    dict_summary['NW']['2010']['max_elev_std'],dict_summary['CW']['2010']['max_elev_std'],
                    dict_summary['SW']['2010']['max_elev_std']]
        
        dplotstd_20112012=[dict_summary['NO']['2011-2012']['max_elev_std'],
                        dict_summary['NW']['2011-2012']['max_elev_std'],dict_summary['CW']['2011-2012']['max_elev_std'],
                        dict_summary['SW']['2011-2012']['max_elev_std']]
        
        dplotstd_20132014=[dict_summary['NO']['2013-2014']['max_elev_std'],
                        dict_summary['NW']['2013-2014']['max_elev_std'],dict_summary['CW']['2013-2014']['max_elev_std'],
                        dict_summary['SW']['2013-2014']['max_elev_std']]
        
        dplotstd_20172018=[dict_summary['NO']['2017-2018']['max_elev_std'],
                        dict_summary['NW']['2017-2018']['max_elev_std'],dict_summary['CW']['2017-2018']['max_elev_std'],
                        dict_summary['SW']['2017-2018']['max_elev_std']]
        
        
        #Stack data for maximum elevation difference calculation
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
        N=4 #Number of regions
        ind= np.arange(N) #Position of regions
                
        axelev.bar(ind, dplot_20022003, width, label='2002-2003',color='#2ECC71', yerr= dplotstd_20022003) #yerr=men_std
        axelev.bar(ind+1*width, dplot_2010, width, label='2010',color='#9ecae1', yerr= dplotstd_2010)
        axelev.bar(ind+2*width, dplot_20112012, width, label='2011-2012',color='#6baed6', yerr= dplotstd_20112012)
        axelev.bar(ind+3*width, dplot_20132014, width, label='2013-2014',color='#3182bd', yerr= dplotstd_20132014)
        axelev.bar(ind+4*width, dplot_20172018, width, label='2017-2018',color='#d73027', yerr= dplotstd_20172018)
        axelev.set_xticks(ind + 2*width)
        axelev.set_xticklabels(labels)
        axelev.set_ylim(1000,2050)
        
        axelev.text(ind[0],np.nanmax(max_elev_diff_NO)+55,str(int(np.round(np.nanmax(max_elev_diff_NO)-np.nanmin(max_elev_diff_NO))))+' m')
        #axelev.text(ind[1],np.nanmax(max_elev_diff_NW)+180,str(int(np.round(np.nanmax(max_elev_diff_NW)-np.nanmin(max_elev_diff_NW))))+' m')
        axelev.text(ind[2],np.nanmax(max_elev_diff_CW)+30,str(int(np.round(np.nanmax(max_elev_diff_CW)-np.nanmin(max_elev_diff_CW))))+' m')
        axelev.text(ind[3],np.nanmax(max_elev_diff_SW)+90,str(int(np.round(np.nanmax(max_elev_diff_SW)-np.nanmin(max_elev_diff_SW))))+' m')
        
        axelev.set_ylabel('Elevation [m]')
        
        #Custom legend myself
        from matplotlib.patches import Patch
        from matplotlib.lines import Line2D
        
        legend_elements = [Patch(facecolor='#2ECC71',label='2002-2003'),
                           Patch(facecolor='#9ecae1',label='2010'),
                           Patch(facecolor='#6baed6',label='2011-2012'),
                           Patch(facecolor='#3182bd',label='2013-2014'),
                           Patch(facecolor='#d73027',label='2017-2018')]#,
                           #Line2D([0], [0], color='k', lw=2, label='Standard deviation around the mean')]
        axelev.legend(handles=legend_elements,loc='upper left')
        plt.legend()
        
        '''
        #Calculate ice slabs inland expansion rate
        table_for_rate_max=np.zeros((4,5))
        table_for_rate_max[:]=np.nan
        
        table_for_rate_median=np.zeros((4,5))
        table_for_rate_median[:]=np.nan
        
        i_rate=0
        for rate_region in list(['NO','NW','CW','SW']):
            j_rate=0
            for rate_year in list(['2002-2003','2010','2011-2012','2013-2014','2017-2018']):
                table_for_rate_max[i_rate,j_rate]=dict_summary[rate_region][rate_year]['max_elev_max']
                table_for_rate_median[i_rate,j_rate]=dict_summary[rate_region][rate_year]['max_elev_median']
                j_rate=j_rate+1
            i_rate=i_rate+1
                
        f, (ax1,ax2,ax3,ax4) = plt.subplots(1,4)
        ax1.plot(np.array([2003,2010,2012,2014,2018]), table_for_rate_max[0,:])
        ax1.scatter(np.array([2003,2010,2012,2014,2018]), table_for_rate_max[0,:])
        ax2.plot(np.array([2003,2010,2012,2014,2018]), table_for_rate_max[1,:])
        ax2.scatter(np.array([2003,2010,2012,2014,2018]), table_for_rate_max[1,:])
        ax3.plot(np.array([2003,2010,2012,2014,2018]), table_for_rate_max[2,:])
        ax3.scatter(np.array([2003,2010,2012,2014,2018]), table_for_rate_max[2,:])
        ax4.plot(np.array([2003,2010,2012,2014,2018]), table_for_rate_max[3,:])
        ax4.scatter(np.array([2003,2010,2012,2014,2018]), table_for_rate_max[3,:])
        plt.show()
        
        f, (ax1,ax2,ax3,ax4) = plt.subplots(1,4)
        ax1.plot(np.array([2003,2010,2012,2014,2018]), table_for_rate_median[0,:])
        ax1.scatter(np.array([2003,2010,2012,2014,2018]), table_for_rate_median[0,:])
        ax2.plot(np.array([2003,2010,2012,2014,2018]), table_for_rate_median[1,:])
        ax2.scatter(np.array([2003,2010,2012,2014,2018]), table_for_rate_median[1,:])
        ax3.plot(np.array([2003,2010,2012,2014,2018]), table_for_rate_median[2,:])
        ax3.scatter(np.array([2003,2010,2012,2014,2018]), table_for_rate_median[2,:])
        ax4.plot(np.array([2003,2010,2012,2014,2018]), table_for_rate_median[3,:])
        ax4.scatter(np.array([2003,2010,2012,2014,2018]), table_for_rate_median[3,:])
        plt.show()
        
        #rates:
        #NO
        print('Rate max NO:',str(np.round((table_for_rate_max[0,-1]-table_for_rate_max[0,0])/(2018-2003),2)),'m/yr')
        #NW
        print('Rate max NW:',str(np.round((table_for_rate_max[1,-1]-table_for_rate_max[1,1])/(2018-2010),2)),'m/yr')
        #CW
        print('Rate max CW:',str(np.round((table_for_rate_max[2,-1]-table_for_rate_max[2,0])/(2018-2003),2)),'m/yr')
        #SW
        print('Rate max SW:',str(np.round((table_for_rate_max[3,-1]-table_for_rate_max[3,0])/(2018-2003),2)),'m/yr')

        '''
        pdb.set_trace()

        # -------------------------------- PANEL B --------------------------------    
    
    if (plot_panelc=='TRUE'):
        
        hull_computation='TRUE'
        likelihood_display='TRUE'
        
        # -------------------------------- PANELS C -------------------------------        
        if (hull_computation=='TRUE'):
            #Panel C
            #Load convex hull mask over which convex hull must be computed
            path_convexhull_masks='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2002_2018/shapefiles/'
            
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
            
            #plt.scatter(df_2010_2018_high[df_2010_2018_high.Track_name.str[:4]=='2010']['lon_3413'],df_2010_2018_high[df_2010_2018_high.Track_name.str[:4]=='2010']['lat_3413'],s=1,color='#3690c0',label='2010-2014 ice slabs')
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
            do_plot='TRUE'
            high_end_summary=concave_hull_computation(df_2010_2018_high,dictionnaries_convexhullmasks,ax1c,do_plot,'high_end')
            do_plot='TRUE'
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
            #Save figure according to different regions
            for region in list(['NE','NO','NW','CW','SW']):
                if (region =='NE'):
                    x0=283000
                    x1=670000
                    y0=-1300000
                    y1=-940000
                    #display data
                    display_panels_c(axNE,NE_rignotetal,x0,x1,y0,y1,flightlines_20022018,df_thickness_likelihood_20102018,crs)
                    #Display region name
                    axNE.text(330000,-1280000,'NE',fontsize=15)
                    '''
                    ################ DISPLAY AREA CHANGE #####################
                    low_end_change=(int((low_end_summary['2017-2018'][region]-low_end_summary['2011-2012'][region])/low_end_summary['2011-2012'][region]*100))
                    high_end_change=(int((high_end_summary['2017-2018'][region]-high_end_summary['2011-2012'][region])/high_end_summary['2011-2012'][region]*100))
                    
                    #Compute and display relative change
                    axNE.text(300000,-1650000,'[+'+str(low_end_change)+' : +'+str(high_end_change)+'] %')
                    ################ DISPLAY AREA CHANGE #####################
                    '''
                elif(region =='NO'):
                    x0=-605000
                    x1=302000
                    y0=-1215000
                    y1=-785000
                    #display data
                    display_panels_c(axNO,NO_rignotetal,x0,x1,y0,y1,flightlines_20022018,df_thickness_likelihood_20102018,crs)
                    #Display region name
                    axNO.text(-90000,-1170000,'NO',fontsize=15)
                    '''
                    ################ DISPLAY AREA CHANGE #####################
                    low_end_change=(int((low_end_summary['2017-2018'][region]-low_end_summary['2011-2012'][region])/low_end_summary['2011-2012'][region]*100))
                    high_end_change=(int((high_end_summary['2017-2018'][region]-high_end_summary['2011-2012'][region])/high_end_summary['2011-2012'][region]*100))
                    
                    #Compute and display relative change
                    axNO.text(-90000,-1150000,'[+'+str(low_end_change)+' : +'+str(high_end_change)+'] %')
                    ################ DISPLAY AREA CHANGE #####################
                    '''
                elif(region =='NW'):
                    x0=-610000
                    x1=-189000
                    y0=-1140000
                    y1=-1985000
                    #display data
                    display_panels_c(axNW,NW_rignotetal,x0,x1,y0,y1,flightlines_20022018,df_thickness_likelihood_20102018,crs)
                    #Display region name
                    axNW.text(-395000,-1400000,'NW',fontsize=15)
                    
                    ################ DISPLAY AREA CHANGE #####################
                    low_end_change=(int((low_end_summary['2017-2018'][region]-low_end_summary['2011-2012'][region])/low_end_summary['2011-2012'][region]*100))
                    high_end_change=(int((high_end_summary['2017-2018'][region]-high_end_summary['2011-2012'][region])/high_end_summary['2011-2012'][region]*100))
                    
                    #Compute and display relative change
                    axNW.text(-400000,-1480000,'[+'+str(low_end_change)+':+'+str(high_end_change)+']%')
                    ################ DISPLAY AREA CHANGE #####################
                    
                elif(region =='CW'):
                    x0=-259000
                    x1=-60500
                    y0=-2385000
                    y1=-1935000
                    #display data
                    display_panels_c(axCW,CW_rignotetal,x0,x1,y0,y1,flightlines_20022018,df_thickness_likelihood_20102018,crs)
                    #Display region name
                    axCW.text(-207500,-2327000,'CW',fontsize=15)
                    
                    ################ DISPLAY AREA CHANGE #####################
                    low_end_change=(int((low_end_summary['2017-2018'][region]-low_end_summary['2011-2012'][region])/low_end_summary['2011-2012'][region]*100))
                    high_end_change=(int((high_end_summary['2017-2018'][region]-high_end_summary['2011-2012'][region])/high_end_summary['2011-2012'][region]*100))
                    
                    #Compute and display relative change
                    axCW.text(-285000,-2360000,'[+'+str(low_end_change)+':+'+str(high_end_change)+']%')
                    ################ DISPLAY AREA CHANGE #####################
                    
                elif(region =='SW'):
                    x0=-265000
                    x1=-55600
                    y0=-2899000
                    y1=-2370000
                    #display data
                    display_panels_c(axSW,SW_rignotetal,x0,x1,y0,y1,flightlines_20022018,df_thickness_likelihood_20102018,crs)
                    #Display region name
                    axSW.text(-167000,-2800000,'SW',fontsize=15)
                    
                    ################ DISPLAY AREA CHANGE #####################
                    low_end_change=(int((low_end_summary['2017-2018'][region]-low_end_summary['2011-2012'][region])/low_end_summary['2011-2012'][region]*100))
                    high_end_change=(int((high_end_summary['2017-2018'][region]-high_end_summary['2011-2012'][region])/high_end_summary['2011-2012'][region]*100))
                    
                    #Compute and display relative change
                    axSW.text(-200000,-2840000,'[+'+str(low_end_change)+':+'+str(high_end_change)+']%')
                    ################ DISPLAY AREA CHANGE #####################
                    
                else:
                    print('Region not known')
            
        # -------------------------------- PANELS C -------------------------------        
    #Force legend of pannel b to be upper left
    axelev.legend(handles=legend_elements,loc='upper left')
    
    #Display panels label
    axmap.text(0, 0.85,'a', ha='center', va='center', transform=axmap.transAxes,fontsize=25)#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot
    axNO.text(0, 0.45,'b', ha='center', va='center', transform=axNO.transAxes,fontsize=25)#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot
    axelev.text(0, 1.05,'c', ha='center', va='center', transform=axelev.transAxes,fontsize=25)#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot
    
    #Maximize plot size
    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()
    
    pdb.set_trace()
    
    #Save figure
    plt.savefig('C:/Users/jullienn/switchdrive/Private/research/RT1/figures/fig1/v4/fig1c_last.png',dpi=1000)

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
import cartopy.crs as ccrs
from matplotlib.lines import Line2D

import seaborn as sns
sns.set_theme(style="whitegrid")
from scalebar import scale_bar
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


#Set fontsize plot
plt.rcParams.update({'font.size': 10})


### -------------------------- Load GrIS DEM ----------------------------- ###
#https://towardsdatascience.com/reading-and-visualizing-geotiff-images-with-python-8dcca7a74510
import rasterio
from rasterio.plot import show

path_GrIS_DEM = r'C:/Users/jullienn/switchdrive/Private/research/backup_Aglaja/working_environment/greenland_topo_data/elevations/greenland_dem_mosaic_100m_v3.0.tif'
GrIS_DEM = rasterio.open(path_GrIS_DEM)

#fig, axs = plt.subplots()
#show(img,ax=axs)
### -------------------------- Load GrIS DEM ----------------------------- ###

### -------------------------- Load shapefiles --------------------------- ###
#from https://gis.stackexchange.com/questions/113799/how-to-read-a-shapefile-in-python
path_IceBridgeArea_Shape='C:/Users/jullienn/switchdrive/Private/research/backup_Aglaja/working_environment/greenland_topo_data/IceBridge Area Shapefiles/IceBridge Area Shapefiles/'
IceBridgeArea_Shape=gpd.read_file(path_IceBridgeArea_Shape+'IceBridgeArea_Shape.shp')

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

#Load Rignot et al., 2016 GrIS mask
path_rignotetal2016_GrIS='C:/Users/jullienn/switchdrive/Private/research/backup_Aglaja/working_environment/greenland_topo_data/GRE_IceSheet_IMBIE2/GRE_IceSheet_IMBIE2/'
GrIS_rignotetal2016=gpd.read_file(path_rignotetal2016_GrIS+'GRE_IceSheet_IMBIE2_v1_EPSG3413.shp',rows=slice(1,2,1)) #the regions are the last rows of the shapefile
GrIS_mask=GrIS_rignotetal2016[GrIS_rignotetal2016.SUBREGION1=='ICE_SHEET']
### -------------------------- Load shapefiles --------------------------- ###

### ---------------------------- Load dataset ---------------------------- ###
#Dictionnaries have already been created, load them
path_df_with_elevation='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2002_2018/final_excel/' 
#Load 2002-2003
f_20022003 = open(path_df_with_elevation+'df_2002_2003_with_elevation_rignotetalregions', "rb")
df_2002_2003 = pickle.load(f_20022003)
f_20022003.close()

#Load 2010-2018 high estimate
f_20102018_high = open(path_df_with_elevation+'high_estimate/df_20102018_with_elevation_high_estimate_rignotetalregions', "rb")
df_2010_2018_high = pickle.load(f_20102018_high)
f_20102018_high.close()

#Load 2010-2018 low estimate
f_20102018_low = open(path_df_with_elevation+'low_estimate/df_20102018_with_elevation_low_estimate_rignotetalregions', "rb")
df_2010_2018_low = pickle.load(f_20102018_low)
f_20102018_low.close()
### ---------------------------- Load dataset ---------------------------- ###

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
path_thickness_likelihood='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2002_2018/final_excel/high_estimate/'
df_thickness_likelihood_20102018 = pd.read_csv(path_thickness_likelihood+'Ice_Layer_Output_Thicknesses_2010_2018_jullienetal2021_high_estimate.csv',delimiter=',',decimal='.')
#Transform miege coordinates from WGS84 to EPSG:3413
points=transformer.transform(np.asarray(df_thickness_likelihood_20102018["lon"]),np.asarray(df_thickness_likelihood_20102018["lat"]))

#Store lat/lon in 3413
df_thickness_likelihood_20102018['lon_3413']=points[0]
df_thickness_likelihood_20102018['lat_3413']=points[1]

#######################################################################
###          Inland expansion of iceslabs from 2002 to 2018         ###
#######################################################################

#Select 2002-2003 with green ice slabs only
df_2002_2003_green=df_2002_2003[df_2002_2003['colorcode_icelens']==1]

#Set the year for plotting
df_2002_2003_green['str_year']=["2002-2003" for x in range(len(df_2002_2003_green))]

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

#Display Fig.1

path_flightlines='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/flightlines/'
flightlines_20022018_load=pd.read_csv(path_flightlines+'flightlines_20022018_GrIS.csv',decimal='.',sep=',')#,low_memory=False)

#Differentiate 2002-2003 VS 2010-2018
flightlines_20022003=flightlines_20022018_load[flightlines_20022018_load.str_year=='2002-2003']
flightlines_20102018=flightlines_20022018_load[flightlines_20022018_load.str_year!='2002-2003']

#Transform the coordinates from WGS84 to EPSG:3413
transformer = Transformer.from_crs("EPSG:4326", "EPSG:3413", always_xy=True)
points=transformer.transform(np.asarray(flightlines_20102018["LON"]),np.asarray(flightlines_20102018["LAT"]))

#Store lat/lon in 3413
flightlines_20102018['lon_3413']=points[0]
flightlines_20102018['lat_3413']=points[1]

#Aggregate 2002-2003 with 2010-2018
flightlines_20022018=flightlines_20022003
flightlines_20022018=flightlines_20022018.append(flightlines_20102018)


###############################################################################
###                  Aggregate data for pannel b display                    ###
###############################################################################
#IV. From here on, work with the different periods separated by strong melting summers.
#    Work thus with 2002-2003 VS 2010 VS 2011-2012 VS 2013-2014 VS 2017-2018
#    Select the absolute low and absolute high of 2002-2003, 2010-2014 and 2017-2018

#Let's create ~10km latitudinal (resp. longitudinal) slices for SW Greenland (resp. NW Greenland)
#and calculate the low and high end in each slice for elevation difference:

#1. Create the latitudinal (resp. longitudinal) slices
############ Is this grid to change?? This is about 10km width but not even whether north or south!
#lat_slices=np.arange(-2808967.3912300025,-1150000,10054.34783)#do not consider southern data
lat_slices=np.linspace(-3000000,-1150000,int((np.abs(-3000000)-np.abs(-1150000))/10000))#original
lon_slices=np.linspace(-600000,650000,int((np.abs(650000)+np.abs(-600000))/10000))#original

#2. Select and store all the data belonging to the lon/lat slices in a dictionnary.
#### ------------------------- 2010-2018 -------------------------------- ####
#   Retreive and store min and max elevation of each slice in a dataframe
#   ----- Latitudinal slices

#Create a dictionnary where to store slices information
dict_lat_slice={}

#Create a dictionnary to store np arrays storing slices min and max elevation for each region
dict_lat_slices_summary={k: {} for k in list(df_2010_2018_high['key_shp'].unique())}

#loop over the regions, create the room for each time period in each region
for region in list(df_2010_2018_high['key_shp'].unique()):
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
    ind_slice=np.logical_and(np.array(df_2010_2018_high['lat_3413']>=low_bound),np.array(df_2010_2018_high['lat_3413']<high_bound))
    df_slice=df_2010_2018_high[ind_slice]
    
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
dict_lon_slices_summary={k: {} for k in list(df_2010_2018_high['key_shp'].unique())}

#loop over the regions, create the room for each time period in each region
for region in list(df_2010_2018_high['key_shp'].unique()):
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
    ind_slice=np.logical_and(np.array(df_2010_2018_high['lon_3413']>=low_bound),np.array(df_2010_2018_high['lon_3413']<high_bound))
    df_slice=df_2010_2018_high[ind_slice]
    
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
    elif (trace[0] == 'may09_03_30_aggregated'):
        #if the firn aquier region, continue
        print('trace is ',str(trace[0]))
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
dict_summary={k: {} for k in list(df_2010_2018_high['key_shp'].unique())}

#Loop over the regions
for region in list(df_2010_2018_high['key_shp'].unique()):
    
    #Continue building the dictionnary
    dict_summary[region]={k: {} for k in list(['2002-2003','2010','2011-2012','2013-2014','2017-2018'])}
    
    #Loop over the 5 time periods
    
    for time_period in list(['2002-2003','2010','2011-2012','2013-2014','2017-2018']):
        dict_summary[region][time_period]={k: {} for k in list(['max_elev_mean','max_elev_median','max_elev_std','max_elev_max'])}
        
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
        dict_summary[region][time_period]['max_elev_max']=np.nanmax(dict_temp[:,1])


###############################################################################
###                  Aggregate data for pannel b display                    ###
###############################################################################

#Append all the dataframes together
df_all=df_2002_2003_green
df_all=df_all.append(df_2010_2018_high)

######################### Keep only data on the GrIS ##########################
# This is from aggregate_20022018_flightlines.py
df_firn_aquifer_all['coords'] = list(zip(df_firn_aquifer_all['lon_3413'],df_firn_aquifer_all['lat_3413']))
df_firn_aquifer_all['coords'] = df_firn_aquifer_all['coords'].apply(Point)
points = gpd.GeoDataFrame(df_firn_aquifer_all, geometry='coords', crs="EPSG:3413")
pointInPolys = gpd.tools.sjoin(points, GrIS_mask, op="within", how='left') #This is from https://www.matecdev.com/posts/point-in-polygon.html
df_firn_aquifer_all_GrIS = points[pointInPolys.SUBREGION1=='ICE_SHEET']
######################### Keep only data on the GrIS ##########################

######################### Keep only data on the GrIS ##########################
# This is from aggregate_20022018_flightlines.py
df_thickness_likelihood_20102018['coords'] = list(zip(df_thickness_likelihood_20102018['lon_3413'],df_thickness_likelihood_20102018['lat_3413']))
df_thickness_likelihood_20102018['coords'] = df_thickness_likelihood_20102018['coords'].apply(Point)
points = gpd.GeoDataFrame(df_thickness_likelihood_20102018, geometry='coords', crs="EPSG:3413")
pointInPolys = gpd.tools.sjoin(points, GrIS_mask, op="within", how='left') #This is from https://www.matecdev.com/posts/point-in-polygon.html
df_thickness_likelihood_20102018_all_GrIS = points[pointInPolys.SUBREGION1=='ICE_SHEET']
######################### Keep only data on the GrIS ##########################

######################### Keep only data on the GrIS ##########################
# This is from aggregate_20022018_flightlines.py
df_all['coords'] = list(zip(df_all['lon_3413'],df_all['lat_3413']))
df_all['coords'] = df_all['coords'].apply(Point)
points = gpd.GeoDataFrame(df_all, geometry='coords', crs="EPSG:3413")
pointInPolys = gpd.tools.sjoin(points, GrIS_mask, op="within", how='left') #This is from https://www.matecdev.com/posts/point-in-polygon.html
df_all_GrIS = points[pointInPolys.SUBREGION1=='ICE_SHEET']
######################### Keep only data on the GrIS ########################## 

'''
df_all_GrIS = df_all
df_firn_aquifer_all_GrIS=df_firn_aquifer_all
df_thickness_likelihood_20102018_all_GrIS=df_thickness_likelihood_20102018
'''

#Create Fig.1
plot_fig1(df_all_GrIS,flightlines_20022018,df_2010_2018_low,df_2010_2018_high,df_firn_aquifer_all_GrIS,df_thickness_likelihood_20102018_all_GrIS,dict_summary)

pdb.set_trace()


#Save 2002-2003 green dataset
#Open filename (same procedure as MacFerrin et al., 2019)
filename_excel_output='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2002_2018/final_excel/2002_2003_green_excel.csv'
fout = open(filename_excel_output, 'w')
header = "lat,lon,tracenum,key_shp,elevation,year\n"
fout.write(header)

tracecount = 0
for lat, lon, tracenum, key, elev, year in zip(np.asarray(df_2002_2003_green['lat_3413']),np.asarray(df_2002_2003_green['lon_3413']),np.asarray(df_2002_2003_green['Track_name']),np.asarray(df_2002_2003_green['key_shp']),np.asarray(df_2002_2003_green['elevation']),np.asarray(df_2002_2003_green['year'])):

    line = "{0},{1},{2},{3},{4},{5}\n".format(lat, lon, tracenum, key, elev, year)
    fout.write(line)
    tracecount += 1
print()

fout.close()


'''
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
#######################################################################
###   Slice plot - Inland expansion of iceslabs from 2002 to 2018   ###
#######################################################################   
'''

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