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


def draw_map(ax_plot,panel_label):
    
    ###################### From Tedstone et al., 2022 #####################
    #from plot_map_decadal_change.py
    # Define the CartoPy CRS object.
    int_crs = ccrs.NorthPolarStereo(central_longitude=-45., true_scale_latitude=70.)
    ###################### From Tedstone et al., 2022 #####################
    
    #Draw plot of GrIS map
    ax_plot.coastlines(edgecolor='black',linewidth=0.075)
    #Display GrIS drainage bassins limits
    GrIS_drainage_bassins.plot(ax=ax_plot,color='none', edgecolor='black',linewidth=0.075)
    #Display region name
    ax_plot.text(NO_rignotetal.centroid.x-125000,NO_rignotetal.centroid.y-150000,np.asarray(NO_rignotetal.SUBREGION1)[0])
    ax_plot.text(NE_rignotetal.centroid.x-150000,NE_rignotetal.centroid.y-100000,np.asarray(NE_rignotetal.SUBREGION1)[0])
    ax_plot.text(SE_rignotetal.centroid.x-100000,SE_rignotetal.centroid.y,np.asarray(SE_rignotetal.SUBREGION1)[0])
    ax_plot.text(SW_rignotetal.centroid.x-150000,SW_rignotetal.centroid.y-170000,np.asarray(SW_rignotetal.SUBREGION1)[0])
    ax_plot.text(CW_rignotetal.centroid.x-100000,CW_rignotetal.centroid.y-60000,np.asarray(CW_rignotetal.SUBREGION1)[0])
    ax_plot.text(NW_rignotetal.centroid.x-25000,NW_rignotetal.centroid.y-150000,np.asarray(NW_rignotetal.SUBREGION1)[0])
    
    # x0, x1, y0, y1
    ax_plot.set_extent([-692338, 916954, -3392187, -627732], crs=int_crs)

    ###################### From Tedstone et al., 2022 #####################
    #from plot_map_decadal_change.py
    gl=ax_plot.gridlines(draw_labels=True, xlocs=[-35, -50], ylocs=[65,75], x_inline=False, y_inline=False,linewidth=0.5,linestyle='dashed')
    #Customize lat labels
    gl.ylabels_right = False
    ax_plot.axis('off')
    
    if (panel_label in list(['a','b','c','d','e'])):
        gl.xlabels_bottom = False
        gl.xlabels_top = False
    else:
        gl.xlabels_top = False
    
    ###################### From Tedstone et al., 2022 #####################
    
    #Add panel labels
    ax_plot.text(-0.1, 0.95, panel_label,zorder=10, ha='center', va='center', transform=ax_plot.transAxes, weight='bold',fontsize=25)#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot
    
    if (panel_label == 'a'):
        #Display scalebar
        scale_bar(ax_plot, (0.745, 0.125), 200, 3,5)# axis, location (x,y), length, linewidth, rotation of text
        #by measuring on the screen, the difference in precision between scalebar and length of transects is about ~200m
        
    return

def display_flightlines(ax_plot,flightlines_to_display,timing):
    #Display flightlines of this time period    
    ax_plot.scatter(flightlines_to_display['lon_3413'],
                    flightlines_to_display['lat_3413'],
                    s=1,marker='.',linewidths=0,c='#d9d9d9',label='Flightlines')
    
    ax_plot.set_title(timing,fontsize=20, weight='bold')
    return

def display_iceslabs(ax_plot,iceslabs_20022018,timing):
    
    #Display 2002-2003 iceslabs
    ax_plot.scatter(iceslabs_20022018[iceslabs_20022018.str_year=='2002-2003']['lon_3413'],
                    iceslabs_20022018[iceslabs_20022018.str_year=='2002-2003']['lat_3413'],
                    s=10,marker='.',color='#8c6bb1',linewidths=0,
                    label='2002-2003 ice slabs')

    if (timing!='2002-2003'):
        #Display iceslabs thickness of the corresponding time period
        lik_blues=ax_plot.scatter(iceslabs_20022018['lon_3413'],
                                  iceslabs_20022018['lat_3413'],
                                  c=iceslabs_20022018['20m_ice_content_m'],
                                  s=10,marker='.',cmap=plt.get_cmap('Blues'),linewidths=0)  
        
    if (timing=='2017-2018'):
        #Inspired from this https://matplotlib.org/stable/gallery/axes_grid1/demo_colorbar_with_inset_locator.html
        axins1 = inset_axes(ax_plot,
                            width="5%",  # width = 50% of parent_bbox width
                            height="100%",  # height : 5%
                            loc='lower left',
                            bbox_to_anchor=(1, 0., 1, 1),
                            bbox_transform=ax_plot.transAxes,
                            borderpad=0)
        
        cbar_blues=plt.colorbar(lik_blues, ax=ax_plot, cax=axins1, shrink=1,orientation='vertical')
        cbar_blues.set_label('Total ice content [m]',fontsize=17)
        cbar_blues.ax.tick_params(labelsize=17)#This is from https://stackoverflow.com/questions/15305737/python-matplotlib-decrease-size-of-colorbar-labels
    return

def plot_pannels_supp(axplot_indiv,axplot_cum,flightlines_20022018,df_firn_aquifer_all,df_all,time_period_function,label_panel_top,label_panel_bottom):
        
    #Generate maps
    draw_map(axplot_indiv,label_panel_top)
    draw_map(axplot_cum,label_panel_bottom)
    
    #Display flightlines
    if (time_period_function=='2010'):
        display_flightlines(axplot_indiv,flightlines_20022018[np.logical_or(flightlines_20022018.str_year==time_period_function,flightlines_20022018.str_year==int(time_period_function))],time_period_function)
    else:
        display_flightlines(axplot_indiv,flightlines_20022018[flightlines_20022018.str_year==time_period_function],time_period_function)
    
    #Display ice slabs
    display_iceslabs(axplot_indiv,df_all[df_all.str_year==time_period_function],time_period_function)
    display_iceslabs(axplot_cum,df_all[df_all.year<=int(time_period_function[-4:])],time_period_function)
    
    return

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
    plot_fig_S6='FALSE'
    plot_standalone_map='TRUE'
    plot_panela='TRUE'
    plot_panelb='TRUE'
    
    ###################### From Tedstone et al., 2022 #####################
    #from plot_map_decadal_change.py
    # Define the CartoPy CRS object.
    crs = ccrs.NorthPolarStereo(central_longitude=-45., true_scale_latitude=70.)
    # This can be converted into a `proj4` string/dict compatible with GeoPandas
    crs_proj4 = crs.proj4_init
    ###################### From Tedstone et al., 2022 #####################
    
    if (plot_fig_S6 == 'TRUE'):
        # -------------------------------- FIG S6 --------------------------------
        plt.rcParams.update({'font.size': 17})
        plt.rcParams["figure.figsize"] = (22,11.3)#from https://pythonguides.com/matplotlib-increase-plot-size/
        fig = plt.figure()
        gs = gridspec.GridSpec(14, 25)
        gs.update(hspace = 2.5)
        gs.update(wspace = 2.5)
        #projection set up from https://stackoverflow.com/questions/33942233/how-do-i-change-matplotlibs-subplot-projection-of-an-existing-axis
        ax1_indiv = plt.subplot(gs[0:7, 0:5],projection=crs)
        ax2_indiv = plt.subplot(gs[0:7, 5:10],projection=crs)
        ax3_indiv = plt.subplot(gs[0:7, 10:15],projection=crs)
        ax4_indiv = plt.subplot(gs[0:7, 15:20],projection=crs)
        ax5_indiv = plt.subplot(gs[0:7, 20:25],projection=crs)
        
        ax1_cum = plt.subplot(gs[7:14, 0:5],projection=crs)
        ax2_cum = plt.subplot(gs[7:14, 5:10],projection=crs)
        ax3_cum = plt.subplot(gs[7:14, 10:15],projection=crs)
        ax4_cum = plt.subplot(gs[7:14, 15:20],projection=crs)
        ax5_cum = plt.subplot(gs[7:14, 20:25],projection=crs)
                
        plot_pannels_supp(ax1_indiv,ax1_cum,flightlines_20022018,df_firn_aquifer_all,df_all,'2002-2003','a','f')
        plot_pannels_supp(ax2_indiv,ax2_cum,flightlines_20022018,df_firn_aquifer_all,df_all,'2010','b','g')
        plot_pannels_supp(ax3_indiv,ax3_cum,flightlines_20022018,df_firn_aquifer_all,df_all,'2011-2012','c','h')
        plot_pannels_supp(ax4_indiv,ax4_cum,flightlines_20022018,df_firn_aquifer_all,df_all,'2013-2014','d','i')
        plot_pannels_supp(ax5_indiv,ax5_cum,flightlines_20022018,df_firn_aquifer_all,df_all,'2017-2018','e','j')
        
        #Add title to 2nd row of plots
        ax1_cum.set_title('2002-2003',fontsize=20, weight='bold')
        ax2_cum.set_title('2002-2010',fontsize=20, weight='bold')
        ax3_cum.set_title('2002-2012',fontsize=20, weight='bold')
        ax4_cum.set_title('2002-2014',fontsize=20, weight='bold')
        ax5_cum.set_title('2002-2018',fontsize=20, weight='bold')
        
        pdb.set_trace()
        
        #Save the figure
        plt.savefig('C:/Users/jullienn/switchdrive/Private/research/RT1/figures/S6/v5/figS6.png',dpi=300,bbox_inches='tight')
        #bbox_inches is from https://stackoverflow.com/questions/32428193/saving-matplotlib-graphs-to-image-as-full-screen)
        # -------------------------------- FIG S6 --------------------------------
    
    if (plot_standalone_map=='TRUE'):
        pdb.set_trace()
        
        #Prepare supp map (old Fig. 1)
        plt.rcParams.update({'font.size': 15})
        plt.rcParams["figure.figsize"] = (22,11.3)#from https://pythonguides.com/matplotlib-increase-plot-size/
        fig = plt.figure()
        
        gs = gridspec.GridSpec(20, 16)
        #projection set up from https://stackoverflow.com/questions/33942233/how-do-i-change-matplotlibs-subplot-projection-of-an-existing-axis
        axmap_standalone = plt.subplot(gs[0:20, 0:16],projection=crs)
        
        #Draw plot of GrIS map
        axmap_standalone.coastlines(edgecolor='black',linewidth=0.075)
        #Display GrIS drainage bassins limits
        GrIS_drainage_bassins.plot(ax=axmap_standalone,color='none', edgecolor='black',linewidth=0.075)
        NO_rignotetal.plot(ax=axmap_standalone,color='none', edgecolor='black',linewidth=0.5)
        NE_rignotetal.plot(ax=axmap_standalone,color='none', edgecolor='black',linewidth=0.5) 
        SE_rignotetal.plot(ax=axmap_standalone,color='none', edgecolor='black',linewidth=0.5) 
        SW_rignotetal.plot(ax=axmap_standalone,color='none', edgecolor='black',linewidth=0.5) 
        CW_rignotetal.plot(ax=axmap_standalone,color='none', edgecolor='black',linewidth=0.5) 
        NW_rignotetal.plot(ax=axmap_standalone,color='none', edgecolor='black',linewidth=0.5)     
        
        #Display scalebar
        scale_bar(axmap_standalone, (0.745, 0.2), 200, 3,5)# axis, location (x,y), length, linewidth, rotation of text
        #by measuring on the screen, the difference in precision between scalebar and length of transects is about ~200m
        
        #Display 2002-2018 flightlines
        axmap_standalone.scatter(flightlines_20022018['lon_3413'],flightlines_20022018['lat_3413'],s=0.1,marker='.',linewidths=0,color='#d9d9d9',label='Flightlines')#,label='2002-2003')
        
        #Display 2010-2014 iceslabs
        axmap_standalone.scatter(df_all[df_all.Track_name.str[:4]=='2010']['lon_3413'],df_all[df_all.Track_name.str[:4]=='2010']['lat_3413'],s=7,marker='.',linewidths=0,color='#4575b4',label='2010-2014 ice slabs')
        axmap_standalone.scatter(df_all[df_all.Track_name.str[:4]=='2011']['lon_3413'],df_all[df_all.Track_name.str[:4]=='2011']['lat_3413'],s=7,marker='.',linewidths=0,color='#4575b4')
        axmap_standalone.scatter(df_all[df_all.Track_name.str[:4]=='2012']['lon_3413'],df_all[df_all.Track_name.str[:4]=='2012']['lat_3413'],s=7,marker='.',linewidths=0,color='#4575b4')
        axmap_standalone.scatter(df_all[df_all.Track_name.str[:4]=='2013']['lon_3413'],df_all[df_all.Track_name.str[:4]=='2013']['lat_3413'],s=7,marker='.',linewidths=0,color='#4575b4')
        axmap_standalone.scatter(df_all[df_all.Track_name.str[:4]=='2014']['lon_3413'],df_all[df_all.Track_name.str[:4]=='2014']['lat_3413'],s=7,marker='.',linewidths=0,color='#4575b4')
        
        #Display 2017-2018 iceslabs
        axmap_standalone.scatter(df_all[df_all.Track_name.str[:4]=='2017']['lon_3413'],df_all[df_all.Track_name.str[:4]=='2017']['lat_3413'],s=3,marker='.',linewidths=0,color='#d73027',label='2017-2018 ice slabs')
        axmap_standalone.scatter(df_all[df_all.Track_name.str[:4]=='2018']['lon_3413'],df_all[df_all.Track_name.str[:4]=='2018']['lat_3413'],s=3,marker='.',linewidths=0,color='#d73027')
                
        #Display 2002-2003 iceslabs
        axmap_standalone.scatter(df_all[df_all.str_year=='2002-2003']['lon_3413'],
                                 df_all[df_all.str_year=='2002-2003']['lat_3413'],
                                 s=12,marker='.',linewidths=0,color='black')
        
        #Display 2002-2003 iceslabs
        axmap_standalone.scatter(df_all[df_all.str_year=='2002-2003']['lon_3413'],
                                 df_all[df_all.str_year=='2002-2003']['lat_3413'],
                                 s=10,marker='.',linewidths=0,color='#ffb300',label='2002-2003 ice slabs')
        
        #Display firn aquifers
        axmap_standalone.scatter(df_firn_aquifer_all['lon_3413'],
                                 df_firn_aquifer_all['lat_3413'],
                                 s=3,marker='.',linewidths=0,color='#238443',label='Firn aquifers')
        
        #Display region name on panel a 
        offset_NO=[-50000,-80000]
        offset_NE=[-50000,-20000]
        offset_SE=[-30000,100000]
        offset_SW=[-40000,-80000]
        offset_CW=[-20000,20000]
        offset_NW=[50000,-50000]
        
        axmap_standalone.text(NO_rignotetal.centroid.x+offset_NO[0],NO_rignotetal.centroid.y+offset_NO[1],np.asarray(NO_rignotetal.SUBREGION1)[0],color='black')
        axmap_standalone.text(NE_rignotetal.centroid.x+offset_NE[0],NE_rignotetal.centroid.y+offset_NE[1],np.asarray(NE_rignotetal.SUBREGION1)[0],color='black')
        axmap_standalone.text(SE_rignotetal.centroid.x+offset_SE[0],SE_rignotetal.centroid.y+offset_SE[1],np.asarray(SE_rignotetal.SUBREGION1)[0],color='black')
        axmap_standalone.text(SW_rignotetal.centroid.x+offset_SW[0],SW_rignotetal.centroid.y+offset_SW[1],np.asarray(SW_rignotetal.SUBREGION1)[0],color='black')
        axmap_standalone.text(CW_rignotetal.centroid.x+offset_CW[0],CW_rignotetal.centroid.y+offset_CW[1],np.asarray(CW_rignotetal.SUBREGION1)[0],color='black')
        axmap_standalone.text(NW_rignotetal.centroid.x+offset_NW[0],NW_rignotetal.centroid.y+offset_NW[1],np.asarray(NW_rignotetal.SUBREGION1)[0],color='black')
                
        # Plot legend. This is from https://stackoverflow.com/questions/24706125/setting-a-fixed-size-for-points-in-legend
        lgnd = axmap_standalone.legend(loc="lower right", scatterpoints=1)
        lgnd.legendHandles[0]._sizes = [30]
        lgnd.legendHandles[1]._sizes = [30]
        lgnd.legendHandles[2]._sizes = [30]
        lgnd.legendHandles[3]._sizes = [30]
        lgnd.legendHandles[4]._sizes = [30]
        
        ###################### From Tedstone et al., 2022 #####################
        #from plot_map_decadal_change.py
        axmap_standalone.set_extent([-634797, 856884, -3345483, -764054], crs=crs)# x0, x1, y0, y1
        gl=axmap_standalone.gridlines(draw_labels=True, xlocs=[-35, -50], ylocs=[65,75], x_inline=False, y_inline=False,linewidth=0.5,linestyle='dashed')
        axmap_standalone.axis('off')
        ###################### From Tedstone et al., 2022 #####################
        
        pdb.set_trace()
        
        plt.savefig('C:/Users/jullienn/switchdrive/Private/research/RT1/figures/fig1/v6/standalone_map.png',dpi=1000,bbox_inches='tight')
        #bbox_inches is from https://stackoverflow.com/questions/32428193/saving-matplotlib-graphs-to-image-as-full-screen)
    
    # --------------------------------- FIG 1 --------------------------------        
    #Prepare Fig. 1
    plt.rcParams.update({'font.size': 15})
    plt.rcParams["figure.figsize"] = (22,11.3)#from https://pythonguides.com/matplotlib-increase-plot-size/
    fig = plt.figure()
    
    gs = gridspec.GridSpec(20, 16)
    gs.update(wspace = 0)
    #gs.update(wspace=0.001)
    #projection set up from https://stackoverflow.com/questions/33942233/how-do-i-change-matplotlibs-subplot-projection-of-an-existing-axis
    axmap = plt.subplot(gs[0:20, 0:11],projection=crs)
    axelev = plt.subplot(gs[0:20, 11:16])

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
        
        axelev.text(ind[0],np.nanmax(max_elev_diff_NO)+55,str(int(np.round(max_elev_diff_NO[4]-max_elev_diff_NO[0])))+' m')
        #axelev.text(ind[1],np.nanmax(max_elev_diff_NW)+180,str(int(np.round(np.nanmax(max_elev_diff_NW)-np.nanmin(max_elev_diff_NW))))+' m')
        axelev.text(ind[2],np.nanmax(max_elev_diff_CW)+30,str(int(np.round(max_elev_diff_CW[4]-max_elev_diff_CW[0])))+' m')
        axelev.text(ind[3],np.nanmax(max_elev_diff_SW)+90,str(int(np.round(max_elev_diff_SW[4]-max_elev_diff_SW[0])))+' m')
        
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
        # -------------------------------- PANEL B --------------------------------    
    
    if (plot_panela=='TRUE'):
        
        hull_computation='FALSE'
        shapefile_display='TRUE'
        
        # -------------------------------- PANEL A --------------------------------
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
        
        if (shapefile_display=='TRUE'):
            
            #Open shapefiles for area change calculations
            #Load high and low estimates ice slabs extent 2010-11-12 and 2010-2018, manually drawn on QGIS
            path_iceslabs_shape='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2002_2018/shapefiles/'
            iceslabs_jullien_highend_20102018_ForCalculations=gpd.read_file(path_iceslabs_shape+'iceslabs_jullien_highend_20102018.shp')
            iceslabs_jullien_highend_2010_11_12_ForCalculations=gpd.read_file(path_iceslabs_shape+'iceslabs_jullien_highend_2010_11_12.shp')
            iceslabs_jullien_lowend_20102018_ForCalculations=gpd.read_file(path_iceslabs_shape+'iceslabs_jullien_lowend_20102018.shp')
            iceslabs_jullien_lowend_2010_11_12_ForCalculations=gpd.read_file(path_iceslabs_shape+'iceslabs_jullien_lowend_2010_11_12.shp')
            
            ### High end ###
            print('--------------- High end ---------------')
            print('Ice slabs extent in 2018:', str(np.round(np.sum(iceslabs_jullien_highend_20102018_ForCalculations.area/(1000000)))), 'km2')
            #Difference in km2. Divide by 1000000 to convert from m2 to km2
            print('Difference 2018 VS 2012:', str(np.round(np.sum(iceslabs_jullien_highend_20102018_ForCalculations.area/(1000000))-np.sum(iceslabs_jullien_highend_2010_11_12_ForCalculations.area/(1000000)),2)), 'km2')
            #Difference in %
            print('Difference 2018 VS 2012:', str((np.sum(iceslabs_jullien_highend_20102018_ForCalculations.area/(1000000))-np.sum(iceslabs_jullien_highend_2010_11_12_ForCalculations.area/(1000000)))/np.sum(iceslabs_jullien_highend_2010_11_12_ForCalculations.area/(1000000))*100),'%')
            ### High end ###

            ### Low end ###
            print('--------------- Low end ---------------')
            print('Ice slabs extent in 2018:', str(np.round(np.sum(iceslabs_jullien_lowend_20102018_ForCalculations.area/(1000000)))), 'km2')
            #Difference in km2. Divide by 1000000 to convert from m2 to km2
            print('Difference 2018 VS 2012:', str(np.round(np.sum(iceslabs_jullien_lowend_20102018_ForCalculations.area/(1000000))-np.sum(iceslabs_jullien_lowend_2010_11_12_ForCalculations.area/(1000000)),2)), 'km2')
            #Difference in %
            print('Difference 2018 VS 2012:', str((np.sum(iceslabs_jullien_lowend_20102018_ForCalculations.area/(1000000))-np.sum(iceslabs_jullien_lowend_2010_11_12_ForCalculations.area/(1000000)))/np.sum(iceslabs_jullien_lowend_2010_11_12_ForCalculations.area/(1000000))*100),'%')
            ### High end ###
                        
            #Draw plot of GrIS map
            axmap.coastlines(edgecolor='black',linewidth=0.075)
            #Display GrIS drainage bassins limits
            GrIS_drainage_bassins.plot(ax=axmap,color='none', edgecolor='black',linewidth=0.075)
            NO_rignotetal.plot(ax=axmap,color='none', edgecolor='black',linewidth=0.5)
            NE_rignotetal.plot(ax=axmap,color='none', edgecolor='black',linewidth=0.5) 
            SE_rignotetal.plot(ax=axmap,color='none', edgecolor='black',linewidth=0.5) 
            SW_rignotetal.plot(ax=axmap,color='none', edgecolor='black',linewidth=0.5) 
            CW_rignotetal.plot(ax=axmap,color='none', edgecolor='black',linewidth=0.5) 
            NW_rignotetal.plot(ax=axmap,color='none', edgecolor='black',linewidth=0.5)     
            axmap.set_extent([-634797, 856884, -3345483, -764054], crs=crs)# x0, x1, y0, y1
            ###################### From Tedstone et al., 2022 #####################
            #from plot_map_decadal_change.py
            gl=axmap.gridlines(draw_labels=True, xlocs=[-35, -50], ylocs=[65,75], x_inline=False, y_inline=False,linewidth=0.5,linestyle='dashed')
            axmap.axis('off')
            ###################### From Tedstone et al., 2022 #####################
            
            #Add panel labels
            axmap.text(0.075, 0.95, 'a',zorder=10, ha='center', va='center', transform=axmap.transAxes, weight='bold',fontsize=25)#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot
            #Display scalebar
            scale_bar(axmap, (0.745, 0.2), 200, 3,5)# axis, location (x,y), length, linewidth, rotation of text
            #by measuring on the screen, the difference in precision between scalebar and length of transects is about ~200m
            
            #Shapefiles
            # --- 2010-2018
            iceslabs_jullien_highend_20102018.plot(ax=axmap,color='#d73027', edgecolor='none',linewidth=0.5)
            
            # --- 2010-11-12
            iceslabs_jullien_highend_2010_11_12.plot(ax=axmap,color='#4575b4', edgecolor='none',linewidth=0.5)
            
            #Flightlines            
            # --- 2013-2014
            axmap.scatter(flightlines_20022018[flightlines_20022018.str_year=='2013-2014']['lon_3413'],
                          flightlines_20022018[flightlines_20022018.str_year=='2013-2014']['lat_3413'],
                          s=0.1,marker='.',linewidths=0,c='#969696',label='flightlines 2013-2014')
            
            # --- 2017-2018
            axmap.scatter(flightlines_20022018[flightlines_20022018.str_year=='2017-2018']['lon_3413'],
                          flightlines_20022018[flightlines_20022018.str_year=='2017-2018']['lat_3413'],
                          s=0.1,marker='.',linewidths=0,c='#969696',label='flightlines 2017-2018')
            
            # --- 2010
            axmap.scatter(flightlines_20022018[np.logical_or(flightlines_20022018.str_year=='2010',flightlines_20022018.str_year==int('2010'))]['lon_3413'],
                          flightlines_20022018[np.logical_or(flightlines_20022018.str_year=='2010',flightlines_20022018.str_year==int('2010'))]['lat_3413'],
                          s=0.1,marker='.',linewidths=0,c='#d9d9d9',label='flightlines 2010')
            # --- 2011-2012
            axmap.scatter(flightlines_20022018[flightlines_20022018.str_year=='2011-2012']['lon_3413'],
                          flightlines_20022018[flightlines_20022018.str_year=='2011-2012']['lat_3413'],
                          s=0.1,marker='.',linewidths=0,c='#d9d9d9',label='flightlines 2011-2012')
            
            #Ice slabs            
            # --- 2002-2003
            axmap.scatter(df_all[df_all.str_year=='2002-2003']['lon_3413'],
                          df_all[df_all.str_year=='2002-2003']['lat_3413'],
                          s=12,marker='.',linewidths=0,color='black',label='2002-2003 ice slabs')
            
            axmap.scatter(df_all[df_all.str_year=='2002-2003']['lon_3413'],
                          df_all[df_all.str_year=='2002-2003']['lat_3413'],
                          s=10,marker='.',linewidths=0,color='#ffb300',label='2002-2003 ice slabs')
            #973de0 purple
            
            #Firn aquifers
            axmap.scatter(df_firn_aquifer_all['lon_3413'],
                          df_firn_aquifer_all['lat_3413'],
                          s=3,marker='.',linewidths=0,color='#238443',label='Firn aquifers')
            
            #Display region name on panel a 
            offset_NO=[-50000,-80000]
            offset_NE=[-50000,-20000]
            offset_SE=[-30000,100000]
            offset_SW=[-40000,-80000]
            offset_CW=[-20000,20000]
            offset_NW=[50000,-50000]
            
            axmap.text(NO_rignotetal.centroid.x+offset_NO[0],NO_rignotetal.centroid.y+offset_NO[1],np.asarray(NO_rignotetal.SUBREGION1)[0],color='black')
            axmap.text(NE_rignotetal.centroid.x+offset_NE[0],NE_rignotetal.centroid.y+offset_NE[1],np.asarray(NE_rignotetal.SUBREGION1)[0],color='black')
            axmap.text(SE_rignotetal.centroid.x+offset_SE[0],SE_rignotetal.centroid.y+offset_SE[1],np.asarray(SE_rignotetal.SUBREGION1)[0],color='black')
            axmap.text(SW_rignotetal.centroid.x+offset_SW[0],SW_rignotetal.centroid.y+offset_SW[1],np.asarray(SW_rignotetal.SUBREGION1)[0],color='black')
            axmap.text(CW_rignotetal.centroid.x+offset_CW[0],CW_rignotetal.centroid.y+offset_CW[1],np.asarray(CW_rignotetal.SUBREGION1)[0],color='black')
            axmap.text(NW_rignotetal.centroid.x+offset_NW[0],NW_rignotetal.centroid.y+offset_NW[1],np.asarray(NW_rignotetal.SUBREGION1)[0],color='black')
            
            #Custom legend myself
            from matplotlib.patches import Patch
            from matplotlib.lines import Line2D
            
            #Custom legend myself            
            legend_elements = [Patch(facecolor='#d73027',label='2010-18 ice slabs extent'),
                               Patch(facecolor='#4575b4',label='2010-12 ice slabs extent'),
                               Line2D([0], [0], marker='.', linestyle='none', label='2002-2003 ice slabs', color='#ffb300'),
                               Line2D([0], [0], marker='.', linestyle='none', label='Firn aquifers', color='#238443'),
                               Line2D([0], [0], color='#969696', lw=2, label='flightlines 2013-14-17-18'),
                               Line2D([0], [0], color='#d9d9d9', lw=2, label='flightlines 2010-11-12')]
            axmap.legend(handles=legend_elements,loc='lower right')
            plt.legend()
            
            #Loop over the region, and extract corresponding total area for 10-11-12 and 10-18 in this region
            for region in list(['NE','NW','CW','SW']):
                #Keep only where name == region
                regional_area_101112=np.sum(iceslabs_jullien_highend_2010_11_12_ForCalculations[iceslabs_jullien_highend_2010_11_12_ForCalculations['region']==region].area/1000000)
                regional_area_1018=np.sum(iceslabs_jullien_highend_20102018_ForCalculations[iceslabs_jullien_highend_20102018_ForCalculations['region']==region].area/1000000)

                #Compute and display area change in % (relative change)
                if (region =='NE'):
                    off_display=[NE_rignotetal.centroid.x+offset_NE[0],NE_rignotetal.centroid.y+offset_NE[1]]
                elif(region =='NO'):
                    off_display=[NO_rignotetal.centroid.x+offset_NO[0],NO_rignotetal.centroid.y+offset_NO[1]]
                elif(region =='NW'):
                    off_display=[NW_rignotetal.centroid.x+offset_NW[0],NW_rignotetal.centroid.y+offset_NW[1]]
                elif(region =='CW'):
                    off_display=[CW_rignotetal.centroid.x+offset_CW[0],CW_rignotetal.centroid.y+offset_CW[1]]
                elif(region =='SW'):
                    off_display=[SW_rignotetal.centroid.x+offset_SW[0],SW_rignotetal.centroid.y+offset_SW[1]]
                else:
                    print('Region not known')
                high_end_change=(regional_area_1018-regional_area_101112)/regional_area_101112*100
                axmap.text(off_display[0]-30000,off_display[1]-80000,'+'+str(int(np.round(high_end_change)))+' %')
    
    # -------------------------------- PANEL A --------------------------------
        
    #Display panels label
    axelev.text(0, 1.05,'b',zorder=10, ha='center', va='center', transform=axelev.transAxes, weight='bold',fontsize=25)#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot
    
    pdb.set_trace()
    
    #Save figure
    plt.savefig('C:/Users/jullienn/switchdrive/Private/research/RT1/figures/fig1/v6/fig1_v5.png',dpi=1000,bbox_inches='tight')
    #bbox_inches is from https://stackoverflow.com/questions/32428193/saving-matplotlib-graphs-to-image-as-full-screen)

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
import rasterio
from rasterio.plot import show
from scalebar import scale_bar

#Set fontsize plot
plt.rcParams.update({'font.size': 10})
### -------------------------- Load GrIS DEM ----------------------------- ###
#https://towardsdatascience.com/reading-and-visualizing-geotiff-images-with-python-8dcca7a74510
path_GrIS_DEM = r'C:/Users/jullienn/switchdrive/Private/research/backup_Aglaja/working_environment/greenland_topo_data/elevations/greenland_dem_mosaic_100m_v3.0.tif'
GrIS_DEM = rasterio.open(path_GrIS_DEM)
### -------------------------- Load GrIS DEM ----------------------------- ###

### -------------------------- Load shapefiles --------------------------- ###
#from https://gis.stackexchange.com/questions/113799/how-to-read-a-shapefile-in-python
path_IceBridgeArea_Shape='C:/Users/jullienn/switchdrive/Private/research/backup_Aglaja/working_environment/greenland_topo_data/IceBridge Area Shapefiles/IceBridge Area Shapefiles/'
IceBridgeArea_Shape=gpd.read_file(path_IceBridgeArea_Shape+'IceBridgeArea_Shape.shp')

#Load Rignot et al., 2016 Greenland drainage bassins
path_rignotetal2016_GrIS_drainage_bassins='C:/Users/jullienn/switchdrive/Private/research/backup_Aglaja/working_environment/greenland_topo_data/GRE_Basins_IMBIE2_v1.3/'
GrIS_drainage_bassins=gpd.read_file(path_rignotetal2016_GrIS_drainage_bassins+'GRE_Basins_IMBIE2_v1.3_EPSG_3413.shp')

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

#Load high estimates ice slabs extent 2010-11-12 and 2010-2018, manually drawn on QGIS
path_iceslabs_shape='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2002_2018/shapefiles/'
iceslabs_jullien_highend_20102018=gpd.read_file(path_iceslabs_shape+'iceslabs_jullien_highend_20102018.shp')
iceslabs_jullien_highend_2010_11_12=gpd.read_file(path_iceslabs_shape+'iceslabs_jullien_highend_2010_11_12.shp')
'''
#Display extent 2010-18 and 2010-11-12
crs = ccrs.NorthPolarStereo(central_longitude=-45., true_scale_latitude=70.) 
fig = plt.figure(figsize=(14,50))
gs = gridspec.GridSpec(7, 25)
gs.update(wspace = 2.5)
ax1 = plt.subplot(gs[0:7, 0:25],projection=crs)
iceslabs_jullien_highend_20102018.plot(ax=ax1,color='blue', edgecolor='black',linewidth=0.5)
iceslabs_jullien_highend_2010_11_12.plot(ax=ax1,color='red', edgecolor='black',linewidth=0.5)
'''
### -------------------------- Load shapefiles --------------------------- ###

### ---------------------------- Load dataset ---------------------------- ###
#Dictionnaries have already been created, load them
path_df_with_elevation='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2002_2018/' 
#Load 2002-2003
f_20022003 = open(path_df_with_elevation+'2002_2003/df_2002_2003_with_elevation_rignotetalregions', "rb")
df_2002_2003 = pickle.load(f_20022003)
f_20022003.close()

#Load cleaned 2010-2018 high estimate
#f_20102018_high = open(path_df_with_elevation+'final_excel/high_estimate/df_20102018_with_elevation_high_estimate_rignotetalregions', "rb")
f_20102018_high = open(path_df_with_elevation+'final_excel/high_estimate/clipped/df_20102018_with_elevation_high_estimate_rignotetalregions_cleaned', "rb")
df_2010_2018_high = pickle.load(f_20102018_high)
f_20102018_high.close()

#Load cleaned 2010-2018 low estimate
#f_20102018_low = open(path_df_with_elevation+'final_excel/low_estimate/df_20102018_with_elevation_low_estimate_rignotetalregions', "rb")
f_20102018_low = open(path_df_with_elevation+'final_excel/low_estimate/clipped/df_20102018_with_elevation_low_estimate_rignotetalregions_cleaned', "rb")
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
path_thickness_likelihood='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2002_2018/final_excel/high_estimate/clipped/'
df_thickness_likelihood_20102018 = pd.read_csv(path_thickness_likelihood+'Ice_Layer_Output_Thicknesses_2010_2018_jullienetal2021_high_estimate_cleaned.csv',delimiter=',',decimal='.')
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
