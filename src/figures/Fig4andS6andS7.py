# -*- coding: utf-8 -*-
"""
Created on Sun Nov 28 12:56:32 2021

@author: jullienn

Code adapted from lenses_thickening_visualisation.py
"""

def compute_distances(eastings,northings):
    #This function is from plot_2002_2003.py, which was originally taken from MacFerrin et al., 2019
    '''Compute the distance (in m here, not km as written originally) of the traces in the file.'''
    # C = sqrt(A^2  + B^2)
    distances = np.power(np.power((eastings[1:] - eastings[:-1]),2) + np.power((northings[1:] - northings[:-1]),2), 0.5)

    #Calculate the cumsum of the distances
    cumsum_distances=np.nancumsum(distances)
    #Seeting the first value of the cumsum to be zero as it is the origin
    return_cumsum_distances=np.zeros(eastings.shape[0])
    return_cumsum_distances[1:eastings.shape[0]]=cumsum_distances

    return return_cumsum_distances


def display_iceslabs_product(dataframe,year,ax_plotting,cmap_year,zorder_indiv,year_background):
    
    #Identify index where time < 20 m
    ind_lower_20m=np.where(dataframe[str(year)]['depth']<20)[0]
    
    #Select x vector and select only where we want data to be displayed
    #Find indexes where within bounds
    indexes_within_bounds=np.logical_and(dataframe[str(year)]['lon_appended']>=-47.4233,dataframe[str(year)]['lon_appended']<=-46.2981)
    #Select only data within bounds
    X=dataframe[str(year)]['lon_appended'][indexes_within_bounds]
    
    #Select only probabilistic data within bounds
    Y=np.arange(0,20,20/dataframe[str(year)]['probabilistic'].shape[0])
    C=dataframe[str(year)]['probabilistic'][:,indexes_within_bounds]

    #Display only where probability is higher or equal than prob
    C_bool=(C>=0.00001).astype(int)
    C_bool_plot=np.zeros([C_bool.shape[0],C_bool.shape[1]])
    C_bool_plot[:]=np.nan
    C_bool_plot[C_bool==1]=1
    
    mask_plot=dataframe[str(year)]['mask'][indexes_within_bounds]
    
    #Create lat/lon vectors for display
    lon3413_plot_int=dataframe[str(year)]['lon_3413'][indexes_within_bounds]
    lat3413_plot_int=dataframe[str(year)]['lat_3413'][indexes_within_bounds]
    
    lon_plot_int=dataframe[str(year)]['lon_appended'][indexes_within_bounds]
    lat_plot_int=dataframe[str(year)]['lat_appended'][indexes_within_bounds]
    
    #Update lat/lon vectors for display with the masks
    lon3413_plot=np.zeros(len(lon3413_plot_int))
    lon3413_plot[:]=np.nan
    lon3413_plot[mask_plot]=lon3413_plot_int[mask_plot]

    lat3413_plot=np.zeros(len(lat3413_plot_int))
    lat3413_plot[:]=np.nan
    lat3413_plot[mask_plot]=lat3413_plot_int[mask_plot]
    
    lon_plot=np.zeros(len(lon_plot_int))
    lon_plot[:]=np.nan
    lon_plot[mask_plot]=lon_plot_int[mask_plot]
    
    lat_plot=np.zeros(len(lat_plot_int))
    lat_plot[:]=np.nan
    lat_plot[mask_plot]=lat_plot_int[mask_plot]
    
    #Calculate distances
    distances_with_start_transect=compute_distances(lon3413_plot,lat3413_plot)
    
    #Display probability
    cb_prob=ax_plotting.pcolor(distances_with_start_transect, Y, C_bool_plot,cmap=cmap_year,alpha=0.8, zorder=zorder_indiv,antialiased=True, linewidth=0.0,label=str(year))
    #for getting rid of mesh lines, this is from https://stackoverflow.com/questions/27092991/white-lines-in-matplotlibs-pcolor
    ax_plotting.invert_yaxis() #Invert the y axis = avoid using flipud.    
    ax_plotting.set_ylim(20,0)
    
    #Add yticks
    ax_plotting.yaxis.tick_left()
    
    #Set xlim
    start_display=10000
    end_display=30000
    ax_plotting.set_xlim(start_display,end_display)
    
    #Display bottom xtick in km instead of m
    xtick_distance=ax_plotting.get_xticks()
    ax_plotting.set_xticks(xtick_distance)
    ax_plotting.set_xticklabels((xtick_distance/1000),fontsize=15)
    

    #Do not dispay x ticks in ax1, ax3
    if (ax_plotting in list([ax1,ax3])):
        #Set xticks labels to empty
        ax_plotting.set_xticklabels([])
        #Set last yticks labels to empty
        ax_plotting.set_yticklabels(['0','5','10','15',''],fontsize=15)
    else:
        #Set fontsize
        ax_plotting.set_yticklabels(['0','5','10','15','20'],fontsize=15)
    
    if (zorder_indiv==1):
        if (ax_plotting==ax5):
            #Custom plot
            ax_plotting.xaxis.tick_bottom()
            ax_plotting.set_xlabel('Distance [km]',fontsize=15)
            
            #Custom legend myself,  line2D from https://stackoverflow.com/questions/39500265/how-to-manually-create-a-legend, marker from https://stackoverflow.com/questions/47391702/how-to-make-a-colored-markers-legend-from-scratch
            legend_elements = [Patch(facecolor=my_pal[year_background],label=str(year_background)),
                               Patch(facecolor=my_pal[year],label=str(year)),
                               Line2D([0], [0], marker='P', linestyle='none', label='KAN_U', color='#b2182b')]
            ax_plotting.legend(handles=legend_elements,loc='lower right',fontsize=15,framealpha=1)
            
        else:
            #Custom legend myself,  line2D from https://stackoverflow.com/questions/39500265/how-to-manually-create-a-legend, marker from https://stackoverflow.com/questions/47391702/how-to-make-a-colored-markers-legend-from-scratch
            legend_elements = [Patch(facecolor=my_pal[year_background],label=str(year_background)),
                               Patch(facecolor=my_pal[year],label=str(year))]
            ax_plotting.legend(handles=legend_elements,loc='lower right',fontsize=15,framealpha=1)     

    elif (zorder_indiv==2):
        
        #Add dashed lines such as in Fig. 3
        ax_plotting.axvline(x=distances_with_start_transect[np.nanargmin(np.abs(np.abs(lon_plot)-np.abs(-47.11)))],zorder=3,linestyle='--',color='k')
        ax_plotting.axvline(x=distances_with_start_transect[np.nanargmin(np.abs(np.abs(lon_plot)-np.abs(-47.023)))],zorder=3,linestyle='--',color='k')
        ax_plotting.axvline(x=distances_with_start_transect[np.nanargmin(np.abs(np.abs(lon_plot)-np.abs(-47.07)))],zorder=3,linestyle='--',color='k',linewidth=1)#Line at km 15.6
        ax_plotting.axvline(x=distances_with_start_transect[np.nanargmin(np.abs(np.abs(lon_plot)-np.abs(-47.0487)))],zorder=3,linestyle='--',color='k',linewidth=1)#Line at km 16.7
        
        #Display KAN_U
        ax_plotting.scatter(distances_with_start_transect[np.nanargmin(np.abs(np.abs(lon_plot)-np.abs(-47.030473)))],0.5,marker='P',s=20,c='#b2182b',zorder=10)
    
    else:
        print('Enter a correct zorder value!')
    
    ax_plotting.grid()

    #display loc on map
    ##### Same procedure of mask appliance for EPSG:32622 coordinates #####
    #Create lat/lon vectors for display
    lon32622_plot_int=dataframe[str(year)]['lon_32622'][indexes_within_bounds]
    lat32622_plot_int=dataframe[str(year)]['lat_32622'][indexes_within_bounds]
    
    #Update lat/lon vectors for display with the masks
    lon32622_plot=np.zeros(len(lon32622_plot_int))
    lon32622_plot[:]=np.nan
    lon32622_plot[mask_plot]=lon32622_plot_int[mask_plot]
    
    lat32622_plot=np.zeros(len(lat32622_plot_int))
    lat32622_plot[:]=np.nan
    lat32622_plot[mask_plot]=lat32622_plot_int[mask_plot]
    ##### Same procedure of mask appliance for EPSG:32622 coordinates #####
        
    #display loc of the trace on zoomed map
    ax_map_region.scatter(lon32622_plot[distances_with_start_transect<=40000],lat32622_plot[distances_with_start_transect<=40000],c='#969696',s=0.1,zorder=10,transform=crs)
    
    ax_map_region.scatter(lon32622_plot[np.logical_and(distances_with_start_transect>=start_display,distances_with_start_transect<=end_display)],
                          lat32622_plot[np.logical_and(distances_with_start_transect>=start_display,distances_with_start_transect<=end_display)],
                          c='k',s=0.1,zorder=11,transform=crs)
    
    ax_map_region.scatter(lon32622_plot[np.logical_and(lon_plot>=-47.11,lon_plot<=-47.023)],
                          lat32622_plot[np.logical_and(lon_plot>=-47.11,lon_plot<=-47.023)],
                          c='#d9d9d9',s=0.1,zorder=12,transform=crs)
    
    '''
    #Same as previously by distance-wise
    ax_map_region.scatter(lon32622_plot[np.logical_and(distances_with_start_transect>=13800,distances_with_start_transect<=17700)],
                          lat32622_plot[np.logical_and(distances_with_start_transect>=13800,distances_with_start_transect<=17700)],
                          c='green',s=0.1,zorder=12,transform=crs)
    '''
    #display loc of the trace on GrIS map
    ax_map_GrIS.scatter(lon3413_plot[distances_with_start_transect<=40000],lat3413_plot[distances_with_start_transect<=40000],c='k',s=5,zorder=10,transform=crs_3413)


def plot_thickness(dictionnary_case_study,dataframe,df_2010_2018_elevation,GrIS_DEM,axt,my_pal):
    #This function is adapted from plot_thickness_evolution from fig_2_paper_iceslabs.py
    
    #Define the probability
    prob=0
    
    ###########################################################################
    ###                    Display thickness evolution                      ###
    ###########################################################################
    
    #Define empty dictionnary for elevation slice definition
    df_for_elev=pd.DataFrame(columns=list(df_2010_2018_elevation.keys()))
    
    #Loop over the years
    for year in dictionnary_case_study.keys():
        if (dictionnary_case_study[year] == 'empty'):
            continue  
        #Select data for the trace
        df_for_elev_temp=df_2010_2018_elevation[df_2010_2018_elevation['Track_name']==dictionnary_case_study[year][0][5:20]+'_'+dictionnary_case_study[year][-1][17:20]]
        #Do not keep where lon<-47.4233 (2012 is the constraining year) and lon> -46.2981 (2018 is the constraining year)
        df_for_elev_temp=df_for_elev_temp[np.logical_and(df_for_elev_temp['lon']>=-47.4233,df_for_elev_temp['lon']<=-46.2981)]
        #2011 data are displayed only from 66.8707 to 67.2
        if (year == 2011):
            df_for_elev_temp=df_for_elev_temp[np.logical_and(df_for_elev_temp['lat']>=66.8707,df_for_elev_temp['lat']<=67.2)]
        #Append data to each other
        df_for_elev=df_for_elev.append(df_for_elev_temp)
     
    #Create an empty df_sampling
    df_sampling=pd.DataFrame(columns=['Track_name','time_period','low_bound', 'high_bound', 'bound_nb', 'mean', 'median', 'q025', 'q075','stddev','rolling_10_median_scatter'])
    
    #Sort df_for_elev from low to high longitude (from west to east)
    df_for_elev_sorted=df_for_elev.sort_values(by=['lon_3413'])
    
    #Create a nan array for storing distances
    df_for_elev_sorted['distances']=np.nan
        
    #Exclude 2010 and 2011 for start of transect definition
    df_for_elev_sorted_transect=df_for_elev_sorted[~np.logical_or(df_for_elev_sorted['year']==2010,df_for_elev_sorted['year']==2011)]
    
    #Store coordinates of the bounds of the transect
    bounds_transect=np.array([[df_for_elev_sorted_transect.iloc[0]['lon_3413'], df_for_elev_sorted_transect.iloc[0]['lat_3413']],
                              [df_for_elev_sorted_transect.iloc[-1]['lon_3413'], df_for_elev_sorted_transect.iloc[-1]['lat_3413']]])
    
    #Compute distance between the westernmost and easternmost point
    bounds_distances=compute_distances(bounds_transect[:,0],bounds_transect[:,1])
    
    #Distance divide every 300m.
    dist_bin_desired=300
    dist_divide=np.arange(np.floor(np.min(bounds_distances)),np.floor(np.max(bounds_distances))+1+dist_bin_desired,dist_bin_desired)
    
    #Define window size for smoothing
    winsize=3
    
    #Calculate distance for every single year
    for indiv_year in np.array([2010,2011,2012,2013,2014,2017,2018]):
        #Extract the indexes of the corresponding year
        ind_indiv_year=np.where(df_for_elev_sorted['year']==indiv_year)
        #Select the corresponding year
        df_trace_year_sorted_for_dist=df_for_elev_sorted.iloc[ind_indiv_year]
                
        if (len(df_trace_year_sorted_for_dist)>0):            
            #Calculate the distance compared to start of transect
            #if 2010 or 2011, calculate the distance along their own path because they are not collocated with 2012-2018
            if (indiv_year in list([2010,2011])):
                
                if (indiv_year==2011):
                    #Sort from low to less low latitude (from south to north)
                    df_trace_year_sorted_for_dist=df_trace_year_sorted_for_dist.sort_values(by=['lat_3413'])
                
                #a. Add the start of the transect for calculating distances
                coordinates_df=[np.asarray(df_trace_year_sorted_for_dist['lon_3413']),
                                np.asarray(df_trace_year_sorted_for_dist['lat_3413'])]
                #b. Calculate the distances
                distances_with_start_transect=compute_distances(coordinates_df[0],coordinates_df[1])
                #c. Store the distances
                df_for_elev_sorted['distances'].iloc[ind_indiv_year]=distances_with_start_transect
                
            else:
                #a. Add the start of the transect for calculating distances
                coordinates_df=[np.append(bounds_transect[0,0],np.asarray(df_trace_year_sorted_for_dist['lon_3413'])),
                                np.append(bounds_transect[0,1],np.asarray(df_trace_year_sorted_for_dist['lat_3413']))]
                #b. Calculate the distances
                distances_with_start_transect=compute_distances(coordinates_df[0],coordinates_df[1])
                #c. Store the distances
                df_for_elev_sorted['distances'].iloc[ind_indiv_year]=distances_with_start_transect[1:]
    
    #Discard here where likelihood<prob
    df_for_elev_sorted=df_for_elev_sorted[df_for_elev_sorted['likelihood']>=prob]
    
    #Define empty list
    app_time_period=[]
    app_low_bound=[]
    app_high_bound=[]
    app_bound_nb=[]
    #app_Track_name=[]
    app_mean=[]
    app_median=[]
    app_q025=[]
    app_q075=[]
    app_stddev=[]
    
    #Loop over every year (2012, 2013, 2014, 2017, 2018) except 2010, 2011
    for indiv_year in range(2012,2019):
        
        #Get data for that specific time period
        df_trace_year_sorted=df_for_elev_sorted[df_for_elev_sorted['year']==indiv_year]
          
        if (len(df_trace_year_sorted)==0):
            #No data in this time period, continue
            continue
        else:
            #Set bound_nb to 0
            bound_nb=0
            #Loop over the dist divide
            for i in range(1,len(dist_divide)):
                
                #Identify low and higher end of the slice
                low_bound=dist_divide[i-1]
                high_bound=dist_divide[i]
        
                #Select all the data belonging to this elev slice
                ind_slice=np.logical_and(np.array(df_trace_year_sorted['distances']>=low_bound),np.array(df_trace_year_sorted['distances']<high_bound))
                df_select=df_trace_year_sorted[ind_slice]
                
                #Append data to each other - general info                
                app_time_period=np.append(app_time_period,np.asarray(indiv_year))
                app_low_bound=np.append(app_low_bound,low_bound)
                app_high_bound=np.append(app_high_bound,high_bound)
                app_bound_nb=np.append(app_bound_nb,str(bound_nb))
                
                if (len(df_select)==0):
                    #Append data to each other - data
                    #app_Track_name=np.append(app_Track_name,np.nan)
                    app_mean=np.append(app_mean,np.nan)
                    app_median=np.append(app_median,np.nan)
                    app_q025=np.append(app_q025,np.nan)
                    app_q075=np.append(app_q075,np.nan)
                    app_stddev=np.append(app_stddev,np.nan)
                else:
                    #Append data to each other -data
                    #app_Track_name=np.append(app_Track_name,np.asarray(df_select['Track_name'].unique()))
                    app_mean=np.append(app_mean,df_select["20m_ice_content_m"].mean())
                    app_median=np.append(app_median,df_select["20m_ice_content_m"].median())
                    app_q025=np.append(app_q025,df_select["20m_ice_content_m"].quantile(q=0.25))
                    app_q075=np.append(app_q075,df_select["20m_ice_content_m"].quantile(q=0.75))
                    app_stddev=np.append(app_stddev,df_select["20m_ice_content_m"].std())

                #Update bound_nb
                bound_nb=bound_nb+1
    
    #Create df_sampling which is the dataframe reuniting all the appended lists
    df_sampling = pd.DataFrame(data={'time_period': app_time_period})
    df_sampling['low_bound']=app_low_bound
    df_sampling['high_bound']=app_high_bound
    df_sampling['bound_nb']=app_bound_nb
    #df_sampling['Track_name']=app_Track_name
    df_sampling['mean']=app_mean
    df_sampling['median']=app_median
    df_sampling['q025']=app_q025
    df_sampling['q075']=app_q075
    df_sampling['stddev']=app_stddev
    df_sampling['rolling_10_median_scatter']=[np.nan]*len(app_mean)
    
    #For ice slabs filling
    summary_filling_year=[]
    summary_filling_data=[]
    summary_filling_low_bound=[]
    summary_filling_high_bound=[]
    
    for time_period in range(2012,2019):
        '''
        if (str(time_period) not in list(['2012','2013','2018'])):
            continue
        '''
        if (len(df_sampling[df_sampling['time_period']==time_period])==0):
            #Empty time period, continue
            continue
        else:
            df_plot=df_sampling[df_sampling['time_period']==time_period]
            
            #Rolling window, size = winsize
            df_plot['rolling_10_median']=df_plot.rolling(winsize, win_type=None,center=True).quantile(quantile=0.5)['median'] 
            df_plot['rolling_10_q025']=df_plot.rolling(winsize, win_type=None,center=True).quantile(quantile=0.5)['q025'] 
            df_plot['rolling_10_q075']=df_plot.rolling(winsize, win_type=None,center=True).quantile(quantile=0.5)['q075'] 
            
            #Where window rolling shows NaN because of NaN surrounding, add raw data
            for i in range(0,len(df_plot)):
                if (np.isnan(np.asarray(df_plot['rolling_10_median'].iloc[i]))):
                    df_plot['rolling_10_median_scatter'].iloc[i]=df_plot['median'].iloc[i]
                    df_plot['rolling_10_q025'].iloc[i]=np.nan#df_plot['q025'].iloc[i]
                    df_plot['rolling_10_q075'].iloc[i]=np.nan#df_plot['q075'].iloc[i]
                
                if (i>0):
                    if (np.isnan(np.asarray(df_plot['rolling_10_median'].iloc[i-1])) and not(np.isnan(np.asarray(df_plot['rolling_10_median'].iloc[i])))):
                        df_plot['rolling_10_median_scatter'].iloc[i]=df_plot['rolling_10_median'].iloc[i]
                        
                if (i<(len(df_plot)-1)):
                    if (np.isnan(np.asarray(df_plot['rolling_10_median'].iloc[i+1])) and not(np.isnan(np.asarray(df_plot['rolling_10_median'].iloc[i])))):
                        df_plot['rolling_10_median_scatter'].iloc[i]=df_plot['rolling_10_median'].iloc[i]
                
            # Plot the median
            axt.plot(df_plot["low_bound"],df_plot["rolling_10_median"],color=my_pal[time_period])
            #Display IQR
            axt.fill_between(df_plot['low_bound'], df_plot['rolling_10_q025'], df_plot['rolling_10_q075'], alpha=0.3,color=my_pal[time_period])
            #Display raw data where moving window do not display
            axt.plot(df_plot['low_bound'], df_plot['rolling_10_median_scatter'],color=my_pal[time_period],alpha=0.5)
            axt.scatter(df_plot['low_bound'], df_plot['rolling_10_median_scatter'],c=my_pal[time_period],marker='.',s=0.5,alpha=1)
            
            #For ice slabs filling
            summary_filling_year=np.append(summary_filling_year,np.ones(len(df_plot['rolling_10_median']))*time_period)
            summary_filling_data=np.append(summary_filling_data,np.asarray(df_plot['rolling_10_median']))
            summary_filling_low_bound=np.append(summary_filling_low_bound,np.asarray(df_plot['low_bound']))
            summary_filling_high_bound=np.append(summary_filling_high_bound,np.asarray(df_plot['high_bound']))

    #Get rid of legend
    #axt.legend_.remove()
    axt.set_xlabel('')
    axt.set_ylabel('')
    
    #Activate ticks x and y label
    axt.yaxis.tick_left()
    axt.xaxis.tick_bottom()
        
    # ---------------------------- Extract elevation ------------------------ #
    #This is from fig_2_paper_iceslabs.py
    #Find closest corresponding elevation. Use only the 2013 transect to ensure we are correctly picking up elevation along a transect
    df_for_elev_sorted_2013=df_for_elev_sorted[df_for_elev_sorted['year']==2013]
    #Define the vectors of latitude and longitude for elevation sampling
    increase_x=7000
    lon_for_elevation_sampling=np.arange(df_for_elev_sorted_2013['lon_3413'].iloc[0],df_for_elev_sorted_2013['lon_3413'].iloc[-1]+increase_x,100)
    #For a and c, the longitude is not constant
    add_yx=(df_for_elev_sorted_2013['lat_3413'].iloc[-1]-df_for_elev_sorted_2013['lat_3413'].iloc[0])/(df_for_elev_sorted_2013['lon_3413'].iloc[-1]-df_for_elev_sorted_2013['lon_3413'].iloc[0])*increase_x
    lat_for_elevation_sampling=np.linspace(df_for_elev_sorted_2013['lat_3413'].iloc[0],df_for_elev_sorted_2013['lat_3413'].iloc[-1]+add_yx,len(lon_for_elevation_sampling))
    
    #Compute distance along this transect
    distances_for_elevation_sampling=compute_distances(lon_for_elevation_sampling,lat_for_elevation_sampling)
        
    #Create the vector for elevations storing
    vect_for_elevation=[]
    
    '''
    #Display on the map to make sure this is correct
    ax8map.scatter(lon_for_elevation_sampling,
                lat_for_elevation_sampling,
                s=5,color='red')
    '''
    #This is from extract_elevation.py
    for i in range(0,len(lon_for_elevation_sampling)):
        #This is from https://gis.stackexchange.com/questions/190423/getting-pixel-values-at-single-point-using-rasterio
        for val in GrIS_DEM.sample([(lon_for_elevation_sampling[i], lat_for_elevation_sampling[i])]):
            #Calculate the corresponding elevation
            vect_for_elevation=np.append(vect_for_elevation,val)
    # ---------------------------- Extract elevation ------------------------ #

    #Place x distance ticks on the top
    axt.xaxis.tick_top()
    
    #Set xlims
    axt.set_xlim(0,40000)
    
    # ---------------------------- Display elevation ------------------------ #
    #This is from fig_2_paper_iceslabs.py
    axt.xaxis.set_ticks(np.arange(0, 41000, 5000)) #ax.xaxis.set_ticks(np.arange(start, end, stepsize))

    #4. Display elevation
    #Store the xticks for the distance
    xtick_distance=axt.get_xticks()
    #Set the xticks
    axt.set_xticks(xtick_distance)
    #Display top xtick in km instead of m
    axt.set_xticklabels((xtick_distance/1000).astype(int))
    
    #This is from https://stackoverflow.com/questions/11244514/modify-tick-label-text
    elevation_display=[np.nan]*len(xtick_distance)
    count=0
    for indiv_dist in xtick_distance:
        if (indiv_dist<0):
            elevation_display[count]=''
        else:
            #Extract index where distance is minimal
            index_closest=np.argmin(np.abs(np.abs(distances_for_elevation_sampling)-np.abs(indiv_dist)))
            #If minimum distance is higher than 1km, store nan. If not, store corresponding elevation
            if (np.abs(np.abs(distances_for_elevation_sampling)-np.abs(indiv_dist))[index_closest] > 1000):
                elevation_display[count]=''
            else:
                elevation_display[count]=np.round(vect_for_elevation[index_closest]).astype(int)
            
        count=count+1
    
    # ---------------------------- Display elevation ------------------------ #    
    #Display elevation on the bottom xticklabels
    #This is from https://stackoverflow.com/questions/19884335/matplotlib-top-bottom-ticks-different "Zaus' reply"
    ax_t = axt.secondary_xaxis('bottom')
    ax_t.set_xticks(xtick_distance)
    ax_t.set_xticklabels(elevation_display)
    ax_t.set_xlabel('Elevation[m]')
    
    #Modify spacing between xticklabels and xticks
    axt.tick_params(pad=1.2)
    ax_t.tick_params(pad=1.2)

    #Use the 2018 dataset for displaying the dashed lines
    lon_for_dashed_lines=df_for_elev_sorted[df_for_elev_sorted['year']==2018]['lon']
    dist_for_dashed_lines=df_for_elev_sorted[df_for_elev_sorted['year']==2018]['distances']
    
    #Display limits of area of focus
    axt.axvline(x=dist_for_dashed_lines.iloc[np.nanargmin(np.abs(np.abs(lon_for_dashed_lines)-np.abs(-47.11)))],zorder=1,linestyle='--',color='k')
    axt.axvline(x=dist_for_dashed_lines.iloc[np.nanargmin(np.abs(np.abs(lon_for_dashed_lines)-np.abs(-47.023)))],zorder=1,linestyle='--',color='k')
    #axt.axvline(x=dist_for_dashed_lines.iloc[np.nanargmin(np.abs(np.abs(lon_for_dashed_lines)-np.abs(-47.07)))],zorder=1,linestyle='--',color='k') #Ice slabs filling, line at km 
    #axt.axvline(x=dist_for_dashed_lines.iloc[np.nanargmin(np.abs(np.abs(lon_for_dashed_lines)-np.abs(-47.0487)))],zorder=1,linestyle='--',color='k') #Ice slabs accretion, line at km x

    #Display KAN_U
    axt.scatter(dist_for_dashed_lines.iloc[np.nanargmin(np.abs(np.abs(lon_for_dashed_lines)-np.abs(-47.030473)))],15.5,marker='P',s=20,c='#b2182b',zorder=10)
    
    '''
    # Hide grid lines, from https://stackoverflow.com/questions/45148704/how-to-hide-axes-and-gridlines-in-matplotlib-python
    axt.grid(False)
    '''
    plt.show()
    
    ###########################################################################
    ###      Extract total columnal ice content inside area of focus        ###
    ###########################################################################
    '''
    #lon coordinates defining in different parts of the transects
    (47.4233,47.11)#0km to first thick dashed line
    (47.11,-47.023)#thick dashed lines
    (-47.023,-46.2981)#right thick dashed line to end of transect
    '''
    coord_parts_transect=np.array([[-47.4233,-47.11],[-47.11,-47.023],[-47.023,-46.2981]])

    count_ice=0
    columnal_sum_studied_case=np.ones((len(np.arange(2009,2022)),len(coord_parts_transect)))*0
    columnal_sum_studied_case[:]=np.nan
    
    for indiv_year in np.arange(2009,2022):
        print(indiv_year)
        if (indiv_year in np.array([2012,2013,2018])):
            #Select data
            df_indiv_year=df_for_elev_sorted[df_for_elev_sorted['year']==indiv_year]
            for count_parts in range(0,len(coord_parts_transect)):
                #Keep only within studied area
                df_studied_case=df_indiv_year[np.logical_and(df_indiv_year['lon']>=coord_parts_transect[count_parts,0],df_indiv_year['lon']<=coord_parts_transect[count_parts,1])]
                #Define the mean delta horizontal dimensions
                delta_horizontal_m = np.mean(np.asarray(df_studied_case['distances'][1:])-np.asarray(df_studied_case['distances'][:-1])) #This is inspired from probabilisitc_iceslabs.py
                #Extract total ice content within this area (in m2 because vertical content [m] * horizontal content [m] #/ distance [m])
                columnal_sum_studied_case[count_ice,count_parts]=np.sum(df_studied_case['20m_ice_content_m']) * delta_horizontal_m #/ (df_studied_case['distances'].iloc[-1]-df_studied_case['distances'].iloc[0]) #if average wanted
        #Update count_ice
        count_ice=count_ice+1
    
    ########################## Ice slabs filling #############################
    #Calculate difference in columnal ice content between 2012 and 2018 for ice slabs filling beneath pre existing one, i.e. between -47.11 and -47.07
    left_end=dist_for_dashed_lines.iloc[np.nanargmin(np.abs(np.abs(lon_for_dashed_lines)-np.abs(-47.11)))]
    right_end=dist_for_dashed_lines.iloc[np.nanargmin(np.abs(np.abs(lon_for_dashed_lines)-np.abs(-47.07)))]
    
    #Create a dataframe
    df_filling = pd.DataFrame(summary_filling_year,columns=['year'])
    df_filling['low_bound']=summary_filling_low_bound
    df_filling['high_bound']=summary_filling_high_bound
    df_filling['ice_content']=summary_filling_data

    #Extract the years
    df_2012_filling=df_filling[df_filling['year']==2012]
    df_2018_filling=df_filling[df_filling['year']==2018]

    #Extract within the bounds
    df_2012_filling_focused=df_2012_filling[np.logical_and(df_2012_filling['low_bound']>=left_end,df_2012_filling['low_bound']<=right_end)]
    df_2018_filling_focused=df_2018_filling[np.logical_and(df_2018_filling['low_bound']>=left_end,df_2018_filling['low_bound']<=right_end)]
    
    '''
    #Display to make sure we extracted th right data
    axt.scatter(df_2012_filling_focused['low_bound'],df_2012_filling_focused['ice_content'],c=my_pal[2012],s=5,alpha=1)
    axt.scatter(df_2018_filling_focused['low_bound'],df_2018_filling_focused['ice_content'],c=my_pal[2018],s=5,alpha=1)
    plt.show()
    '''
    '''
    #Calculate the difference
    diff_20182012_filling=np.asarray(df_2018_filling_focused['ice_content'])-np.asarray(df_2012_filling_focused['ice_content'])
    print('ice filling - median: ',np.median(diff_20182012_filling))
    print('ice filling - mean: ',np.mean(diff_20182012_filling))
    '''
    ########################## Ice slabs filling #############################

    ########################## Ice slabs accretion #############################
    #Calculate difference in columnal ice content between 2012 and 2018 for ice slabs accretion on top pre existing lens/slab, i.e. between -47.0487 and -47.023
    left_end=dist_for_dashed_lines.iloc[np.nanargmin(np.abs(np.abs(lon_for_dashed_lines)-np.abs(-47.0487)))]
    right_end=dist_for_dashed_lines.iloc[np.nanargmin(np.abs(np.abs(lon_for_dashed_lines)-np.abs(-47.023)))]
    
    #Extract within the bounds
    df_2012_accretion_focused=df_2012_filling[np.logical_and(df_2012_filling['low_bound']>=left_end,df_2012_filling['low_bound']<=right_end)]
    df_2018_accretion_focused=df_2018_filling[np.logical_and(df_2018_filling['low_bound']>=left_end,df_2018_filling['low_bound']<=right_end)]
    
    '''
    #Display to make sure we extracted th right data
    axt.scatter(df_2012_accretion_focused['low_bound'],df_2012_accretion_focused['ice_content'],c=my_pal[2012],s=5,alpha=1)
    axt.scatter(df_2018_accretion_focused['low_bound'],df_2018_accretion_focused['ice_content'],c=my_pal[2018],s=5,alpha=1)
    plt.show()
    '''
    '''
    #Calculate the difference
    diff_20182012_accretion=np.asarray(df_2018_accretion_focused['ice_content'])-np.asarray(df_2012_accretion_focused['ice_content'])
    print('ice accretion - median: ',np.median(diff_20182012_accretion))
    print('ice accretion - mean: ',np.mean(diff_20182012_accretion))
    '''
    ########################## Ice slabs accretion #############################
    
    ######################### Spatial variability #############################   
    count_ice_sp=0
    columnal_sum_spatial_variability=np.zeros(len(np.arange(2013,2015)))
    columnal_sum_spatial_variability[:]=np.nan
        
    for indiv_year in np.arange(2013,2015):
        print(indiv_year)
        #Select data
        df_indiv_year=df_for_elev_sorted[df_for_elev_sorted['year']==indiv_year]
        #Keep only within studied area
        df_studied_case=df_indiv_year[np.logical_and(df_indiv_year['lon']>=-47.11,df_indiv_year['lon']<=-47.07)]
        #Define the mean delta horizontal dimensions
        delta_horizontal_m = np.mean(np.asarray(df_studied_case['distances'][1:])-np.asarray(df_studied_case['distances'][:-1])) #This is inspired from probabilisitc_iceslabs.py
        #Extract total ice content within this area (in m2 because vertical content [m] * horizontal content [m] #/ distance [m])
        columnal_sum_spatial_variability[count_ice_sp]=np.sum(df_studied_case['20m_ice_content_m']) * delta_horizontal_m #/ (df_studied_case['distances'].iloc[-1]-df_studied_case['distances'].iloc[0]) #if average wanted
        #Update count_ice_sp
        count_ice_sp=count_ice_sp+1
    
    print('Spatial variability difference (2013-2014)/2014*100: ',(columnal_sum_spatial_variability[0]-columnal_sum_spatial_variability[1])/columnal_sum_spatial_variability[1]*100)
    ######################### Spatial variability #############################

    ###########################################################################
    ###      Extract total columnal ice content inside area of focus        ###
    ###########################################################################  
    
    print('End plotting thickness data')
    
    ###########################################################################
    ###                    Display thickness evolution                      ###
    ###########################################################################
    
    ###########################################################################
    ###                           Display radargrams                        ###
    ###########################################################################
    #Loop over the years
    for year in dictionnary_case_study.keys():
        #Pick up the corresponding datetrack
        date_track=dictionnary_case_study[year][0][5:20]+'_'+dictionnary_case_study[year][-1][17:20]
        
        #Reset depths to 0
        dataframe[str(year)]['depth']=dataframe[str(year)]['depth']-dataframe[str(year)]['depth'][0]
        
        #Select radar slice
        depth_corrected_file=dataframe[str(year)]['radar']
        
        #Identify index where time < 20 m
        ind_lower_20m=np.where(dataframe[str(year)]['depth']<20)[0]
        depth_corrected_20m=depth_corrected_file[ind_lower_20m,:]
        
        #Identify axis for plotting
        if (year==2010):
            ax_plotting=ax1r
            label_for_map=u'\u03B1'
            #Activate ticks xlabel
            ax_plotting.xaxis.tick_bottom()
            #Modify spacing between xticklabels and xticks
            ax_plotting.tick_params(pad=1.2)
            #Define pannel label
            casestudy_nb='a'
        elif (year==2011):
            ax_plotting=ax2r
            label_for_map=u'\u03B2'
            #Activate ticks xlabel
            ax_plotting.xaxis.tick_bottom()
            #Modify spacing between xticklabels and xticks
            ax_plotting.tick_params(pad=1.2)
            #Define pannel label
            casestudy_nb='b'
        elif (year==2012):
            ax_plotting=ax3r
            label_for_map=u'\u03B3'
            #Define pannel label
            casestudy_nb='c'
        elif (year==2013):
            ax_plotting=ax4r
            ax4r.set_xlabel('Depth [m]')
            label_for_map=u'\u03B3'
            #Define pannel label
            casestudy_nb='d'
        elif (year==2014):
            ax_plotting=ax5r
            label_for_map=u'\u03B4'
            #Adapt xticklabels, from https://stackoverflow.com/questions/43673884/change-x-axis-ticks-to-custom-strings
            #Define pannel label
            casestudy_nb='e'
        elif (year==2017):
            ax_plotting=ax6r
            label_for_map=u'\u03B4'
            #Define pannel label
            casestudy_nb='f'
        elif (year==2018):
            ax_plotting=ax7r
            label_for_map=u'\u03B3'
            #Activate ticks xlabel
            ax_plotting.xaxis.tick_bottom()
            #Define pannel label
            casestudy_nb='g'
        else:
            print('Year not existing')
                
        #Select x vector and select only where we want data to be displayed
        #Find indexes where within bounds
        if (year==2011):
            #2011 lat>= 66.8707 and lat <67.2
            indexes_within_bounds=np.logical_and(dataframe[str(year)]['lat_appended']>=66.8707,dataframe[str(year)]['lat_appended']<=67.2)
            #Select only data within bounds
            X=dataframe[str(year)]['lat_appended'][indexes_within_bounds]
        else:
            #2010, 2012-2018 lon>-47.4233 and lon<-46.2981
            indexes_within_bounds=np.logical_and(dataframe[str(year)]['lon_appended']>=-47.4233,dataframe[str(year)]['lon_appended']<=-46.2981)
            #Select only data within bounds
            X=dataframe[str(year)]['lon_appended'][indexes_within_bounds]
                
        #Select only radar data within bounds
        Y_data=np.arange(0,100,100/dataframe[str(year)]['radar'].shape[0])
        C_data=dataframe[str(year)]['radar'][:,indexes_within_bounds]
        
        #Keep only the first 20m of radar data
        C_data=C_data[Y_data<20,:]
        Y_data=Y_data[Y_data<20]

        #Select only probabilistic data within bounds
        Y=np.arange(0,20,20/dataframe[str(year)]['probabilistic'].shape[0])
        C=dataframe[str(year)]['probabilistic'][:,indexes_within_bounds]
    
        #Display only where probability is higher or equal than prob
        C_bool=(C>=0.00001).astype(int)
        C_bool_plot=np.zeros([C_bool.shape[0],C_bool.shape[1]])
        C_bool_plot[:]=np.nan
        C_bool_plot[C_bool==1]=1
        
        mask_plot=dataframe[str(year)]['mask'][indexes_within_bounds]
        
        #Create lat/lon vectors for display
        lon_plot_int=dataframe[str(year)]['lon_appended'][indexes_within_bounds]
        lat_plot_int=dataframe[str(year)]['lat_appended'][indexes_within_bounds]
        
        lon3413_plot_int=dataframe[str(year)]['lon_3413'][indexes_within_bounds]
        lat3413_plot_int=dataframe[str(year)]['lat_3413'][indexes_within_bounds]
        
        #Update lat/lon vectors for display with the masks
        lon_plot=np.zeros(len(lon_plot_int))
        lon_plot[:]=np.nan
        lon_plot[mask_plot]=lon_plot_int[mask_plot]
        
        lat_plot=np.zeros(len(lat_plot_int))
        lat_plot[:]=np.nan
        lat_plot[mask_plot]=lat_plot_int[mask_plot]
        
        lon3413_plot=np.zeros(len(lon3413_plot_int))
        lon3413_plot[:]=np.nan
        lon3413_plot[mask_plot]=lon3413_plot_int[mask_plot]
        
        lat3413_plot=np.zeros(len(lat3413_plot_int))
        lat3413_plot[:]=np.nan
        lat3413_plot[mask_plot]=lat3413_plot_int[mask_plot]
        
        #Calculate distances
        distances_with_start_transect=compute_distances(lon3413_plot,lat3413_plot)
        
        #Display radargram
        cb=ax_plotting.pcolor(distances_with_start_transect, Y_data, C_data,cmap=plt.get_cmap('gray'),zorder=-2)#,norm=divnorm)
        ax_plotting.invert_yaxis() #Invert the y axis = avoid using flipud.    
        ax_plotting.set_ylim(20,0)
        
        #Set yticks
        ax_plotting.set_yticks([0,5,10,15,20])

        #Display probability
        cb_prob=ax_plotting.pcolor(distances_with_start_transect, Y, C_bool_plot,cmap=plt.get_cmap('autumn'),zorder=-1,alpha=0.1, antialiased=True, linewidth=0.0)
        #for getting rid of mesh lines, this is from https://stackoverflow.com/questions/27092991/white-lines-in-matplotlibs-pcolor
        
        #Display bottom xtick in km instead of m
        xtick_distance=np.arange(0, 41000, 5000)
        ax_plotting.set_xticks(xtick_distance) #ax.xaxis.set_ticks(np.arange(start, end, stepsize)) #This is from fig_2_paper_iceslabs.py
        ax_plotting.set_xticklabels((xtick_distance/1000).astype(int))
        
        #Display limits of area of focus
        if (str(year) in list(['2012','2013','2018'])):
            ax_plotting.axvline(x=distances_with_start_transect[np.nanargmin(np.abs(np.abs(lon_plot)-np.abs(-47.11)))],zorder=1,linestyle='--',color='k')
            ax_plotting.axvline(x=distances_with_start_transect[np.nanargmin(np.abs(np.abs(lon_plot)-np.abs(-47.023)))],zorder=1,linestyle='--',color='k')
            
            #Ice slabs filling
            print(year)
            ax_plotting.axvline(x=distances_with_start_transect[np.nanargmin(np.abs(np.abs(lon_plot)-np.abs(-47.07)))],zorder=1,linestyle='--',color='k',linewidth=1)#Line at km 15.6
            print('dist filling: ',distances_with_start_transect[np.nanargmin(np.abs(np.abs(lon_plot)-np.abs(-47.07)))])
            
            ##Ice slabs accretion
            ax_plotting.axvline(x=distances_with_start_transect[np.nanargmin(np.abs(np.abs(lon_plot)-np.abs(-47.0487)))],zorder=1,linestyle='--',color='k',linewidth=1)#Line at km 16.7
            print('dist accretion: ',distances_with_start_transect[np.nanargmin(np.abs(np.abs(lon_plot)-np.abs(-47.0487)))])
            
            #Full transect
            print('dist full end: ',distances_with_start_transect[np.nanargmin(np.abs(np.abs(lon_plot)-np.abs(-47.023)))])
            print('dist full start: ',distances_with_start_transect[np.nanargmin(np.abs(np.abs(lon_plot)-np.abs(-47.11)))])

        ###########################################################################
        ###                           Display radargrams                        ###
        ###########################################################################
                
        #Do not dispay x ticks in 2012, 2013, 2014, 2017
        if (str(year) in list(['2012','2013','2014','2017'])):
            #Fix xticks
            ax_plotting.set_xticks(ax_plotting.get_xticks())
            #Set xticks labels to empty
            ax_plotting.set_xticklabels([])
            #Do not display 20m yticks labels
            ax_plotting.set_yticklabels(['0','5','10','15',''])
        
        ###########################################################################
        ###                       Display data localisation                     ###
        ###########################################################################
        
        if (year==2011):
            ax_plotting.set_xlim(0,40000)
            #Display KAN_U
            ax_plotting.scatter(distances_with_start_transect[np.nanargmin(np.abs(np.abs(lat_plot)-np.abs(67.000425)))],1,s=20,marker='P',c='#b2182b',zorder=10)
        else:
            ax_plotting.set_xlim(0,40000)
            #Display KAN_U
            ax_plotting.scatter(distances_with_start_transect[np.nanargmin(np.abs(np.abs(lon_plot)-np.abs(-47.030473)))],1,s=20,marker='P',c='#b2182b',zorder=10)
        
        ##### Same procedure of mask appliance for EPSG:32622 coordinates #####
        #Create lat/lon vectors for display
        lon32622_plot_int=dataframe[str(year)]['lon_32622'][indexes_within_bounds]
        lat32622_plot_int=dataframe[str(year)]['lat_32622'][indexes_within_bounds]
        
        #Update lat/lon vectors for display with the masks
        lon32622_plot=np.zeros(len(lon32622_plot_int))
        lon32622_plot[:]=np.nan
        lon32622_plot[mask_plot]=lon32622_plot_int[mask_plot]
        
        lat32622_plot=np.zeros(len(lat32622_plot_int))
        lat32622_plot[:]=np.nan
        lat32622_plot[mask_plot]=lat32622_plot_int[mask_plot]
        ##### Same procedure of mask appliance for EPSG:32622 coordinates #####
                
        #display loc on map
        ax8map.scatter(lon32622_plot[distances_with_start_transect<=40000],lat32622_plot[distances_with_start_transect<=40000],c='k',s=0.1,zorder=10,transform=crs)
        
        #Add year on radargram
        ax_plotting.text(0.9535, 0.775,str(year)+', '+label_for_map, color='black',zorder=10, ha='center', va='center', transform=ax_plotting.transAxes,fontsize=15,weight='bold')#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot
        ax_plotting.text(0.9535, 0.775,str(year)+', '+label_for_map, color=my_pal[year],zorder=10, ha='center', va='center', transform=ax_plotting.transAxes,fontsize=15,weight='bold')#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot
        
        #Add pannel label
        ax_plotting.text(0.01, 0.825,casestudy_nb,ha='center', va='center', transform=ax_plotting.transAxes,fontsize=20,zorder=10,weight='bold')#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot
        
        #Activate ticks ylabel
        ax_plotting.yaxis.tick_left()
        ###########################################################################
        ###                       Display data localisation                     ###
        ###########################################################################
        
        ########################## Firn replenishment #########################
        if (str(year) in list(['2010','2011'])):
            continue
        else:
            #Extract top of ice slabs depth
            indexes_top_ice=np.logical_and(dataframe[str(year)]['lon_appended'][indexes_within_bounds]>=-47.0487,dataframe[str(year)]['lon_appended'][indexes_within_bounds]<=-47.023)
            prob_top_ice=dataframe[str(year)]['probabilistic'][:,indexes_within_bounds]
            prob_top_ice=prob_top_ice[:,indexes_top_ice]
            
            top_of_ice=np.zeros((1,prob_top_ice.shape[1]))
            for indiv_col in range (0,prob_top_ice.shape[1]):
                top_of_ice[0,indiv_col]=np.min(np.where(prob_top_ice[:,indiv_col]>0))
            
            top_of_ice=top_of_ice.astype(int)
            #Extract depth of corresponding top of ice
            depth_top_of_ice=dataframe[str(year)]['depth'][top_of_ice[0]]
            
            print(str(year),' top of ice is: ',str(np.mean(depth_top_of_ice)),' +/- ',str(np.std(depth_top_of_ice)))
            print(str(year),' distance_top ice: ',str(distances_with_start_transect[indexes_top_ice][-1]-distances_with_start_transect[indexes_top_ice][0]))
        ########################## Firn replenishment #########################


        ########################## Burrial #########################
        if (str(year) in list(['2010','2011'])):
            continue
        else:
            #Extract bottom of ice slabs depth
            indexes_bottom_ice=np.logical_and(dataframe[str(year)]['lon_appended'][indexes_within_bounds]>=-47.0487,dataframe[str(year)]['lon_appended'][indexes_within_bounds]<=-47.023)
            prob_bottom_ice=dataframe[str(year)]['probabilistic'][:,indexes_within_bounds]
            prob_bottom_ice=prob_bottom_ice[:,indexes_bottom_ice]
            
            bottom_of_ice=np.zeros((1,prob_bottom_ice.shape[1]))
            for indiv_col in range (0,prob_bottom_ice.shape[1]):
                bottom_of_ice[0,indiv_col]=np.max(np.where(prob_bottom_ice[:,indiv_col]>0))
            
            bottom_of_ice=bottom_of_ice.astype(int)
            #Extract depth of corresponding bottom of ice
            depth_bottom_of_ice=dataframe[str(year)]['depth'][bottom_of_ice[0]]
            '''
            #Display where we extracted the values
            ax_plotting.axvline(x=distances_with_start_transect[np.nanargmin(np.abs(np.abs(lon_plot)-np.abs(-47.0487)))],zorder=1,linestyle='--',color='red',linewidth=1)#Line at km 17.7
            ax_plotting.axvline(x=distances_with_start_transect[np.nanargmin(np.abs(np.abs(lon_plot)-np.abs(-47.023)))],zorder=1,linestyle='--',color='red',linewidth=1)#Line at km 16.5
            '''
            print(str(year),' bottom of ice is: ',str(np.mean(depth_bottom_of_ice)),' +/- ',str(np.std(depth_bottom_of_ice)))
            print(str(year),' bottom ice: ',str(distances_with_start_transect[indexes_bottom_ice][-1]-distances_with_start_transect[indexes_bottom_ice][0]))
        ########################## Burrial #########################
    
        ################## Ice slabs thick right at KAN_U #####################
        if (str(year) in list(['2010','2011'])):
            continue
        else:
            #Extract top of ice slabs depth
            indexes_KAN_U=np.logical_and(dataframe[str(year)]['lon_appended'][indexes_within_bounds]>=-47.0329,dataframe[str(year)]['lon_appended'][indexes_within_bounds]<=-47.030473)
            prob_KAN_U=dataframe[str(year)]['probabilistic'][:,indexes_within_bounds]
            prob_KAN_U=prob_KAN_U[:,indexes_KAN_U]
            
            #Define the mean delta vertical dimensions
            delta_vertical_m = np.mean(np.asarray(dataframe[str(year)]['depth'][1:])-np.asarray(dataframe[str(year)]['depth'][:-1])) #This is inspired from probabilisitc_iceslabs.py
            
            KAN_U_ice=np.zeros((1,prob_KAN_U.shape[1]))
            for indiv_col in range (0,prob_KAN_U.shape[1]):
                KAN_U_ice[0,indiv_col]=np.sum(prob_KAN_U[:,indiv_col]>0)*delta_vertical_m
            
            KAN_U_ice_mean=np.mean(KAN_U_ice)
            '''
            #Display where we extracted the values
            ax_plotting.axvline(x=distances_with_start_transect[np.nanargmin(np.abs(np.abs(lon_plot)-np.abs(-47.0329)))],zorder=1,linestyle='--',color='red',linewidth=1)#Line at km 15.6
            ax_plotting.axvline(x=distances_with_start_transect[np.nanargmin(np.abs(np.abs(lon_plot)-np.abs(-47.030473)))],zorder=1,linestyle='--',color='red',linewidth=1)#Line at km 15.6
            '''
            print(str(year),' mean ice right at KAN_U: ',str(KAN_U_ice_mean))
            print(str(year),' distance ice right KAN_U: ',str(distances_with_start_transect[indexes_KAN_U][-1]-distances_with_start_transect[indexes_KAN_U][0]))
        ################## Ice slabs thick right at KAN_U #####################
        
        ####################### Top ice accretion #############################
        if (str(year) in list(['2010','2011','2014','2017'])):
            continue
        else:
            #Extract sector of interest
            indexes_iceacc=np.logical_and(dataframe[str(year)]['lon_appended'][indexes_within_bounds]>=-47.0487,dataframe[str(year)]['lon_appended'][indexes_within_bounds]<=-47.023)
            prob_iceacc=dataframe[str(year)]['probabilistic'][:,indexes_within_bounds]
            prob_iceacc=prob_iceacc[:,indexes_iceacc]
            
            #Define the mean delta vertical dimensions
            delta_vertical_m = np.mean(np.asarray(dataframe[str(year)]['depth'][1:])-np.asarray(dataframe[str(year)]['depth'][:-1])) #This is inspired from probabilisitc_iceslabs.py
            
            iceacc=np.zeros((1,prob_iceacc.shape[1]))
            slab_bottom=np.zeros((1,prob_iceacc.shape[1]))

            for indiv_col in range (0,prob_iceacc.shape[1]):
                #Iceacc
                iceacc[0,indiv_col]=np.sum(prob_iceacc[:,indiv_col]>0)*delta_vertical_m
                #Ice slab bottom depth
                slab_bottom[0,indiv_col]=dataframe[str(year)]['depth'][np.where(prob_iceacc[:,indiv_col]!=0)[0][-1]][0]
            '''
            #Display where we extracted the values
            ax_plotting.axvline(x=distances_with_start_transect[np.nanargmin(np.abs(np.abs(lon_plot)-np.abs(-47.0487)))],zorder=1,linestyle='--',color='red',linewidth=1)#Line at km 17.7
            ax_plotting.axvline(x=distances_with_start_transect[np.nanargmin(np.abs(np.abs(lon_plot)-np.abs(-47.023)))],zorder=1,linestyle='--',color='red',linewidth=1)#Line at km 16.5
            '''
            print('-----> To use:',str(year),' mean ice accretion: ',str(np.mean(iceacc)))
            print('-----> To use:',str(year),' median ice accretion: ',str(np.median(iceacc)))
            print('-----> To use:',str(year),' mean ice bottom depth: ',str(np.mean(slab_bottom)))
            print('-----> To use:',str(year),' median ice bottom depth: ',str(np.median(slab_bottom)))
            print(str(year),' distance ice accretion: ',str(distances_with_start_transect[indexes_iceacc][-1]-distances_with_start_transect[indexes_iceacc][0]))
        ####################### Top ice accretion #############################
        
        ######################## In depth filling #############################
        if (str(year) in list(['2010','2011','2014','2017'])):
            continue
        else:
            #Extract sector of interest
            indexes_icefill=np.logical_and(dataframe[str(year)]['lon_appended'][indexes_within_bounds]>=-47.11,dataframe[str(year)]['lon_appended'][indexes_within_bounds]<=-47.07)
            prob_icefill=dataframe[str(year)]['probabilistic'][:,indexes_within_bounds]
            prob_icefill=prob_icefill[:,indexes_icefill]
            
            #Define the mean delta vertical dimensions
            delta_vertical_m = np.mean(np.asarray(dataframe[str(year)]['depth'][1:])-np.asarray(dataframe[str(year)]['depth'][:-1])) #This is inspired from probabilisitc_iceslabs.py
            
            icefill=np.zeros((1,prob_icefill.shape[1]))
            for indiv_col in range (0,prob_icefill.shape[1]):
                icefill[0,indiv_col]=np.sum(prob_icefill[:,indiv_col]>0)*delta_vertical_m
            
            '''
            #Display where we extracted the values
            ax_plotting.axvline(x=distances_with_start_transect[np.nanargmin(np.abs(np.abs(lon_plot)-np.abs(-47.0329)))],zorder=1,linestyle='--',color='red',linewidth=1)#Line at km 15.6
            ax_plotting.axvline(x=distances_with_start_transect[np.nanargmin(np.abs(np.abs(lon_plot)-np.abs(-47.030473)))],zorder=1,linestyle='--',color='red',linewidth=1)#Line at km 15.6
            '''
            print('-----> To use:',str(year),' mean ice filling: ',str(np.mean(icefill)))
            print('-----> To use:',str(year),' median ice filling: ',str(np.median(icefill)))
            print(str(year),' distance ice filling: ',str(distances_with_start_transect[indexes_iceacc][-1]-distances_with_start_transect[indexes_iceacc][0]))
        ######################## In depth filling #############################
        
        ################ Near surface ice layer thickness #####################
        if (str(year) in list(['2013'])):
            
            #Extract sector of interest
            prob_nsil=dataframe[str(year)]['probabilistic']
            prob_nsil=prob_nsil[:,1800:dataframe[str(year)]['probabilistic'].shape[1]]
            
            #Define the mean delta vertical dimensions
            delta_vertical_m = np.mean(np.asarray(dataframe[str(year)]['depth'][1:])-np.asarray(dataframe[str(year)]['depth'][:-1])) #This is inspired from probabilisitc_iceslabs.py
            
            nsil=np.zeros((1,prob_nsil.shape[1]))
            for indiv_col in range (0,prob_nsil.shape[1]):
                nsil[0,indiv_col]=np.sum(prob_nsil[:,indiv_col]>0)*delta_vertical_m
            
            #Get rid of zeros
            nsil[nsil==0]=np.nan
            
            '''
            #Display where we extracted the values
            ax_plotting.axvline(x=distances_with_start_transect[1800-665],zorder=1,linestyle='--',color='red',linewidth=1)#Line at km 15.6
            '''
            
            print(str(year),' mean near-surface ice layer thickness: ',str(np.nanmean(nsil)))
            print(str(year),' median near-surface ice layer thickness: ',str(np.nanmedian(nsil)))
        ################ Near surface ice layer thickness #####################
        
    return np.min(df_for_elev['elevation']),np.max(df_for_elev['elevation']),columnal_sum_studied_case




import pickle
import scipy.io
import numpy as np
import pdb
import h5py
import matplotlib.pyplot as plt
import pandas as pd
import geopandas as gpd
import matplotlib.gridspec as gridspec
from pyproj import Transformer
import seaborn as sns
sns.set_theme(style="whitegrid")
import cartopy.crs as ccrs
from pyproj import Transformer
import rasterio
import matplotlib.colors #Customize colormap from https://stackoverflow.com/questions/43197474/how-to-customize-the-colorbar-in-python
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import matplotlib.patches as patches
from scalebar import scale_bar

#Define palette as a function of time
#This is from https://www.python-graph-gallery.com/33-control-colors-of-boxplot-seaborn
my_pal = {2010: "k", 2011: "k", 2012: "#80cdc1", 2013: "#542788", 2014: "#d9a764", 2017:"#2166ac",2018:"#b2182b"}

### -------------------------- Load shapefiles ---------------------------- ###
#Load Rignot et al., 2016 Greenland drainage bassins
path_rignotetal2016_GrIS_drainage_bassins='C:/Users/jullienn/switchdrive/Private/research/backup_Aglaja/working_environment/greenland_topo_data/GRE_Basins_IMBIE2_v1.3/'
GrIS_drainage_bassins=gpd.read_file(path_rignotetal2016_GrIS_drainage_bassins+'GRE_Basins_IMBIE2_v1.3_EPSG_3413.shp')
#Extract indiv regions and create related indiv shapefiles
NW_rignotetal=GrIS_drainage_bassins[GrIS_drainage_bassins.SUBREGION1=='NW']
CW_rignotetal=GrIS_drainage_bassins[GrIS_drainage_bassins.SUBREGION1=='CW']
SW_rignotetal=GrIS_drainage_bassins[GrIS_drainage_bassins.SUBREGION1=='SW']
NO_rignotetal=GrIS_drainage_bassins[GrIS_drainage_bassins.SUBREGION1=='NO']
NE_rignotetal=GrIS_drainage_bassins[GrIS_drainage_bassins.SUBREGION1=='NE']
SE_rignotetal=GrIS_drainage_bassins[GrIS_drainage_bassins.SUBREGION1=='SE']
### -------------------------- Load shapefiles ---------------------------- ###

### -------------------------- Load GrIS DEM ------------------------------ ###
#This is from extract_elevation.py
#https://towardsdatascience.com/reading-and-visualizing-geotiff-images-with-python-8dcca7a74510
import rasterio
from rasterio.plot import show

path_GrIS_DEM = r'C:/Users/jullienn/switchdrive/Private/research/backup_Aglaja/working_environment/greenland_topo_data/elevations/greenland_dem_mosaic_100m_v3.0.tif'
GrIS_DEM = rasterio.open(path_GrIS_DEM)
### -------------------------- Load GrIS DEM ------------------------------ ###

### ----------------- Load satellite image in August 2021 ----------------- ###
from pyproj import CRS
import rioxarray as rxr
path_satellite='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/satellite_image/'
#Load satelite data for display
sat_image = rxr.open_rasterio(path_satellite+'T22WFV_20210823T145759_B2348.tif',
                              masked=True).squeeze() #No need to reproject satelite image
#Define extents based on the satelite image
extent_image = [np.asarray(sat_image.x[0]), np.asarray(sat_image.x[-1]), np.asarray(sat_image.y[-1]), np.asarray(sat_image.y[0])]#[west limit, east limit., south limit, north limit]
### ----------------- Load satellite image in August 2021 ----------------- ###

### -------------------- Define data to load and paths -------------------- ###
#Best continuity between overlap through 2012-2018
investigation_year={2010:['Data_20100515_01_007.mat','Data_20100515_01_008.mat','Data_20100515_01_009.mat'],
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

#Define paths
path_data='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/'
path_depth_corrected='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2002_2018/ii_out_from_iceslabs_processing_jullien.py/pickles/'
path_mask='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2002_2018/i_out_from_IceBridgeGPR_Manager_v2.py/pickles_and_images/Boolean Array Picklefiles/'
path_probabilistic='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2002_2018/iii_out_from_probabilistic_iceslabs.py/pickles/'
### -------------------- Define data to load and paths -------------------- ###


### ----------------------- Load and organise data ------------------------ ###
#Define transformer for coordinates transform from "EPSG:4326" to "EPSG:3413"
transformer = Transformer.from_crs("EPSG:4326", "EPSG:3413", always_xy=True)

#Compute the speed (Modified Robin speed):
# self.C / (1.0 + (coefficient*density_kg_m3/1000.0))
v= 299792458 / (1.0 + (0.734*0.873/1000.0))

dataframe={}
for single_year in investigation_year.keys():
    print(single_year)
        
    #If no data, continue
    if (investigation_year[single_year]=='empty'):
        print('No data for year '+str(single_year)+', continue')
        continue
    
    ###1. Load the depth_corrected and probabilistic files:
    start_date_track=investigation_year[single_year][0]
    end_date_track=investigation_year[single_year][-1]
    date_track=start_date_track[5:20]+'_'+end_date_track[17:20]
    
    '''
    if (single_year==2011):
        date_track='20110408_01_087_103'
    if (single_year==2018):
        date_track='20180421_01_004_007'
    #pdb.set_trace()
    '''
    
    filename_depth_corrected=date_track+'_Depth_Corrected_surf_removal_100m.pickle'
    filename_mask=date_track+'_mask.pickle'
    filename_probabilistic=date_track+'_probability_iceslabs_presence_after_DF.pickle'
    
    #Open files
    f_depth_corrected = open(path_depth_corrected+filename_depth_corrected, "rb")
    radar = pickle.load(f_depth_corrected)
    f_depth_corrected.close()
    
    f_mask = open(path_mask+filename_mask, "rb")
    mask = pickle.load(f_mask)
    f_mask.close()
    
    f_probabilistic = open(path_probabilistic+filename_probabilistic, "rb")
    probabilistic_file = pickle.load(f_probabilistic)
    f_probabilistic.close() 
    
    ###2. Load the latitude and longitude
    lat_appended=[]
    lon_appended=[]
    
    for indiv_file_load in investigation_year[single_year]:
        print(indiv_file_load)
        
        #Create the path
        path_raw_data=path_data+str(single_year)+'_Greenland_P3/CSARP_qlook/'+indiv_file_load[5:16]+'/'
        
        #Load data
        if (single_year>=2014):
            
            fdata_filename = h5py.File(path_raw_data+indiv_file_load)
            lat_filename=fdata_filename['Latitude'][:,:]
            lon_filename=fdata_filename['Longitude'][:,:]
            time_filename=fdata_filename['Time'][:,:]
            
        else:
            fdata_filename = scipy.io.loadmat(path_raw_data+indiv_file_load)
            lat_filename = fdata_filename['Latitude']
            lon_filename = fdata_filename['Longitude']
            time_filename = fdata_filename['Time']
            
        #Append data
        lat_appended=np.append(lat_appended,lat_filename)
        lon_appended=np.append(lon_appended,lon_filename)
        
    #Check whether the data are acquired ascending or descending elevation wise.
    #I choose the ascending format. For the data that are descending, reverse them
    #To have ascending data, the longitude should increase
    #(-48 = low elevations, -46 = higher elevations)
    
    if (np.sum(np.diff(lon_appended))<0):
        #It is going toward lower elevations, thus flip left-right
        #(or up-down) all the data!
        
        lat_appended=np.flipud(lat_appended)
        lon_appended=np.flipud(lon_appended)
        radar=np.fliplr(radar)
        mask=np.flipud(mask)
        probabilistic_file=np.fliplr(probabilistic_file)
    
    #Define transformer for coordinates transform from "EPSG:4326" to "EPSG:3413"
    transformer = Transformer.from_crs("EPSG:4326", "EPSG:3413", always_xy=True)
    #Transform the coordinated from WGS84 to EPSG:3413
    #Example from: https://pyproj4.github.io/pyproj/stable/examples.html
    points=transformer.transform(np.array(lon_appended),np.array(lat_appended))
    lon_3413=points[0]
    lat_3413=points[1]
    
    #Define transformer for coordinates transform from "EPSG:4326" to "EPSG:32622" (=UTM 22N)
    transformer = Transformer.from_crs("EPSG:4326", "EPSG:32622", always_xy=True)
    #Transform the coordinated from WGS84 to EPSG:3413
    #Example from: https://pyproj4.github.io/pyproj/stable/examples.html
    points=transformer.transform(np.array(lon_appended),np.array(lat_appended))
    lon_32622=points[0]
    lat_32622=points[1]
    
    #Calculate the depth from the time
    #########################################################################
    # From plot_2002_2003.py - BEGIN
    #########################################################################
    depth_check = v * time_filename / 2.0
    
    #If 2014, transpose the vector
    if (str(single_year)>='2014'):
        depth_check=np.transpose(depth_check)
    
    #Reset times to zero! This is from IceBridgeGPR_Manager_v2.py
    if (depth_check[10]<0):
        #depth_check[10] so that I am sure that the whole vector is negative and
        #not the first as can be for some date were the proccessing is working
        depth_check=depth_check+abs(depth_check[0])
        depth = depth_check
    else:
        depth = depth_check
    
    if (str(single_year) in list(['2011','2012','2014','2017','2018'])):
        if (depth_check[10]>1):
            #depth_check[10] so that I am sure that the whole vector is largely positive and
            #not the first as can be for some date were the proccessing is working
            depth_check=depth_check-abs(depth_check[0])
            depth = depth_check
        
    #Store reunited lat/lon, slice output and mask in a dictionnary:
    dataframe[str(single_year)]={'lat_appended':lat_appended,
                                 'lon_appended':lon_appended,
                                 'lat_3413':lat_3413,
                                 'lon_3413':lon_3413,
                                 'lat_32622':lat_32622,
                                 'lon_32622':lon_32622,
                                 'depth':depth,
                                 'radar':radar,
                                 'mask':mask,
                                 'probabilistic':probabilistic_file}

#Load 2010-2018 elevation dataset
path_df_with_elevation='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2002_2018/final_excel/high_estimate/' 
f_20102018 = open(path_df_with_elevation+'df_20102018_with_elevation_high_estimate_rignotetalregions', "rb")
df_2010_2018_elevation = pickle.load(f_20102018)
f_20102018.close()

#Define transformer for coordinates transform from "EPSG:4326" to "EPSG:32622" (=UTM 22N)
transformer = Transformer.from_crs("EPSG:4326", "EPSG:32622", always_xy=True)

#Transform the coordinated from "EPSG:4326" to "EPSG:32622"
#Example from: https://pyproj4.github.io/pyproj/stable/examples.html
points=transformer.transform(np.array(df_2010_2018_elevation['lon']),np.array(np.array(df_2010_2018_elevation['lat'])))
df_2010_2018_elevation['lon_32622']=points[0]
df_2010_2018_elevation['lat_32622']=points[1]

# Define CRS of the satelite image as being the base for plotting - from FigS9,py
###################### From Tedstone et al., 2022 #####################
#from plot_map_decadal_change.py
# Define the CartoPy CRS object.
crs = ccrs.UTM(zone=22)
# This can be converted into a `proj4` string/dict compatible with GeoPandas
crs_proj4 = crs.proj4_init
###################### From Tedstone et al., 2022 #####################
### ----------------------- Load and organise data ------------------------ ###

###############################################################################
############## Fig. S6 showing radargrams and ice slabs product ###############
###############################################################################
### --------------------------- Prepare figure ---------------------------- ###
plt.rcParams.update({'font.size': 12})
#plt.rcParams["figure.figsize"] = (3.54,5.79)#GRL limits in insches
plt.rcParams["figure.figsize"] = (10.63,17.36)#from https://pythonguides.com/matplotlib-increase-plot-size/

fig = plt.figure()
gs = gridspec.GridSpec(42, 10)  
gs.update(wspace=0.1)
gs.update(hspace=0.1)

ax1r = plt.subplot(gs[0:3, 0:10])
ax2r = plt.subplot(gs[4:7, 0:10])
ax3r = plt.subplot(gs[8:11, 0:10])
ax4r = plt.subplot(gs[11:14, 0:10])
ax5r = plt.subplot(gs[14:17, 0:10])
ax6r = plt.subplot(gs[17:20, 0:10])
ax7r = plt.subplot(gs[20:23, 0:10])
ax11t = plt.subplot(gs[25:31, 0:10])
ax8map = plt.subplot(gs[33:42, 0:5],projection=crs)

# Define the CartoPy CRS object for whole Greenland map display
###################### From Tedstone et al., 2022 #####################
#from plot_map_decadal_change.py
# Define the CartoPy CRS object.
crs_3413 = ccrs.NorthPolarStereo(central_longitude=-45., true_scale_latitude=70.)
###################### From Tedstone et al., 2022 #####################
ax8map_GrIS = plt.subplot(gs[33:42, 5:10],projection=crs_3413)
### --------------------------- Prepare figure ---------------------------- ###

### -------------------------- Display GrIS map  -------------------------- ###
ax8map_GrIS.coastlines(edgecolor='black',linewidth=0.075)
#Display GrIS drainage bassins limits
GrIS_drainage_bassins.plot(ax=ax8map_GrIS,color='none', edgecolor='black',linewidth=0.075)
#Display region name
ax8map_GrIS.text(NO_rignotetal.centroid.x-200000,NO_rignotetal.centroid.y-40000,np.asarray(NO_rignotetal.SUBREGION1)[0])
ax8map_GrIS.text(NE_rignotetal.centroid.x-200000,NE_rignotetal.centroid.y+20000,np.asarray(NE_rignotetal.SUBREGION1)[0])
ax8map_GrIS.text(SE_rignotetal.centroid.x-100000,SE_rignotetal.centroid.y,np.asarray(SE_rignotetal.SUBREGION1)[0])
ax8map_GrIS.text(SW_rignotetal.centroid.x-250000,SW_rignotetal.centroid.y-120000,np.asarray(SW_rignotetal.SUBREGION1)[0])
ax8map_GrIS.text(CW_rignotetal.centroid.x-200000,CW_rignotetal.centroid.y-100000,np.asarray(CW_rignotetal.SUBREGION1)[0])
ax8map_GrIS.text(NW_rignotetal.centroid.x-100000,NW_rignotetal.centroid.y-150000,np.asarray(NW_rignotetal.SUBREGION1)[0])

###################### From Tedstone et al., 2022 #####################
#from plot_map_decadal_change.py
gl=ax8map_GrIS.gridlines(draw_labels=True, xlocs=[-20,-30,-40,-50,-60,-70], ylocs=[60,65,70,75,80], x_inline=False, y_inline=False,linewidth=0.5,linestyle='dashed')
#Customize lat labels
gl.ylabels_right = False
gl.xlabels_top = False
ax8map_GrIS.axis('off')
###################### From Tedstone et al., 2022 #####################

#Display scalebar
scale_bar(ax8map_GrIS, (0.725, 0.125), 200, 3,5)# axis, location (x,y), length, linewidth, rotation of text
#by measuring on the screen, the difference in precision between scalebar and length of transects is about ~200m
ax8map_GrIS.text(-0.45, 0.95,'j',ha='center', va='center', transform=ax8map_GrIS.transAxes,fontsize=20,zorder=10,weight='bold')#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot
### -------------------------- Display GrIS map  -------------------------- ###

### ------------------ Display map of zoom on transect   ------------------ ###
#This section of displaying sat data was coding using tips from
#https://www.earthdatascience.org/courses/use-data-open-source-python/intro-raster-data-python/raster-data-processing/reproject-raster/
#https://towardsdatascience.com/visualizing-satellite-data-using-matplotlib-and-cartopy-8274acb07b84
#Display satelite image
ax8map.imshow(sat_image[3,:,:], extent=extent_image, transform=crs, origin='upper', cmap='Blues_r',zorder=1) #NIR
#Display scalebar
scale_bar(ax8map, (0.8, 0.175), 10, 3,5)# axis, location (x,y), length, linewidth, rotation of text
### ------------------ Display map of zoom on transect   ------------------ ###

### ---------- Display radargrams, ice slabs and thickness plot  ---------- ###
min_elev,max_elev,columnal_sum_studied_case=plot_thickness(investigation_year,dataframe,df_2010_2018_elevation,GrIS_DEM,ax11t,my_pal)
### ---------- Display radargrams, ice slabs and thickness plot  ---------- ###

### ---------------------- Finalise and polish plot ----------------------- ###
ax11t.set_ylabel('Column ice thickness [m]')
#Add pannel label
ax11t.text(0.01, 0.05,'h',ha='center', va='center', transform=ax11t.transAxes,weight='bold',fontsize=20)#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot

#Finalize radargrams plot
ax7r.set_xlabel('Distance [km]')
ax7r.xaxis.labelpad = -1#this is from https://stackoverflow.com/questions/6406368/matplotlib-move-x-axis-label-downwards-but-not-x-axis-ticks

ax4r.set_ylabel('Depth [m]')
#Define transformer for coordinates transform from "EPSG:3413" to "EPSG:32622"
transformer = Transformer.from_crs("EPSG:3413", "EPSG:32622",always_xy=True)

#Display flightlines correspondance with year
coord_u03B2=transformer.transform([-94000],[-2505000])
ax8map.text(coord_u03B2[0][0],coord_u03B2[1][0],u'\u03B2')

coord_u03B1=transformer.transform([-87990],[-2562500])
ax8map.text(coord_u03B1[0][0],coord_u03B1[1][0],u'\u03B1')

coord_u03B3=transformer.transform([-79280],[-2522000])
ax8map.text(coord_u03B3[0][0],coord_u03B3[1][0],u'\u03B3')

coord_uu03B4=transformer.transform([-79280],[-2534500])
ax8map.text(coord_uu03B4[0][0],coord_uu03B4[1][0],s=u'\u03B4')

#Show KAN_U
#Calculate coordinates in "EPSG:32622"
#Define transformer for coordinates transform from "EPSG:4326" to "EPSG:32622"
transformer = Transformer.from_crs("EPSG:4326", "EPSG:32622",always_xy=True)
KAN_U_coord=transformer.transform([-47.0253],[67.0003])
ax8map.scatter(KAN_U_coord[0][0],KAN_U_coord[1][0],s=15,c='#b2182b',label='KAN_U',zorder=10,transform=crs)
ax8map.legend()
#Add pannel label
ax8map.text(-0.2,0.95,'i',ha='center', va='center',transform=ax8map.transAxes, fontsize=20,weight='bold')#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot
ax8map.set_xlim(644858, 699388)
ax8map.set_ylim(7396158, 7459711)

#Display rectangle around datalocation in Greenland map this is from Fig. 3.py   
#Define transformer for coordinates transform from "EPSG:4326" to "EPSG:32622"
transformer = Transformer.from_crs("EPSG:32622", "EPSG:3413",always_xy=True)

coord_origin=transformer.transform([ax8map.get_xlim()[0]],[ax8map.get_ylim()[0]])
coord_topright=transformer.transform([ax8map.get_xlim()[1]],[ax8map.get_ylim()[1]])

#This is from https://stackoverflow.com/questions/37435369/matplotlib-how-to-draw-a-rectangle-on-image
rect = patches.Rectangle((coord_origin[0][0],coord_origin[1][0]),
                         np.abs(coord_origin[0][0]-coord_topright[0][0]),
                         np.abs(coord_origin[1][0]-coord_topright[1][0]),
                         angle=5, linewidth=1, edgecolor='red', facecolor='none')
ax8map_GrIS.add_patch(rect)

###################### From Tedstone et al., 2022 #####################
#from plot_map_decadal_change.py
# x0, x1, y0, y1
gl=ax8map.gridlines(draw_labels=True, xlocs=[-46.5, -47, -47.5], ylocs=[67,66.75], x_inline=False, y_inline=False,linewidth=0.5,linestyle='dashed')
#Customize lat labels
gl.ylabels_right = False
gl.xlabels_top = False
#Set frame in red
plt.setp(ax8map.spines.values(), color='red')#this is from https://stackoverflow.com/questions/7778954/elegantly-changing-the-color-of-a-plot-frame-in-matplotlib
###################### From Tedstone et al., 2022 #####################

#Display legend on thickness plot
#Custom legend myself,  line2D from https://stackoverflow.com/questions/39500265/how-to-manually-create-a-legend, marker from https://stackoverflow.com/questions/47391702/how-to-make-a-colored-markers-legend-from-scratch
legend_elements = [Line2D([0], [0],color=my_pal[2012],label='2012'),
                   Line2D([0], [0],color=my_pal[2013],label='2013'),
                   Line2D([0], [0],color=my_pal[2014],label='2014'),
                   Line2D([0], [0],color=my_pal[2017],label='2017'),
                   Line2D([0], [0],color=my_pal[2018],label='2018'),
                   Line2D([0], [0], marker='P', linestyle='none', label='KAN_U', color='#b2182b')]
                   
ax11t.legend(handles=legend_elements,loc='upper right',fontsize=12)
### ---------------------- Finalise and polish plot ----------------------- ###
'''
#Save figure
plt.savefig('C:/Users/jullienn/switchdrive/Private/research/RT1/figures/fig3/v10/figS6.png',dpi=300,bbox_inches='tight')
#bbox_inches is from https://stackoverflow.com/questions/32428193/saving-matplotlib-graphs-to-image-as-full-screen)
'''

#Calcul ice content change from 0 to 13.8 km: (columnal_sum_studied_case[9,0]-columnal_sum_studied_case[3,0])/columnal_sum_studied_case[3,0]*100
#Calcul ice content change from 13.8 to 17.7 km: (columnal_sum_studied_case[9,1]-columnal_sum_studied_case[3,1])/columnal_sum_studied_case[3,1]*100
#Calcul ice content change from 17.7km to the end of transect: (columnal_sum_studied_case[9,2]-columnal_sum_studied_case[3,2])/columnal_sum_studied_case[3,2]*100
###############################################################################
############## Fig. S6 showing radargrams and ice slabs product ###############
###############################################################################

pdb.set_trace()

###############################################################################
###################### Fig. 4 showing ice slabs product #######################
###############################################################################
#Prepare plot
plt.rcParams.update({'font.size': 15})
plt.rcParams["figure.figsize"] = (22,11.3)#from https://pythonguides.com/matplotlib-increase-plot-size/
fig = plt.figure()
gs = gridspec.GridSpec(14, 14)
gs.update(wspace=0)
gs.update(hspace=0.1)

ax1 = plt.subplot(gs[5:8, 0:14])
ax3 = plt.subplot(gs[8:11, 0:14])
ax5 = plt.subplot(gs[11:14, 0:14])
ax_map_region = plt.subplot(gs[0:4, 0:12],projection=crs)
ax_map_GrIS = plt.subplot(gs[0:4, 12:14],projection=crs_3413)

#Reset depths to 0
for year in np.asarray([2012,2013,2014,2017,2018]):
    dataframe[str(year)]['depth']=dataframe[str(year)]['depth']-dataframe[str(year)]['depth'][0]

#Display ice slabs products
display_iceslabs_product(dataframe,2012,ax1,matplotlib.colors.ListedColormap([my_pal[2012]]),2,'empty')
display_iceslabs_product(dataframe,2013,ax1,matplotlib.colors.ListedColormap([my_pal[2013]]),1,2012)

display_iceslabs_product(dataframe,2014,ax3,matplotlib.colors.ListedColormap([my_pal[2014]]),2,'empty')
display_iceslabs_product(dataframe,2017,ax3,matplotlib.colors.ListedColormap([my_pal[2017]]),1,2014)

display_iceslabs_product(dataframe,2012,ax5,matplotlib.colors.ListedColormap([my_pal[2012]]),2,'empty')
display_iceslabs_product(dataframe,2018,ax5,matplotlib.colors.ListedColormap([my_pal[2018]]),1,2012)

ax1.set_ylabel('Depth [m]',fontsize=15)
ax3.set_ylabel('Depth [m]',fontsize=15)
ax5.set_ylabel('Depth [m]',fontsize=15)

#Display satelite image
ax_map_region.imshow(sat_image[3,:,:], extent=extent_image, transform=crs, origin='upper', cmap='Blues_r',zorder=1) #NIR
'''
path_2016='T22WFV_20160824T145915/S2A_MSIL1C_20160824T145912_N0204_R125_T22WFV_20160824T145915.SAFE/GRANULE/L1C_T22WFV_A006128_20160824T145915/IMG_DATA/'
#Load 2016 satelite data
sat_image = rxr.open_rasterio(path_satellite+path_2016+'T22WFV_20160824T145912_B08.jp2',
                              masked=True).squeeze() #No need to reproject satelite image
#Define extents based on the satelite image
extent_image = [np.asarray(sat_image.x[0]), np.asarray(sat_image.x[-1]), np.asarray(sat_image.y[-1]), np.asarray(sat_image.y[0])]#[west limit, east limit., south limit, north limit]
ax_map_region.imshow(sat_image, extent=extent_image, transform=crs,cmap='Blues_r', origin='upper',zorder=1) #NIR
'''

#Define path where DEM is found
path_DEM='C:/Users/jullienn/switchdrive/Private/research/backup_Aglaja/working_environment/greenland_topo_data/elevations/'
#Display contours - this is fropm Fig.2.py
GrIS_DEM_display_SW_EPSG32622 = rxr.open_rasterio(path_DEM+'SW_zoom/greenland_dem_mosaic_100m_v3.0_SW_EPSG32622.tif',
                              masked=True).squeeze()

extent_DEM_SW_EPSG32622 = [np.asarray(GrIS_DEM_display_SW_EPSG32622.x[0]), np.asarray(GrIS_DEM_display_SW_EPSG32622.x[-1]), np.asarray(GrIS_DEM_display_SW_EPSG32622.y[-1]), np.asarray(GrIS_DEM_display_SW_EPSG32622.y[0])]
cont=ax_map_region.contour(GrIS_DEM_display_SW_EPSG32622[:,:], levels=np.arange(1800,2100,100), extent=extent_DEM_SW_EPSG32622, transform=crs, origin='upper', colors=['#8c510a'],linewidth=0.05,zorder=2)

#Display KAN_U
#Define transformer for coordinates transform from "EPSG:4326" to "EPSG:32622"
transformer = Transformer.from_crs("EPSG:4326", "EPSG:32622",always_xy=True)
KAN_U_coord=transformer.transform([-47.0253],[67.0003])

ax_map_region.scatter(KAN_U_coord[0][0],KAN_U_coord[1][0],s=35,c='#b2182b',label='KAN_U',zorder=10,transform=crs)

#set x and y limits
ax_map_region.set_xlim(649090, 701160)
ax_map_region.set_ylim(7428626, 7437212)

#Add elevation contour values
ax_map_region.text(0.22, -0.11,'1800', ha='center', va='center', rotation=90,transform=ax_map_region.transAxes,fontsize=15,color='#8c510a')#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot
ax_map_region.text(0.525, -0.11,'1900', ha='center', va='center', rotation=90,transform=ax_map_region.transAxes,fontsize=15,color='#8c510a')#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot
ax_map_region.text(0.83, -0.11,'2000', ha='center', va='center', rotation=90,transform=ax_map_region.transAxes,fontsize=15,color='#8c510a')#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot

#Custom legend myself,  line2D from https://stackoverflow.com/questions/39500265/how-to-manually-create-a-legend, marker from https://stackoverflow.com/questions/47391702/how-to-make-a-colored-markers-legend-from-scratch
legend_elements = [Line2D([0], [0], label='Complete transect (0-40 km)', color='#969696'),
                   Line2D([0], [0], label='Displayed transect (10-30 km)', color='black'),
                   Line2D([0], [0], label='~4 km-long transect (13.8-17.7 km)', color='#d9d9d9'),
                   Line2D([0], [0], label='Elevation', color='#8c510a')]#,
                   #Line2D([0], [0], marker='o', linestyle='none', label='KAN_U', color='#b2182b')]
ax_map_region.legend(handles=legend_elements,loc='lower right',fontsize=15)

###################### From Tedstone et al., 2022 #####################
#from plot_map_decadal_change.py
gl=ax_map_region.gridlines(draw_labels=True, xlocs=[-47,-47.5,-46.5], ylocs=[66.95,67], x_inline=False, y_inline=False,linewidth=0.5,linestyle='dashed')
#Customize lat labels
gl.ylabels_right = False
gl.xlabels_top = False
#ax_map_region.axis('off')
###################### From Tedstone et al., 2022 #####################

ax_map_region.axis('off')
#Display yeas on region map
ax_map_region.text(KAN_U_coord[0][0]+16600,KAN_U_coord[1][0]-1050,'2012',color=my_pal[2012],rotation=5,weight='bold')
ax_map_region.text(KAN_U_coord[0][0]+18250,KAN_U_coord[1][0]-900,',',color='black',rotation=5,weight='bold')
ax_map_region.text(KAN_U_coord[0][0]+18600,KAN_U_coord[1][0]-900,'2013',color=my_pal[2013],rotation=5,weight='bold')
ax_map_region.text(KAN_U_coord[0][0]+20200,KAN_U_coord[1][0]-740,',',color='black',rotation=5,weight='bold')
ax_map_region.text(KAN_U_coord[0][0]+20550,KAN_U_coord[1][0]-700,'2018',color=my_pal[2018],rotation=5,weight='bold')

ax_map_region.text(KAN_U_coord[0][0]+18750,KAN_U_coord[1][0]-2600,'2014',color=my_pal[2014],rotation=5,weight='bold')
ax_map_region.text(KAN_U_coord[0][0]+20400,KAN_U_coord[1][0]-2430,',',color='black',rotation=5,weight='bold')
ax_map_region.text(KAN_U_coord[0][0]+20750,KAN_U_coord[1][0]-2430,'2017',color=my_pal[2017],rotation=5,weight='bold')

#Draw plot of GrIS map
ax_map_GrIS.coastlines(edgecolor='black',linewidth=0.075)
#Display GrIS drainage bassins limits
GrIS_drainage_bassins.plot(ax=ax_map_GrIS,color='none', edgecolor='black',linewidth=0.075)
#Display region name
ax_map_GrIS.text(NO_rignotetal.centroid.x-200000,NO_rignotetal.centroid.y-40000,np.asarray(NO_rignotetal.SUBREGION1)[0])
ax_map_GrIS.text(NE_rignotetal.centroid.x-200000,NE_rignotetal.centroid.y+20000,np.asarray(NE_rignotetal.SUBREGION1)[0])
ax_map_GrIS.text(SE_rignotetal.centroid.x-100000,SE_rignotetal.centroid.y,np.asarray(SE_rignotetal.SUBREGION1)[0])
ax_map_GrIS.text(SW_rignotetal.centroid.x-325000,SW_rignotetal.centroid.y-120000,np.asarray(SW_rignotetal.SUBREGION1)[0])
ax_map_GrIS.text(CW_rignotetal.centroid.x-200000,CW_rignotetal.centroid.y-100000,np.asarray(CW_rignotetal.SUBREGION1)[0])
ax_map_GrIS.text(NW_rignotetal.centroid.x-200000,NW_rignotetal.centroid.y-150000,np.asarray(NW_rignotetal.SUBREGION1)[0])
ax_map_GrIS.axis('off')
ax_map_GrIS.set_xlim(-693308, 912076)
ax_map_GrIS.set_ylim(-3440844, -541728)

#Display scalebar
scale_bar(ax_map_region, (0.55, 0.175), 10, 3,5)# axis, location (x,y), length, linewidth, rotation of text
scale_bar(ax_map_GrIS, (0.745, 0.125), 200, 3,5)# axis, location (x,y), length, linewidth, rotation of text
#by measuring on the screen, the difference in precision between scalebar and length of transects is about ~200m

#Add vertical and horizontal arrows to indicate lateral movement of the ice and burrial rate
#Lateral ice motion at KAN_U from Sept 2008 to Sept 2013 = 52.26 +/- 0.01m/year (Doyle et al., 2014)
ax5.arrow(27500,5,-52*6,0,color='black',head_width=1,head_length=100,length_includes_head=True)
ax5.text(25050,5.75,'312 m (52 $m\cdot y^{-1}$)')

#Burrial rate from Spring 2013 to Spring 2017 = ~1.7 (roughly measured with ruler on the Fig. S2a in Rennermalm et al., (2021)).
#This value is probably a bit overestimated because of rounding. However should be okay because we do not capture 2018
ax5.arrow(27500,5,0,1.7,color='black',head_width=100,head_length=1,length_includes_head=True)
ax5.text(27200,9.5,'~1.7 m')

#Add panel label
ax_map_region.text(0.015, 0.85,'a',ha='center', va='center', transform=ax_map_region.transAxes,fontsize=25,zorder=10)#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot
ax1.text(0.0105, 0.85,'   ',backgroundcolor='white',ha='center', va='center', transform=ax1.transAxes,fontsize=17,zorder=10)#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of
ax1.text(0.01, 0.85,'b',ha='center', va='center', transform=ax1.transAxes,fontsize=25,zorder=10)#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot
ax3.text(0.0105, 0.85,'   ',backgroundcolor='white',ha='center', va='center', transform=ax3.transAxes,fontsize=17,zorder=10)#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of
ax3.text(0.01, 0.85,'c',ha='center', va='center', transform=ax3.transAxes,fontsize=25,zorder=10)#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot
ax5.text(0.0105, 0.85,'   ',backgroundcolor='white',ha='center', va='center', transform=ax5.transAxes,fontsize=17,zorder=10)#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of
ax5.text(0.01, 0.85,'d',ha='center', va='center', transform=ax5.transAxes,fontsize=25,zorder=10)#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot

plt.show()

'''
#Save figure
plt.savefig('C:/Users/jullienn/switchdrive/Private/research/RT1/figures/fig3/v10/fig4.png',dpi=300,bbox_inches='tight')
'''
#bbox_inches is from https://stackoverflow.com/questions/32428193/saving-matplotlib-graphs-to-image-as-full-screen
###############################################################################
###################### Fig. 4 showing ice slabs product #######################
###############################################################################

pdb.set_trace()

###############################################################################
############### Supp Fig. 6 PDH and total columnal ice content ################
###############################################################################
fig = plt.figure()
gs = gridspec.GridSpec(5, 10)
gs.update(wspace=0.1)
gs.update(hspace=0.1)
ax10m = plt.subplot(gs[0:5, 0:10])

#Load KAN_U data
path_KAN_U_data='C:/Users/jullienn/switchdrive/Private/research/RT1/KAN_U_data/'
df_KAN_U_csv = pd.read_csv(path_KAN_U_data+'KAN_U_hourly_v03_csv.csv',sep=';',decimal=',',header=0,na_values=-999)

#Calculate PDH
df_KAN_U_csv['PDH']=(df_KAN_U_csv['AirTemperature(C)']>0).astype(int)
df_KAN_U_csv['PDH_temperature']=df_KAN_U_csv[df_KAN_U_csv['AirTemperature(C)']>0]['AirTemperature(C)']

#5. show total cumulative melt
ax = sns.barplot(x="Year", y="PDH_temperature", data=df_KAN_U_csv,palette=['black'],ax=ax10m,estimator=sum,ci=None,alpha=0.8)

#This is from https://stackoverflow.com/questions/14762181/adding-a-y-axis-label-to-secondary-y-axis-in-matplotlib
ax10m_second = ax10m.twinx()
ax10m_second.bar(np.arange(0,13)-0.5,columnal_sum_studied_case[:,1],width=0.15,color='indianred')
ax10m_second.yaxis.set_label_position("right")
ax10m_second.yaxis.tick_right()
ax10m_second.set_ylabel('Total ice content [$m^2$]')
ax10m_second.set_xlim(0,8.6)

ax10m.set_ylabel('cPDH [°C]')
ax10m.set_xlabel('Year')
ax10m.set_xlim(-0.5,8.6) #From 2009 to 2017
#Activate ticks xlabel
ax10m.xaxis.tick_bottom()

#Remove the grids
ax10m.grid(False)
ax10m_second.grid(False)

#Adapt colors
ax10m.spines['left'].set_color('black') #from https://stackoverflow.com/questions/1982770/matplotlib-changing-the-color-of-an-axis
ax10m_second.spines['right'].set_color('indianred') #from https://stackoverflow.com/questions/1982770/matplotlib-changing-the-color-of-an-axis
ax10m_second.tick_params(axis='y', colors='indianred') #from https://stackoverflow.com/questions/1982770/matplotlib-changing-the-color-of-an-axis
ax10m_second.set_ylabel('Total ice content [$m^2$]',color='indianred')
pdb.set_trace()
'''
ax10m.legend_.remove()
plt.show()
'''
pdb.set_trace()

'''
#Save figure
plt.savefig('C:/Users/jullienn/switchdrive/Private/research/RT1/figures/fig3/v9/fig4.png',dpi=300)
'''
###############################################################################
############### Supp Fig. 6 PDH and total columnal ice content ################
###############################################################################