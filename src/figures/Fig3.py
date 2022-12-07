# -*- coding: utf-8 -*-
"""
Created on Sun Dec 19 12:14:06 2021

@author: jullienn
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


def plot_thickness_evolution(dictionnary_case_study,df_2010_2018_csv,df_2010_2018_elevation,DEM_for_elevation,ax1,axt,custom_angle,offset_x,offset_y,casestudy_nb):
    
    #Define empty dictionnary for elevation slice definition
    df_for_elev=pd.DataFrame(columns=list(df_2010_2018_elevation.keys()))
    
    #Loop over the years
    for year in dictionnary_case_study.keys():
        if (dictionnary_case_study[year] == 'empty'):
            continue  
        #Select data for the trace
        df_for_elev_temp=df_2010_2018_elevation[df_2010_2018_elevation['Track_name']==dictionnary_case_study[year][0][5:20]+'_'+dictionnary_case_study[year][-1][17:20]]
        
        #If panel, b, then start only at lon=-530298
        if (casestudy_nb=='b'):
            #Do not keep where lon_3413 < -530298 because not monotoneously elevation increase
            df_for_elev_temp=df_for_elev_temp[df_for_elev_temp['lon_3413']>=-530298]
            
        #If panel c, then start only at lon=-47.9337
        if (casestudy_nb=='c'):
            #Do not keep where lon < -47.9337 because bare ice
            df_for_elev_temp=df_for_elev_temp[df_for_elev_temp['lon']>=-47.9337]      
        
        if (casestudy_nb=='d'):
            #If panel, d, then start only at lon=-47.8226
            #Do not keep where lon < -47.8226 because bare ice
            df_for_elev_temp=df_for_elev_temp[df_for_elev_temp['lon']>=-47.8226]
            
        if (casestudy_nb=='e'):
            #If panel, e, then start only at lon=-47.4233
            #Do not keep where lon < -47.4233 because bare ice
            df_for_elev_temp=df_for_elev_temp[df_for_elev_temp['lon']>=-47.4233]    
        
        if (casestudy_nb=='f'):
            #If panel, e, then start only at lon=-48.2106
            #Do not keep where lon < -48.2106 because bare ice
            df_for_elev_temp=df_for_elev_temp[df_for_elev_temp['lon']>=-48.2106]
        
        #Append data to each other
        df_for_elev=df_for_elev.append(df_for_elev_temp)
                
        #Display data
        ax1.scatter(df_for_elev_temp['lon_3413'],
                    df_for_elev_temp['lat_3413'],
                    s=0.1,color='#737373')
    
    #TO UNCOMMENT FOR NEW GENERATION!!!!!
    '''
    #Save pandas dataframe into excel
    path_transects='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2002_2018/final_excel/transects_Fig2/'
    df_for_elev.to_csv(path_transects+'Fig2_transect_'+str(casestudy_nb)+'.csv',index=False)
    '''
    
    #Display rectangle around data
    x=(np.min(df_for_elev.lon_3413)-offset_x)
    y=(np.min(df_for_elev.lat_3413)-offset_y)
    width=10000
    height=np.sqrt(np.power(abs(np.min(df_for_elev.lon_3413))-abs(np.max(df_for_elev.lon_3413)),2)+np.power(abs(np.min(df_for_elev.lat_3413))-abs(np.max(df_for_elev.lat_3413)),2))+2*offset_x
    #This is from https://stackoverflow.com/questions/37435369/matplotlib-how-to-draw-a-rectangle-on-image
    # Create a Rectangle patch
    rect = patches.Rectangle((x,y),width,height, angle=custom_angle, linewidth=1, edgecolor='blue', facecolor='none')
    # Add the patch to the Axes
    ax1.add_patch(rect)
            
    if (casestudy_nb=='a'):
        #Add number of case study on fig localisation
        ax1.text(x+30000,y-40000,casestudy_nb,color='r',weight='bold',fontsize=20)
    elif (casestudy_nb=='c'):
        #Add number of case study on fig localisation
        ax1.text(x-15000,y+30000,casestudy_nb,color='r',weight='bold',fontsize=20)
    else:
        #Add number of case study on fig localisation
        ax1.text(x-45000,y-20000,casestudy_nb,color='r',weight='bold',fontsize=20)
    
    #Add number of case study on fig localisation    
    axt.text(0.99, 0.8,casestudy_nb, ha='center', va='center', transform=axt.transAxes,fontsize=20)#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot
    
    #Define palette for time periods
    #This is from https://www.python-graph-gallery.com/33-control-colors-of-boxplot-seaborn
    my_pal = {'2010': "#fdd49e", '2011-2012': "#fc8d59", '2013-2014':"#d7301f",'2017-2018':"#7f0000"}
        
    #Create an empty df_sampling
    df_sampling=pd.DataFrame(columns=['Track_name','time_period','low_bound', 'high_bound', 'bound_nb', 'mean', 'median', 'q025', 'q075','stddev','rolling_10_median_scatter'])
        
    #Sort df_for_elev from low to high longitude (from west to east)
    df_for_elev_sorted=df_for_elev.sort_values(by=['lon_3413'])
    
    #Create a nan array for storing distances
    df_for_elev_sorted['distances']=np.nan
    
    #Store coordinates of the bounds of the transect
    bounds_transect=np.array([[df_for_elev_sorted.iloc[0]['lon_3413'], df_for_elev_sorted.iloc[0]['lat_3413']],
                              [df_for_elev_sorted.iloc[-1]['lon_3413'], df_for_elev_sorted.iloc[-1]['lat_3413']]])
    
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
            #a. Add the start of the transect for calculating distances
            coordinates_df=[np.append(bounds_transect[0,0],np.asarray(df_trace_year_sorted_for_dist['lon_3413'])),
                            np.append(bounds_transect[0,1],np.asarray(df_trace_year_sorted_for_dist['lat_3413']))]
            #b. Calculate the distances
            distances_with_start_transect=compute_distances(coordinates_df[0],coordinates_df[1])
            #c. Store the distances
            df_for_elev_sorted['distances'].iloc[ind_indiv_year]=distances_with_start_transect[1:]
        
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
    
    #Loop over the different time periods (2010, 2011-2012, 2013-2014, 2017-2018)
    for time_period in list(['2010','2011-2012','2013-2014','2017-2018']):
        
        #Get data for that specific time period
        if (time_period == '2010'):
            df_trace_year_sorted=df_for_elev_sorted[df_for_elev_sorted['year']==2010]
        elif (time_period == '2011-2012'):
            df_trace_year_sorted=df_for_elev_sorted[(df_for_elev_sorted['year']>=2011) & (df_for_elev_sorted['year']<=2012)]
        elif (time_period == '2013-2014'):
            df_trace_year_sorted=df_for_elev_sorted[(df_for_elev_sorted['year']>=2013) & (df_for_elev_sorted['year']<=2014)]
        elif (time_period == '2017-2018'):
            df_trace_year_sorted=df_for_elev_sorted[(df_for_elev_sorted['year']>=2017) & (df_for_elev_sorted['year']<=2018)]
        else:
            print('Time period not known, break')
            break
        
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
                app_time_period=np.append(app_time_period,np.asarray(time_period))
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
    
    for time_period in list(['2010','2011-2012','2013-2014','2017-2018']):
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
    
    #Get rid of legend
    #axt.legend_.remove()
    axt.set_xlabel('')
    axt.set_ylabel('')
    
    #Activate ticks x and y label
    axt.yaxis.tick_left()
    axt.xaxis.tick_bottom()
    
    #Set fontsiez of y ticks labels to 15
    axt.set_yticklabels(['0','5','10','15'],fontsize=15)
    
    #4. Display elevation
    #Set xlims
    axt.set_xlim(0,70000)
    #Set xticks, this is from https://stackoverflow.com/questions/12608788/changing-the-tick-frequency-on-x-or-y-axis-in-matplotlib
    #start, end = axt.get_xlim()
    if (casestudy_nb=='a'):
        axt.xaxis.set_ticks(np.arange(0, 16000, 5000)) #ax.xaxis.set_ticks(np.arange(start, end, stepsize))
    elif (casestudy_nb=='b'):
        axt.xaxis.set_ticks(np.arange(0, 41000, 5000)) #ax.xaxis.set_ticks(np.arange(start, end, stepsize))
    elif (casestudy_nb=='c'):
        axt.xaxis.set_ticks(np.arange(0, 71000, 5000)) #ax.xaxis.set_ticks(np.arange(start, end, stepsize))
    elif (casestudy_nb=='d'):
        axt.xaxis.set_ticks(np.arange(0, 51000, 5000)) #ax.xaxis.set_ticks(np.arange(start, end, stepsize))
    elif (casestudy_nb=='e'):
        axt.xaxis.set_ticks(np.arange(0, 36000, 5000)) #ax.xaxis.set_ticks(np.arange(start, end, stepsize))
    elif (casestudy_nb=='f'):
        axt.xaxis.set_ticks(np.arange(0, 51000, 5000)) #ax.xaxis.set_ticks(np.arange(start, end, stepsize))
    
    #Store the xticks for the distance
    xtick_distance=axt.get_xticks()
    #Set the xticks
    axt.set_xticks(xtick_distance)
    
    # ---------------------------- Extract elevation ------------------------ #
    #Define the vectors of latitude and longitude for elevation sampling
    increase_x=5000
    lon_for_elevation_sampling=np.arange(df_for_elev_sorted['lon_3413'].iloc[0],df_for_elev_sorted['lon_3413'].iloc[-1]+increase_x,1)
    #For a and c, the longitude is not constant
    add_yx=(df_for_elev_sorted['lat_3413'].iloc[-1]-df_for_elev_sorted['lat_3413'].iloc[0])/(df_for_elev_sorted['lon_3413'].iloc[-1]-df_for_elev_sorted['lon_3413'].iloc[0])*increase_x
    lat_for_elevation_sampling=np.linspace(df_for_elev_sorted['lat_3413'].iloc[0],df_for_elev_sorted['lat_3413'].iloc[-1]+add_yx,len(lon_for_elevation_sampling))
    
    #Compute distance along this transect
    distances_for_elevation_sampling=compute_distances(lon_for_elevation_sampling,lat_for_elevation_sampling)

    #Create the vector for elevations storing
    vect_for_elevation=[]
    '''
    #Display on the map to make sure this is correct
    ax1.scatter(lon_for_elevation_sampling,
                lat_for_elevation_sampling,
                s=0.1,color='red')
    '''
    #This is from extract_elevation.py
    for i in range(0,len(lon_for_elevation_sampling)):
        #This is from https://gis.stackexchange.com/questions/190423/getting-pixel-values-at-single-point-using-rasterio
        for val in GrIS_DEM.sample([(lon_for_elevation_sampling[i], lat_for_elevation_sampling[i])]):
            #Calculate the corresponding elevation
            vect_for_elevation=np.append(vect_for_elevation,val)
    # ---------------------------- Extract elevation ------------------------ #
    
    # ---------------------------- Display elevation ------------------------ #
    #prepare vector of elevation steps
    step_elev=50
    steps_elev=np.arange(math.ceil(vect_for_elevation[0]/50)*50,vect_for_elevation[-1],step_elev)
    #round to closest integer larger or equal number is from https://stackoverflow.com/questions/2356501/how-do-you-round-up-a-number
    
    if (casestudy_nb in list(['a'])):
        #delete the last steps_elev
        steps_elev=steps_elev[0:-1]
    
    #Loop over elevation every 50m, and find the closest distance in km to set top x ticks
    distance_display=np.ones((len(steps_elev),2))*np.nan
    count=0
    for indiv_elev in steps_elev:
        #Extract index where distance is minimal
        index_closest=np.argmin(np.abs(np.abs(vect_for_elevation)-np.abs(indiv_elev)))
        #Store corresponding distance and difference in distance
        distance_display[count,0]=np.round(distances_for_elevation_sampling[index_closest]).astype(int)#distance
        distance_display[count,1]=np.abs(np.abs(vect_for_elevation)-np.abs(indiv_elev))[index_closest]#corresponding difference of elevation between elevation vector and desired elevation
        #Update count
        count=count+1
        
    #Display elevation on the top xticklabels
    #This is from https://stackoverflow.com/questions/19884335/matplotlib-top-bottom-ticks-different "Zaus' reply"
    ax_t = axt.secondary_xaxis('top',color='#8c510a')
    #set xticks distances according to rounded elevations steps
    ax_t.set_xticks(distance_display[:,0])
    #Display elevations
    ax_t.set_xticklabels(steps_elev.astype(int),fontsize=15,color='#8c510a')
    
    #Get rid of elevation = 1250m
    if (casestudy_nb=='a'):
        tick_lab_a=ax_t.get_xticklabels()
        tick_lab_a[1]._text=''
        ax_t.set_xticklabels(tick_lab_a)
    # ---------------------------- Display elevation ------------------------ #
    
    #Display bottom xtick in km instead of m
    axt.set_xticklabels((xtick_distance/1000).astype(int),fontsize=15)
    
    #Modify spacing between xticklabels and xticks
    axt.tick_params(pad=1.2)
    ax_t.tick_params(pad=1.2)
    
    plt.show()
    print('End plotting fig 2')
    
    #pdb.set_trace()
    
    #Calculate summary statistics
    if (casestudy_nb=='a'):
        #Affine distances selection!
        summary_statistics_calculations(df_for_elev_sorted,0,15,list(['2011-2012','2017-2018']),my_pal)
    elif (casestudy_nb=='b'):
        #Affine distances selection!
        summary_statistics_calculations(df_for_elev_sorted,0,14,list(['2010','2017-2018']),my_pal)
    elif (casestudy_nb=='c'):
        #Affine distances selection!
        summary_statistics_calculations(df_for_elev_sorted,39,42,list(['2010','2017-2018']),my_pal)
    elif (casestudy_nb=='d'):
        #Affine distances selection!
        summary_statistics_calculations(df_for_elev_sorted,20.5,37.5,list(['2010','2017-2018']),my_pal)
    elif (casestudy_nb=='e'):
        #Affine distances selection!
        summary_statistics_calculations(df_for_elev_sorted,0,18,list(['2011-2012','2017-2018']),my_pal)
    elif (casestudy_nb=='f'):
        #Affine distances selection!
        summary_statistics_calculations(df_for_elev_sorted,22.5,32,list(['2010','2017-2018']),my_pal)
    else:
        print('Error!')
        pdb.set_trace()
    
    return 


def summary_statistics_calculations(df_casestudy,end_well_developped,end_thickening,time_periods,col_palette):
    
    # create plot
    fig = plt.figure()
    gs = gridspec.GridSpec(9, 9)
    gs.update(wspace=0.25)
    gs.update(wspace=0.25)
    ax1_dist = plt.subplot(gs[0:9, 0:3])
    ax2_dist = plt.subplot(gs[0:9, 3:6])
    ax3_dist = plt.subplot(gs[0:9, 6:9])
    
    #Loop over the different time periods (2010, 2011-2012, 2013-2014, 2017-2018)
    for time_period in time_periods:
        
        #Get data for that specific time period
        if (time_period == '2010'):
            df_casestudy_TimePeriod=df_casestudy[df_casestudy['year']==2010]
        elif (time_period == '2011-2012'):
            df_casestudy_TimePeriod=df_casestudy[(df_casestudy['year']>=2011) & (df_casestudy['year']<=2012)]
        elif (time_period == '2013-2014'):
            df_casestudy_TimePeriod=df_casestudy[(df_casestudy['year']>=2013) & (df_casestudy['year']<=2014)]
        elif (time_period == '2017-2018'):
            df_casestudy_TimePeriod=df_casestudy[(df_casestudy['year']>=2017) & (df_casestudy['year']<=2018)]
        else:
            print('Time period not known, break')
            break
        
        if (len(df_casestudy_TimePeriod)==0):
            #No data in this time period, continue
            continue
        else:
            print('Calculating summary statistics')
            #Select data belonging to the well-developped section
            ax1_dist.hist(df_casestudy_TimePeriod[df_casestudy_TimePeriod['distances']<end_well_developped*1000]['20m_ice_content_m'],bins=np.arange(0,17,1),color=col_palette[time_period])
            #Select data belonging to the thickening section
            ax2_dist.hist(df_casestudy_TimePeriod[np.logical_and(df_casestudy_TimePeriod['distances']>=end_well_developped*1000,
                                                                 df_casestudy_TimePeriod['distances']<end_thickening*1000)]
                                                                 ['20m_ice_content_m'],bins=np.arange(0,17,1),color=col_palette[time_period])
            #Select data belonging to the well-developped section
            ax3_dist.hist(df_casestudy_TimePeriod[df_casestudy_TimePeriod['distances']>=end_thickening*1000]['20m_ice_content_m'],bins=np.arange(0,17,1),color=col_palette[time_period])
            
            pdb.set_trace()
    
    
            #from https://stackoverflow.com/questions/7805552/fitting-a-histogram-with-python and https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.norm.html
            from scipy.stats import norm
            import matplotlib.mlab as mlab
            #from https://towardsdatascience.com/pandas-df-to-numpy-array-c1a9e7d8585f
            data_plot=df_casestudy_TimePeriod[np.logical_and(df_casestudy_TimePeriod['distances']>=end_well_developped*1000,
                                                             df_casestudy_TimePeriod['distances']<end_thickening*1000)]['20m_ice_content_m'].to_numpy().astype(float)
            # best fit of data
            (mu, sigma) = norm.fit(data_plot)
            
            # the histogram of the data
            plt.hist(data_plot, 60, density=True, facecolor='green', alpha=0.75)
            
            # add a 'best fit' line
            y = norm.pdf( np.arange(0,17,1), mu, sigma)
            l = plt.plot(np.arange(0,17,1), y, 'r--', linewidth=2)
            
            #plot
            plt.xlabel('Smarts')
            plt.ylabel('Probability')
            plt.title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=%.3f,\ \sigma=%.3f$' %(mu, sigma))
            plt.grid(True)
            
            plt.show()
 
            from scipy.stats import norm
            fig, ax = plt.subplots(1, 1)
            

            ax.hist(data_plot, density=True, histtype='stepfilled', alpha=0.2)
            ax.plot(data_plot,norm.pdf(data_plot),'r-', lw=5, alpha=0.6, label='norm pdf')
            
            
            
    
    
    
    
    return



###     This is from iceslabs_20102018_thickening_analysis.py       ###

#Import librairies
import datetime
from scipy import spatial
import pandas as pd
from pyproj import Transformer
import numpy as np
import pdb
import matplotlib.pyplot as plt
import geopandas as gpd
import pickle
import matplotlib.patches as patches
import matplotlib as mpl
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
import seaborn as sns
sns.set_theme(style="whitegrid")
from scipy import signal
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
from scalebar import scale_bar
import math

#Set fontsize plot
plt.rcParams.update({'font.size': 10})

### -------------------------- Load shapefiles --------------------------- ###
#Load Rignot et al., 2016 Greenland drainage bassins
path_rignotetal2016_GrIS_drainage_bassins='C:/Users/jullienn/switchdrive/Private/research/backup_Aglaja/working_environment/greenland_topo_data/GRE_Basins_IMBIE2_v1.3/'
GrIS_drainage_bassins=gpd.read_file(path_rignotetal2016_GrIS_drainage_bassins+'GRE_Basins_IMBIE2_v1.3_EPSG_3413.shp') #the regions are the last rows of the shapefile

#Extract indiv regions and create related indiv shapefiles
NO_rignotetal=GrIS_drainage_bassins[GrIS_drainage_bassins.SUBREGION1=='NO']
NE_rignotetal=GrIS_drainage_bassins[GrIS_drainage_bassins.SUBREGION1=='NE']
SE_rignotetal=GrIS_drainage_bassins[GrIS_drainage_bassins.SUBREGION1=='SE']
SW_rignotetal=GrIS_drainage_bassins[GrIS_drainage_bassins.SUBREGION1=='SW']
CW_rignotetal=GrIS_drainage_bassins[GrIS_drainage_bassins.SUBREGION1=='CW']
NW_rignotetal=GrIS_drainage_bassins[GrIS_drainage_bassins.SUBREGION1=='NW']
### -------------------------- Load shapefiles --------------------------- ###

path='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2002_2018/final_excel/dataset_for_Fig3/clipped/'
#Load all 2010-2018 data
df_2010_2018_csv = pd.read_csv(path+'Ice_Layer_Output_Thicknesses_2010_2018_jullienetal2021_Fig3_high_estimate_cleaned.csv',delimiter=',',decimal='.')
#Transform the coordinated from WGS84 to EPSG:3413
transformer = Transformer.from_crs("EPSG:4326", "EPSG:3413", always_xy=True)
points=transformer.transform(np.asarray(df_2010_2018_csv["lon"]),np.asarray(df_2010_2018_csv["lat"]))

#Store lat/lon in 3413
df_2010_2018_csv['lon_3413']=points[0]
df_2010_2018_csv['lat_3413']=points[1]

### -------------------------- Load GrIS DEM ----------------------------- ###
#This is from extract_elevation.py
#https://towardsdatascience.com/reading-and-visualizing-geotiff-images-with-python-8dcca7a74510
import rasterio
from rasterio.plot import show

path_GrIS_DEM = r'C:/Users/jullienn/switchdrive/Private/research/backup_Aglaja/working_environment/greenland_topo_data/elevations/greenland_dem_mosaic_100m_v3.0.tif'
GrIS_DEM = rasterio.open(path_GrIS_DEM)
### -------------------------- Load GrIS DEM ----------------------------- ###
#Plot 2010, 2011, 2012, 2013, 2014 ,2017 2018, select overlapping case study: use clean and clear ice slabs tramsects

#CHOOSE LOC1, loc 2, loc 3, loc 15, loc 23

loc1={2010:['Data_20100507_01_008.mat','Data_20100507_01_009.mat','Data_20100507_01_010.mat'],
      2011:['Data_20110426_01_009.mat','Data_20110426_01_010.mat','Data_20110426_01_011.mat'],
      2012:'empty',
      2013:'empty',
      2014:['Data_20140421_01_009.mat','Data_20140421_01_010.mat','Data_20140421_01_011.mat','Data_20140421_01_012.mat','Data_20140421_01_013.mat'],
      2017:['Data_20170424_01_008.mat','Data_20170424_01_009.mat','Data_20170424_01_010.mat','Data_20170424_01_011.mat','Data_20170424_01_012.mat','Data_20170424_01_013.mat','Data_20170424_01_014.mat'],
      2018:'empty'}

loc2={2010:['Data_20100513_01_001.mat','Data_20100513_01_002.mat'],
      2011:['Data_20110411_01_116.mat','Data_20110411_01_117.mat','Data_20110411_01_118.mat'],
      2012:['Data_20120428_01_125.mat','Data_20120428_01_126.mat'],
      2013:'empty',
      2014:['Data_20140408_11_024.mat','Data_20140408_11_025.mat','Data_20140408_11_026.mat'],
      2017:['Data_20170508_02_165.mat','Data_20170508_02_166.mat','Data_20170508_02_167.mat','Data_20170508_02_168.mat','Data_20170508_02_169.mat','Data_20170508_02_170.mat','Data_20170508_02_171.mat'],
      2018:'empty'}

#This one is collocated with FS1, 2, 3.
loc3={2010:'empty',
      2011:'empty',
      2012:['Data_20120423_01_137.mat','Data_20120423_01_138.mat'],
      2013:['Data_20130409_01_010.mat','Data_20130409_01_011.mat','Data_20130409_01_012.mat'],
      2014:'empty',
      2017:'empty',
      2018:['Data_20180421_01_004.mat','Data_20180421_01_005.mat','Data_20180421_01_006.mat','Data_20180421_01_007.mat']}


#in 2017 overestimation of ice content
loc4={2010:['Data_20100512_04_073.mat','Data_20100512_04_074.mat'],
      2011:'empty',
      2012:'empty',
      2013:'empty',
      2014:'empty',
      2017:['Data_20170421_01_171.mat','Data_20170421_01_172.mat','Data_20170421_01_173.mat','Data_20170421_01_174.mat'],
      2018:['Data_20180425_01_166.mat','Data_20180425_01_167.mat','Data_20180425_01_168.mat','Data_20180425_01_169.mat']}

loc6={2010:'empty',
      2011:['Data_20110516_01_009.mat','Data_20110516_01_010.mat'],
      2012:'empty',
      2013:['Data_20130402_01_008.mat'],
      2014:'empty',
      2017:['Data_20170412_01_075.mat','Data_20170412_01_076.mat'],
      2018:'empty'}

loc8={2010:['Data_20100517_02_001.mat','Data_20100517_02_002.mat'],
      2011:['Data_20110502_01_171.mat'],
      2012:['Data_20120516_01_002.mat'],
      2013:['Data_20130419_01_004.mat','Data_20130419_01_005.mat'],
      2014:['Data_20140507_03_007.mat','Data_20140507_03_008.mat'], #test with 20140514_02_087_089 and 20140515_02_173_175 also
      2017:['Data_20170417_01_171.mat','Data_20170417_01_172.mat','Data_20170417_01_173.mat','Data_20170417_01_174.mat'],
      2018:'empty'}

loc9={2010:['Data_20100508_01_114.mat','Data_20100508_01_115.mat'],
      2011:['Data_20110419_01_008.mat','Data_20110419_01_009.mat','Data_20110419_01_010.mat'],
      2012:['Data_20120418_01_129.mat','Data_20120418_01_130.mat','Data_20120418_01_131.mat'],
      2013:['Data_20130405_01_165.mat','Data_20130405_01_166.mat','Data_20130405_01_167.mat'],
      2014:['Data_20140424_01_002.mat','Data_20140424_01_003.mat','Data_20140424_01_004.mat'],
      2017:['Data_20170422_01_168.mat','Data_20170422_01_169.mat','Data_20170422_01_170.mat','Data_20170422_01_171.mat'],
      2018:['Data_20180427_01_170.mat','Data_20180427_01_171.mat','Data_20180427_01_172.mat']}

###################### From Tedstone et al., 2022 #####################
#from plot_map_decadal_change.py
# Define the CartoPy CRS object.
crs = ccrs.NorthPolarStereo(central_longitude=-45., true_scale_latitude=70.)

# This can be converted into a `proj4` string/dict compatible with GeoPandas
crs_proj4 = crs.proj4_init
###################### From Tedstone et al., 2022 #####################

plt.rcParams.update({'font.size': 15})
plt.rcParams["figure.figsize"] = (22,11.3)#from https://pythonguides.com/matplotlib-increase-plot-size/
fig = plt.figure()
gs = gridspec.GridSpec(39, 20)
gs.update(wspace=0.5)
#gs.update(wspace=0.001)
#projection set up from https://stackoverflow.com/questions/33942233/how-do-i-change-matplotlibs-subplot-projection-of-an-existing-axis
ax1 = plt.subplot(gs[0:34, 0:3],projection=crs)
ax_legend = plt.subplot(gs[34:39, 0:3])

ax2t = plt.subplot(gs[0:4, 4:20])
ax3t = plt.subplot(gs[7:11, 4:20])
ax4t = plt.subplot(gs[14:18, 4:20])
ax5t = plt.subplot(gs[21:25, 4:20])
ax6t = plt.subplot(gs[28:32, 4:20])
ax7t = plt.subplot(gs[35:39, 4:20])

#Draw plot of GrIS map
ax1.coastlines(edgecolor='black',linewidth=0.075)
#Display GrIS drainage bassins limits
GrIS_drainage_bassins.plot(ax=ax1,color='none', edgecolor='black',linewidth=0.075)
#Display region name
ax1.text(SW_rignotetal.centroid.x-80000,SW_rignotetal.centroid.y-80000,np.asarray(SW_rignotetal.SUBREGION1)[0])
ax1.text(CW_rignotetal.centroid.x-175000,CW_rignotetal.centroid.y-40000,np.asarray(CW_rignotetal.SUBREGION1)[0])
ax1.text(NW_rignotetal.centroid.x-50000,NW_rignotetal.centroid.y+20000,np.asarray(NW_rignotetal.SUBREGION1)[0])
ax1.text(NO_rignotetal.centroid.x-40000,NO_rignotetal.centroid.y-230000,np.asarray(NO_rignotetal.SUBREGION1)[0])

#Load 2010-2018 elevation dataset
f_20102018 = open(path+'df_20102018_with_elevation_Fig3_high_estimate_rignotetalregions_cleaned', "rb")
df_2010_2018_elevation = pickle.load(f_20102018)
f_20102018.close()

#Where ice content is higher than 16m, replace the ice content by 16!
df_2010_2018_elevation.loc[df_2010_2018_elevation['20m_ice_content_m']>16,'20m_ice_content_m']=16

#Plot data
plot_thickness_evolution(loc6,df_2010_2018_csv,df_2010_2018_elevation,GrIS_DEM,ax1,ax2t,custom_angle=-120,offset_x=7000,offset_y=-18000,casestudy_nb='a')

plot_thickness_evolution(loc8,df_2010_2018_csv,df_2010_2018_elevation,GrIS_DEM,ax1,ax3t,custom_angle=-90,offset_x=10000,offset_y=-5000,casestudy_nb='b')
#previousl b was loc8
plot_thickness_evolution(loc1,df_2010_2018_csv,df_2010_2018_elevation,GrIS_DEM,ax1,ax4t,custom_angle=-52,offset_x=10000,offset_y=1000,casestudy_nb='c')

plot_thickness_evolution(loc9,df_2010_2018_csv,df_2010_2018_elevation,GrIS_DEM,ax1,ax5t,custom_angle=-90,offset_x=10000,offset_y=-5000,casestudy_nb='d')

plot_thickness_evolution(loc3,df_2010_2018_csv,df_2010_2018_elevation,GrIS_DEM,ax1,ax6t,custom_angle=-90,offset_x=10000,offset_y=-5000,casestudy_nb='e')

plot_thickness_evolution(loc2,df_2010_2018_csv,df_2010_2018_elevation,GrIS_DEM,ax1,ax7t,custom_angle=-90,offset_x=10000,offset_y=-5000,casestudy_nb='f')

###################### From Tedstone et al., 2022 #####################
#from plot_map_decadal_change.py
# x0, x1, y0, y1
ax1.set_extent([-615320, -37273, -2909407, -1243910], crs=crs)
gl=ax1.gridlines(draw_labels=True, xlocs=[-40,-50,-60], ylocs=[65,70,75,80], x_inline=False, y_inline=False,linewidth=0.5,linestyle='dashed')
#Customize lat labels
gl.ylabels_right = False
gl.xlabels_bottom = False
ax1.axis('off')
###################### From Tedstone et al., 2022 #####################

#Display scalebar
scale_bar(ax1, (0.60, 0.05), 200, 3,-12)# axis, location (x,y), length, linewidth, rotation of text
#by measuring on the screen, the difference in precision between scalebar and length of transects is about ~200m

#Add shading where merging of diconsitnuous slabs
ax5t.axvspan(23100, 24600, facecolor='gray', alpha=0.3)
ax7t.axvspan(24600, 26700, facecolor='gray', alpha=0.3)

#Display distance as Elevation [m]
ax5t.set_ylabel('Column ice thickness [m]',fontsize=15)
ax7t.set_xlabel('Distance [km]',fontsize=15)
ax2t.xaxis.set_label_position("top")
ax2t.set_xlabel('Elevation [m]',fontsize=15,color='#8c510a')

#Custom legend myself
legend_elements = [Patch(facecolor='#fdd49e',label='2010'),
                   Patch(facecolor='#fc8d59',label='2011-2012'),
                   Patch(facecolor='#d7301f',label='2013-2014'),
                   Patch(facecolor='#7f0000',label='2017-2018')]

ax_legend.legend(handles=legend_elements,fontsize=15)
plt.legend()

#Get rid of axis in legend axis
ax_legend.axis('off')
plt.show()

pdb.set_trace()

#Save the figure
plt.savefig('C:/Users/jullienn/switchdrive/Private/research/RT1/figures/fig2/v10/fig2.png',dpi=300,bbox_inches='tight')

