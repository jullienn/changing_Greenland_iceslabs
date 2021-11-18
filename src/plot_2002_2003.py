# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 17:18:07 2021

@author: JullienN
"""
#################### Define functions for surface picking ###################

##############################################################################
############# Define kernel function for surface identification ##############
##############################################################################
#_gaussian function taken from IceBridgeGPR_Manager_v2.py
# Define a quick guassian function to scale the cutoff mask above
def _gaussian(x,mu,sigma):
    return np.exp(-np.power((x-mu)/sigma, 2.)/2.)
    
#This function have been taken from 'IceBridgeGPR_Manager_v2.py
def kernel_function(traces_input,suggested_pixel):
    
    traces = traces_input
    #Do not take the log10 of traces because 'data have been detrented in the log domain' according to John Paden's email, so I guess they are already log10!
    #traces = np.log10(traces)
    
    # We do not have the original indicies to use as a starter so we use our suggestion for surface picking start
    
    # 3) Perform surface pick crawling threshold behavior mask (assume a step-change analysis [goes from weak->strong at surface], and continuity of surface in original file.)
    # Create a step-change mask to optimze where the returns transition from "dark" to "bright"
    MASK_RADIUS = 50
    vertical_span_mask = np.empty([MASK_RADIUS*2,], dtype=np.float)
    vertical_span_mask[:MASK_RADIUS] = -1.0
    vertical_span_mask[MASK_RADIUS:] = +3.0
    
    vertical_span_mask = vertical_span_mask * _gaussian(np.arange(vertical_span_mask.shape[0]),mu=(MASK_RADIUS-5),sigma=(float(MASK_RADIUS)/3.0))
    
    # Expand the shape to handle array broadcasting below
    vertical_span_mask.shape = vertical_span_mask.shape[0], 1
    
    # This is the vertical window size of the extent of the search.  Should be bigger than any jump from one surface pixel to the next.
    MASK_SEARCH_RADIUS = 40
    
    improved_indices = np.zeros(traces.shape[1], dtype='int64')
    #pdb.set_trace()
    #traces.shape[1] indeed correspond to the horizontal distance
    
    # Start at the left with the hand-picked "suggested surface pick" in the ICEBRIDGE_SURFACE_PICK_SUGGESTIONS_FILE as starting point
    
    last_best_index = suggested_pixel
    
    #pdb.set_trace()
    # A template graph to use, just have to add in the center vertical index at each point and go from there.
    search_indices_template = np.sum(np.indices((vertical_span_mask.shape[0], 2*MASK_SEARCH_RADIUS)),axis=0) - MASK_SEARCH_RADIUS - MASK_RADIUS
    for i in range(traces.shape[1]):
        
        # Create an array of indices spanning the top-to-bottom of the MASK_SEARCH_RADIUS, and fanning out MASK_RADIUS above and below that point.
        search_indices = search_indices_template + last_best_index
        # Handle overflow indices if below zero or above max (shouldn't generally happen)... just assign to the top or bottom pixel
        search_indices[search_indices < 0] = 0
        search_indices[search_indices >= traces.shape[0]] = traces.shape[0]-1
        
        bestfit_sum = np.sum(traces[:,i][search_indices] * vertical_span_mask, axis=0)
            
        assert bestfit_sum.shape[0] == 2*MASK_SEARCH_RADIUS
        
        # Get the best fit (with the highest value from the transformation fit)
        last_best_index = search_indices[MASK_RADIUS,np.argmax(bestfit_sum)]
        improved_indices[i] = last_best_index
        
    #If there are pixels with particularly strong echo that are being erroneously
    #picked up as the surface, erase most the little "jump" artifacts in
    #the surface picker.
    improved_indices = _get_rid_of_false_surface_jumps(improved_indices)
        
    #I do not use any mask so I think I shouldn't need to use that:
    ###### Must re-expand the surface indices to account for masked values (filled w/ nan)
    ##### improved_indices_expanded = self._refill_array(improved_indices, surface_maskname)
    
    #pdb.set_trace()
    return improved_indices
##############################################################################
############# Define kernel function for surface identification ##############
##############################################################################

##############################################################################
################## Define functions for radar slice picking ##################
##############################################################################
def _radar_slice_indices_above_and_below(meters_cutoff_above, meters_cutoff_below,depths):
    #pdb.set_trace()

    delta_distance = np.mean(depths[1:] - depths[:-1])
    idx_above = int(np.round(float(meters_cutoff_above) / delta_distance))
    # Add one to the index below to include that last pixel when array-slicing
    idx_below = int(np.round(float(meters_cutoff_below) / delta_distance)) + 1
    
    return idx_above, idx_below
    
def _return_radar_slice_given_surface(traces,
                                      depths,
                                      surface_indices,
                                      meters_cutoff_above,
                                      meters_cutoff_below):
    '''From this radar track, return a "slice" of the image above and below the surface by
    (meters_cutoff_above, meters_cutoff_below), respectively.
    
    Return value:
    A ((idx_below+idx_above), numtraces]-sized array of trace sample values.
    '''
    idx_above, idx_below = _radar_slice_indices_above_and_below(meters_cutoff_above, meters_cutoff_below,depths)
    
    output_traces = np.empty((idx_above + idx_below, traces.shape[1]), dtype=traces.dtype)
    bottom_indices = np.zeros(shape=(1,traces.shape[1]))
    
    for i,s in enumerate(surface_indices):
        try:
            output_traces[:,i] = traces[(s-idx_above):(s+idx_below), i]
            bottom_indices[0,i]=(s+idx_below)
        except ValueError:
            # If the surf_i is too close to one end of the array or the other, it extends beyond the edge of the array and breaks.
            if s < idx_above:
                start, end = None, idx_above+idx_below
            elif s > (traces.shape[0] - idx_below):
                start, end = traces.shape[0] - (idx_above + idx_below), None
            else:
                # SHouldn't get here.
                print(i, s, traces.shape)
                assert False
            output_traces[:,i] = traces[start:end, i]
            bottom_indices[0,i]=end 
    return output_traces, bottom_indices


def _get_rid_of_false_surface_jumps(surface_indices):
    '''Some of the 2011 files especially, have strong echos that are errantly being picked up as the surface.  Find these big "jumps", and get rid of them.  Use the suggested surface instead.'''
    improved_surface = surface_indices.copy()
    
    jumps = improved_surface[1:] - improved_surface[:-1]
    # Substitute any large jumps with brightest pixel in a window of original surface.  Do this until large jumps either go away or have all been corrected to original surface.
    for i in range(len(jumps)):
        
        # Slope windowsize = number of pixels we use to average the previous slope.
        slope_windowsize = 10
        if i < slope_windowsize:
            continue
        mean_slope = np.mean(np.array(jumps[i-slope_windowsize:i], dtype=np.float))
    
        # Find the difference of this slope from the last five stops
        difference_from_mean_slope = jumps[i] - mean_slope
        # Ignore if it's jumped less than 3 from the mean recent slope, or less than 50% greater than the mean slope at this time.
        if (difference_from_mean_slope < 5) or (difference_from_mean_slope < (1.5*mean_slope)):
            continue
    
        # tune settings
        jump_lookahead = 20 # Number of pixels to look ahead to see if we can find a matching down-jump
        if i+jump_lookahead > len(jumps):
            jump_lookahead = len(jumps) - i
    
        # This is how close the surface on the "other side" of the jump must be to the original slope to be considered for it.
        jump_magnitude_threshold = 1.10
    
        # See if we can find a point in the near future that would approximate the current slope.
        slopes_ahead = np.cumsum(jumps[i:i+jump_lookahead]) / np.arange(1,jump_lookahead+1)
        opposite_match = np.argmax(slopes_ahead <= (mean_slope * jump_magnitude_threshold))
        
        if opposite_match > 0:
            # We found a match, onward!
            opposite_match_index = i + opposite_match
            for j in range(i+1,opposite_match_index+1):
                improved_surface[j] = np.round(improved_surface[i] + float(improved_surface[opposite_match_index+1] - improved_surface[i])*(j-i)/(opposite_match_index+1-i))    
            # now recompute jumps
            jumps = improved_surface[1:] - improved_surface[:-1]
            continue
    
        # IF THE ABOVE DIDN'T WORK, TRY THE 'JUMP' TECHNIQUE, SEEING WHETHER AN ANOMALOUS 'JUMP' IS COUNTERBALANCED BY AN
        # OPPOSITE AND (APPROXIMATELY) EQUAL JUMP IN THE OPPOSITE DIRECTION.
        # Don't worry about any trends less than 12 pixels.  Hills do that.
        jump = jumps[i]
        if abs(jump) < 5:
            continue
    
        # tune settings
        jump_lookahead = 50 # Number of pixels to look ahead to see if we can find a matching down-jump
        jump_magnitude_threshold = 0.50 # What fraction of the original jump the new jump has to be (in the opposite direction) to qualify.
    
        # see if we can find a jump in the near-future that crosses this threshold in the other direction.  If so, we've found our counter-part
        if jump < 0:
            opposite_jump_index = np.argmax((jumps[i:i+jump_lookahead]) > (-jump*jump_magnitude_threshold))
        elif jump > 0:
            opposite_jump_index = np.argmax((jumps[i:i+jump_lookahead]) < (-jump*jump_magnitude_threshold))
    
        if opposite_jump_index > 0:
            opposite_jump_index += i
        else: # If we didn't find a partner opposite offset, skip and move along.
            continue
    
        # Linearly interpolate, get to the closest pixel
        try:
            for j in range(i+1,opposite_jump_index+1):
                improved_surface[j] = np.round(improved_surface[i] + float(improved_surface[opposite_jump_index+1] - improved_surface[i])*(j-i)/(opposite_jump_index+1-i))
        except IndexError:
            print("i", i, "j", j, "opposite_jump_index", opposite_jump_index, improved_surface.shape, jumps.shape)
            # Break the program here.
            100/0
    
        # now recompute jumps
        jumps = improved_surface[1:] - improved_surface[:-1]
        continue
    return improved_surface
##############################################################################
################## Define functions for radar slice picking ##################
##############################################################################
    
#################### Define functions for surface picking ###################

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


def plot_radar_slice(ax_map,ax_plot,ax_elevation,ax_nb,path_radar_slice,lines,folder_year,folder_day,indiv_file,technique,xls_icelenses,trafic_light,elevation_dictionnary):
    
    #Define the uppermost and lowermost limits
    meters_cutoff_above=0
    meters_cutoff_below=30
    
    dt = 2.034489716724874e-09 #Timestep for 2002/2003 traces
    t0 = 0; # Unknown so set to zero
    #Compute the speed (Modified Robin speed):
    # self.C / (1.0 + (coefficient*density_kg_m3/1000.0))
    v= 299792458 / (1.0 + (0.734*0.873/1000.0))
    
    if (folder_day=='jun04'):
        #Open the file and read it
        fdata= scipy.io.loadmat(path_radar_slice)
        #Select radar echogram, lat and lon
        radar_echo=fdata['data']
        lat=fdata['latitude']
        lon=fdata['longitude']
    else:
        #Open the file and read it
        f_agg = open(path_radar_slice, "rb")
        radar_data = pickle.load(f_agg)
        f_agg.close()
                                                
        #Select radar echogram, lat and lon
        radar_echo=radar_data['radar_echogram']
        
        latlontime=radar_data['latlontime']
        lat=latlontime['lat_gps']
        lon=latlontime['lon_gps']
        
        #Remove zeros in lat/lon
        lat.replace(0, np.nan, inplace=True)
        lon.replace(0, np.nan, inplace=True)
    
    #pdb.set_trace()
    #Transform the longitudes. The longitudes are ~46 whereas they should be ~-46! So add a '-' in front of lon
    lon=-lon
                        
    #Transform the coordinated from WGS84 to EPSG:3413
    #Example from: https://pyproj4.github.io/pyproj/stable/examples.html
    transformer = Transformer.from_crs("EPSG:4326", "EPSG:3413", always_xy=True)
    points=transformer.transform(np.array(lon),np.array(lat))
    
    #Reset the lat_3413 and lon_3413 to empty vectors.
    lon_3413=[]
    lat_3413=[]
    
    lon_3413=points[0]
    lat_3413=points[1]
    
    #pdb.set_trace()
    #Define the dates that need reversed display
    list_reverse_agg=['may12_03_36_aggregated','may14_03_51_aggregated',
                       'may13_03_29_aggregated','may30_02_51_aggregated',
                       'may24_02_25_aggregated','may15_03_37_aggregated',
                       'may11_03_29_aggregated']
        
    list_reverse_mat=['jun04_02proc_52.mat','jun04_02proc_53.mat']
    
    #Display on the map where is this track
    ax_map.scatter(lon_3413, lat_3413,s=1,facecolors='black', edgecolors='none')
    
    #1. Compute the vertical resolution
    #a. Time computation according to John Paden's email.
    Nt = radar_echo.shape[0]
    Time = t0 + dt*np.arange(1,Nt+1)
    #b. Calculate the depth:
    #self.SAMPLE_DEPTHS = self.radar_speed_m_s * self.SAMPLE_TIMES / 2.0
    depths = v * Time / 2.0
    
    # Load the suggested pixel for the specific date
    for date_pix in lines:
        if (folder_day=='jun04'):
            if (date_pix.partition(" ")[0]==str(indiv_file.replace(".mat",""))):
                suggested_pixel=int(date_pix.partition(" ")[2])
                #If it has found its suggested pixel, leave the loop
                continue   
        else:
            if (date_pix.partition(" ")[0]==str(indiv_file.replace("_aggregated",""))):
                suggested_pixel=int(date_pix.partition(" ")[2])
                #If it has found its suggested pixel, leave the loop
                continue

    surface_indices=kernel_function(radar_echo, suggested_pixel)
    
    #Get our slice (30 meters as currently set)
    radar_slice, bottom_indices = _return_radar_slice_given_surface(radar_echo,
                                                                    depths,
                                                                    surface_indices,
                                                                    meters_cutoff_above=meters_cutoff_above,
                                                                    meters_cutoff_below=meters_cutoff_below)
    
    #The range have already been computed, plot the data:
    if (technique=='perc_25_75'):
        if (folder_year=='2002'):
            perc_lower_end=-0.08318485583215623
            perc_upper_end=0.09414986209628376
        elif (folder_year=='2003'):
            perc_lower_end=-0.08488332785270308
            perc_upper_end=0.09654050592743407
    elif (technique=='perc_5_95'):
        if (folder_year=='2002'):
            perc_lower_end=-0.2870889087496134
            perc_upper_end=0.3138722799744009
        elif (folder_year=='2003'):
            perc_lower_end=-0.31927843730229416
            perc_upper_end=0.3682849426401127
    elif (technique=='perc_05_995'):
        if (folder_year=='2002'):
            perc_lower_end=-2.1488917418616134
            perc_upper_end=2.650167679823621
        elif (folder_year=='2003'):
            perc_lower_end=-1.661495950494564
            perc_upper_end=1.9431298622848088
    elif (technique=='perc_2p5_97p5'):
        if (folder_year=='2002'):
            perc_lower_end=-0.5709792307554173
            perc_upper_end=0.7082634842114803
        elif (folder_year=='2003'):
            perc_lower_end=-0.6061610403154447
            perc_upper_end=0.7572821079440079      
    
    #Generate the pick for vertical distance display
    ticks_yplot=np.arange(0,radar_slice.shape[0],20)
    
    #Plot the radar slice
    cb2=ax_plot.pcolor(radar_slice,cmap=plt.get_cmap('gray'))#,norm=divnorm)
    ax_plot.set_ylim(0,radar_slice.shape[0])
    ax_plot.invert_yaxis() #Invert the y axis = avoid using flipud.
    ax_plot.set_aspect('equal') # X scale matches Y scale
    #ax_plot.set_xlabel('Horizontal distance')
    
    #Colorbar custom
    cb2.set_clim(perc_lower_end,perc_upper_end)
    #cbar2=fig.colorbar(cb2, ax=[ax_plot], location='left')
    #cbar2.set_label('Signal strength')
    
    #Set the y ticks
    ax_plot.set_yticks(ticks_yplot) 
    ax_plot.set_yticklabels(np.round(depths[ticks_yplot]))
    
    #Set the x ticks
    #remove xtick
    #ax_plot.set_xticks([])
    
    #Distance from start of the trace
    if (ax_nb==2):
        ax_plot.set_title('Ablation zone',fontsize=10)
        ax_plot.set_ylabel('Depth [m]')
        ax_plot.set_xlabel('Distance [km]')
        
        letter_elev='a'
        #Display the correspondance between radar slice and radar location on the map
        ax_plot.text(0,99,letter_elev,color='black',fontsize=20)
        #Display on the map the letter corresponding to the radar slice
        ax_map.text(np.nanmedian(lon_3413)-10000,np.nanmedian(lat_3413),'a',color='black',fontsize=15)
    
    elif (ax_nb==3):
        ax_plot.set_title('Percolation zone - ice lenses',fontsize=10)
        
        letter_elev='b'
        #Display the correspondance between radar slice and radar location on the map
        ax_plot.text(1000,99,letter_elev,color='black',fontsize=20)
        #Display on the map the letter corresponding to the radar slice
        ax_map.text(np.nanmedian(lon_3413)-5000,np.nanmedian(lat_3413)+5000,'b',color='black',fontsize=15)
    
    elif (ax_nb==4):
        ax_plot.set_title('Percolation zone - ice slabs',fontsize=10)
        
        letter_elev='c'
        #Display the correspondance between radar slice and radar location on the map
        ax_plot.text(1000,99,letter_elev,color='black',fontsize=20)
        #Display on the map the letter corresponding to the radar slice
        ax_map.text(np.nanmedian(lon_3413),np.nanmedian(lat_3413)+4000,'c',color='black',fontsize=15)
    
    elif (ax_nb==5):
        ax_plot.set_title('Dry snow zone',fontsize=10)
        
        letter_elev='d'
        #Display the correspondance between radar slice and radar location on the map
        ax_plot.text(1000,99,letter_elev,color='black',fontsize=20)
        #Display on the map the letter corresponding to the radar slice
        ax_map.text(np.nanmedian(lon_3413)-11000,np.nanmedian(lat_3413),'d',color='black',fontsize=15)
        
    #Display the ice lenses identification:
    #pdb.set_trace()

    if (indiv_file in list(xls_icelenses.keys())):
        print(indiv_file+' hold ice lens!')
        #This file have ice lenses in it: read the data:
        df_temp=xls_icelenses[indiv_file]
        df_colnames = list(df_temp.keys())
        x_loc=[]
        
        #Trafic light information
        df_trafic_light=trafic_light[indiv_file]
        df_colnames_trafic_light = list(df_trafic_light.keys())
        
        for i in range (0,int(len(df_colnames)),2):
            x_vect=df_temp[df_colnames[i]]
            y_vect=df_temp[df_colnames[i+1]]
            #Load trafic light color
            trafic_light_indiv_color=df_colnames_trafic_light[i]
            #Define the color in which to display the ice lens
            if (trafic_light_indiv_color[0:3]=='gre'):
                color_to_display='#00441b'
            elif (trafic_light_indiv_color[0:3]=='ora'):
                color_to_display='#fed976'
            elif (trafic_light_indiv_color[0:3]=='red'):
                color_to_display='#c9662c'
            elif (trafic_light_indiv_color[0:3]=='pur'):
                color_to_display='purple'
            else:
                print('The color is not known!')
            
            #Display ice lens
            ax_plot.plot(x_vect,y_vect,color=color_to_display,linestyle='dashed',linewidth=0.5)
    
    #Load the elevation profile
    elevation_vector=elevation_dictionnary[folder_year][folder_day][indiv_file]
    
    #Transpose if june 04
    if (folder_day=='jun04'):
        lat_3413=np.transpose(lat_3413)
        lon_3413=np.transpose(lon_3413)
    
    #Calculate the distances (in m)
    distances=compute_distances(lon_3413,lat_3413)
    
    #Convert distances from m to km
    distances=distances/1000
        
    #Order the radar track from down to up if needed      
    if (indiv_file in list(list_reverse_agg)):
        ax_plot.set_xlim(radar_slice.shape[1],0)
        #plot the reversed elevation profile
        ax_elevation.plot(np.arange(0,len(elevation_vector)),np.flipud(elevation_vector),color='black')
        #Reverse the distances vector:
        distances=np.flipud(distances)
        
    elif (indiv_file in list(list_reverse_mat)):
        ax_plot.set_xlim(radar_slice.shape[1],0)
        #plot the the reversed elevation profile
        ax_elevation.plot(np.arange(0,len(elevation_vector)),np.flipud(elevation_vector),color='black')
        #Reverse the distances vector:
        distances=np.flipud(distances)
    else:
        #plot the elevation profile
        ax_elevation.plot(np.arange(0,len(elevation_vector)),elevation_vector,color='black')
    
    #Generate the pick for horizontal distance display
    ticks_xplot=np.arange(0,distances.shape[0]+1,100)
    #Plot also the last index
    ticks_xplot[-1]=distances.shape[0]-1
    #Set x ticks
    ax_plot.set_xticks(ticks_xplot) 
    #Display the distances from the origin as being the x label
    ax_plot.set_xticklabels(np.round(distances[ticks_xplot]))
    
    #Define the elevation to plot the letter
    if (ax_nb==2):
        #Display which elevation line correspond to which trace
        ax_elevation.text(950,elevation_vector[950]-150,letter_elev)
    elif (ax_nb==3):
        #Display which elevation line correspond to which trace
        ax_elevation.text(950,elevation_vector[50]-150,letter_elev)
    elif (ax_nb==4):
        #Display which elevation line correspond to which trace
        ax_elevation.text(950,elevation_vector[50]-150,letter_elev)
    elif (ax_nb==5):
        #Display which elevation line correspond to which trace
        ax_elevation.text(950,elevation_vector[50]-150,letter_elev)

    return

def compute_distances(eastings,northings):
    #This part of code is from MacFerrin et al., 2019
    '''Compute the distance (in m here, not km as written originally) of the traces in the file.'''
    # C = sqrt(A^2  + B^2)
    distances = np.power(np.power((eastings[1:] - eastings[:-1]),2) + np.power((northings[1:] - northings[:-1]),2), 0.5)
    #Calculate the cumsum of the distances
    cumsum_distances=np.nancumsum(distances)
    #Seeting the first value of the cumsum to be zero as it is the origin
    return_cumsum_distances=np.zeros(eastings.shape[0])
    return_cumsum_distances[1:eastings.shape[0]]=cumsum_distances
    
    return return_cumsum_distances

def plot_radar_slice_with_thickness(ax_map,ax_elevation,ax_plot,path_radar_slice,lines,folder_year,folder_day,indiv_file,technique,xls_icelenses,trafic_light,elevation_dictionnary,icelens_information):
    
    #Define color code for trafic light plotting in elevation plot
    cmap_elevation = ListedColormap(['#c9662c', '#fed976', '#238b45', 'purple'])
    norm_elevation = BoundaryNorm([-1.5, -0.5, 0.5, 1.5, 2.5], cmap_elevation.N)
    
    #Define the uppermost and lowermost limits
    meters_cutoff_above=0
    meters_cutoff_below=30
    
    dt = 2.034489716724874e-09 #Timestep for 2002/2003 traces
    t0 = 0; # Unknown so set to zero
    #Compute the speed (Modified Robin speed):
    # self.C / (1.0 + (coefficient*density_kg_m3/1000.0))
    v= 299792458 / (1.0 + (0.734*0.873/1000.0))
    
    if (folder_day=='jun04'):
        #Open the file and read it
        fdata= scipy.io.loadmat(path_radar_slice)
        #Select radar echogram, lat and lon
        radar_echo=fdata['data']
        lat=fdata['latitude']
        lon=fdata['longitude']
    else:
        #Open the file and read it
        f_agg = open(path_radar_slice, "rb")
        radar_data = pickle.load(f_agg)
        f_agg.close()
                                                
        #Select radar echogram, lat and lon
        radar_echo=radar_data['radar_echogram']
        
        latlontime=radar_data['latlontime']
        lat=latlontime['lat_gps']
        lon=latlontime['lon_gps']
        
        #Remove zeros in lat/lon
        lat.replace(0, np.nan, inplace=True)
        lon.replace(0, np.nan, inplace=True)
    
    #Transform the longitudes. The longitudes are ~46 whereas they should be ~-46! So add a '-' in front of lon
    lon=-lon
                        
    #Transform the coordinated from WGS84 to EPSG:3413
    #Example from: https://pyproj4.github.io/pyproj/stable/examples.html
    transformer = Transformer.from_crs("EPSG:4326", "EPSG:3413", always_xy=True)
    points=transformer.transform(np.array(lon),np.array(lat))
    
    #Reset the lat_3413 and lon_3413 to empty vectors.
    lon_3413=[]
    lat_3413=[]
    
    lon_3413=points[0]
    lat_3413=points[1]
    
    #Transpose coordinates if june 04
    if (folder_day=='jun04'):
        lat_3413=np.transpose(lat_3413)
        lon_3413=np.transpose(lon_3413)
    
    #pdb.set_trace()
    #Define the dates that need reversed display
    list_reverse_agg=[]
    list_reverse_mat=[]
    #list_reverse_agg=['may12_03_36_aggregated','may14_03_51_aggregated',
    #                   'may13_03_29_aggregated','may30_02_51_aggregated',
    #                   'may24_02_25_aggregated','may15_03_37_aggregated',
    #                   'may11_03_29_aggregated']
        
    #list_reverse_mat=['jun04_02proc_52.mat','jun04_02proc_53.mat']
    
    #Display on the map where is this track
    ax_map.scatter(lon_3413, lat_3413,s=5,facecolors='black', edgecolors='none')
    
    #pdb.set_trace()
    
    #Load deepest ice lenses information
    deepest_icelenses=icelens_information[indiv_file]
    #Retrieve the index where deepest data are present
    index_deepest_data_present=~(np.isnan(np.asarray(deepest_icelenses['x'])))
    #Display the depth of the deepest ice lens in the map
    cb_depth=ax_map.scatter(lon_3413[index_deepest_data_present], lat_3413[index_deepest_data_present],c=np.asarray(deepest_icelenses['deepest_depth'])[index_deepest_data_present],cmap=discrete_cmap(10,'RdYlGn_r'),s=5, edgecolors='none')
    cbar_depth=fig.colorbar(cb_depth, ax=[ax_map], location='right')
    cbar_depth.set_label('Ice lens maximum depth [m]')
    cb_depth.set_clim(0,20)
    
    #Display the start of the track
    ax_map.scatter(lon_3413[0],lat_3413[0],c='m',s=5, edgecolors='none')
    
    #Zoom on the trace on the map plot
    ax_map.set_xlim(np.nanmedian(lon_3413)-75000, np.nanmedian(lon_3413)+75000)
    ax_map.set_ylim(np.nanmedian(lat_3413)-75000, np.nanmedian(lat_3413)+75000)
    
    #1. Compute the vertical resolution
    #a. Time computation according to John Paden's email.
    Nt = radar_echo.shape[0]
    Time = t0 + dt*np.arange(1,Nt+1)
    #b. Calculate the depth:
    #self.SAMPLE_DEPTHS = self.radar_speed_m_s * self.SAMPLE_TIMES / 2.0
    depths = v * Time / 2.0
    
    # Load the suggested pixel for the specific date
    for date_pix in lines:
        if (folder_day=='jun04'):
            if (date_pix.partition(" ")[0]==str(indiv_file.replace(".mat",""))):
                suggested_pixel=int(date_pix.partition(" ")[2])
                #If it has found its suggested pixel, leave the loop
                continue   
        else:
            if (date_pix.partition(" ")[0]==str(indiv_file.replace("_aggregated",""))):
                suggested_pixel=int(date_pix.partition(" ")[2])
                #If it has found its suggested pixel, leave the loop
                continue

    surface_indices=kernel_function(radar_echo, suggested_pixel)
    
    #Get our slice (30 meters as currently set)
    radar_slice, bottom_indices = _return_radar_slice_given_surface(radar_echo,
                                                                    depths,
                                                                    surface_indices,
                                                                    meters_cutoff_above=meters_cutoff_above,
                                                                    meters_cutoff_below=meters_cutoff_below)
    
    #The range have already been computed, plot the data:
    if (technique=='perc_25_75'):
        if (folder_year=='2002'):
            perc_lower_end=-0.08318485583215623
            perc_upper_end=0.09414986209628376
        elif (folder_year=='2003'):
            perc_lower_end=-0.08488332785270308
            perc_upper_end=0.09654050592743407
    elif (technique=='perc_5_95'):
        if (folder_year=='2002'):
            perc_lower_end=-0.2870889087496134
            perc_upper_end=0.3138722799744009
        elif (folder_year=='2003'):
            perc_lower_end=-0.31927843730229416
            perc_upper_end=0.3682849426401127
    elif (technique=='perc_05_995'):
        if (folder_year=='2002'):
            perc_lower_end=-2.1488917418616134
            perc_upper_end=2.650167679823621
        elif (folder_year=='2003'):
            perc_lower_end=-1.661495950494564
            perc_upper_end=1.9431298622848088
    elif (technique=='perc_2p5_97p5'):
        if (folder_year=='2002'):
            perc_lower_end=-0.5709792307554173
            perc_upper_end=0.7082634842114803
        elif (folder_year=='2003'):
            perc_lower_end=-0.6061610403154447
            perc_upper_end=0.7572821079440079      
    
    #Generate the pick for vertical distance display
    ticks_yplot=np.arange(0,radar_slice.shape[0],20)
    
    #Plot the radar slice
    cb2=ax_plot.pcolor(radar_slice,cmap=plt.get_cmap('gray'))#,norm=divnorm)
    ax_plot.set_ylim(0,radar_slice.shape[0])
    ax_plot.invert_yaxis() #Invert the y axis = avoid using flipud.
    ax_plot.set_aspect('equal') # X scale matches Y scale
    #ax_plot.set_xlabel('Horizontal distance')
    
    #Colorbar custom
    cb2.set_clim(perc_lower_end,perc_upper_end)
    #cbar2=fig.colorbar(cb2, ax=[ax_plot], location='left')
    #cbar2.set_label('Signal strength')
    
    #Set the y ticks
    ax_plot.set_yticks(ticks_yplot) 
    ax_plot.set_yticklabels(np.round(depths[ticks_yplot]))
    
    #Set the x ticks
    #remove xtick
    #ax_plot.set_xticks([])
    
    #Distance from start of the trace
    ax_plot.set_title('Radar slice')
    ax_plot.set_ylabel('Depth [m]')
    ax_plot.set_xlabel('Distance [km]')
    
    #Display the ice lenses identification:
    #pdb.set_trace()

    if (indiv_file in list(xls_icelenses.keys())):
        print(indiv_file+' hold ice lens!')
        #This file have ice lenses in it: read the data:
        df_temp=xls_icelenses[indiv_file]
        df_colnames = list(df_temp.keys())
        x_loc=[]
        
        #Trafic light information
        df_trafic_light=trafic_light[indiv_file]
        df_colnames_trafic_light = list(df_trafic_light.keys())
        
        for i in range (0,int(len(df_colnames)),2):
            x_vect=df_temp[df_colnames[i]]
            y_vect=df_temp[df_colnames[i+1]]
            #Load trafic light color
            trafic_light_indiv_color=df_colnames_trafic_light[i]
            #Define the color in which to display the ice lens
            if (trafic_light_indiv_color[0:3]=='gre'):
                color_to_display='#00441b'
            elif (trafic_light_indiv_color[0:3]=='ora'):
                color_to_display='#fed976'
            elif (trafic_light_indiv_color[0:3]=='red'):
                color_to_display='#c9662c'
            elif (trafic_light_indiv_color[0:3]=='pur'):
                color_to_display='purple'
            else:
                print('The color is not known!')
            
            #Display ice lens color on radar slice
            ax_plot.plot(x_vect,y_vect,color=color_to_display,linestyle='dashed',linewidth=0.5)
            #Display ice lens color on elevation plot
    
    #pdb.set_trace()
    #Display the deepest identifies icelens
    ax_plot.scatter(np.asarray(deepest_icelenses['x']),np.asarray(deepest_icelenses['deepest_depth_index']),color='red',s=1)

    #Load the elevation profile
    elevation_vector=elevation_dictionnary[folder_year][folder_day][indiv_file]
    
    #Store the color code
    color_code_all=np.asarray(deepest_icelenses['deepest_depth_color'])
    #Create the vector for color code display
    elevation_color=np.zeros(elevation_vector.shape[0])
    #Fill in the elevation_color vector with NaNs
    elevation_color[elevation_color==0]=np.nan
    #Fill in the elevation_color vector with the color code vector
    elevation_color[~np.isnan(color_code_all)]=elevation_vector[~np.isnan(color_code_all)]
    
    #Plot the elevation profile
    ax_elevation.plot(np.arange(0,len(elevation_vector)),elevation_vector,color='black')
    #Plot the elevation profile with the color code were ice lenses
    #ax_elevation.plot(np.arange(0,len(elevation_vector)),elevation_color,cmap=cmap_elevation,norm=norm_elevation)

    ax_elevation.scatter(np.arange(0,len(elevation_vector)),elevation_color,c=color_code_all,cmap=cmap_elevation,norm=norm_elevation,s=10)
    ax_elevation.set_title('Trace elevation profile')
    #pdb.set_trace()

    #Calculate the distances (in m)
    distances=compute_distances(lon_3413,lat_3413)
    
    #Convert distances from m to km
    distances=distances/1000
        
    #Order the radar track from down to up if needed      
    if (indiv_file in list(list_reverse_agg)):
        #Reverse xlim in radar slice plot
        ax_plot.set_xlim(radar_slice.shape[1],0)
        #Reverse xlim in elevation plot
        ax_elevation.set_xlim(radar_slice.shape[1],0)
        #Reverse the distances vector:
        distances=np.flipud(distances)
        
    elif (indiv_file in list(list_reverse_mat)):
        #Reverse xlim in radar slice plot
        ax_plot.set_xlim(radar_slice.shape[1],0)
        #Reverse xlim in elevation plot
        ax_elevation.set_xlim(radar_slice.shape[1],0)
        #Reverse the distances vector:
        distances=np.flipud(distances)
    
    #Generate the pick for horizontal distance display
    ticks_xplot=np.arange(0,distances.shape[0]+1,100)
    #Plot also the last index
    ticks_xplot[-1]=distances.shape[0]-1
    #Set x ticks
    ax_plot.set_xticks(ticks_xplot) 
    #Display the distances from the origin as being the x label
    ax_plot.set_xticklabels(np.round(distances[ticks_xplot]))
    
    #Change xticks for elevation display
    #Set xlim
    ax_elevation.set_xlim(0,radar_slice.shape[1])
    #Set x ticks
    ax_elevation.set_xticks(ticks_xplot) 
    #Display the distances from the origin as being the x label
    ax_elevation.set_xticklabels(np.round(distances[ticks_xplot]))
    #Display the x and ylabel
    ax_elevation.set_xlabel('Distance [km]')
    ax_elevation.set_ylabel('Elevation [m]')
    
    #pdb.set_trace()
    return

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

#Import packages
import rasterio
from rasterio.plot import show
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pandas as pd
from os import listdir
from os.path import isfile, join
import pickle
from pysheds.grid import Grid
import pdb
import numpy as np
from pyproj import Transformer
import matplotlib.gridspec as gridspec
import scipy.io
from osgeo import gdal

from matplotlib.colors import ListedColormap, BoundaryNorm

technique='perc_2p5_97p5'
making_down_to_up='FALSE'
plt.rcParams.update({'font.size': 7})
########################## Load GrIS elevation ##########################
#Open the DEM
grid = Grid.from_raster("C:/Users/jullienn/switchdrive/Private/research/backup_Aglaja/working_environment/greenland_topo_data/elevations/greenland_dem_mosaic_100m_v3.0.tif",data_name='dem')
#Minnor slicing on borders to enhance colorbars
elevDem=grid.dem[:-1,:-1]              
#Scale the colormap
divnorm = mcolors.DivergingNorm(vmin=0, vcenter=1250, vmax=2500)
########################## Load GrIS elevation ##########################

################# Load 2002-2003 flightlines coordinates ################
path_data='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/icelens_identification'

#Open the file and read it
f_flightlines = open(path_data+'/metadata_coord_2002_2003', "rb")
all_2002_3_flightlines = pickle.load(f_flightlines)
f_flightlines.close()
################# Load 2002-2003 flightlines coordinates ################

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

################### Load 2010-2014 ice slabs location ##################
#Load the data
filename_MacFerrin= 'C:/Users/jullienn/switchdrive/Private/research/backup_Aglaja/working_environment/iceslabs_MacFerrin/icelens_identification/indiv_traces_icelenses/MacFerrin_etal2019_iceslabs.xlsx'
#Read MacFerrin data thanks to https://stackoverflow.com/questions/65254535/xlrd-biffh-xlrderror-excel-xlsx-file-not-supported
df_MacFerrin = pd.read_excel(filename_MacFerrin,engine='openpyxl')

#Transform the coordinated from WGS84 to EPSG:3413
#Example from: https://pyproj4.github.io/pyproj/stable/examples.html
transformer = Transformer.from_crs("EPSG:4326", "EPSG:3413", always_xy=True)
points=transformer.transform(np.array(df_MacFerrin.lon),np.array(df_MacFerrin.lat))

lon_3413_MacFerrin=points[0]
lat_3413_MacFerrin=points[1]
################### Load 2010-2014 ice slabs location ##################

################### Load 2017-2018 ice slabs location ##################
#Load the data
filename_Jullien= 'C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2010_2018/excel_spatial_aggreation_and_other/final_excel/prob00/Ice_Layer_Output_Thicknesses_2010_2018_jullienetal2021_prob00.csv'
#Read Jullien data thanks to https://stackoverflow.com/questions/65254535/xlrd-biffh-xlrderror-excel-xlsx-file-not-supported
df_jullien = pd.read_csv(filename_Jullien,delimiter=',',decimal='.')

#Keep only 2017-2018 data
year=df_jullien.Track_name.str[:4]
year_int=year.astype(str).astype(int)

index_2017=(df_jullien.Track_name.str[:4] == '2017')
index_2018=(df_jullien.Track_name.str[:4] == '2018')

df_2017_2018=df_jullien[index_2017]
df_2017_2018=df_2017_2018.append(df_jullien[index_2018])

#Transform the coordinated from WGS84 to EPSG:3413
#Example from: https://pyproj4.github.io/pyproj/stable/examples.html
transformer = Transformer.from_crs("EPSG:4326", "EPSG:3413", always_xy=True)
points=transformer.transform(np.array(df_2017_2018.lon),np.array(df_2017_2018.lat))

lon_3413_2017_2018=points[0]
lat_3413_2017_2018=points[1]
################### Load 2017-2018 ice slabs location ##################

################################### Plot ##################################
#Prepare plot
fig = plt.figure(figsize=(19,10))
fig.suptitle('2002-2003 ice lenses and ice slabs mapping SW Greenland')
gs = gridspec.GridSpec(10, 20)
gs.update(wspace=0.1)
gs.update(wspace=0.001)
ax1 = plt.subplot(gs[0:10, 10:20])
ax2 = plt.subplot(gs[0:2, 0:10])
ax3 = plt.subplot(gs[2:4, 0:10])
ax4 = plt.subplot(gs[4:6, 0:10])
ax5 = plt.subplot(gs[6:8, 0:10])
ax6 = plt.subplot(gs[8:10, 0:10])

#Display elevation
cb1=ax1.imshow(elevDem, extent=grid.extent,cmap=discrete_cmap(10,'cubehelix_r'),alpha=0.5,norm=divnorm)
#cbar1=fig.colorbar(cb1, ax=[ax1], location='left')
#ax1.grid()
ax1.set_title('Ice lenses and slabs location')

#Plot all the 2010-2014 icelenses
ax1.scatter(lon_3413_MacFerrin, lat_3413_MacFerrin,s=1,facecolors='cornflowerblue', edgecolors='none')
#ax1.scatter(lon_3413_MacFerrin, lat_3413_MacFerrin,color='red',marker='o',alpha=0.2)

#Plot all the 2002-2003 flightlines
ax1.scatter(lon_all, lat_all,s=1,facecolors='lightgrey', edgecolors='none',alpha=0.1)
################################### Plot ##################################

#Open several files to display on top of the map

#display may12_03_36, may12_02_1, may09_03_1, jun04_02_53, may14_03_7
#dry snow zone:may18_02_16 or 18_02_12
# east greenland: jun04_proc02_7.mat

#Open, read and close the file of suggested surface picks
f = open('C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/Exclusion_folder/txt/SURFACE_STARTING_PICKS_Suggestions_2002_2003.txt','r')
lines = [line.strip() for line in f.readlines() if len(line.strip()) > 0]
f.close()

#Open and read the excel file having the ice lenses/slabs in it
filename_icelenses='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/icelens_identification/indiv_traces_icelenses/icelenses_22022020.xls'
xls_icelenses = pd.read_excel(filename_icelenses, sheet_name=None,header=2)
trafic_light=pd.read_excel(filename_icelenses, sheet_name=None,header=1)

#Specify the general path name
path_radar_data='C:/Users/jullienn/switchdrive/Private/research/backup_Aglaja/working_environment/iceslabs_MacFerrin/data'

#Plot date 1
folder_year='2003'
folder_day='may11'
indiv_file='may11_03_1_aggregated' #From down to up: OK!
ax_nb=2
path_radar_slice=path_radar_data+'/'+folder_year+'/'+folder_day+'/'+indiv_file
plot_radar_slice(ax1,ax2,ax6,ax_nb,path_radar_slice,lines,folder_year,folder_day,indiv_file,technique,xls_icelenses,trafic_light,elevation_dictionnary)

#pdb.set_trace()
#Plot date 2
folder_year='2002'
folder_day='jun04'
indiv_file='jun04_02proc_53.mat' #From up to down: need reversing! Already done, OK!
ax_nb=3
path_radar_slice=path_radar_data+'/'+folder_year+'/'+folder_day+'/'+indiv_file
plot_radar_slice(ax1,ax3,ax6,ax_nb,path_radar_slice,lines,folder_year,folder_day,indiv_file,technique,xls_icelenses,trafic_light,elevation_dictionnary)

#Plot date 3
folder_year='2003'
folder_day='may12'
indiv_file='may12_03_36_aggregated' #From up to down: need reversing! Already fone, OK!
ax_nb=4
path_radar_slice=path_radar_data+'/'+folder_year+'/'+folder_day+'/'+indiv_file
plot_radar_slice(ax1,ax4,ax6,ax_nb,path_radar_slice,lines,folder_year,folder_day,indiv_file,technique,xls_icelenses,trafic_light,elevation_dictionnary)

#pdb.set_trace()
#Plot date 4
folder_year='2003'
folder_day='may11'
indiv_file='may11_03_29_aggregated' #High elevation, no need: OK!
ax_nb=5
path_radar_slice=path_radar_data+'/'+folder_year+'/'+folder_day+'/'+indiv_file
plot_radar_slice(ax1,ax5,ax6,ax_nb,path_radar_slice,lines,folder_year,folder_day,indiv_file,technique,xls_icelenses,trafic_light,elevation_dictionnary)

#Plot all the 2002-2003 icelenses according to their confidence color 
#1. Red
ax1.scatter(df_2002_2003[df_2002_2003['colorcode_icelens']==-1]['lon_3413'],df_2002_2003[df_2002_2003['colorcode_icelens']==-1]['lat_3413'],s=1,facecolors='#c9662c', edgecolors='none')
#2. Orange
ax1.scatter(df_2002_2003[df_2002_2003['colorcode_icelens']==0]['lon_3413'],df_2002_2003[df_2002_2003['colorcode_icelens']==0]['lat_3413'],s=1,facecolors='#fed976', edgecolors='none')
#3. Green
ax1.scatter(df_2002_2003[df_2002_2003['colorcode_icelens']==1]['lon_3413'],df_2002_2003[df_2002_2003['colorcode_icelens']==1]['lat_3413'],s=1,facecolors='#238b45', edgecolors='none')
#Purple
ax1.scatter(df_2002_2003[df_2002_2003['colorcode_icelens']==2]['lon_3413'],df_2002_2003[df_2002_2003['colorcode_icelens']==2]['lat_3413'],s=1,facecolors='purple', edgecolors='none')

#Zoom on SW Greenland
ax1.set_xlim(-380100,106800)
ax1.set_ylim(-2810000,-2215200)

#Custom ylabel
ax6.set_ylim(950,2600)
start_ytick_elev, end_ytick_elev = ax6.get_ylim()
ax6.yaxis.set_ticks(np.arange(start_ytick_elev, end_ytick_elev, 250))
ax6.set_ylabel('Elevation [m]')

#Custom xlabel
ticks_xplot_elev=np.arange(0,1001,100)
ticks_xplot_elev[-1]=999
ax6.set_xticks(ticks_xplot_elev)
ax6.set_xticklabels([])
ax6.set_xlim(0,1000)
ax6.grid()

##Save the figure
#fig_name=[]
#fig_name='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/icelens_identification/indiv_traces_icelenses/2002_3_SWGr_icelenses.png'
#plt.savefig(fig_name,dpi=1000)
print('Done with SW Greenland plot')

#Plot the whole GrIS 2002-2003 radar tracks
#Prepare plot
fig, (ax1) = plt.subplots(1, 1)#, gridspec_kw={'width_ratios': [1, 3]})
fig.suptitle('2002-2003 ice lensing overview')

#Display elevation
cb=ax1.imshow(elevDem, extent=grid.extent,cmap=discrete_cmap(10,'cubehelix_r'),norm=divnorm,alpha=0.5)
cbar=fig.colorbar(cb, ax=[ax1], location='left')
cbar.set_label('Elevation [m]')

#Plot all the 2010-2014 icelenses
ax1.scatter(lon_3413_MacFerrin, lat_3413_MacFerrin,s=0.1,facecolors='cornflowerblue', edgecolors='none', label='2010-2014 ice slabs mapping by MacFerrin et al., 2019')

#Plot all the 2002-2003 flightlines
ax1.scatter(lon_all, lat_all,s=0.1,facecolors='lightgrey', edgecolors='none',alpha=0.1, label='2002-2003 flight lines')

#Plot all the 2002-2003 icelenses according to their condifence color
#1. Red
ax1.scatter(df_2002_2003[df_2002_2003['colorcode_icelens']==-1]['lon_3413'],df_2002_2003[df_2002_2003['colorcode_icelens']==-1]['lat_3413'],s=5,facecolors='#c9662c', edgecolors='none',label='High confidence')
#2. Orange
ax1.scatter(df_2002_2003[df_2002_2003['colorcode_icelens']==0]['lon_3413'],df_2002_2003[df_2002_2003['colorcode_icelens']==0]['lat_3413'],s=5,facecolors='#fed976', edgecolors='none',label='Medium confidence')
#3. Green
ax1.scatter(df_2002_2003[df_2002_2003['colorcode_icelens']==1]['lon_3413'],df_2002_2003[df_2002_2003['colorcode_icelens']==1]['lat_3413'],s=5,facecolors='#238b45', edgecolors='none',label='Low confidence')
#Purple
ax1.scatter(df_2002_2003[df_2002_2003['colorcode_icelens']==2]['lon_3413'],df_2002_2003[df_2002_2003['colorcode_icelens']==2]['lat_3413'],s=5,facecolors='purple', edgecolors='none',label='Bright layer - very low confidence')

#Correct zoom
ax1.set_xlim(-650000,900000)
ax1.set_ylim(-3360000,-650000)
plt.legend(loc='lower right')
plt.show()

##Save the figure
#fig_name=[]
#fig_name='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/icelens_identification/indiv_traces_icelenses/whole_GrIS_2002_3.png'
#plt.savefig(fig_name,dpi=1000)

print('Done with whole GrIS plot')
pdb.set_trace()
#######################################################################
###                 Identification of deepest ice lenses            ###
#######################################################################
         
#Time calculation variables
dt = 2.034489716724874e-09 #Timestep for 2002/2003 traces
t0 = 0; # Unknown so set to zero
#Compute the speed (Modified Robin speed):
# self.C / (1.0 + (coefficient*density_kg_m3/1000.0))
v= 299792458 / (1.0 + (0.734*0.873/1000.0))

#Path radar data:
path_radar_data= 'C:/Users/jullienn/switchdrive/Private/research/backup_Aglaja/working_environment/iceslabs_MacFerrin/data'

#Create the dictionary to save ice lens information
icelens_information={k: {} for k in list(xls_icelenses.keys())}

#Depth of ice lenses: use the variable 'xsl_icelenses'
for indiv_file in list(xls_icelenses.keys()):
    print(indiv_file)
    
    ####################################################################
    ###    Load the data of interest to retrieve the depth vector    ###
    
    #Define the folder_year
    if (indiv_file[6:8]=='02'):
        folder_year='2002'
    elif (indiv_file[6:8]=='03'):
        folder_year='2003'
    else:
        print('Problem in data!')
    
    path_radar_load=path_radar_data+'/'+folder_year+'/'+indiv_file[0:5]

    if (indiv_file[0:5]=='jun04'):
        fdata= scipy.io.loadmat(path_radar_load+'/'+indiv_file)
        #Select radar echogram
        radar_echo=fdata['data']
    else:
        #Open the file and read it
        f_agg = open(path_radar_load+'/'+indiv_file, "rb")
        data = pickle.load(f_agg)
        f_agg.close()
        
        #Select radar echogram
        radar_echo=data['radar_echogram']

    #Compute the vertical resolution
    Nt = radar_echo.shape[0]
    Time = t0 + dt*np.arange(1,Nt+1)
    #self.SAMPLE_DEPTHS = self.radar_speed_m_s * self.SAMPLE_TIMES / 2.0
    depths = v * Time / 2.0 
    
    ###    Load the data of interest to retrieve the depth vector    ###
    ####################################################################
    
    #Identify the index corresponding to the deepest depth (i.e. 20m)
    deepest_index=np.where(depths>20)[0][0]
    
    #Load ice lenses identifications
    df_temp=xls_icelenses[indiv_file]
    df_colnames = list(df_temp.keys())
    
    #pdb.set_trace()
    
    #Load trafic light information
    df_trafic_light=trafic_light[indiv_file]
    df_colnames_trafic_light = list(df_trafic_light.keys())
    
    #Get the storage dataframe ready
    x_empty=np.zeros(radar_echo.shape[1])
    x_empty[:]=np.nan
    
    df_icelenses_information=pd.DataFrame({'x':x_empty,
                                           'deepest_depth_index':x_empty,
                                           'deepest_depth':x_empty,
                                           'deepest_depth_color':x_empty})
    
    #Empty gathering vectors
    x_all=[]
    y_all=[]
    x_color=[]
    
    for i in range (0,int(len(df_colnames)),2):
        #Load x and y
        x_vect=df_temp[df_colnames[i]]
        y_vect=df_temp[df_colnames[i+1]]
        
        #Load trafic light color
        trafic_light_invid_color=df_colnames_trafic_light[i]
        
        #Floor the horizontal pixels
        x_floored=[]
        x_floored=np.floor(x_vect)
        
        #Floor the vertical pixels
        y_floored=[]
        y_floored=np.floor(y_vect)
        
        #Define the color code in which to display the ice lens
        if (trafic_light_invid_color[0:3]=='gre'):
            color_code=1
        elif (trafic_light_invid_color[0:3]=='ora'):
            color_code=0
        elif (trafic_light_invid_color[0:3]=='red'):
            color_code=-1
        elif (trafic_light_invid_color[0:3]=='pur'):
            color_code=2
        else:
            print('The color is not known!')       
        
        #Append x and y to gather vectors
        x_all=np.append(x_all,x_floored)
        y_all=np.append(y_all,y_floored)
        
        #Append the color codes
        x_color=np.append(x_color,np.ones(x_floored.shape[0])*color_code)
    
    #Remove the nans
    #pdb.set_trace()
    x_color=x_color[~np.isnan(x_all)]
    x_all=x_all[~np.isnan(x_all)]
    y_all=y_all[~np.isnan(y_all)]
    
    #Find the index of unique horizontal pixel
    x_all_unique, idx = np.unique(x_all, return_index=True)
    
    #Set the x_all_unique as integer
    x_all_unique=x_all_unique.astype(int)
    
    #We can have several minimums for one horizontal pixel
    for i in range (0,len(x_all_unique),1):
        #pdb.set_trace()
        
        #Find all the index having the same horizontal pixel value
        index_element_search=np.where(x_all == x_all_unique[i])[0]
        
        #Select all the correponding vertical pixel values
        y_index=y_all[index_element_search]
        #Select all the corresponding xcolor
        x_color_index=x_color[index_element_search]
        
        #For y_index > deepest index, store the deepest index (correpond to 20m depth)
        y_index[y_index>deepest_index]=deepest_index
        
        #Keep the deepest one
        deepest_pixel_index=np.nanmax(y_index)
        #The corresponding location is given by y_index.argmax()
        #Associate the deepest pixel with its color code
        x_color_value=x_color_index[y_index.argmax()]
        
        #If nan, continue and do not store anything
        if np.isnan(deepest_pixel_index):
            print('Identified ice lens all below 20m deep, continue')
            continue
        
        #Retrieve the corresponding deepest depth
        deepest_depth=depths[deepest_pixel_index.astype(int)]
        
        ########################## old method ###########################
        ##Retreive the depth. What is stored in y_all are the y coordinates, which
        ##correspond to the index!
        #depths_retrieval=depths[y_all[index_element_search].astype(int)]
        #
        ##Remove the depths deeper than 20m
        #depths_retrieval[depths_retrieval>20]=np.nan
        #
        ##Retrieve the deepest depth between 0 and 20m deep:
        #deepest_depth=np.nanmax(depths_retrieval)
        ########################## old method ###########################

        #store the information in the dataframe
        df_icelenses_information['x'][x_all_unique[i]]=x_all_unique[i]
        df_icelenses_information['deepest_depth_index'][x_all_unique[i]]=deepest_pixel_index
        df_icelenses_information['deepest_depth'][x_all_unique[i]]=deepest_depth
        df_icelenses_information['deepest_depth_color'][x_all_unique[i]]=x_color_value
    
    #Moving window averaging
    for j in range(4,len(df_icelenses_information['deepest_depth_index'])-5,1):
        if (np.isnan(df_icelenses_information['deepest_depth_index'][j])):
            continue
        
        #Define the moving window as considering -4 and +4 around it, without considering j
        moving_window_temp=df_icelenses_information['deepest_depth_index'][j-4:j+5]
        
        #if (indiv_file=='jun04_02proc_4.mat'):
        #    pdb.set_trace()
        
        if (np.sum(np.isnan(np.asarray(list(moving_window_temp))).astype(int))==8):
            #Pixel of interest is surrounded by NaNs.
            #Let's first try to not bother about it
            continue
        
        #pdb.set_trace()
        #1. Depth index
        #Get rid of the dependance with df_icelenses_information
        moving_window=np.asarray(list(moving_window_temp))
        #Removing the jth element of interest
        moving_window[4]=np.nan
        moving_average=np.nanmean(moving_window)        
        moving_std=np.nanstd(moving_window)
        
        #2. Depth
        #Get rid of the dependance with df_icelenses_information
        window_depth=np.asarray(list(df_icelenses_information['deepest_depth'][j-4:j+5]))
        #Removing the jth element of interest
        window_depth[4]=np.nan
        moving_average_depth=np.nanmean(window_depth)
        
        #3. Color code
        #Get rid of the dependance with df_icelenses_information
        window_color=np.asarray(list(df_icelenses_information['deepest_depth_color'][j-4:j+5]))
        #Removing the jth element of interest
        window_color[4]=np.nan
        moving_average_color=np.round(np.nanmean(window_color))
        
        pixel_studied=df_icelenses_information['deepest_depth_index'][j]

        if ((pixel_studied>(moving_average+2*moving_std)) or (pixel_studied<(moving_average-2*moving_std))):
            #The jumped index is recalculated as a function of is neighboors
            #pdb.set_trace()
            df_icelenses_information['deepest_depth_index'][j]=np.round(moving_average).astype(int)
            df_icelenses_information['deepest_depth'][j]=moving_average_depth
            df_icelenses_information['deepest_depth_color'][j]=moving_average_color
    #Save the dataframe into a dictionnary
    icelens_information[indiv_file]=df_icelenses_information

#Display all the traces with the corresponding depth on the map
################################### Plot ##################################
print('Save indiv files with deepest ice lenses identification')
#Plot all the dates:
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
                print('Plot the deppest ice lenses')
                #Define the path of radar data
                path_radar_slice=path_radar_data+'/'+year+'/'+days+'/'+indiv_file
                
                #Prepare plot
                fig = plt.figure(figsize=(19,10))
                fig.suptitle(indiv_file)
                gs = gridspec.GridSpec(10, 20)
                gs.update(wspace=0.1)
                gs.update(wspace=0.001)
                ax1 = plt.subplot(gs[0:6, 0:10])
                ax2 = plt.subplot(gs[0:6, 12:20])
                ax3 = plt.subplot(gs[6:10, 0:20])
                
                #Display elevation
                cb1=ax1.imshow(elevDem, extent=grid.extent,cmap=discrete_cmap(10,'cubehelix_r'),alpha=0.5,norm=divnorm)
                ax1.set_title('Ice lenses and slabs location')
                cbar1=fig.colorbar(cb1, ax=[ax1], location='left')
                cbar1.set_label('Elevation [m]')
                
                #Plot all the 2010-2014 icelenses
                ax1.scatter(lon_3413_MacFerrin, lat_3413_MacFerrin,s=1,facecolors='cornflowerblue', edgecolors='none')
                #ax1.scatter(lon_3413_MacFerrin, lat_3413_MacFerrin,color='red',marker='o',alpha=0.2)
                
                #Plot all the 2002-2003 flightlines
                ax1.scatter(lon_all, lat_all,s=1,facecolors='lightgrey', edgecolors='none',alpha=0.1)
                
    
                #Plot the individual results
                plot_radar_slice_with_thickness(ax1,ax2,ax3,path_radar_slice,lines,year,days,indiv_file,technique,xls_icelenses,trafic_light,elevation_dictionnary,icelens_information)
                #pdb.set_trace()
                
                ##Save the figure
                #fig_name=[]
                #fig_name='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/icelens_identification/indiv_traces_icelenses/deepest_lenses'+indiv_file+'.png'
                #plt.savefig(fig_name,dpi=1000)
                plt.close()
                print('Done with deepest',indiv_file)
                
#######################################################################
###                 Identification of deepest ice lenses            ###
#######################################################################
