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


def load_2002_2003_radargram(path_radar_slice,lines,folder_year,folder_day,indiv_file):
    #this function loads the individual 2002-2003 radargram
    
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
    #Transpose if june 04
    if (folder_day=='jun04'):
        lat_3413=np.transpose(lat_3413)
        lon_3413=np.transpose(lon_3413)
    #Calculate the distances (in m)
    distances=compute_distances(lon_3413,lat_3413)
    #Convert distances from m to km
    distances=distances/1000
    
    #Save lat, lon and radar slice from 0 to 30m deep into a dictionnary
    dict_returned={'lat_3413':lat_3413,'lon_3413':lon_3413,'radar_slice_0_30m':radar_slice,'distances':distances,'depths':depths}
    
    return dict_returned


def plot_radar_slice(ax_map,ax_plot,ax_nb,path_radar_slice,lines,folder_year,folder_day,indiv_file,technique,xls_icelenses,trafic_light,elevation_dictionnary):
        
    #load radargram
    radargram_data = load_2002_2003_radargram(path_radar_slice,lines,folder_year,folder_day,indiv_file)
        
    #Display on the map where is this track
    ax_map.scatter(radargram_data['lon_3413'], radargram_data['lat_3413'],s=1,facecolors='black', edgecolors='black')
    
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
    ticks_yplot=np.linspace(0,radargram_data['radar_slice_0_30m'].shape[0],4).astype(int)
    
    #Plot the radar slice
    cb2=ax_plot.pcolor(radargram_data['radar_slice_0_30m'],cmap=plt.get_cmap('gray'))
    ax_plot.set_ylim(0,radargram_data['radar_slice_0_30m'].shape[0])
    ax_plot.invert_yaxis() #Invert the y axis = avoid using flipud.
    #Colorbar custom
    cb2.set_clim(perc_lower_end,perc_upper_end)
    #Set the y ticks
    ax_plot.set_yticks(ticks_yplot) 
    ax_plot.set_yticklabels(np.round(radargram_data['depths'][ticks_yplot]).astype(int))
    #Set ylabel
    ax_plot.set_ylabel('Depth [m]')
        
    #Distance from start of the trace
    if (ax_nb==2):
        letter_elev='b'
        #Display on the map the letter corresponding to the radar slice
        ax_map.text(np.nanmedian(radargram_data['lon_3413'])-15000,np.nanmedian(radargram_data['lat_3413']),letter_elev,color='black',fontsize=20)
    
    elif (ax_nb==3):
        #ax_plot.set_title('Percolation zone - ice lenses',fontsize=10)
        letter_elev='c'
        #Display on the map the letter corresponding to the radar slice
        ax_map.text(np.nanmedian(radargram_data['lon_3413'])-30000,np.nanmedian(radargram_data['lat_3413'])-10000,letter_elev,color='black',fontsize=20)
    
    elif (ax_nb==4):
        #ax_plot.set_title('Percolation zone - ice slabs',fontsize=10)
        letter_elev='d'
        #Display on the map the letter corresponding to the radar slice
        ax_map.text(np.nanmedian(radargram_data['lon_3413'])-15000,np.nanmedian(radargram_data['lat_3413'])+5000,letter_elev,color='black',fontsize=20)
    
    elif (ax_nb==5):
        #ax_plot.set_title('Dry snow zone',fontsize=10)
        ax_plot.set_xlabel('Distance [km]')
        letter_elev='e'
        #Display on the map the letter corresponding to the radar slice
        ax_map.text(np.nanmedian(radargram_data['lon_3413']),np.nanmedian(radargram_data['lat_3413'])-37000,letter_elev,color='black',fontsize=20)
        
        #Display colorbar. This is from FigS1.py
        cbar_depth=fig.colorbar(cb2, cax=axc)#aspect is from https://stackoverflow.com/questions/33443334/how-to-decrease-colorbar-width-in-matplotlib
        cbar_depth.set_label('Radar signal strength [dB]')
        axc.yaxis.labelpad = -5#this is from https://stackoverflow.com/questions/6406368/matplotlib-move-x-axis-label-downwards-but-not-x-axis-ticks
    
    #Display the correspondance between radar slice and radar location on the map
    ax_plot.text(0.0172, 0.8975,'   ',backgroundcolor='white',ha='center', va='center', transform=ax_plot.transAxes,fontsize=17)#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of
    ax_plot.text(0.0172, 0.8975,letter_elev, ha='center', va='center', transform=ax_plot.transAxes,fontsize=25)#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot
    
    #Display the ice lenses identification:
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
                color_to_display='#ffb300'#'#00441b'
            elif (trafic_light_indiv_color[0:3]=='ora'):
                color_to_display='#fed976'
                continue
            elif (trafic_light_indiv_color[0:3]=='red'):
                color_to_display='#c9662c'
                continue
            elif (trafic_light_indiv_color[0:3]=='pur'):
                color_to_display='purple'
                continue
            else:
                print('The color is not known!')
            #Display ice lens
            ax_plot.plot(x_vect,y_vect,color=color_to_display,linestyle='dashed',linewidth=2.5)
    
    #Load the elevation profile
    elevation_vector=elevation_dictionnary[folder_year][folder_day][indiv_file]
  
    #Define the dates that need reversed display
    list_reverse_agg_mat=['may12_03_36_aggregated','may14_03_51_aggregated',
                       'may13_03_29_aggregated','may30_02_51_aggregated',
                       'may24_02_25_aggregated','may15_03_37_aggregated',
                       'may11_03_29_aggregated','jun04_02proc_52.mat',
                       'jun04_02proc_53.mat']
        
    #Order the radar track from down to up if needed      
    if (indiv_file in list(list_reverse_agg_mat)):
        ax_plot.set_xlim(radargram_data['radar_slice_0_30m'].shape[1],0)
        #Reverse the distances vector:
        distances=np.flipud(radargram_data['distances'])
    else:
        distances=radargram_data['distances']
    
    #Generate the pick for horizontal distance display
    ticks_xplot=np.arange(0,distances.shape[0]+1,100)
    #Plot also the last index
    ticks_xplot[-1]=distances.shape[0]-1
    #Set x ticks
    ax_plot.set_xticks(ticks_xplot) 
    #Display the distances from the origin as being the x label
    ax_plot.set_xticklabels(np.round(distances[ticks_xplot]).astype(int))
    
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
    
    return



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
from matplotlib.patches import Patch
import matplotlib.patches as patches
from matplotlib.lines import Line2D
from scalebar import scale_bar

save_2002_2003_data='FALSE'
technique='perc_2p5_97p5'
making_down_to_up='FALSE'

################# Load 2002-2003 flightlines coordinates ################
path_data='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/icelens_identification'

#Open the file and read it
f_flightlines = open(path_data+'/metadata_coord_2002_2003', "rb")
all_2002_3_flightlines = pickle.load(f_flightlines)
f_flightlines.close()
################# Load 2002-2003 flightlines coordinates ################

############################ Load GrIS elevation ##############################
#This is from fig3_paper_iceslabs.py
#Open and display satelite image behind map
from pyproj import CRS
import rioxarray as rxr
import cartopy.crs as ccrs
import geopandas as gpd
#This section of displaying sat data was coding using tips from
#https://www.earthdatascience.org/courses/use-data-open-source-python/intro-raster-data-python/raster-data-processing/reproject-raster/
#https://towardsdatascience.com/visualizing-satellite-data-using-matplotlib-and-cartopy-8274acb07b84

path_DEM='C:/Users/jullienn/switchdrive/Private/research/backup_Aglaja/working_environment/greenland_topo_data/elevations/'
#Load DEM for display
GrIS_DEM_display = rxr.open_rasterio(path_DEM+'greenland_dem_mosaic_100m_v3.0.tif',
                              masked=True).squeeze()
#Load DEM for elevation pick up
GrIS_DEM = rasterio.open(path_DEM+'greenland_dem_mosaic_100m_v3.0.tif')

###################### From Tedstone et al., 2022 #####################
#from plot_map_decadal_change.py
# Define the CartoPy CRS object.
crs = ccrs.NorthPolarStereo(central_longitude=-45., true_scale_latitude=70.)
# This can be converted into a `proj4` string/dict compatible with GeoPandas
crs_proj4 = crs.proj4_init
###################### From Tedstone et al., 2022 #####################

#Define extents
#ease_extent = [west limit, east limit., south limit, north limit]
extent_DEM = [np.asarray(GrIS_DEM_display.x[0]), np.asarray(GrIS_DEM_display.x[-1]), np.asarray(GrIS_DEM_display.y[-1]), np.asarray(GrIS_DEM_display.y[0])]

'''
plt.figure(figsize=(14,10))
ax = plt.axes(projection=crs)
#ax.set_extent(ease_extent, crs=crs)
ax.imshow(GrIS_DEM_display[:,:], extent=extent_DEM, transform=crs,cmap='gray', origin='upper') #NIR
ax.gridlines(color='gray', linestyle='--')
ax.coastlines()
#ax.set_xlim(extent_image[0],extent_image[1])
#ax.set_ylim(extent_image[3],extent_image[2])
plt.tight_layout()
'''
############################ Load GrIS elevation ##############################

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
                lat_elev=[]
                lon_elev=[]
                if (days=='jun04'):
                    #Append all the flightlines
                    lat_all=np.append(lat_all,all_2002_3_flightlines[year][days][indiv_file][0][0])
                    lon_all=np.append(lon_all,all_2002_3_flightlines[year][days][indiv_file][1][0])
                    
                    #For elevation extraction:
                    lat_elev=np.transpose(all_2002_3_flightlines[year][days][indiv_file][0][0])
                    lon_elev=np.transpose(all_2002_3_flightlines[year][days][indiv_file][1][0])
                else:
                    #Append all the flightlines
                    lat_all=np.append(lat_all,all_2002_3_flightlines[year][days][indiv_file][0])
                    lon_all=np.append(lon_all,all_2002_3_flightlines[year][days][indiv_file][1])
                    
                    #For elevation extraction:
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
                        #This is adpated from fi3_paper_iceslabs.py, originally from extract_elevation.py
                        #This is from https://gis.stackexchange.com/questions/190423/getting-pixel-values-at-single-point-using-rasterio
                        for val in GrIS_DEM.sample([(indiv_coord[0], indiv_coord[1])]):
                            #Calculate the corresponding elevation
                            elev_indiv_file=np.append(elev_indiv_file,val)
                
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

######## Load 2010-2018 ice slabs location from Jullien et al., 2022 ##########
#Load the data
filename_Jullien= 'C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2002_2018/final_excel/high_estimate/clipped/Ice_Layer_Output_Thicknesses_2010_2018_jullienetal2023_high_estimate_cleaned.csv'
#Read Jullien data thanks to https://stackoverflow.com/questions/65254535/xlrd-biffh-xlrderror-excel-xlsx-file-not-supported
df_20102018 = pd.read_csv(filename_Jullien,delimiter=',',decimal='.')

#Transform the coordinated from WGS84 to EPSG:3413
#Example from: https://pyproj4.github.io/pyproj/stable/examples.html
transformer = Transformer.from_crs("EPSG:4326", "EPSG:3413", always_xy=True)
points=transformer.transform(np.array(df_20102018.lon),np.array(df_20102018.lat))

lon_3413_20102018=points[0]
lat_3413_20102018=points[1]
######## Load 2010-2018 ice slabs location from Jullien et al., 2022 ##########

### -------------------------- Load shapefiles --------------------------- ###
#From fig3_paper_iceslabs.py
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

#Load high estimates ice slabs extent 2010-2018, manually drawn on QGIS
path_iceslabs_shape='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2002_2018/shapefiles/'
iceslabs_jullien_highend_20102018=gpd.read_file(path_iceslabs_shape+'iceslabs_jullien_highend_20102018.shp')

### -------------------------- Load shapefiles --------------------------- ###
#Set fontsize
plt.rcParams.update({'font.size': 15}) #from https://stackoverflow.com/questions/3899980/how-to-change-the-font-size-on-a-matplotlib-plot

################################### Plot ##################################
#Prepare plot
plt.rcParams["figure.figsize"] = (22,11.3)#from https://pythonguides.com/matplotlib-increase-plot-size/
fig = plt.figure()
#fig.suptitle('2002-2003 ice lenses and ice slabs mapping SW Greenland')
gs = gridspec.GridSpec(8, 200)
gs.update(hspace=0.5)
gs.update(wspace=0.5)

'''
ax1 = plt.subplot(gs[0:5, 5:55],projection=crs)
ax_InsetMap = plt.subplot(gs[5:8, 0:45],projection=crs)
'''
ax1 = plt.subplot(gs[0:8, 5:55],projection=crs)
ax2 = plt.subplot(gs[0:2, 65:195])
ax3 = plt.subplot(gs[2:4, 65:195])
ax4 = plt.subplot(gs[4:6, 65:195])
ax5 = plt.subplot(gs[6:8, 65:195])
axc = plt.subplot(gs[0:8, 197:200])

#Load DEM clipped over the SW
GrIS_DEM_display_SW = rxr.open_rasterio(path_DEM+'SW_zoom/greenland_dem_mosaic_100m_v3.0_SW.tif',
                              masked=True).squeeze()
extent_DEM_SW = [np.asarray(GrIS_DEM_display_SW.x[0]), np.asarray(GrIS_DEM_display_SW.x[-1]), np.asarray(GrIS_DEM_display_SW.y[-1]), np.asarray(GrIS_DEM_display_SW.y[0])]

#Display elevation
#cb1=ax1.imshow(GrIS_DEM_display[:,:], extent=extent_DEM, transform=crs, origin='upper', cmap='Blues_r',zorder=1,alpha=0.5)
#cb1=ax1.imshow(GrIS_DEM_display_SW[:,:], extent=extent_DEM_SW, transform=crs, origin='upper', cmap='Blues_r',zorder=1,alpha=0.5)
#cbar1=fig.colorbar(cb1, ax=[ax1], location='right')

#Display SW and CW regions
SW_rignotetal.plot(ax=ax1,color='white', edgecolor='black',linewidth=0.5) 
CW_rignotetal.plot(ax=ax1,color='white', edgecolor='black',linewidth=0.5) 

#Add region names
ax1.text(30000,-2369000,'CW',color='black')
ax1.text(-12830,-2600000,'SW',color='black')
ax1.text(82500,-2609000,'SE',color='black')

#Display contours
cont=ax1.contour(GrIS_DEM_display_SW[:,:], levels=np.arange(1000,3750,250), extent=extent_DEM_SW, transform=crs, origin='upper', colors=['#8c510a'],linewidth=0.1)

#Add elevation contour values
ax1.text(-0.05, 0.42,'1000', ha='center', va='center', transform=ax1.transAxes,fontsize=10,color='#8c510a')#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot
ax1.text(-0.05, 0.39,'1250', ha='center', va='center', transform=ax1.transAxes,fontsize=10,color='#8c510a')#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot
ax1.text(-0.05, 0.235,'1500', ha='center', va='center', transform=ax1.transAxes,fontsize=10,color='#8c510a')#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot
ax1.text(-0.05, 0.18,'1500', ha='center', va='center', transform=ax1.transAxes,fontsize=10,color='#8c510a')#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot
ax1.text(0.055, -0.02,'1750', ha='center', va='center', rotation=90,transform=ax1.transAxes,fontsize=10,color='#8c510a')#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot
ax1.text(0.12, -0.02,'2000', ha='center', va='center', rotation=90,transform=ax1.transAxes,fontsize=10,color='#8c510a')#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot
ax1.text(0.25, -0.02,'2250', ha='center', va='center', rotation=90,transform=ax1.transAxes,fontsize=10,color='#8c510a')#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot
ax1.text(0.39, -0.02,'2500', ha='center', va='center', rotation=90,transform=ax1.transAxes,fontsize=10,color='#8c510a')#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot
ax1.text(1.06, 0.455,'2500', ha='center', va='center',transform=ax1.transAxes,fontsize=10,color='#8c510a')#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot

#Zoom over SW Greenland
ease_extent = [-201529, 106696, -2889749, -2160162]
ax1.set_extent(ease_extent, crs=crs) 

#Shapefiles
# --- 2010-2018
iceslabs_jullien_highend_20102018.plot(ax=ax1,color='#d73027', edgecolor='none',linewidth=0.5)

'''
#Plot all the 2010-2018 ice slabs
ax1.scatter(lon_3413_20102018, lat_3413_20102018,s=1,facecolors='#0570b0', edgecolors='none')
'''
#Plot all the 2002-2003 flightlines
ax1.scatter(lon_all, lat_all,s=1,facecolors='#969696', edgecolors='none',alpha=0.1)
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
path_radar_data='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data'

#Plot date 1
folder_year='2003'
folder_day='may11'
indiv_file='may11_03_1_aggregated' #From down to up: OK!
ax_nb=2
path_radar_slice=path_radar_data+'/'+folder_year+'/'+folder_day+'/'+indiv_file
plot_radar_slice(ax1,ax2,ax_nb,path_radar_slice,lines,folder_year,folder_day,indiv_file,technique,xls_icelenses,trafic_light,elevation_dictionnary)

#Plot date 2
folder_year='2002'
folder_day='jun04'
indiv_file='jun04_02proc_53.mat' #From up to down: need reversing! Already done, OK!
ax_nb=3
path_radar_slice=path_radar_data+'/'+folder_year+'/'+folder_day+'/'+indiv_file
plot_radar_slice(ax1,ax3,ax_nb,path_radar_slice,lines,folder_year,folder_day,indiv_file,technique,xls_icelenses,trafic_light,elevation_dictionnary)

#Plot date 3
folder_year='2003'
folder_day='may12'
indiv_file='may12_03_36_aggregated' #From up to down: need reversing! Already fone, OK!
ax_nb=4
path_radar_slice=path_radar_data+'/'+folder_year+'/'+folder_day+'/'+indiv_file
plot_radar_slice(ax1,ax4,ax_nb,path_radar_slice,lines,folder_year,folder_day,indiv_file,technique,xls_icelenses,trafic_light,elevation_dictionnary)

#Plot date 4
folder_year='2003'
folder_day='may11'
indiv_file='may11_03_29_aggregated' #High elevation, no need: OK!
ax_nb=5
path_radar_slice=path_radar_data+'/'+folder_year+'/'+folder_day+'/'+indiv_file
plot_radar_slice(ax1,ax5,ax_nb,path_radar_slice,lines,folder_year,folder_day,indiv_file,technique,xls_icelenses,trafic_light,elevation_dictionnary)

#Plot all the 2002-2003 icelenses according to their confidence color 
'''
#1. Red
ax1.scatter(df_2002_2003[df_2002_2003['colorcode_icelens']==-1]['lon_3413'],df_2002_2003[df_2002_2003['colorcode_icelens']==-1]['lat_3413'],s=2,facecolors='#c9662c', edgecolors='#c9662c')
#2. Orange
ax1.scatter(df_2002_2003[df_2002_2003['colorcode_icelens']==0]['lon_3413'],df_2002_2003[df_2002_2003['colorcode_icelens']==0]['lat_3413'],s=2,facecolors='#fed976', edgecolors='#fed976')
'''
#3. Green
ax1.scatter(df_2002_2003[df_2002_2003['colorcode_icelens']==1]['lon_3413'],df_2002_2003[df_2002_2003['colorcode_icelens']==1]['lat_3413'],s=2,facecolors='#ffb300', edgecolors='#ffb300')
'''
#Purple
ax1.scatter(df_2002_2003[df_2002_2003['colorcode_icelens']==2]['lon_3413'],df_2002_2003[df_2002_2003['colorcode_icelens']==2]['lat_3413'],s=2,facecolors='purple', edgecolors='purple')
'''

#Display lat/lon lines in map
gl=ax1.gridlines(draw_labels=True, xlocs=[-42, -44, -46, -48, -50], ylocs=[62, 63, 64, 65, 66, 67, 68, 69, 70], x_inline=False, y_inline=False,linewidth=0.5,linestyle='dashed')
#Customize lat and lon labels
gl.ylabels_right = False
gl.xlabels_bottom = False

#Custom legend myself,  line2D from https://stackoverflow.com/questions/39500265/how-to-manually-create-a-legend, marker from https://stackoverflow.com/questions/47391702/how-to-make-a-colored-markers-legend-from-scratch
legend_elements = [Line2D([0], [0], label='Ice sheet regional divide', color='black', linewidth=0.5),
                   Line2D([0], [0], label='Elevation contours', color='#8c510a'),
                   Line2D([0], [0], label='2002-2003 flightlines', color='#969696'),
                   Line2D([0], [0], label='Transects of interest', color='k', linewidth=3),
                   Patch(facecolor='#d73027',label='2010-2018 ice slabs'),
                   Line2D([0], [0], label='2002-2003 ice layers and slabs', color='#ffb300', linewidth=3)]
#Add legend
ax1.legend(handles=legend_elements,loc='upper left',fontsize=12,framealpha=1,bbox_to_anchor=(-0.019, 1.008))
plt.legend()
axc.legend_.remove()

#Display the map panel label
ax1.text(-0.06, 1.025,'a', ha='center', va='center', transform=ax1.transAxes,fontsize=25)#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot
'''
#Display GrIS inset map
ax_InsetMap.coastlines(edgecolor='black',linewidth=0.075)
#Display GrIS drainage bassins limits
GrIS_drainage_bassins.plot(ax=ax_InsetMap,color='none', edgecolor='black',linewidth=0.075)
#Display region name
ax_InsetMap.text(NO_rignotetal.centroid.x-200000,NO_rignotetal.centroid.y-80000,np.asarray(NO_rignotetal.SUBREGION1)[0])
ax_InsetMap.text(NE_rignotetal.centroid.x-200000,NE_rignotetal.centroid.y+20000,np.asarray(NE_rignotetal.SUBREGION1)[0])
ax_InsetMap.text(SE_rignotetal.centroid.x-100000,SE_rignotetal.centroid.y,np.asarray(SE_rignotetal.SUBREGION1)[0])
ax_InsetMap.text(SW_rignotetal.centroid.x-185000,SW_rignotetal.centroid.y-120000,np.asarray(SW_rignotetal.SUBREGION1)[0])
ax_InsetMap.text(CW_rignotetal.centroid.x-200000,CW_rignotetal.centroid.y-100000,np.asarray(CW_rignotetal.SUBREGION1)[0])
ax_InsetMap.text(NW_rignotetal.centroid.x-150000,NW_rignotetal.centroid.y-150000,np.asarray(NW_rignotetal.SUBREGION1)[0])

#Display rectangle around datalocation - this is from Fig. 3.py   
#Extract corner coordinates
coord_origin=[ax1.get_xlim()[0],ax1.get_ylim()[0]]
coord_topright=[ax1.get_xlim()[1],ax1.get_ylim()[1]]
#This is from https://stackoverflow.com/questions/37435369/matplotlib-how-to-draw-a-rectangle-on-image
# Create a Rectangle patch
rect = patches.Rectangle((coord_origin[0],coord_origin[1]),
                         np.abs(coord_origin[0]-coord_topright[0]),
                         np.abs(coord_origin[1]-coord_topright[1]),
                         angle=0, linewidth=1, edgecolor='black', facecolor='none')
# Add the patch to the Axes
ax_InsetMap.add_patch(rect)
ax_InsetMap.axis('on')

#Display scalebar
scale_bar(ax_InsetMap, (0.745, 0.125), 200, 3,5)# axis, location (x,y), length, linewidth, rotation of text
'''
scale_bar(ax1, (0.32, 0.0425), 50, 3,0)# axis, location (x,y), length, linewidth, rotation of text
#by measuring on the screen, the difference in precision between scalebar and length of transects is about ~200m

plt.show()
pdb.set_trace()

#Save the figure
fig_name=[]
#fig_name='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/icelens_identification/indiv_traces_icelenses/2002_3_SWGr_icelenses.png'
fig_name='C:/Users/jullienn/switchdrive/Private/research/RT1/figures/S1/v6/Fig_S1.png'
plt.savefig(fig_name,dpi=300,bbox_inches='tight') #bbox_inches is from https://stackoverflow.com/questions/32428193/saving-matplotlib-graphs-to-image-as-full-screen
print('Done with SW Greenland plot')

######################## Save 2002-2003 radargram data ########################

if (save_2002_2003_data=='TRUE'):
    #1 Loop over all 2002-2003 data
    for year in list(icelens_2002_3_flightlines.keys()):
        for days in list(icelens_2002_3_flightlines[year].keys()):
            for indiv_file in list(icelens_2002_3_flightlines[year][days].keys()):
                print(indiv_file)
                if (indiv_file[0:7]=='quality'):
                    print('Quality file, continue')
                    continue
                else:
                    #2. Load radargrams
                    #Define path
                    path_radar_slice=path_radar_data+'/'+year+'/'+days+'/'+indiv_file
                    #load radargram
                    radargram_data = load_2002_2003_radargram(path_radar_slice,lines,year,days,indiv_file)
                    #reset depths so that it matches with radagrams 0-30m dimensions
                    radargram_data['depths']=radargram_data['depths'][0:radargram_data['radar_slice_0_30m'].shape[0]]              
                    
                    #briefly check the data being generated
                    fig = plt.figure(figsize=(15,7))
                    ax1 = plt.subplot()
                    #Plot the radar slice
                    cb=ax1.pcolor(radargram_data['radar_slice_0_30m'],cmap=plt.get_cmap('gray'))
                    ax1.set_ylim(0,radargram_data['radar_slice_0_30m'].shape[0])
                    ax1.invert_yaxis() #Invert the y axis = avoid using flipud.
                    #Colorbar custom
                    if (year=='2002'):
                        perc_lower_end=-0.5709792307554173
                        perc_upper_end=0.7082634842114803
                    elif (year=='2003'):
                        perc_lower_end=-0.6061610403154447
                        perc_upper_end=0.7572821079440079     
                    cb.set_clim(perc_lower_end,perc_upper_end)
                    #Set the y ticks
                    ticks_yplot=np.linspace(0,radargram_data['radar_slice_0_30m'].shape[0]-1,4).astype(int)
    
                    ax1.set_yticks(ticks_yplot) 
                    ax1.set_yticklabels(np.round(radargram_data['depths'][ticks_yplot]).astype(int))
                    #Set ylabel
                    ax1.set_ylabel('Depth [m]')
                    #Generate the pick for horizontal distance display
                    ticks_xplot=np.linspace(0,radargram_data['distances'].shape[0],10).astype(int)
                    #Plot also the last index
                    ticks_xplot[-1]=radargram_data['distances'].shape[0]-1
                    #Set x ticks
                    ax1.set_xticks(ticks_xplot) 
                    #Display the distances from the origin as being the x label
                    ax1.set_xticklabels(np.round(radargram_data['distances'][ticks_xplot]).astype(int))
                    ax1.set_xlabel('Distance [km]')
                    ax1.set_title(indiv_file)
                    plt.show()
                                        
                    #3. Save as pickle
                    path_save_radargrams='C:/Users/jullienn/switchdrive/Private/research/RT1/final_dataset_2002_2018/2002_2003/radargram_data/'
                    filename_tosave=path_save_radargrams+year+'/L1_'+indiv_file.split('.')[0]+'.pickle'
                    
                    outfile= open(filename_tosave, "wb" )
                    pickle.dump(radargram_data,outfile)
                    outfile.close()
                    
                    #Save figure
                    plt.savefig(path_save_radargrams+'fig_check/'+indiv_file.split('.')[0]+'.png',dpi=300,bbox_inches='tight')
                    plt.close()
                    
                    #bbox_inches is from https://stackoverflow.com/questions/32428193/saving-matplotlib-graphs-to-image-as-full-screen)      
######################## Save 2002-2003 radargram data ########################

pdb.set_trace()

#Plot the whole GrIS 2002-2003 radar tracks
#Prepare plot
fig = plt.figure(figsize=(19,10))
gs = gridspec.GridSpec(10, 20)
ax1 = plt.subplot(gs[0:10, 0:20],projection=crs)

#Display coastlines
ax1.coastlines(edgecolor='black',linewidth=0.75)

#Display regions
NO_rignotetal.plot(ax=ax1,color='white', edgecolor='#081d58',linewidth=0.5) 
NE_rignotetal.plot(ax=ax1,color='white', edgecolor='#081d58',linewidth=0.5) 
SE_rignotetal.plot(ax=ax1,color='white', edgecolor='#081d58',linewidth=0.5) 
CW_rignotetal.plot(ax=ax1,color='white', edgecolor='#081d58',linewidth=0.5) 
SW_rignotetal.plot(ax=ax1,color='white', edgecolor='#081d58',linewidth=0.5) 
NW_rignotetal.plot(ax=ax1,color='white', edgecolor='#081d58',linewidth=0.5) 

#Display region name
ax1.text(NO_rignotetal.centroid.x,NO_rignotetal.centroid.y+20000,np.asarray(NO_rignotetal.SUBREGION1)[0])
ax1.text(NE_rignotetal.centroid.x,NE_rignotetal.centroid.y+20000,np.asarray(NE_rignotetal.SUBREGION1)[0])
ax1.text(SE_rignotetal.centroid.x,SE_rignotetal.centroid.y+20000,np.asarray(SE_rignotetal.SUBREGION1)[0])
ax1.text(SW_rignotetal.centroid.x,SW_rignotetal.centroid.y+20000,np.asarray(SW_rignotetal.SUBREGION1)[0])
ax1.text(CW_rignotetal.centroid.x,CW_rignotetal.centroid.y+20000,np.asarray(CW_rignotetal.SUBREGION1)[0])
ax1.text(NW_rignotetal.centroid.x,NW_rignotetal.centroid.y+20000,np.asarray(NW_rignotetal.SUBREGION1)[0])

#Plot all the 2002-2003 flightlines
ax1.scatter(lon_all, lat_all,s=0.1,facecolors='lightgrey', edgecolors='none',alpha=0.1, label='2002-2003 flight lines')

#Plot all the 2010-2018 ice slabs
ax1.scatter(lon_3413_20102018, lat_3413_20102018,s=1,facecolors='#0570b0', label='2010-18 ice slabs')

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
plt.legend(loc='lower right', fontsize=9)

###################### From Tedstone et al., 2022 #####################
#from plot_map_decadal_change.py
# Define the CartoPy CRS object.
crs = ccrs.NorthPolarStereo(central_longitude=-45., true_scale_latitude=70.)
# This can be converted into a `proj4` string/dict compatible with GeoPandas
crs_proj4 = crs.proj4_init
# x0, x1, y0, y1
ax1.set_extent([-692338, 916954, -3392187, -627732], crs=crs)
gl=ax1.gridlines(draw_labels=True, xlocs=[-50, -35], ylocs=[65, 75], x_inline=False, y_inline=False, color='#969696',linewidth=0.5)

#Customize lat labels
ax1.axis('off')
#scalebar.scale_bar(ax, (0, 0), 300, zorder=200)
###################### From Tedstone et al., 2022 #####################

plt.show()
pdb.set_trace()

#Save the figure
fig_name=[]
fig_name='C:/Users/jullienn/switchdrive/Private/research/RT1/figures/S1/v3/Fig_20022003_icelenses.png'
plt.savefig(fig_name,dpi=300)

print('Done with whole GrIS plot')
pdb.set_trace()

### --------------------------- BELOW: NOT USED --------------------------- ### 
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
