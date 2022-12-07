# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 17:55:47 2022

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

import pickle
import numpy as np
import matplotlib.pyplot as plt
import pdb
from PIL import Image
import matplotlib.gridspec as gridspec
from pyproj import Transformer

#define path
data_path='C:/Users/jullienn/switchdrive/Private/research/RT1/masking_iceslabs/quantiles_threshold_application/'

#Open dryfirn and ice slabs distribution
f_dryfirn = open(data_path+'referece_dry_firn_distrib.pickle', "rb")
dryfirn_distrib = pickle.load(f_dryfirn)
f_dryfirn.close()

f_iceslabs = open(data_path+'referece_iceslabs_distrib.pickle', "rb")
iceslabs_distrib = pickle.load(f_iceslabs)
f_iceslabs.close()

#Open dataframe
f_dataframe = open(data_path+'20130409_01_010_012_dataframeforS2.pickle', "rb")
dataframe = pickle.load(f_dataframe)
f_dataframe.close()

#Open quantile 0.65
f_quant065 = open(data_path+'20130409_01_010_012_SG1_cutoff_0.65_threshold_350.pickle', "rb")
quant065 = pickle.load(f_quant065)
f_quant065.close()

#Open quantile 0.79
f_quant079 = open(data_path+'20130409_01_010_012_SG1_cutoff_0.79_threshold_350.pickle', "rb")
quant079 = pickle.load(f_quant079)
f_quant079.close()

#Open conservative mask
path_mask='C:/Users/jullienn/switchdrive/Private/research/RT1/masking_iceslabs/'
boolean_mask = Image.open(path_mask+'binary_conservative_mask_20130409_01_010_012_XDEPTHCORRECT_AFTER.png').convert("L")
arr_boolean_mask = np.asarray(boolean_mask)

#Keep only the first 20m of the 30m boolean mask
ind_20m=np.where(dataframe['depth']<20)[0]
arr_boolean_mask_20m=arr_boolean_mask[ind_20m,:]
# Convert the image mask into a boolean
arr_boolean_mask_20m=~arr_boolean_mask_20m
arr_boolean_mask_20m[arr_boolean_mask_20m==255]=1

#Finalize the mask
final_mask=np.zeros((arr_boolean_mask_20m.shape[0],arr_boolean_mask_20m.shape[1]))
final_mask[:]=np.nan
final_mask[arr_boolean_mask_20m==1]=arr_boolean_mask_20m[arr_boolean_mask_20m==1] 

#Finalize the quantile 0.65
final_quant065=np.zeros((quant065.shape[0],quant065.shape[1]))
final_quant065[:]=np.nan
final_quant065[quant065==1]=quant065[quant065==1]
#Finalize the quantile 0.79
final_quant079=np.zeros((quant079.shape[0],quant079.shape[1]))
final_quant079[:]=np.nan
final_quant079[quant079==1]=quant079[quant079==1]

#Define transformer for coordinates transform from "EPSG:4326" to "EPSG:3413"
transformer = Transformer.from_crs("EPSG:4326", "EPSG:3413", always_xy=True)

#Transform the coordinated from WGS84 to EPSG:3413
#Example from: https://pyproj4.github.io/pyproj/stable/examples.html
points=transformer.transform(np.array(dataframe['lon_appended']),np.array(dataframe['lat_appended']))
lon_3413=points[0]
lat_3413=points[1]

lon3413_plot=np.zeros(len(lon_3413))
lon3413_plot[:]=np.nan
lon3413_plot[dataframe['mask'][:,0]]=lon_3413[dataframe['mask'][:,0]]

lat3413_plot=np.zeros(len(lat_3413))
lat3413_plot[:]=np.nan
lat3413_plot[dataframe['mask'][:,0]]=lat_3413[dataframe['mask'][:,0]]

#Calculate distances
distances_with_start_transect=compute_distances(lon3413_plot,lat3413_plot)

#Create the figure
plt.rcParams.update({'font.size': 15})
fig = plt.figure(figsize=(8,21))
gs = gridspec.GridSpec(14, 51)
gs.update(wspace=0.001)
ax2 = plt.subplot(gs[0:1, 0:50])
ax3 = plt.subplot(gs[1:2, 0:50])
ax4 = plt.subplot(gs[2:3, 0:50])
ax5 = plt.subplot(gs[3:4, 0:50])
axc = plt.subplot(gs[0:4, 50:51])

#Display depth corrected radargram
ax2.pcolor(distances_with_start_transect, dataframe['depth'][ind_20m], dataframe['depth_corrected_after_surf_removal_without_norm'][ind_20m,:],cmap=plt.get_cmap('gray'),zorder=-2)#,norm=divnorm)
ax2.invert_yaxis() #Invert the y axis = avoid using flipud.    
ax2.set_ylim(20,0)
plt.setp(ax2.get_xticklabels(), visible=False)
ax2.set_yticks([0,10,20])
ax2.set_yticklabels(['0','10',''])
ax2.text(0.01, 0.75,'a',ha='center', va='center', transform=ax2.transAxes,fontsize=15,zorder=10,weight='bold',color='white')#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot
ax2.set_aspect(dataframe['depth_corrected_after_surf_removal_without_norm'][ind_20m,:].shape[1]/dataframe['depth_corrected_after_surf_removal_without_norm'][ind_20m,:].shape[0])

#Display depth corrected radargram and manual mask over it
ax3.pcolor(distances_with_start_transect, dataframe['depth'][ind_20m], dataframe['depth_corrected_after_surf_removal_without_norm'][ind_20m,:],cmap=plt.get_cmap('gray'),zorder=-2)#,norm=divnorm)
ax3.invert_yaxis() #Invert the y axis = avoid using flipud.    
ax3.set_ylim(20,0)
plt.setp(ax3.get_xticklabels(), visible=False)
ax3.set_yticks([0,10,20])
ax3.set_yticklabels(['0','10',''])
ax3.pcolor(distances_with_start_transect, dataframe['depth'][ind_20m], final_mask,cmap=plt.get_cmap('gray'),zorder=0)#,norm=divnorm)
ax3.text(0.01, 0.75,'b',ha='center', va='center', transform=ax3.transAxes,fontsize=15,zorder=10,weight='bold',color='white')#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot
ax3.set_aspect(dataframe['depth_corrected_after_surf_removal_without_norm'][ind_20m,:].shape[1]/dataframe['depth_corrected_after_surf_removal_without_norm'][ind_20m,:].shape[0])


#Display depth corrected radargram and quantile 0.65 over it
ax4.pcolor(distances_with_start_transect, dataframe['depth'][ind_20m], dataframe['depth_corrected_after_surf_removal_without_norm'][ind_20m,:],cmap=plt.get_cmap('gray'),zorder=-2)#,norm=divnorm)
ax4.invert_yaxis() #Invert the y axis = avoid using flipud.    
ax4.set_ylim(20,0)
ax4.set_yticks([0,10,20])
plt.setp(ax4.get_xticklabels(), visible=False)
ax4.set_yticklabels(['0','10',''])
ax4.pcolor(distances_with_start_transect, dataframe['depth'][ind_20m], final_quant065,cmap=plt.get_cmap('gray'),zorder=0)#,norm=divnorm)
ax4.text(0.01, 0.75,'c',ha='center', va='center', transform=ax4.transAxes,fontsize=15,zorder=10,weight='bold',color='white')#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot
ax4.set_aspect(dataframe['depth_corrected_after_surf_removal_without_norm'][ind_20m,:].shape[1]/dataframe['depth_corrected_after_surf_removal_without_norm'][ind_20m,:].shape[0])


#Display depth corrected radargram and quantile 0.79 over it
cb=ax5.pcolor(distances_with_start_transect, dataframe['depth'][ind_20m], dataframe['depth_corrected_after_surf_removal_without_norm'][ind_20m,:],cmap=plt.get_cmap('gray'),zorder=-2)#,norm=divnorm)
ax5.invert_yaxis() #Invert the y axis = avoid using flipud.    
ax5.set_ylim(20,0)
ax5.set_ylabel('Depth [m]')
ax5.set_yticks([0,10,20])
ax5.pcolor(distances_with_start_transect, dataframe['depth'][ind_20m], final_quant079,cmap=plt.get_cmap('gray'),zorder=0)#,norm=divnorm)
ax5.text(0.01, 0.75,'d',ha='center', va='center', transform=ax5.transAxes,fontsize=15,zorder=10,weight='bold',color='white')#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot
ax5.set_aspect(dataframe['depth_corrected_after_surf_removal_without_norm'][ind_20m,:].shape[1]/dataframe['depth_corrected_after_surf_removal_without_norm'][ind_20m,:].shape[0])

#Set distance
low_xlim=ax5.get_xlim()[0]
high_xlim=ax5.get_xlim()[1]
#Display bottom xtick in km instead of m
xtick_distance=ax5.get_xticks()
ax5.set_xticks(xtick_distance)
ax5.set_xticklabels((xtick_distance/1000).astype(int))
ax5.set_xlim(low_xlim,high_xlim)
ax5.set_xlabel('Distance [km]')
ax5.set_aspect(dataframe['depth_corrected_after_surf_removal_without_norm'][ind_20m,:].shape[1]/dataframe['depth_corrected_after_surf_removal_without_norm'][ind_20m,:].shape[0])


#Display colorbar. This is from FigS1.py
cbar_depth=fig.colorbar(cb, cax=axc)
cbar_depth.set_label('Radar signal strength [dB]')
plt.show()

pdb.set_trace()
'''
#Save figure
plt.savefig('C:/Users/jullienn/switchdrive/Private/research/RT1/figures/S4and5/v6/figS4_radargrams.png',dpi=300,bbox_inches='tight')
#bbox_inches is from https://stackoverflow.com/questions/32428193/saving-matplotlib-graphs-to-image-as-full-screen
'''
#Create the figure
fig,ax1 = plt.subplots()
plt.rcParams.update({'font.size': 15})

#Display histograms
ax1.hist(iceslabs_distrib,bins=500,density=True,label='In-situ ice content')
ax1.hist(dryfirn_distrib,bins=500,density=True,alpha=0.2,label='Porous firn')
ax1.legend()
ax1.set_xlabel('Radar signal strength [dB]',fontsize=15)
ax1.set_ylabel('Probability density [ ]',fontsize=15)
ax1.tick_params(axis='both', which='major', labelsize=15)

#Define the desired quantiles
desired_quantiles=np.arange(0.65,0.80,0.01)
#Define quantiles for investigation of accuracy
quantile_investigation=np.quantile(iceslabs_distrib,desired_quantiles)
#Add low and high quantile as dashed lines

ax1.axvline(x=quantile_investigation[0],linestyle='--',color='k')
ax1.axvline(x=quantile_investigation[-1],linestyle='--',color='k')
ax1.text(quantile_investigation[0]-0.025, 0.85, 'quantile 0.65', rotation=90, va='center',fontsize=15)
ax1.text(quantile_investigation[-1]+0.005, 0.85, 'quantile 0.79', rotation=90, va='center',fontsize=15)

ax1.set_xlim(-0.5,0.5)
'''
#Save figure
plt.savefig('C:/Users/jullienn/switchdrive/Private/research/RT1/figures/S4and5/v6/figS5_distribution.png',dpi=300,bbox_inches='tight')
#bbox_inches is from https://stackoverflow.com/questions/32428193/saving-matplotlib-graphs-to-image-as-full-screen
'''
pdb.set_trace()

#Display contingency table between mask, quantile 0.65 and quantile 0.79
from sklearn.metrics.cluster import contingency_matrix
#Replace NaNs by 0
final_mask[np.isnan(final_mask)]=0
final_quant065[np.isnan(final_quant065)]=0
final_quant079[np.isnan(final_quant079)]=0

#contingency matrix from https://scikit-learn.org/stable/modules/generated/sklearn.metrics.cluster.contingency_matrix.html#sklearn-metrics-cluster-contingency-matrix
contingency_matrix(final_mask,final_quant065) #contingency_matrix(labels_true, labels_pred)
contingency_matrix(final_mask,final_quant079) #contingency_matrix(labels_true, labels_pred)


#Display for each quantile the ice slabs identification
desired_quantiles=np.round(np.arange(0.60,0.91,0.01),2)

for single_quantile in desired_quantiles:
    print(single_quantile)
    
    #Open quantile
    f_indiv_quant = open(data_path+'20130409_01_010_012_SG1_cutoff_'+str(single_quantile)+'_threshold_350.pickle', "rb")
    indiv_quant = pickle.load(f_indiv_quant)
    f_indiv_quant.close()
    
    #Finalize the quantile
    final_indiv_quant=np.zeros((indiv_quant.shape[0],indiv_quant.shape[1]))
    final_indiv_quant[:]=np.nan
    final_indiv_quant[indiv_quant==1]=indiv_quant[indiv_quant==1]
        
    fig, (ax1, ax2) = plt.subplots(2)
    #Display depth corrected radargram and manual mask over it
    ax1.pcolor(distances_with_start_transect, dataframe['depth'][ind_20m], dataframe['depth_corrected_after_surf_removal_without_norm'][ind_20m,:],cmap=plt.get_cmap('gray'),zorder=-2)#,norm=divnorm)
    ax1.invert_yaxis() #Invert the y axis = avoid using flipud.    
    ax1.set_ylim(20,0)
    plt.setp(ax3.get_xticklabels(), visible=False)
    ax1.pcolor(distances_with_start_transect, dataframe['depth'][ind_20m], final_mask,cmap=plt.get_cmap('gray'),zorder=0)#,norm=divnorm)
    ax1.set_aspect(100)
    
    #Display depth corrected radargram and indiv_quantile over it
    ax2.pcolor(distances_with_start_transect, dataframe['depth'][ind_20m], dataframe['depth_corrected_after_surf_removal_without_norm'][ind_20m,:],cmap=plt.get_cmap('gray'),zorder=-2)#,norm=divnorm)
    ax2.invert_yaxis() #Invert the y axis = avoid using flipud.    
    ax2.set_ylim(20,0)
    plt.setp(ax4.get_xticklabels(), visible=False)
    ax2.pcolor(distances_with_start_transect, dataframe['depth'][ind_20m], final_indiv_quant,cmap=plt.get_cmap('gray'),zorder=0)#,norm=divnorm)
    ax2.text(0.01, 0.75,'d',ha='center', va='center', transform=ax4.transAxes,fontsize=15,zorder=10,weight='bold',color='white')#This is from https://pretagteam.com/question/putting-text-in-top-left-corner-of-matplotlib-plot
    ax2.set_aspect(100)
    
    #Save figure
    plt.savefig('C:/Users/jullienn/switchdrive/Private/research/RT1/figures/S4and5/v4/figS4_'+str(single_quantile)+'.png',dpi=300)
