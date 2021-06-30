# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 17:25:03 2021

@author: JullienN
"""

from IPython.display import Image
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from PIL import Image
from os import listdir
from os.path import isfile, join
import os
import pdb

##############################################################################
################### Define function for ice lenses logging ###################
##############################################################################
#This function if adapted from https://stackoverflow.com/questions/37363755/python-mouse-click-coordinates-as-simply-as-possible
def onclick(event):
    #This functions print and save the x and y coordinates in pixels!
    print(event.xdata, event.ydata)
    #Fill in the file to log on the information
    #filename_flog='C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/2002_2003_radar_slice/flog_icelenses_alldates.txt'
    #f_log = open(filename_flog, "a")
    #f_log.write(str(round(event.xdata,2))+','+str(round(event.ydata,2))+'\n')
    #f_log.close() #Close the quality assessment file when we’re done!
##############################################################################
################### Define function for ice lenses logging ###################
##############################################################################


#Define the working environment
path= 'C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/exported/'

f = open('C:/Users/jullienn/Documents/working_environment/iceslabs_MacFerrin/data/Exclusion_folder/data_2017_toberun.txt','r')
dates_surf_2018 = [line.strip() for line in f.readlines() if len(line.strip()) > 0]
f.close()

for indiv_file in list(dates_surf_2018):
    print(indiv_file)
    #Open only the roll_corrected files of interest
    #filename=indiv_file+'_XDEPTHCORRECT_AFTER.png'
    
    filename=indiv_file+'_0m_30m_BESTFIT_V1.png'
    
    #Open and plot image
    ROLLCORRECT__BEFORE = Image.open(path+filename).convert("L")
    arr_ROLLCORRECT__BEFORE = np.asarray(ROLLCORRECT__BEFORE)
    #Prepare plot$
    
    fig, (ax1) = plt.subplots(1, 1)#, gridspec_kw={'width_ratios': [1, 3]})
    fig.suptitle(indiv_file)
    ax1.imshow(arr_ROLLCORRECT__BEFORE, cmap='gray', vmin=0, vmax=255)
    fig.canvas.mpl_connect('button_press_event', onclick)
    plt.show()
    
    pdb.set_trace()

