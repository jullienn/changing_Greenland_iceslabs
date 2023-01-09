# Greenland Ice Sheet Ice Slab Expansion and Thickening

## Introduction

All the codes used for data preparation, processing, analysis and figure generation are present in this folder, as well as related helper files.
These codes were written by Nicolas Jullien between 2020-2022.
The way the codes and folders are now organized in this code repository were not tested, and some parts may fail. Furthermore, path adaptations are mandatory before running any code.
Some codes (that we may or may not have modified) were directly taken from MacFerrin et al., (2019). Some sections in other codes use functions developped by MacFerrin et al., (2019), that may or not have been modified by us.
Feel free to contact the author for any questions (nicolas.jullien@unifr.ch).

The datasets used in this paper, as well as exported data are archived in the data repository related to the paper, accessible at https://doi.org/10.5281/zenodo.7505426.

If reproduction of 2010-2018 ice slabs is desired, here is the workflow to follow:
1. Run 'IceBridgeGPR_Manager_v2.py', for each year
2. Run 'iceslabs_processing_jullien.py'
3. Run 'probabilistic_ice_slabs.py'

If e.g. Fig. 2 generation is desired:
4. Run extract_elevation.py'
5. Run Fig2andS7andS8andS12.py

## Codes
### Data processing

#### General

* 'Downdload_AR_data.py': Download accumulation radar data. May require a text file as input giving a list of data to download. All the dates with individual files processed in this study are listed in 'datetrack_20102018.txt'.

* 'aggregate_20022018_flightlines.py': Flighlines aggregation from 2002 to 2018. Note: the file 'metadata_coor_2002_2003' is an intermediate product and is not provided in the data repository.

* 'extract_2010_flightlines.py': Flighlines extractions for 2010

* 'extract_elevation.py': Extract elevation of ice slabs and its region from Rignot and Mouginot (2012), for Ice_Layer_Output_Thicknesses_2010_2018_jullienetal2023_highend/lowend.csv and metadata_coord_icelens_2002_2003_26022020.pickle using Greenland DEM 100m (elevation above ellipsoid WGS84).
The resulting outputs are df_20102018_with_elevation_highend/lowend_rignotetalregions.pickle and df_2002_2003_with_elevation_rignotetalregions.pickle, but are not provided in the data repository because they are intermediate products.

* 'Identification_SurfacePixel_start.py': Help in manual identification of the surface pixel at the start of a track.

#### 2002-2003

* 'Retreive_coordinates_AR_2002_03.py': Matching between coordinates in the master files and individual radargrams. Each individual has a quality file associated with the matching. Here is the correspondance between the headers and their signification (from 'Retreive_Coordinates_AR_2002_2003.py'):
#date: date of the file assessed
#B_match_dftimearr_dftimearrdec: correspondance in df_final between begin 'timearr' and begin 'timearr_dec' -> 1 is matching, 0 is no match
#E_match_dftimearr_dftimearrdec: correspondance in df_final between end 'timearr' and end 'timearr_dec' -> 1 is matching, 0 is no match
#B_match_dftimearrdec_filetimearr: correspondance between begin 'timearr_dec' in df_final and begin 'timearr' of df_file_being_read -> 1 is matching, 0 is no match
#E_match_dftimearrdec_filetimearr: correspondance between end 'timearr_dec' in df_final and end 'timearr' of df_file_being_read for which we have echogram data (defined by the size of filftin)-> 1 is matching, 0 is no match
#length_match_df_file: length of df_final and file_being_read['filtfin'] should be the same -> 1 is yes, 0 is no
#B0to2_df_timearrdec_df_secondsgps: the first 3 rows of floor(df_final['timearr_dec']) and df_final['seconds_gps'] must be identical. If they are, the value 3 should be stored
#Em1tom3_df_timearrdec_df_secondsgps: the last 3 rows of floor(df_final['timearr_dec']) and df_final['seconds_gps'] must be identical. If they are, the value 3 should be stored
#file_being_read[filtfin].shape[1]: size of file_being_read[filtfin].shape[1]
#df_final.shape[0]: size of df_final.shape[0]

* 'Manual_SurfacePicking_WhereFail.py': Where the automatic procedure of surface identification failed, the surface was picked manually running this function. The resulting files ('surf_MMMDD_YY_ID_aggregated.txt') storing the surface pixel at any trace of the track was created, and used later on for data display.

* 'Signal_Bounds_2002_03.py': Derivation of the 2.5th and 97.5th percentiles of the distribution of the 30m thick potentially holding ice slabs radargram for each year in 2002-03, used for 2002-03 data display.

* 'Data_visualisation_2002_2003_AR.py': In this code, the visualising of the whole 2002-03 dataset can be done under many different ways.

* 'picking_icelenses_2002_2003.py': In this code, manually identification of ice lenses in 2002-03 data was performed by clicking on the figures and saving the resulting coordinates in pixel space.
Output file storing this information: icelenses_22022020.csv

#### 2010-2018

* 'rename_2012.py': Used to rename all the 2012 .mat files for data processing. 

* 'refine_location_2017_2018.py': Once all the 2017-2018 data were downloaded, this code was run to shorten the number of file where to look for ice slabs in 2017-18. For each individual track, it was checked if one point was contained within the iceslabs area shapefile from MacFerrin et al., (2019). If yes, the corresponding track was kept.
Then, the previous and following individual track of the corresponding suite of individual tracks was manually added to the file so that they are processed to make sure any inland expansion of ice slabs can be catched.

* 'initial_data_selection_2017_2018.py': Helper function to select 2017-18 data. Similar as 'refine_location_2017_2018.py'.

* 'delete_AR.py': When tracks are not listed in a text file (e.g. data_2017_toberun.txt), they are deleted in the folder where data are stored. The file 'data_2017_toberun.txt' can be created from 'datetrack_20102018.txt'.

* 'Investigation_time_issue.py': Allowed to find the issue related to the fail in the processing of some 2012-2013 data.

The following files are from MacFerrin et al., (2019). They were updated from Python 2.7 to Python 3. Some modifications were made on them, especially on 'IceBridgeGPR_Manager_v2.py'. For further details on these codes, please refer to: https://github.com/mmacferrin/Greenland_Ice_Slabs
* 'GPR_FileData.py': Code for filenames definition and organisation.
* 'IceBridgeGPR_Manager_v2.py': Code for processing 2010-2017 accumulation radar data (reading data, surface identification, obivous exclusions appliance, correction for the roll of the aricraft). Note that some fuctions that were developped in this code by MacFerrin et al., (2019) were also used in the processing of the 2010-2018 data at a later stage, and for the 2002-2003 data processing (sometimes with some adaptations). Some modifications were made to this code, e.g. :
- Sometimes, the start of the time variable was not zero (either too positive, or negative). Because the time variable increases linearly, we reset the start of the time variable to zero when it was not the case, increasing towards positive values.
- Some little modifications were made from times to times, probably because of the moving of the code from Python 2.7 to Python 3, or maybe because the processing chain is a bit different compared to how the code was used by MacFerrin et al., (2019). For example, line 2911-2924, we added to load in .self the boolean traces after they have been generated. This was not done previously and self.FNAME_ice_lenses_picklefile was empty. It was apparently not needed in MacFerrin's version, although this is essential here as self.TRACES_boolean_ice_layers is used to extract the ice content to be stored in the csv file.
* 'IceBridgeRadarFlightlineDB.py': Code for generating the 2010-2018 accumuluation radar database, later used by 'IceBridgeGPR_Manager_v2.py'. The automatic conversion of this file from Python 2.7 to Python 3 did not work. This file was successfully manually updated and several functions were commented.
* 'InSituGPR_Manager.py': Code used for the definition of the speed of the radar signal in icy firn. The rest of the code was commented as not useful for our study.

* 'identification_exclusions.py': Helper function to identify exclusions in 2010-2018 accumulation radar data. This function was partly used. The non-used parts use some text files which are not provided because they were not used in the final processing of the data.

* 'Quantiles_threshold_bestRange.py': Code used to identify the best range of quantile as threshold to differentiate between porous firn and ice from the appliance of the manual mask defined on the transect 20130409_01_010_012.

* 'iceslabs_processing_jullien.py: Identification of ice with the appliance of the 15 quantiles as threshold to differentiate between porous firn and ice in 2010-2018 data.
This code reads as input the roll corrected tracks from 'IceBridgeGPR_Manager_v2.py' from MacFerrin et al., (2019), and performs surface extraction, depth correction without normalisation of the signal.
For transects having too little signal variation, a rescaling of the signal was performed by setting the minimum and maximum signal strenght value to be the 5th and 95th percentiles from the distribution of the remaining transects.
Then for each transect: appliance independently of the 15 quantile to differentiate between porous firn and ice, appliance of smoothing and skrink and growth algorithm to the ice product.
Outputs are for each transects fifteen different ice product (not provided in the data repository). 

* 'probabilistic_ice_slabs.py': Computation of final ice slabs in 2010-2018 data. Essentially, this code reads for each transect the fifteen (i.e. one for each quantile) ice slab identification products, aggregate them together and generate a pickle and an image file of ice likelihood.
In a cell, if ice was identified in every quantile, its corresponding likelihood is 1; if it was identified only in the last quantile, its likelihood is 1/15.
Exclusion of areas where porous firn was erroneously identified as ice was performed at the end, but before generating the individual files.
Uses as inputs the list of dates to process 'datetrack_20102018.txt' and the porous firn areas to exclude 'dry_firn_removal.txt'.

The outputs are as such: 'YYYYMMDD_ID_START_END_probability_iceslabs_presence_after_DF.png/pickle'

Finally, the dataset storing all the data together is generated and names 'Ice_Layer_Output_Thicknesses_2010_2018_jullienetal2023_low/high_estimate.csv'

### Figures

All the figures can be regenerated using the corresponding code.
Individual pickle files are needed sometimes (e.g. in Fig4andS6andS9.py, FigS3and_S4.py), but they are not provided as they are an intermediate product. Instead, corresponding individual radargrams (as .png) are avalaible in the data repository. Feel free to contact the author is the individual pickle files are desired.  

## Datasets

Our datasets are available at https://doi.org/10.5281/zenodo.7505426

### Inputs

#### 2002-2003

Raw 2002-2003 data can be downloaded on the CReSIS data repository: https://data.cresis.ku.edu/data/accum/old_format/
A matching between the individual radargrams and the master GPS file for each folder was done based on a common timestamp (this is done in 'Retreive_coordinates_AR_2002_03.py').
Path to data : MatchedRadargrams_20022003/raw_radargrams in data repository.

* 'SURFACE_STARTING_PICKS_Suggestions_2002_2003.txt': File storing the start surface pixel of each track in 2002_03. Path: data/HelperFiles_20022003 in code repository.
* 'surf_picking_working_boolean_20022003_after_improvement.txt': File indicating whether an individual file of surface pixel for a specific track must be loaded for surface picking. Path: data/HelperFiles_20022003 in code repository.
* 'surf_MMMDD_YY_ID_aggregated.txt': Individual file storing the surface pixels of the whole track after manual foricing. Path: data/HelperFiles_20022003 in code repository.

#### 2010-2018:

Manual mask differentiating porous firn and ice defined on the track 20130409_01_010_012_XDEPTHCORRECT_AFTER.png, mask being binary_conservative_mask_20130409_01_010_012_XDEPTHCORRECT_AFTER.png
Path: AssociatedFiles/Mask20130409_01_010_012 in data repository.

### Outputs

#### 2002-2003

* '2002_2003_green_excel.csv': Ice layers and slabs coordinates identified in the subsurface (from 0-30m deep).
Path: IceSlabs_20022018/20022003 in data repository.

* 'ice_lenses_identification_MMDD_YY_ID_aggregated.png': Images, radargrams showing the ice slabs identification together with their traffic light classification code.
Path: IndividualIceSlabsIdentification_20022018/20022003 in data repository.

* 'metadata_coord_icelens_2002_2003_26022020.pickle': Ice layers and slabs coordinates identified in the subsurface (from 0-30m deep) in 2002-2003 together with their traffic light classification code. We associated a traffic light classification related to the confidence of mapping at the moment of identification (green = high confidence, orange = medium confidence, red = low confidence, purple = near-surface refreezing layer). This traffic light classification is specified by the number store in the third column (-1 = red, 0 = orange, 1 = green, 2 = purple). Note that we kept only green identifications in the analysis. This file was generated running ice_lenses_identification_final.py
Path: AssociatedFiles in data repository.

* 'metadata_coord_2002_2003.pickle': All 2002-2003 flightlines.
Not provided, can be build running ice_lenses_identification_final.py. 2002-2003 flightlines are present in the file 'flightlines_20022018_GrIS.csv' (see further down).

#### 2010-2018
* 'referece_dry_firn_distrib.pickle': Porous firn radar signal strenght distribution after depth correction of the 20130409_01_010_012 transect after appliance of the manual mask. Path: AssociatedFiles/Mask20130409_01_010_012 in data repository.
* 'referece_iceslabs_distrib.pickle': Ice radar signal strenght distribution after depth correction of the 20130409_01_010_012 transect after appliance of the manual mask. Path: AssociatedFiles/Mask20130409_01_010_012 in data repository.
* 'quantile_file_0.65_0.79_onfulltrace.txt': Value of quantiles (from 0.65 to 0.79) of the ice radar signal strength distribution extracted from the appliance of the manual mask on the 20130409_01_010_012 transect after depth correction. Path: /data in the code repository.

#### Flightlines 2002-2018

'flightlines_20022018_GrIS.csv': All the flightlines intersecting with the GrIS (mask from Rignot and Mouginot (2012)) in 2002, 2003, 2010, 2011, 2012, 2013, 2014, 2017, 2018.
Path: Flightlines_20022018 in data repository.

## License
This code repository is under the GNU General Public License v3.0. If it is used in publication (academic or otherwise), we request the following paper be cited:
N. Jullien, A. J. Tedstone, H. Machguth, N. B. Karlsson, V. Helm (2023) "Greenland Ice Slab Expansion and Thickening,". Submitted to Geophysical Research Letters (in review),

and MacFerrin et al., (2019) when applicable:

MacFerrin, M., Machguth, H., van As, D., Charalampidis, C., Stevens, C. M., Heilig, A., Vandecrux, B., Langen, P. L., Mottram, R., Fettweis, X., Van den Broeke, M. R. , Pfeffer, W.T., Moussavi, M., Abdalati, W. (2019) "Rapid expansion of Greenland's low- permeability ice slabs". Nature.

