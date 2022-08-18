# Greenland Ice Sheet ice slab expansion and thickening

## Introduction

All the codes used for data preparation, processing, analysis and figure generation are present in this folder, as well as related text files and logbooks.
The way the codes and folders are now organized were not tested. Path adaptations are mandatory before running any code!
Some codes comes directly from MacFerrin et al., (2019), while some sections in other codes use function developped by MacFerrin et al., (2019), that may or not have been modified by us.
Feel free to contact the author for any questions (nicolas.jullien@unifr.ch).

If reproduction of 2010-2018 ice slabs is desired, here is the workflow to follow:
1. Run 'IceBridgeGPR_Manager_v2.py', for each year
2. Run 'iceslabs_processing_jullien.py'
3. Run 'probabilistic_ice_slabs.py'

If e.g. Fig. 1 generation is desired:
4. Run extract_elevation.py'
5. Run Fig1.py

## Codes
### Data processing

#### General

'Downdload_AR_data.py': Download accumulation radar data. May require a text file as input giving a list of data to download. All the dates with individual files processed in this study are listed in 'datetrack_20102018.txt'.

'aggregate_20022018_flightlines.py': Flighlines aggregation from 2002 to 2018. Note: the file 'metadata_coor_2002_2003' is an intermediate product and is not provided in the data repository.
'extract_2010_flightlines.py': Flighlines extractions for 2010

'extract_elevation.py': Extract elevation of ice slabs and its region from Rignot and Mouginot et al., (2012), for Ice_Layer_Output_Thicknesses_2010_2018_jullienetal2021_highend/lowend.csv and metadata_coord_icelens_2002_2003_26022020.pickle using Greenland DEM 100m (elevation above ellipsoid WGS84).
The resulting outputs are df_20102018_with_elevation_highend/lowend_rignotetalregions.pickle and df_2002_2003_with_elevation_rignotetalregions.pickle, but are not provided in the data repository because they are intermediate products.

'Identification_SurfacePixel_start.py': Help in manual identification of the surface pixel at the start of a track.

#### 2002-2003

'Retreive_coordinates_AR_2002_03.py': Matching between coordinates in the master files and individual radargrams. Each individual has a quality file associated with the matching. Here is the correspondance between the headers and their signification (from 'Retreive_Coordinates_AR_2002_2003.py'):
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

'Manual_SurfacePicking_WhereFail.py': Where the automatic procedure of surface identification failed, the surface was picked manually running this function. The resulting files ('surf_MMMDD_YY_ID_aggregated.txt') storing the surface pixel at any trace of the track was created, and used later on for data display.

'Signal_Bounds_2002_03.py': Derivation of the 2.5th and 97.5th percentiles of the distribution of the 30m thick potentially holding ice slabs radar tracks for each year in 2002-03, used for 2002-03 data display.

'Data_visualisation_2002_2003_AR.py': In this code, the visualising of the whole 2002-03 dataset can be done under many different ways.

'picking_icelenses_2002_2003.py': In this code, manually identification of ice lenses in 2002-03 data was performed by clicking on the figures and saving the resulting coordinates in pixel space.
Output file storing this information: icelenses_22022020.csv

#### 2010-2018

'rename_2012.py': Used to rename all the 2012 .mat files for data processing. 

'refine_location_2017_2018.py': Once all the 2017-2018 data were downloaded, this code was run to shorten the number of file where to look for ice slabs in 2017-18. For each individual track, it was checked if one point was contained within the iceslabs are shapefile from MacFerrin et al., (2019). If yes, the corresponding track was kept.
Then, the previous and following individual track of the corresponding suite of individual tracks was manually added to the file so that they are processed to make sure any inland expansion of ice slabs were previously not identified can be catched.

'initial_data_selection_2017_2018.py': Helper function to select 2017-18 data. Similar as 'refine_location_2017_2018.py'.

'delete_AR.py': When tracks are not listed in a text file (e.g. data_2017_toberun.txt), they are deleted in the folder where data are stored. The file 'data_2017_toberun.txt' can be created from 'datetrack_20102018.txt'.

'Investigation_time_issue.py': Allowed to find the issue related to the fail in the processing of some 2012-2013 data.

The following files are from MacFerrin et al., (2019). They were updated from Python 2 to Python 3. Some modifications were made on them, especially on 'IceBridgeGPR_Manager_v2.py'. For further details on these codes, please refer to: https://github.com/mmacferrin/Greenland_Ice_Slabs
'GPR_FileData.py'
'IceBridgeGPR_Manager_v2.py'
'IceBridgeRadarFlightlineDB.py'
'InSituGPR_Manager.py'

'identification_exclusions.py': Helper function to identify exclusions in 2010-2018 accumulation radar data. This function was partly used. The non-used parts use some text files which are not provided because they were not used in the final processing of the data.

'Quantiles_threshold_bestRange.py': Code used to identify the best range of quantile as threshold to differentiate between porous firn and ice from the appliance of the manual mask defined on the track 20130409_01_010_012.

'iceslabs_processing_jullien.py: Identification of ice with the appliance of the 15 quantiles as threshold to differentiate between porous firn and ice in 2010-2018 data.
This code reads as input the roll corrected tracks from 'IceBridgeGPR_Manager_v2.py' from MacFerrin et al., (2019), and performs surface extraction, depth correction without normalisation of the signal.
For tracks having too little signal variation, a rescaling of the signal was performed by setting the minimum and maximum signal strenght value to be the 5th and 95th percentiles from the distribution of the remaining tracks.
Then for each track: appliance independently of the 15 quantile to differentiate between porous firn and ice, appliance of smoothing and skrink and growth algorithm to the ice product.
Outputs are for each track 15 ice product (not provided in the data repository). 

'probabilistic_ice_slabs.py': Computation of final ice slabs in 2010-2018 data. Essentially, this code reads for each track the 15 ice slab identification products for each quantile, aggregate them together and generate a pickle and an image file of ice likelihood.
In a cell, if ice was identified in every quantile, its corresponding likelihood is 1; if it was identified only in the last quantile, its likelihood is 1/15.
Exclusion of areas (uses the file) where porous firn was erroneously identified as ice was performed at the end, before generating the individual files.
Uses as inputs the list of dates to process 'datetrack_20102018.txt' and the porous firn areas to exclude 'dry_firn_removal.txt'

The outputs are as such: 'YYYYMMDD_ID_START_END_probability_iceslabs_presence_after_DF.png/pickle'

The dataset storing all the data together is generated and names 'Ice_Layer_Output_Thicknesses_2010_2018_jullienetal2021_low/high_estimate.csv'

### Figures

All the figures can be regenerated using the corresponding code. Fig S5 can be generated using the corresponding section in Fig1.py. Fig 4, S6 and S7 are generated using the code Fig4andS6andS7.py, Fig S3 and S4 are both generated using the code FigS3andS4.py.
Individual pickle files are needed sometimes (e.g. in Fig4andS6andS7.py, FigS3and_S4.py), but they are not provided as they are an intermediate product. Instead, corresponding individual radargrams are avalaible in the data repository. Feel free to contact the author is the individual pickle files are desired.  

## Datasets

Our datasets are available at https://zenodo.org/record/6981378.

###Inputs

####2002-2003

Raw 2002-2003 data can be downloaded on the CReSIS data repository: https://data.cresis.ku.edu/data/accum/old_format/
A matching between the individual radargrams and the master GPS file for each folder was done based on a common timestamp (this is done in 'Retreive_coordinates_AR_2002_03.py'). The resulting matched data are found in: EnterFolderNameInDataRepository
Path: MatchedRadargrams_20022003.zip in data repository

'SURFACE_STARTING_PICKS_Suggestions_2002_2003.txt': File storing the start surface pixel of each track in 2002_03. Path: data/HelperFiles_20022003 in code repository.
'surf_picking_working_boolean_20022003_after_improvement.txt': File indicating whether an individual file of surface pixel for a specific track must be loaded for surface picking. Path: data/HelperFiles_20022003 in code repository.
'surf_MMMDD_YY_ID_aggregated.txt': Individual file storing the surface pixels of the whole track after manual foricing. Path: data/HelperFiles_20022003 in code repository.

####2010-2018:

Manual mask differentiating porous firn and ice defined on the track 20130409_01_010_012_XDEPTHCORRECT_AFTER.png, mask being binary_conservative_mask_20130409_01_010_012_XDEPTHCORRECT_AFTER.png
Path: AssociatedFiles.zip in data repository

###Outputs

####2002-2003
'2002_2003_green_excel.csv': Ice layers and slabs coordinates identified in the subsurface (from 0-30m deep).
Path: IceSlabs_20022018.zip in data repository

'ice_lenses_identification_MMDD_YY_ID_aggregated.png': Images, radargrams showing the ice slabs identification together with their traffic light classification code. The traffic light classification related to the confidence of mapping at the moment of identification is specified on the 2nd line of each sheet (green = high confidence, orange = medium confidence, red = low confidence, purple = near-surface refreezing layer). We kept only green identifications in the analysis.
Path: IndividualIceSlabsIdentification_20022018.zip in data repository

'metadata_coord_icelens_2002_2003_26022020.pickle': Ice layers and slabs coordinates identified in the subsurface (from 0-30m deep) together with their traffic light classification code. Can be build in ice_lenses_identification_final.py
Path: AssociatedFiles.zip in data repository

'metadata_coord_2002_2003.pickle': All 2002-2003 flightlines.
Path: Not provided, can be build in ice_lenses_identification_final.py.

####2010-2018
'referece_dry_firn_distrib.pickle': Porous firn radar signal strenght distribution after depth correction of the 20130409_01_010_012 track after appliance of the manual mask
'referece_iceslabs_distrib.pickle': Ice radar signal strenght distribution after depth correction of the 20130409_01_010_012 track after appliance of the manual mask
'quantile_file_0.65_0.79_onfulltrace.txt': Value of quantiles (from 0.65 to 0.79) of the ice radar signal strength distribution extracted from the appliance of the manual mask on the 20130409_01_010_012 track after depth correction
Path: AssociatedFiles.zip in data repository

#### Flightlines 2002-2018

'flightlines_20022018_GrIS.csv': All the flightlines intersecting with the GrIS (mask from Rignot et al., 2012) in 2002, 2003, 2010, 2011, 2012, 2013, 2014, 2017, 2018.
Path: Flightlines_20022018.zip in data repository

------------------------- To adapt still ------------------------- 
!! Add description of the modifications I have made !!

# Modifications of the files by Nicolas Jullien, starting on July 29th 2020

 July 29th 2020:
1. Moved from version 2.7 to python 3 the follwing codes:
	- FirnCore_Manager.py
	- GPR_Coordinates.py
	- GPR_Data.py
	- IceBridgeGPR_Manager_v2.py
	- InSituGPR_Manager.py
	- RCM_FileData.py
	- RCM_Manager.py
2. The file GPR_FileData.py have been into the changing version process 
but anything changed as no difference were spottable by the process
3. The file IceBridgeRadarFlightlineDB.py have been through the changing version process but it did not work. This file have been manually updated and several functions have been commented. This version have been tested and it worked in generating the database for the test data subset.
4. The filenames that are saved in 'self.FILENAMES', line 747 in the function _read_metadata() in the file 'IceBridgeGPR_Manager_v2.py' is now fixed.
5. Working in the branch 'debugging_with_subset' to debug the code working only with a subset of the data
# Greenland_Ice_Slabs
Code associated with detecting and modeling ice slabs in Greenland firn.

Note: This code was written between 2015-2019 by Michael J. MacFerrin. Committed to GitHub on June 20, 2019, upon acceptance of manuscript. The code served its purpose when it was written an initially executed (in stages as the project progressed). It has not been exhaustively tested for functionality, there may be parts of it that fail if used in ways that it was not intended, or on files that do not fit the exact assumptions of when it was written. If you would like to use the code but are running into problems, feel free to contact the author at michael.macferrin@colorado.edu.

Some datasets used in the paper are archived at the paper's [FigShare Repository](https://figshare.com/account/home#/projects/47690)

## License
This code is distributable and reusable under the GNU General Public License v3.0. It may be used or distributed freely. If it is used in a publication (academic or otherwise), we request that the following paper be cited:

MacFerrin, M., Machguth, H., van As, D., Charalampidis, C., Stevens, C. M., Heilig, A., 
Vandecrux, B., Langen, P. L., Mottram, R., Fettweis, X., Van den Broeke, M. R. , 
Pfeffer, W.T., Moussavi, M., Abdalati, W. (2019) "Rapid expansion of Greenland's low-
permeability ice slabs". Nature.

## Files (in 'src')
* **GPR_Data.py** : Contains file locations for data files used by the in-situ and IceBridge processing scripts.
* **GPR_Coordinates.py** : Used for processing and resampling GPR Coordinates (.cor) files.
* **GPR_Data.py** : Code used for processing steps on in-situ GPR data. NOTE: Much of the functionality in this script was deprecated and has since been replaced by "InSituGPR_Manager.py." However, some of this code was still used in final processing of the GPR data, so the file remains.)
* **InSituGPR_Manager.py** : Code for processing in-situ GPR data in Greenland for this manuscript.
* **IceBridgeGPR_Manager_v2.py** : Code for processing IceBridge airborne radar data in this manuscript. (Note: v1 code has all been deprecated, is not used in the final results, and is not included here.)
* **RCM_FileData.py** : Contains file locations for data files used by the RCM processing code (see below).
* **RCM_Manager.py** : Code for processing regional climate model (RCM) data for use in this manuscript.


------------------------- To adapt still ------------------------- 