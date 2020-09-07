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
