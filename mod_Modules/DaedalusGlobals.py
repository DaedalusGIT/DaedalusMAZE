'''
This file contains some global variables which the modules can use.
You can acess it after you type at your code: import SourceCode.DaedalusGlobals as DaedalusGlobals
In case you need some variable to be added here, please contact the Daedalus team.
'''



# Absolute path where the files with the orbit data reside.
Orbit_Files_Path  = "/home/NAS/Data_Files/OrbitData/"
CSV_Files_Path = "/home/NAS/Data_Files/OrbitData/" 
NetCDF_Files_Path = "/home/NAS/Data_Files/OrbitData/"

# Absolute parent path where all data files are stored
AllData_Files_Path = "/home/NAS/"

# Absolute path where the TIEGCM data files are stored
TIEGCM_Files_Path = "/home/NAS/TIEGCM_DATA/"
TIEGCM_Files_Path_2 = "/home/NAS/TIEGCM_DATA_2/TIEGCM_Lifetime_2015_to_2018_JH_QD/"

# Absolute path where the WACCAM-X data files are stored
WACCAMX_Files_Path = "/home/NAS/WACCAM_X/"


# Absolute path where the files with Geomagnetic Kp-Ap Indices are stored
GeomagneticIndices_Files_Path = "/home/NAS/Data_Files/GeomagneticIndices/"

# Absolute path where the files with Geomagnetic Kp-Ap Indices are stored
CoverageResults_Files_Path = "/home/NAS/Data_Files/CoverageResults/"


# Absolute path where the filenames with calculated values reside.
Processed_Files_Path = "/home/NAS/Data_Files/ModelOutputs/" 

# Absolute path where temporary filenames 
Temporary_Files_Path = "/home/NAS/Data_Files/Temp/" 


# Absolute path where the filenames with the interpolated values reside.
Interpolated_Files_Path = "/home/NAS/Data_Files/ModelOutputs/A_Interpolated/" 

# Absolute path where the filenames with the Synthetic Truth values reside.
# That is after Sub-Grid Variability is applied on the data.
SyntheticTruth_Files_Path = "/home/NAS/Data_Files/ModelOutputs/B_SyntheticTruth/" 

# Absolute path where the filenames with the Stimuli values reside.
# That is after Wake-Effects and Charging-Effects are applied on the data.
Stimuli_Files_Path = "/home/NAS/Data_Files/ModelOutputs/C_Stimuli/" 

# Absolute path where the filenames with the Synthetic Science values reside.
# That is after the results of the satellite instruments are applied on the data.
SyntheticScience_Files_Path = "/home/NAS/Data_Files/ModelOutputs/D_SyntheticScience/" 
