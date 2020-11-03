# Daedalus_Interpolator_V2.0
# ********************************************************************************
# This program reads Daedalus' orbit (LLA) and TIEGCM's data and then interpolates
# TIGCM's variables (scalars) to Daedalus' orbit. In charge of the interpolation are
# a set of FORTRAN subroutines wrapped under f2py3. There are 2 different versions of
# those subroutines, Serial_Version and Parallel_Version. The Parallel_Version is
# using the OpenMP interface for multithreading.Interpolation Method: Trilinear 3D
# Instructions for FORTRAN subroutines compilation

# Serial_Version
# No Optimization: "f2py3 -c Interpolation_Serial_Version.f95 -m interpolation"

# Parallel_Version
# No Optimization: "f2py3 -c --f90flags="-fopenmp " -lgomp  Interpolation_Parallel_Version.f95 -m interpolation"
# -O1 Level Optimization: "f2py3 -c --f90flags="-fopenmp " -lgomp --opt='-O1' Interpolation_Parallel_Version.f95 -m interpolation"
# -O2 Level Optimization: "f2py3 -c --f90flags="-fopenmp " -lgomp --opt='-O2' Interpolation_Parallel_Version.f95 -m interpolation"
# -O3 Level Optimization: "f2py3 -c --f90flags="-fopenmp " -lgomp --opt='-O3' Interpolation_Parallel_Version.f95 -m interpolation"
# -ffast Level Optimization: "f2py3 -c --f90flags="-fopenmp " -lgomp --opt='-ffast' Interpolation_Parallel_Version.f95 -m interpolation"
# KEP 5/6/19


# Usage Instructions:
# 1) Compile either parallel or serial version selecting proper options from above
# 2) Run python code
# Compilation of fortran subroutines needs to take place once.
# ********************************************************************************



            #  ____                 _       _
            # |  _ \  __ _  ___  __| | __ _| |_   _ ___
            # | | | |/ _` |/ _ \/ _` |/ _` | | | | / __|
            # | |_| | (_| |  __/ (_| | (_| | | |_| \__ \
            # |____/ \__,_|\___|\__,_|\__,_|_|\__,_|___/



from netCDF4 import Dataset
import matplotlib.pyplot as plt
import pandas as pd
import SourceCode.ModulesSourceCode.Interpolator.interpolation as interpolation 
import numpy as np
import pandas as pd
from pandas import DataFrame
from datetime import datetime
import os

def Interpolator(model,model_data_file,orbit_file,save,VAR):
# model-->string:: name of model eg "TIEGCM"
# model_data_file--> string:: model data files in data folder netCDF file provided from modellers
# orbit_file-->string :: orbit filename in Time-Lat-Lon-Alt format as porvided in Jupyter
# save--> Logical:: if true saves interpolated values to directory
# VAR--> string:: variable to inteprolate, must be  one of the variables included in the model data
#Outputs:: Plots+ 1 csv file


    # ******************************************************************************

    # Interpolation Site Min and Max Altitude Along Track Depending on Model and User Preference
    min_alt=130
    max_alt=470                            

    #============== Read File + Orbitt and Manipulate Data=======
    
    TIEGCM=Dataset("../../NAS/TIEGCM_DATA/"+model_data_file)
    df = pd.read_csv(orbit_file)

    # =================================================
    
#   Get Orbit Data from Dataset
    daed_lat_temp = df["Lat_GEOD(deg)"]
    daed_lon_temp = df["Lon_GEOD(deg)"]
    daed_alt_temp = df["Height_WGS84 (km)"]
    daed_time_temp = df["Epoch(UTCG)"]


    # Sample Orbit and Find Positions Below max alt!
    counter=0
    for i in range(0,len(daed_alt_temp)):
        if (daed_alt_temp[i] < max_alt and daed_alt_temp[i] > min_alt):
            counter=counter+1
            
    # ******************************************************************************

    

    # ******************************************************************************
   
    # Allocate Arrays
    daed_lat=np.zeros((counter))
    daed_lon=np.zeros((counter))
    daed_alt=np.zeros((counter))
    daed_time=np.zeros((counter),dtype=datetime)


    #Keep Data in MIN MAX Range
    counter=0
    for i in range(0,len(daed_alt_temp)):
        if (daed_alt_temp[i] < max_alt and daed_alt_temp[i] > min_alt):
            daed_time[counter]=daed_time_temp[i]
            daed_lat[counter]=daed_lat_temp[i]
            daed_lon[counter]=daed_lon_temp[i]
            daed_alt[counter]=daed_alt_temp[i]
            counter=counter+1

    #Transfrom Lon to Match TIEGCM (-180-180 deg)
    for i in range (0, len(daed_alt)):
            daed_lon[i]=daed_lon[i]-180
    
    # ******************************************************************************


    # ******************************************************************************
    
    # Get data from Model Dataset
    grid_lat=TIEGCM.variables['lat'][:]
    grid_lon=TIEGCM.variables['lon'][:]
    grid_lev=TIEGCM.variables['ilev'][:]
    grid_time=TIEGCM.variables['time'][:]
    zg=TIEGCM.variables['ZG'][:]

    # Select Variable to Interpolate
    ne=TIEGCM.variables[VAR][:]
    
    # For surface plotting during debugging
    surf_Z=ne[1,10,:,:]     
    
    # ******************************************************************************


    
    # ******************************************************************************
    
    #Find model's temporal resolution
    model_dt=(grid_time[1]-grid_time[0])*60   #in seconds
    orbit_dt=1/16                             #in seconds
    #print("Timestep of Model Data is ",model_dt, "seconds")
    #print("Timestep of Orbit Data is ",orbit_dt, "seconds")
    
    # ******************************************************************************


    
    # ******************************************************************************
    
    #Allocate a matrix for the interpolated values
    int_data=np.zeros(((counter)))
    print("Interpolation Starting...",counter, "Positions to Interpolate")
    #*********************Call Interpolation Subroutines***************************
    interpolation.daed_interp(grid_lat,grid_lon,grid_lev,daed_lat,daed_lon,daed_alt,zg,model_dt,orbit_dt,ne,int_data)
    # *****************************************************************************
    print("Interpolation Done!")

    # ******************************************************************************
    


    # ******************************************************************************

    #Export Data to NAS
    if save==True:
        directory = "../../NAS/Data_Files/ModelsOutput/Interpolation/"

        Exports={'Time (UTCG)':daed_time,'Lat (deg)':daed_lat,'Lon (deg)':daed_lon,
                        'Alt (km)':daed_alt,'Interpolated_Data':int_data }
        
        df = DataFrame(Exports, columns= ['Time (UTCG)', 'Lat (deg)','Lon (deg)','Alt (km)','Interpolated_Data'])

        
        
#     input_path=os.path.dirname(orbit_file)
        export_name="../../NAS/Data_Files/ModelsOutput/Interpolation/"+os.path.splitext(os.path.basename(orbit_file))[0]+"_"+model+ "_"+ VAR+".csv"
            
        export_csv = df.to_csv (export_name, index = None, header=True)
#       print( "Output Saved!","Path-->", directory+orbit_file+"_" + model+ "_"  +VAR + '.csv')
    # ******************************************************************************


    
    # ******************************************************************************
    # Plotting
    xaxis=np.arange(0,len(int_data))
    plt.figure(1)
    plt.scatter(xaxis/(1*60),int_data,c=daed_alt,s=1)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.title("TIEGCM:" +VAR+" Along Daedalus' Orbit")
    plt.xlabel("Time [Minutes]")
    plt.ylabel("Eletron Density [$cm^{-3}$]")
    clb = plt.colorbar()
    clb.set_label('Altitude [km]')
    plt.grid("True")


    # Surface Plot
    # plt.figure(2)
    # plt.imshow(surf_Z,cmap="rainbow")
    # plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    # plt.colorbar()


    plt.show()
    # ******************************************************************************

    return(export_name)



# ******************************************************************************

# Example Call
# model_data_file='tiegcm_s_24900.nc'
# orbit_file="DAED_ORB_Evt0_LLA_Per120_Lat00_Srt01Hz_Msc"
# save=True
# VAR="NE"
# model="TIEGCM"
# DaedalusInterpolator(model,model_data_file,orbit_file,save,VAR)