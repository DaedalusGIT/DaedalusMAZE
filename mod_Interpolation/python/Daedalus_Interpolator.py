# Daedalus_Interpolator_V5.0
# ******************************************************************************
# This program reads Daedalus' orbit (LLA) and TIEGCM's data and then interpolates
# TIEGCM's variables (scalars) to Daedalus' orbit. In charge of the interpolation are
# a set of FORTRAN subroutines wrapped under f2py3. There are 2 different versions of
# those subroutines, Serial_Version and Parallel_Version. The Parallel_Version is
# using the OpenMP interface for multithreading.Interpolation Method: Trilinear 3D
# Instructions for FORTRAN subroutines compilation
# Update: Python Interpolation Routines both serial and parallel
# Serial_Version
# No Optimization: "f2py3 -c Interpolation_Serial.f90 -m interpolation"

# Parallel_Version
# No Optimization: "f2py3 -c --f90flags="-fopenmp " -lgomp  Interpolation_Parallel.f90 -m interpolation"
# -O1 Level Optimization: "f2py3 -c --f90flags="-fopenmp " -lgomp --opt='-O1' Interpolation_Parallel.f90 -m interpolation_OpenMP"
# -O2 Level Optimization: "f2py3 -c --f90flags="-fopenmp " -lgomp --opt='-O2' Interpolation_Parallel.f90 -m interpolation_OpenMP"
# -O3 Level Optimization: "f2py3 -c --f90flags="-fopenmp " -lgomp --opt='-O3' Interpolation_Parallel.f90 -m interpolation_OpenMP"
# -ffast Level Optimization: "f2py3 -c --f90flags="-fopenmp " -lgomp --opt='-ffast' Interpolation_Parallel.f90 -m interpolation_OpenMP"
# KEP 5/6/19


# Usage Instructions:
# 1) Compile either parallel or serial version selecting proper options from above
# 2) Run python code
# Compilation of FORTRAN subroutines needs to take place once.
# ******************************************************************************



            #  ____                 _       _
            # |  _ \  __ _  ___  __| | __ _| |_   _ ___
            # | | | |/ _` |/ _ \/ _` |/ _` | | | | / __|
            # | |_| | (_| |  __/ (_| | (_| | | |_| \__ \
            # |____/ \__,_|\___|\__,_|\__,_|_|\__,_|___/



from netCDF4 import Dataset
import matplotlib.pyplot as plt
import pandas as pd
import interpolation
import interpolation_OpenMP
import numpy as np
import pandas as pd
from pandas import DataFrame
from datetime import datetime
import time
import os
import multiprocessing
from functools import partial


def local(dim1,y,x):
    local_pos=0
    for i in range(0,len(x)-1):
        if y >= x[i] and y < x[i+1]:
            local_pos=i
            break
    return (local_pos)

def Interpolate_Parallel(grid_lat,grid_lon,grid_lev,daed_lat,daed_lon,daed_alt,zg,model_dt,orbit_dt,ne,i):
    counter=1
    deltaphi= np.abs(grid_lon[2]-grid_lon[1])
    deltatheta=np.abs(grid_lat[2]-grid_lat[1])

    if i >-1:


        phi_local=local(len(grid_lon),daed_lon[i],grid_lon)
        theta_local=local(len(grid_lat),daed_lat[i],grid_lat)

        if daed_lon[i] >= 177.5 :
            phi_local=len(grid_lon)-2

        if daed_lon[i] <= -177.5:
            phi_local=0

        if daed_lat[i] >= 87.5:
            theta_local=len(grid_lat)-2

        if daed_lat[i] <= -87.5:
            theta_local=0


        alts=zg[counter,:,theta_local,phi_local]/100000
        r_local=local(len(alts),daed_alt[i],alts)

        deltarho=alts[r_local+1]-alts[r_local]

        dx=np.abs(((daed_alt[i]-alts[r_local])/deltarho))
        dy=np.abs(((daed_lat[i]-grid_lat[theta_local])/deltatheta))
        dz=np.abs(((daed_lon[i]-grid_lon[phi_local])/deltaphi))

        w1=(1-dx)*(1-dy)*(1-dz)
        w2=(dx)*(1-dy)*(1-dz)
        w3=(1-dx)*(dy)*(1-dz)
        w4=(dx)*(dy)*(1-dz)
        w5=(1-dx)*(1-dy)*(dz)
        w6=(dx)*(1-dy)*(dz)
        w7=(1-dx)*(dy)*(dz)
        w8=(dx)*(dy)*(dz)

        m=0.0
        m=       ne[counter,r_local,theta_local,phi_local]*w1
        m=m+  ne[counter,r_local+1,theta_local,phi_local]*w2
        m=m+  ne[counter,r_local,theta_local+1,phi_local]*w3
        m=m+  ne[counter,r_local+1,theta_local+1,phi_local]*w4
        m=m+  ne[counter,r_local,theta_local,phi_local+1]*w5
        m=m+  ne[counter,r_local+1,theta_local,phi_local+1]*w6
        m=m+  ne[counter,r_local,theta_local+1,phi_local+1]*w7
        m=m+  ne[counter,r_local+1,theta_local+1,phi_local+1]*w8

    return (m)

def Interpolate_Serial(grid_lat,grid_lon,grid_lev,daed_lat,daed_lon,daed_alt,zg,model_dt,orbit_dt,ne):
    counter=1
    deltaphi= np.abs(grid_lon[2]-grid_lon[1])
    deltatheta=np.abs(grid_lat[2]-grid_lat[1])
    m=np.zeros((len(daed_alt)))
    for i in range(0,len(daed_alt)):



        phi_local=local(len(grid_lon),daed_lon[i],grid_lon)
        theta_local=local(len(grid_lat),daed_lat[i],grid_lat)

        if daed_lon[i] >= 177.5 :
            phi_local=len(grid_lon)-2

        if daed_lon[i] <= -177.5:
            phi_local=0

        if daed_lat[i] >= 87.5:
            theta_local=len(grid_lat)-2

        if daed_lat[i] <= -87.5:
            theta_local=0


        alts=zg[counter,:,theta_local,phi_local]/100000
        r_local=local(len(alts),daed_alt[i],alts)

        deltarho=alts[r_local+1]-alts[r_local]

        dx=np.abs(((daed_alt[i]-alts[r_local])/deltarho))
        dy=np.abs(((daed_lat[i]-grid_lat[theta_local])/deltatheta))
        dz=np.abs(((daed_lon[i]-grid_lon[phi_local])/deltaphi))

        w1=(1-dx)*(1-dy)*(1-dz)
        w2=(dx)*(1-dy)*(1-dz)
        w3=(1-dx)*(dy)*(1-dz)
        w4=(dx)*(dy)*(1-dz)
        w5=(1-dx)*(1-dy)*(dz)
        w6=(dx)*(1-dy)*(dz)
        w7=(1-dx)*(dy)*(dz)
        w8=(dx)*(dy)*(dz)

        m[i]=0.0
        m[i]=       ne[counter,r_local,theta_local,phi_local]*w1
        m[i]=m[i]+  ne[counter,r_local+1,theta_local,phi_local]*w2
        m[i]=m[i]+  ne[counter,r_local,theta_local+1,phi_local]*w3
        m[i]=m[i]+  ne[counter,r_local+1,theta_local+1,phi_local]*w4
        m[i]=m[i]+  ne[counter,r_local,theta_local,phi_local+1]*w5
        m[i]=m[i]+  ne[counter,r_local+1,theta_local,phi_local+1]*w6
        m[i]=m[i]+  ne[counter,r_local,theta_local+1,phi_local+1]*w7
        m[i]=m[i]+  ne[counter,r_local+1,theta_local+1,phi_local+1]*w8

    return (m)


def DaedalusInterpolator(model,model_data_file,orbit_file,save,VAR,Parallel,F90):
# directory-->string:: output directory
# model_data_file,model_data_secondary_file--> strings:: model data files in data folder netCDF file provided from modellers
# orbit_file-->string :: orbit filename in Time-Lat-Lon-Alt format as porvided in Jupyter
# max_alt-->float:: max altitude for the interpolation--depends on model and user preference
# save--> Logical:: if true saves interpolated values to directory
# VAR--> string:: variable to inteprolate, must be  one of the variables included in the model data
#Outputs:: Plots+ 1 csv file





    #***************Read Model Data and Orbit Files*****************************
    # Interpolation Min and Max Altitude Along Track Depending on Model and User Preference
    min_alt=140
    max_alt=450

#     TIEGCM=Dataset(model_data_file)
#     df = pd.read_csv(orbit_file)
    TIEGCM=Dataset("../../../NAS/TIEGCM_DATA/"+model_data_file)
    df = pd.read_csv("../../../NAS/Data_Files/OrbitData/"+orbit_file)

    # Get data from Orbit
    daed_lat_temp = df["Lat_GEOD(deg)"]
    daed_lon_temp = df["Lon_GEOD(deg)"]
    daed_alt_temp = df["Height_WGS84 (km)"]
    daed_time_temp = df["Epoch(UTCG)"]

    # Get data from Model
    grid_lat=TIEGCM.variables['lat'][:]
    grid_lon=TIEGCM.variables['lon'][:]
    grid_lev=TIEGCM.variables['ilev'][:]
    grid_time=TIEGCM.variables['time'][:]
    zg=TIEGCM.variables['ZG'][:]

    # Select Variable to Interpolate
    ne=TIEGCM.variables[VAR][:]

    #Find model's temporal resolution
    model_dt=(grid_time[1]-grid_time[0])*60   #in seconds
    orbit_dt=1/16                             #in seconds

    # Sample Orbit and Find Positions Below 450km
    counter=0
    for i in range(0,len(daed_alt_temp)):
        if (daed_alt_temp[i] < max_alt and daed_alt_temp[i] > min_alt):
            counter=counter+1
    # **************************************************************************



    # ***************Allocate Arrays********************************************
    daed_lat=np.zeros((counter))
    daed_lon=np.zeros((counter))
    daed_alt=np.zeros((counter))
    daed_time=np.zeros((counter),dtype=datetime)
    index=[None]*counter
    int_final=[None]*len(daed_alt_temp)

    #Keep Data in MIN MAX Range For the Interpolation
    counter=0
    for i in range(0,len(daed_alt_temp)):
        if (daed_alt_temp[i] < max_alt and daed_alt_temp[i] > min_alt):
            daed_time[counter]=daed_time_temp[i]
            daed_lat[counter]=daed_lat_temp[i]
            daed_lon[counter]=daed_lon_temp[i]
            daed_alt[counter]=daed_alt_temp[i]
            index[counter]=i          #keep indices for merging data
            counter=counter+1

    #Transfrom Lon to Match TIEGCM (Orbit Longitudes range from 0 to 360 deg)
    for i in range (0, len(daed_alt)):
            daed_lon[i]=daed_lon[i]-180
    # **************************************************************************



    # **************************************************************************
    #Allocate a matrix for the interpolated values
    int_data=np.zeros(((counter)))


    print("Interpolation Starting...",counter, "Positions to Interpolate")
    #======================Call Interpolation Subroutines=======================
    t1=time.time()
    # ********************FORTRAN Subroutines***********************************
    if F90 == True:
        if Parallel==True:
            interpolation_OpenMP.daed_interp(grid_lat,grid_lon,grid_lev,daed_lat,daed_lon,daed_alt,zg,model_dt,orbit_dt,ne,int_data)
        else:
            interpolation.daed_interp(grid_lat,grid_lon,grid_lev,daed_lat,daed_lon,daed_alt,zg,model_dt,orbit_dt,ne,int_data)
    else:
    # **************************************************************************


    # ********************Python Subroutines***********************************
        if Parallel==True:
            if __name__=="__main__":
                pool=multiprocessing.Pool()
                func=partial(Interpolate_Parallel,grid_lat,grid_lon,grid_lev,daed_lat,daed_lon,daed_alt,zg,model_dt,orbit_dt,ne)
                i=range(len(daed_alt))
                int_data=pool.map(func,i)
                pool.close()
                pool.join()
        else:
            int_data=Interpolate_Serial(grid_lat,grid_lon,grid_lev,daed_lat,daed_lon,daed_alt,zg,model_dt,orbit_dt,ne)
    # **************************************************************************


    t2=time.time()
    # ==========================================================================
    print("Interpolation Finished in ",str(t2-t1),"s")


    # Merge Interpolated Values to Full Orbit Array (For printing full orbit)
    for i in range(0,len(int_data)):
        int_final[index[i]]=int_data[i]
    # **************************************************************************


     # ********************Export Data to DaedalusNAS****************************
    if save==True:
        directory = "../../../NAS/Data_Files/ModelsOutput/Interpolation/"

        Exports={'Time (UTCG)':daed_time,'Lat (deg)':daed_lat,'Lon (deg)':daed_lon,
                        'Alt (km)':daed_alt,'Interpolated_Data':int_data }

        df = DataFrame(Exports, columns= ['Time (UTCG)', 'Lat (deg)','Lon (deg)','Alt (km)','Interpolated_Data'])

        export_name=directory +os.path.splitext(os.path.basename(orbit_file))[0]+"_"+model+ "_"+ VAR+".csv"

        export_csv = df.to_csv (export_name, index = None, header=True)
    # **************************************************************************


    # **************************Plotting****************************************
    xaxis=np.arange(0,len(int_final)) #time axis based on Daedalus Sampling rate
    plt.figure(2)
    plt.scatter(xaxis/(1*60),int_final,c=daed_alt_temp,s=1)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.title("TIEGCM:" +VAR+" Along Daedalus' Orbit")
    plt.xlabel("Time [Minutes]")
    plt.ylabel("Eletron Density [$cm^{-3}$]")
    clb = plt.colorbar()
    plt.xlim(0,xaxis[-1]/(1*60))
    clb.set_label('Altitude [km]')
    plt.grid("True")

    plt.show()
    # **************************************************************************
    return(export_name)



# ******************************************************************************
# # # Example Call
# model_data_file="../Data/1.nc"
# orbit_file="../Data/2.csv"
# save=False
# VAR="NE"
# Parallel=True
# F90=True
# Parallel=False
# DaedalusInterpolator(model_data_file,orbit_file,save,VAR,Parallel,F90)
