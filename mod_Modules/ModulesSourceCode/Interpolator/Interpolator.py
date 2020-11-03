# Daedalus_Interpolator_V6.0
# ******************************************************************************
# This program reads Daedalus' orbit (LLA) and TIEGCM's data and then interpolates
# TIEGCM's variables (scalars) to Daedalus' orbit. In charge of the interpolation are
# a set of FORTRAN subroutines wrapped under f2py3. There are 2 different versions of
# those subroutines, Serial_Version and Parallel_Version. The Parallel_Version is
# using the OpenMP interface for multithreading.Interpolation Method: Trilinear 3D
# Instructions for FORTRAN subroutines compilation
# Update: Python I_lat),daed_lnterpolation Routines both serial and parallel
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
import SourceCode.DaedalusGlobals as DaedalusGlobals
# import interpolation
# import interpolation_OpenMP
import numpy as np
import pandas as pd
from pandas import DataFrame
from datetime import datetime
import time
import os
import multiprocessing
from functools import partial
from shutil import copyfile
import os


def local(dim1,y,x):
    local_pos=0
    for i in range(0,len(x)-1):
        if y >= x[i] and y < x[i+1]:
            local_pos=i
            break
    return (local_pos)


def geod_lat2geo_lat(phi):
    # Taken from Panagiotis Pirnaris

    # calculate geocentric latitude from geodetic latitude
    # according to WGS 84
    a = 6378137  # meter semi major axis of earth
    f = 1 / 298.257  # flattening
    b = a - f * a  # semi minor axis
    e = ((a ** 2 - b ** 2) ** (1 / 2)) / a
    phi_rad = np.deg2rad(phi)
    geo_lat = np.arctan((1 - e ** 2) * np.tan(phi_rad))
    geo_lat = np.rad2deg(geo_lat)
    return geo_lat  # in degrees

def Interpolate_Parallel(grid_lat,grid_lon,grid_lev,daed_lat,daed_lon,daed_alt,zg,model_dt,orbit_dt,ne,i):
    counter=1
    deltaphi= np.abs(grid_lon[2]-grid_lon[1])
    deltatheta=np.abs(grid_lat[2]-grid_lat[1])

    if i >-1:


        phi_local=local(len(grid_lon),daed_lon[i],grid_lon)
        theta_local=local(len(grid_lat),daed_lat[i],grid_lat)

        if daed_lon[i] >= 177.5 :
            phi_local=len(grid_lon)-1

        if daed_lon[i] <= -177.5:
            phi_local=0

        if daed_lat[i] >= 87.5:
            theta_local=len(grid_lat)-1

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
    Re=6378137.0/1000
    arc_theta=2*np.pi*Re*(deltatheta/360)
    arc_phi=2*np.pi*Re*(deltaphi/360)

    m=np.zeros((len(daed_alt)))

    #orbit block number of orbit timesteps to cover one model timestep
    # orbit at 1 Hz model at 20 min so we need 60*20 seconds to jump
    # a model timestep
    block=60 
    timecvrd=0


    for i in range(0,len(daed_alt)):

        if i%block  == 0 and counter<22:
            counter=counter +1
            timecvrd=0


        phi_local=local(len(grid_lon),daed_lon[i],grid_lon)
        theta_local=local(len(grid_lat),daed_lat[i],grid_lat)

        if daed_lon[i] >= 177.5 :
            phi_local=len(grid_lon)-1

        if daed_lon[i] <= -177.5:
            phi_local=1

        if daed_lat[i] >= 87.5:
            theta_local=len(grid_lat)-1

        if daed_lat[i] <= -87.5:
            theta_local=1

        if i==0:
            alts=zg[counter,:,theta_local,phi_local]/100000
        


        r_local=local(len(alts),daed_alt[i],alts)

        deltarho=alts[r_local+1]-alts[r_local]

        dx=(((daed_alt[i]-alts[r_local])/deltarho))
        dy=(((daed_lat[i]-grid_lat[theta_local])/deltatheta))
        dz=(((daed_lon[i]-grid_lon[phi_local])/deltaphi))



        w1=np.abs((1-dx)*(1-dy)*(1-dz))
        w2=np.abs((dx)*(1-dy)*(1-dz))
        w3=np.abs((1-dx)*(dy)*(1-dz))
        w4=np.abs((dx)*(dy)*(1-dz))
        w5=np.abs((1-dx)*(1-dy)*(dz))
        w6=np.abs((dx)*(1-dy)*(dz))
        w7=np.abs((1-dx)*(dy)*(dz))
        w8=np.abs((dx)*(dy)*(dz))

        m[i]=0.0
        m[i]=       ne[counter,r_local,theta_local,phi_local]*w1
        m[i]=m[i]+  ne[counter,r_local+1,theta_local,phi_local]*w2
        m[i]=m[i]+  ne[counter,r_local,theta_local+1,phi_local]*w3
        m[i]=m[i]+  ne[counter,r_local+1,theta_local+1,phi_local]*w4
        m[i]=m[i]+  ne[counter,r_local,theta_local,phi_local+1]*w5
        m[i]=m[i]+  ne[counter,r_local+1,theta_local,phi_local+1]*w6
        m[i]=m[i]+  ne[counter,r_local,theta_local+1,phi_local+1]*w7
        m[i]=m[i]+  ne[counter,r_local+1,theta_local+1,phi_local+1]*w8




        #Temporal Interpolation
        r_local=local(len(alts),daed_alt[i],alts)

        deltarho=alts[r_local+1]-alts[r_local]

        dx=(((daed_alt[i]-alts[r_local])/deltarho))
        dy=(((daed_lat[i]-grid_lat[theta_local])/deltatheta))
        dz=(((daed_lon[i]-grid_lon[phi_local])/deltaphi))
        
        w1=np.abs((1-dx)*(1-dy)*(1-dz))
        w2=np.abs((dx)*(1-dy)*(1-dz))
        w3=np.abs((1-dx)*(dy)*(1-dz))
        w4=np.abs((dx)*(dy)*(1-dz))
        w5=np.abs((1-dx)*(1-dy)*(dz))
        w6=np.abs((dx)*(1-dy)*(dz))
        w7=np.abs((1-dx)*(dy)*(dz))
        w8=np.abs((dx)*(dy)*(dz))
        
        m2=0.0
        m2=       ne[counter+1,r_local,theta_local,phi_local]*w1
        m2=m2+  ne[counter+1,r_local+1,theta_local,phi_local]*w2
        m2=m2+  ne[counter+1,r_local,theta_local+1,phi_local]*w3
        m2=m2+  ne[counter+1,r_local+1,theta_local+1,phi_local]*w4
        m2=m2+  ne[counter+1,r_local,theta_local,phi_local+1]*w5
        m2=m2+  ne[counter+1,r_local+1,theta_local,phi_local+1]*w6
        m2=m2+  ne[counter+1,r_local,theta_local+1,phi_local+1]*w7
        m2=m2+  ne[counter+1,r_local+1,theta_local+1,phi_local+1]*w8

        timecvrd=timecvrd+1

        t_w1=(block-timecvrd)/block
        t_w2=(timecvrd)/block

        m[i]=m[i]*t_w1+m2*t_w2

    return (m)





def Interpolator(model_data_file, orbit_file):
# model data file--> netcdf to read model data from
# orbit file--> netcdf to read orbit from 
# save--> Logical:: if true saves interpolated values to directory
# VAR--> string array:: variables to inteprolate, must be  one of the variables included in the model data
#Outputs:: Plots+ 1 csv file

    if orbit_file=="":
        return "" # <<<<

    save=True
    VAR=["NE","DEN","HE","NO","O1","O2","OP","TN","UI_ExB","VI_ExB","WI_ExB","UN","VN"]
    Parallel=False
    F90=False

    #***************Read Model Data and Orbit Files*****************************
    # Interpolation Min and Max Altitude Along Track Depending on Model and User Preference
    min_alt=110      #min altitude to interpolate
    max_alt=450      #max altitude to interpolate      
    srt=1            # daedalus sampling rate

    # Copy the orbit file from the OrbitData folder to ModelsOutput folder
    # TODO: if file exists just open it for append and check the Calculated attribute
    orbit_path = orbit_file[ 0 : orbit_file.rfind("/")+1]
    orbit_name = orbit_file[ orbit_file.rfind("/")+1 : -3]
    result_filename = DaedalusGlobals.Interpolated_Files_Path + orbit_name + "_Interpolated" + ".nc"
    if os.path.isfile(result_filename)==False  or  "Surface" in result_filename:
        copyfile(orbit_file, result_filename)
    
    # Get data from Orbit
    resultCDF=Dataset(result_filename, "a")
    
    if resultCDF.Calculated == "yes":
        print ("Skipping calculation because values are valid for" , result_filename)
        return result_filename # <<<<
    
    daed_lat_temp = resultCDF.variables['lat'][:]
    daed_lon_temp = resultCDF.variables['lon'][:]
    daed_alt_temp = resultCDF.variables['altitude'][:]
    daed_time_temp = resultCDF.variables['time'][:]

    # Get data from Model
    TIEGCM=Dataset(model_data_file)
    grid_lat=TIEGCM.variables['lat'][:]
    grid_lon=TIEGCM.variables['lon'][:]
    grid_lev=TIEGCM.variables['ilev'][:]
    grid_time=TIEGCM.variables['time'][:]
    zg=TIEGCM.variables['ZG'][:]


    #Find model's temporal resolution
    model_dt=(grid_time[1]-grid_time[0])*60   #in seconds
    orbit_dt=1/srt                            #in seconds

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
            daed_lat[i]=  geod_lat2geo_lat(daed_lat[i])


    # **************************************************************************



    # **************************************************************************
    #Allocate a matrix for the interpolated values
    int_data=np.zeros(((counter)))
    # **************************************************************************


    #======================Call Interpolation Subroutines=======================
    print("Interpolation Starting...",counter*len(VAR), "Positions to Interpolate")
    t1=time.time()
    exports=0
    for jj in range(0,len(VAR)):
        # Select Variable to Interpolate
        ne=TIEGCM.variables[VAR[jj]][:]

        # ********************FORTRAN Subroutines***********************************
        if F90 == True:
            if Parallel==True:
                print(VAR[jj])
                interpolation_OpenMP.daed_interp(grid_lat,grid_lon,grid_lev,daed_lat,daed_lon,daed_alt,zg,model_dt,orbit_dt,ne,int_data)
            else:
                print(VAR[jj])
                interpolation.daed_interp(grid_lat,grid_lon,grid_lev,daed_lat,daed_lon,daed_alt,zg,model_dt,orbit_dt,ne,int_data)
        else:
        # ********************Python Subroutines***********************************
            if Parallel==True:
                print(VAR[jj])
                if __name__=="__main__":
                    pool=multiprocessing.Pool(2)
                    func=partial(Interpolate_Parallel,grid_lat,grid_lon,grid_lev,daed_lat,daed_lon,daed_alt,zg,model_dt,orbit_dt,ne)
                    i=range(len(daed_alt))
                    int_data=pool.map(func,i)
                    pool.close()
                    pool.join()
            else:
                print(VAR[jj])
                int_data=Interpolate_Serial(grid_lat,grid_lon,grid_lev,daed_lat,daed_lon,daed_alt,zg,model_dt,orbit_dt,ne)
        # **************************************************************************
    #===============================================================================

        # Merge Interpolated Values to Full Orbit Array (For printing full orbit)
        for i in range(0,len(int_data)):
            int_final[index[i]]=int_data[i]
        # **************************************************************************

        if save==True:
            resultCDF.variables[VAR[jj]][:]= int_final
            # **********************************************************************§

    t2=time.time()
    resultCDF.Calculated = "yes"
    resultCDF.EditTime   = str(datetime.now())
    resultCDF.close()    
    
    print("Interpolation Finished in ",str(t2-t1),"s")
    
    return result_filename



# ******************************************************************************
# #Example Call
# model_data_file="../Data/weimer.nc"
# orbit_file="../Data/DAED_ORB_Evt0_LLA_Per120_Lat00_Srt01Hz_Msc.nc"
# save=True
# VAR=["NE","DEN","HE","NO","O1","O2","OP","TN","UI_ExB","VI_ExB","WI_ExB","UN","VN"]
# Parallel=False
# F90=False
# DaedalusInterpolator(model_data_file,orbit_file,save,VAR,Parallel,F90)
