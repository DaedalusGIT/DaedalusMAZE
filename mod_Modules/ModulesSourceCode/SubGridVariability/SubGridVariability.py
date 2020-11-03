# Sub_Grid_Variability_Demo_(V1.0)
# ********************************************************************************
#  Earth is mapped onto a rectrangular grid (NxM) with dimension N=144 and M=72
# (LAT-LON). The computational domain is one lattice taken from the full domain.
# Then the lattice is refined by adding extra gridpoints turning it into a new
# AxB mesh grid. The artificial noise is then analytically calculated at every
# grid point in the lattice. Surface and spectral plots of the noise are
# evaluated. Daedalus is then sampling the noise passing through the lattice.
# The noise can be modified by means of wavenumber and amplitude in each dimension.
# KEP 5/6/19
# ********************************************************************************


            #  ____                 _       _
            # |  _ \  __ _  ___  __| | __ _| |_   _ ___
            # | | | |/ _` |/ _ \/ _` |/ _` | | | | / __|
            # | |_| | (_| |  __/ (_| | (_| | | |_| \__ \
            # |____/ \__,_|\___|\__,_|\__,_|_|\__,_|___/


import numpy as np
import matplotlib.colors
import math
import matplotlib.pyplot as plt
import csv
import pandas as pd
from pandas import DataFrame
import pickle
import random
import numpy.fft as fft
import scipy.fftpack as syfp
import cProfile




def noise(w_x,w_y,A_x,A_y,x,y):
    #Inputs
    # A_x,A_y--->Amplitude of noise Kind:: real
    # w_x,w_y--->Radial Frequency of noise Kind:: real
    #Outputs
    # Kind:: real , Target main
    return A_x*np.sin((w_x)*x*(np.pi/180))+A_y*np.sin(w_y*y*(np.pi/180))

def SubGridVariability( NoiseWx, NoiseWy, NoiseAx, NoiseAy, filename, ValueName ):
#Inputs
#NoiseAx,NoiseAy: Amplitude of noise Kind:: real
# NoiseWx,NoiseWy: Wavenumber of noise Kind:: real
# Filename: Character, field to read orbit from
# ValueName: variable from model data to use and add SGV. Corresponds to CSV Column : Character
#Outputs
# Plots..

    if filename is None:
        return None
    if len(filename) == 0:
        return ""

# ******************************************************************************************
    w_x=NoiseWx
    w_y=NoiseWy
    A_x=NoiseAx
    A_y=NoiseAy
# ******************************************************************************************


# ******************************************************************************************
    # Calculate number of gridpoints in each dimension considering LLA co-ordinates
    Lx=2*np.pi                              #Domain length in x dimension (lat)
    Ly=np.pi                                #Domain length in y dimension (long)
    deltaphi=2.5                            #Anglular spatial step
    Nx=int((Lx*(180/np.pi))/deltaphi)       #Number of gridpoints in x dimension
    Ny=int((Ly*(180/np.pi))/deltaphi )      #Number of gridpoints in y dimension
# ******************************************************************************************


# ******************************************************************************************

    # Input Daedalus' Orbit Data with Interpolated Values
    file=filename
    outfile=file[:-4]+"_"+ValueName+"_"+"SGV"        # build name of outfile from input file + variable

    Data=pd.read_csv(filename)
    
    try:
        daed_time=Data["Time (UTCG)"]
        daed_lat=Data["Lat (deg)"]
        daed_lon=Data["Lon (deg)"]
        daed_alt=Data["Alt (km)"]
    except:
        daed_time=Data["Epoch(UTCG)"]
        daed_lat=Data["Lat_GEOD(deg)"]
        daed_lon=Data["Lon_GEOD(deg)"]
        daed_alt=Data["Height_WGS84 (km)"]
  
    if len(ValueName) > 0:
        daed_int_data=Data[ValueName]
    else:
        daed_int_data=Data[ Data.columns[4] ]

    # Spatial Step in azimuth direction
    deltaphi=daed_lon[3]-daed_lon[2]
# ******************************************************************************************


# ******************************************************************************************
    # Fist Spectral Analysis of Interpolated Data without the added variability
    plt.figure(1)
    plt.xlabel("Frequency")
    plt.ylabel("Power")
    plt.grid(True)
    plt.title("Power Spectrum Analysis of Measured Data w/o Sub-Grid Variability")
    sig=daed_int_data
    sig=[0 if math.isnan(x) else x for x in sig]
    N=len(sig)
    sampled_points=N
    sig_fft = np.fft.fft(sig)/N*2
    fft_freq = np.fft.fftfreq(sampled_points,deltaphi/360)
    plt.plot(fft_freq[:N//2], 20*np.log10(abs((sig_fft[:N//2]))))

    # Plot of Measured Data
    plt.figure(3)
    plt.xlabel("Measurements")
    plt.ylabel("Electron Density $cm-3$")
    plt.grid(True)
    plt.title("IRI: Electron Density without Sub-Grid Variability")
    plt.plot(daed_int_data)
# ******************************************************************************************


# ******************************************************************************************
    # Add Sub Grid Variability

    # Allocate Array
    daed_int_data_sgv=np.zeros((len(daed_lon)))

    # Add noise to Interpolated Data to get Sub Grid Variability
    for i in range(0,len(daed_lon)):
        daed_int_data_sgv[i]=daed_int_data[i]+noise(w_x,w_y,A_x,A_y,daed_lat[i],daed_lon[i])
# ******************************************************************************************


# ******************************************************************************************
    # Second Spectral Analysis of Interpolated Data with the added variability
    plt.figure(2)
    plt.xlabel("Frequency")
    plt.ylabel("Power")
    plt.title("Power Spectrum Analysis of Measured Data with Sub-Grid Variability")
    N=len(daed_int_data_sgv)
    sampled_points=N
    sig=daed_int_data_sgv
    sig=[0 if math.isnan(x) else x for x in sig]
    sig_fft = np.fft.fft(sig)/N*2
    fft_freq =np.fft.fftfreq(sampled_points,deltaphi/360)
    plt.plot(fft_freq[:N//2], 20*np.log10(abs((sig_fft[:N//2]))))
    plt.grid(True)

    # Plot of Measured Data with Sub-Grid Variability
    plt.figure(4)
    plt.xlabel("Measurements")
    plt.ylabel("Electron Density (cm-3)")
    plt.grid(True)
    plt.title("IRI: Electron Density with Sub-Grid Variability")
    plt.plot(daed_int_data_sgv)
# ******************************************************************************************


# ******************************************************************************************
# Computational Domain Build for Visualization
    M=np.zeros((Nx,Ny))
    x=0
    # #Create computational domain from lattice and sample noise
    for i in range(0,Nx):
        y=0
        for j in range(0,Ny):
            x=x+deltaphi
            y=j+deltaphi
            M[i,j]=noise(w_x,w_y,A_x,A_y,x,y)


# Surface Plot of Added Variability
    M=np.transpose(M,(1,0))
    plt.figure(5)
    plt.title("Variability Visualiztion")
    plt.xlabel("Lon")
    plt.ylabel("Lat")
    plt.imshow(M,cmap="rainbow")
    plt.show()
# ******************************************************************************************


# ******************************************************************************************
# Export Sub Grid Variability to CSV
    Exports={'Time (UTCG)':daed_time,'Lat (deg)':daed_lat,'Lon (deg)':daed_lon,
                'Alt (km)':daed_alt,'Variable_SGV':daed_int_data_sgv }


    df = DataFrame(Exports, columns= ['Time (UTCG)', 'Lat (deg)','Lon (deg)','Alt (km)','Variable_SGV'])
  #  export_csv = df.to_csv (outfile + '.csv', index = None, header=True)
# ******************************************************************************************

    return outfile + '.csv'

# # Example Call and  Pre Defined Values
# NoiseWx=100                 #Wavenumber in x dimension
# NoiseWy=50                  #Wavenumber in y dimension
# NoiseAx=2000/2              #Amplitude of Variability in x dimension
# NoiseAy=2000/2              #Amplitude of Variability in x dimension
# ValueName="ne_iri16_cm-3"   #Value to operate on
# filename="IRI.csv"          #Filename to Read Model Data
# Daedalus_SGV(NoiseWx,NoiseWy,NoiseAx,NoiseAy,filename,ValueName)