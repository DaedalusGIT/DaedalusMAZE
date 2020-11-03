
# This script calculates the motional Electric field VxB due to spacecraft motion in orbit.
# Inputs are velocity of spacecraft in the orbit and magnetic field using IGRF model

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas import DataFrame


# Plotting Function
# x--> variable to be plotted 1D array
# plt_serial--> integer :: plot number
# label_y--> string:: y label
# title--> string :: title of plot
def plot_along_orbit(x,plt_serial,label_y,title):
    plt.figure(plt_serial)
    plt.plot(x)
    plt.title(title)
    plt.xlabel("Measurements")
    plt.ylabel(label_y)
    plt.grid(False)
    plt.show()

# ******************************************************************************************
# Main
# file_IGRF -->String:: file to read Magnetic Field and Daedalus Velocities(INPUT)
# outfile--> string::filename to save plasma products(OUTPUT)

def main(file_IGRF,outfile,Save_Products,var_to_plot,Plot_Products,label_y,title):

#    Orbit_Data=pd.read_csv("../../../NAS/Data_Files/ModelOutputs/IGRF/"+ file_IGRF)
    Orbit_Data=pd.read_csv("./data/"+ file_IGRF)
    daed_time=Orbit_Data["Epoch(UTCG)"]
    daed_Bx=Orbit_Data["Bx_Tesla"]
    daed_By=Orbit_Data["By_Tesla"]
    daed_Bz=Orbit_Data["Bz_Tesla"]
    
    # Read Velocity Along Orbit
#    Ram_Data = pd.read_csv("../../../NAS/Data_Files/OrbitData/"+file_Velocities)
    Ram_Data = pd.read_csv("./data/"+file_Velocities)
    V_ram = Ram_Data["VMag(km/s)"]*1e3      # in m/s
    daed_vx = Ram_Data["RamX_GSE"]*V_ram    # in m/s
    daed_vy = Ram_Data["RamY_GSE"]*V_ram    # in m/s
    daed_vz = Ram_Data["RamZ_GSE"]*V_ram    # in m/s
    X_GSE = Ram_Data["X_GSE(km)"]     # in km 
    Y_GSE = Ram_Data["Y_GSE(km)"]     # in km
    Z_GSE = Ram_Data["Z_GSE(km)"]     # in km

  
    # ******************************************************************************************
    # Calculating Motional Electric Field Using VxB 
    # ******************************************************************************************
    Ex = -(daed_vy*daed_Bz - daed_vz*daed_By)  # in V/m  
    Ey = -(daed_vz*daed_Bx - daed_vx*daed_Bz)  # in V/m
    Ez = -(daed_vx*daed_By - daed_vy*daed_Bx)  # in V/m

    Etot = np.sqrt(Ex*Ex + Ey*Ey + Ez*Ez)         # in V/m

#     plt.plot(Ex)
#     plt.ylabel("Ex $(V/m)$")
#     plt.title('Motional Ex')
#     plt.savefig('MEx_DAED_ORB_Evt0_LLA_Per120_Lat80_Srt01Hz.png') 
    
#     plt.plot(Ey)
#     plt.ylabel("Ey $(V/m)$")
#     plt.title('Motional Ey')
#     plt.savefig('MEy_DAED_ORB_Evt0_LLA_Per120_Lat80_Srt01Hz.png') 
    

#     plt.plot(Ez)
#     plt.ylabel("Ez $(V/m)$")
#     plt.title('Motional Ez')
#     plt.savefig('MEz_DAED_ORB_Evt0_LLA_Per120_Lat80_Srt01Hz.png') 
    

#     plt.plot(daed_vx/1e3, color='b', label='vx')
#     plt.plot(daed_vy/1e3, color='g', label='vy')
#     plt.plot(daed_vz/1e3, color='r', label='vz')
#     plt.legend(loc=0)
#     plt.ylabel('V ram in GSE (km/s)')
#     plt.savefig('VramGSE.png')    

    plt.plot(Ex, label='Ex')
    plt.plot(Ey, label='Ey')    
    plt.plot(Ez, label='Ez')  
    plt.ylabel("E $(V/m)$")    
    plt.title('Motional E in GSE') # close to GSE not quite right yet
    plt.legend(loc=0)
    plt.savefig('ME_DAED_ORB_Evt0_LLA_Per120_Lat80_Srt01Hz.png') 

#     # ******************************************************************************************
#   Transforming E from GSE to spacecraft frame
#   Along is opposite of ram direction
#   Nadir/Earth is toward Earth center or r in spherical coordinates
#   Across is the cross product of Along x Nadir

    Ealong = -(Ex*daed_vx + Ey*daed_vy + Ez*daed_vz)/V_ram   # Component of E along spacecraft -E.Vram/|Vram| and it must be zero
    plt.plot(Ealong, label='Ealong')   
#     if var_to_plot== "":
#         plot_along_orbit(Ex  ,1  ,"Ex $(V/m)$","Motional Ex")
#         plot_along_orbit(Ey  ,2  ,"Ey $(V/m)$","Motional Ey")
#         plot_along_orbit(Ez  ,3  ,"Ez $(V/m)$","Motional Ez")
#         plot_along_orbit(Etot  ,4  ,"Et $(V/m)$","Motional E total")
#     else:

#         plot_along_orbit(Orbit_Data[var_to_plot],1 ,label_y,title)
#     # ******************************************************************************************

#     if Save_Products==True:
#         # Export Products for PIC initialization to CSV
#         Exports={'Time (UTCG)':daed_time,'Ex (V/m)':Ex,'Ey (V/m)':Ey,
#                     'Ez (V/m)':Ez,'Etot (V/m)':Etot}



#         df = DataFrame(Exports, columns= ['Time (UTCG)', 'Ex (V/m)','Ey (V/m)','Ez (V/m)','Etot (V/m)'])
#         export_csv = df.to_csv ("data/"+ outfile + '.csv', index = None, header=True)
#     # ******************************************************************************************

    return

## Example Call
from motionalE import *
file_IGRF = "DAED_ORB_Evt0_LLA_Per120_Lat80_Srt01Hz_Msc_IGRF_all.csv"
file_Velocities = "DAED_ORB_Evt0_PTG_Per120_Lat80_Srt01Hz_Msc.csv"
outfile="MotionalE"
Save_Products=False
Plot_Products=False
# Plot_All=False
var_to_plot=""
label_y="mV/m"
title="Motional E Along Orbit in GSE"
main(file_IGRF,outfile,Save_Products,var_to_plot,Plot_Products,label_y,title)