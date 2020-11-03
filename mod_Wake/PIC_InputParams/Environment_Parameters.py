# ===============Plasma Parameters nd Products================================!
# THIS PROGRAMS READS THE MODEL OUTPUTS FOR A GIVEN DAEDALUS' ORBIT AND
# PLOTS OUTPUT VARIABLES. IT ALSO CALCULATES PRODUCTS FROM THE PLASMA PARAMETERS
# TO BE USED IN PARTICLE IN CELL SIMULATIONS FOR DAEDALUS.
# KEP 2019
# ============================================================================!
#
#                    ____                 _       _
#                   |  _ \  __ _  ___  __| | __ _| |_   _ ___
#                   | | | |/ _` |/ _ \/ _` |/ _` | | | | / __|
#                   | |_| | (_| |  __/ (_| | (_| | | |_| \__ \
#                   |____/ \__,_|\___|\__,_|\__,_|_|\__,_|___/




import matplotlib.pyplot as plt
import numpy as np
from scipy.io import FortranFile
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
from pandas import DataFrame
import os



def plot_along_orbit(x,plt_serial,label_y,title,daed_alt,lim1,lim2,dt):
# Plotting Function
# x--> variable to be plotted 1D array
# plt_serial--> integer :: plot number
# label_y--> string:: y label
# title--> string :: title of plot

    plt.figure(plt_serial)
    xaxis=np.arange(0,len(daed_alt))
    plt.plot(xaxis/dt,x)      #line plot
    # plt.scatter(xaxis/dt,x,c=daed_alt,s=0.8)
    # clb = plt.colorbar()
    # clb.set_label('Altitude [km]')
    plt.title(title)
    plt.xlabel("Time [Minutes]")
    plt.ylabel(label_y)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.grid(True)
    plt.xlim(lim1,lim2)
    # plt.show()



# ******************************************************************************************
# Main
# file_IRI -->String:: file to read IRI Data (INPUT)
# file_Velocities -->String:: file to read Daedalus Velocities(INPUT)
# outfile--> string::filename to save plasma products(OUTPUT)
# Save_Plasma_Products-->Logical"" save plasma products to outfile
# var_to_plot--> String:: specifies which variale to plot. if empty ("") it plots all available variables
# Plot_Products--> Logical:: plots plasma all plasma products calculated
# label_y,title-->String:: if var_to_plot is given then the user specifies the y axis label and title of plot
#                     i.e. for electron density : title="Electron Density", label_y="NE $cm^-3$"
# OTHER OUTPUTS:: plots
def main(file_IRI,file_MSISE00,file_IGRF,file_Velocities,outfile,Save_Plasma_Products,var_to_plot,Plot_Products,label_y,title):

    # **********************************************************
    # Read Files
    IRI_Data=pd.read_csv("data/"+file_IRI)
    Orbit_Data=pd.read_csv("data/"+file_Velocities)
    MSISE00_Data=pd.read_csv("data/"+file_MSISE00)
    IGRF_Data=pd.read_csv("data/"+file_IGRF)

    #Daedalus Stats
    daed_time=IRI_Data["Epoch(UTCG)"]
    daed_lat=IRI_Data["Lat_GEOD(deg)"]
    daed_lon=IRI_Data["Lon_GEOD(deg)"]
    daed_alt=IRI_Data["Height_WGS84 (km)"]

    vmag=Orbit_Data["VMag(km/s)"]
    daed_vx=Orbit_Data["RamX_GSE"]*vmag*1000
    daed_vy=Orbit_Data["RamY_GSE"]*vmag*1000
    daed_vz=Orbit_Data["RamZ_GSE"]*vmag*1000

    #Get IRI Outputs
    NE=IRI_Data["ne_iri16_cm-3"]
    Te=IRI_Data["Te_iri16_K"]
    Ti=IRI_Data["Ti_iri16_K"]
    Op=IRI_Data["O+_iri16_cm-3"]
    O2p=IRI_Data["O2+_iri16_cm-3"]
    NO=IRI_Data["NO+_iri16_cm-3"]


    #Get MSISE00 Outputs
    Tn=MSISE00_Data["Tn_msise00_K"]

    #Get IGRF Outputs
    Bt=IGRF_Data["B_Tesla"]
    Bx=IGRF_Data["Bx_Tesla"]
    By=IGRF_Data["By_Tesla"]
    Bz=IGRF_Data["Bz_Tesla"]
    Dst=IGRF_Data["DST_index_nT"]

    # *********************************************************



    # *********************************************
    #OLD FILE FORMAT
    #Get Orbit Data
    # daed_time=IRI_Data["Time (UTCG)"]
    # daed_lat=IRI_Data["Lat (deg)"]
    # daed_lon=IRI_Data["Lon (deg)"]
    # daed_alt=IRI_Data["Alt (km)"]

    #Get Velocity Along Orbit
    # daed_vx=Orbit_Data["vx (km/sec)"]
    # daed_vy=Orbit_Data["vy (km/sec)"]
    # daed_vz=Orbit_Data["vz (km/sec)"]
    # *********************************************

    # ******************************************************************************************


    # ******************************************************************************************

    #Sampling Rate
    srt=1 #Hz
    dt=60*srt

    #Find Orbit Position of MIN DST
    min=Dst[0]
    for i in range(0,len(Dst)):
        if Dst[i] < min :
            min=Dst[i]
            pos=i

    #Find apogee
    apogee=daed_alt[0]
    for i in range(0,len(daed_alt)):
        if daed_alt[i]>apogee:
            apogee=daed_alt[i]

    #Find Closest previous apogee to Min DST
    for i in range(pos,1,-1):
        if (daed_alt[i]<daed_alt[i-1] and daed_alt[i-2]<daed_alt[i-1]):
            apg1=i
            break

    #Find Closest next apogee to Min DST
    for i in range(pos,len(daed_alt)):
        if (daed_alt[i]<daed_alt[i+1] and daed_alt[i+2]<daed_alt[i+1]):
            apg2=i
            break


    # Convert apg1 and apg2 to minutes flying
    apg1_t=np.floor(apg1/dt)
    apg2_t=np.ceil(apg2/dt)


    # Plasma Parameters Calculations

    # Constants
    sound_speed=334       #m/s
    cc=299792458          #speed of light in [m/s]
    e0=8.85418e-12        #permitivitty of free space
    q_e=1.60217e-19       #electron charge [C]
    k_b=1.3806e-23        #boltzman constant
    m_e=9.10938356e-31    #electron mass  [kg]
    m_i=1.67e-27
    Op_mass=16*m_i
    O2p_mass=32*m_i
    NO_mass=30*m_i


    # Convert Densities to m^-3
    NE=NE*1e6
    Op=Op*1e6
    O2p=O2p*1e6
    NO=NO*1e6


    # Ti/Te t_ratio
    t_ratio=Ti/Te

    #Debye Length
    L_D=np.sqrt((e0*k_b*Te)/(NE*q_e*q_e))



    # Total Ion Density
    Ions_total=Op+O2p+NO

    # Weights fot Ion Density
    w1=Op/Ions_total
    w2=O2p/Ions_total
    w3=NO/Ions_total

    # Average Ion Mass from weights
    avg_mi=(w1*Op_mass+w2*O2p_mass+w3*NO_mass)
    # M=w1*Op_mass+w2*O2p_mass+w3*NO_mass
    gamma=1.4       #adiabatic constant
    Z=(Op+O2p+NO)/NE    #plasma charge state

    # Sound Speed
#     cs=np.sqrt(gamma*Z*k_b*Te/avg_mi)
    cs=np.sqrt(k_b*Te/m_i)
    # cs=np.sqrt((Te+0.7*Ti)/avg_mi )

    # Mach number
    sat_speed=np.sqrt(daed_vx*daed_vx+daed_vy*daed_vy+daed_vz*daed_vz)  #m/s
    mach_number=sat_speed/cs

    # Plasma Frequency rad/s
    wp=np.sqrt(((NE)*q_e*q_e)/(e0*m_e))  #rad/s
    # wp=np.sqrt((4*np.pi*NE*1e-6/q_e*q_e)) #rad/s

    # Plasma Frequency Hz
    fp=(1/2*np.pi)*np.sqrt(((NE)*q_e*q_e)/(e0*m_e))  #Hz


    #Plasma skin Depth
    skin_depth=(cc/(np.sqrt(2)*fp))

    # ******************************************************************************************



    # ******************************************************************************************
    # Plotting Model Data
    if var_to_plot== "":
        plot_along_orbit(NE  ,1  ,"NE [$m^-3$]","Electron Density NE",daed_alt,apg1_t,apg2_t,dt)
        plot_along_orbit(Te  ,2  ,"Te [$K$]","Electron Temperature ",daed_alt,apg1_t,apg2_t,dt)
        plot_along_orbit(Ti  ,3  ,"Ti [$K$]","Ion Temperature ",daed_alt,apg1_t,apg2_t,dt)
        plot_along_orbit(Tn  ,4  ,"Tn [$K$]","Neutral Temperature ",daed_alt,apg1_t,apg2_t,dt)
        plot_along_orbit(Op  ,5  ,"OP [$m^-3$]","Ion Oxygen Density OP",daed_alt,apg1_t,apg2_t,dt)
        plot_along_orbit(O2p ,6  ,"O2P [$m^-3$]","Ion Oxygen 2 Density O2P",daed_alt,apg1_t,apg2_t,dt)
        plot_along_orbit(NO  ,7  ,"NO [$m^-3$]","Nitric Oxide Density NO",daed_alt,apg1_t,apg2_t,dt)
        plot_along_orbit(Bt  ,8  ,"Bt [T]","Total B Field",daed_alt,apg1_t,apg2_t,dt)
        plot_along_orbit(Dst  ,9  ,"DST  [nT]","DST Index",daed_alt,apg1_t,apg2_t,dt)


    else:

        plot_along_orbit(IRI_Data[var_to_plot],1 ,label_y,title,daed_alt,apg1_t,apg2_t,apg1_t,apg2_t,dt)

    if Plot_Products==True:
        plot_along_orbit(t_ratio,10 ,"Ti/Te","T Ratio Along Orbit",daed_alt,apg1_t,apg2_t,dt)
        plot_along_orbit(L_D,11 ,"$L_{d}$ [m]","Debye Length Along Orbit",daed_alt,apg1_t,apg2_t,dt)
        plot_along_orbit(mach_number,12 ,"M","Mach Number Along Orbit",daed_alt,apg1_t,apg2_t,dt)
        plot_along_orbit(fp,13 ,"$f_{p}$ [Hz]","Plasma Frequency Along Orbit",daed_alt,apg1_t,apg2_t,dt)
        plot_along_orbit(skin_depth  ,14  ,"$c_{i}$ [m]","Plasma Skin Depth - Electron Inertial Length",daed_alt,apg1_t,apg2_t,dt)
        plot_along_orbit(cs ,15,"$C_{s}$ [m/s]","Speed of Sound Along Orbit",daed_alt,apg1_t,apg2_t,dt)
        plot_along_orbit(daed_alt ,16,"Height WGS84 [km]","Altitude",daed_alt,apg1_t,apg2_t,dt)
        plot_along_orbit(daed_lat ,17,"Lat GEOD [deg]","Latitude",daed_alt,apg1_t,apg2_t,dt)
        # plot_along_orbit(cs ,16,"$C_{s}$ [m/s]","Speed of Sound Along Orbit",daed_alt)
        # plot_along_orbit(cs ,14,"$C_{s}$ [m/s]","Speed of Sound Along Orbit",daed_alt)


    # ******************************************************************************************



    # ******************************************************************************************
    # if not os.path.exists("output"):
    #     os.makedirs("output")


    if Save_Plasma_Products==True:
        # Export Plasma Products for PIC initialization to CSV
        Exports={"Epoch(UTCG)":daed_time,"Lat_GEOD(deg)":daed_lat,"Lon_GEOD(deg)":daed_lon,
                    "Height_WGS84 (km)":daed_alt,"RamX_GSE":daed_vx,"RamY_GSE":daed_vy,"RamZ_GSE":daed_vz,
                    'Ti/Te':t_ratio,'Debye Length (m)':L_D,'Plasma Frequency (rad/sec)':wp,
                    'Mach Number':mach_number }



        df = DataFrame(Exports, columns= ["Epoch(UTCG)", "Lat_GEOD(deg)","Lon_GEOD(deg)","Height_WGS84 (km)","RamX_GSE",
                        "RamY_GSE","RamZ_GSE",'Ti/Te','Debye Length (m)','Plasma Frequency (rad/sec)','Mach Number'])
        export_csv = df.to_csv ("data/"+outfile +"_PLASMA_PRODUCTS"+ '.csv', index = None, header=True)
    # ******************************************************************************************
    plt.show()
    return







# #****************************** Example Call************************************
# #Model Output Filename
# # *****************************EVT0*********************************************
# # file_IRI="Evt0/DAED_ORB_Evt0_LLA_Per120_Lat00_Srt01Hz_Msc_IRI16_all.csv"
# # file_MSISE00="Evt0/DAED_ORB_Evt0_LLA_Per120_Lat00_Srt01Hz_Msc_MSISE00_all.csv"
# # file_IGRF="Evt0/DAED_ORB_Evt0_LLA_Per120_Lat00_Srt01Hz_Msc_IGRF_all.csv"
# # # File with Daedalus Velocity
# # file_Velocities="Evt0/DAED_ORB_Evt0_PTG_Per120_Lat00_Srt01Hz_Msc.csv"
# # # Filename of Outputs for this code
# # outfile="DAED_ORB_Evt0_LLA_Per120_Lat00_Srt01Hz_Msc"
# # ******************************************************************************

# # *****************************EVT1*********************************************
# file_IRI="Evt1/DAED_ORB_Evt1_LLA_Per120_Lat80_Srt01Hz_Msc_IRI16_all.csv"
# file_MSISE00="Evt1/DAED_ORB_Evt1_LLA_Per120_Lat80_Srt01Hz_Msc_MSISE00_all.csv"
# file_IGRF="Evt1/DAED_ORB_Evt1_LLA_Per120_Lat80_Srt01Hz_Msc_IGRF_all2.csv"
# # File with Daedalus Velocity
# file_Velocities="Evt1/DAED_ORB_Evt1_PTG_Per120_Lat80_Srt01Hz_Msc.csv"
# # Filename of Outputs for this code
# outfile="DAED_ORB_Evt1_LLA_Per120_Lat80_Srt01Hz_Msc"
# # *******************************************************************************


# # Toggles
# Save_Plasma_Products=True
# Plot_Products=True
# # Plot_All=False

# # Specifiers Model Variable +label for plot + title for plot. If var to plot ="" then plots all.
# # var_to_plot="ne_iri16_cm-3"
# var_to_plot=""
# label_y="NE [$cm^{-3}$]"
# title="IRI: NE Density Along Orbit"

# # Call
# main(file_IRI,file_MSISE00,file_IGRF,file_Velocities,outfile,Save_Plasma_Products,var_to_plot,Plot_Products,label_y,title)
# # **********************************************************************
