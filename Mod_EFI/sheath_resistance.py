# This script calculates the sheath resistance in orbit.
# Inputs are plasma parameters along the orbit
# Spacecraft potential needs to be calculated for I = 0




import matplotlib.pyplot as plt
import numpy as np
from scipy.io import FortranFile
from mpl_toolkits.mplot3d import Axes3D
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
    plt.grid(True)
    plt.show()



# ******************************************************************************************
# Main
# file_IRI -->String:: file to read IRI Data (INPUT)
# file_Velocities -->String:: file to read Daedalus Velocities(INPUT)
# outfile--> string::filename to save plasma products(OUTPUT)
# Save_Products-->Logical"" save products to outfile
# var_to_plot--> String:: specifies which variale to plot. if empty ("") plots all available variables
# Plot_Products--> Logical:: plots all products calculated
# label_y,title-->String:: if var_to_plot is given then the user specifies the y axis label and title of plot
#                     i.e. for electron density : title="Electron Density", label_y="NE $cm^-3$"
# OTHER OUTPUTS:: plots
def main(file_IRI,file_Velocities,outfile,Save_Products,var_to_plot,Plot_Products,label_y,title):

    Model_Data = pd.read_csv("../../../DataFiles/OrbitData/"+file_IRI)
    daed_time = Model_Data["Epoch(UTCG)"]
    daed_lat = Model_Data["Lat_GEOD(deg)"]
    daed_lon = Model_Data["Lon_GEOD(deg)"]
    daed_alt = Model_Data["Height_WGS84 (km)"]
    n_e = Model_Data["ne_iri16_cm-3"]*1e6    # in 1/m3
    n_Op = Model_Data["O+_iri16_cm-3"]*1e6   # in 1/m3
    n_i = n_e-n_Op
    #n_O2p = Model_Data["O2+_iri16_cm-3"]*1e6 # in 1/m3
    #n_NOp = Model_Data["NO+_iri16_cm-3"]*1e6  # in 1/m3
    T_e = Model_Data["Te_iri16_K"]       # in K
    T_i = Model_Data["Ti_iri16_K"]       # in K
    T_Op = T_i
    #T_O2p = T_i
    #T_NOp = T_i

    #=============
    # Constants
    #=============

    #=============
    # Boom Dimensions
    #=============
    length = 5.0     # in m
    radius = 0.06     # 6 cm

    A_full = 4. * np.pi * radius * radius   # Fully Sphere
    A_ram  = np.pi * radius * radius        # Cross-sectional area

    #Lambda_D = np.sqrt( e0 * Kb * T_e / (n_e * qe * qe))


    e0 = 8.85e-12    # permitivity of free space
    qe = 1.6e-19     # fundamental charge in Coloumb
    Kb = 1.38e-23    # Boltzman constant
    me = 9.11e-31    # electron mass in kg
    mi = 1.67e-27    # ion mass in kg
    mOp = 16.0 * mi   # Oxygen mass in kg
    #mO2p = 16.0 * mi   # Oxygen mass in kg
    #mNOp = 16.0 * mi   # Oxygen and Nitrogen mass in kg
    RE = 6371.0         # in km
 
    Ram_Data = pd.read_csv("../../../DataFiles/OrbitData/"+file_Velocities)
    V_ram = Ram_Data["VMag(km/s)"]*1e3 # in m/s
    RamX = Ram_Data["RamX_GSE"]*V_ram   # in m/s
    RamY = Ram_Data["RamY_GSE"]*V_ram    # in m/s
    RamZ = Ram_Data["RamZ_GSE"]*V_ram    # in m/s
    X_GSE = Ram_Data["X_GSE(km)"]     # in km 
    Y_GSE = Ram_Data["Y_GSE(km)"]     # in km
    Z_GSE = Ram_Data["Z_GSE(km)"]     # in km

    npoints = daed_time.shape[0]
    RS = np.zeros(npoints)
    V_I0 = np.zeros(npoints)
    I_V0 = np.zeros(npoints)
    Gain = np.zeros(npoints)
    FRC = np.zeros(npoints)
    GDC = np.zeros(npoints)
    GAC = np.zeros(npoints)


    RB = 100e6 # Base Resistance in Ohms
    CB = 100e-12 # Capacitance in Farads
    CS = 4*np.pi*e0*radius
    #omega = 1e2 # in Hz
    for k in range(1616,5416):
        print('k =', k)
         # ******************************************************************************************

        #================
        # For different values of spacecraft potential 
        #===============

        Vsc = np.linspace(-2.0,2.0,200)  #between -2 and 2 Volts

        # ******************************************************************************************
        #=============
        # Currents
        #=============

        #=============
        # Photoelectron Current
        #=============

        alpha = 0.05
        V1 = 2.7   # volts empirical photoeletron energy
        V2 = 10.0  #volts empirical photoeletron energy
        jphoto0 = 20.0e-6 # 20 to 60 microA/m^2 at 1AU

        Jphoto = np.zeros(200)

        Jphoto[np.where(Vsc >=0)] = jphoto0*( (1.0-alpha)*np.exp(-Vsc[np.where(Vsc >=0)]/V1) + alpha*np.exp(-Vsc[np.where(Vsc >=0)]/V2) )  # Vsc >= 0
        Jphoto[np.where(Vsc <0)] = jphoto0  # Vsc < 0 
        A_photo = A_ram
        I_photo     = Jphoto*A_photo

       # Finding shadows in orbit X_GSE < 0 AND np.sqrt(Y_GSE*Y_GSE + Z_GSE*Z_GSE) > 1 RE -> I_photo = 0 in shadow
        # This goes into the for loop along the orbit.

        if ( (X_GSE[k] < 0) and (np.sqrt(Y_GSE[k]*Y_GSE[k] + Z_GSE[k]*Z_GSE[k])> RE) ):
            I_photo = 0.0


        #=============
        # Thermal currents -- thermal velocity sqrt(kb Tj/mj), j = i, e, o
        #=============

        Jthe0     = np.sqrt( Kb * T_e[k] / me ) * n_e[k] * qe   
        Jthi0     = np.sqrt( Kb * T_i[k] / mi ) * n_i[k] * qe
        Jtho0     = np.sqrt( Kb * T_Op[k] / mOp ) * n_Op[k] * qe

        #================
        # Thermal currents for positive spacecraft potential Vsc >=0
        #================
        I_e_thermal = np.zeros(200)
        I_i_thermal = np.zeros(200)
        I_o_thermal = np.zeros(200)


        I_e_thermal[np.where(Vsc >=0)]  = A_full * Jthe0 * np.sqrt(1.0 + qe*Vsc[np.where(Vsc >=0)] / (Kb * T_e[k]) )
        I_i_thermal[np.where(Vsc >=0)]  = A_full * Jthi0 * np.exp( -qe * Vsc[np.where(Vsc >=0)] / (Kb * T_i[k]) )
        I_o_thermal[np.where(Vsc >=0)]  = A_full * Jtho0 * np.exp( -qe * Vsc[np.where(Vsc >=0)] / (Kb * T_Op[k]) )


        #================
        # Thermal currents for negative spacecraft potential Vsc < 0
        #================

        I_e_thermal[np.where(Vsc <0)] = A_full * Jthe0 * np.exp( qe * Vsc[np.where(Vsc <0)] / (Kb * T_e[k]) )
        I_i_thermal[np.where(Vsc <0)] = A_full * Jthi0 * np.sqrt(1. - qe*Vsc[np.where(Vsc <0)] / (Kb * T_i[k]) )
        I_o_thermal[np.where(Vsc <0)] = A_full * Jtho0 * np.sqrt(1. - qe*Vsc[np.where(Vsc <0)] / (Kb * T_Op[k]) )




        #=============
        # Ram Current V_ram is total ram velocity
        #=============

        I_e_ram     = A_ram * V_ram[k] * qe * n_e[k]
        I_i_ram     = A_ram * V_ram[k] * qe * n_i[k]
        I_o_ram     = A_ram * V_ram[k] * qe * n_Op[k]

        #================
        # Total current
        #================

        I_total     = -(I_photo - I_e_thermal + I_i_thermal + I_o_thermal + I_i_ram + I_o_ram - I_e_ram)


        R = 1./np.gradient(I_total, Vsc) 

        plt.plot(Vsc, I_total)
        plt.plot(Vsc, R)

        IL0_ind = np.asarray(I_total <=0).nonzero() # make it better
        I0_ind = IL0_ind[0][-1]

        RS[k] = R[I0_ind]
        V_I0[k] = Vsc[I0_ind]
        I_V0[k] = I_total[100]
        
        nn = 4000
        omega = np.linspace(0.01,4e6,nn)
        ZB = RB / (1. + 1j*omega*CB*RB)
        ZS = RS[k] / ( 1.+ 1j*omega*CS*RS[k])
        G = np.zeros(nn)
        G = ZB /(ZB+ZS)
        plt.semilogx(omega,G.real)
        plt.ylabel("Gain")
        plt.xlabel("Frequency")
        #plt.savefig('Gain.png')    
        gind = np.where((G.real < 0.55) & (G.real > 0.48))
        FRC[k] = omega[gind[0][0]]
        GDC[k] = G.real[2]
        GAC[k] = G.real[3900]
    # end of loop

    plt.plot(RS/1e6)
    plt.ylabel("Resistance (Mega Ohm)")

    plt.plot(V_I0)
    plt.plot(I_V0)

    plt.semilogy(FRC/(2*np.pi*1e3))
    plt.ylabel("R-C cross over Frequency (kHz)")
    plt.savefig('FRC.png') 
    
    plt.plot(np.sqrt(X_GSE*X_GSE+Y_GSE*Y_GSE+Z_GSE*Z_GSE)/RE)
    plt.plot(X_GSE/RE)
    
    
#     # ******************************************************************************************
#     # Plotting Model Data
#     if var_to_plot== "":
#         plot_along_orbit(NE  ,1  ,"NE $cm^-3$","Electron Density NE")
#         plot_along_orbit(Te  ,2  ,"Te $K$","Electron Temperature K")
#         plot_along_orbit(Ti  ,3  ,"Ti $K$","Ion Temperature K")
#         plot_along_orbit(Op  ,4  ,"OP $cm^-3$","Ion Oxygen Density OP")
#         plot_along_orbit(O2p ,5  ,"O2P $cm^-3$","Ion Oxygen 2 Density O2P")
#         plot_along_orbit(NO  ,6  ,"NO $cm^-3$","Nitric Oxide Density NO")
#     else:

#         plot_along_orbit(Model_Data[var_to_plot],1 ,label_y,title)

#     if Plot_Products==True:
#         plot_along_orbit(t_ratio,7 ,"Ti/Te","T Ratio")
#         plot_along_orbit(L_D,8 ,"Ld m","Debye Length Along Orbit")
#         plot_along_orbit(mach_number,9 ,"Mach Number","Mach Number (km/sec)")
#         plot_along_orbit(wp,10 ,"$V_{p}$" ,"Plasma Frequency")

#     # ******************************************************************************************



#     # ******************************************************************************************

    if Save_Products==True:
        # Export Products for PIC initialization to CSV
        Exports={'Time (UTCG)':daed_time,'Lat (deg)':daed_lat,'Lon (deg)':daed_lon,'Alt (km)':daed_alt,
                 'RSheath (Ohm)':RS,'FRC (Hz)':FRC,'GDC':GDC,'GAC':GAC,'V_I0':V_I0}



        df = DataFrame(Exports, columns= ['Time (UTCG)', 'Lat (deg)','Lon (deg)','Alt (km)','RSheath (Ohm)','FRC (Hz)',
                                          'GDC','GAC','V_I0 (volt)'])
        export_csv = df.to_csv (outfile + '.csv', index = None, header=True)
    # ******************************************************************************************




    return

## Example Call
#file_IRI="DAED_ORB_Evt0_LLA_Per120_Lat00_Srt01Hz_Msc_IRI16_all.csv"
#file_Velocities="DAED_ORB_Evt0_PTG_Per120_Lat00_Srt01Hz_Msc.csv"

file_IRI='DAED_ORB_Evt0_LLA_Per120_Lat80_Srt01Hz_Msc_IRI16_all.csv'
file_Velocities="DAED_ORB_Evt0_PTG_Per120_Lat80_Srt01Hz_Msc.csv"

outfile="Gain_Products"
Save_Products=True
Plot_Products=True
# Plot_All=False
var_to_plot="NO+_iri16_cm-3"
label_y="NO cm^-3"
title="IRI: Nitric Oxide Density Along Orbit"
main(file_IRI,file_Velocities,outfile,Save_Products,var_to_plot,Plot_Products,label_y,title)


