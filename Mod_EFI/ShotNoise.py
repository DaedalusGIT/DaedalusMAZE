# This script calculates the shot noise in orbit.
# Inputs are T_e, n_e, Zs, Vsc, along the orbit



file_IRI='DAED_ORB_Evt0_LLA_Per120_Lat80_Srt01Hz_Msc_IRI16_all.csv'
file_Velocities="DAED_ORB_Evt0_PTG_Per120_Lat80_Srt01Hz_Msc.csv"
file_Gain="Gain_Products.csv"


import matplotlib.pyplot as plt
import numpy as np
from scipy.io import FortranFile
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
from pandas import DataFrame


#

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
def main(file_IRI,file_Velocities,file_Gain,outfile,Save_Products,var_to_plot,Plot_Products,label_y,title):

    #Model_Data = pd.read_csv("../../../DataFiles/OrbitData/"+file_IRI)
    Model_Data = pd.read_csv("./data/"+file_IRI)
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
    A0 = 1.0         # Spacecraft charging correction factor (A focusing factor), set to 1 to disable
 
    #Ram_Data = pd.read_csv("../../../DataFiles/OrbitData/"+file_Velocities)
    Ram_Data = pd.read_csv("./data/"+file_Velocities)
    V_ram = Ram_Data["VMag(km/s)"]*1e3 # in m/s
    RamX = Ram_Data["RamX_GSE"]*V_ram   # in m/s
    RamY = Ram_Data["RamY_GSE"]*V_ram    # in m/s
    RamZ = Ram_Data["RamZ_GSE"]*V_ram    # in m/s
    X_GSE = Ram_Data["X_GSE(km)"]     # in km 
    Y_GSE = Ram_Data["Y_GSE(km)"]     # in km
    Z_GSE = Ram_Data["Z_GSE(km)"]     # in km

    Gain_Data = pd.read_csv("./data/"+file_Gain)
    RS= Gain_Data["RSheath (Ohm)"] # in ohm
    V_I0 = Gain_Data["V_I0 (volt)"]   # in volt



    RB = 100e6 # Base Resistance in Ohms
    CB = 100e-12 # Capacitance in Farads
    CS = 4*np.pi*e0*radius

    #omega = 1e2 # in Hz
    for k in range(1616,5416):
        print('k =', k)
         # ******************************************************************************************
        
        nn = 50000
        omega = np.linspace(0.01,5e5,nn) #0.01 to 5e5 Hz frequency range
        #ZB = RB / (1. + 1j*omega*CB*RB)
        ZS = RS[k] / ( 1.+ 1j*omega*CS*RS[k])
        vth_e = np.sqrt(Kb * T_e[k]/ me)  # in m/s

        # Current density
        j_e = n_e[k] * vth_e / 4.   # why divide by 4?
        
        N_impact = j_e*A_full*A0
        
        ZS = RS[k] / ( 1.+ 1j*omega*CS*RS[k])
        ZS_sq = ZS*np.conj(ZS)
        
        # shot noise voltage due to electrons in Volt^2/Hz
        V_sn_sq = 2.* qe**2 * N_impact * ZS_sq * np.exp((qe*V_I0[k])/(Kb*T_e[k]))
        # check qtn for probes Meyer-Vernet 1983 paper    
        plt.loglog(omega,V_sn_sq)
        plt.ylabel("SN e- (in volt^2/Hz) ")
        plt.xlim(1,5e5)
        plt.xlabel('Frequency (Hz)')
        plt.title('SN e- (in volt^2/Hz) Rs=%08.2f' % (RS[k]))
        plt.savefig("./Figs/"+"shotnoise%04d.png" % (k))
        plt.clf()
        
        
        
        # end of loop

    plt.semilogy(FRC/(2*np.pi*1e3))
    plt.ylabel("R-C cross over Frequency (kHz)")
    plt.savefig('FRC.png') 
    
    plt.plot(np.sqrt(X_GSE*X_GSE+Y_GSE*Y_GSE+Z_GSE*Z_GSE)/RE)
    plt.plot(X_GSE/RE)
        

#     # ******************************************************************************************

    if Save_Products==True:
        # Export Products for PIC initialization to CSV
        Exports={'Time (UTCG)':daed_time,'Lat (deg)':daed_lat,'Lon (deg)':daed_lon,
                    'Alt (km)':daed_alt,'RSheath (Ohm)':RS,'FRC (Hz)':FRC,'GDC':GDC,'GAC':GAC}



        df = DataFrame(Exports, columns= ['Time (UTCG)', 'Lat (deg)','Lon (deg)','Alt (km)','RSheath (Ohm)','FRC (Hz)',
                                          'GDC','GAC'])
        export_csv = df.to_csv (outfile + '.csv', index = None, header=True)
    # ******************************************************************************************




    return


outfile="Gain_Products"
Save_Products=True
Plot_Products=True
# Plot_All=False
var_to_plot="NO+_iri16_cm-3"
label_y="NO cm^-3"
title="IRI: Nitric Oxide Density Along Orbit"
main(file_IRI,file_Velocities,outfile,Save_Products,var_to_plot,Plot_Products,label_y,title)



