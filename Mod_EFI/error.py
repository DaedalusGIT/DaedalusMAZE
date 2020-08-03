# error propagation in E field measurements
# calculate for each component and differnt boom length

# Esig = dV*Gain/Leff - EMotion
# Emotion=-VscxB independent from boom length
# Gain is independent from boom length
# dV depends on shot noise, preamplifier, analog an digital electronics
# calculate Gain, Emotion for one orbit as input



import matplotlib.pyplot as plt
import numpy as np
from scipy.io import FortranFile
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
from pandas import DataFrame
from scipy.interpolate import interp1d
from scipy import signal

file_IRI = 'DAED_ORB_Evt0_LLA_Per120_Lat80_Srt01Hz_Msc_IRI16_all.csv'
file_Velocities= "DAED_ORB_Evt0_PTG_Per120_Lat80_Srt01Hz_Msc.csv"
file_IGRF = "DAED_ORB_Evt0_LLA_Per120_Lat80_Srt01Hz_Msc_IGRF_all.csv"
file_Gain = "Gain_Products.csv"
file_bend = "Tipdeflections.xlsx"
# getting magnetic field data along the orbit
Orbit_Data=pd.read_csv("./data/"+ file_IGRF)
Bx=Orbit_Data["Bx_Tesla"]
By=Orbit_Data["By_Tesla"]
Bz=Orbit_Data["Bz_Tesla"]
 

# getting plasma parameters along the orbit
Model_Data = pd.read_csv("./data/"+file_IRI)
time = Model_Data["Epoch(UTCG)"]
lat = Model_Data["Lat_GEOD(deg)"]
lon = Model_Data["Lon_GEOD(deg)"]
alt = Model_Data["Height_WGS84 (km)"]
ne = Model_Data["ne_iri16_cm-3"]*1e6    # in 1/m3
nOp = Model_Data["O+_iri16_cm-3"]*1e6   # in 1/m3
ni = ne-nOp

Te = Model_Data["Te_iri16_K"]       # in K
Ti = Model_Data["Ti_iri16_K"]       # in K
TOp = Ti

# getting spacecraft velocity along the orbit
Ram_Data = pd.read_csv("./data/"+file_Velocities)
V_ram = Ram_Data["VMag(km/s)"]*1e3 # in m/s
Vx = Ram_Data["RamX_GSE"]*V_ram   # in m/s
Vy = Ram_Data["RamY_GSE"]*V_ram    # in m/s
Vz = Ram_Data["RamZ_GSE"]*V_ram    # in m/s
X_GSE = Ram_Data["X_GSE(km)"]     # in km 
Y_GSE = Ram_Data["Y_GSE(km)"]     # in km
Z_GSE = Ram_Data["Z_GSE(km)"]     # in km

# getting Gain parameters along the orbit
Gain_Data = pd.read_csv("./data/"+file_Gain)
RS = Gain_Data["RSheath (Ohm)"] # in Ohm
#FRC = Gain_Data["FRC (Hz)"] # in Hz
VI0 = Gain_Data["V_I0 (volt)"] # in Volts where Itot=0
GDC = Gain_Data["GDC"] # Gain in DC Frequency
GAC = Gain_Data["GAC"] # Gain in AC Frequency

Bend_data = pd.read_excel("./data/"+file_bend)
alt_bend = Bend_data["Altitude(km)"]
anglex= Bend_data["6m_angle_orthogonally_mounted"] #*0.0254 # in meter
angley= Bend_data["6m-Zstar_angle"] #*0.0254  # in meter
anglez= Bend_data["6m+Zstar_angle"] #*0.0254  # in meter
# anglex= Bend_data["4m_angle_orthogonally_mounted"] #*0.0254 # in meter
# angley= Bend_data["4m-Zstar_angle"] #*0.0254  # in meter
# anglez= Bend_data["4m+Zstar_angle"] #*0.0254  # in meter


anglex_fix = np.zeros(7080)
angley_fix = np.zeros(7080)
anglez_fix = np.zeros(7080)

anglex_fix = np.where(alt>552,anglex[84],0)
anglex_interpol = signal.resample_poly(anglex,9,1)
anglex_fix[3505:4251]=anglex_interpol[0:746]
anglex_fix[2763:3504]=anglex_interpol[741:0:-1]
anglex_fix[3490:3520] = anglex[0]

angley_fix = np.where(alt>552,angley[84],0)
angley_interpol = signal.resample_poly(angley,9,1)
angley_fix[3505:4251]=angley_interpol[0:746]
angley_fix[2763:3504]=angley_interpol[741:0:-1]
angley_fix[3490:3520] = angley[0]

anglez_fix = np.where(alt>552,anglez[84],0)
anglez_interpol = signal.resample_poly(anglez,9,1)
anglez_fix[3505:4251]=anglez_interpol[0:746]
anglez_fix[2763:3504]=anglez_interpol[741:0:-1]
anglez_fix[3490:3520] = anglez[0]



FRC = 2000. #10.

# constants
#length = 6.0     # in m (4,5,6,7.5,10)
radius = 0.03     # 3 cm

Afull = 4. * np.pi * radius * radius   # Fully Sphere
Aram  = np.pi * radius * radius        # Cross-sectional area

e0 = 8.85e-12    # permitivity of free space
qe = 1.6e-19     # fundamental charge in Coloumb
Kb = 1.38e-23    # Boltzman constant
me = 9.11e-31    # electron mass in kg
mi = 1.67e-27    # ion mass in kg
mOp = 16.0 * mi   # Oxygen mass in kg
RE = 6371.0         # in km


RB = 100e6 # Base Resistance in Ohms
CB = 100e-12 # Capacitance in Farads
CS = 4*np.pi*e0*radius
ZB = RB / (1. + 1j*FRC*CB*RB)
ZS = RS / ( 1.+ 1j*FRC*CS*RS)

# finding vhat, nhat and chat (Ram, Nadir, Cross track in spacecraft frame)

e1=np.array([-0.88260019,0.15159485,-0.44501226]) #perp to orbital plane constant in all orbit
e2=np.array([0.40270657,-0.24463839,-0.88203145])
e3=np.array([0.24257851,0.95769048,-0.15486965])


#vhat = -(Vx i^ + Vy j^ + Vz k^)/Vram along sc velocity or anti-ram
# nhat = chat X vhat
# chat = e1 from inertia tensor
vhat = np.zeros((3,7080))
nhat = np.zeros((3,7080))
chat = np.zeros((3,7080))

vhat = -1*np.vstack((Vx/V_ram,Vy/V_ram,Vz/V_ram))

chat[0,:]=e1[0]
chat[1,:]=e1[1]
chat[2,:]=e1[2]

nhatx = chat[1,:]*vhat[2,:] - chat[2,:]*vhat[1,:] 
nhaty = chat[2,:]*vhat[0,:] - chat[0,:]*vhat[2,:] 
nhatz = chat[0,:]*vhat[1,:] - chat[1,:]*vhat[0,:] 



Vv = Vx*vhat[0,:] + Vy*vhat[1,:] + Vz*vhat[2,:]
Vn = Vx*nhat[0,:] + Vy*nhat[1,:] + Vz*nhat[2,:]
Vc = Vx*chat[0,:] + Vy*chat[1,:] + Vz*chat[2,:]

Bv = Bx*vhat[0,:] + By*vhat[1,:] + Bz*vhat[2,:]
Bn = Bx*nhat[0,:] + By*nhat[1,:] + Bz*nhat[2,:]
Bc = Bx*chat[0,:] + By*chat[1,:] + Bz*chat[2,:]

# ******************************************************************************************
# Calculating Motional Electric Field Using VxB in SC frame
# ******************************************************************************************

# Ev = E.vhat = 0, En = E.nhat, Ec = E.chat

# Ev = Emx*vhat[0,:] + Emy*vhat[1,:] + Emz*vhat[2,:]
# En = Emx*nhat[0,:] + Emy*nhat[1,:] + Emz*nhat[2,:]
# Ec = Emx*chat[0,:] + Emy*chat[1,:] + Emz*chat[2,:]
Bv=Bx
Bn=By
Bc=Bz


Ev = -(Vn*Bc - Vc*Bn)  # in V/m  
En = -(Vc*Bv - Vv*Bc)  # in V/m
Ec = -(Vv*Bn - Vn*Bv)  # in V/m

# Rotate Velocity and B to boom frame
Vel_12 = Vv*np.cos(45*np.pi/180.)
Vel_34 = Vv*np.sin(45*np.pi/180.)*np.cos(45*np.pi/180.)
Vel_56 = Vv*np.sin(45*np.pi/180.)*np.cos(45*np.pi/180.)

B_12 = Bv*np.cos(45*np.pi/180.)-Bn*np.cos(45*np.pi/180.)
B_34 = Bv*np.sin(45*np.pi/180.)*np.cos(45*np.pi/180.)+Bn*np.cos(45*np.pi/180.)*np.cos(45*np.pi/180.)-Bc*np.cos(45*np.pi/180.)
B_56 = Bv*np.sin(45*np.pi/180.)*np.cos(45*np.pi/180.)+Bn*np.cos(45*np.pi/180.)*np.cos(45*np.pi/180.)+Bc*np.cos(45*np.pi/180.)


# error terms

# 0.1% error in spacecraft velocity
dVv = 0.001*Vel_12   # in m/s
dVn = 0.001*Vel_34   # in m/s
dVc = 0.001*Vel_56   # in m/s

# 5nT error in Magnetic field
dBv = 5 * 10 ** (-9) # in Tesla
dBn = 5 * 10 ** (-9) # in Tesla
dBc = 5 * 10 ** (-9) # in Tesla

# error in plasma parameters
dTi = 0.1*Ti
dTe = 0.1*Te
dTOp = 0.1*TOp
dni = 0.05*ni
dne = 0.05*ne
dnOp = 0.05*nOp


# Motional E error terms
dEmxerr = np.sqrt( (B_56*dVn)**2 + (Vel_34*dBc)**2 + (B_34*dVc)**2 + (Vel_56*dBn)**2 )
dEmyerr = np.sqrt( (B_12*dVc)**2 + (Vel_56*dBv)**2 + (B_56*dVv)**2 + (Vel_12*dBc)**2 )
dEmzerr = np.sqrt( (B_12*dVn)**2 + (Vel_34*dBv)**2 + (B_34*dVv)**2 + (Vel_12*dBn)**2 )
dEmerr = np.sqrt(dEmxerr*dEmxerr+dEmyerr*dEmyerr+dEmzerr*dEmzerr)


# Gain error Vsc<0 in this case
Ithi = Afull*np.sqrt( Kb * Ti / mi ) * ni * qe * np.sqrt(1. - qe*VI0 / (Kb * Ti) )
IthOp = Afull*np.sqrt( Kb * TOp / mOp ) * nOp * qe * np.sqrt(1. - qe*VI0 / (Kb * TOp) )
Ithe = Afull*np.sqrt( Kb * Te / me ) * ne * qe * np.exp( qe * VI0 / (Kb * Te) )

dIthidV = -qe * Ithi / (Kb*Ti - qe*VI0)
dIthOpdV = -qe * IthOp / (Kb*TOp - qe*VI0)
dIthedV = qe * Ithe / (Kb*Te)

dRsdn = (dIthidV*dni/ni)**2 + (dIthOpdV*dnOp/nOp)**2 + (dIthedV*dne/ne)**2


dRsdT = (Kb*dIthidV*dTi/np.sqrt(Kb*Ti-qe*VI0))**2 + (Kb*dIthOpdV*dTOp/np.sqrt(Kb*TOp-qe*VI0))**2 + ((-qe*VI0/(Kb*Te)-0.5)*dIthedV*dTe/Te)**2


dRs = RS*np.sqrt( dRsdn + dRsdT )


dGerr = (-(ZB*ZS**2)/((RS**2)*(ZB+ZS)**2)*dRs)**2
dGerr_sq = dGerr*np.conj(dGerr)

# error from electronics (analog, digital, preamplifier) based on Maven LF receiver
dVPDA = 9.3808315e-07  #in V/m

vth_e = np.sqrt(Kb * Te/ me)  # in m/s
j_e = ne * vth_e / 4.   # why divide by 4?
N_impact = j_e*Afull
ZS_sq = ZS*np.conj(ZS)

# shot noise voltage due to electrons in Volt/sqrt(Hz) Hz = 16
Vsn = np.sqrt(2.* qe**2 * N_impact * ZS_sq * np.exp((qe*VI0)/(Kb*Te)))*4



dVsn = np.sqrt( (0.5*Vsn*dne/ne)**2 + (ZS_sq*Vsn*dRs/RS**3)**2 + ((0.25+0.5*qe*VI0/(Kb*Te))*Vsn*dTe/Te)**2)


# calculate Eperp and Epar for each boom pair
# 
E12par = En*np.cos(45*np.pi/180.)
E12perp = np.sqrt(Ec**2 + (En*np.sin(45*np.pi/180.))**2)

E34par = En*(np.cos(45*np.pi/180.)**2)-Ec*np.cos(45*np.pi/180.)
E34perp = En*np.cos(45*np.pi/180.)*np.sin(45*np.pi/180.)+Ec*np.sin(45*np.pi/180.)

E56par = En*(np.cos(45*np.pi/180.)**2)+Ec*np.cos(45*np.pi/180.)
E56perp = En*np.cos(45*np.pi/180.)*np.sin(45*np.pi/180.)-Ec*np.sin(45*np.pi/180.)



# bending angle of the booms
#theta = 2.*np.pi/180.

# for 4m and 6m boom
length=6.
thetax= anglex_fix*np.pi/180. #np.arctan(lx_fix/length)*180./np.pi
thetazm= angley_fix*np.pi/180. #np.arctan(ly_fix/length)*180./np.pi
thetazp= anglez_fix*np.pi/180. #np.arctan(lz_fix/length)*180./np.pi

dL = 0.0

Leff = length*(np.cos(thetazm)+ np.cos(thetazp))

# V is motional field*Leff    #+ natural field*Leff  set to zero 

dV12 = E12par*Leff  #+50e-3*Lxeff
dV34 = E34par*Leff  #+50e-3*Lyeff
dV56 = E56par*Leff  #+50e-3*Lzeff

dV12_par_err = length*E12par*(2.-np.cos(thetazm)-np.cos(thetazp))
dV34_par_err = length*E34par*(2.-np.cos(thetazm)-np.cos(thetazp))
dV56_par_err = length*E56par*(2.-np.cos(thetazm)-np.cos(thetazp))
dV12_perp_err = length*E12perp*(np.sin(thetazm)-np.sin(thetazp))
dV34_perp_err = length*E34perp*(np.sin(thetazm)-np.sin(thetazp))
dV56_perp_err = length*E56perp*(np.sin(thetazm)-np.sin(thetazp))



# FIXME add dV_par_err  and dV_perp_err
dVerr12 = np.sqrt(dVPDA**2+dVsn**2+dV12_par_err**2+dV12_perp_err**2 )
dVerr34 = np.sqrt(dVPDA**2+dVsn**2+dV34_par_err**2+dV34_perp_err**2 )
dVerr56 = np.sqrt(dVPDA**2+dVsn**2+dV56_par_err**2+dV56_perp_err**2 )




# total error

#first = dV*GDC/Leff*np.sqrt((dVerr/VI0)**2 +np.sqrt(np.real(dGerr_sq/GDC))**2 + (0.5*dL/length)**2)
first12 = dV12*GDC/Leff*np.sqrt((dVerr12/dV12)**2 + np.sqrt(np.real(dGerr_sq/GDC))**2 + (0.5*dL/length)**2)
first34 = dV34*GDC/Leff*np.sqrt((dVerr34/dV34)**2 + np.sqrt(np.real(dGerr_sq/GDC))**2 + (0.5*dL/length)**2)
first56 = dV56*GDC/Leff*np.sqrt((dVerr56/dV56)**2 + np.sqrt(np.real(dGerr_sq/GDC))**2 + (0.5*dL/length)**2)

#dEerr = np.sqrt((first)**2 + dEmerr**2)
dE12err = np.sqrt((first12)**2 + dEmxerr**2)
dE34err = np.sqrt((first34)**2 + dEmyerr**2)
dE56err = np.sqrt((first56)**2 + dEmzerr**2)

fig, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel('orbit')
ax1.set_ylabel('Bending angle (degrees)', color=color)
ax1.plot(thetazp*180./np.pi, color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('Altitude (Km)', color=color)  # we already handled the x-label with ax1
ax2.plot(alt, color=color)
ax2.tick_params(axis='y', labelcolor=color)

plt.title('Boom length = 6m')
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.savefig('thetax6.png')



fig, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel('orbit')
ax1.set_ylabel('E12 Error (mV/m)', color=color)
ax1.plot(dE12err*1e3, color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('Altitude (Km)', color=color)  # we already handled the x-label with ax1
ax2.plot(alt, color=color)
ax2.tick_params(axis='y', labelcolor=color)

plt.title('Boom length = 6m')
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.savefig('XErrorinorbit6.png')

fig, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel('orbit')
ax1.set_ylabel('E34 Error (mV/m)', color=color)
ax1.plot(dE34err*1e3, color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('Altitude (Km)', color=color)  # we already handled the x-label with ax1
ax2.plot(alt, color=color)
ax2.tick_params(axis='y', labelcolor=color)

plt.title('Boom length = 6m')
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.savefig('YErrorinorbit6.png')


fig, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel('orbit')
ax1.set_ylabel('E56 Error (mV/m)', color=color)
ax1.plot(dE56err*1e3, color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('Altitude (Km)', color=color)  # we already handled the x-label with ax1
ax2.plot(alt, color=color)
ax2.tick_params(axis='y', labelcolor=color)

plt.title('Boom length = 6m')
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.savefig('ZErrorinorbit6.png')



