
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas import DataFrame

file_Velocities = "DAED_ORB_Evt0_PTG_Per120_Lat80_Srt01Hz_Msc.csv"

Ram_Data = pd.read_csv("./data/"+file_Velocities)
V_ram = Ram_Data["VMag(km/s)"]*1e3      # in m/s
vx = Ram_Data["RamX_GSE"]*V_ram    # in m/s
vy = Ram_Data["RamY_GSE"]*V_ram    # in m/s
vz = Ram_Data["RamZ_GSE"]*V_ram    # in m/s

RE = 6371.0
X = Ram_Data["X_GSE(km)"] #/RE    # in km 
Y = Ram_Data["Y_GSE(km)"] #/RE     # in km
Z = Ram_Data["Z_GSE(km)"] #/RE     # in km


#Apply Inertia tensor analysis
Ixx = np.sum( Z**2 + Y**2)
Ixy = -np.sum( X * Y)
Ixz = -np.sum( X * Z)

Iyx = Ixy
Iyy = np.sum( Z**2 + X**2)
Iyz = -np.sum( Y * Z)

Izx = Ixz
Izy = Iyz
Izz = np.sum( Y**2 + X**2)

tensor = np.zeros((3,3))

tensor[0,0] = Ixx
tensor[1,0] = Ixy
tensor[2,0] = Ixz

tensor[0,1] = Iyx
tensor[1,1] = Iyy
tensor[2,1] = Iyz

tensor[0,2] = Izx
tensor[1,2] = Izy
tensor[2,2] = Izz

#Normalize Interia Tensor
inertia = tensor / np.sum(tensor)

eig_val,eig_vec = np.linalg.eig(inertia)

e1=np.array([-0.88260019,0.15159485,-0.44501226]) #perp to orbital plane
e2=np.array([0.40270657,-0.24463839,-0.88203145])
e3=np.array([0.24257851,0.95769048,-0.15486965])


pos = np.array([X,Y,Z])
for k in range(0,7080):
    perp[k]=np.dot(pos[:,k],e1)
    
