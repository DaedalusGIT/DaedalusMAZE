#!/usr/bin/bash
import numpy as np 
from netCDF4 import Dataset
import sys
from tqdm import tqdm





class varArray:
    def __init__(self,name=None):
        self.name=name
        self.data=[]

    def append(self,list1):
        self.data.extend(list1)
    
    def makeArray(self):
        a=np.array(self.data)
        return a


# Get files
files=[]
for i in range(1,len(sys.argv)):
    print(sys.argv[i])
    files.append(sys.argv[i])

# Get variable names
varNames=[]
nc=Dataset(files[0])
varNames+=nc.variables.keys()
nc.close()

# Init varArrays
variables={}
for i,name in enumerate(varNames):
    variables[name]=varArray(name)



# Merge objects
for i,file in enumerate(files):
    nc=Dataset(file)
    for j,var in enumerate(tqdm(varNames)):
        # print(var)
        list1=nc.variables[var][:]
        variables[var].append(list1)
    nc.close()



# Save to 1 file
ncout= Dataset("Merged.nc",'w',format='NETCDF4')
ncout.createDimension("time",len(variables["time"].data))

for name in variables:
    print(name)
    data1 = ncout.createVariable(name,'f4',"time")
    data1[:] = variables[name].makeArray()
   
ncout.close()   
