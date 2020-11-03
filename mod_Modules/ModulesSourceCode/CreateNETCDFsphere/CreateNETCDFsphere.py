

from netCDF4 import Dataset  # https://unidata.github.io/netcdf4-python/netCDF4/index.html
from datetime import datetime
import numpy as np
import shutil 
import SourceCode.DaedalusGlobals as DaedalusGlobals


'''
Creates a NetCDF file which describes positions at a spherical surface around the Earth.
The file has the same structure as the Daedalus NetCDF orbit files. 
The values for Latitude  range between 87.5 and -87.5 with step <LatitudeStep>
The values for Longitude range between -180 and   175 with step <LongitudeStep>
Time and Altitude are fixed values for the whole surface.
A tiny offset value is added to all Lat,Lon,Time values so that Panoply can plot the data as trajectory.
RETURNS: the filename which has been created
'''
def CreateNETCDFsphere( fixedDatetimeString, fixedAltitude, LatitudeStep, LongitudeStep ):
    
    # if no date is provided then return an empty string
    if fixedDatetimeString == "":
        return "" # <<<<
    
    # convert the Time-string to Unix Timestamp
    aTimeString = fixedDatetimeString[:-6] + " " + "+0000" # indicate explicity that this is UTC time
    TimeObj = datetime.strptime( aTimeString, "%b %d %Y %H:%M:%S.%f %z" ) 
    FixedTimestamp = float( datetime.timestamp(TimeObj) ) * 1000  # timestamp in milliseconds

    # copy the standard surface file to a temporary location, where we will process it
    shutil.copyfile( DaedalusGlobals.Temporary_Files_Path+"SurfaceTemplate.nc", DaedalusGlobals.Temporary_Files_Path+"Surface.nc" )

    # open file
    CDFroot =  Dataset( DaedalusGlobals.Temporary_Files_Path+"Surface.nc", "a")

    # update general attributes
    CDFroot.Type = "surface"
    CDFroot.Calculated = "no" # reset state
    CDFroot.ResetTime = str(datetime.now()) # log

    # fill with pre-fabricated data for all Latitude-Lonitude pairs according to their respective steps
    n = 0;
    for i in np.arange(87.5, -92.5, -5.0):
        for j in np.arange(-180.0,  180.0, 5.0):
            CDFroot.variables["lat"][n] = i + 0.00000001*n # add a tiny offset, so that Panoply can plot the data as trajectory
            CDFroot.variables["lon"][n] = j + 0.00000001*n # add a tiny offset, so that Panoply can plot the data as trajectory
            #print( round(float(CDFroot.variables["lat"][n]),2), "    ",  round(float(CDFroot.variables["lon"][n]),2) )
            CDFroot.variables["time"][n]     = FixedTimestamp + n # add a tiny offset, so that Panoply can plot the data as trajectory
            CDFroot.variables["altitude"][n] = fixedAltitude
            CDFroot.variables["ALFA"][n]     = np.random.random_sample() # assign random values for test purposes
            n = n + 1

    # close file
    CDFroot.close() 
    
    
       
    return  DaedalusGlobals.Temporary_Files_Path+"Surface.nc"
