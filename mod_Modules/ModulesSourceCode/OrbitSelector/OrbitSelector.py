"""
Constructs and returns orbit full-path filenames of csv and NetCDF formats from the orbit properties.
The filename includes all parameters required to select the appropriate orbit.
All filenames  have the same length and same number of parameters. 
In case OrbitFilename and EvtXY are empty then an empty string is returned.
FILE FORMAT: DAED_ORB_EVTXY_TYP_PERYYY_LATZZ_SRTQQHz_XSC.csv
  Parameter OrbitFilename:
    If it is not empty then the rest arguments are ignored. 
    If it contains slashes then is is assumed it is a full path name and it is returned as it is.
    If it does not contain slashes then the orbits path is added.
    
  Parameter EvtXY values:
    EVTXS    X Event, Single Orbit
    EVTXA    X Event, All Orbit
    EVT1Y    1st Event: St Patrick’s day event [17 Mar 2015 – 20 Mar 2015]
    EVT2Y    2nd Event ...  ...
	
  Parameter TYP	values:
    LLA    Time,Latitude Longitude Altitude 
    VEL    Time,VmagVxVyVz
    PTG    Time, X_GSE, Y_GSE, Z_GSE, RamX_GSE, RamY_GSE, RamZ_GSE

  Parameter PerYYY values:
    PER120, PER150    Perigee Altitude at 120km or 150km

  Parameter LatZZ values:
    LAT00, LAT40, LAT80    Perigee Latitude at 0°, 40° or 80°
	
  Parameter SRXXHZ values:
    SRT16Hz, SRT01Hz    Sampling rate 16Hz or 1HZ
	
  Parameter SC values:
    MSC    Mother Spacecraft
    SSC    Sub-Spacecraft
"""
CSV_FILES_PATH = "/home/NAS/Data_Files/OrbitData/" # absolute path where the CSV filenames with the orbit data reside.
NETCDF_FILES_PATH = "/home/NAS/Data_Files/OrbitData/" # absolute path where the NetCDF filenames with the orbit data reside.

def OrbitSelector( OrbitFilename, EvtXY, TYP, PerYYY, LatZZ, SRXXHZ, SC ):
    filename_csv = ""
    filename_cdf = ""
	
    if len(OrbitFilename)==0  and  len(EvtXY)==0 :
        filename_csv = ""
        filename_cdf = ""
    elif len(OrbitFilename) > 0:
        if "/" in OrbitFilename: 
            filename_csv = OrbitFilename
            filename_cdf = OrbitFilename
        else:
            filename_csv = CSV_FILES_PATH + OrbitFilename
            filename_cdf = NETCDF_FILES_PATH + OrbitFilename
    else:
        filename_csv = CSV_FILES_PATH    + "DAED_ORB_" + EvtXY + "_" + TYP + "_" + PerYYY + "_" + LatZZ + "_" + SRXXHZ + "_" + SC + ".csv"
        filename_cdf = NETCDF_FILES_PATH + "DAED_ORB_" + EvtXY + "_" + TYP + "_" + PerYYY + "_" + LatZZ + "_" + SRXXHZ + "_" + SC + ".nc"
	
    return filename_csv, filename_cdf


