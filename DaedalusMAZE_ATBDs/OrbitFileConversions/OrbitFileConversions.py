import sys
sys.path.insert(1, '../CoverageCalculator/')
import Conversions as Conversions

import pandas as pd
import apexpy as ap
import numpy as np
from datetime import datetime
from dateutil.relativedelta import relativedelta
import csv


'''
reads a csv orbit file with format:
    Time (UTCG)
    Lat (deg)
    Lon (deg)
    Alt (km)
and creates a new csv orbit file with the same fields plus the extra fields:
    Daedalus.Magnetic Latitude
    Daedalus.Magnetic Longitude
    Daedalus.MLT
The values are calculated for modified apex at 90 km - used by TIEGCM as well
'''
def AddMagneticCoordinates_ModifiedApex90( sourceFilename, resultFilename ):
    CSVfileOUT = open(resultFilename, 'w')
    CSVwriter  = csv.writer(CSVfileOUT, delimiter=',')
    CSVwriter.writerow( ["Time (UTCG)", "Lat (deg)", "Lon (deg)", "Alt (km)", "Daedalus.Magnetic Latitude", "Daedalus.Magnetic Longitude", "Daedalus.MLT"] )
    with open( sourceFilename ) as CSVfileIN:        
        CSVreader = csv.reader( CSVfileIN )
        # locate the column numnbers of interest inside the csv file
        CSVheader = next( CSVreader )
        Time_idx     = CSVheader.index( "Time (UTCG)" ) #CSVheader.index( "Daedalus.EpochText" )
        Lat_idx      = CSVheader.index( "Lat (deg)" ) #CSVheader.index( "Daedalus.Latitude" )
        Lon_idx      = CSVheader.index( "Lon (deg)" ) #CSVheader.index( "Daedalus.Longitude" )
        Alt_idx      = CSVheader.index( "Alt (km)" ) #CSVheader.index( "Daedalus.Longitude" )
        # read the orbit file
        for row in CSVreader: # for each satellite position
            resultItems = list()
            # add the standard fields to the result file
            resultItems.append( row[Time_idx] )
            resultItems.append( row[Lat_idx] ) # Latitude is geodetic inside the orbit file
            resultItems.append( row[Lon_idx] ) 
            resultItems.append( row[Alt_idx] )
            # Calculate the extra fields
            try: # read time for the current orbit position 
                current_time = datetime.strptime(row[Time_idx], '%d %b %Y %H:%M:%S.%f')
            except:
                try:
                    current_time = datetime.strptime(row[Time_idx], '%b %d %Y %H:%M:%S.%f')
                except:
                    print( "ERROR - Wrong time format:", row[Time_idx] )
                    return
            # take care of time so that it is compatible with igrf
            #if current_time.year > 2024:  
            #    time_for_igrf = current_time - relativedelta(years=13)
            #else:
            #    time_for_igrf = current_time
            time_for_igrf = current_time
            MagneticLatitude, MagneticLongitude, MagneticLocalTime = Conversions.getMagneticProperties( time_for_igrf, float(row[Lat_idx]), float(row[Lon_idx]) )
            # add the extra fields to the result file
            resultItems.append( MagneticLatitude )
            resultItems.append( MagneticLongitude )
            resultItems.append( MagneticLocalTime )
            # write it
            CSVwriter.writerow( resultItems )
    # clean up
    CSVfileOUT.close()

    
'''
reads a csv orbit file with format:
    Time (UTCG)
    Lat (deg)
    Lon (deg)
    Alt (km)
and creates a new csv orbit file with the same fields plus the extra fields:
    Daedalus.Magnetic Latitude
    Daedalus.Magnetic Longitude
    Daedalus.MLT
The values are calculated for quasi-dipole approximation
'''
def AddMagneticCoordinates_QuasiDipole( sourceFilename, resultFilename ):
    input_file = pd.read_csv( sourceFilename ) 
    geod_lat = np.array(input_file["Lat (deg)"])
    geod_lon = np.array(input_file["Lon (deg)"])
    geod_alt = np.array(input_file["Alt (km)"])
    orbit_time = np.array(input_file["Time (UTCG)"])

    #DAY_OF_YEAR = np.array(input_file["Spacecraft1.DayOfYear"])

    re_geod_lat = []
    re_geod_lon = []
    re_geod_alt = []

    #apexpy needs at most 4 decimal digits precision
    for i in range (0,len(geod_alt)):
        re_geod_lat.append("{:.4f}".format(geod_lat[i]))
        if geod_lon[i] >180:
             geod_lon[i] = geod_lon[i] - 360 
        re_geod_lon.append("{:.4f}".format(geod_lon[i]))
        re_geod_alt.append("{:.4f}".format(geod_alt[i])) 

    geod_lat = np.array(re_geod_lat,dtype=float)
    geod_lon = np.array(re_geod_lon,dtype=float)

    MAG_LAT = []
    MAG_LON = []
    MLT = []
    #NEW_DOY = []
    new_time_string = []
    new_alt = []
    step=1 #refers to in how many mins we take sample from the orbit
    for i in range (0,len(geod_alt)):
        Time_value_String = orbit_time[i]
        # convert to datetime object from text....Jan 01 2015 00:00:01.000000000 to Jan 01 2015 00:00:01.0
        date_time,microseconds=Time_value_String.split('.')
        rounding = len(microseconds[:3])
        divisor=10**rounding
        new_micros=int(round(int(microseconds)/divisor,1))
        Time_value_String=date_time+'.'+str(new_micros)
        new_time_string.append(Time_value_String)
        dt_object = datetime.strptime(Time_value_String, '%d %b %Y %H:%M:%S.%f')
        k=ap.Apex(date=dt_object)
        # converting to magnetic coordinates
        mag_lat , mag_lon = k.geo2qd(geod_lat[i],geod_lon[i],re_geod_alt[i])
        MAG_LAT.append(mag_lat)
        MAG_LON.append(mag_lon)
        mlt = k.mlon2mlt(mag_lon,dt_object)
        MLT.append(mlt)
        #NEW_DOY.append(DAY_OF_YEAR[i])
        new_alt.append(re_geod_alt[i])
        #print (MAG_LON[i],MAG_LAT[i])
        
    # creating new columns for dataframe named = input_file
    input_file['Daedalus.Magnetic Latitude']=MAG_LAT
    input_file['Daedalus.Magnetic Longitude']=MAG_LON
    input_file['Daedalus.MLT']=MLT

    input_file.to_csv( resultFilename, index=None )

#new_file_dict = {"Epoch(UTCG)":new_time_string,'DAY_OF_YEAR':NEW_DOY, 'Lat_MAG (deg)':MAG_LAT, 'Lon_MAG (deg)':MAG_LON, 'Alt_MAG (km)':new_alt,'Magnetic Local Time (hours)':MLT}
#df = pd.DataFrame(new_file_dict)
#df.to_csv('Deadalus Orbits/DAED_ORB_LIFETIME_93deg_QD_'+str(step)+'mins_step.csv', index=False)
########################################################################################

