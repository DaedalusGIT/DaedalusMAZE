import csv
import numpy as np
import SourceCode.DaedalusGlobals as DaedalusGlobals


'''
Creates a CSV file with the format [Time, Lat, Lon, Alt]
Time and Alt are fixed values given as arguments and Lat and Lon are produced for the whole globe accodring to steps specified.
CSVfilename: string, the csv filename to be written (overwrite mode)
fixedDatetimeString: string, the date-time fixed text with format as example: "17 Mar 2015 00:12:00.000". Fixed value for every csv entry. If this is an empty string then no file is created and empty string and zeros are returned.
fixedAltitude: float positive, kilometers above earth. Fixed value for every csv entry.
LatitudeStep: float positive, Latiude values range between 87.5 and -88.5. The step is how Lat should be incremented at each iteration.
LongitudeStep: float positive, Longitude values range between -180 and 180.  The step is how Lon should be incremented at each iteration.
RETURNS: the filename it has written (same as the argument)
         the altitude (same as the argument)
         the number of lines written
'''
def CreateCSV_Sphere( CSVfilename, fixedDatetimeString, fixedAltitude, LatitudeStep, LongitudeStep ):
   
    if len(fixedDatetimeString) == 0:
        return "", 0, 0 # <<<<

    linesWritten = 0 # count how many lines you write in the CSV in order to inform the user at the end.
    
    if CSVfilename is None  or  len(CSVfilename)==0:
        CSVfilename = DaedalusGlobals.Temporary_Files_Path + "sphere.csv"
    

    # open the new file as csv
    CSVfile = open( CSVfilename, mode='w')
    
    # write header of csv
    HeaderWriter = csv.writer(CSVfile, delimiter=',', quotechar='\"', quoting=csv.QUOTE_ALL, lineterminator="\n" )
    HeaderWriter.writerow(['Time (UTCG)', 'Lat (deg)', 'Lon (deg)', 'Alt (km)'])
    linesWritten = linesWritten + 1
    
    # write data of csv
    DataWriter = csv.writer(CSVfile, delimiter=',', quotechar='\"', quoting=csv.QUOTE_NONE, lineterminator="\n" )
    for lat in np.arange(87.5,  -88.5, -LatitudeStep):
        for lon in np.arange(-180.0,  180.0, LongitudeStep):
            DataWriter.writerow( [fixedDatetimeString, lat, lon, fixedAltitude])
            linesWritten = linesWritten + 1
    
    # finish
    CSVfile.close()
    return  CSVfilename, fixedAltitude,  linesWritten

# Example of Usage:
#s = CreateCSV_Sphere("sphere.csv", "17 Mar 2015 00:12:00.000", 120, 5, 5)
#print( str(s) + " lines written." )
