# This file contains source code for the Data Managent tab of DaedalusMaze GUI


from netCDF4 import Dataset  # https://unidata.github.io/netcdf4-python/netCDF4/index.html
from netCDF4 import MFDataset
from netCDF4 import stringtochar
from datetime import datetime
import SourceCode.DaedalusGlobals as DaedalusGlobals
import hashlib
import ipywidgets as widgets
import traceback
import numpy as np
import csv
import os


AuthPassword = "" # global - Necessary to remember password from the widgets.Password() widget.
PasswordBox = widgets.Password( ) # global
SpecificFilename = widgets.Text( )  # global

# Constructs the GUI which allows the user manage the NetCDF files.
def ConstructDataManagementPanel():
    # the main panel
    DataManagementPanel = widgets.VBox() 
    # for separating gui elements
    SeparatorLabel = widgets.HTML(value="<hr>")
    # for explaining the usage of this panel
    UsageLabel = widgets.HTML(value="<p style='line-height:1.4'><b>Usage:</b><br>Here you can select one or several orbit files which will be tranformed from csv to NetCDF files (.nc). If the NetCDF file already exists then it is overwritten. To start the process you have to type the correct password and press the 1st blue button. The 2nd blue button alters the state of the NetCDF file, making Calculated='no' so that all values will be recalculated during a simulation run.</p>")
    # for selecting orbits
    AllOrbitsLabel = widgets.HTML(value="<b>Select one or several orbits:</b>")
    AllOrbitsLabel.layout.margin = '20px 0px 0px 0px'
    AllOrbitFilenames = list()
    for f in os.listdir( DaedalusGlobals.CSV_Files_Path ):
        if f.endswith("csv"):
            AllOrbitFilenames.append( f[:-4] )
    global AllOrbitsSelector
    AllOrbitsSelector = widgets.SelectMultiple( options=AllOrbitFilenames, rows=8, )
    AllOrbitsSelector.layout.min_width='400px'
    # for entering the authorization password 
    PasswordLabel = widgets.HTML(value="<b>Type the authorization password:</b>")
    PasswordBox = widgets.Password()
    PasswordBox.observe(PasswordBox_onValueChange, names='value')
    # for selecting a specific file by typing its full-path name (in order to reset its state)
    SpecificFilenameLabel = widgets.HTML(value="<b>Type a specific full-path filename to reset:<b>")
    SpecificFilenameLabel.layout.margin = '20px 0px 0px 0px'
    # for selecting type of file according to the data stored inside it (in order to reset its state)
    AllStagesLabel = widgets.HTML(value="<b>Or select one/several stages of calculation to reset:</b>")
    AllStagesLabel.layout.margin = '20px 0px 0px 0px'
    AllStages = list()
    for f in os.listdir(DaedalusGlobals.Processed_Files_Path):
        if os.path.isdir(DaedalusGlobals.Processed_Files_Path+f):  
            AllStages.append( f )
    global AllStagesSelector
    AllStagesSelector = widgets.SelectMultiple( options=AllStages, rows=5, )
    AllStagesSelector.layout.min_width='400px'
    # construct the action buttons
    FillNETCDFbutton = widgets.Button(description='Create Orbit NetCDF files', button_style='primary')
    FillNETCDFbutton.layout.margin = '20px 0px 0px 0px'
    FillNETCDFbutton.layout.min_width='220px'
    FillNETCDFbutton.on_click( FillNETCDF_button_clicked )
    ResetNETCDFbutton = widgets.Button(description='Reset state of NetCDF files', button_style='primary')
    ResetNETCDFbutton.layout.margin = '20px 0px 0px 0px'
    ResetNETCDFbutton.layout.min_width='220px'
    ResetNETCDFbutton.on_click( ResetNETCDFbutton_clicked )
    
    # Put all together
    DataManagementPanel.children = [ UsageLabel, SeparatorLabel, PasswordLabel, PasswordBox, AllOrbitsLabel, AllOrbitsSelector, FillNETCDFbutton, SeparatorLabel, SpecificFilenameLabel, SpecificFilename, AllStagesLabel, AllStagesSelector, ResetNETCDFbutton ]
    return DataManagementPanel
    
# Event listener, handles the "Reset state of NetCDF files" button click
# It alters the flag attribute 'Calculated' inside the NetCDF file.
# When this attribute is "yes" then a simulation module may not bother recalculating all values
# and can just pass-along the file to the next module.
# Here, if the user decides that the already-calculated values are obsolete, she can change the flag to "no", 
# so that all values will re-calculated.
def ResetNETCDFbutton_clicked( b ):
    global AuthPassword
    global AllOrbitsSelector
    global SpecificFilename
    if len( SpecificFilename.value ) > 0:
        if os.path.isfile( SpecificFilename.value ):
            print( "Reseting " + SpecificFilename.value )
            ResetState_of_NETCDF_file( SpecificFilename.value ) 
        else:
            print( "File", SpecificFilename.value, "does not exist." )
    elif len(AllOrbitsSelector.value) <= 0:
        print( "Please select one or more orbit files." )
    elif len(AllStagesSelector.value) <= 0:
        print( "Please select one or more stages of calculation." )
    elif hashlib.sha224(AuthPassword.encode('utf-8')).hexdigest() == "658716d9784e508836425ed823de4440268f351f68fbf910c3b81696":
        n = 0
        for orbitname in AllOrbitsSelector.value:
            for stagename in AllStagesSelector.value:
                try:
                    for f in os.listdir(DaedalusGlobals.Processed_Files_Path + stagename + "/"):
                        if f.startswith( orbitname ):
                            fullname = DaedalusGlobals.Processed_Files_Path + stagename + "/" + orbitName + ".nc"
                            print( "Reseting " + fullname )
                            ResetState_of_NETCDF_file( fullname )
                            n = n + 1
                except Exception as err:
                    CDFroot.close() # so that the file will not be locked in case of error
                    print("ERROR: ", err)
                    traceback.print_exc()
        print( n, "files were reseted." ) 
    else:
        print( "Wrong Password.", hashlib.sha224(AuthPassword.encode('utf-8')).hexdigest() )
    
    
         
    
# Resets the state of a Daedalus NetCDF file.
# Arguments:
def ResetState_of_NETCDF_file( FullpathFilename):
    CDFroot =  Dataset( FullpathFilename, "a") # open file
    CDFroot.Calculated = "no" # reset state
    CDFroot.ResetTime = str(datetime.now()) # log 
    CDFroot.close() # close file
    
    
    
# Event listener, handles the "Add/Update orbits to NetCDF file" button click
# Allows the user to select one or many orbit files in csv format and batch convert them into NetCDF format.
def FillNETCDF_button_clicked( b ):
    global AuthPassword
    global AllOrbitsSelector
    if len(AllOrbitsSelector.value) <= 0:
        print( "Please select one or more orbit files." )
    elif hashlib.sha224(AuthPassword.encode('utf-8')).hexdigest() == "658716d9784e508836425ed823de4440268f351f68fbf910c3b81696":
        for orbitname in AllOrbitsSelector.value:
            try: 
                #### Create the NetCDF for the orbit
                print( "Creating " + DaedalusGlobals.NetCDF_Files_Path+orbitname+".nc" )
                CDFroot =  Dataset( DaedalusGlobals.NetCDF_Files_Path+orbitname+".nc", 'w') # open file
                CreateCDFstructure( CDFroot ) # create the structure
                FillCDF_withCSVdata( CDFroot, DaedalusGlobals.CSV_Files_Path+orbitname+".csv" ) # fill with data from the csv file
                CDFroot.close() # close file
                #### Create the NetCDF template file for the surface as well, so that it's structure is always up-to-date
                CDFroot =  Dataset( DaedalusGlobals.Temporary_Files_Path + "SurfaceTemplate.nc", 'w') # open file
                CreateCDFstructure( CDFroot ) # create the structure
                CDFroot.close() # close file
            except Exception as err:
                CDFroot.close() # so that the file will not be locked in case of error
                print("ERROR: ", err)
                traceback.print_exc()
    else:
        print( "Wrong Password." )

# Necessary to remember password from the widgets.Password() widget.
# Because this widget's "value" attribute is always empty.
def PasswordBox_onValueChange( change ):
    global AuthPassword
    AuthPassword = change['new']
        

################################################################################################
################################################################################################
################################################################################################

# Prints information about a NetCDF file. 
# Arguments: 
#    nc_fid: a reference to the NetCDF file, as returned from Dataset command.
def printCDFinfo( nc_fid ):
    print( "GENERAL INFO\n===========" )
    print( "Format: ", nc_fid.file_format )
    print( "Groups: ", nc_fid.groups )
    print( "\nKEYS\n===========" )
    print ( nc_fid.dimensions.keys() )
    print( "\nVARIABLES\n===========" )
    for i in nc_fid.variables:
        print(i, "\t", nc_fid.variables[i].dtype, end =" ")
        try:
            print("\t",nc_fid.variables[i].units, end =" ")
        except:
            print("\t","???", end =" ")
        try:
            print ( "\t","[" + nc_fid.variables[i].long_name + "]")
        except:
            print ( "\t","[" + "xxx" + "]")
        #print( nc_fid.variables[i].__dict__ ) # display all attributes of a variable    
    print( "\nDIMENSIONS\n===========" )        
    for dimobj in nc_fid.dimensions.values():
        print(dimobj)
        

# Creates an empty NetCDF file, assigning it an appropriate structure with all variables which will be used.
# Arguments: 
#    CDFroot: a reference to the NetCDF file, as returned from Dataset command.        
def CreateCDFstructure( CDFroot ):
    # create attributes
    CDFroot.Type = "orbit"
    CDFroot.Calculated = "no"
    CDFroot.CreationTime = str(datetime.now())
    CDFroot.ResetTime    = str(datetime.now())
    CDFroot.EditTime     = str(datetime.now())
    # create dimensions
    CDFroot.createDimension("altitude", None)
    CDFroot.createDimension("time",     None)
    CDFroot.createDimension("lat",      None)
    CDFroot.createDimension("lon",      None)
    # create variables
    CDFroot.createVariable("time",     "u8", ("time",))
    CDFroot.createVariable("lat",      "f8", ("time",))
    CDFroot.createVariable("lon",      "f8", ("time",))
    CDFroot.createVariable("altitude", "f8", ("time",))
    V = CDFroot.createVariable("ALFA", "f8", ("time",))
    V.description = "Aurora Characteristic Energy"
    V.units = "keV"
    V = CDFroot.createVariable("BARM", "f8", ("time",))
    V.description = "Mean Molecular Weight"
    V.units = "g/mol"
    V = CDFroot.createVariable("DEN", "f8", ("time",))
    V.description = "Total Density"
    V.units = "g/cm^3"
    V = CDFroot.createVariable("EFLUX", "f8", ("time",))
    V.description = "Aurora Energy Flux"
    V.units = "erg/cm^2/s"
    V = CDFroot.createVariable("HE", "f8", ("time",))
    V.description = "Helium"
    V.units = "mmr"
    V = CDFroot.createVariable("LAMDA_HAL", "f8", ("time",))
    V.description = "Hall Ion Drag Coefficient"
    V.units = "1/s"
    V = CDFroot.createVariable("LAMDA_PED", "f8", ("time",))
    V.description = "Pedersen Ion Drag Coefficient"
    V.units = "1/s"    
    V = CDFroot.createVariable("N4S", "f8", ("time",))
    V.description = ""
    V.units = "mmr"
    V = CDFroot.createVariable("NE", "f8", ("time",))
    V.description = "Electron Density"
    V.units = "cm^-3"
    V = CDFroot.createVariable("NO", "f8", ("time",))
    V.description = "Nitric Oxide"
    V.units = "mmr"
    V = CDFroot.createVariable("O1", "f8", ("time",))
    V.description = "Atomic Oxygen"
    V.units = "mmr"
    V = CDFroot.createVariable("O2", "f8", ("time",))
    V.description = "Molecular Oxygen"
    V.units = "mmr"
    V = CDFroot.createVariable("OP", "f8", ("time",))
    V.description = "O+ ION Density"
    V.units = "cm^-3"
    V = CDFroot.createVariable("PEDERSEN", "f8", ("time",))
    V.description = ""
    V.units = ""    
    V = CDFroot.createVariable("POTEN", "f8", ("time",))
    V.description = "Electric Potential"
    V.units = "Volts"
    V = CDFroot.createVariable("QAMIE", "f8", ("time",))
    V.description = ""
    V.units = ""    
    V = CDFroot.createVariable("QAURORA", "f8", ("time",))
    V.description = ""
    V.units = "cm^-3 s^-1"
    V = CDFroot.createVariable("QWIND", "f8", ("time",))
    V.description = ""
    V.units = ""    
    V = CDFroot.createVariable("sdtide", "f8", ("time",))
    V.description = "amplitudes and phases of semidiurnal tide"
    V.units = ""
    V = CDFroot.createVariable("TBAR", "f8", ("time",))
    V.description = ""
    V.units = ""    
    V = CDFroot.createVariable("Te", "f8", ("time",))
    V.description = "Electron Temperature"
    V.units = "Kelvin"
    V = CDFroot.createVariable("Ti", "f8", ("time",))
    V.description = "Ion Temperature"
    V.units = "Kelvin"
    V = CDFroot.createVariable("TEC", "f8", ("time",))
    V.description = "Total Electron Content"
    V.units = "1/cm^2"
    V = CDFroot.createVariable("TN", "f8", ("time",))
    V.description = "Neutral Temperature"
    V.units = "Kelvin"
    V = CDFroot.createVariable("UI_ExB", "f8", ("time",))
    V.description = "Zonal ExB Velocity"
    V.units = "cm/s"
    V = CDFroot.createVariable("VI_ExB", "f8", ("time",))
    V.description = " Meridional ExB Velocity"
    V.units = "cm/s"
    V = CDFroot.createVariable("VN", "f8", ("time",))
    V.description = "Neutral Meridional Wind (+North)"
    V.units = "cm/s"
    V = CDFroot.createVariable("UN", "f8", ("time",))
    V.description = "Neutral Zonal Wind (+East)"
    V.units = "cm/s"
    V = CDFroot.createVariable("WI_ExB", "f8", ("time",))
    V.description = "Vertical ExB Velocity"
    V.units = "cm/s"
    V = CDFroot.createVariable("WN", "f8", ("time",))
    V.description = "Neutral Vertical Wind (plus up)"
    V.units = "cm/s"
    V = CDFroot.createVariable("Z", "f8", ("time",))
    V.description = "Geopotential Height"
    V.units = "cm"
    V = CDFroot.createVariable("ZG", "f8", ("time",))
    V.description = "Geometric Height ZG"
    V.units = "cm"
    V = CDFroot.createVariable("ZMAG", "f8", ("time",))
    V.description = "ZMAG from pdynamo"
    V.units = "km"
    V = CDFroot.createVariable("ZGMID", "f8", ("time",))
    V.description = "Geometric Height at midpoints"
    V.units = "cm"
    
    

# Reads the contents of a CSV file and initializes a NetCDF file with its data.
# The data are time, lat, lon, elevation quadruples.
# Arguments: 
#    CDFroot: a reference to the NetCDF file, as returned from Dataset command.        
#    CSVfilename: the filename of the CSV file from which data will be incorporated into the NetCDF file.
def FillCDF_withCSVdata( CDFroot, CSVfilename ):
    CSVreader = csv.reader( open( CSVfilename ) ) # read the csv file. format: (time lat lon alt value)
    next( CSVreader ) # ignore the csv header
    csvData = [r[0:4] for r in CSVreader ] # read all data from the CSV file into a 2D array
    TimesAsStrings = np.array( csvData )[:,0] # parse Times from the CSV file
    Latitudes = np.array( csvData )[:,1].astype(np.float) # parse Latitudes from the CSV file
    Longitudes = np.array( csvData )[:,2].astype(np.float) # parse Longitudes from the CSV file
    Elevations = np.array( csvData )[:,3].astype(np.float) # parse Elevations from the CSV file
    # convert the time-strings of the CSV file to UTC milliseconds
    Times = list()
    for aTimeString in TimesAsStrings: # format in CSV file: 17 Mar 2015 00:00:00.000
        aTimeString = aTimeString[:-6] + " " + "+0000" # indicate explicity that this is UTC time ---- TO SHOW MY LOCAL TIME ZONE CORRECTION: from time import gmtime, strftime\n strftime("%z", gmtime())
        TimeObj = datetime.strptime( aTimeString, "%b %d %Y %H:%M:%S.%f %z" ) 
        Times.append( float( datetime.timestamp(TimeObj) ) * 1000 ) # store timestamp as milliseconds
    # store the parsed data into the NetCDF file
    CDFroot.variables["time"][:] = Times
    CDFroot.variables["lat"][:] = Latitudes
    CDFroot.variables["lon"][:] = Longitudes
    CDFroot.variables["altitude"][:] = Elevations
    CDFroot.variables["ALFA"][:] = np.random.random( len(Times) ) # assign random values for test purposes


'''
04/11/2019
This function creates a NetCDF file for a sphere surface around the Earth, 
with a structure which can be parsed and plotted correctly by Panoply.
We decided not to use this structure because it would require different methods for acessing data of the NetCDF file,
depending on wether it would be an orbit or a surface file. 
Instead we created a separate module Create_NetCDF_surface which creates an orbit-like structured NetCDF and updates
some of its attributes and variables, like Type and altitude.
def Create_NetCDF_surface( CDFroot, OrbitTemplate_CDFroot ):
    # create attributes
    CDFroot.Type = "surface"
    CDFroot.Calculated = "no"
    CDFroot.CreationTime = str(datetime.now())
    CDFroot.ResetTime = str(datetime.now())
    CDFroot.Time = "Jan 01 2028 00:00:01.000000000"
    CDFroot.Timestamp = ""
    CDFroot.Elevation = "400" # meters
    # create dimensions
    CDFroot.createDimension("lat", 36)
    CDFroot.createDimension("lon", 72)
    # create lat & lon variables
    CDFroot.createVariable("lat","f8",("lat",))
    CDFroot.createVariable("lon","f8",("lon",))
    # fill lat & lon variables
    CDFroot.variables["lat"][:] = np.arange(87.5, -92.5, -5)
    CDFroot.variables["lon"][:] = np.arange( -180,  180, 5)
    # create and fill all the rest variables
    for v in OrbitTemplate_CDFroot.variables:
        if v!="lat" and v!="lon" and v!="time" and v!="altitude":
            CDFroot.createVariable(v,"f4",("lat","lon",))
            CDFroot.variables[v][:, :] = np.random.random((36, 72)) # assign random values for test purposes
            CDFroot.variables[v][9, 21] = 2 # highlight New York for test purposes
# USAGE:
#### Create the NetCDF for the surface as well ( a sphere around Earth, with the same structure as the orbit file )
orbitTemplateCDFroot = Dataset(DaedalusGlobals.NetCDF_Files_Path+orbitname+".nc", 'r')
SurfaceCDFroot = Dataset(DaedalusGlobals.Temporary_Files_Path+"Surface.nc", 'w')
Create_NetCDF_surface( SurfaceCDFroot, orbitTemplateCDFroot )
orbitTemplateCDFroot.close()
SurfaceCDFroot.close()  
'''            

