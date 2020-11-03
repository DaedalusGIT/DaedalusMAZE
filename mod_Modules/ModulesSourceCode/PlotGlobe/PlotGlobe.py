'''
The source code of this file is based on https://plot.ly/~empet/14813/heatmap-plot-on-a-spherical-map/#/
The main function of this file is PlotGlobe. The others are ancillary.

Notes
-----
We can add these plot types to Plotables:
    'area', 'bar', 'barpolar', 'box', 'candlestick', 'carpet', 'choropleth', 'cone', 'contour', 'contourcarpet', 
    'funnel', 'funnelarea', 'heatmap', 'heatmapgl', 'histogram', 'histogram2d', 'histogram2dcontour', 'isosurface', 
    'mesh3d', 'ohlc', 'parcats', 'parcoords', 'pie', 'pointcloud', 'sankey', 'scatter', 'scatter3d', 'scattercarpet', 
    'scattergeo','scattergl', 'scattermapbox', 'scatterpolar','scatterpolargl', 'scatterternary', 'splom', 'streamtube', 
    'sunburst', 'surface', 'table', 'violin', 'volume', 'waterfall'
Plotly colorscale names:
    'aggrnyl', 'agsunset', 'algae', 'amp', 'armyrose', 'balance', 'blackbody', 'bluered', 'blues', 'blugrn', 'bluyl', 'brbg',
    'brwnyl', 'bugn', 'bupu', 'burg', 'burgyl', 'cividis', 'curl', 'darkmint', 'deep', 'delta', 'dense', 'earth', 'edge', 'electric',
    'emrld', 'fall', 'geyser', 'gnbu', 'gray', 'greens', 'greys', 'haline', 'hot', 'hsv', 'ice', 'icefire', 'inferno', 'jet',
    'magenta', 'magma', 'matter', 'mint', 'mrybm', 'mygbm', 'oranges', 'orrd', 'oryel', 'peach', 'phase', 'picnic', 'pinkyl', 'piyg',
    'plasma', 'plotly3', 'portland', 'prgn', 'pubu', 'pubugn', 'puor', 'purd', 'purp', 'purples', 'purpor', 'rainbow', 'rdbu', 'rdgy',
    'rdpu', 'rdylbu', 'rdylgn', 'redor', 'reds', 'solar', 'spectral', 'speed', 'sunset', 'sunsetdark', 'teal', 'tealgrn', 'tealrose',
    'tempo', 'temps', 'thermal', 'tropic', 'turbid', 'twilight', 'viridis', 'ylgn', 'ylgnbu', 'ylorbr', 'ylorrd'

'''

from numpy import pi, sin, cos
import time
from netCDF4 import Dataset  # https://unidata.github.io/netcdf4-python/netCDF4/index.html
from netCDF4 import MFDataset
import plotly
import chart_studio.plotly as py # import plotly.plotly as py
import csv
import numpy as np 
import math
from scipy.io import netcdf  
from mpl_toolkits.basemap import Basemap
import warnings
import random

EarthRadius = 6378.137 # km - global, its value will not change through out the code
#EarthCircumference = 40007.863 # km
#pi = 3.14159265359

LatStep  =  5 # global - its value will change after parsing the surface file. It will be declared as global inside the fnctions.
LonStep  =  5 # global - its value will change after parsing the surface file. It will be declared as global inside the fnctions.
theSurfaceOpacity = 0.90

theSurfaceTimestep      = 1 # the timestep where the surface data to be plotted belong. This is parameter of the NetCDF file.
theSurfacePressureLevel = 0 # the pressure level where the surface data to be plotted belong. This is parameter of the NetCDF file.
theVectorsTimestep      = 1 # the timestep where the vector data to be plotted belong. This is parameter of the NetCDF file.
theVectorsPressureLevel = 0 # the pressure level where the vector data to be plotted belong. This is parameter of the NetCDF file.


# from plot.ly - convert degrees to radians
def degree2radians(degree):    
    return degree*pi/180

# from plot.ly - maps the points of coords (lon, lat) to points onto the sphere of radius radius

def mapping_map_to_sphere(lon, lat, radius=1):  
    lon=np.array(lon, dtype=np.float64)
    lat=np.array(lat, dtype=np.float64)
    lon=degree2radians(lon)
    lat=degree2radians(lat)
    xs=radius*cos(lon)*cos(lat)
    ys=radius*sin(lon)*cos(lat)
    zs=radius*sin(lat)
    return xs, ys, zs

def mapping_map_to_sphere_Azimuthal_equidistant_projection(lon, lat, radius=1):  
    lon=np.array(lon, dtype=np.float64)
    lat=np.array(lat, dtype=np.float64)
    lon=degree2radians(lon)
    lat=degree2radians(lat)
    Ro    = radius * ( pi/2 - lat )
    Theta = lon
    xs=Ro*cos(Theta)
    ys=Ro*sin(Theta)
    zs=radius*sin(lat)
    return xs, ys, zs


def mapping_map_to_blob(lon, lat, radius):  
    lon=np.array(lon, dtype=np.float64)
    lat=np.array(lat, dtype=np.float64)
    lon=degree2radians(lon)
    lat=degree2radians(lat)
    
    print( ">>> ", radius.shape, lon.shape, lat.shape )
    xs=radius*cos(lon)*cos(lat)
    ys=radius*sin(lon)*cos(lat)
    zs=radius*sin(lat)
    
    '''
    xs = np.zeros( lon.shape )
    ys = np.zeros( lon.shape )
    zs = np.zeros( lon.shape )
    for i in range(len(lon)):
        for j in range(len(lon[i])):
            xs[i][j] = (1+cos(lat[i][j]/5)/2)* radius*cos(lon[i][j])*cos(lat[i][j])
            ys[i][j] = (1+cos(lat[i][j]/5)/2)* radius*sin(lon[i][j])*cos(lat[i][j])
            zs[i][j] = (1+cos(lat[i][j]**2))* radius*sin(lat[i][j])
    '''
    return xs, ys, zs


# Make shortcut to Basemap object, not specifying projection type
m = Basemap() 

# from plot.ly - Functions converting coastline/country polygons to lon/lat traces
# pos arg 1. (poly_paths): paths to polygons
# pos arg 2. (N_poly): number of polygon to convert
def polygons_to_traces(poly_paths, N_poly):
    # init. plotting list
    lons=[]
    lats=[]
    for i_poly in range(N_poly):
        poly_path = poly_paths[i_poly]
        # get the Basemap coordinates of each segment
        coords_cc = np.array(
            [(vertex[0],vertex[1]) for (vertex,code) in poly_path.iter_segments(simplify=False)]
        )
        # convert coordinates to lon/lat by 'inverting' the Basemap projection
        lon_cc, lat_cc = m(coords_cc[:,0],coords_cc[:,1], inverse=True)
        lats.extend(lat_cc.tolist()+[None]) 
        lons.extend(lon_cc.tolist()+[None])
    return lons, lats

# from plot.ly - Function generating coastline lon/lat 
def get_coastline_traces():
    poly_paths = m.drawcoastlines().get_paths() # coastline polygon paths
    N_poly = 91  # use only the 91st biggest coastlines (i.e. no rivers)
    cc_lons, cc_lats= polygons_to_traces(poly_paths, N_poly)
    return cc_lons, cc_lats

# from plot.ly - Function generating country lon/lat 
def get_country_traces():
    poly_paths = m.drawcountries().get_paths() # country polygon paths
    N_poly = len(poly_paths)  # use all countries
    country_lons, country_lats= polygons_to_traces(poly_paths, N_poly)
    return country_lons, country_lats

###############

'''
Creates a 3D plot of an earth globe. Can plot a sphere surface and/or a satellite orbit. 
The surface and the orbit are colored according to the selected variable values in the data files. 
The data files can be of CSV o NetCDF4 format.
The same or different variables can be selected for plotting at the surface and the orbit.
The same or different colorscales can be selected for the surface and the orbit.
The same or different value range can be selected for the surface and the orbit by assigning the same or different colorbar titles.
Valid Colorscale Names: ‘Blackbody’, ‘Bluered’, ‘Blues’, ‘Earth’, ‘Electric’, ‘Greens’, ‘Greys’, ‘Hot’, ‘Jet’, ‘Picnic’, ‘Portland’, ‘Rainbow’, ‘RdBu’, ‘Reds’, ‘Viridis’, ‘YlGnBu’, ‘YlOrRd’
ARGUMENTS:
    SurfaceFilename: The file which contains data for a surface. If empty string then no surface will be plotted. 
                     CSV format: Time,Lat,Lon,Alt,value. Contains the data for the sphere surface. 
                     NetCDF format: see Daedalus User Manual
    SurfaceVariableToPlot: The name of the variable to read from the surface data file and plot upon the globe.
    SurfaceColorbarTitle: A title to display above the colorbar which refers to the surface data
    SurfaceColorscaleName: The name of the Colorscale to use for the surface data. 
                           In case an empty string is passed then a default HeatMap colorscale will be applied.
                           In case None is passed then all points will be black irrespective of value.
    OrbitFilename: counterpart of the corresponding argument for the Surface data
    OrbitVariableToPlot: counterpart of the corresponding argument for the Surface data 
    OrbitColorbarTitle: counterpart of the corresponding argument for the Surface data 
    OrbitColorscaleName: counterpart of the corresponding argument for the Surface data
    PlotTitle: A title which is displayed on the top of the plot.
OPTIONAL ARGUMENTS
    VectorsFilename: The NetCDF file which contains data for a vector plot. A cone will be displayed for each lat-lon pair.
    VectorsVariablesToPlot: a string of 3 NetCDF variables separated by comma. They represent the x,y,z components.
    VectorsColorbarTitle: A title to display above the colorbar which refers to the vector data
    VectorsColorscaleName: The name of the Colorscale to use for the vector data. 
    SurfaceOpacity: a float between 0 and 1 which defines how opaque or transparent the surface will be displayed
    LogScale: True or False. If True then log10(x) will be applied to all data values.
    VectorConeSize: how big the vector cones will be displayed
    VectorConeOpacity: a float between 0 and 1 which defines how opaque or transparent the vector cones will be displayed
    SurfaceMultiplicationFactor: the surface data values will be multiplied by this value before plotting
    OrbitMultiplicationFactor: the orbit data values will be multiplied by this value before plotting
    VectorsMultiplicationFactor: the vectors data values will be multiplied by this value before plotting
    OneSizeVectorCones: if True then all cones will have the same size, but different color according to the vector's length
    SurfaceTimestep: the timestep where the surface data to be plotted belong. This is parameter of the NetCDF file.
    SurfacePressureLevel: the pressure level where the surface data to be plotted belong. This is parameter of the NetCDF file.
    VectorsTimestep: the timestep where the vector data to be plotted belong. This is parameter of the NetCDF file.
    VectorsPressureLevel: the pressure level where the vector data to be plotted belong. This is parameter of the NetCDF file.
RETURNS: a string containing general information about the Data
'''
def PlotGlobe( SurfaceFilename, SurfaceVariableToPlot, SurfaceColorbarTitle, SurfaceColorscaleName, OrbitFilename, OrbitVariableToPlot, OrbitColorbarTitle, OrbitColorscaleName, PlotTitle, VectorsFilename="", VectorsVariablesToPlot="", VectorsColorbarTitle="", VectorsColorscaleName="Jet", SurfaceOpacity=0.90, LogScale=False, VectorConeSize=44, VectorConeOpacity=1, SurfaceMultiplicationFactor=1, OrbitMultiplicationFactor=1, VectorsMultiplicationFactor=1, OneSizeVectorCones=False, SurfaceTimestep= 1, SurfacePressureLevel=0, VectorsTimestep=1, VectorsPressureLevel=0 ):
    global theSurfaceTimestep
    global theSurfacePressureLevel
    global theVectorsTimestep
    global theVectorsPressureLevel
    
    result = ""
    startSecs = time.time()
    
    theSurfaceTimestep      = SurfaceTimestep
    theSurfacePressureLevel = SurfacePressureLevel
    theVectorsTimestep      = VectorsTimestep
    theVectorsPressureLevel = VectorsPressureLevel

    global theSurfaceOpacity
    theSurfaceOpacity = SurfaceOpacity
    
    
    # define colorscales for ploting
    colorscaleWhite=[[0.0, '#ffffff'], [1.0, '#ffffff']]
    colorscaleBlack=[[0.0, '#666666'], [1.0, '#666666']]
    if SurfaceColorscaleName is None:
        SurfaceColorscaleName = colorscaleBlack
    elif len( SurfaceColorscaleName ) == 0: 
        SurfaceColorscaleName=[[0.0, '#313695'], [0.07692307692307693, '#3a67af'], [0.15384615384615385, '#5994c5'], [0.23076923076923078, '#84bbd8'], [0.3076923076923077, '#afdbea'], [0.38461538461538464, '#d8eff5'], [0.46153846153846156, '#d6ffe1'], [0.5384615384615384, '#fef4ac'], [0.6153846153846154, '#fed987'], [0.6923076923076923, '#fdb264'], [0.7692307692307693, '#f78249'], [0.8461538461538461, '#e75435'], [0.9230769230769231, '#cc2727'], [1.0, '#a50026']]
    if OrbitColorscaleName is None:
        OrbitColorscaleName = colorscaleBlack
    elif len( OrbitColorscaleName ) == 0: 
        OrbitColorscaleName=[[0.0, '#313695'], [0.07692307692307693, '#3a67af'], [0.15384615384615385, '#5994c5'], [0.23076923076923078, '#84bbd8'], [0.3076923076923077, '#afdbea'], [0.38461538461538464, '#d8eff5'], [0.46153846153846156, '#d6ffe1'], [0.5384615384615384, '#fef4ac'], [0.6153846153846154, '#fed987'], [0.6923076923076923, '#fdb264'], [0.7692307692307693, '#f78249'], [0.8461538461538461, '#e75435'], [0.9230769230769231, '#cc2727'], [1.0, '#a50026']]
    if VectorsColorscaleName is None:
        VectorsColorscaleName = colorscaleBlack
    elif len( VectorsColorscaleName ) == 0: 
        VectorsColorscaleName=[[0.0, '#313695'], [0.07692307692307693, '#3a67af'], [0.15384615384615385, '#5994c5'], [0.23076923076923078, '#84bbd8'], [0.3076923076923077, '#afdbea'], [0.38461538461538464, '#d8eff5'], [0.46153846153846156, '#d6ffe1'], [0.5384615384615384, '#fef4ac'], [0.6153846153846154, '#fed987'], [0.6923076923076923, '#fdb264'], [0.7692307692307693, '#f78249'], [0.8461538461538461, '#e75435'], [0.9230769230769231, '#cc2727'], [1.0, '#a50026']]
        
        
    #### create the earth surface
    # construct the values for longitude and latitude
    global LatStep
    global LonStep
    if len(SurfaceFilename) > 0:        
        LatStep, LonStep = CalculateLatLonSteps( SurfaceFilename )
    lat = np.arange(   87.5,  -88.5, -LatStep )
    lon = np.arange( -180.0,  180.0,  LonStep )
    # To ensure color continuity we extend the lon list with [180] (its last value was lon[-1]=177.5). In this way we can identify lon=-180 with lon=180.
    # We do the same with latitudes. We extend both directions with [90] and [-90] in order to have values for the poles.
    clons = np.array( lon.tolist() + [180]        , dtype=np.float64)
    clats = np.array( [90] + lat.tolist() + [-90] , dtype=np.float64)
    clons, clats=np.meshgrid(clons, clats)
    earthX, earthY, earthZ = mapping_map_to_sphere(clons, clats, radius=EarthRadius)
    earthS = np.zeros(clons.shape, dtype=np.float64)
    EarthSurface=dict(type='surface',
        x = earthX, y = earthY, z = earthZ,
        colorscale = colorscaleWhite, surfacecolor = earthS,
        showscale = False, cmin=-20, cmax=20,
    )

    # Get list of of coastline, country, and state lon/lat and concatenate them
    coastline_lons, coastline_lats = get_coastline_traces()
    country_lons, country_lats = get_country_traces()
    boundaries_lons = coastline_lons+[None]+country_lons
    boundaries_lats = coastline_lats+[None]+country_lats
    boundaries_dataX, boundaries_dataY, boundaries_dataZ = mapping_map_to_sphere(boundaries_lons, boundaries_lats, radius = EarthRadius+0.1) # radius is slightly greater than 1 to ensure lines visibility
    # create the visual element of the world map
    Boundaries=dict(type='scatter3d',
        x=boundaries_dataX, y=boundaries_dataY, z=boundaries_dataZ,
        mode='lines', name= "", showlegend=False, line=dict(color='black', width=1)
    )    
    
    # define the layout of the plot
    theLayout = dict(
        title = PlotTitle,
        width=1200, height=1000,
        scene = dict(
            xaxis = dict( zeroline=False ), yaxis = dict( zeroline=False ), zaxis = dict( zeroline=False ),
            #xaxis_type="log", yaxis_type="log", zaxis_type="log",
            aspectratio=dict(x=1, y=1,z=1), camera=dict(eye=dict(x=1.20, y=1.20, z=1.20))
        )
    )
    
    #### read the surface Data file and find the Altitude of the surface described by the data file
    if len(SurfaceFilename) > 0:
        Altitude = CalculateAltitudeFromData( SurfaceFilename )
        
    #### read the Data file into an Array
    if len(SurfaceFilename) > 0:
        SurfaceData =  SurfaceFile_to_array( SurfaceFilename, SurfaceVariableToPlot ) 
        if SurfaceMultiplicationFactor != 1:
            SurfaceData = SurfaceData * SurfaceMultiplicationFactor
        if LogScale:
            for i in range(len(SurfaceData)):
                for j in range(len(SurfaceData[i])):
                    if SurfaceData[i][j] == 0:
                        SurfaceData[i][j] = SurfaceData[i][j]
                    elif SurfaceData[i][j] > 0:
                        SurfaceData[i][j] = math.log( SurfaceData[i][j], 10 )
                    else:
                        SurfaceData[i][j] = -1 * math.log( abs(SurfaceData[i][j]), 10 ) 
    else:
        SurfaceData = None
        
    
    #### read the Orbit Data file into an Array
    if len(OrbitFilename) > 0:
        OrbitData =  OrbitFile_to_array( OrbitFilename, OrbitVariableToPlot )
        if OrbitMultiplicationFactor != 1:
            OrbitData = OrbitData * OrbitMultiplicationFactor    
        if LogScale  and  (OrbitData is not None):
            for i in range(len(OrbitData)):
                if OrbitData[i][3] == 0:
                    OrbitData[i][3] = OrbitData[i][3]
                elif OrbitData[i][3] > 0:
                    OrbitData[i][3] = math.log( OrbitData[i][3], 10 )
                else:
                    OrbitData[i][3] = -1 * math.log( abs(OrbitData[i][3]), 10 )
    else:
        OrbitData = None


    #### Calculate min / max values for surface and orbit data 
    # calculate min/max for the surface data
    if len(SurfaceFilename) > 0:
        SurfaceMin = np.amin(SurfaceData)
        SurfaceMax = np.amax(SurfaceData)
    # calculate min/max for the orbit data
    if len(OrbitFilename) > 0:
        OrbitMin = np.amin( OrbitData[:,3] )
        OrbitMax = np.amax( OrbitData[:,3] )
    # combine if necessary
    if len(SurfaceFilename) > 0  and  len(OrbitFilename) > 0:
        if SurfaceVariableToPlot==OrbitVariableToPlot and SurfaceColorbarTitle==OrbitColorbarTitle: # user is ploting the same variable, so use a global min/max
            SurfaceMin = min( SurfaceMin, OrbitMin) 
            SurfaceMax = max( SurfaceMax, OrbitMax)
            OrbitMin = SurfaceMin
            OrbitMax = SurfaceMax
           
    # add all visual elements in a list and assign them to a figure
    Plotables = list()
    Plotables.append( EarthSurface )
    Plotables.append( Boundaries )
    if len(SurfaceFilename) > 0:
        Plotables.append( CreatePlotable_Surface(SurfaceData, Altitude, SurfaceMin, SurfaceMax, SurfaceColorscaleName, SurfaceColorbarTitle) )        
    if len(OrbitFilename) > 0:
        Plotables.append( CreatePlotable_Orbit(OrbitData, OrbitMin, OrbitMax, OrbitColorscaleName, OrbitColorbarTitle) )
        
    # For vectors
    if len(VectorsFilename) > 0:
        # construct Altitudes, Longtitudes, Altitudes 1D arrays
        V_Lats = list()
        V_Lons = list()
        V_Alts = list()
        VectorX = list()
        VectorY = list()
        VectorZ = list()
        latidx = 0
        lonidx = 0
        # parse variable names as given by user
        varXname, varYname, varZname = VectorsVariablesToPlot.split( ",", 2 )
        # open NectCDF file and read the variables
        CDFroot = Dataset( VectorsFilename, 'r' )        
        for i in np.arange(-88.75, 90, 2.5): # for i in np.arange(87.5, -92.5, -5.0):
            lonidx = 0
            a = 6378137; 
            e = 2.71828;
            aLat = degree2radians( i )
            N = EarthRadius #N = EarthRadius / math.sqrt(abs(1 - e*e * math.sin(aLat)*math.sin(aLat)));
            for j in np.arange(-180.0,  177.5, 2.5): # for j in np.arange(-180.0,  180.0, 5.0):
                aLon = degree2radians( j )
                V_Lats.append( (N+300) * math.cos(aLat) * math.cos(aLon) )
                V_Lons.append( (N+300) * math.cos(aLat) * math.sin(aLon) )
                V_Alts.append( (N+300) * math.sin(aLat) ) # V_Alts.append( ((1-e*e)*N+500) * math.sin(aLat) )
                
                # For Neutral Winds "UN,VN,WN", the ENU east,north,up coordinates have to be transformed to ECEF x,y,z coordinates
                if VectorsVariablesToPlot == "UN,VN,WN":
                    tmpX, tmpY, tmpZ = ENU_to_ECEF( aLat, aLon, CDFroot.variables["UN"][theVectorsTimestep,theVectorsPressureLevel,latidx,lonidx], CDFroot.variables["VN"][theVectorsTimestep,theVectorsPressureLevel,latidx,lonidx], CDFroot.variables["WN"][theVectorsTimestep,theVectorsPressureLevel,latidx,lonidx] )
                else:
                    tmpX = CDFroot.variables[varXname][theVectorsTimestep,theVectorsPressureLevel,latidx,lonidx]
                    tmpY = CDFroot.variables[varYname][theVectorsTimestep,theVectorsPressureLevel,latidx,lonidx]
                    tmpZ = CDFroot.variables[varZname][theVectorsTimestep,theVectorsPressureLevel,latidx,lonidx]
                    
                VectorX.append( tmpX )
                VectorY.append( tmpY )
                VectorZ.append( tmpZ )
                    
                lonidx = lonidx + 1
            latidx = latidx + 1
        CDFroot.close()
        VectorX = np.asarray(VectorX, dtype=np.float64)
        VectorY = np.asarray(VectorY, dtype=np.float64)
        VectorZ = np.asarray(VectorZ, dtype=np.float64)
        
        #print( len(VectorX), len(VectorY), len(VectorZ), len(V_Lats), len(V_Lons), len(V_Alts) )
        #np.set_printoptions(suppress=False)
        #np.set_printoptions(precision=12)
        
        # throw away some cones so that the plot becomes more clear
        #step = 3
        #V_Lats  = V_Lats[::step]
        #V_Lons  = V_Lons[::step]
        #V_Alts  = V_Alts[::step]
        #VectorX = VectorX[::step]
        #VectorY = VectorY[::step]
        #VectorZ = VectorZ[::step]
        
        if VectorsMultiplicationFactor != 1:
            VectorX = VectorX * VectorsMultiplicationFactor
            VectorY = VectorY * VectorsMultiplicationFactor
            VectorZ = VectorZ * VectorsMultiplicationFactor
        if LogScale:
            for i in range(len(VectorX)):
                if VectorX[i] == 0:
                    VectorX[i] = VectorX[i]
                elif VectorX[i] > 0:
                    VectorX[i] = math.log( VectorX[i], 10 )
                else:
                    VectorX[i] = -1 * math.log( abs(VectorX[i]), 10 ) 
            for i in range(len(VectorY)):
                if VectorY[i] == 0:
                    VectorY[i] = VectorY[i]
                elif VectorY[i] > 0:
                    VectorY[i] = math.log( VectorY[i], 10 )
                else:
                    VectorX[i] = -1 * math.log( abs(VectorX[i]), 10 ) 
            for i in range(len(VectorZ)):
                if VectorZ[i] == 0:
                    VectorZ[i] = VectorZ[i]
                elif VectorZ[i] > 0:
                    VectorZ[i] = math.log( VectorZ[i], 10 )
                else:
                    VectorZ[i] = -1 * math.log( abs(VectorZ[i]), 10 ) 
                    

        
        if OneSizeVectorCones: # zoro
            '''
            # TODO: calculate VectorsMin VectorsMax
            VectorsMin = 0
            VectorsMax = 100
            for i in range( 0, len(VectorX)-1, 1 ):
                #print( i, V_Lats[i], V_Lons[i],  V_Alts[i], " ~~ ", VectorX[i], VectorY[i], VectorZ[i] )
                SingleVectorCone = dict(type='cone', x=[V_Lats[i]],  y=[V_Lons[i]],  z=[V_Alts[i]],  u=[VectorX[i]], v=[VectorY[i]], w=[VectorZ[i]],
                    cmin = VectorsMin, cmax = VectorsMax,
                    sizemode='absolute', sizeref=VectorConeSize, opacity=VectorConeOpacity, colorscale=VectorsColorscaleName, anchor='tip'
                )
                Plotables.append( SingleVectorCone )
                '''
            
            VectorCones = dict(type='cone', x=V_Lats,  y=V_Lons,  z=V_Alts,   u=VectorX, v=VectorY, w=VectorZ,
                sizemode='scaled', sizeref=VectorConeSize, opacity=VectorConeOpacity, colorscale=VectorsColorscaleName, anchor='tip',
                colorbar=dict(thickness=20,len=0.75,ticklen=4,title=dict(text=VectorsColorbarTitle,side="top"),xanchor="center",x=0) 
            )        
            Plotables.append( VectorCones )            
        else :
            VectorCones = dict(type='cone', x=V_Lats,  y=V_Lons,  z=V_Alts,   u=VectorX, v=VectorY, w=VectorZ,
                sizemode='scaled', sizeref=VectorConeSize, opacity=VectorConeOpacity, colorscale=VectorsColorscaleName, anchor='tip',
                colorbar=dict(thickness=20,len=0.75,ticklen=4,title=dict(text=VectorsColorbarTitle,side="top"),xanchor="center",x=0) 
            )        
            Plotables.append( VectorCones )
    

    # plot all
    fig = dict( data=Plotables, layout=theLayout )
    plotly.offline.init_notebook_mode(connected=True)
    plotly.offline.iplot(fig)
    
    
    finishSecs = time.time()
    # be verbose
    #result = result + "Lon shape: " + str(lon.shape) + "  Lat shape: " + str(lat.shape) + "  Data shape: " + str(Data.shape) + "\n"
    if len(SurfaceFilename) > 0: result = result + "Surface Min: " + str(SurfaceMin) + "  Surface Max: " + str(SurfaceMax) + "\n"
    if len(OrbitFilename) > 0: result = result + "Orbit Min: " + str(OrbitMin) + "  Orbit Max: " + str(OrbitMax) + "\n"
    if len(OrbitFilename) > 0: result = result + "Number of Orbit positions: "  + str(len(OrbitData)) + "\n"    
        
        
    ########### Plot a Polar Chart of the surface-data and  orbit-data                   ###########
    ########### info from https://en.wikipedia.org/wiki/Azimuthal_equidistant_projection ###########
    PolarPlotables = list()

    # for the Surface
    if len(SurfaceFilename) > 0:
        SurfacePolarColors = SurfaceData.flatten()
        LATs = list()
        LONs = list()
        all_lats = np.arange(   -88.5,  87.5,  LatStep )
        all_lons = np.arange( -180.0,  180.0,  LonStep )
        for i in range( len(all_lats) ):
            for j in range( len(all_lons) ):
                LATs.append( all_lats[i] )
                LONs.append( all_lons[j] )
        Surface_Lats = np.array( LATs , dtype=np.float64)
        Surface_Lons = np.array( LONs , dtype=np.float64)
        Surface_Ro    = np.zeros( len(Surface_Lons) )
        Surface_Theta = np.zeros( len(Surface_Lons) )
        for i in range( len( Surface_Lons ) ):
            Surface_Ro[i]    = ( pi/2 - Surface_Lats[i] )
            Surface_Theta[i] = Surface_Lons[i] - 90
            #if  Surface_Lats[i]+33.75<5  and  Surface_Lons[i]+70<5:
            #    SurfacePolarColors[i] = 1000
            #    SurfaceMax = 1000
        SurfacePolarScatter = dict( 
            type = "scatterpolar", mode = "markers", 
            r = Surface_Ro,  theta = Surface_Theta,
            marker = dict( size=10, color = SurfacePolarColors, colorscale=SurfaceColorscaleName, opacity = 0.80, 
                          cmin = SurfaceMin, cmax = SurfaceMax,
                colorbar=dict(thickness=20, len=0.75, ticklen=4, title=dict(text=SurfaceColorbarTitle,side="top"), xanchor="center", x=0) # x must be between [-2, 3]
            ),
        ) 
        
        #PolarPlotables.append( SurfacePolarScatter )

        PolarPlotables.append( CreatePlotable_PolarSurface (SurfaceData, SurfaceMin, SurfaceMax, SurfaceColorscaleName, SurfaceColorbarTitle) )
    
    # for theEarth
    if len(SurfaceFilename) > 0  or  len(OrbitFilename) > 0 :
        Boundaries_Ro    = np.zeros( len(boundaries_lons) )
        Boundaries_Theta = np.zeros( len(boundaries_lons) )
        for i in range( len( boundaries_lons ) ):
            if boundaries_lats[i] is None:
                Boundaries_Ro[i] = 0
            else:
                Boundaries_Ro[i]    = ( pi/2 - boundaries_lats[i] )
            if boundaries_lons[i] is None:
                Boundaries_Theta[i] = 0
            else:
                Boundaries_Theta[i] = boundaries_lons[i] - 90
        
        # Get list of of coastline, country, and state lon/lat and concatenate them 
        coastline_lons, coastline_lats = get_coastline_traces()
        country_lons, country_lats = get_country_traces()   
        boundaries_lons = coastline_lons+[None]+country_lons
        boundaries_lats = coastline_lats+[None]+country_lats
        # keep only 60-90 degrees
        for i in range( len(boundaries_lats) ):
            if (boundaries_lats[i] is not None) and boundaries_lats[i] < 55:
                boundaries_lats[i] = None
        # adjust latitudes for polar plot
        for i in range( len(boundaries_lats) ):
            if boundaries_lats[i] is None:
                boundaries_lats[i] = None
            elif boundaries_lats[i] >= 0:
                boundaries_lats[i] = 45 + boundaries_lats[i] / 2
            elif boundaries_lats[i] < 0:
                boundaries_lats[i] = 45 + boundaries_lats[i]/2
        boundaries_dataX, boundaries_dataY, boundaries_dataZ = mapping_map_to_sphere_Azimuthal_equidistant_projection(boundaries_lons, boundaries_lats, radius=EarthRadius)
                
        XXX = boundaries_dataX 
        YYY = boundaries_dataY 
        ZZZ = np.zeros(len(boundaries_dataZ))
            
        BoundariesPolarScatter3D = dict( 
            type = "scatter3d", mode = "lines", opacity = 1, showlegend=False, marker = dict( size=2, color="black"),
            x=XXX, y=YYY, z=ZZZ,
        )          
        PolarPlotables.append( BoundariesPolarScatter3D )
    
    # for the lines
    if len(SurfaceFilename) > 0:    
        line_lons = np.arange(-180, 190, 10).tolist()
        line_lats = np.full( len(line_lons), 0 ).tolist()
        for i in range( 0, int(EarthRadius/6), int(EarthRadius/18) ):
            LinesX, LinesY, LinesZ = mapping_map_to_sphere_Azimuthal_equidistant_projection( line_lons, line_lats, radius=i )
            LinesZ = np.zeros( LinesZ.shape )
            Lines = dict( 
                type = "scatter3d", mode = "lines", opacity = 0.7, marker = dict( size=2, color="gray"), showlegend=False,
                x=LinesX, y=LinesY, z=LinesZ,
            )
            PolarPlotables.append( Lines )
        ## add a label for each line indicating the latitude number
        labelsXpos = np.arange( 140+555, 2000, 555 )
        labels = dict( 
            type = "scatter3d", mode = "text", opacity = 0.7, showlegend=False,
            x=labelsXpos, y=[0,0,0], z=[0,0,0], text=["80", "70", "60"],
        )
        PolarPlotables.append( labels )
            
    # for the Orbit  
    if len(OrbitFilename) > 0: 
        Orbit_Lats  = OrbitData[:,0]
        Orbit_Lons  = OrbitData[:,1]
        # keep only 60-90 degrees
        for i in range( len(Orbit_Lats) ):
            if Orbit_Lats[i] < 55:
                Orbit_Lats[i] = None
        for i in range( len(Orbit_Lats) ):
            if Orbit_Lats[i] is None:
                Orbit_Lats[i] = None
            elif Orbit_Lats[i] >= 0:
                Orbit_Lats[i] = 45 + Orbit_Lats[i] / 2
            elif Orbit_Lats[i] < 0:
                Orbit_Lats[i] = 45 + Orbit_Lats[i]/2

        XXX, YYY, ZZZ = mapping_map_to_sphere_Azimuthal_equidistant_projection(Orbit_Lons, Orbit_Lats, radius=EarthRadius)
        ZZZ = np.zeros(len(ZZZ))
        OrbitPolarScatter3D = dict(
            type = "scatter3d", mode = "markers", opacity = 1, showlegend=False, x=XXX, y=YYY, z=ZZZ,
            marker = dict( size=2, color=OrbitData[:,2], cmin = OrbitMin, cmax = OrbitMax, colorscale=OrbitColorscaleName,
                colorbar=dict(thickness=20, len=0.75, ticklen=4, title=dict(text=OrbitColorbarTitle,side="top"), xanchor="center", x=0) # x must be between [-2, 3]                          
            ),
        )          
        PolarPlotables.append( OrbitPolarScatter3D )
    
    # for the Orbit
    '''
    if len(OrbitFilename) > 0:    
        Orbit_Lats  = OrbitData[:,0]
        Orbit_Lons  = OrbitData[:,1]
        Orbit_Ro    = np.zeros( len(Orbit_Lats) )
        Orbit_Theta = np.zeros( len(Orbit_Lats) )
        OrbitPolarColors = OrbitData[:,2]
        for i in range( len( OrbitData ) ):
            Orbit_Ro[i]    = ( pi/2 - Orbit_Lats[i] )
            Orbit_Theta[i] = Orbit_Lons[i] - 90
        OrbitPolarScatter = dict( 
            type = "scatterpolar", mode = "markers", 
            r = Orbit_Ro,  theta = Orbit_Theta,
            marker = dict( size=4, color = OrbitPolarColors, colorscale = "Jet", opacity = 0.80, 
                          cmin = OrbitMin, cmax = OrbitMax,
                colorbar=dict(thickness=20, len=0.75, ticklen=4, title=dict(text=SurfaceColorbarTitle,side="top"), xanchor="center", x=0) # x must be between [-2, 3]
            ),
        )
        #PolarPlotables.append( OrbitPolarScatter )
    '''
        
    #### Plot the polar chart 
    if len(PolarPlotables) > 0:
        fig = dict( data=PolarPlotables, layout=theLayout )
        plotly.offline.init_notebook_mode(connected=True)
        plotly.offline.iplot(fig)
        
    result = result + "Plot creation duration: " + str(finishSecs-startSecs) + " seconds.\n"
    return result;


###################################################################################################
###################################################################################################
###################################################################################################

#### parse the data file and find the Altitude of the surface described
# Arguments:
#    SurfaceFilename: the filename containing the data for the surface. Can be in csv or netCDF format.
# Returns: the Altitude (float)
def CalculateAltitudeFromData( SurfaceFilename ):
    result = 0
    if len(SurfaceFilename) == 0:
        result = 0
    elif SurfaceFilename.endswith( ".csv" ):
        with open(SurfaceFilename) as CSVfile:
            CSVreader = csv.reader( CSVfile )
            next( CSVreader ) # ignore the csv header
            result = float(   next( CSVreader )[3]   )
    else:
        CDFroot = Dataset( SurfaceFilename, 'r' )
        try:
            result = CDFroot.variables["altitude"][0]
            if result==0: 
                result = 2
        except:
            try:
                # convert pressure levels to altitudes
                result = EarthRadius + CDFroot.variables["ZG"][theSurfaceTimestep,theSurfacePressureLevel,:,:] /100000# Geometric Height in cm. float ZG(time=17, ilev=29, lat=36, lon=72)
                result[10][40] = result[10][40]+800 # Athens
                if len(result) > 36: # TODO explain or correct
                    result = result[0::2, 0::2]
                print( "Converted pressure levels to altitudes." )
                result = np.mean( result ) - EarthRadius # TODO: make it work as blob
            except:
                result = 100
        CDFroot.close()
    ##
    return result


# Calculates the step between succesive values of Latitudes and Longitudes as stored inside the surface file.
# The steps are absolute values, always positive.
# Different step values lead to smaller or bigger arrays which hold the surface data and display the plot.
# This function enables the automation of such an adjustment.
# Arguments:
#    SurfaceFilename: the filename containing the data for the surface. Can be in csv or netCDF format.
# Returns 2 floats: the LatitudeStep and the LongitudeStep
def CalculateLatLonSteps( SurfaceFilename ):
    theLatStep = -5
    theLonStep =  5
    if len(SurfaceFilename) == 0:
        return
    elif SurfaceFilename.endswith( ".csv" ):
        with open(SurfaceFilename) as CSVfile:
            CSVreader = csv.reader( CSVfile ) # read the csv file. format: (time lat lon alt value)
            next( CSVreader ) # ignore header
            # read data
            RawData =  [r[1:] for r in CSVreader] # read the rest file into a 2D list. Ignore 1st column (time).
            RawData = np.array( RawData, dtype=np.float64 ) # convert to 2D numpy array
            # calculate Latitude step
            Lats = np.unique( RawData[:,0:1] )
            theLatStep = abs( Lats[0] - Lats[1] )
            # calculate Longitude step
            Lons = np.unique( RawData[:,1:2] )
            theLonStep = abs( Lons[0] - Lons[1] )
    else: 
        CDFroot = Dataset( SurfaceFilename, 'r' )
        if "(unlimited)" in str(CDFroot.dimensions["lat"]):
            # calculate Latitude step
            Lats = CDFroot.variables["lat"][:]
            for i in range(0, len(Lats)):
                Lats[i] = round( Lats[i], 2 )
            Lats = np.unique( Lats )
            theLatStep = abs( Lats[0] - Lats[1] )
            # calculate Longitude step
            Lons = CDFroot.variables["lon"][:]
            for i in range(0, len(Lons)):
                Lons[i] = round( Lons[i], 2 )
            Lons = np.unique( Lons )
            theLonStep = abs( Lons[0] - Lons[1] )            
        else:
            # calculate Latitude step
            Lats = CDFroot.variables["lat"][:]
            theLatStep = abs( Lats[0] - Lats[1] )
            # calculate Longitude step
            Lons = CDFroot.variables["lon"][:]
            theLonStep = abs( Lons[0] - Lons[1] )
        CDFroot.close()
        theLatStep = 5 # TODO explain or correct
        theLonStep = 5 # TODO explain or correct
    return theLatStep, theLonStep




#### read the values of the surface data file into a 2D Array
# Arguments:
#    SurfaceFilename: the filename containing the data for the surface. Can be in csv or netCDF format.
#    VariableName: the name of the variable which will be plotted. The file may contain values for several variables.
# Returns: a 2D array containing the data to be plotted. The item Data[0,0] represents the value at Lat=87.5 and Lon=-180 etc.
def SurfaceFile_to_array( SurfaceFilename, VariableName ):
    global LatStep
    global LonStep
    if len(SurfaceFilename) == 0:
        DataGrid = None
    elif SurfaceFilename.endswith( ".csv" ):
        with open(SurfaceFilename) as CSVfile:
            CSVreader = csv.reader( open( SurfaceFilename ) ) # read the csv file. format: (time lat lon alt value)
            # read csv header to resolve the column number which corresponds to VariableName
            VariableColumn = 3 # use the 1st column containing variable values as default
            HeaderTitles = next( CSVreader )
            for i in range(0, len(HeaderTitles)):
                if HeaderTitles[i] == VariableName:
                    VariableColumn = i - 1 # we subtract by 1 because we will throw away the time column (the 1st column)
            # read the rest file into a 2D array. 
            RawData =  [r[1:] for r in CSVreader] # Ignore 1st column (time). -> get lat,lon,elevation,values
            RawData = np.array( RawData, dtype=np.float64 ) # convert to 2D numpy array
            # construct a 2D array like a grid wirh Lat and Lon axis with the variable's values inside
            Lat = np.arange(   87.5,  -88.5,  -LatStep )
            Lon = np.arange( -180.0,  180.0,   LonStep )
            DataGrid = np.zeros( (len(Lat), len(Lon)), dtype=np.float64 )
            for i in range(0, len(Lat)):
                for j in range(0, len(Lon)):
                    for x in range(0, len(RawData)):
                        if( RawData[x,0] == Lat[i]  and  RawData[x,1] == Lon[j] ):
                            DataGrid[i, j] = RawData[x, VariableColumn]
                            break
    else: 
        CDFroot = Dataset( SurfaceFilename, 'r' )
        if "(unlimited)" in str(CDFroot.dimensions["lat"]): # surface's file format is the same as the orbit's
            # read the NetCDF file into a 2D array with columns lat,lon,elevation,values
            RawData = np.vstack( (CDFroot.variables["lat"][:], CDFroot.variables["lon"][:], CDFroot.variables["altitude"][:], CDFroot.variables[VariableName][:] ) ).T
            # Round lat and lon columns to remove small variations which trick Panoply to see different values and plot the surface
            for i in range(0, len(RawData)):
                RawData[i, 0] = round( RawData[i, 0], 2 )
                RawData[i, 1] = round( RawData[i, 1], 2 )
            # construct a 2D array like a grid wirh Lat and Lon axis with the variable's values inside
            Lat = np.arange(   87.5,  -88.5,  -LatStep )
            Lon = np.arange( -180.0,  180.0,   LonStep )
            DataGrid = np.zeros( (len(Lat), len(Lon)), dtype=np.float64 )
            for i in range(0, len(Lat)):
                for j in range(0, len(Lon)):
                    for x in range(0, len(RawData)):
                        if( RawData[x,0] == Lat[i]  and  RawData[x,1] == Lon[j] ):
                            DataGrid[i, j] = RawData[x, 3]
                            break
        else: # the Panoply-oriented format for surfaces is used
            DataGrid = np.asarray(CDFroot.variables[VariableName][theSurfaceTimestep,theSurfacePressureLevel,:,:] , dtype=np.float64)
            if len(DataGrid) > 36: # TODO explain or correct
                DataGrid = DataGrid[0::2, 0::2]
        ##
        CDFroot.close()
    ##
    return DataGrid



#### read the values of the orbit data file into a 1D Array
#    SurfaceFilename: the filename containing the data for the surface. Can be in csv or netCDF format.
#    VariableName: the name of the variable which will be plotted. The file may contain values for several variables.
# Returns: a 2D array containing the data to be plotted, with columns Lat,Lon,Elevation,Value
def OrbitFile_to_array( OrbitFilename, VariableName ):
    if len(OrbitFilename) == 0:
        OrbitData = None
    elif OrbitFilename.endswith( ".csv" ):    
        orbitCSVreader = csv.reader( open( OrbitFilename ) )
        next( orbitCSVreader ) # ignore the csv header
        OrbitData =  [r[1:5] for r in orbitCSVreader if len(r[4])>0] # read the rest CSV file at once into a 2D list. Ignore first column which is the time and take only the first value in case there are more. Ignore lines without value, as well.
        OrbitData = np.array( OrbitData, dtype=np.float64 ) # convert to 1D numpy array
    else:
        CDFroot = Dataset( OrbitFilename, 'r' )
        OrbitData = np.vstack( (CDFroot.variables["lat"][:], CDFroot.variables["lon"][:], CDFroot.variables["altitude"][:], CDFroot.variables[VariableName][:] ) ).T
        CDFroot.close()
    ##
    return OrbitData

###################################################################################################


# Creates a dictionary containing all information needed to plot a spherical surface around the earth
# with data from the surface-data-file.
# Arguments:
#     Data: the data in 2D array to be plotted. The item Data[0,0] represents the value at Lat=87.5 and Lon=-180. This array can be cosntructed by the function SurfaceFile_to_array from a csv or netCDF file.
#     Altitude: the altitude of the sphere above sea level
#     MinValue: the minimum value of the data. Needed because there may be a global minimum calculated upon both the orbit and the surface data.
#     MaxValue: the maximum value of the data
#     ColorScale: the colorscale name. This changes the colors which will be displayed. The colors will correspond to the values and will ne in range between the min and max values
#     ColorbarTitle: the title which will be shows above the color-bar. This bar displays the colors and their corresponding values 
def CreatePlotable_Surface( Data, Altitude, MinValue, MaxValue , ColorScale, ColorbarTitle):
    # construct the values for longitude and latitude
    global LatStep
    global LonStep
    lat = np.arange(   87.5,  -88.5,  -LatStep ) # lat = np.arange(   87.5,  -88.5,  int(-180 / len(Data   ) ) )
    lon = np.arange( -180.0,  180.0,   LonStep ) # lon = np.arange( -180.0,  180.0,  int( 360 / len(Data[0]) ) ) 
    
    #### SHIFTING
    # Shift 'lon' from [0,360] to [-180,180]
    tmp_lon = np.array([lon[n]-360 if l>=180 else lon[n] for n,l in enumerate(lon)])  # => [0,180]U[-180,2.5]
    i_east, = np.where(tmp_lon>=0)  # indices of east lon
    i_west, = np.where(tmp_lon<0)   # indices of west lon
    lon = np.hstack((tmp_lon[i_west], tmp_lon[i_east]))  # stack the 2 halves
    # Correspondingly, shift the Data array
    Data_ground = np.array(Data)
    try:  # TODO explain or correct
        Data = np.hstack((Data_ground[:,i_west], Data_ground[:,i_east]))
    except:
        print( "Data NOT shifted." )
    
    # To ensure color continuity we extend the lon list with [180] (its last value was lon[-1]=177.5). In this way we can identify lon=-180 with lon=180.
    # We do the same with latitudes. We extend both directions with [90] and [-90] in order to have values for the poles.
    clons = np.array( lon.tolist() + [180]        , dtype=np.float64)
    clats = np.array( [90] + lat.tolist() + [-90] , dtype=np.float64)
    clons, clats=np.meshgrid(clons, clats)
    # Map the meshgrids clons, clats onto the sphere
    if isinstance(Altitude, np.ndarray): 
        cAltitude = np.array( Altitude )
        cAltitude[0,0] = cAltitude[0,0] + 500
        cAltitude = np.vstack(( cAltitude, cAltitude[35] ))
        cAltitude = np.vstack(( cAltitude, cAltitude[ 0] ))

        Altitude = np.hstack((cAltitude[:,i_west], cAltitude[:,i_east]))
        
        tmp_array = (Altitude[:,0]+Altitude[:,71])/2
        tmp_array = np.reshape( tmp_array, (38,1) )
        Altitude = np.hstack(( Altitude, tmp_array ))

        print( "cAltitude.shape =", cAltitude.shape )
        print( "Altitude.shape =", Altitude.shape )
        print( "Data.shape =", Data.shape )
        XS, YS, ZS = mapping_map_to_blob(clons, clats, Altitude)
    else:
        XS, YS, ZS = mapping_map_to_sphere(clons, clats, radius=EarthRadius+Altitude)
    # for color continuity
    nrows, ncolumns=clons.shape
    DATA = np.zeros(clons.shape, dtype=np.float64)
    DATA[1:nrows-1, :ncolumns-1] = np.copy(np.array(Data,  dtype=np.float64)) # ignore the extended values
    DATA[1:nrows-1,  ncolumns-1] = np.copy(Data[:, 0]) # ignore the extended values
    
    # Create a sphere above earth which will be colored in accordance with the data
    DataSphere=dict(type='surface', x = XS,  y = YS,  z = ZS,
        colorscale = ColorScale, surfacecolor = DATA, opacity = theSurfaceOpacity,
        cmin = MinValue, cmax = MaxValue,
        colorbar=dict(thickness=20, len=0.75, ticklen=4, title=ColorbarTitle),
    )
    #if len(ColorbarTitle) == 0 : DataSphere["showscale"] = False

    return DataSphere
    

###########################################################################################################

def CreatePlotable_PolarSurface( Data, MinValue, MaxValue, ColorScale, ColorbarTitle):
    # construct the values for longitude and latitude
    global LatStep
    global LonStep
    lat = np.arange(   87.5,  -88.5,  -LatStep ) # lat = np.arange(   87.5,  -88.5,  int(-180 / len(Data   ) ) )
    lon = np.arange( -180.0,  180.0,   LonStep ) # lon = np.arange( -180.0,  180.0,  int( 360 / len(Data[0]) ) ) 
    
    #### SHIFTING
    # Shift 'lon' from [0,360] to [-180,180]
    tmp_lon = np.array([lon[n]-360 if l>=180 else lon[n] for n,l in enumerate(lon)])  # => [0,180]U[-180,2.5]
    i_east, = np.where(tmp_lon>=0)  # indices of east lon
    i_west, = np.where(tmp_lon<0)   # indices of west lon
    lon = np.hstack((tmp_lon[i_west], tmp_lon[i_east]))  # stack the 2 halves
    # Correspondingly, shift the Data array
    Data_ground = np.array(Data)
    try:  # TODO explain or correct
        Data = np.hstack((Data_ground[:,i_west], Data_ground[:,i_east]))
    except:
        print( "Data NOT shifted." )
        
    # To ensure color continuity we extend the lon list with [180] (its last value was lon[-1]=177.5). In this way we can identify lon=-180 with lon=180.
    # We do the same with latitudes. We extend both directions with [90] and [-90] in order to have values for the poles.
    clons = np.array( lon.tolist() + [180]        , dtype=np.float64)
    clats = np.array( [90] + lat.tolist() + [-90] , dtype=np.float64)
    

    # keep only 60-90 degrees
    for i in range( len(clats) ):
        if clats[i] < 55:
            clats[i] = None
    # manipulate for polar plot    
    for i in range( len(clats) ):
        if clats[i] is None:
            clats[i] = None
        elif clats[i] >= 0:
            clats[i] = 45 + clats[i] / 2
        elif clats[i] < 0:
            clats[i] = 45 + clats[i]/2
            

    # create 2d mesh from the 1D lats and lons
    clons, clats = np.meshgrid(clons, clats)            
            
    # Map the meshgrids clons, clats onto the sphere
    XS, YS, ZS = mapping_map_to_sphere_Azimuthal_equidistant_projection(clons, clats, radius=EarthRadius)
    ZS = np.full( ZS.shape, 0 )
            
    # for color continuity
    nrows, ncolumns=clons.shape
    DATA = np.zeros(clons.shape, dtype=np.float64)
    DATA[1:nrows-1, :ncolumns-1] = np.copy(np.array(Data,  dtype=np.float64)) # ignore the extended values
    DATA[1:nrows-1,  ncolumns-1] = np.copy(Data[:, 0]) # ignore the extended values
    
    # flatten to polar plain
    ##ZS = np.full( ZS.shape, 0 )
    # Create the mesh in polar coordinates and compute corresponding Z.
    ##r = np.linspace(0, 90, 73) #r = np.linspace(0, 90, 73)
    ##p = np.linspace(-np.pi, np.pi, 39) # p = np.linspace(0, 2*np.pi, 38)
    ##R, P = np.meshgrid(r, p)
    # Express the mesh in the cartesian system.
    ##XS, YS = R*np.cos(P), R*np.sin(P)
    
    # Create a polar surface colored in accordance with the data
    DataSphere=dict(type='surface', x = XS,  y = YS,  z = ZS,
        colorscale = ColorScale, surfacecolor = DATA, opacity = 0.9,
        cmin = MinValue, cmax = MaxValue,
        colorbar=dict(thickness=20, len=0.75, ticklen=4, title=ColorbarTitle),
    )
    #if len(ColorbarTitle) == 0 : DataSphere["showscale"] = False

    return DataSphere



# Creates a dictionary containing all information needed to Plotly in order to plot a spherical surface around the earth
# with data from the surface-data-file.
# Arguments:
#     Data: the data in 2D array to be plotted. The item Data[0,0] represents the value at Lat=87.5 and Lon=-180. This array can be cosntructed by the function SurfaceFile_to_array from a csv or netCDF file.
#     MinValue: the minimum value of the data. Needed because there may be a global minimum calculated upon both the orbit and the surface data.
#     MaxValue: the maximum value of the data
#     ColorScale: the colorscale name. This changes the colors which will be displayed. The colors will correspond to the values and will ne in range between the min and max values
#     ColorbarTitle: the title which will be shows above the color-bar. This bar displays the colors and their corresponding values 
def CreatePlotable_Orbit( OrbitData, MinValue, MaxValue, ColorScale, ColorbarTitle ):    
    # process orbit-data in order to plot them in relation with the globe: Convert lat,lon,alt to x,y,z
    for i in range(0, len(OrbitData)): OrbitData[i, 2] = OrbitData[i, 2] + EarthRadius # form radius from earth center using altitude
    tmplon = degree2radians( OrbitData[:,1] )  
    tmplat = degree2radians( OrbitData[:,0] )    
    tmprad = OrbitData[:,2]
    OrbitX = tmprad * cos(tmplon) * cos(tmplat)
    OrbitY = tmprad * sin(tmplon) * cos(tmplat)
    OrbitZ = tmprad * sin(tmplat)
    # create the orbit scatter points for the plot
    OrbitScatter = dict( type = "scatter3d", mode = "markers", x = OrbitX,  y = OrbitY,  z = OrbitZ, showlegend = False, 
        marker = dict( 
            size=3, color = OrbitData[:,3], colorscale = ColorScale,  
            cmin = MinValue, cmax = MaxValue,
            colorbar=dict(thickness=20, len=0.75, ticklen=4, title=dict(text=ColorbarTitle,side="top"), xanchor="center", x=0) # x must be between [-2, 3]
        ),
    )
    return OrbitScatter




'''
Converts the ENU east,north,up coordinates to ECEF x,y,z coordinates
https://en.wikipedia.org/wiki/Geographic_coordinate_conversion#From_ENU_to_ECEF
'''
def ENU_to_ECEF( LatRadians, LonRadians, EastSpeed, NorthSpeed, UpSpeed ):
    #F = Lat * pi/180 # consider (90-Lat)*pi/180   
    #L = Lon * pi/180
    F = LatRadians
    L = LonRadians
    Ux = -sin(L)*EastSpeed + -sin(F)*cos(L)*NorthSpeed + cos(F)*cos(L)*UpSpeed
    Uy =  cos(L)*EastSpeed + -sin(F)*sin(L)*NorthSpeed + cos(F)*sin(L)*UpSpeed
    Uz =       0*EastSpeed +         cos(F)*NorthSpeed +        sin(F)*UpSpeed
    return Ux, Uy, Uz
    


'''
18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35
    lat = np.arange(   87.5,  -88.5, -LatStep ) 36
    lon = np.arange( -180.0,  180.0,  LonStep ) 72

FOR TESTING
Place		Lat	Lon
---------------------------
Zero Point	0		0			SurfaceData[18][36] = 50000 # Zero Point
Sydney 		-34		151			SurfaceData[22][66] = 50000 # Sydney
Santiago	-33		-70			SurfaceData[24][22] = 50000 # Santiago
Tokyo		35		140			SurfaceData[10][64] = 50000 # Tokyo
NY			41.25	-75			SurfaceData[9][21]  = 50000 # NY
Athens		38		23			SurfaceData[10][40] = 50000 # Athens
'''



        
