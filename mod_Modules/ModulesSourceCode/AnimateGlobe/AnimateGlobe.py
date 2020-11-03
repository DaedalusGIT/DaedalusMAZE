'''
The source code of this file is based on 
    https://plot.ly/~empet/14813/heatmap-plot-on-a-spherical-map/#/
    https://plot.ly/~rreusser/62/animating-the-contour-plot-of-a-wave-ux/#/
    https://community.plot.ly/t/3d-surface-animation-using-two-types-data/13385/9
The main function of this file is AnimateGlobe. The others are ancillary.

Notes
-----
We can add these plot types to Plotables:
    'area', 'bar', 'barpolar', 'box', 'candlestick', 'carpet', 'choropleth', 'cone', 'contour', 'contourcarpet', 
    'funnel', 'funnelarea', 'heatmap', 'heatmapgl', 'histogram', 'histogram2d', 'histogram2dcontour', 'isosurface', 
    'mesh3d', 'ohlc', 'parcats', 'parcoords', 'pie', 'pointcloud', 'sankey', 'scatter', 'scatter3d', 'scattercarpet', 
    'scattergeo','scattergl', 'scattermapbox', 'scatterpolar','scatterpolargl', 'scatterternary', 'splom', 'streamtube', 
    'sunburst', 'surface', 'table', 'violin', 'volume', 'waterfall'

'''

from numpy import pi, sin, cos
import time
import plotly
import plotly.graph_objs as go
# import plotly.plotly as py
import chart_studio.plotly as py
import csv
import numpy as np           
from scipy.io import netcdf  
from mpl_toolkits.basemap import Basemap
import warnings

EarthRadius = 6378.137 # km

#convert degrees to radians
def degree2radians(degree):    
    return degree*pi/180

# maps the points of coords (lon, lat) to points onto the sphere of radius radius
def mapping_map_to_sphere(lon, lat, radius=1):  
    lon=np.array(lon, dtype=np.float64)
    lat=np.array(lat, dtype=np.float64)
    lon=degree2radians(lon)
    lat=degree2radians(lat)
    xs=radius*cos(lon)*cos(lat)
    ys=radius*sin(lon)*cos(lat)
    zs=radius*sin(lat)
    return xs, ys, zs

# Make shortcut to Basemap object, not specifying projection type
m = Basemap() 

# Functions converting coastline/country polygons to lon/lat traces
''' 
    pos arg 1. (poly_paths): paths to polygons
    pos arg 2. (N_poly): number of polygon to convert
'''
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

# Function generating coastline lon/lat 
def get_coastline_traces():
    poly_paths = m.drawcoastlines().get_paths() # coastline polygon paths
    N_poly = 91  # use only the 91st biggest coastlines (i.e. no rivers)
    cc_lons, cc_lats= polygons_to_traces(poly_paths, N_poly)
    return cc_lons, cc_lats

# Function generating country lon/lat 
def get_country_traces():
    poly_paths = m.drawcountries().get_paths() # country polygon paths
    N_poly = len(poly_paths)  # use all countries
    country_lons, country_lats= polygons_to_traces(poly_paths, N_poly)
    return country_lons, country_lats

###############

'''
Creates a 3D plot of an earth globe, a sphere surface and a satellite orbit. The surface and the orbit are colored according to the data values in the CSV files.
  DataCSVfilename: format: Time,Lat,Lon,Alt,value. Contains the data for the sphere surface. If empty then no surface will be plotted.
  OrbitDataCSVfilename: format: Time,Lat,Lon,Alt,value. Contains the data for the orbit. If empty then no orbit will be plotted.
  PlotTitle: the title of the plot. It will be displayed at the top of the globe.
  ColorbarTitle: the title of the colorbar beside the globe. Usually contains measuremnt units
  ColorscaleName: valid values are: ‘Blackbody’, ‘Bluered’, ‘Blues’, ‘Earth’, ‘Electric’, ‘Greens’, ‘Greys’, ‘Hot’, ‘Jet’, ‘Picnic’, ‘Portland’, ‘Rainbow’, ‘RdBu’, ‘Reds’, ‘Viridis’, ‘YlGnBu’, ‘YlOrRd’
                  In case an empty string is passed then a default HeatMap colorscale will be applied.
                  In case None is passed then all points will be black irrespective of value.
  RETURNS: a string containing information about the Data
'''
def AnimateGlobe( DataCSVfilenames, OrbitDataCSVfilename, PlotTitle, ColorbarTitle, ColorscaleName ):
    result = ""
    Plotables = list()
    AllSurfaces = list()
    startSecs = time.time()
    
    # construct FrameNames to display at the animation slider
    FrameNames = list()
    for s in DataCSVfilenames:
        FrameNames.append( s[s.rfind('/')+1 : len(s)-4] )
    
    # create coloscales for ploting
    global colorscale
    colorscaleWhite=[[0.0, '#ffffff'], [1.0, '#ffffff']]
    colorscaleBlack=[[0.0, '#000000'], [1.0, '#000000']]
    if ColorscaleName is None:
        colorscale = colorscaleBlack
    elif len( ColorscaleName ) > 0: 
        colorscale = ColorscaleName
    else:
        colorscale=[[0.0, '#313695'], [0.07692307692307693, '#3a67af'], [0.15384615384615385, '#5994c5'], [0.23076923076923078, '#84bbd8'], [0.3076923076923077, '#afdbea'], [0.38461538461538464, '#d8eff5'], [0.46153846153846156, '#d6ffe1'], [0.5384615384615384, '#fef4ac'], [0.6153846153846154, '#fed987'], [0.6923076923076923, '#fdb264'], [0.7692307692307693, '#f78249'], [0.8461538461538461, '#e75435'], [0.9230769230769231, '#cc2727'], [1.0, '#a50026']]

    #### create the earth surface
    # construct the values for longitude and latitude
    lon = np.arange( -180.0,  180.0,  5 )
    lat = np.arange(   87.5,  -88.5, -5 )
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
    
    # define the general layout of the plot
    theLayout = dict(
        title = PlotTitle,
        scene = dict(
            xaxis = dict( zeroline=False, autorange=False, range=[-7500, 7500] ),
            yaxis = dict( zeroline=False, autorange=False, range=[-7500, 7500] ), 
            zaxis = dict( zeroline=False, autorange=False, range=[-7500, 7500] ),
            aspectratio=dict(x=1, y=1,z=1), camera=dict(eye=dict(x=1.20, y=1.20, z=1.20))
        )
    )
    

    for f in DataCSVfilenames:
        #### parse DataCSVfilename and find the Altitude of the surface described by the csv file
        CSVreader = csv.reader( open( f ) )
        next( CSVreader ) # ignore the csv header
        Altitude = float(   next( CSVreader )[3]   )
        
        #### read DataCSVfilename into a 2D Array
        Data =  CSV_to_Datagrid( f )         

        
        #### Calculate min / max values of all plotable data 
        # for the spere surface around earth
        if len(f) > 0:
            DataMin = np.amin(Data)
            DataMax = np.amax(Data)
        # for the orbit
        if len(OrbitDataCSVfilename) > 0:
            orbitCSVreader = csv.reader( open( OrbitDataCSVfilename ) )
            next( orbitCSVreader ) # ignore the csv header
            OrbitData =  [r[1:5] for r in orbitCSVreader if len(r[4])>0] # read the rest CSV file at once into a 2D list. Ignore first column which is the time and take only the first value in case there are more. Ignore lines without value, as well.
            OrbitData = np.array( OrbitData, dtype=np.float64 ) # convert to 2D numpy array
            OrbitMin = np.amin( OrbitData[:,3] )
            OrbitMax = np.amax( OrbitData[:,3] )
        # Verdict
        if len(f) > 0  and  len(OrbitDataCSVfilename) > 0:
            GeneralMin = min( DataMin, OrbitMin)
            GeneralMax = max( DataMax, OrbitMax)
        elif len(f) > 0:
            GeneralMin = DataMin
            GeneralMax = DataMax
        elif len(OrbitDataCSVfilename) > 0:
            GeneralMin = OrbitMin
            GeneralMax = OrbitMax    

        AllSurfaces.append( CreatePlotable_Surface( Data, Altitude, GeneralMin, GeneralMax, colorscale, ColorbarTitle ) )        
        
    
    # add all visual elements in a list and assign them to a figure
    '''
    Plotables.append( AllSurfaces[0] )
    Plotables.append( EarthSurface )
    Plotables.append( Boundaries )
    if len(OrbitDataCSVfilename) > 0:
        Plotables.append( CreatePlotable_Orbit( OrbitDataCSVfilename, GeneralMin, GeneralMax, colorscale, "" ) )
    fig = go.Figure( data=Plotables, layout=theLayout  ) #fig = dict( data=Plotables, layout=theLayout )
    '''
    
    
    # make figure
    fig_dict = {
        "data": [EarthSurface,Boundaries, AllSurfaces[0]], "layout": theLayout, "frames": []
        #"data": [ AllSurfaces[0], AllSurfaces[1], AllSurfaces[2]], "layout": theLayout, "frames": []
    }

    # sliders
    fig_dict["layout"]["sliders"] = {
        "args": [ "transition", {"duration": 1000} ],
        "plotlycommand": "animate", "values": FrameNames, "visible": True
    }
    
    fig_dict["layout"]["updatemenus"] = [{        
        'buttons': [
            {
                'args': [None, {'frame': {'duration': 1000, 'redraw': True},'fromcurrent': True, 'transition': {'duration': 1000}}],
                'label': 'Play', 'method': 'animate'
            },
            {
                'args': [[None], {'frame': {'duration': 0, 'redraw': True}, 'mode': 'immediate', 'transition': {'duration': 0}}],
                'label': 'Pause', 'method': 'animate'
            }
        ],
        'direction': 'left', 'pad': {'r': 10, 't': 87}, 'showactive': False, 'type': 'buttons',
        'x': 0, 'xanchor': 'right',  'y': 0, 'yanchor': 'top'
    } ]
    
    # create frames
    for i in range( len(AllSurfaces) ):
        frame = { "data": [AllSurfaces[i], EarthSurface, Boundaries ], "name": FrameNames[i] }
        fig_dict["frames"].append( frame )
        

    

    sliders_dict = {
        "active": 0, "yanchor": "top", "xanchor": "left",
        "currentvalue": { "font": {"size": 20}, "prefix": "frame: ", "visible": True, "xanchor": "right" },
        "transition": {"duration": 1000},
        "pad": {"b": 10, "t": 50},  "len": 0.9, "x": 0.1, "y": 0,  "steps": []
    }
    for i in range( len(AllSurfaces) ):
        slider_step = {"method": 'animate',
                       "args": [
                        [FrameNames[i]],
                        {"frame": {"duration": 1000, "redraw": True}}
                      ],
                      "label": FrameNames[i],
                      }
        sliders_dict["steps"].append(slider_step)
    
    fig_dict["layout"]["sliders"] = [sliders_dict]
    
    
    ###
    fig = go.Figure(fig_dict)

    
    # plot all
    plotly.offline.init_notebook_mode(connected=True)
    plotly.offline.iplot(fig)    
   
    
    finishSecs = time.time()
    # be verbose
    #result = result + "Lon shape: " + str(lon.shape) + "  Lat shape: " + str(lat.shape) + "  Data shape: " + str(Data.shape) + "\n"
    #if len(DataCSVfilename) > 0: result = result + "Data Min: " + str(DataMin) + "  Data Max: " + str(DataMax) + "\n"
    if len(OrbitDataCSVfilename) > 0: result = result + "Orbit Min: " + str(OrbitMin) + "  Orbit Max: " + str(OrbitMax) + "\n"
    result = result + "General Min: " + str(GeneralMin) + "  General Max: " + str(GeneralMax) + "\n"
    if len(OrbitDataCSVfilename) > 0: result = result + "Number of Orbit positions: "  + str(len(OrbitData)) + "\n"    
    result = result + "Plot creation duration: " + str(finishSecs-startSecs) + " seconds.\n"
    
    print("ready")

    #for i in range( len(AllSurfaces) ):
    #    print( i )
    #    fig["data"][0]["surfacecolor"] = AllSurfaces[i]["surfacecolor"]
    #    time.sleep(4)

    
    return result;


###################################################################################################
###################################################################################################
###################################################################################################


def CSV_to_Datagrid( InputFilenameCSV ):    
    ## read the csv file. format: (time lat lon alt value)
    CSVreader = csv.reader( open( InputFilenameCSV ) )
    # ignore the csv header
    next( CSVreader ) 
    # read the rest file into a 2D list. Ignore 1st column (time) and take only the first value. Ignore lines without value.
    RawData =  [r[1:5] for r in CSVreader if len(r[4])>0]
    # convert to 2D numpy array
    RawData = np.array( RawData, dtype=np.float64 ) 
    
    # TODO: Look at the data and find out the step for lat and lon values
    LatStep  = 5
    LonStep  = 5
    
    Lat = np.arange(   87.5,  -88.5,  -1*LatStep )
    Lon = np.arange( -180.0,  180.0,     LonStep )
    DataGrid = np.zeros( (len(Lat), len(Lon)), dtype=np.float64 )
    
    for i in range(0, len(Lat)):
        for j in range(0, len(Lon)):
            for x in range(0, len(RawData)):
                if( RawData[x,0] == Lat[i]  and  RawData[x,1] == Lon[j] ):
                    DataGrid[i, j] = RawData[x, 3]
    return DataGrid

###################################################################################################

def CreatePlotable_Surface( Data, Altitude, GeneralMin, GeneralMax , ColorScale, ColorbarTitle):
    # assign/read data which will be plotted
    #Data = np.array( list(csv.reader(open(  DataCSVfilename  ))), dtype=np.float64)  # make a 2D array out of the rest csv file

    # construct the values for longitude and latitude
    lon = np.arange( -180.0,  180.0,  int( 360 / len(Data[0]) ) )
    lat = np.arange(   87.5,  -88.5,  int(-180 / len(Data   ) ) )
    
    #### SHIFTING
    # Shift 'lon' from [0,360] to [-180,180]
    tmp_lon = np.array([lon[n]-360 if l>=180 else lon[n] for n,l in enumerate(lon)])  # => [0,180]U[-180,2.5]
    i_east, = np.where(tmp_lon>=0)  # indices of east lon
    i_west, = np.where(tmp_lon<0)   # indices of west lon
    lon = np.hstack((tmp_lon[i_west], tmp_lon[i_east]))  # stack the 2 halves
    # Correspondingly, shift the Data array
    Data_ground = np.array(Data)
    Data = np.hstack((Data_ground[:,i_west], Data_ground[:,i_east]))
    
    # cut out a pizza slice
    #for i in range (0,16):
    #    for j in range (0,16):
    #        Data[i,j] = float('nan')
    
    # To ensure color continuity we extend the lon list with [180] (its last value was lon[-1]=177.5). In this way we can identify lon=-180 with lon=180.
    # We do the same with latitudes. We extend both directions with [90] and [-90] in order to have values for the poles.
    clons = np.array( lon.tolist() + [180]        , dtype=np.float64)
    clats = np.array( [90] + lat.tolist() + [-90] , dtype=np.float64)
    clons, clats=np.meshgrid(clons, clats)
    # Map the meshgrids clons, clats onto the sphere
    XS, YS, ZS = mapping_map_to_sphere(clons, clats, radius=EarthRadius+Altitude)
    # for color continuity
    nrows, ncolumns=clons.shape
    DATA = np.zeros(clons.shape, dtype=np.float64)
    DATA[1:nrows-1, :ncolumns-1] = np.copy(np.array(Data,  dtype=np.float64)) # ignore the extended values
    DATA[1:nrows-1,  ncolumns-1] = np.copy(Data[:, 0]) # ignore the extended values

    # Create a sphere above earth which will be colored in accordance with the data
    DataSphere=dict(type='surface', x = XS,  y = YS,  z = ZS,
        colorscale = ColorScale, surfacecolor = DATA, opacity = 0.90,
        cmin = GeneralMin, cmax = GeneralMax,
        colorbar=dict(thickness=20, len=0.75, ticklen=4, title=ColorbarTitle),
    )
    if len(ColorbarTitle) > 0 : DataSphere["showscale"] = False

    return DataSphere
    

###########################################################################################################


def CreatePlotable_Orbit( OrbitDataCSVfilename, GeneralMin, GeneralMax, ColorScale, ColorbarTitle ):
    ## read the csv file. format: (time lat lon alt value)
    orbitCSVreader = csv.reader( open( OrbitDataCSVfilename ) )
    # ignore the csv header
    next( orbitCSVreader ) 
    # read the rest file into a 2D list. Ignore 1st column (time) and take only the first value. Ignore lines without value.
    OrbitData =  [r[1:5] for r in orbitCSVreader if len(r[4])>0]
    # convert to numpy array
    OrbitData = np.array( OrbitData, dtype=np.float64 ) # convert to 2D numpy array
    
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
        marker = dict( size=3, color = OrbitData[:,3], colorscale = ColorScale,  cmin = GeneralMin, cmax = GeneralMax )
    )
    
    #if len(ColorbarTitle) > 0 :
    #else:
    
    return OrbitScatter


# Example of Usage:
'''
AllSurfaces = [
                "./files/IRI100.csv", 
                "./files/IRI150.csv",     
                "./files/IRI200.csv",     
                "./files/IRI250.csv", 
                "./files/IRI300.csv",
                "./files/IRI350.csv",     
                "./files/IRI400.csv",     
                "./files/IRI450.csv",     
                "./files/IRI500.csv",     
              ]
s = AnimateGlobe.AnimateGlobe( AllSurfaces, "", "Electron Density currents", "cm^-3", "Jet" )
print(s)
'''