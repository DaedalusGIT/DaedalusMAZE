import sys
sys.path.insert(1, '../../SourceCode/')
import DaedalusGlobals as DaedalusGlobals

import numpy as np
from datetime import datetime
import plotly
import plotly.graph_objs as go
import chart_studio.plotly as py 
import random
import math
import time

from scipy import signal
from matplotlib import pyplot as plt
import plotly.offline as offline
plotly.io.orca.config.executable = '/home/maze/anaconda3/bin/orca'
from scipy import fftpack

from numba import vectorize
from numba import guvectorize
from numba import jit
from numba import cuda



R_Earth = 6371 # the Earth radius in km
LatMin =  50   # the min latitude for which the gaussian functions will be created
LatMax =  90   # the max latitude for which the gaussian functions will be created
LonMin = -180  # the min longitude for which the gaussian functions will be created
LonMax =  180  # the max longitude for which the gaussian functions will be created
AmpMax = 200   # the maximum Amplitude for which the gaussian functions will be created
ScaleLat = 100 # km
FuncParams = list() # a list of lists. each element is a set of attributes for a Gaussian function which represents a spike

# Construct random functions of the type:
#     f(LAT,LON) = Amplitude * math.e ^ ( - ((LAT+slideLat)^2) / (2*sigmaLat^2)  -  ((LON+slideLon)^2) / (2*sigmaLon^2) )
# and stores their attributes in FuncParams.
# Each Gaussian function represents a spike in 2d world of lat-lon coordinates.
# All spikes together construct the Sub-Grid-Variability
def initGaussianFunctions():
    global FuncParams, col0, col1, col2, col3, col4
    for i in range( 0, 10000 ): # 600 or 90000 for world
        Amplitude = random.random() * 2 *AmpMax - AmpMax # spike height
        sigmaLat  = 57.295 * np.random.normal( 0.5, 0.18 ) *  ScaleLat / (R_Earth * 2.3548) # could be gaussian np.random.normal( 0.5, 0.18 )
        slideLat  = random.random() * 180 - 90   # spike position at Latitudes - they will be contained by LatitudeModulationNorth
        sigmaLon  = 2 * sigmaLat #*   np.rad2deg(math.sin( np.deg2rad(sigmaLat-90) ))
        slideLon  = random.random() * (LonMax-LonMin) - (LonMax-LonMin)/2  # spike position at Longitudes
        FuncParams.append( [Amplitude, sigmaLat, slideLat, sigmaLon, slideLon] )

    for i in range( 0, 10000 ): # 3000
        Amplitude = (random.random() * 2 *AmpMax - AmpMax )# spike height
        sigmaLat  = 57.295 * np.random.normal( 0.5, 0.18 ) *  (ScaleLat/10) / (R_Earth * 2.3548) # could be gaussian np.random.normal( 0.5, 0.18 )
        slideLat  = random.random() * 180 - 90   # spike position at Latitudes - they will be contained by LatitudeModulationNorth
        sigmaLon  = 2 * sigmaLat # *    np.rad2deg(math.sin( np.deg2rad(sigmaLat-90) ))
        slideLon  = random.random() * (LonMax-LonMin) - (LonMax-LonMin)/2  # spike position at Longitudes
        FuncParams.append( [Amplitude, sigmaLat, slideLat, sigmaLon, slideLon] )

    for i in range( 0, 10000 ): #1000
        Amplitude = (random.random() * 2 *AmpMax - AmpMax) # spike height
        sigmaLat  = 57.295 * np.random.normal( 0.5, 0.18 ) *  (ScaleLat/100) / (R_Earth * 2.3548) # could be gaussian np.random.normal( 0.5, 0.18 )
        slideLat  = random.random() * 180 - 90   # spike position at Latitudes - they will be contained by LatitudeModulationNorth
        sigmaLon  = 2 * sigmaLat #*  np.rad2deg(math.sin( np.deg2rad(sigmaLat-90) ))
        slideLon  = random.random() * (LonMax-LonMin) - (LonMax-LonMin)/2  # spike position at Longitudes
        FuncParams.append( [Amplitude, sigmaLat, slideLat, sigmaLon, slideLon] )
    # convert to numpy arrays    
    np_FuncParams = np.array(FuncParams, dtype="float32")
    col0 = np.ascontiguousarray(np_FuncParams[:,0], dtype="float32")
    col1 = np.ascontiguousarray(np_FuncParams[:,1], dtype="float32")
    col2 = np.ascontiguousarray(np_FuncParams[:,2], dtype="float32")
    col3 = np.ascontiguousarray(np_FuncParams[:,3], dtype="float32")
    col4 = np.ascontiguousarray(np_FuncParams[:,4], dtype="float32")
    # copy to cuda graphics card
    col0 = cuda.to_device( col0 )
    col1 = cuda.to_device( col1 )
    col2 = cuda.to_device( col2 )
    col3 = cuda.to_device( col3 )
    col4 = cuda.to_device( col4 )
    print( len(FuncParams), "Gaussian functions initialized.")

# NOT USED
@cuda.jit
def GaussianJIT( LON, LAT,    Amplitude, sigmaLat, slideLat, sigmaLon, slideLon,   out ):
    num_of_blocks  = cuda.gridDim.x   # number of blocks in the grid
    num_of_threads = cuda.blockDim.x  # number of threads per block
    BlockIDX  = cuda.blockIdx.x  # this is the unique block ID within the 1D grid
    ThreadIDX = cuda.threadIdx.x # this is the unique thread ID within a 1D block
    ThreadStep = int(len(Amplitude)/num_of_threads)+1
    result = 0.0
    LatMin =  50
    LatMax =  90    
    sigma_mod =  (LatMax - LatMin)/2.3548
    #if abs(LAT) < abs(LatMin)-5: return 0  # <<<<
    
    for idx in range( ThreadIDX*ThreadStep,  ThreadIDX*ThreadStep + ThreadStep ):
        # main function - for North Altitudes
        LatitudeModulationNorth =  -((LAT-LatMax)**2) / (2*sigma_mod**2) 
        F1 = ((LAT-slideLat[idx])**2)/(2*sigmaLat[idx]**2)
        F2 = (2*sigmaLon[idx]**2)
        result += Amplitude[idx] * math.e ** ( - F1 - ((LON-slideLon[idx])**2)/F2  + LatitudeModulationNorth )
    
    # for symmetry - for South Altitudes
    LatitudeModulationSouth =  -((LAT+LatMax)**2) / (2*sigma_mod**2) 
    result += Amplitude[idx] * math.e ** ( - F1 - ((LON-slideLon[idx])**2)/F2  + LatitudeModulationSouth )
            
    # respect the boundary conditions: spikes at the edge of Longitudes have to be applied to the other edge as well so that the grid becomes a cylinder
    #if LON > 167:
    #    result += Amplitude[idx] * math.e ** ( - F1 - ((LON-360-slideLon[idx])**2)/F2  + LatitudeModulationNorth )
    #    result += Amplitude[idx] * math.e ** ( - F1 - ((LON-360-slideLon[idx])**2)/F2  + LatitudeModulationSouth )
    #elif LON < -167:
    #    result += Amplitude[idx] * math.e ** ( - F1 - ((LON+360-slideLon[idx])**2)/F2  + LatitudeModulationNorth )
    #    result += Amplitude[idx] * math.e ** ( - F1 - ((LON+360-slideLon[idx])**2)/F2  + LatitudeModulationSouth )
    ##
    out[idx] = result
'''
### WITH CUDA-JIT
startSecs = time.time()
out = np.zeros(  [len(col0)]  )
blocks_per_grid = (1)
threads_per_block = (1024)
for i in range( 0, len(Lons) ):
    for j in range( 0, len(Lats) ):
        LON = Lons[i]
        LAT = Lats[j]
        if LAT < LatMin -15: continue # <<<<        
        GaussianJIT[blocks_per_grid, threads_per_block](LON, LAT, col0, col1, col2, col3, col4, out)
        Data[i][j] = np.sum(out)
finishSecs = time.time()
print( "TEST Sub-Grid-Variability for world finished in" , finishSecs-startSecs, "sec")
'''    

    
    
    
@vectorize(['float32(float32, float32, float32, float32, float32, float32, float32)'], target='cuda')
def Gaussian( LON, LAT, Amplitude, sigmaLat, slideLat, sigmaLon, slideLon ):
    result = 0.0
    sigma_mod =  (LatMax - LatMin)/2.3548
    #if abs(LAT) < abs(LatMin)-5: return 0  # <<<<
    
    # main function - for North Altitudes
    LatitudeModulationNorth =  -((LAT-LatMax)**2) / (2*sigma_mod**2) 
    F1 = ((LAT-slideLat)**2)/(2*sigmaLat**2)
    F2 = (2*sigmaLon**2)
    result += Amplitude * math.e ** ( - F1 - ((LON-slideLon)**2)/F2  + LatitudeModulationNorth )
    
    # for symmetry - for South Altitudes
    LatitudeModulationSouth =  -((LAT+LatMax)**2) / (2*sigma_mod**2) 
    result += Amplitude * math.e ** ( - F1 - ((LON-slideLon)**2)/F2  + LatitudeModulationSouth )
            
    # respect the boundary conditions: spikes at the edge of Longitudes have to be applied to the other edge as well so that the grid becomes a cylinder
    if LON > 160:
        result += Amplitude * math.e ** ( - F1 - ((LON-360-slideLon)**2)/F2  + LatitudeModulationNorth )
        result += Amplitude * math.e ** ( - F1 - ((LON-360-slideLon)**2)/F2  + LatitudeModulationSouth )
    elif LON < -160:
        result += Amplitude * math.e ** ( - F1 - ((LON+360-slideLon)**2)/F2  + LatitudeModulationNorth )
        result += Amplitude * math.e ** ( - F1 - ((LON+360-slideLon)**2)/F2  + LatitudeModulationSouth )
    ##
    return result

# Calculate Sub-Grid-Variability for world #################################################################
def SGVpreview_forWorlds():
    if len(FuncParams) == 0: initGaussianFunctions()
    startSecs = time.time()
    # init
    Lons = list( range( -180, 180, 2 ) )
    Lats = list( range(  -90,  90, 2 ) )
    Zaxis = np.zeros( (len(Lons), len(Lats)) ).tolist()
    Data  = np.zeros( (len(Lons), len(Lats)) ).tolist()
    Lons = np.array(Lons, dtype="float32")
    Lats = np.array(Lats, dtype="float32")
    # calc
    for i in range( 0, len(Lons) ):
        for j in range( 0, len(Lats) ):
            LON = Lons[i]
            LAT = Lats[j]
            AllResults = Gaussian( LON, LAT, col0, col1, col2, col3, col4 )
            Data[i][j] = np.sum(AllResults)
    finishSecs = time.time()
    print( "Sub-Grid-Variability for world finished in" , finishSecs-startSecs, "sec")
    # PLOT 
    # define the layout of the plot
    theLayout = dict( title="Sub Grid Variability (Lat-Lon World)", width=1200, height=800, scene = dict( xaxis = dict(title="Latitude",dtick=10), yaxis = dict(title="Longitude",dtick=20), zaxis = dict( title="" ), ) )
    # plot 3d world    
    surface=go.Surface(x=Lats, y=Lons, z=Data,   colorscale = "jet", surfacecolor = Data, colorbar=dict(thickness=20, len=0.75, ticklen=4, title="mV/m"), )
    fig = dict( data=[surface], layout=theLayout )
    plotly.offline.init_notebook_mode(connected=True)
    plotly.offline.iplot(fig)
    # plot projection on XY plane
    surface=go.Surface(x=Lats, y=Lons, z=Zaxis,   colorscale = "jet", surfacecolor = Data, colorbar=dict(thickness=20, len=0.75, ticklen=4, title="mV/m"), )
    fig = dict( data=[surface], layout=theLayout )
    plotly.offline.init_notebook_mode(connected=True)
    plotly.offline.iplot(fig)
    # plot projections on XZ and YZ planes
    surface=go.Surface(x=Lats, y=np.full(len(Lons),-180).tolist(), z=Data,   colorscale = "jet", surfacecolor = Data, colorbar=dict(thickness=20, len=0.75, ticklen=4, title="mV/m"), )
    surface2=go.Surface(x=np.full(len(Lats),-90).tolist(), y=Lons, z=Data,   colorscale = "jet", surfacecolor = Data, colorbar=dict(thickness=20, len=0.75, ticklen=4, title="mV/m"), )
    fig = dict( data=[surface,surface2], layout=theLayout )
    plotly.offline.init_notebook_mode(connected=True)
    plotly.offline.iplot(fig)
       


    
    



# Calculate Sub-Grid-Variability for a certain Longitude #################################################################
def SGVpreview_forCertainLongitude( aLongitude ):
    if len(FuncParams) == 0: initGaussianFunctions()
    startSecs = time.time()
    SectionData = list()
    LON = aLongitude
    myRange = np.arange( 50, 90, 0.01 )
    for LAT in myRange:
        AllResults = Gaussian( LON, LAT, col0, col1, col2, col3, col4 )
        value = np.sum(AllResults)
        SectionData.append( value )
    finishSecs = time.time()    
    print( "Sub-Grid-Variability for LON =",aLongitude,"finished in" , finishSecs-startSecs, "sec.", "Points:", len(SectionData) )
    # plot the Sub-Grid-Variability for a certain Longitude as line-plot
    fig = go.Figure(data=go.Scatter(x=list(myRange), y=SectionData)) #fig = go.Figure(data=go.Scatter( mode='markers+lines', x=list(myRange), y=SectionData ,marker=dict( colorscale = "jet", color=SectionData, size=2 )) )
    fig.update_layout( title="values along Lon=0", xaxis_title="Latitudes", yaxis_title="mV/m" )
    plotly.offline.init_notebook_mode(connected=True)
    plotly.offline.iplot(fig) 
    # plot frequences
    
    SectionData = np.array(SectionData, dtype=np.float32)
    c = 0.001
    fs=1/c
    X = fftpack.fft(SectionData)
    freqs = fftpack.fftfreq(len(SectionData)) * fs
    fig = go.Figure(data=go.Scatter(x=list(freqs), y=list(10*np.log(np.abs(X)))))
    fig.update_layout( title="values along Lon=0", xaxis_title="Frequencies", yaxis_title="??" )
    plotly.offline.init_notebook_mode(connected=True)
    plotly.offline.iplot(fig)         
    
    
    
'''
def GG( LAT, LON ):
    result = 0.0
    LatMin =  50
    LatMax =  90
    for G in FuncParams:
        slideLat  = G[2]
        slideLon  = G[4]
        if abs(LAT-slideLat) > 6  and  abs(LON-slideLon) > 6: continue  # <<<<
        Amplitude = G[0]
        sigmaLat  = G[1]
        sigmaLon  = G[3]
        # main function - for North Altitudes
        sigma_mod =  (LatMax - LatMin)/2.3548
        LatitudeModulationNorth =  -((LAT-LatMax)**2) / (2*sigma_mod**2) 
        F1 = ((LAT-slideLat)**2)/(2*sigmaLat**2)
        F2 = (2*sigmaLon**2)
        result += Amplitude * math.e ** ( - F1 - ((LON-slideLon)**2)/F2  + LatitudeModulationNorth )
        # for symmetry - for South Altitudes
        LatitudeModulationSouth =  -((LAT+LatMax)**2) / (2*sigma_mod**2) 
        result += Amplitude * math.e ** ( - F1 - ((LON-slideLon)**2)/F2  + LatitudeModulationSouth )
        #
        return result
startSecs = time.time()
SectionData = list()
LON = 0          
myRange = np.arange( 50, 90, 0.01 )
for LAT in myRange:
    value = Gaussian( LON, LAT, col0, col1, col2, col3, col4 )
    SectionData.append( value )
finishSecs = time.time()    
print( "Sub-Grid-Variability (ver2) for LON=0 finished in" , finishSecs-startSecs, "sec.", "Points:", len(SectionData) )
# plot the Sub-Grid-Variability for a certain Longitude as line-plot
fig = go.Figure(data=go.Scatter(x=list(myRange), y=SectionData))
fig.update_layout( title="values along Lon=0", xaxis_title="Latitudes", yaxis_title="mV/m" )
plotly.offline.init_notebook_mode(connected=True)
plotly.offline.iplot(fig)           
'''

  

# Creates an image for every longitude and stores it in a folder in order to become an animation
# You can make an animated gif using $> convert -delay 10 -loop 0 *.png myimage.gif
# You can make an mp4 movie using    $> ffmpeg -framerate 10 -i lon%03d.png -c:v libx264 -r 10 -pix_fmt yuv420p out.mp4
# You can batch  convert PNG tp JPG using  $> mogrify -format jpg *.png
# You can convert PNG tp JPG using         $> convert image.png image.jpg
def createSliceImages_forEveryLongitude():
    print("Start creating images for all Longitudes")
    startSecs = time.time()
    myRange = np.arange( 50, 90, 0.01 )
    for LON in range(-180, 180, 1):
        if LON%20==0: print( "working Longitude", LON )
        SectionData = list()
        for LAT in myRange:
            AllResults = Gaussian( LON, LAT, col0, col1, col2, col3, col4 )
            value = np.sum(AllResults)
            SectionData.append( value )
        fig = go.Figure(data=go.Scatter(x=list(myRange), y=SectionData))
        fig.update_layout( title="Sub-Grid_Variability along Longitude "+"{:03.0f}".format(LON), xaxis_title="Latitudes", yaxis_title="mV/m" )
        fig.update_yaxes(range=[-180, 180])
        img_filename = DaedalusGlobals.CoverageResults_Files_Path + "SGVimages/lon"+"{:03.0f}".format(LON+181)+".png"
        try:
            fig.write_image( img_filename )    
        except:
            print("Plotly-Orca failed to save", img_filename)
    finishSecs = time.time()    
    print( "Sub-Grid-Variability for each Longitude finished in" , finishSecs-startSecs, "sec.")

    
    
def calculateSGV_forSinglePoint( Longitude, Latitude, MaxAmplitude=200, Scale=100 ):
    global AmpMax, ScaleLat
    if len(FuncParams)==0  or AmpMax!=MaxAmplitude or ScaleLat!=Scale:
        AmpMax = MaxAmplitude
        ScaleLat = Scale
        initGaussianFunctions()
    return np.sum ( Gaussian( Longitude, Latitude, col0, col1, col2, col3, col4 ) )

# EXAMPLES OF HOW TO USE THIS MODULE:
#print( "Test: SGV for single point (11.97,17.28) is", calculateSGV_forSinglePoint(11.97,17.28) )
#SGVpreview_forWorlds()
#SGVpreview_forCertainLongitude(0)
#createSliceImages_forEveryLongitude():
