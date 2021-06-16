# local imports
import Data as D #from Data import *
from scicolorscales import *
import scicolorscales

# system imports
import math
import datetime
import numpy as np
import plotly
import chart_studio.plotly as py 
import plotly.graph_objects as go
import plotly.figure_factory as ff
from plotly.subplots import make_subplots

MyColors = ["#217ca3", "#e29930", "#919636", "#af1c1c", "#e7552c", "#1b4b5a", "#e4535e", "#aebd38", "#ffbb00", "#2c7873"]
VariableMinValue =  999999999
VariableMaxValue = -999999999

def plotDistribution( VariableName, Buckets, SuperTitle="", ColorMap=vik, OrbitBuckets=None ):        
    if VariableName == "Joule Heating":
        MultiplicationFactor = 10**8 
        new_units = "10^-8 W/m3"
        
    fig_log = make_subplots(rows=1, cols=len(D.KPsequence), shared_yaxes=True, horizontal_spacing=0.015)
    fig_lin = make_subplots(rows=1, cols=len(D.KPsequence), shared_yaxes=True, horizontal_spacing=0.015)
    includeInPlot( fig_log, fig_lin, 1,   VariableName, Buckets, SuperTitle, ColorMap)
    if OrbitBuckets is not None: includeInPlot( fig_log, fig_lin, 2,   VariableName, OrbitBuckets, SuperTitle, ColorMap)
        
    # update layout
    #fig_log.update_layout( annotations=BinAnnotations )
    #fig_log.update_layout( shapes=FigureShapes )
    fig_log.update_layout( title=SuperTitle + "<br>" + " (" + new_units + ")", width=4800/4, height=1800/4, legend_orientation="h", legend= {'itemsizing': 'constant'}) 
    fig_log.update_yaxes(title="Altitude (km)", row=1, col=1)
    #fig_lin.update_layout( annotations=BinAnnotations )
    #fig_lin.update_layout( shapes=FigureShapes )
    fig_lin.update_layout( title=SuperTitle + "<br>" + " (" + new_units + ")", width=4800/4, height=1800/4, legend_orientation="h", legend= {'itemsizing': 'constant'}) 
    fig_lin.update_yaxes(title="Altitude (km)", row=1, col=1)
    # increase font size
    fig_log.update_xaxes( tickfont=dict(size=52/4) )
    fig_log.update_yaxes( tickfont=dict(size=52/4) )
    fig_lin.update_xaxes( tickfont=dict(size=52/4) )
    fig_lin.update_yaxes( tickfont=dict(size=52/4) )
    
    # ======== plot log scale
    for i in range(0, len(D.KPsequence)):
        fig_log.update_yaxes(range=[D.ALT_min, D.ALT_max+D.ALT_distance_of_a_bucket], row=1, col=i+1)
        fig_log.update_xaxes(type="log", row=1, col=i+1 )
    plotly.offline.init_notebook_mode(connected=True)
    plotly.offline.iplot(fig_log)    
    
    # ======== plot linear scale
    for i in range(0, len(D.KPsequence)):
        fig_lin.update_yaxes(range=[D.ALT_min, D.ALT_max+D.ALT_distance_of_a_bucket], row=1, col=i+1)
        fig_lin.update_xaxes(range=[VariableMinValue, VariableMaxValue], row=1, col=i+1 )
    plotly.offline.init_notebook_mode(connected=True)
    plotly.offline.iplot(fig_lin)
        
        
        
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

def includeInPlot( fig_log, fig_lin, TraceNum,   VariableName, Buckets, SuperTitle="", ColorMap=vik):  
    global VariableMinValue, VariableMaxValue
    if VariableName == "Joule Heating":
        MultiplicationFactor = 10**8 
        new_units = "10^-8 W/m3"
    
    # calculate min max variable values
    for aKP in D.KPsequence:
        for anALT in D.ALTsequence:
            for aMagLat in D.MLATsequence:
                for aMLT in D.MLTsequence:
                    if Buckets[(aKP, anALT, aMagLat, aMLT, "Minimum")] < VariableMinValue: 
                        VariableMinValue = Buckets[(aKP, anALT, aMagLat, aMLT, "Minimum")]
                    if Buckets[(aKP, anALT, aMagLat, aMLT, "Maximum")] > VariableMaxValue: 
                        VariableMaxValue = Buckets[(aKP, anALT, aMagLat, aMLT, "Maximum")]
    print("Min variable value =", VariableMinValue, "  Max variable value =",  VariableMaxValue)
        
    if TraceNum == 1:
        LineType = "solid"
        LineFade = 0
    else:
        LineType = "dot"
        LineFade = 0.5
    
    #### Plot 
    for aKP in D.KPsequence:
        for anALT in D.ALTsequence:
            currentColor = MyColors[ D.KPsequence.index(aKP) ]
    
            localmin = Buckets[aKP, anALT, D.MLAT_min, D.MLT_min, "Minimum"]
            localmax = Buckets[aKP, anALT, D.MLAT_min, D.MLT_min, "Maximum"]   
            slots    = len(Buckets[aKP, anALT, D.MLAT_min, D.MLT_min, "Distribution"])
            bucket_starts = np.arange(localmin, localmax, (localmax-localmin)/slots ) 
            bucket_values = Buckets[aKP, anALT, D.MLAT_min, D.MLT_min, "Distribution"]
            
            # relocate points to fit into their correspondent altitude range
            if len(Buckets[aKP, anALT, D.MLAT_min, D.MLT_min, "Distribution"]) > 0:
                local_max_distr_counts = max(Buckets[aKP, anALT, D.MLAT_min, D.MLT_min, "Distribution"])
                for i in range(0, len(bucket_values)):
                    bucket_values[i] /= local_max_distr_counts
                    bucket_values[i] = anALT + bucket_values[i]*D.ALT_distance_of_a_bucket
    
            fig_log.add_trace( go.Scatter(x=bucket_starts, y=bucket_values, mode='lines', line=dict(color=currentColor,width=12/4,dash=LineType), opacity=1-LineFade, showlegend=False), row=1, col=D.KPsequence.index(aKP)+1)
            fig_lin.add_trace( go.Scatter(x=bucket_starts, y=bucket_values, mode='lines', line=dict(color=currentColor,width=12/4,dash=LineType), opacity=1-LineFade, showlegend=False), row=1, col=D.KPsequence.index(aKP)+1)
    
            # **** add visuals for the median line
            Percentile50 = Buckets[aKP, anALT, D.MLAT_min, D.MLT_min, "Percentile50"]
            fig_log.add_trace( go.Scatter(x=[Percentile50,Percentile50], y=[anALT,anALT+D.ALT_distance_of_a_bucket], mode='lines', line=dict(color=currentColor,width=12/4,dash=LineType), opacity=1-LineFade, showlegend=False), row=1, col=D.KPsequence.index(aKP)+1)
            fig_lin.add_trace( go.Scatter(x=[Percentile50,Percentile50], y=[anALT,anALT+D.ALT_distance_of_a_bucket], mode='lines', line=dict(color=currentColor,width=12/4,dash=LineType), opacity=1-LineFade, showlegend=False), row=1, col=D.KPsequence.index(aKP)+1)
            # **** add visuals for standard deviation
            Variance = Buckets[aKP, anALT, D.MLAT_min, D.MLT_min, "Variance"]
            fig_log.add_trace( go.Scatter(x=[Percentile50-(Variance)**(1/2)/2,Percentile50+(Variance)**(1/2)/2], y=[anALT+(D.ALT_distance_of_a_bucket)/2, anALT+(D.ALT_distance_of_a_bucket)/2], mode='lines', line=dict(color=currentColor,width=8/4,dash=LineType), opacity=1-LineFade, showlegend=False), row=1, col=D.KPsequence.index(aKP)+1)            
            fig_lin.add_trace( go.Scatter(x=[Percentile50-(Variance)**(1/2)/2,Percentile50+(Variance)**(1/2)/2], y=[anALT+(D.ALT_distance_of_a_bucket)/2, anALT+(D.ALT_distance_of_a_bucket)/2], mode='lines', line=dict(color=currentColor,width=8/4,dash=LineType), opacity=1-LineFade, showlegend=False), row=1, col=D.KPsequence.index(aKP)+1)            
    
    
    
    