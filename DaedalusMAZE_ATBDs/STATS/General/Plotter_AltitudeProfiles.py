import Data as D

import numpy as np
import plotly
import chart_studio.plotly as py 
import plotly.graph_objects as go
import plotly.figure_factory as ff
from plotly.subplots import make_subplots

def plotAltProfiles( VariableName, Buckets, SuperTitle="" ):
    if VariableName == "Joule Heating":
        x_axes_range=[0, 4]
        MultiplicationFactor = 10**8 
        new_units = "10^-8 W/m3"
    else:
        x_axes_range=[0, 0.1]
        MultiplicationFactor = 10**3 
        new_units = "mS/m"
    
    # alter visibleALTsequence so that the point is displayed in the middle of the sub-bin
    visibleALTsequence = D.ALTsequence.copy()
    for i in range(1, len(visibleALTsequence)-1):
        visibleALTsequence[i] += D.ALT_distance_of_a_bucket/2
    visibleALTsequence[0] = D.ALTsequence[0]
    visibleALTsequence[-1] = D.ALTsequence[-1] + D.ALT_distance_of_a_bucket
    
    # plot
    Color10 = '#c4dfe6'
    Color25 = '#a1d6e2'
    Color50 = '#1995ad'
    Color75 = '#a1d6e2'
    Color90 = '#c4dfe6'
    
    # construct the column MLT titles #("0-3", "3-6", "6-9", "9-12", "12-15", "15-18", "18-21", "21-24")
    ColumnTitles = list()
    
    for i in range(0, len(D.MLTsequence)):
        MLTfrom = int(D.MLTsequence[i])
        if MLTfrom > 24: MLTfrom -=24
        MLTto = int(D.MLTsequence[i]+D.MLT_duration_of_a_bucket)
        if MLTto > 24: MLTto -=24
        ColumnTitles.append( "MLT " + str(MLTfrom) + "-"  + str(MLTto) )
    # define secondary y-axis at the right of the plot
    mySpecs = list()
    for row in range(0, len(D.KPsequence)):
        mySpecs.append( list() )
        for col in range(0, len(D.MLTsequence)):
            mySpecs[row].append( {"secondary_y": True} )

    #make plot
    if D.MLAT_min<=65: print("\nWARNING: LOW LATITUDES MAY HAVE ALMOST ZERO VALUES LEADING TO EMPTY ALTITUDE-PROFILES PLOT.\n")
    fig = make_subplots(rows=len(D.KPsequence), cols=len(D.MLTsequence), shared_xaxes=True, shared_yaxes=True, vertical_spacing=0.035, horizontal_spacing=0.02, subplot_titles=ColumnTitles, specs=mySpecs)
    for aKP in D.KPsequence:
        for aMLT in D.MLTsequence:
            #Means = list()
            Percentiles10 = list()
            Percentiles25 = list()
            Percentiles50 = list()            
            Percentiles75 = list()
            Percentiles90 = list()
            hits  = 0
            
            # compute percentiles
            for anALT in D.ALTsequence:
                Percentiles10.append( Buckets[aKP, anALT, D.MLAT_min, aMLT, "Percentile10"] * MultiplicationFactor )
                Percentiles25.append( Buckets[aKP, anALT, D.MLAT_min, aMLT, "Percentile25"] * MultiplicationFactor )
                Percentiles50.append( Buckets[aKP, anALT, D.MLAT_min, aMLT, "Percentile50"] * MultiplicationFactor )
                Percentiles75.append( Buckets[aKP, anALT, D.MLAT_min, aMLT, "Percentile75"] * MultiplicationFactor )
                Percentiles90.append( Buckets[aKP, anALT, D.MLAT_min, aMLT, "Percentile90"] * MultiplicationFactor )
            
            fig.add_trace( go.Scatter(x=[0]*len(visibleALTsequence), y=visibleALTsequence, mode='lines', fill='tonexty', fillcolor=Color10, line=dict(color='gray',width=1,), showlegend=False), row=D.KPsequence.index(aKP)+1, col=D.MLTsequence.index(aMLT)+1 )
            fig.add_trace( go.Scatter(x=Percentiles10, y=visibleALTsequence, mode='lines', fill='tonexty', fillcolor=Color10, line=dict(color='gray',width=1,), showlegend=False), row=D.KPsequence.index(aKP)+1, col=D.MLTsequence.index(aMLT)+1 )
            fig.add_trace( go.Scatter(x=Percentiles25, y=visibleALTsequence, mode='lines', fill='tonexty', fillcolor=Color25, line=dict(color='gray',width=1,), showlegend=False), row=D.KPsequence.index(aKP)+1, col=D.MLTsequence.index(aMLT)+1 )
            fig.add_trace( go.Scatter(x=Percentiles50, y=visibleALTsequence, mode='lines', fill='tonexty', fillcolor=Color50, line=dict(color='black',width=2,), showlegend=False), row=D.KPsequence.index(aKP)+1, col=D.MLTsequence.index(aMLT)+1 )
            # plot mean
            #fig.add_trace( go.Scatter(x=Means, y=visibleALTsequence, mode='lines', fill='tonexty', fillcolor='black', line=dict(color='black',width=1,), showlegend=False), row=KPsequence.index(aKP)+1, col=MLTsequence.index(aMLT)+1 )
            # plot percentiles
            fig.add_trace( go.Scatter(x=Percentiles75, y=visibleALTsequence, mode='lines', fill='tonexty', fillcolor=Color75, line=dict(color='gray',width=1,), showlegend=False), row=D.KPsequence.index(aKP)+1, col=D.MLTsequence.index(aMLT)+1 )
            fig.add_trace( go.Scatter(x=Percentiles90, y=visibleALTsequence, mode='lines', fill='tonexty', fillcolor=Color90, line=dict(color='gray',width=1,), showlegend=False), row=D.KPsequence.index(aKP)+1, col=D.MLTsequence.index(aMLT)+1,  )
            # add a trace in order to display secondary y-axis at the right
            fig.add_trace( go.Scatter(x=[-1000], y=[-1000], showlegend=False), row=D.KPsequence.index(aKP)+1, col=D.MLTsequence.index(aMLT)+1, secondary_y=True )
            
    # display legends
    fig.add_trace( go.Scatter(name='10th Perc.', x=Percentiles10, y=visibleALTsequence, mode='lines', fill='tonexty', fillcolor=Color10, line=dict(color='gray',width=1,), showlegend=True), row=D.KPsequence.index(aKP)+1, col=D.MLTsequence.index(aMLT)+1 )
    fig.add_trace( go.Scatter(name='25th Perc.', x=Percentiles25, y=visibleALTsequence, mode='lines', fill='tonexty', fillcolor=Color25, line=dict(color='gray',width=1,), showlegend=True), row=D.KPsequence.index(aKP)+1, col=D.MLTsequence.index(aMLT)+1 )
    fig.add_trace( go.Scatter(name='50th Perc.', x=Percentiles50, y=visibleALTsequence, mode='lines', fill='tonexty', fillcolor=Color50, line=dict(color='black',width=2,), showlegend=True), row=D.KPsequence.index(aKP)+1, col=D.MLTsequence.index(aMLT)+1 )
    #fig.add_trace( go.Scatter(name='Mean value', x=Means, y=visibleALTsequence, mode='lines', fill='tonexty', fillcolor='#5cc5ef', line=dict(color='black',width=1,), showlegend=True), row=KPsequence.index(aKP)+1, col=MLTsequence.index(aMLT)+1 )            
    fig.add_trace( go.Scatter(name='75th Perc.', x=Percentiles75, y=visibleALTsequence, mode='lines', fill='tonexty', fillcolor=Color75, line=dict(color='gray',width=1,), showlegend=True), row=D.KPsequence.index(aKP)+1, col=D.MLTsequence.index(aMLT)+1 )
    fig.add_trace( go.Scatter(name='90th Perc.', x=Percentiles90, y=visibleALTsequence, mode='lines', fill='tonexty', fillcolor=Color90, line=dict(color='gray',width=1,), showlegend=True), row=D.KPsequence.index(aKP)+1, col=D.MLTsequence.index(aMLT)+1 )

    fig.update_xaxes( range=x_axes_range, row=1, col=1)
    fig.update_xaxes( range=x_axes_range, row=1, col=2)
    fig.update_xaxes( range=x_axes_range, row=1, col=3)
    fig.update_xaxes( range=x_axes_range, row=1, col=4)
    
    fig.update_xaxes( range=x_axes_range, row=2, col=1)
    fig.update_xaxes( range=x_axes_range, row=2, col=2)
    fig.update_xaxes( range=x_axes_range, row=2, col=3)
    fig.update_xaxes( range=x_axes_range, row=2, col=4)
    
    fig.update_xaxes( range=x_axes_range, row=3, col=1)
    fig.update_xaxes( range=x_axes_range, row=3, col=2)
    fig.update_xaxes( range=x_axes_range, row=3, col=3)
    fig.update_xaxes( range=x_axes_range, row=3, col=4)
    
    for aKP in D.KPsequence:
        fig.update_yaxes( title_text="Altitude (km)", row=D.KPsequence.index(aKP)+1, col=1, side='left', secondary_y=False)
        row_title = "Kp " + str(aKP) + " - "
        if aKP == 0:
            row_title +=  "2"
        elif aKP == 2:
            row_title +=  "4"
        else:
            row_title +=  "9"
        fig.update_yaxes( title_text=row_title, row=D.KPsequence.index(aKP)+1, col=len(D.MLTsequence),  side='right', secondary_y=True, showticklabels=False )
        for aMLT in D.MLTsequence:
            fig.update_yaxes( row=D.KPsequence.index(aKP)+1, col=D.MLTsequence.index(aMLT)+1, secondary_y=True, showticklabels=False )
    #fig.update_xaxes( range=x_axes_range )
    fig.update_yaxes( range=[80, 150], dtick=10 )  
    fig.update_layout( title = SuperTitle + "<br>" + "Altitude Profiles of " + VariableName + " (" + new_units + ")",
                       width=400+len(D.MLTsequence)*250, height=200+200*len(D.KPsequence), showlegend=True, legend_orientation="h", legend_y=-0.04) 

    
    plotly.offline.init_notebook_mode(connected=True)
    plotly.offline.iplot(fig) 

    # plot more zoom versions
    '''
    new_x_axes_range = [x * (2/3) for x in x_axes_range]
    fig.update_xaxes( range=new_x_axes_range )
    plotly.offline.iplot(fig) 
    new_x_axes_range = [x * (1/2) for x in x_axes_range]
    fig.update_xaxes( range=new_x_axes_range )
    plotly.offline.iplot(fig) 
    new_x_axes_range = [x * (3/2) for x in x_axes_range]
    fig.update_xaxes( range=new_x_axes_range )
    plotly.offline.iplot(fig) 
    new_x_axes_range = [x * (2.5) for x in x_axes_range]
    fig.update_xaxes( range=new_x_axes_range )
    plotly.offline.iplot(fig) 
    new_x_axes_range = [x * (10) for x in x_axes_range]
    fig.update_xaxes( range=new_x_axes_range )
    plotly.offline.iplot(fig) 
    '''
