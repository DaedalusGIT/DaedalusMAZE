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


'''
Creates color-spread plots with several sub-plots. Color represents the mean value of the variable at the certain space-time area.
 Arguments:
   VariableName: string, the name of the variable to be plotted
   Buckets: a dictionary containing the data as described in Data.py
   LogScale: if true then the mean values colors will be plotted in logarithmic scale
   SuperTitle: a title to be added at the top of the plot
   ColorMap: the name of the colormap ot use for ploting
             values: http://www.fabiocrameri.ch/colourmaps.php: acton, bamako, batlow, berlin, bilbao, broc, buda, cork, davos, devon, grayC, hawaii, imola, lajolla, lapaz, lisbon, nuuk, oleron , oslo, roma, tofino, tokyo, turku, vik - romaO, brocO, corkO, vikO 
             values: plotly: Blackbody, Bluered, Blues, Earth, Electric, Greens, Greys, Hot, Jet, Picnic, Portland, Rainbow, RdBu, Reds, Viridis, YlGnBu, YlOrRd
 Returns: -
'''
def plotColorSpread_perKpRange( VariableName, Buckets, LogScale=True, SuperTitle="", ColorMap=vik ):    
    print( "Plot creation started", datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S") )
    #global ALTsequence
    MultiplicationFactor = 1
        
    # resize the Altitude range:
    #oldALTsequence  = D.ALTsequence.copy()
    #D.ALTsequence   = list( range(100, 180, int(D.ALT_distance_of_a_bucket)) )
    
    # construct the column titles 
    ColumnTitles = list()    
    for i in range(0, len(D.ALTsequence)):
        ColumnTitles.append( "<b>" + str(D.ALTsequence[i]) + "-"  + str(D.ALTsequence[i]+D.ALT_distance_of_a_bucket) + "km" + "</b>")
        
    #make plot
    HitsStr = ""
    fig = make_subplots(rows=len(D.KPsequence), cols=len(D.ALTsequence), shared_xaxes=True, shared_yaxes=True, vertical_spacing=0.035, horizontal_spacing=0.015, subplot_titles=ColumnTitles)
    
    # bundle data, min and max values
    allMeans_min = allMeans_logscale_min = 999999
    allMeans_max = allMeans_logscale_max = -99999
    for aKP in D.KPsequence:
        for anALT in D.ALTsequence:
            Means = np.zeros( ( len(D.MLATsequence), len(D.MLTsequence)) )
            hits = 0
            for aMLT in D.MLTsequence:
                for aMagLat in D.MLATsequence:
                    hits += float( Buckets[(aKP, anALT, aMagLat, aMLT, "Len")]  )
                    i = D.MLATsequence.index(aMagLat)
                    j = D.MLTsequence.index(aMLT)
                    if Buckets[(aKP, anALT, aMagLat, aMLT, "Len")] > 0: 
                        Means[i, j] = Buckets[(aKP, anALT, aMagLat, aMLT, "Sum")] / Buckets[(aKP, anALT, aMagLat, aMLT, "Len")] 
                        
            #print( "Kp = ", aKP, "ALT =", anALT, "   Hits =", hits)
            
            # change units
            Means *= MultiplicationFactor
            
            # logScale
            Means_logscale = np.zeros( Means.shape )
            for i in range(0, len(Means)):
                for j in range(0, len(Means[i])):
                    if Means[i][j] > 0:
                        Means_logscale[i][j] = np.log10(Means[i][j])
                    else:
                        Means_logscale[i][j] = None

            # min-max
            Means_min =  999999
            Means_max = -999999
            for i in range(0, len(Means)):
                for j in range(0, len(Means[i])):
                    if Means[i][j] is not None and Means[i][j]!=0 and Means_min>Means[i][j]: Means_min = Means[i][j]
                    if Means[i][j] is not None and Means[i][j]!=0 and Means_max<Means[i][j]: Means_max = Means[i][j]
            Means_logscale_min = np.nanmin(Means_logscale)
            Means_logscale_max = np.nanmax(Means_logscale)
            if Means_logscale_min==float("-inf"): Means_logscale_min = 0
            if Means_logscale_max==float("-inf"): Means_logscale_max = 0
            if Means_min is not None and allMeans_min > Means_min: allMeans_min = Means_min
            if Means_max is not None and allMeans_max < Means_max: allMeans_max = Means_max
            if allMeans_logscale_min > Means_logscale_min: allMeans_logscale_min = Means_logscale_min
            if allMeans_logscale_max < Means_logscale_max: allMeans_logscale_max = Means_logscale_max
                
            #Means.append( Means[-1] )
            #Means_logscale.append( Means_logscale[-1] )
            Means = np.hstack((Means, np.tile(Means[:, [-1]], 1)))
            Means = np.vstack([Means, Means[-1]])
            Means_logscale = np.hstack((Means_logscale, np.tile(Means_logscale[:, [-1]], 1)))
            Means_logscale = np.vstack([Means_logscale, Means_logscale[-1]])
            MLTseq_for_plot = D.MLTsequence + [D.MLTsequence[-1]+D.MLT_duration_of_a_bucket]
            MagLatseq_for_plot = D.MLATsequence + [D.MLATsequence[-1]+D.MLAT_degrees_of_a_bucket]

            # force all min/max equal to tiegcm's
            #allMeans_logscale_min = -1.801815958043412 
            #allMeans_logscale_max = 1.1062579615522905 
            #allMeans_min = 0.0
            #allMeans_max = 12.771972111394563
    
            # plot heatmap
            if LogScale:
                fig.add_trace( go.Heatmap(z=Means_logscale.tolist(), x=MLTseq_for_plot, y=MagLatseq_for_plot, zsmooth='best', showlegend=False, coloraxis="coloraxis1"), row=D.KPsequence.index(aKP)+1, col=D.ALTsequence.index(anALT)+1,  )
            else:
                fig.add_trace( go.Heatmap(z=Means.tolist(), x=MLTseq_for_plot, y=MagLatseq_for_plot, zsmooth='best', showlegend=False, coloraxis="coloraxis1"), row=D.KPsequence.index(aKP)+1, col=D.ALTsequence.index(anALT)+1,  )

    print("* * * * * * * * * * * * * * * * * * * * * * * * * *  *")
    print("allMeans_logscale_min", "allMeans_logscale_max", "allMeans_min", "allMeans_max" )
    print(allMeans_logscale_min, allMeans_logscale_max, allMeans_min, allMeans_max )
    print("* * * * * * * * * * * * * * * * * * * * * * * * * *  *")
    
    fig.update_layout(coloraxis=dict(colorscale=ColorMap), showlegend=False) #fig.update_traces(zmin=0.07687949e-02, zmax=3.07687949e-01, selector=dict(type="heatmap"))
    # display titles
    fig.update_yaxes( title_text="<b>" + "Kp 0-3" + "</b>" + "<br><br>" + "Magnetic Latitude (deg)", row=1, col=1, side='left', secondary_y=False)
    fig.update_yaxes( title_text="<b>" + "Kp 3-9" + "</b>" + "<br><br>" + "Magnetic Latitude (deg)", row=2, col=1, side='left', secondary_y=False)
    for anAlt in D.ALTsequence: fig.update_xaxes( title_text="MLT (hours)", row=len(D.KPsequence), col=D.ALTsequence.index(anAlt)+1)
        
    # Set the same min/max for all figures
    fig.update_traces(zmin=allMeans_min, zmax=allMeans_max)
    # tick values at the color bar
    if LogScale:
        my_Tickvals    = np.linspace(allMeans_min, allMeans_max, 5, endpoint=True)
        my_logTickvals = list()
        my_Ticktexts   = list()
        for t in range( 0, len(my_Tickvals) ):
            try:
                my_logTickvals.append( math.log10(my_Tickvals[t]) )
                my_Ticktexts.append( "{:.3e}".format(my_Tickvals[t]) )                
            except Exception as e:
                #print(e)
                pass
        fig.update_layout(coloraxis_colorbar=dict( title="Log scale",  tickvals=my_logTickvals,  ticktext=my_Ticktexts, ))
    #
    fig.update_yaxes( range=[D.MLAT_min,  D.MLAT_max], dtick=D.MLAT_degrees_of_a_bucket )
    fig.update_xaxes( range=[D.MLT_min, D.MLT_max], dtick=D.MLT_duration_of_a_bucket )
    # font
    fig.update_layout(font_family="Helvetica",)
    # Set title
    mainTitle = SuperTitle
    fig.update_layout( title = mainTitle, width=400+len(D.ALTsequence)*160, height=250+200*len(D.KPsequence), showlegend=True, legend_orientation="h", legend_y=-0.04) 
    plotly.offline.init_notebook_mode(connected=True)
    plotly.offline.iplot(fig)
        
    # resize the Altitude range to its initial values
    #D.ALTsequence = oldALTsequence


    
    
    
    
    
    
    
    
    
    
    
def plotColorSpread_MLTvsMAGLAT_slice( VariableName, Buckets, LogScale=True, SuperTitle="", ColorMap=vik ):    
    print( "Plot creation started", datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S") )
    #global ALTsequence
    MultiplicationFactor = 1
        
    # resize the Altitude range:
    oldALTsequence = D.ALTsequence.copy()
    D.ALTsequence    = [120]
    
    # construct the column titles 
    ColumnTitles = list()    
    for i in range(0, len(D.ALTsequence)):
        ColumnTitles.append( "<b>Altitude " + str(D.ALTsequence[i]) + "-"  + str(D.ALTsequence[i]+D.ALT_distance_of_a_bucket) + "km" + "</b>")
        
    #make plot
    HitsStr = ""
    fig = make_subplots(rows=1, cols=1, shared_xaxes=True, shared_yaxes=True, vertical_spacing=0.035, horizontal_spacing=0.015, subplot_titles=ColumnTitles)
    
    # bundle data, min and max values
    allMeans_min = allMeans_logscale_min = 999999
    allMeans_max = allMeans_logscale_max = -99999
    for anALT in D.ALTsequence:
        Means = np.zeros( ( len(D.MLATsequence), len(D.MLTsequence)) )
        hits = 0
        for aMLT in D.MLTsequence:
            for aMagLat in D.MLATsequence:
                i = D.MLATsequence.index(aMagLat)
                j = D.MLTsequence.index(aMLT)
                for aKP in D.KPsequence:
                    hits += float( Buckets[(aKP, anALT, aMagLat, aMLT, "Len")]  )
                    bucket_sum = 0
                    bucket_len = 0
                    if Buckets[(aKP, anALT, aMagLat, aMLT, "Len")] > 0: 
                        bucket_sum += Buckets[(aKP, anALT, aMagLat, aMLT, "Sum")]
                        bucket_len += Buckets[(aKP, anALT, aMagLat, aMLT, "Len")]
                        Means[i, j] = bucket_sum / bucket_len
                        
        #print( "Kp = ", aKP, "ALT =", anALT, "   Hits =", hits)
            
        # change units
        Means *= MultiplicationFactor
            
        # logScale
        Means_logscale = np.zeros( Means.shape )
        for i in range(0, len(Means)):
            for j in range(0, len(Means[i])):
                if Means[i][j] > 0:
                    Means_logscale[i][j] = np.log10(Means[i][j])
                else:
                    Means_logscale[i][j] = None

        # min-max
        Means_min =  999999
        Means_max = -999999
        for i in range(0, len(Means)):
            for j in range(0, len(Means[i])):
                if Means[i][j] is not None and Means[i][j]!=0 and Means_min>Means[i][j]: Means_min = Means[i][j]
                if Means[i][j] is not None and Means[i][j]!=0 and Means_max<Means[i][j]: Means_max = Means[i][j]
        Means_logscale_min = np.nanmin(Means_logscale)
        Means_logscale_max = np.nanmax(Means_logscale)
        if Means_logscale_min==float("-inf"): Means_logscale_min = 0
        if Means_logscale_max==float("-inf"): Means_logscale_max = 0
        if Means_min is not None and allMeans_min > Means_min: allMeans_min = Means_min
        if Means_max is not None and allMeans_max < Means_max: allMeans_max = Means_max
        if allMeans_logscale_min > Means_logscale_min: allMeans_logscale_min = Means_logscale_min
        if allMeans_logscale_max < Means_logscale_max: allMeans_logscale_max = Means_logscale_max                

        # force all min/max equal to tiegcm's
        #allMeans_logscale_min = -1.801815958043412 
        #allMeans_logscale_max = 1.1062579615522905 
        #allMeans_min = 0.0
        #allMeans_max = 12.771972111394563
    
        # plot heatmap
        if LogScale:
            fig.add_trace( go.Heatmap(z=Means_logscale.tolist(), x=D.MLTsequence, y=D.MLATsequence, showlegend=False, coloraxis="coloraxis1"), row=1, col=D.ALTsequence.index(anALT)+1,  )
        else:
            fig.add_trace( go.Heatmap(z=Means.tolist(), x=D.MLTsequence, y=D.MLATsequence, showlegend=False, coloraxis="coloraxis1"), row=1, col=D.ALTsequence.index(anALT)+1,  )

    print("* * * * * * * * * * * * * * * * * * * * * * * * * *  *")
    print("allMeans_logscale_min", "allMeans_logscale_max", "allMeans_min", "allMeans_max" )
    print(allMeans_logscale_min, allMeans_logscale_max, allMeans_min, allMeans_max )
    print("* * * * * * * * * * * * * * * * * * * * * * * * * *  *")
    
    fig.update_layout(coloraxis=dict(colorscale=ColorMap), showlegend=False) #fig.update_traces(zmin=0.07687949e-02, zmax=3.07687949e-01, selector=dict(type="heatmap"))
    # display titles
    fig.update_yaxes( title_text="<b>" + "Kp 0-9" + "</b>" + "<br><br>" + "Magnetic Latitude (deg)", row=1, col=1, side='left', secondary_y=False)
    for anAlt in D.ALTsequence: fig.update_xaxes( title_text="MLT (hours)", row=1, col=D.ALTsequence.index(anAlt)+1)
        
    # Set the same min/max for all figures
    fig.update_traces(zmin=allMeans_min, zmax=allMeans_max)
    # tick values at the color bar
    if LogScale:
        my_Tickvals    = np.linspace(allMeans_min, allMeans_max, 5, endpoint=True)
        my_logTickvals = list()
        my_Ticktexts   = list()
        for t in range( 0, len(my_Tickvals) ):
            try:
                my_logTickvals.append( math.log10(my_Tickvals[t]) )
                my_Ticktexts.append( "{:.3e}".format(my_Tickvals[t]) )                
            except Exception as e:
                #print(e)
                pass
        fig.update_layout(coloraxis_colorbar=dict( title="Log scale",  tickvals=my_logTickvals,  ticktext=my_Ticktexts, ))
    #
    fig.update_yaxes( range=[D.MLAT_min,  D.MLAT_max], dtick=30 )
    fig.update_xaxes( range=[D.MLT_min, D.MLT_max], dtick=6 )
    # font
    fig.update_layout(font_family="Helvetica",)
    # Set title
    mainTitle = SuperTitle
    fig.update_layout( title = mainTitle, width=1200, height=800, showlegend=True, legend_orientation="h", legend_y=-0.04) 
    plotly.offline.init_notebook_mode(connected=True)
    plotly.offline.iplot(fig)
        
    # resize the Altitude range to its initial values
    D.ALTsequence = oldALTsequence
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
def plotColorSpread_AltAxes( VariableName, Buckets, LogScale=True, SuperTitle="", ColorMap=lajolla):
    print( "Plot creation started", datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S") )
    MultiplicationFactor = 1
        
    # construct the column titles 
    ColumnTitles = ["MLT 23:00-01:00", "MLT 11:00-13:00"]
        
    #make plot
    HitsStr = ""
    fig = make_subplots(rows=1, cols=2, shared_xaxes=True, shared_yaxes=True, vertical_spacing=0.03, horizontal_spacing=0.03, subplot_titles=ColumnTitles)
    
    # bundle data, min and max values
    allMeans_min = allMeans_logscale_min = 999999
    allMeans_max = allMeans_logscale_max = -99999
    
    MLTmoments = [ [0,23], [11,12] ]
    for m in range(0, len(MLTmoments)):
        # compute a mean for each MagLat-Altitude pair
        Sums  = np.zeros( ( len(D.ALTsequence), len(D.MLATsequence)) )
        Lens  = np.zeros( ( len(D.ALTsequence), len(D.MLATsequence)) )
        Means = np.zeros( ( len(D.ALTsequence), len(D.MLATsequence)) )
        #print( D.KPsequence)
        #print( D.ALTsequence)
        #print( D.MLATsequence )
        #print( D.MLTsequence )
        for aMagLat in D.MLATsequence:
            for anALT in D.ALTsequence:
                i = D.ALTsequence.index(anALT)
                j = D.MLATsequence.index(aMagLat)
                Sums[i, j] += Buckets[(0, anALT, aMagLat, MLTmoments[m][0], "Sum")] 
                Sums[i, j] += Buckets[(3, anALT, aMagLat, MLTmoments[m][0], "Sum")] 
                Sums[i, j] += Buckets[(0, anALT, aMagLat, MLTmoments[m][1], "Sum")] 
                Sums[i, j] += Buckets[(3, anALT, aMagLat, MLTmoments[m][1], "Sum")] 
                Lens[i, j] += Buckets[(0, anALT, aMagLat, MLTmoments[m][0], "Len")] 
                Lens[i, j] += Buckets[(3, anALT, aMagLat, MLTmoments[m][0], "Len")] 
                Lens[i, j] += Buckets[(0, anALT, aMagLat, MLTmoments[m][1], "Len")] 
                Lens[i, j] += Buckets[(3, anALT, aMagLat, MLTmoments[m][1], "Len")] 
                if Lens[i, j] > 0:
                    Means[i, j] = Sums[i, j] / Lens[i, j]
            
        # change units
        Means *= MultiplicationFactor
            
        # logScale
        Means_logscale = np.zeros( Means.shape )
        for i in range(0, len(Means)):
            for j in range(0, len(Means[i])):
                if Means[i][j] > 0:
                    Means_logscale[i][j] = np.log10(Means[i][j])
                else:
                    Means_logscale[i][j] = None

        # min-max
        Means_min =  999999
        Means_max = -999999
        for i in range(0, len(Means)):
            for j in range(0, len(Means[i])):
                if Means[i][j] is not None and Means[i][j]!=0 and Means_min>Means[i][j]: Means_min = Means[i][j]
                if Means[i][j] is not None and Means[i][j]!=0 and Means_max<Means[i][j]: Means_max = Means[i][j]
        Means_logscale_min = np.nanmin(Means_logscale)
        Means_logscale_max = np.nanmax(Means_logscale)
        if Means_logscale_min==float("-inf"): Means_logscale_min = 0
        if Means_logscale_max==float("-inf"): Means_logscale_max = 0
        if Means_min is not None and allMeans_min > Means_min: allMeans_min = Means_min
        if Means_max is not None and allMeans_max < Means_max: allMeans_max = Means_max
        if allMeans_logscale_min > Means_logscale_min: allMeans_logscale_min = Means_logscale_min
        if allMeans_logscale_max < Means_logscale_max: allMeans_logscale_max = Means_logscale_max                

        # force all min/max equal to tiegcm's
        #allMeans_logscale_min = -1.801815958043412 
        #allMeans_logscale_max = 1.1062579615522905 
        #allMeans_min = 0.0
        #allMeans_max = 12.771972111394563
    
        # plot heatmap # zsmooth='best',
        if LogScale:
            fig.add_trace( go.Heatmap(z=Means_logscale.tolist(), x=D.MLATsequence, y=D.ALTsequence,  showlegend=False, coloraxis="coloraxis1"), row=1, col=m+1,  )
        else:
            fig.add_trace( go.Heatmap(z=Means.tolist(), x=D.MLATsequence, y=D.ALTsequence, showlegend=False, coloraxis="coloraxis1"), row=1, col=m+1,  )

    print("* * * * * * * * * * * * * * * * * * * * * * * * * *  *")
    print("allMeans_logscale_min", "allMeans_logscale_max", "allMeans_min", "allMeans_max" )
    print(allMeans_logscale_min, allMeans_logscale_max, allMeans_min, allMeans_max )
    print("* * * * * * * * * * * * * * * * * * * * * * * * * *  *")
    
    fig.update_layout(coloraxis=dict(colorscale=ColorMap), showlegend=False)
    # display titles
    fig.update_yaxes( title_text="Altitude (Km)", row=1, col=1, side='left', secondary_y=False)
    fig.update_xaxes( title_text="Magnetic Latitude (degrees)", row=1, col=1)
        
    # Set the same min/max for all figures
    fig.update_traces(zmin=allMeans_min, zmax=allMeans_max)
    # tick values at the color bar
    if LogScale:
        my_Tickvals    = np.linspace(allMeans_min, allMeans_max, 5, endpoint=True)
        my_logTickvals = list()
        my_Ticktexts   = list()
        for t in range( 0, len(my_Tickvals) ):
            try:
                my_logTickvals.append( math.log10(my_Tickvals[t]) )
                my_Ticktexts.append( "{:.3e}".format(my_Tickvals[t]) )                
            except Exception as e:
                #print(e)
                pass
        fig.update_layout(coloraxis_colorbar=dict( title="Log scale",  tickvals=my_logTickvals,  ticktext=my_Ticktexts, ))
    #
    fig.update_yaxes( range=[D.Alt_min,  D.Alt_max], dtick=50 )
    fig.update_xaxes( range=[D.MLAT_min,  D.MLAT_max], dtick=30)
    # font
    fig.update_layout(font_family="Helvetica",)
    # Set title
    mainTitle = SuperTitle
    fig.update_layout( title = mainTitle, width=1800, height=740, showlegend=True, legend_orientation="h", legend_y=-0.04) 
    plotly.offline.init_notebook_mode(connected=True)
    plotly.offline.iplot(fig)
    print("Done")    
    