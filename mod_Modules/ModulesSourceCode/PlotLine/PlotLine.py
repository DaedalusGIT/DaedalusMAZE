
import plotly
import plotly.plotly as py
import plotly.graph_objs as go

"""
Plots X-Y pair values connected by a line.
  Xvalues: 1D array 
  Yvalues: 1D array. Must be the same size as Xvalues array.
  TitleMain: the main title of the graph
  TitleXaxis: the title for the X axis
  TitleYaxis: the title for the Y axis
"""
def PlotLine( Xvalues, Yvalues, TitleMain, TitleXaxis, TitleYaxis ):
    '''
    For future may add these arguments: 
        inside trace: name, for the name of the line, especialy in case there are more lines
                      mode, 'lines+markers', 'lines', 'markers'
                      line, example:     line = dict(color = ('rgb(22, 96, 167)'), width = 4, dash = 'dot')
  
    ''' 
    trace = go.Scatter( x = Xvalues, y = Yvalues )
    
    layout = dict( title = TitleMain,
                   xaxis = dict(title = TitleXaxis),
                   yaxis = dict(title = TitleYaxis), )
    
    data = [trace]
    fig = dict(data=data, layout=layout)
    plotly.offline.init_notebook_mode(connected=True)
    plotly.offline.iplot(fig)
    return 


