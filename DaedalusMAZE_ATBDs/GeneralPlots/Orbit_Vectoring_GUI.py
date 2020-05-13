# -*- coding: utf-8 -*-
"""
Created on Tue May 12 22:11:26 2020

@author: Stelios
"""

from ipywidgets import*
import ipywidgets as widgets


def run_run(x): 
    output_file_v1(datetime_value,lat_value,lon_value,min_alt_value,max_alt_value,int(vector_step_value),ap_ind_value,f10_7_ind_value)

style = {'description_width': '150px'}
layout1stcolumn = {'width': '500px'}
layout2ndcolumn= {'width': '250px'}
layout3rdcolumn= {'width': '100px'}
button_layout={'width': '150px'}

#lat widget
lat=widgets.BoundedFloatText(
    value=-90,
    min=-90,
    max=90,
    step=0.5,
    description='Lat_GEOD(deg):',
    description_tooltip='Insert the latitude for the output file',
    disabled=False,
    layout=layout1stcolumn,
    style=style,   
)


lon=widgets.BoundedFloatText(
    value=0,
    min=-180,
    max=180,
    step=0.5,
    description='Lon_GEOD(deg):',
    description_tooltip='Insert the longitude for the output file',
    disabled=False,
    layout=layout1stcolumn,
    style=style,
)

min_alt=widgets.BoundedIntText(
    value=0,
    min=100,
    max=10000,
    step=5,
    description='Min_Altitude(km):',
    description_tooltip='Insert the minimum altitude for the altitude range',
    disabled=False,
    layout=layout1stcolumn,
    style=style,
    
)

max_alt=widgets.BoundedIntText(
    value=0,
    min=1000,
    max=10000,
    step=5,
    description='Max_Altitude(km):',
    description_tooltip='Insert the maximum altitude for the altitude range',
    disabled=False,
    layout=layout1stcolumn,
    style=style,
)

vector_step=widgets.Text(
    value='10',
    placeholder='Enter the vector step',
    description='Vector_step:',
    description_tooltip='Insert the step of the Vectors. The higher the value the lower the vectors to be shown on plot',
    disabled=False,
    layout=layout1stcolumn,
    style=style,
)

filename=widgets.Text(
    placeholder='Path_Orbit_Vector_Data',
    description_tooltip='Enter the full path with the filename, like this: /../DaedalusMAZE_ATBDs/GeneralPlots/Orbit_Vectoring_GUI.csv ',
    description='Orbit_Filename:',
    disabled=False,
    layout=layout1stcolumn,
    style=style,
)

run_button=widgets.Button(
    value=False,
    description='Plot Orbit_Vectors',
    disabled=False,
    button_style='success', # 'success', 'info', 'warning', 'danger' or ''
    tooltip='Press to create the orbit plot with vectors',
    icon='check',
    layout=button_layout,
    style=style,
    
)

#lat widget value
lat_value=lat.value
#lon widget value
lon_value=lon.value
#min_alt widget value
min_alt_value=min_alt.value
#max_alt widget value
max_alt_value=max_alt.value
#vector_step widget value
vector_step_value=vector_step.value
#datetime widget value
filename_value=filename.value

y1=widgets.VBox([filename,lon,lat,min_alt,max_alt,vector_step])
y3=widgets.VBox([run_button])

run_button.on_click(run_run)
box=widgets.HBox([y1,y3])
display(box)
