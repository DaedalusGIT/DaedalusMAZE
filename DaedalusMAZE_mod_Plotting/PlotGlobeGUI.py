import sys
sys.path.insert(1, '../SourceCode/ModulesSourceCode/PlotGlobe')
import PlotGlobe as PlotGlobe
sys.path.insert(1, '../SourceCode/')
import DaedalusGlobals as DaedalusGlobals


from netCDF4 import Dataset
import ipywidgets as w
import os
import glob

# A list of available colorscales
ColorScales = [
    'aggrnyl', 'agsunset', 'algae', 'amp', 'armyrose', 'balance', 'blackbody', 'bluered', 'blues', 'blugrn', 'bluyl', 'brbg',
    'brwnyl', 'bugn', 'bupu', 'burg', 'burgyl', 'cividis', 'curl', 'darkmint', 'deep', 'delta', 'dense', 'earth', 'edge', 'electric',
    'emrld', 'fall', 'geyser', 'gnbu', 'gray', 'greens', 'greys', 'haline', 'hot', 'hsv', 'ice', 'icefire', 'inferno', 'jet',
    'magenta', 'magma', 'matter', 'mint', 'mrybm', 'mygbm', 'oranges', 'orrd', 'oryel', 'peach', 'phase', 'picnic', 'pinkyl', 'piyg',
    'plasma', 'plotly3', 'portland', 'prgn', 'pubu', 'pubugn', 'puor', 'purd', 'purp', 'purples', 'purpor', 'rainbow', 'rdbu', 'rdgy',
    'rdpu', 'rdylbu', 'rdylgn', 'redor', 'reds', 'solar', 'spectral', 'speed', 'sunset', 'sunsetdark', 'teal', 'tealgrn', 'tealrose',
    'tempo', 'temps', 'thermal', 'tropic', 'turbid', 'twilight', 'viridis', 'ylgn', 'ylgnbu', 'ylorbr', 'ylorrd'
]


# returns a list of all the variables and their properties in the NetCDF file
def Get_NetCDF_VariablesInfo( NetCDF_filename ):
    NetCDF = Dataset(NetCDF_filename, 'r')
    ##
    AllVariables = [""]
    for V in NetCDF.variables:
        Vdtype = Vunits = Vlongname = ""
        try:
            Vdtype = "(" + NetCDF.variables[V].dtype + ") "
            Vdtype = Vdtype.replace("() ","")
        except:
            pass
        try:
            Vunits = "[" + NetCDF.variables[V].units + "] "
            Vunits = Vunits.replace("[] ","")
        except:
            pass
        try:
            Vlongname = "(" +NetCDF.variables[V].long_name + ") "
            Vlongname = Vlongname.replace("() ","")
        except:
            pass
        AllVariables.append( V + "   " + Vlongname +  Vunits + Vdtype )
    NetCDF.close()
    return AllVariables

################################################################################################################
####################################### Event Listeners ########################################################
################################################################################################################
def Orbit_Filename_onChange(change):
    if change['type']=='change' and change['name']=='value' and len(change['new'])>0:
        # load available variables of the NetCDF file and propose them to user
        filename = DaedalusGlobals.Orbit_Files_Path + "DAED_ORB_"
        filename += Orbit_EventDropdown.value + "_"
        filename += Orbit_DataFormatDropdown.value + "_"
        filename += Orbit_PerigeeDropdown.value + "_"
        filename += Orbit_LatitudeDropdown.value + "_"
        filename += Orbit_SamplingDropdown.value + "_"
        filename += Orbit_SpacecraftDropdown.value
        filename += ".nc"
        try:
            AllVariables = Get_NetCDF_VariablesInfo( filename )
            Orbit_VariableDropdown.options = AllVariables
        except:
            pass
            

def Surface_TypeDropdown_onChange(change):
    if change['type']=='change' and change['name']=='value' and len(change['new'])>0:
        # display available event folder names at the dropdown menu
        AllEvents = glob.glob( DaedalusGlobals.AllData_Files_Path+str(change['new'])+"/"+str(change['new'])[0:4]+"*" )
        for i in range(0, len(AllEvents)):
            AllEvents[i] = AllEvents[i][ AllEvents[i].rfind('/')+1 : ]                
        AllEvents.sort()
        Surface_EventDropdown.options = AllEvents
def Surface_EventDropdown_onChange(change):
    if change['type']=='change' and change['name']=='value' and len(change['new'])>0:
        # display available data file names at the dropdown menu
        AllDataFiles = glob.glob( DaedalusGlobals.AllData_Files_Path +Surface_TypeDropdown.value+"/"+str(change['new'])+"/*.nc" )
        for i in range(0, len(AllDataFiles)):
            AllDataFiles[i] = AllDataFiles[i][ AllDataFiles[i].rfind('/')+1 : ]
        AllDataFiles.sort()
        Surface_DatafileDropdown.options = AllDataFiles
def Surface_DatafileDropdown_onChange(change):
    if change['type']=='change' and change['name']=='value' and len(change['new'])>0:
        # load available variables of the NetCDF file and propose them to user
        AllVariables = Get_NetCDF_VariablesInfo( DaedalusGlobals.AllData_Files_Path +Surface_TypeDropdown.value+"/"+Surface_EventDropdown.value+"/"+str(change['new']) )
        Surface_VariableDropdown.options = AllVariables
    
def Vectors_TypeDropdown_onChange(change):
    if change['type']=='change' and change['name']=='value' and len(change['new'])>0:
        # display available event folder names at the dropdown menu
        AllEvents = glob.glob( DaedalusGlobals.AllData_Files_Path+str(change['new'])+"/"+str(change['new'])[0:4]+"*" )
        for i in range(0, len(AllEvents)):
            AllEvents[i] = AllEvents[i][ AllEvents[i].rfind('/')+1 : ]                
        AllEvents.sort()
        Vectors_EventDropdown.options = AllEvents
def Vectors_EventDropdown_onChange(change):
    if change['type']=='change' and change['name']=='value' and len(change['new'])>0:
        # display available data file names at the dropdown menu
        AllDataFiles = glob.glob( DaedalusGlobals.AllData_Files_Path +Vectors_TypeDropdown.value+"/"+str(change['new'])+"/*.nc" )
        for i in range(0, len(AllDataFiles)):
            AllDataFiles[i] = AllDataFiles[i][ AllDataFiles[i].rfind('/')+1 : ]        
        AllDataFiles.sort()
        Vectors_DatafileDropdown.options = AllDataFiles
def Vectors_DatafileDropdown_onChange(change):
    if change['type']=='change' and change['name']=='value' and len(change['new'])>0:
        # load available variables of the NetCDF file and propose them to user
        AllVariables = Get_NetCDF_VariablesInfo( DaedalusGlobals.AllData_Files_Path +Vectors_TypeDropdown.value+"/"+Vectors_EventDropdown.value+"/"+str(change['new']) )
        Vectors_VarXDropdown.options = AllVariables
        Vectors_VarYDropdown.options = AllVariables
        Vectors_VarZDropdown.options = AllVariables
                
def Orbit_Var_onChange(change):
    if Orbit_VariableDropdown.value is not None and len(Orbit_VariableDropdown.value)>0:
        filename = "DAED_ORB_"
        filename += Orbit_EventDropdown.value + "_"
        filename += Orbit_DataFormatDropdown.value + "_"
        filename += Orbit_PerigeeDropdown.value + "_"
        filename += Orbit_LatitudeDropdown.value + "_"
        filename += Orbit_SamplingDropdown.value + "_"
        filename += Orbit_SpacecraftDropdown.value
        varname = Orbit_VariableDropdown.value[:Orbit_VariableDropdown.value.find('(')-1].strip()
        OrbitAccordion.set_title(0, "Orbit: " +filename + " [" + varname + "]" )
    else:
        OrbitAccordion.set_title(0, "Click here to plot an Orbit" )
def Surface_Var_onChange(change):
    if Surface_VariableDropdown.value is not None and len(Surface_VariableDropdown.value)>0:
        varname = Surface_VariableDropdown.value[:Surface_VariableDropdown.value.find('(')-1].strip()
        SurfaceAccordion.set_title(0, "Surface: " + Surface_DatafileDropdown.value + " [" + varname + "]" )
        # auto fill units, pressure level limits, timestep limits
        NetCDF_filename = DaedalusGlobals.AllData_Files_Path +Surface_TypeDropdown.value+"/"+Surface_EventDropdown.value+"/"+Surface_DatafileDropdown.value
        NetCDF = Dataset(NetCDF_filename, 'r')
        Surface_ColorbartitleText.value = NetCDF.variables[varname].units
        if len(NetCDF.variables[varname].shape) >= 2:
            Surface_TimestepText.max = NetCDF.variables[varname].shape[0] - 1
            Surface_PressurelevelText.max = NetCDF.variables[varname].shape[1] - 1 
        NetCDF.close()        
    else:
        SurfaceAccordion.set_title(0, "Click here to plot a variable onto a Surface above the Earth" )
def Vectors_Var_onChange(change):
    VarSelected = False
    if Vectors_VarXDropdown.value is not None and len(Vectors_VarXDropdown.value)>0:
        if Vectors_VarYDropdown.value is not None and len(Vectors_VarYDropdown.value)>0:
            if Vectors_VarZDropdown.value is not None and len(Vectors_VarZDropdown.value)>0:
                VarSelected = True
    if VarSelected:
        s = "Vectors: " + Vectors_DatafileDropdown.value + " ["
        s +=  Vectors_VarXDropdown.value[:Vectors_VarXDropdown.value.find('(')-1].strip() + ", "
        s +=  Vectors_VarYDropdown.value[:Vectors_VarYDropdown.value.find('(')-1].strip() + ", "
        s +=  Vectors_VarZDropdown.value[:Vectors_VarZDropdown.value.find('(')-1].strip() + "]"
        VectorsAccordion.set_title(0, s )
        # auto fill units, pressure level limits, timestep limits
        varname = Vectors_VarXDropdown.value[:Vectors_VarXDropdown.value.find('(')-1].strip()
        NetCDF_filename = DaedalusGlobals.AllData_Files_Path +Vectors_TypeDropdown.value+"/"+Vectors_EventDropdown.value+"/"+Vectors_DatafileDropdown.value
        NetCDF = Dataset(NetCDF_filename, 'r')
        Vectors_ColorbartitleText.value = NetCDF.variables[varname].units
        if len(NetCDF.variables[varname].shape) >= 2:
            Vectors_TimestepText.max = NetCDF.variables[varname].shape[0] - 1
            Vectors_PressurelevelText.max = NetCDF.variables[varname].shape[1] - 1
        NetCDF.close()         
    else:
        VectorsAccordion.set_title(0, "Click here to plot a Vector field" )

def plot_button_clicked( b ):
    SurfaceFilename = SurfaceColorbarTitle = SurfaceColorscaleName = SurfaceMultiplicationFactor = ""
    SurfaceOpacity  = SurfaceTimestep      = SurfacePressureLevel  = ""
    OrbitFilename   = OrbitColorbarTitle   = OrbitColorscaleName   = OrbitMultiplicationFactor   = ""
    VectorsFilename = VectorsColorbarTitle = VectorsColorscaleName = VectorsMultiplicationFactor = ""
    VectorConeSize  = VectorConeOpacity    = VectorsTimestep       = VectorsPressureLevel        = ""
    PlotTitle = Preferences_Plottitle.value
    LogScale = Preferences_Logscale.value
    OneSizeVectorCones = False
    ####
    SurfaceVariableToPlot = Surface_VariableDropdown.value[:Surface_VariableDropdown.value.find('(')-1].strip()
    if( len(SurfaceVariableToPlot) > 0 ):
        SurfaceFilename  = DaedalusGlobals.AllData_Files_Path + Surface_TypeDropdown.value + "/" + Surface_EventDropdown.value
        SurfaceFilename += "/" + Surface_DatafileDropdown.value
        SurfaceColorbarTitle = Surface_ColorbartitleText.value
        SurfaceColorscaleName = Surface_ColorscaleText.value
        SurfaceOpacity = Surface_OpacityText.value
        SurfaceMultiplicationFactor = Surface_MultiplyText.value
        SurfaceTimestep = Surface_TimestepText.value
        SurfacePressureLevel = Surface_PressurelevelText.value
    ####
    OrbitVariableToPlot = Orbit_VariableDropdown.value[:Orbit_VariableDropdown.value.find('(')-1].strip()
    if( len(OrbitVariableToPlot) > 0 ):
        OrbitFilename = DaedalusGlobals.Orbit_Files_Path +"DAED_ORB_"
        OrbitFilename += Orbit_EventDropdown.value + "_"
        OrbitFilename += Orbit_DataFormatDropdown.value + "_"
        OrbitFilename += Orbit_PerigeeDropdown.value + "_"
        OrbitFilename += Orbit_LatitudeDropdown.value + "_"
        OrbitFilename += Orbit_SamplingDropdown.value + "_"
        OrbitFilename += Orbit_SpacecraftDropdown.value + ".nc"
        OrbitColorbarTitle = Orbit_ColorbartitleText.value
        OrbitColorscaleName = Orbit_ColorscaleText.value
        OrbitMultiplicationFactor = Orbit_MultiplyText.value
    ####
    VectorsVariablesToPlot  = Vectors_VarXDropdown.value[:Vectors_VarXDropdown.value.find('(')-1].strip() + ","
    VectorsVariablesToPlot += Vectors_VarYDropdown.value[:Vectors_VarYDropdown.value.find('(')-1].strip() + ","
    VectorsVariablesToPlot += Vectors_VarZDropdown.value[:Vectors_VarZDropdown.value.find('(')-1].strip()
    if( len(VectorsVariablesToPlot) > 2 ):
        VectorsFilename  = DaedalusGlobals.AllData_Files_Path + Vectors_TypeDropdown.value + "/" + Vectors_EventDropdown.value
        VectorsFilename += "/" + Vectors_DatafileDropdown.value 
        VectorsColorbarTitle = Vectors_ColorbartitleText.value
        VectorsColorscaleName = Vectors_ColorscaleText.value
        VectorConeSize = Vectors_ConesizeText.value
        VectorConeOpacity = Vectors_OpacityText.value
        VectorsMultiplicationFactor = Vectors_MultiplyText.value
        VectorsTimestep = Vectors_TimestepText.value
        VectorsPressureLevel = Vectors_PressurelevelText.value
    ####
    print(VectorsFilename)
    print(VectorsVariablesToPlot)
    print(VectorConeSize)
    print(VectorConeOpacity)
    print(VectorsTimestep)
    print(VectorsPressureLevel)
    PlotGlobe.PlotGlobe( SurfaceFilename, SurfaceVariableToPlot, SurfaceColorbarTitle, SurfaceColorscaleName, OrbitFilename, OrbitVariableToPlot, OrbitColorbarTitle, OrbitColorscaleName, PlotTitle, VectorsFilename, VectorsVariablesToPlot, VectorsColorbarTitle, VectorsColorscaleName, SurfaceOpacity, LogScale, VectorConeSize, VectorConeOpacity, SurfaceMultiplicationFactor, OrbitMultiplicationFactor, VectorsMultiplicationFactor, OneSizeVectorCones, SurfaceTimestep, SurfacePressureLevel, VectorsTimestep, VectorsPressureLevel )
################################################################################################################
######################################## GUI Creation ##########################################################
################################################################################################################        
        

WholeGUI = w.VBox() # the top level visual element
style1 = {'description_width':'100px'}
style2 = {'description_width':'60px'}
layout1 = {'width':'200px'}
layout2 = {'width':'280px'}

# ---------------- construct visual elements for the selection of an ORBIT file ----------------
Orbit_Column1 = w.VBox( layout={'margin':'0px 65px 0px 0px'} )
Orbit_Column2 = w.VBox( layout={'margin':'0px 65px 0px 0px'} )
Orbit_Column3 = w.VBox()
OrbitSelectionPanel = w.HBox()

Orbit_EventDropdown = w.Dropdown( options=['', 'Evt0', 'Evt1', 'Evt2', 'Evt3', 'Evt4', 'Evt5', 'Evt6', 'Evt7', 'Evt8', 'Evt9'], description='Event:', style=style1, layout=layout1)
Orbit_DataFormatDropdown = w.Dropdown( options=['LLA', 'VEL', 'PTG'], description='Data format:', style=style1, layout=layout1)
Orbit_PerigeeDropdown = w.Dropdown( options=['Per120', 'Per150'], description='Perigee Altitude:', style=style1, layout=layout1)
Orbit_LatitudeDropdown = w.Dropdown( options=['Lat00', 'Lat40', 'Lat80'],  description='Perigee Latitude:', style=style1, layout=layout1)
Orbit_SamplingDropdown = w.Dropdown( options=['Srt16Hz', 'Srt01Hz'],  description='Sampling rate:', style=style1, layout=layout1)
Orbit_SpacecraftDropdown = w.Dropdown( options=['Msc', 'Ssc'],  description='Spacecraft:', style=style1, layout=layout1)
Orbit_EventDropdown.observe( Orbit_Filename_onChange )
Orbit_DataFormatDropdown.observe( Orbit_Filename_onChange )
Orbit_PerigeeDropdown.observe( Orbit_Filename_onChange )
Orbit_LatitudeDropdown.observe( Orbit_Filename_onChange )
Orbit_SamplingDropdown.observe( Orbit_Filename_onChange )
Orbit_SpacecraftDropdown.observe( Orbit_Filename_onChange )
Orbit_Column1.children += (w.HTML(value="<div align='left'><b><font color='LightSeaGreen'>Select an orbit:</b></div>"),)
Orbit_Column1.children += (Orbit_EventDropdown,)
Orbit_Column1.children += (Orbit_DataFormatDropdown,)
Orbit_Column1.children += (Orbit_PerigeeDropdown,)
Orbit_Column1.children += (Orbit_LatitudeDropdown,)
Orbit_Column1.children += (Orbit_SamplingDropdown,)
Orbit_Column1.children += (Orbit_SpacecraftDropdown,)

Orbit_VariableDropdown = w.Dropdown( value='', options=[''], description='', style={'description_width':'50px'}, layout=layout1)
Orbit_VariableDropdown.observe( Orbit_Var_onChange )
Orbit_MultiplyText = w.FloatText(value=1, description='Multiply by', layout=layout1)
Orbit_Column2.children += (w.HTML(value="<div align='left'><b><font color='LightSeaGreen'>Select a variable:</b></div>"),)
Orbit_Column2.children += (Orbit_VariableDropdown,)
Orbit_Column2.children += (Orbit_MultiplyText,)

Orbit_ColorbartitleText = w.Text(value="", description='Units:', layout=layout1)
Orbit_ColorscaleText = w.Dropdown( options=ColorScales, value='jet', description='Colorscale:', layout=layout1 )
Orbit_Column3.children += (w.HTML(value="<div align='left'><b><font color='LightSeaGreen'>Select Properties:</b></div>"),)
Orbit_Column3.children += (Orbit_ColorbartitleText,)
Orbit_Column3.children += (Orbit_ColorscaleText,)

OrbitSelectionPanel.children = (Orbit_Column1, Orbit_Column2, Orbit_Column3,)

# ---------------- construct visual elements for the selection of a SURFACE file ----------------
SurfaceSelectionPanel = w.VBox()
Surface_Column1 = w.VBox( layout={'margin':'0px 65px 0px 0px'} )
Surface_Column2 = w.VBox( layout={'margin':'0px 65px 0px 0px'} )
Surface_Column3 = w.VBox( layout={'margin':'0px 65px 0px 0px'} )
SurfaceSelectionPanel = w.HBox()

Surface_TypeDropdown = w.Dropdown(options=['', 'TIEGCM_DATA', 'WACCAM_X'], description='File Type:', style=style2, layout=layout2)
Surface_TypeDropdown.observe( Surface_TypeDropdown_onChange )
Surface_EventDropdown = w.Dropdown(options=[''], description='Event:', style=style2, layout=layout2)
Surface_EventDropdown.observe( Surface_EventDropdown_onChange )
Surface_DatafileDropdown = w.Dropdown(options=[''], description='Data File:', style=style2, layout=layout2)
Surface_DatafileDropdown.observe( Surface_DatafileDropdown_onChange )
Surface_Column1.children += (w.HTML(value="<div align='left'><b><font color='LightSeaGreen'>Select a Data File:</b></div>"),)
Surface_Column1.children += (Surface_TypeDropdown,)
Surface_Column1.children += (Surface_EventDropdown,)
Surface_Column1.children += (Surface_DatafileDropdown,)

Surface_VariableDropdown = w.Dropdown( value='', options=[''], description='', style={'description_width':'50px'}, layout=layout1)
Surface_VariableDropdown.observe( Surface_Var_onChange )
Surface_TimestepText = w.BoundedIntText( value=0, min=0, description='Timestep:', layout=layout1)
Surface_PressurelevelText = w.BoundedIntText( value=0, min=0, description='Pres. Level:', layout=layout1)
Surface_MultiplyText = w.FloatText(value=1, description='Multiply by:', layout=layout1)
Surface_Column2.children += (w.HTML(value="<div align='left'><b><font color='LightSeaGreen'>Select a variable:</b></div>"),)
Surface_Column2.children += (Surface_VariableDropdown,)
Surface_Column2.children += (Surface_PressurelevelText,)
Surface_Column2.children += (Surface_TimestepText,)
Surface_Column2.children += (Surface_MultiplyText,)

Surface_ColorbartitleText = w.Text(value="", description='Units:', layout=layout1)
Surface_ColorscaleText = w.Dropdown( options=ColorScales, value='jet', description='Colorscale:', layout=layout1 )
Surface_OpacityText = w.BoundedFloatText(value=0.9, min=0, max=1, step=0.05, description='Opacity:', layout=layout1)
Surface_Column3.children += (w.HTML(value="<div align='left'><b><font color='LightSeaGreen'>Select Properties:</b></div>"),)
Surface_Column3.children += (Surface_ColorbartitleText,)
Surface_Column3.children += (Surface_ColorscaleText,)
Surface_Column3.children += (Surface_OpacityText,)

SurfaceSelectionPanel.children = (Surface_Column1, Surface_Column2, Surface_Column3,)

# ---------------- construct visual elements for the selection of a VECTORS file ----------------
VectorsSelectionPanel = w.VBox()
Vectors_Column1 = w.VBox( layout={'margin':'0px 65px 0px 0px'} )
Vectors_Column2 = w.VBox( layout={'margin':'0px 65px 0px 0px'} )
Vectors_Column3 = w.VBox( layout={'margin':'0px 65px 0px 0px'} )
VectorsSelectionPanel = w.HBox()

Vectors_TypeDropdown = w.Dropdown(options=['', 'TIEGCM_DATA', 'WACCAM_X'], description='File Type:', style=style2, layout=layout2)
Vectors_TypeDropdown.observe( Vectors_TypeDropdown_onChange )
Vectors_EventDropdown = w.Dropdown(options=[''], description='Event:', style=style2, layout=layout2)
Vectors_EventDropdown.observe( Vectors_EventDropdown_onChange )
Vectors_DatafileDropdown = w.Dropdown(options=[''], description='Data File:', style=style2, layout=layout2)
Vectors_DatafileDropdown.observe( Vectors_DatafileDropdown_onChange )
Vectors_Column1.children += (w.HTML(value="<div align='left'><b><font color='LightSeaGreen'>Select a Data File:</b></div>"),)
Vectors_Column1.children += (Vectors_TypeDropdown,)
Vectors_Column1.children += (Vectors_EventDropdown,)
Vectors_Column1.children += (Vectors_DatafileDropdown,)

Vectors_VarXDropdown = w.Dropdown( value='', options=[''], description='X:', style={'description_width':'20px'}, layout=layout1)
Vectors_VarYDropdown = w.Dropdown( value='', options=[''], description='Y:', style={'description_width':'20px'}, layout=layout1)
Vectors_VarZDropdown = w.Dropdown( value='', options=[''], description='Z:', style={'description_width':'20px'}, layout=layout1)
Vectors_VarXDropdown.observe( Vectors_Var_onChange )
Vectors_VarYDropdown.observe( Vectors_Var_onChange )
Vectors_VarZDropdown.observe( Vectors_Var_onChange )
Vectors_TimestepText = w.BoundedIntText( value=0, min=0, description='Timestep:', layout=layout1)
Vectors_PressurelevelText = w.BoundedIntText( value=0, min=0, description='Pres. Level:', layout=layout1)
Vectors_MultiplyText = w.FloatText(value=1, description='Multiply by:', layout=layout1)
Vectors_Column2.children += (w.HTML(value="<div align='left'><b><font color='LightSeaGreen'>Select a variable:</b></div>"),)
Vectors_Column2.children += (Vectors_VarXDropdown,Vectors_VarYDropdown,Vectors_VarZDropdown,)
Vectors_Column2.children += (Vectors_PressurelevelText,)
Vectors_Column2.children += (Vectors_TimestepText,)
Vectors_Column2.children += (Vectors_MultiplyText,)

Vectors_ColorbartitleText = w.Text(value="", description='Units:', layout=layout1)
Vectors_ColorscaleText = w.Dropdown( options=ColorScales, value='jet', description='Colorscale:', layout=layout1 )
Vectors_OpacityText = w.BoundedFloatText(value=0.9, min=0, max=1, step=0.05, description='Opacity:', layout=layout1)
Vectors_ConesizeText = w.FloatText(value=44, description='Cone size:', layout=layout1)
Vectors_Column3.children += (w.HTML(value="<div align='left'><b><font color='LightSeaGreen'>Select Properties:</b></div>"),)
Vectors_Column3.children += (Vectors_ColorbartitleText,)
Vectors_Column3.children += (Vectors_ColorscaleText,)
Vectors_Column3.children += (Vectors_OpacityText,)
Vectors_Column3.children += (Vectors_ConesizeText,)

VectorsSelectionPanel.children = (Vectors_Column1, Vectors_Column2, Vectors_Column3,)

# ---------------- construct visual elements for the selection of PLOT PREFERENCES ----------------
PreferencesPanel = w.VBox()

Preferences_Plottitle = w.Text(value="", description='Plot Title:', layout=layout1)
Preferences_Logscale = w.Checkbox( value=False, description='Log Scale', layout=layout1)
PreferencesPanel.children += (Preferences_Plottitle,)
PreferencesPanel.children += (Preferences_Logscale,)




# ---------------- construct the top level visual elements ----------------
OrbitAccordion = w.Accordion( children=[OrbitSelectionPanel], selected_index=None )
OrbitAccordion.set_title(0, "Click here to plot an Orbit" )
SurfaceAccordion = w.Accordion( children=[SurfaceSelectionPanel], selected_index=None )
SurfaceAccordion.set_title(0, 'Click here to plot a variable onto a Surface above the Earth')
VectorsAccordion = w.Accordion( children=[VectorsSelectionPanel], selected_index=None )
VectorsAccordion.set_title(0, 'Click here to plot a Vector field')
PreferencesAccordion = w.Accordion( children=[PreferencesPanel], selected_index=None )
PreferencesAccordion.set_title(0, 'Click here to select the Plot Preferences')

Plot_Btn = w.Button (description='P l o t',tooltip="Click here to plot the above selected data",)
Plot_Btn.layout.min_width = '140px'
Plot_Btn.layout.min_height = '44px'
Plot_Btn.style.font_weight = 'bold'
Plot_Btn.style.button_color = 'MediumAquamarine'
Plot_Btn.on_click(plot_button_clicked)

WholeGUI.children += (OrbitAccordion,)
WholeGUI.children += (SurfaceAccordion,)
WholeGUI.children += (VectorsAccordion,)
WholeGUI.children += (PreferencesAccordion,)
WholeGUI.children += (Plot_Btn,)


def display():
    return WholeGUI






