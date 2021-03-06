import ipywidgets as widgets


def Construct_IGRF():
	#Create Containers
	IGRF_Panel = widgets.Box()
	IGRF_Panel.layout.overflow_x = 'scroll'
	IGRF_Panel.layout.flex_flow = 'row'
	IGRF_Panel.layout.display = 'flex'
	########
	## GUI code for module 'OrbitSelector'
	########
	# Create widgets for module's inputs
	InputsPanel = widgets.VBox()
	InputsPanel.layout.min_width = '330px'
	global WIDGET_OrbitSelector_Filename
	WIDGET_OrbitSelector_Filename = widgets.Text()
	WIDGET_OrbitSelector_Filename.description = 'Filename'
	WIDGET_OrbitSelector_Filename.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_OrbitSelector_Filename,)
	global WIDGET_OrbitSelector_EvtXY
	WIDGET_OrbitSelector_EvtXY = widgets.Dropdown( options=['', 'Evt0', 'Evt1', 'Evt2', 'Evt3', 'Evt4', 'Evt5', 'Evt6', 'Evt7', 'Evt8', 'Evt9'], description='Variable')
	WIDGET_OrbitSelector_EvtXY.description = 'EvtXY'
	WIDGET_OrbitSelector_EvtXY.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_OrbitSelector_EvtXY,)
	global WIDGET_OrbitSelector_TYP
	WIDGET_OrbitSelector_TYP = widgets.Dropdown( options=['LLA', 'VEL', 'PTG'], description='Variable')
	WIDGET_OrbitSelector_TYP.description = 'TYP'
	WIDGET_OrbitSelector_TYP.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_OrbitSelector_TYP,)
	global WIDGET_OrbitSelector_PerYYY
	WIDGET_OrbitSelector_PerYYY = widgets.Dropdown( options=['Per120', 'Per150'], description='Variable')
	WIDGET_OrbitSelector_PerYYY.description = 'PerYYY'
	WIDGET_OrbitSelector_PerYYY.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_OrbitSelector_PerYYY,)
	global WIDGET_OrbitSelector_LatZZ
	WIDGET_OrbitSelector_LatZZ = widgets.Dropdown( options=['Lat00', 'Lat40', 'Lat80'],  description='Variable')
	WIDGET_OrbitSelector_LatZZ.description = 'LatZZ'
	WIDGET_OrbitSelector_LatZZ.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_OrbitSelector_LatZZ,)
	global WIDGET_OrbitSelector_SRXXHZ
	WIDGET_OrbitSelector_SRXXHZ = widgets.Dropdown( options=['Srt16Hz', 'Srt01Hz'],  description='Variable')
	WIDGET_OrbitSelector_SRXXHZ.description = 'SRXXHZ'
	WIDGET_OrbitSelector_SRXXHZ.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_OrbitSelector_SRXXHZ,)
	global WIDGET_OrbitSelector_SC
	WIDGET_OrbitSelector_SC = widgets.Dropdown( options=['Msc', 'Ssc'],  description='Variable')
	WIDGET_OrbitSelector_SC.description = 'SC'
	WIDGET_OrbitSelector_SC.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_OrbitSelector_SC,)
	IGRF_Panel.children += (InputsPanel,)
	# Create widget for the moudle black-box body
	OrbitSelector_Btn = widgets.Button (
		description='OrbitSelector',
		tooltip="Constructs and returns an orbit full-path filename from the orbit properties.\nThe filename includes all parameters required to select the appropriate orbit.\nAll filenames  have the same length and same number of parameters. \nIn case OrbitFilename and EvtXY are empty then an empty string is returned.\nFILE FORMAT: DAED_ORB_EVTXY_TYP_PERYYY_LATZZ_SRTQQHz_XSC.csv\n  Parameter OrbitFilename:\n    If it is not empty then the rest arguments are ignored. \n    If it contains slashes then is is assumed it is a full path name and it is returned as it is.\n    If it does not contain slashes then the orbits path is added.\n    \n  Parameter EvtXY values:\n    EVTXS    X Event, Single Orbit\n    EVTXA    X Event, All Orbit\n    EVT1Y    1st Event: St Patrick’s day event [17 Mar 2015 – 20 Mar 2015]\n    EVT2Y    2nd Event ...  ...\n \n  Parameter TYP values:\n    LLA    Time,Latitude Longitude Altitude \n    VEL    Time,VmagVxVyVz\n    PTG    Time, X_GSE, Y_GSE, Z_GSE, RamX_GSE, RamY_GSE, RamZ_GSE\n  Parameter PerYYY values:\n    PER120, PER150    Perigee Altitude at 120km or 150km\n  Parameter LatZZ values:\n    LAT00, LAT40, LAT80    Perigee Latitude at 0°, 40° or 80°\n \n  Parameter SRXXHZ values:\n    SRT16Hz, SRT01Hz    Sampling rate 16Hz or 1HZ\n \n  Parameter SC values:\n    MSC    Mother Spacecraft\n    SSC    Sub-Spacecraft\n",
	)
	OrbitSelector_Btn.layout.min_width = '200px'
	OrbitSelector_Btn.style.button_color = 'YellowGreen'
	IGRF_Panel.children += (OrbitSelector_Btn,)
	# Create widgets for module's outputs
	OutputsPanel = widgets.VBox()
	OutputsPanel.layout.min_width = '300px'
	WIDGET_OrbitSelector_OrbitCSVfilename = widgets.Label(value='  --> OrbitCSVfilename  ')
	WIDGET_OrbitSelector_OrbitCSVfilename.layout.border = '1px dashed green'
	WIDGET_OrbitSelector_OrbitCSVfilename.layout.margin ='0px 40px 0px 0px' 
	WIDGET_OrbitSelector_OrbitCSVfilename.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_OrbitSelector_OrbitCSVfilename,)
	WIDGET_OrbitSelector_OrbitNetCDFfilename = widgets.Label(value='  --> OrbitNetCDFfilename  ')
	WIDGET_OrbitSelector_OrbitNetCDFfilename.layout.border = '1px dashed green'
	WIDGET_OrbitSelector_OrbitNetCDFfilename.layout.margin ='0px 40px 0px 0px' 
	WIDGET_OrbitSelector_OrbitNetCDFfilename.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_OrbitSelector_OrbitNetCDFfilename,)
	IGRF_Panel.children += (OutputsPanel,)
	########
	## GUI code for module 'IGRF'
	########
	# Create widgets for module's inputs
	InputsPanel = widgets.VBox()
	InputsPanel.layout.min_width = '330px'
	WIDGET_IGRF_CSVfilename_forOrbit = widgets.Label(value='  OrbitSelector.OrbitCSVfilename --> CSVfilename_forOrbit  ')
	WIDGET_IGRF_CSVfilename_forOrbit.layout.border = '1px dashed blue'
	WIDGET_IGRF_CSVfilename_forOrbit.layout.padding = '0px 10px 0px 10px'
	InputsPanel.children += (WIDGET_IGRF_CSVfilename_forOrbit,)
	global WIDGET_IGRF_Variable_forOrbit
	WIDGET_IGRF_Variable_forOrbit = widgets.Dropdown( options=['', 'B', 'Bx', 'By', 'Bz'], value='', description='Variable')
	WIDGET_IGRF_Variable_forOrbit.description = 'Variable_forOrbit'
	WIDGET_IGRF_Variable_forOrbit.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_IGRF_Variable_forOrbit,)
	IGRF_Panel.children += (InputsPanel,)
	# Create widget for the moudle black-box body
	IGRF_Btn = widgets.Button (
		description='IGRF',
		tooltip="International Geomagnetic Reference Field empirical model (IGRF)\nfilename: the CSV file to read. Format: Time,Alt,Lon,Lat. The file defines the coordinates for which the IGRF calculation will take place.\nParameter: selects which values to calculate. Possible values are:\n    ''   for all\n    'B'  for magnetic field\n    'Bx' for magnetic field component at x axis\n    'By' for magnetic field component at y axis\n    'Bz' for magnetic field component at z axis\nReturns: the name of the CSV file which contains the calculated values. Format: Time,Alt,Lon,Lat,B,Bx,By,Bz",
	)
	IGRF_Btn.layout.min_width = '200px'
	IGRF_Btn.style.button_color = 'gold'
	IGRF_Panel.children += (IGRF_Btn,)
	# Create widgets for module's outputs
	OutputsPanel = widgets.VBox()
	OutputsPanel.layout.min_width = '300px'
	WIDGET_IGRF_OrbitResultCSV = widgets.Label(value='  --> OrbitResultCSV  ')
	WIDGET_IGRF_OrbitResultCSV.layout.border = '1px dashed green'
	WIDGET_IGRF_OrbitResultCSV.layout.margin ='0px 40px 0px 0px' 
	WIDGET_IGRF_OrbitResultCSV.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_IGRF_OrbitResultCSV,)
	IGRF_Panel.children += (OutputsPanel,)
	########
	## GUI code for module 'CreateCSV_Sphere'
	########
	# Create widgets for module's inputs
	InputsPanel = widgets.VBox()
	InputsPanel.layout.min_width = '330px'
	global WIDGET_CreateCSV_Sphere_CSVfilename
	WIDGET_CreateCSV_Sphere_CSVfilename = widgets.Text(value="")
	WIDGET_CreateCSV_Sphere_CSVfilename.description = 'CSVfilename'
	WIDGET_CreateCSV_Sphere_CSVfilename.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_CreateCSV_Sphere_CSVfilename,)
	global WIDGET_CreateCSV_Sphere_fixedDatetimeString
	WIDGET_CreateCSV_Sphere_fixedDatetimeString = widgets.Text(value="Mar 17 2015 00:50:00.000000000")
	WIDGET_CreateCSV_Sphere_fixedDatetimeString.description = 'fixedDatetimeString'
	WIDGET_CreateCSV_Sphere_fixedDatetimeString.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_CreateCSV_Sphere_fixedDatetimeString,)
	global WIDGET_CreateCSV_Sphere_fixedAltitude
	WIDGET_CreateCSV_Sphere_fixedAltitude = widgets.FloatText(value=120)
	WIDGET_CreateCSV_Sphere_fixedAltitude.description = 'fixedAltitude'
	WIDGET_CreateCSV_Sphere_fixedAltitude.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_CreateCSV_Sphere_fixedAltitude,)
	global WIDGET_CreateCSV_Sphere_LatitudeStep
	WIDGET_CreateCSV_Sphere_LatitudeStep = widgets.FloatText(value=5)
	WIDGET_CreateCSV_Sphere_LatitudeStep.description = 'LatitudeStep'
	WIDGET_CreateCSV_Sphere_LatitudeStep.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_CreateCSV_Sphere_LatitudeStep,)
	global WIDGET_CreateCSV_Sphere_LongitudeStep
	WIDGET_CreateCSV_Sphere_LongitudeStep = widgets.FloatText(value=5)
	WIDGET_CreateCSV_Sphere_LongitudeStep.description = 'LongitudeStep'
	WIDGET_CreateCSV_Sphere_LongitudeStep.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_CreateCSV_Sphere_LongitudeStep,)
	IGRF_Panel.children += (InputsPanel,)
	# Create widget for the moudle black-box body
	CreateCSV_Sphere_Btn = widgets.Button (
		description='CreateCSV_Sphere',
		tooltip="Creates a CSV file with the format [Time, Lat, Lon, Alt]\nTime and Alt are fixed values given as arguments and Lat and Lon are produced for the whole globe accodring to steps specified.\nCSVfilename: string, the csv filename to be written (overwrite mode)\nfixedDatetimeString: string, the date-time fixed text with format as example: Mar 17 2015 00:12:00.000. Fixed value for every csv entry.\nfixedAltitude: float positive, kilometers above earth. Fixed value for every csv entry.\nLatitudeStep: float positive, Latiude values range between 87.5 and -88.5. The step is how Lat should be incremented at each iteration.\nLongitudeStep: float positive, Longitude values range between -180 and 180.  The step is how Lon should be incremented at each iteration.\nRETURNS: the filename it has written (same as the argument)\n         the altitude (same as the argument)\n         the number of lines written",
	)
	CreateCSV_Sphere_Btn.layout.min_width = '200px'
	CreateCSV_Sphere_Btn.style.button_color = 'gold'
	IGRF_Panel.children += (CreateCSV_Sphere_Btn,)
	# Create widgets for module's outputs
	OutputsPanel = widgets.VBox()
	OutputsPanel.layout.min_width = '300px'
	WIDGET_CreateCSV_Sphere_theCSVfilename = widgets.Label(value='  --> theCSVfilename  ')
	WIDGET_CreateCSV_Sphere_theCSVfilename.layout.border = '1px dashed green'
	WIDGET_CreateCSV_Sphere_theCSVfilename.layout.margin ='0px 40px 0px 0px' 
	WIDGET_CreateCSV_Sphere_theCSVfilename.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_CreateCSV_Sphere_theCSVfilename,)
	WIDGET_CreateCSV_Sphere_theAltitude = widgets.Label(value='  --> theAltitude  ')
	WIDGET_CreateCSV_Sphere_theAltitude.layout.border = '1px dashed green'
	WIDGET_CreateCSV_Sphere_theAltitude.layout.margin ='0px 40px 0px 0px' 
	WIDGET_CreateCSV_Sphere_theAltitude.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_CreateCSV_Sphere_theAltitude,)
	WIDGET_CreateCSV_Sphere_NumberOfLinesWritten = widgets.Label(value='  --> NumberOfLinesWritten  ')
	WIDGET_CreateCSV_Sphere_NumberOfLinesWritten.layout.border = '1px dashed green'
	WIDGET_CreateCSV_Sphere_NumberOfLinesWritten.layout.margin ='0px 40px 0px 0px' 
	WIDGET_CreateCSV_Sphere_NumberOfLinesWritten.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_CreateCSV_Sphere_NumberOfLinesWritten,)
	IGRF_Panel.children += (OutputsPanel,)
	########
	## GUI code for module 'IGRF'
	########
	# Create widgets for module's inputs
	InputsPanel = widgets.VBox()
	InputsPanel.layout.min_width = '330px'
	WIDGET_IGRF_CSVfilename_forMeshgrid = widgets.Label(value='  CreateCSV_Sphere.theCSVfilename --> CSVfilename_forMeshgrid  ')
	WIDGET_IGRF_CSVfilename_forMeshgrid.layout.border = '1px dashed blue'
	WIDGET_IGRF_CSVfilename_forMeshgrid.layout.padding = '0px 10px 0px 10px'
	InputsPanel.children += (WIDGET_IGRF_CSVfilename_forMeshgrid,)
	global WIDGET_IGRF_Variable_forMeshgrid
	WIDGET_IGRF_Variable_forMeshgrid = widgets.Dropdown( options=['', 'B', 'Bx', 'By', 'Bz'], value='', description='Variable')
	WIDGET_IGRF_Variable_forMeshgrid.description = 'Variable_forMeshgrid'
	WIDGET_IGRF_Variable_forMeshgrid.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_IGRF_Variable_forMeshgrid,)
	IGRF_Panel.children += (InputsPanel,)
	# Create widget for the moudle black-box body
	IGRF_Btn = widgets.Button (
		description='IGRF',
		tooltip="International Geomagnetic Reference Field empirical model (IGRF)\nfilename: the CSV file to read. Format: Time,Alt,Lon,Lat. The file defines the coordinates for which the IGRF calculation will take place.\nParameter: selects which values to calculate. Possible values are:\n    ''   for all\n    'B'  for magnetic field\n    'Bx' for magnetic field component at x axis\n    'By' for magnetic field component at y axis\n    'Bz' for magnetic field component at z axis\nReturns: the name of the CSV file which contains the calculated values. Format: Time,Alt,Lon,Lat,B,Bx,By,Bz",
	)
	IGRF_Btn.layout.min_width = '200px'
	IGRF_Btn.style.button_color = 'gold'
	IGRF_Panel.children += (IGRF_Btn,)
	# Create widgets for module's outputs
	OutputsPanel = widgets.VBox()
	OutputsPanel.layout.min_width = '300px'
	WIDGET_IGRF_MeshgridResultCSV = widgets.Label(value='  --> MeshgridResultCSV  ')
	WIDGET_IGRF_MeshgridResultCSV.layout.border = '1px dashed green'
	WIDGET_IGRF_MeshgridResultCSV.layout.margin ='0px 40px 0px 0px' 
	WIDGET_IGRF_MeshgridResultCSV.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_IGRF_MeshgridResultCSV,)
	IGRF_Panel.children += (OutputsPanel,)
	########
	## GUI code for module 'PlotGlobe'
	########
	# Create widgets for module's inputs
	InputsPanel = widgets.VBox()
	InputsPanel.layout.min_width = '330px'
	WIDGET_PlotGlobe_SurfaceFilename = widgets.Label(value='  IGRF.MeshgridResultCSV --> SurfaceFilename  ')
	WIDGET_PlotGlobe_SurfaceFilename.layout.border = '1px dashed blue'
	WIDGET_PlotGlobe_SurfaceFilename.layout.padding = '0px 10px 0px 10px'
	InputsPanel.children += (WIDGET_PlotGlobe_SurfaceFilename,)
	global WIDGET_PlotGlobe_SurfaceVariableToPlot
	WIDGET_PlotGlobe_SurfaceVariableToPlot = widgets.Text(value='ALFA')
	WIDGET_PlotGlobe_SurfaceVariableToPlot.description = 'SurfaceVariableToPlot'
	WIDGET_PlotGlobe_SurfaceVariableToPlot.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_PlotGlobe_SurfaceVariableToPlot,)
	global WIDGET_PlotGlobe_SurfaceColorbarTitle
	WIDGET_PlotGlobe_SurfaceColorbarTitle = widgets.Text()
	WIDGET_PlotGlobe_SurfaceColorbarTitle.description = 'SurfaceColorbarTitle'
	WIDGET_PlotGlobe_SurfaceColorbarTitle.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_PlotGlobe_SurfaceColorbarTitle,)
	global WIDGET_PlotGlobe_SurfaceColorscaleName
	WIDGET_PlotGlobe_SurfaceColorscaleName = widgets.Dropdown( options=['', 'None', 'Blackbody', 'Bluered', 'Blues', 'Earth', 'Electric', 'Greens', 'Greys', 'Hot', 'Jet', 'Picnic', 'Portland', 'Rainbow', 'RdBu', 'Reds', 'Viridis', 'YlGnBu', 'YlOrRd'], value='Jet', description='Variable')
	WIDGET_PlotGlobe_SurfaceColorscaleName.description = 'SurfaceColorscaleName'
	WIDGET_PlotGlobe_SurfaceColorscaleName.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_PlotGlobe_SurfaceColorscaleName,)
	WIDGET_PlotGlobe_OrbitFilename = widgets.Label(value='  IGRF.OrbitResultCSV --> OrbitFilename  ')
	WIDGET_PlotGlobe_OrbitFilename.layout.border = '1px dashed blue'
	WIDGET_PlotGlobe_OrbitFilename.layout.padding = '0px 10px 0px 10px'
	InputsPanel.children += (WIDGET_PlotGlobe_OrbitFilename,)
	global WIDGET_PlotGlobe_OrbitVariableToPlot
	WIDGET_PlotGlobe_OrbitVariableToPlot = widgets.Text(value='ALFA')
	WIDGET_PlotGlobe_OrbitVariableToPlot.description = 'OrbitVariableToPlot'
	WIDGET_PlotGlobe_OrbitVariableToPlot.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_PlotGlobe_OrbitVariableToPlot,)
	global WIDGET_PlotGlobe_OrbitColorbarTitle
	WIDGET_PlotGlobe_OrbitColorbarTitle = widgets.Text()
	WIDGET_PlotGlobe_OrbitColorbarTitle.description = 'OrbitColorbarTitle'
	WIDGET_PlotGlobe_OrbitColorbarTitle.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_PlotGlobe_OrbitColorbarTitle,)
	global WIDGET_PlotGlobe_OrbitColorscaleName
	WIDGET_PlotGlobe_OrbitColorscaleName = widgets.Dropdown( options=['', 'None', 'Blackbody', 'Bluered', 'Blues', 'Earth', 'Electric', 'Greens', 'Greys', 'Hot', 'Jet', 'Picnic', 'Portland', 'Rainbow', 'RdBu', 'Reds', 'Viridis', 'YlGnBu', 'YlOrRd'], value='Jet', description='Variable')
	WIDGET_PlotGlobe_OrbitColorscaleName.description = 'OrbitColorscaleName'
	WIDGET_PlotGlobe_OrbitColorscaleName.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_PlotGlobe_OrbitColorscaleName,)
	global WIDGET_PlotGlobe_PlotTitle
	WIDGET_PlotGlobe_PlotTitle = widgets.Text()
	WIDGET_PlotGlobe_PlotTitle.description = 'PlotTitle'
	WIDGET_PlotGlobe_PlotTitle.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_PlotGlobe_PlotTitle,)
	IGRF_Panel.children += (InputsPanel,)
	# Create widget for the moudle black-box body
	PlotGlobe_Btn = widgets.Button (
		description='PlotGlobe',
		tooltip="Creates a 3D plot of an earth globe. Can plot a sphere surface and/or a satellite orbit. \nThe surface and the orbit are colored according to the selected variable values in the data files. \nThe data files can be of CSV o NetCDF4 format.\nThe same or different variables can be selected for plotting at the surface and the orbit.\nThe same or different colorscales can be selected for the surface and the orbit.\nThe same or different value range can be selected for the surface and the orbit by assigning the same or different colorbar titles.\nValid Colorscale Names: ‘Blackbody’, ‘Bluered’, ‘Blues’, ‘Earth’, ‘Electric’, ‘Greens’, ‘Greys’, ‘Hot’, ‘Jet’, ‘Picnic’, ‘Portland’, ‘Rainbow’, ‘RdBu’, ‘Reds’, ‘Viridis’, ‘YlGnBu’, ‘YlOrRd’\nARGUMENTS:\n  SurfaceFilename: The file which contains data for a surface. If empty string then no surface will be plotted. \n                   CSV format: Time,Lat,Lon,Alt,value. Contains the data for the sphere surface. \n                   NetCDF format: see Daedalus User Manual\n  SurfaceVariableToPlot: The name of the variable to read from the surface data file and plot upon the globe.\n  SurfaceColorbarTitle: A title to display above the colorbar which refers to the surface data\n  SurfaceColorscaleName: The name of the Colorscale to use for the surface data. \n                         In case an empty string is passed then a default HeatMap colorscale will be applied.\n                         In case None is passed then all points will be black irrespective of value.\n  OrbitFilename: counterpart of the corresponding argument for the Surface data\n  OrbitVariableToPlot: counterpart of the corresponding argument for the Surface data \n  OrbitColorbarTitle: counterpart of the corresponding argument for the Surface data \n  OrbitColorscaleName: counterpart of the corresponding argument for the Surface data\n  PlotTitle: A title which is displayed on the top of the plot.\nRETURNS: a string containing information about the Data\n",
	)
	PlotGlobe_Btn.layout.min_width = '200px'
	PlotGlobe_Btn.style.button_color = 'deeppink'
	IGRF_Panel.children += (PlotGlobe_Btn,)
	# Create widgets for module's outputs
	OutputsPanel = widgets.VBox()
	OutputsPanel.layout.min_width = '300px'
	IGRF_Panel.children += (OutputsPanel,)
	return IGRF_Panel


def Construct_Interpolator_NetCDF():
	#Create Containers
	Interpolator_NetCDF_Panel = widgets.Box()
	Interpolator_NetCDF_Panel.layout.overflow_x = 'scroll'
	Interpolator_NetCDF_Panel.layout.flex_flow = 'row'
	Interpolator_NetCDF_Panel.layout.display = 'flex'
	########
	## GUI code for module 'OrbitSelector'
	########
	# Create widgets for module's inputs
	InputsPanel = widgets.VBox()
	InputsPanel.layout.min_width = '330px'
	global WIDGET_OrbitSelector_Filename
	WIDGET_OrbitSelector_Filename = widgets.Text("")
	WIDGET_OrbitSelector_Filename.description = 'Filename'
	WIDGET_OrbitSelector_Filename.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_OrbitSelector_Filename,)
	global WIDGET_OrbitSelector_EvtXY
	WIDGET_OrbitSelector_EvtXY = widgets.Dropdown( options=['', 'Evt0', 'Evt1', 'Evt2', 'Evt3', 'Evt4', 'Evt5', 'Evt6', 'Evt7', 'Evt8', 'Evt9'], value='Evt0', description='Variable')
	WIDGET_OrbitSelector_EvtXY.description = 'EvtXY'
	WIDGET_OrbitSelector_EvtXY.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_OrbitSelector_EvtXY,)
	global WIDGET_OrbitSelector_TYP
	WIDGET_OrbitSelector_TYP = widgets.Dropdown( options=['LLA', 'VEL', 'PTG'], description='Variable')
	WIDGET_OrbitSelector_TYP.description = 'TYP'
	WIDGET_OrbitSelector_TYP.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_OrbitSelector_TYP,)
	global WIDGET_OrbitSelector_PerYYY
	WIDGET_OrbitSelector_PerYYY = widgets.Dropdown( options=['Per120', 'Per150'], description='Variable')
	WIDGET_OrbitSelector_PerYYY.description = 'PerYYY'
	WIDGET_OrbitSelector_PerYYY.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_OrbitSelector_PerYYY,)
	global WIDGET_OrbitSelector_LatZZ
	WIDGET_OrbitSelector_LatZZ = widgets.Dropdown( options=['Lat00', 'Lat40', 'Lat80'],  description='Variable')
	WIDGET_OrbitSelector_LatZZ.description = 'LatZZ'
	WIDGET_OrbitSelector_LatZZ.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_OrbitSelector_LatZZ,)
	global WIDGET_OrbitSelector_SRXXHZ
	WIDGET_OrbitSelector_SRXXHZ = widgets.Dropdown( options=['Srt16Hz', 'Srt01Hz'],  value='Srt01Hz', description='Variable')
	WIDGET_OrbitSelector_SRXXHZ.description = 'SRXXHZ'
	WIDGET_OrbitSelector_SRXXHZ.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_OrbitSelector_SRXXHZ,)
	global WIDGET_OrbitSelector_SC
	WIDGET_OrbitSelector_SC = widgets.Dropdown( options=['Msc', 'Ssc'],  description='Variable')
	WIDGET_OrbitSelector_SC.description = 'SC'
	WIDGET_OrbitSelector_SC.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_OrbitSelector_SC,)
	Interpolator_NetCDF_Panel.children += (InputsPanel,)
	# Create widget for the moudle black-box body
	OrbitSelector_Btn = widgets.Button (
		description='OrbitSelector',
		tooltip="Constructs and returns an orbit full-path filename from the orbit properties.\nThe filename includes all parameters required to select the appropriate orbit.\nAll filenames  have the same length and same number of parameters. \nIn case OrbitFilename and EvtXY are empty then an empty string is returned.\nFILE FORMAT: DAED_ORB_EVTXY_TYP_PERYYY_LATZZ_SRTQQHz_XSC.csv\n  Parameter OrbitFilename:\n    If it is not empty then the rest arguments are ignored. \n    If it contains slashes then is is assumed it is a full path name and it is returned as it is.\n    If it does not contain slashes then the orbits path is added.\n    \n  Parameter EvtXY values:\n    EVTXS    X Event, Single Orbit\n    EVTXA    X Event, All Orbit\n    EVT1Y    1st Event: St Patrick’s day event [17 Mar 2015 – 20 Mar 2015]\n    EVT2Y    2nd Event ...  ...\n \n  Parameter TYP values:\n    LLA    Time,Latitude Longitude Altitude \n    VEL    Time,VmagVxVyVz\n    PTG    Time, X_GSE, Y_GSE, Z_GSE, RamX_GSE, RamY_GSE, RamZ_GSE\n  Parameter PerYYY values:\n    PER120, PER150    Perigee Altitude at 120km or 150km\n  Parameter LatZZ values:\n    LAT00, LAT40, LAT80    Perigee Latitude at 0°, 40° or 80°\n \n  Parameter SRXXHZ values:\n    SRT16Hz, SRT01Hz    Sampling rate 16Hz or 1HZ\n \n  Parameter SC values:\n    MSC    Mother Spacecraft\n    SSC    Sub-Spacecraft\n",
	)
	OrbitSelector_Btn.layout.min_width = '200px'
	OrbitSelector_Btn.style.button_color = 'YellowGreen'
	Interpolator_NetCDF_Panel.children += (OrbitSelector_Btn,)
	# Create widgets for module's outputs
	OutputsPanel = widgets.VBox()
	OutputsPanel.layout.min_width = '300px'
	WIDGET_OrbitSelector_OrbitCSVfilename = widgets.Label(value='  --> OrbitCSVfilename  ')
	WIDGET_OrbitSelector_OrbitCSVfilename.layout.border = '1px dashed green'
	WIDGET_OrbitSelector_OrbitCSVfilename.layout.margin ='0px 40px 0px 0px' 
	WIDGET_OrbitSelector_OrbitCSVfilename.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_OrbitSelector_OrbitCSVfilename,)
	WIDGET_OrbitSelector_OrbitNetCDFfilename = widgets.Label(value='  --> OrbitNetCDFfilename  ')
	WIDGET_OrbitSelector_OrbitNetCDFfilename.layout.border = '1px dashed green'
	WIDGET_OrbitSelector_OrbitNetCDFfilename.layout.margin ='0px 40px 0px 0px' 
	WIDGET_OrbitSelector_OrbitNetCDFfilename.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_OrbitSelector_OrbitNetCDFfilename,)
	Interpolator_NetCDF_Panel.children += (OutputsPanel,)
	########
	## GUI code for module 'Interpolator'
	########
	# Create widgets for module's inputs
	InputsPanel = widgets.VBox()
	InputsPanel.layout.min_width = '330px'
	global WIDGET_Interpolator_ModelFile_forOrbit
	WIDGET_Interpolator_ModelFile_forOrbit = widgets.Text(value="/home/NAS/TIEGCM_DATA/TIEGCM_EVT2_2017_Sept_S/tiegcm_s_24900.nc")
	WIDGET_Interpolator_ModelFile_forOrbit.description = 'ModelFile_forOrbit'
	WIDGET_Interpolator_ModelFile_forOrbit.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_Interpolator_ModelFile_forOrbit,)
	WIDGET_Interpolator_OrbitFilename = widgets.Label(value='  OrbitSelector.OrbitNetCDFfilename --> OrbitFilename  ')
	WIDGET_Interpolator_OrbitFilename.layout.border = '1px dashed blue'
	WIDGET_Interpolator_OrbitFilename.layout.padding = '0px 10px 0px 10px'
	InputsPanel.children += (WIDGET_Interpolator_OrbitFilename,)
	Interpolator_NetCDF_Panel.children += (InputsPanel,)
	# Create widget for the moudle black-box body
	Interpolator_Btn = widgets.Button (
		description='Interpolator',
		tooltip="model-->string:: name of model eg TIEGCM\n model_data_file--> string:: model data files stored on DaedalusNAS\n orbit_file-->string :: orbit filename in Time-Lat-Lon-Alt format stored on DaedalusNAS\n save--> Logical:: if true saves interpolated values to directory\n VAR--> string:: variable to inteprolate, must be  one of the variables included in the model data\nOutputs:: Plots+ 1 csv file stored on DaedalusNAS in ModelOutputs/Interpolator",
	)
	Interpolator_Btn.layout.min_width = '200px'
	Interpolator_Btn.style.button_color = 'gold'
	Interpolator_NetCDF_Panel.children += (Interpolator_Btn,)
	# Create widgets for module's outputs
	OutputsPanel = widgets.VBox()
	OutputsPanel.layout.min_width = '300px'
	WIDGET_Interpolator_InterpolatedOrbit = widgets.Label(value='  --> InterpolatedOrbit  ')
	WIDGET_Interpolator_InterpolatedOrbit.layout.border = '1px dashed green'
	WIDGET_Interpolator_InterpolatedOrbit.layout.margin ='0px 40px 0px 0px' 
	WIDGET_Interpolator_InterpolatedOrbit.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_Interpolator_InterpolatedOrbit,)
	Interpolator_NetCDF_Panel.children += (OutputsPanel,)
	########
	## GUI code for module 'CreateNETCDFsphere'
	########
	# Create widgets for module's inputs
	InputsPanel = widgets.VBox()
	InputsPanel.layout.min_width = '330px'
	global WIDGET_CreateNETCDFsphere_fixedDatetimeString
	WIDGET_CreateNETCDFsphere_fixedDatetimeString = widgets.Text(value="Mar 17 2015 00:50:00.000000000")
	WIDGET_CreateNETCDFsphere_fixedDatetimeString.description = 'fixedDatetimeString'
	WIDGET_CreateNETCDFsphere_fixedDatetimeString.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_CreateNETCDFsphere_fixedDatetimeString,)
	global WIDGET_CreateNETCDFsphere_fixedAltitude
	WIDGET_CreateNETCDFsphere_fixedAltitude = widgets.FloatText(value=120)
	WIDGET_CreateNETCDFsphere_fixedAltitude.description = 'fixedAltitude'
	WIDGET_CreateNETCDFsphere_fixedAltitude.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_CreateNETCDFsphere_fixedAltitude,)
	global WIDGET_CreateNETCDFsphere_LatitudeStep
	WIDGET_CreateNETCDFsphere_LatitudeStep = widgets.FloatText(value=5)
	WIDGET_CreateNETCDFsphere_LatitudeStep.description = 'LatitudeStep'
	WIDGET_CreateNETCDFsphere_LatitudeStep.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_CreateNETCDFsphere_LatitudeStep,)
	global WIDGET_CreateNETCDFsphere_LongitudeStep
	WIDGET_CreateNETCDFsphere_LongitudeStep = widgets.FloatText(value=5)
	WIDGET_CreateNETCDFsphere_LongitudeStep.description = 'LongitudeStep'
	WIDGET_CreateNETCDFsphere_LongitudeStep.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_CreateNETCDFsphere_LongitudeStep,)
	Interpolator_NetCDF_Panel.children += (InputsPanel,)
	# Create widget for the moudle black-box body
	CreateNETCDFsphere_Btn = widgets.Button (
		description='CreateNETCDFsphere',
		tooltip="Creates a CSV file with the format [Time, Lat, Lon, Alt]\nTime and Alt are fixed values given as arguments and Lat and Lon are produced for the whole globe accodring to steps specified.\nCSVfilename: string, the csv filename to be written (overwrite mode)\nfixedDatetimeString: string, the date-time fixed text with format as example: Mar 17 2015 00:12:00.000. Fixed value for every csv entry.\nfixedAltitude: float positive, kilometers above earth. Fixed value for every csv entry.\nLatitudeStep: float positive, Latiude values range between 87.5 and -88.5. The step is how Lat should be incremented at each iteration.\nLongitudeStep: float positive, Longitude values range between -180 and 180.  The step is how Lon should be incremented at each iteration.\nRETURNS: the filename it has written (same as the argument)\n         the altitude (same as the argument)\n         the number of lines written",
	)
	CreateNETCDFsphere_Btn.layout.min_width = '200px'
	CreateNETCDFsphere_Btn.style.button_color = 'gold'
	Interpolator_NetCDF_Panel.children += (CreateNETCDFsphere_Btn,)
	# Create widgets for module's outputs
	OutputsPanel = widgets.VBox()
	OutputsPanel.layout.min_width = '300px'
	WIDGET_CreateNETCDFsphere_theSurfaceFilename = widgets.Label(value='  --> theSurfaceFilename  ')
	WIDGET_CreateNETCDFsphere_theSurfaceFilename.layout.border = '1px dashed green'
	WIDGET_CreateNETCDFsphere_theSurfaceFilename.layout.margin ='0px 40px 0px 0px' 
	WIDGET_CreateNETCDFsphere_theSurfaceFilename.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_CreateNETCDFsphere_theSurfaceFilename,)
	Interpolator_NetCDF_Panel.children += (OutputsPanel,)
	########
	## GUI code for module 'Interpolator'
	########
	# Create widgets for module's inputs
	InputsPanel = widgets.VBox()
	InputsPanel.layout.min_width = '330px'
	global WIDGET_Interpolator_ModelFile_forSurface
	WIDGET_Interpolator_ModelFile_forSurface = widgets.Text(value="/home/NAS/TIEGCM_DATA/TIEGCM_EVT2_2017_Sept_S/tiegcm_s_24900.nc")
	WIDGET_Interpolator_ModelFile_forSurface.description = 'ModelFile_forSurface'
	WIDGET_Interpolator_ModelFile_forSurface.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_Interpolator_ModelFile_forSurface,)
	WIDGET_Interpolator_SurfaceFile = widgets.Label(value='  CreateNETCDFsphere.theSurfaceFilename --> SurfaceFile  ')
	WIDGET_Interpolator_SurfaceFile.layout.border = '1px dashed blue'
	WIDGET_Interpolator_SurfaceFile.layout.padding = '0px 10px 0px 10px'
	InputsPanel.children += (WIDGET_Interpolator_SurfaceFile,)
	Interpolator_NetCDF_Panel.children += (InputsPanel,)
	# Create widget for the moudle black-box body
	Interpolator_Btn = widgets.Button (
		description='Interpolator',
		tooltip="model-->string:: name of model eg TIEGCM\n model_data_file--> string:: model data files stored on DaedalusNAS\n orbit_file-->string :: orbit filename in Time-Lat-Lon-Alt format stored on DaedalusNAS\n save--> Logical:: if true saves interpolated values to directory\n VAR--> string:: variable to inteprolate, must be  one of the variables included in the model data\nOutputs:: Plots+ 1 csv file stored on DaedalusNAS in ModelOutputs/Interpolator",
	)
	Interpolator_Btn.layout.min_width = '200px'
	Interpolator_Btn.style.button_color = 'gold'
	Interpolator_NetCDF_Panel.children += (Interpolator_Btn,)
	# Create widgets for module's outputs
	OutputsPanel = widgets.VBox()
	OutputsPanel.layout.min_width = '300px'
	WIDGET_Interpolator_InterpolatedSurface = widgets.Label(value='  --> InterpolatedSurface  ')
	WIDGET_Interpolator_InterpolatedSurface.layout.border = '1px dashed green'
	WIDGET_Interpolator_InterpolatedSurface.layout.margin ='0px 40px 0px 0px' 
	WIDGET_Interpolator_InterpolatedSurface.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_Interpolator_InterpolatedSurface,)
	Interpolator_NetCDF_Panel.children += (OutputsPanel,)
	########
	## GUI code for module 'PlotGlobe'
	########
	# Create widgets for module's inputs
	InputsPanel = widgets.VBox()
	InputsPanel.layout.min_width = '330px'
	WIDGET_PlotGlobe_SurfaceFilename = widgets.Label(value='  Interpolator.InterpolatedSurface --> SurfaceFilename  ')
	WIDGET_PlotGlobe_SurfaceFilename.layout.border = '1px dashed blue'
	WIDGET_PlotGlobe_SurfaceFilename.layout.padding = '0px 10px 0px 10px'
	InputsPanel.children += (WIDGET_PlotGlobe_SurfaceFilename,)
	global WIDGET_PlotGlobe_SurfaceVariableToPlot
	WIDGET_PlotGlobe_SurfaceVariableToPlot = widgets.Text(value='NE')
	WIDGET_PlotGlobe_SurfaceVariableToPlot.description = 'SurfaceVariableToPlot'
	WIDGET_PlotGlobe_SurfaceVariableToPlot.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_PlotGlobe_SurfaceVariableToPlot,)
	global WIDGET_PlotGlobe_SurfaceColorbarTitle
	WIDGET_PlotGlobe_SurfaceColorbarTitle = widgets.Text()
	WIDGET_PlotGlobe_SurfaceColorbarTitle.description = 'SurfaceColorbarTitle'
	WIDGET_PlotGlobe_SurfaceColorbarTitle.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_PlotGlobe_SurfaceColorbarTitle,)
	global WIDGET_PlotGlobe_SurfaceColorscaleName
	WIDGET_PlotGlobe_SurfaceColorscaleName = widgets.Dropdown( options=['', 'None', 'Blackbody', 'Bluered', 'Blues', 'Earth', 'Electric', 'Greens', 'Greys', 'Hot', 'Jet', 'Picnic', 'Portland', 'Rainbow', 'RdBu', 'Reds', 'Viridis', 'YlGnBu', 'YlOrRd'], value='Jet', description='Variable')
	WIDGET_PlotGlobe_SurfaceColorscaleName.description = 'SurfaceColorscaleName'
	WIDGET_PlotGlobe_SurfaceColorscaleName.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_PlotGlobe_SurfaceColorscaleName,)
	WIDGET_PlotGlobe_OrbitFilename = widgets.Label(value='  Interpolator.InterpolatedOrbit --> OrbitFilename  ')
	WIDGET_PlotGlobe_OrbitFilename.layout.border = '1px dashed blue'
	WIDGET_PlotGlobe_OrbitFilename.layout.padding = '0px 10px 0px 10px'
	InputsPanel.children += (WIDGET_PlotGlobe_OrbitFilename,)
	global WIDGET_PlotGlobe_OrbitVariableToPlot
	WIDGET_PlotGlobe_OrbitVariableToPlot = widgets.Text(value='NE')
	WIDGET_PlotGlobe_OrbitVariableToPlot.description = 'OrbitVariableToPlot'
	WIDGET_PlotGlobe_OrbitVariableToPlot.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_PlotGlobe_OrbitVariableToPlot,)
	global WIDGET_PlotGlobe_OrbitColorbarTitle
	WIDGET_PlotGlobe_OrbitColorbarTitle = widgets.Text()
	WIDGET_PlotGlobe_OrbitColorbarTitle.description = 'OrbitColorbarTitle'
	WIDGET_PlotGlobe_OrbitColorbarTitle.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_PlotGlobe_OrbitColorbarTitle,)
	global WIDGET_PlotGlobe_OrbitColorscaleName
	WIDGET_PlotGlobe_OrbitColorscaleName = widgets.Dropdown( options=['', 'None', 'Blackbody', 'Bluered', 'Blues', 'Earth', 'Electric', 'Greens', 'Greys', 'Hot', 'Jet', 'Picnic', 'Portland', 'Rainbow', 'RdBu', 'Reds', 'Viridis', 'YlGnBu', 'YlOrRd'], value='Jet', description='Variable')
	WIDGET_PlotGlobe_OrbitColorscaleName.description = 'OrbitColorscaleName'
	WIDGET_PlotGlobe_OrbitColorscaleName.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_PlotGlobe_OrbitColorscaleName,)
	global WIDGET_PlotGlobe_PlotTitle
	WIDGET_PlotGlobe_PlotTitle = widgets.Text()
	WIDGET_PlotGlobe_PlotTitle.description = 'PlotTitle'
	WIDGET_PlotGlobe_PlotTitle.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_PlotGlobe_PlotTitle,)
	Interpolator_NetCDF_Panel.children += (InputsPanel,)
	# Create widget for the moudle black-box body
	PlotGlobe_Btn = widgets.Button (
		description='PlotGlobe',
		tooltip="Creates a 3D plot of an earth globe. Can plot a sphere surface and/or a satellite orbit. \nThe surface and the orbit are colored according to the selected variable values in the data files. \nThe data files can be of CSV o NetCDF4 format.\nThe same or different variables can be selected for plotting at the surface and the orbit.\nThe same or different colorscales can be selected for the surface and the orbit.\nThe same or different value range can be selected for the surface and the orbit by assigning the same or different colorbar titles.\nValid Colorscale Names: ‘Blackbody’, ‘Bluered’, ‘Blues’, ‘Earth’, ‘Electric’, ‘Greens’, ‘Greys’, ‘Hot’, ‘Jet’, ‘Picnic’, ‘Portland’, ‘Rainbow’, ‘RdBu’, ‘Reds’, ‘Viridis’, ‘YlGnBu’, ‘YlOrRd’\nARGUMENTS:\n  SurfaceFilename: The file which contains data for a surface. If empty string then no surface will be plotted. \n                   CSV format: Time,Lat,Lon,Alt,value. Contains the data for the sphere surface. \n                   NetCDF format: see Daedalus User Manual\n  SurfaceVariableToPlot: The name of the variable to read from the surface data file and plot upon the globe.\n  SurfaceColorbarTitle: A title to display above the colorbar which refers to the surface data\n  SurfaceColorscaleName: The name of the Colorscale to use for the surface data. \n                         In case an empty string is passed then a default HeatMap colorscale will be applied.\n                         In case None is passed then all points will be black irrespective of value.\n  OrbitFilename: counterpart of the corresponding argument for the Surface data\n  OrbitVariableToPlot: counterpart of the corresponding argument for the Surface data \n  OrbitColorbarTitle: counterpart of the corresponding argument for the Surface data \n  OrbitColorscaleName: counterpart of the corresponding argument for the Surface data\n  PlotTitle: A title which is displayed on the top of the plot.\nRETURNS: a string containing information about the Data\n",
	)
	PlotGlobe_Btn.layout.min_width = '200px'
	PlotGlobe_Btn.style.button_color = 'deeppink'
	Interpolator_NetCDF_Panel.children += (PlotGlobe_Btn,)
	# Create widgets for module's outputs
	OutputsPanel = widgets.VBox()
	OutputsPanel.layout.min_width = '300px'
	Interpolator_NetCDF_Panel.children += (OutputsPanel,)
	return Interpolator_NetCDF_Panel


def Construct_ExecuteTest1():
	#Create Containers
	ExecuteTest1_Panel = widgets.Box()
	ExecuteTest1_Panel.layout.overflow_x = 'scroll'
	ExecuteTest1_Panel.layout.flex_flow = 'row'
	ExecuteTest1_Panel.layout.display = 'flex'
	########
	## GUI code for module 'Test1'
	########
	# Create widgets for module's inputs
	InputsPanel = widgets.VBox()
	InputsPanel.layout.min_width = '330px'
	ExecuteTest1_Panel.children += (InputsPanel,)
	# Create widget for the moudle black-box body
	Test1_Btn = widgets.Button (
		description='Test1',
		tooltip="",
	)
	Test1_Btn.layout.min_width = '200px'
	Test1_Btn.style.button_color = 'gold'
	ExecuteTest1_Panel.children += (Test1_Btn,)
	# Create widgets for module's outputs
	OutputsPanel = widgets.VBox()
	OutputsPanel.layout.min_width = '300px'
	WIDGET_Test1_txt = widgets.Label(value='  --> txt  ')
	WIDGET_Test1_txt.layout.border = '1px dashed green'
	WIDGET_Test1_txt.layout.margin ='0px 40px 0px 0px' 
	WIDGET_Test1_txt.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_Test1_txt,)
	ExecuteTest1_Panel.children += (OutputsPanel,)
	return ExecuteTest1_Panel


def Construct_HWM14():
	#Create Containers
	HWM14_Panel = widgets.Box()
	HWM14_Panel.layout.overflow_x = 'scroll'
	HWM14_Panel.layout.flex_flow = 'row'
	HWM14_Panel.layout.display = 'flex'
	########
	## GUI code for module 'OrbitSelector'
	########
	# Create widgets for module's inputs
	InputsPanel = widgets.VBox()
	InputsPanel.layout.min_width = '330px'
	global WIDGET_OrbitSelector_Filename
	WIDGET_OrbitSelector_Filename = widgets.Text()
	WIDGET_OrbitSelector_Filename.description = 'Filename'
	WIDGET_OrbitSelector_Filename.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_OrbitSelector_Filename,)
	global WIDGET_OrbitSelector_EvtXY
	WIDGET_OrbitSelector_EvtXY = widgets.Dropdown( options=['', 'Evt0', 'Evt1', 'Evt2', 'Evt3', 'Evt4', 'Evt5', 'Evt6', 'Evt7', 'Evt8', 'Evt9'], description='Variable')
	WIDGET_OrbitSelector_EvtXY.description = 'EvtXY'
	WIDGET_OrbitSelector_EvtXY.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_OrbitSelector_EvtXY,)
	global WIDGET_OrbitSelector_TYP
	WIDGET_OrbitSelector_TYP = widgets.Dropdown( options=['LLA', 'VEL', 'PTG'], description='Variable')
	WIDGET_OrbitSelector_TYP.description = 'TYP'
	WIDGET_OrbitSelector_TYP.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_OrbitSelector_TYP,)
	global WIDGET_OrbitSelector_PerYYY
	WIDGET_OrbitSelector_PerYYY = widgets.Dropdown( options=['Per120', 'Per150'], description='Variable')
	WIDGET_OrbitSelector_PerYYY.description = 'PerYYY'
	WIDGET_OrbitSelector_PerYYY.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_OrbitSelector_PerYYY,)
	global WIDGET_OrbitSelector_LatZZ
	WIDGET_OrbitSelector_LatZZ = widgets.Dropdown( options=['Lat00', 'Lat40', 'Lat80'],  description='Variable')
	WIDGET_OrbitSelector_LatZZ.description = 'LatZZ'
	WIDGET_OrbitSelector_LatZZ.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_OrbitSelector_LatZZ,)
	global WIDGET_OrbitSelector_SRXXHZ
	WIDGET_OrbitSelector_SRXXHZ = widgets.Dropdown( options=['Srt16Hz', 'Srt01Hz'],  description='Variable')
	WIDGET_OrbitSelector_SRXXHZ.description = 'SRXXHZ'
	WIDGET_OrbitSelector_SRXXHZ.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_OrbitSelector_SRXXHZ,)
	global WIDGET_OrbitSelector_SC
	WIDGET_OrbitSelector_SC = widgets.Dropdown( options=['Msc', 'Ssc'],  description='Variable')
	WIDGET_OrbitSelector_SC.description = 'SC'
	WIDGET_OrbitSelector_SC.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_OrbitSelector_SC,)
	HWM14_Panel.children += (InputsPanel,)
	# Create widget for the moudle black-box body
	OrbitSelector_Btn = widgets.Button (
		description='OrbitSelector',
		tooltip="Constructs and returns an orbit full-path filename from the orbit properties.\nThe filename includes all parameters required to select the appropriate orbit.\nAll filenames  have the same length and same number of parameters. \nIn case OrbitFilename and EvtXY are empty then an empty string is returned.\nFILE FORMAT: DAED_ORB_EVTXY_TYP_PERYYY_LATZZ_SRTQQHz_XSC.csv\n  Parameter OrbitFilename:\n    If it is not empty then the rest arguments are ignored. \n    If it contains slashes then is is assumed it is a full path name and it is returned as it is.\n    If it does not contain slashes then the orbits path is added.\n    \n  Parameter EvtXY values:\n    EVTXS    X Event, Single Orbit\n    EVTXA    X Event, All Orbit\n    EVT1Y    1st Event: St Patrick’s day event [17 Mar 2015 – 20 Mar 2015]\n    EVT2Y    2nd Event ...  ...\n \n  Parameter TYP values:\n    LLA    Time,Latitude Longitude Altitude \n    VEL    Time,VmagVxVyVz\n    PTG    Time, X_GSE, Y_GSE, Z_GSE, RamX_GSE, RamY_GSE, RamZ_GSE\n  Parameter PerYYY values:\n    PER120, PER150    Perigee Altitude at 120km or 150km\n  Parameter LatZZ values:\n    LAT00, LAT40, LAT80    Perigee Latitude at 0°, 40° or 80°\n \n  Parameter SRXXHZ values:\n    SRT16Hz, SRT01Hz    Sampling rate 16Hz or 1HZ\n \n  Parameter SC values:\n    MSC    Mother Spacecraft\n    SSC    Sub-Spacecraft\n",
	)
	OrbitSelector_Btn.layout.min_width = '200px'
	OrbitSelector_Btn.style.button_color = 'YellowGreen'
	HWM14_Panel.children += (OrbitSelector_Btn,)
	# Create widgets for module's outputs
	OutputsPanel = widgets.VBox()
	OutputsPanel.layout.min_width = '300px'
	WIDGET_OrbitSelector_OrbitCSVfilename = widgets.Label(value='  --> OrbitCSVfilename  ')
	WIDGET_OrbitSelector_OrbitCSVfilename.layout.border = '1px dashed green'
	WIDGET_OrbitSelector_OrbitCSVfilename.layout.margin ='0px 40px 0px 0px' 
	WIDGET_OrbitSelector_OrbitCSVfilename.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_OrbitSelector_OrbitCSVfilename,)
	WIDGET_OrbitSelector_OrbitNetCDFfilename = widgets.Label(value='  --> OrbitNetCDFfilename  ')
	WIDGET_OrbitSelector_OrbitNetCDFfilename.layout.border = '1px dashed green'
	WIDGET_OrbitSelector_OrbitNetCDFfilename.layout.margin ='0px 40px 0px 0px' 
	WIDGET_OrbitSelector_OrbitNetCDFfilename.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_OrbitSelector_OrbitNetCDFfilename,)
	HWM14_Panel.children += (OutputsPanel,)
	########
	## GUI code for module 'HWM14'
	########
	# Create widgets for module's inputs
	InputsPanel = widgets.VBox()
	InputsPanel.layout.min_width = '330px'
	WIDGET_HWM14_CSVfilename_forOrbit = widgets.Label(value='  OrbitSelector.OrbitCSVfilename --> CSVfilename_forOrbit  ')
	WIDGET_HWM14_CSVfilename_forOrbit.layout.border = '1px dashed blue'
	WIDGET_HWM14_CSVfilename_forOrbit.layout.padding = '0px 10px 0px 10px'
	InputsPanel.children += (WIDGET_HWM14_CSVfilename_forOrbit,)
	global WIDGET_HWM14_Variable_forOrbit
	WIDGET_HWM14_Variable_forOrbit = widgets.Dropdown( options=['', 'v', 'u'], value='', description='Variable')
	WIDGET_HWM14_Variable_forOrbit.description = 'Variable_forOrbit'
	WIDGET_HWM14_Variable_forOrbit.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_HWM14_Variable_forOrbit,)
	HWM14_Panel.children += (InputsPanel,)
	# Create widget for the moudle black-box body
	HWM14_Btn = widgets.Button (
		description='HWM14',
		tooltip="Horizontal Wind Model empirical (HWM) \nfilename: the CSV file to read. Format: Time,Alt,Lon,Lat. The file defines the coordinates for which the HWM calculation will take place.\nParameter: selects which values to calculate. Possible values are:\n    ''  for all\n    'v' for meridional wind (northward) \n    'u' for zonal wind (eastward)  \nReturns: the name of the CSV file which contains the calculated values. Format: Time,Alt,Lon,Lat,v,u \n",
	)
	HWM14_Btn.layout.min_width = '200px'
	HWM14_Btn.style.button_color = 'gold'
	HWM14_Panel.children += (HWM14_Btn,)
	# Create widgets for module's outputs
	OutputsPanel = widgets.VBox()
	OutputsPanel.layout.min_width = '300px'
	WIDGET_HWM14_OrbitResultCSV = widgets.Label(value='  --> OrbitResultCSV  ')
	WIDGET_HWM14_OrbitResultCSV.layout.border = '1px dashed green'
	WIDGET_HWM14_OrbitResultCSV.layout.margin ='0px 40px 0px 0px' 
	WIDGET_HWM14_OrbitResultCSV.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_HWM14_OrbitResultCSV,)
	HWM14_Panel.children += (OutputsPanel,)
	########
	## GUI code for module 'CreateCSV_Sphere'
	########
	# Create widgets for module's inputs
	InputsPanel = widgets.VBox()
	InputsPanel.layout.min_width = '330px'
	global WIDGET_CreateCSV_Sphere_CSVfilename
	WIDGET_CreateCSV_Sphere_CSVfilename = widgets.Text(value="")
	WIDGET_CreateCSV_Sphere_CSVfilename.description = 'CSVfilename'
	WIDGET_CreateCSV_Sphere_CSVfilename.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_CreateCSV_Sphere_CSVfilename,)
	global WIDGET_CreateCSV_Sphere_fixedDatetimeString
	WIDGET_CreateCSV_Sphere_fixedDatetimeString = widgets.Text(value="Mar 17 2015 00:50:00.000000000")
	WIDGET_CreateCSV_Sphere_fixedDatetimeString.description = 'fixedDatetimeString'
	WIDGET_CreateCSV_Sphere_fixedDatetimeString.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_CreateCSV_Sphere_fixedDatetimeString,)
	global WIDGET_CreateCSV_Sphere_fixedAltitude
	WIDGET_CreateCSV_Sphere_fixedAltitude = widgets.FloatText(value=120)
	WIDGET_CreateCSV_Sphere_fixedAltitude.description = 'fixedAltitude'
	WIDGET_CreateCSV_Sphere_fixedAltitude.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_CreateCSV_Sphere_fixedAltitude,)
	global WIDGET_CreateCSV_Sphere_LatitudeStep
	WIDGET_CreateCSV_Sphere_LatitudeStep = widgets.FloatText(value=5)
	WIDGET_CreateCSV_Sphere_LatitudeStep.description = 'LatitudeStep'
	WIDGET_CreateCSV_Sphere_LatitudeStep.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_CreateCSV_Sphere_LatitudeStep,)
	global WIDGET_CreateCSV_Sphere_LongitudeStep
	WIDGET_CreateCSV_Sphere_LongitudeStep = widgets.FloatText(value=5)
	WIDGET_CreateCSV_Sphere_LongitudeStep.description = 'LongitudeStep'
	WIDGET_CreateCSV_Sphere_LongitudeStep.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_CreateCSV_Sphere_LongitudeStep,)
	HWM14_Panel.children += (InputsPanel,)
	# Create widget for the moudle black-box body
	CreateCSV_Sphere_Btn = widgets.Button (
		description='CreateCSV_Sphere',
		tooltip="Creates a CSV file with the format [Time, Lat, Lon, Alt]\nTime and Alt are fixed values given as arguments and Lat and Lon are produced for the whole globe accodring to steps specified.\nCSVfilename: string, the csv filename to be written (overwrite mode)\nfixedDatetimeString: string, the date-time fixed text with format as example: Mar 17 2015 00:12:00.000. Fixed value for every csv entry.\nfixedAltitude: float positive, kilometers above earth. Fixed value for every csv entry.\nLatitudeStep: float positive, Latiude values range between 87.5 and -88.5. The step is how Lat should be incremented at each iteration.\nLongitudeStep: float positive, Longitude values range between -180 and 180.  The step is how Lon should be incremented at each iteration.\nRETURNS: the filename it has written (same as the argument)\n         the altitude (same as the argument)\n         the number of lines written",
	)
	CreateCSV_Sphere_Btn.layout.min_width = '200px'
	CreateCSV_Sphere_Btn.style.button_color = 'gold'
	HWM14_Panel.children += (CreateCSV_Sphere_Btn,)
	# Create widgets for module's outputs
	OutputsPanel = widgets.VBox()
	OutputsPanel.layout.min_width = '300px'
	WIDGET_CreateCSV_Sphere_theCSVfilename = widgets.Label(value='  --> theCSVfilename  ')
	WIDGET_CreateCSV_Sphere_theCSVfilename.layout.border = '1px dashed green'
	WIDGET_CreateCSV_Sphere_theCSVfilename.layout.margin ='0px 40px 0px 0px' 
	WIDGET_CreateCSV_Sphere_theCSVfilename.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_CreateCSV_Sphere_theCSVfilename,)
	WIDGET_CreateCSV_Sphere_theAltitude = widgets.Label(value='  --> theAltitude  ')
	WIDGET_CreateCSV_Sphere_theAltitude.layout.border = '1px dashed green'
	WIDGET_CreateCSV_Sphere_theAltitude.layout.margin ='0px 40px 0px 0px' 
	WIDGET_CreateCSV_Sphere_theAltitude.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_CreateCSV_Sphere_theAltitude,)
	WIDGET_CreateCSV_Sphere_NumberOfLinesWritten = widgets.Label(value='  --> NumberOfLinesWritten  ')
	WIDGET_CreateCSV_Sphere_NumberOfLinesWritten.layout.border = '1px dashed green'
	WIDGET_CreateCSV_Sphere_NumberOfLinesWritten.layout.margin ='0px 40px 0px 0px' 
	WIDGET_CreateCSV_Sphere_NumberOfLinesWritten.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_CreateCSV_Sphere_NumberOfLinesWritten,)
	HWM14_Panel.children += (OutputsPanel,)
	########
	## GUI code for module 'HWM14'
	########
	# Create widgets for module's inputs
	InputsPanel = widgets.VBox()
	InputsPanel.layout.min_width = '330px'
	WIDGET_HWM14_CSVfilename_forMeshgrid = widgets.Label(value='  CreateCSV_Sphere.theCSVfilename --> CSVfilename_forMeshgrid  ')
	WIDGET_HWM14_CSVfilename_forMeshgrid.layout.border = '1px dashed blue'
	WIDGET_HWM14_CSVfilename_forMeshgrid.layout.padding = '0px 10px 0px 10px'
	InputsPanel.children += (WIDGET_HWM14_CSVfilename_forMeshgrid,)
	global WIDGET_HWM14_Variable_forMeshgrid
	WIDGET_HWM14_Variable_forMeshgrid = widgets.Dropdown( options=['', 'v', 'u'], value='', description='Variable')
	WIDGET_HWM14_Variable_forMeshgrid.description = 'Variable_forMeshgrid'
	WIDGET_HWM14_Variable_forMeshgrid.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_HWM14_Variable_forMeshgrid,)
	HWM14_Panel.children += (InputsPanel,)
	# Create widget for the moudle black-box body
	HWM14_Btn = widgets.Button (
		description='HWM14',
		tooltip="Horizontal Wind Model empirical (HWM) \nfilename: the CSV file to read. Format: Time,Alt,Lon,Lat. The file defines the coordinates for which the HWM calculation will take place.\nParameter: selects which values to calculate. Possible values are:\n    ''  for all\n    'v' for meridional wind (northward) \n    'u' for zonal wind (eastward)  \nReturns: the name of the CSV file which contains the calculated values. Format: Time,Alt,Lon,Lat,v,u \n",
	)
	HWM14_Btn.layout.min_width = '200px'
	HWM14_Btn.style.button_color = 'gold'
	HWM14_Panel.children += (HWM14_Btn,)
	# Create widgets for module's outputs
	OutputsPanel = widgets.VBox()
	OutputsPanel.layout.min_width = '300px'
	WIDGET_HWM14_MeshgridResultCSV = widgets.Label(value='  --> MeshgridResultCSV  ')
	WIDGET_HWM14_MeshgridResultCSV.layout.border = '1px dashed green'
	WIDGET_HWM14_MeshgridResultCSV.layout.margin ='0px 40px 0px 0px' 
	WIDGET_HWM14_MeshgridResultCSV.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_HWM14_MeshgridResultCSV,)
	HWM14_Panel.children += (OutputsPanel,)
	########
	## GUI code for module 'PlotGlobe'
	########
	# Create widgets for module's inputs
	InputsPanel = widgets.VBox()
	InputsPanel.layout.min_width = '330px'
	WIDGET_PlotGlobe_SurfaceFilename = widgets.Label(value='  HWM14.MeshgridResultCSV --> SurfaceFilename  ')
	WIDGET_PlotGlobe_SurfaceFilename.layout.border = '1px dashed blue'
	WIDGET_PlotGlobe_SurfaceFilename.layout.padding = '0px 10px 0px 10px'
	InputsPanel.children += (WIDGET_PlotGlobe_SurfaceFilename,)
	global WIDGET_PlotGlobe_SurfaceVariableToPlot
	WIDGET_PlotGlobe_SurfaceVariableToPlot = widgets.Text(value='ALFA')
	WIDGET_PlotGlobe_SurfaceVariableToPlot.description = 'SurfaceVariableToPlot'
	WIDGET_PlotGlobe_SurfaceVariableToPlot.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_PlotGlobe_SurfaceVariableToPlot,)
	global WIDGET_PlotGlobe_SurfaceColorbarTitle
	WIDGET_PlotGlobe_SurfaceColorbarTitle = widgets.Text()
	WIDGET_PlotGlobe_SurfaceColorbarTitle.description = 'SurfaceColorbarTitle'
	WIDGET_PlotGlobe_SurfaceColorbarTitle.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_PlotGlobe_SurfaceColorbarTitle,)
	global WIDGET_PlotGlobe_SurfaceColorscaleName
	WIDGET_PlotGlobe_SurfaceColorscaleName = widgets.Dropdown( options=['', 'None', 'Blackbody', 'Bluered', 'Blues', 'Earth', 'Electric', 'Greens', 'Greys', 'Hot', 'Jet', 'Picnic', 'Portland', 'Rainbow', 'RdBu', 'Reds', 'Viridis', 'YlGnBu', 'YlOrRd'], value='Jet', description='Variable')
	WIDGET_PlotGlobe_SurfaceColorscaleName.description = 'SurfaceColorscaleName'
	WIDGET_PlotGlobe_SurfaceColorscaleName.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_PlotGlobe_SurfaceColorscaleName,)
	WIDGET_PlotGlobe_OrbitFilename = widgets.Label(value='  HWM14.OrbitResultCSV --> OrbitFilename  ')
	WIDGET_PlotGlobe_OrbitFilename.layout.border = '1px dashed blue'
	WIDGET_PlotGlobe_OrbitFilename.layout.padding = '0px 10px 0px 10px'
	InputsPanel.children += (WIDGET_PlotGlobe_OrbitFilename,)
	global WIDGET_PlotGlobe_OrbitVariableToPlot
	WIDGET_PlotGlobe_OrbitVariableToPlot = widgets.Text(value='ALFA')
	WIDGET_PlotGlobe_OrbitVariableToPlot.description = 'OrbitVariableToPlot'
	WIDGET_PlotGlobe_OrbitVariableToPlot.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_PlotGlobe_OrbitVariableToPlot,)
	global WIDGET_PlotGlobe_OrbitColorbarTitle
	WIDGET_PlotGlobe_OrbitColorbarTitle = widgets.Text()
	WIDGET_PlotGlobe_OrbitColorbarTitle.description = 'OrbitColorbarTitle'
	WIDGET_PlotGlobe_OrbitColorbarTitle.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_PlotGlobe_OrbitColorbarTitle,)
	global WIDGET_PlotGlobe_OrbitColorscaleName
	WIDGET_PlotGlobe_OrbitColorscaleName = widgets.Dropdown( options=['', 'None', 'Blackbody', 'Bluered', 'Blues', 'Earth', 'Electric', 'Greens', 'Greys', 'Hot', 'Jet', 'Picnic', 'Portland', 'Rainbow', 'RdBu', 'Reds', 'Viridis', 'YlGnBu', 'YlOrRd'], value='Jet', description='Variable')
	WIDGET_PlotGlobe_OrbitColorscaleName.description = 'OrbitColorscaleName'
	WIDGET_PlotGlobe_OrbitColorscaleName.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_PlotGlobe_OrbitColorscaleName,)
	global WIDGET_PlotGlobe_PlotTitle
	WIDGET_PlotGlobe_PlotTitle = widgets.Text()
	WIDGET_PlotGlobe_PlotTitle.description = 'PlotTitle'
	WIDGET_PlotGlobe_PlotTitle.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_PlotGlobe_PlotTitle,)
	HWM14_Panel.children += (InputsPanel,)
	# Create widget for the moudle black-box body
	PlotGlobe_Btn = widgets.Button (
		description='PlotGlobe',
		tooltip="Creates a 3D plot of an earth globe. Can plot a sphere surface and/or a satellite orbit. \nThe surface and the orbit are colored according to the selected variable values in the data files. \nThe data files can be of CSV o NetCDF4 format.\nThe same or different variables can be selected for plotting at the surface and the orbit.\nThe same or different colorscales can be selected for the surface and the orbit.\nThe same or different value range can be selected for the surface and the orbit by assigning the same or different colorbar titles.\nValid Colorscale Names: ‘Blackbody’, ‘Bluered’, ‘Blues’, ‘Earth’, ‘Electric’, ‘Greens’, ‘Greys’, ‘Hot’, ‘Jet’, ‘Picnic’, ‘Portland’, ‘Rainbow’, ‘RdBu’, ‘Reds’, ‘Viridis’, ‘YlGnBu’, ‘YlOrRd’\nARGUMENTS:\n  SurfaceFilename: The file which contains data for a surface. If empty string then no surface will be plotted. \n                   CSV format: Time,Lat,Lon,Alt,value. Contains the data for the sphere surface. \n                   NetCDF format: see Daedalus User Manual\n  SurfaceVariableToPlot: The name of the variable to read from the surface data file and plot upon the globe.\n  SurfaceColorbarTitle: A title to display above the colorbar which refers to the surface data\n  SurfaceColorscaleName: The name of the Colorscale to use for the surface data. \n                         In case an empty string is passed then a default HeatMap colorscale will be applied.\n                         In case None is passed then all points will be black irrespective of value.\n  OrbitFilename: counterpart of the corresponding argument for the Surface data\n  OrbitVariableToPlot: counterpart of the corresponding argument for the Surface data \n  OrbitColorbarTitle: counterpart of the corresponding argument for the Surface data \n  OrbitColorscaleName: counterpart of the corresponding argument for the Surface data\n  PlotTitle: A title which is displayed on the top of the plot.\nRETURNS: a string containing information about the Data\n",
	)
	PlotGlobe_Btn.layout.min_width = '200px'
	PlotGlobe_Btn.style.button_color = 'deeppink'
	HWM14_Panel.children += (PlotGlobe_Btn,)
	# Create widgets for module's outputs
	OutputsPanel = widgets.VBox()
	OutputsPanel.layout.min_width = '300px'
	HWM14_Panel.children += (OutputsPanel,)
	return HWM14_Panel


def Construct_IRI16():
	#Create Containers
	IRI16_Panel = widgets.Box()
	IRI16_Panel.layout.overflow_x = 'scroll'
	IRI16_Panel.layout.flex_flow = 'row'
	IRI16_Panel.layout.display = 'flex'
	########
	## GUI code for module 'OrbitSelector'
	########
	# Create widgets for module's inputs
	InputsPanel = widgets.VBox()
	InputsPanel.layout.min_width = '330px'
	global WIDGET_OrbitSelector_Filename
	WIDGET_OrbitSelector_Filename = widgets.Text()
	WIDGET_OrbitSelector_Filename.description = 'Filename'
	WIDGET_OrbitSelector_Filename.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_OrbitSelector_Filename,)
	global WIDGET_OrbitSelector_EvtXY
	WIDGET_OrbitSelector_EvtXY = widgets.Dropdown( options=['', 'Evt0', 'Evt1', 'Evt2', 'Evt3', 'Evt4', 'Evt5', 'Evt6', 'Evt7', 'Evt8', 'Evt9'], description='Variable')
	WIDGET_OrbitSelector_EvtXY.description = 'EvtXY'
	WIDGET_OrbitSelector_EvtXY.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_OrbitSelector_EvtXY,)
	global WIDGET_OrbitSelector_TYP
	WIDGET_OrbitSelector_TYP = widgets.Dropdown( options=['LLA', 'VEL', 'PTG'], description='Variable')
	WIDGET_OrbitSelector_TYP.description = 'TYP'
	WIDGET_OrbitSelector_TYP.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_OrbitSelector_TYP,)
	global WIDGET_OrbitSelector_PerYYY
	WIDGET_OrbitSelector_PerYYY = widgets.Dropdown( options=['Per120', 'Per150'], description='Variable')
	WIDGET_OrbitSelector_PerYYY.description = 'PerYYY'
	WIDGET_OrbitSelector_PerYYY.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_OrbitSelector_PerYYY,)
	global WIDGET_OrbitSelector_LatZZ
	WIDGET_OrbitSelector_LatZZ = widgets.Dropdown( options=['Lat00', 'Lat40', 'Lat80'],  description='Variable')
	WIDGET_OrbitSelector_LatZZ.description = 'LatZZ'
	WIDGET_OrbitSelector_LatZZ.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_OrbitSelector_LatZZ,)
	global WIDGET_OrbitSelector_SRXXHZ
	WIDGET_OrbitSelector_SRXXHZ = widgets.Dropdown( options=['Srt16Hz', 'Srt01Hz'],  description='Variable')
	WIDGET_OrbitSelector_SRXXHZ.description = 'SRXXHZ'
	WIDGET_OrbitSelector_SRXXHZ.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_OrbitSelector_SRXXHZ,)
	global WIDGET_OrbitSelector_SC
	WIDGET_OrbitSelector_SC = widgets.Dropdown( options=['Msc', 'Ssc'],  description='Variable')
	WIDGET_OrbitSelector_SC.description = 'SC'
	WIDGET_OrbitSelector_SC.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_OrbitSelector_SC,)
	IRI16_Panel.children += (InputsPanel,)
	# Create widget for the moudle black-box body
	OrbitSelector_Btn = widgets.Button (
		description='OrbitSelector',
		tooltip="Constructs and returns an orbit full-path filename from the orbit properties.\nThe filename includes all parameters required to select the appropriate orbit.\nAll filenames  have the same length and same number of parameters. \nIn case OrbitFilename and EvtXY are empty then an empty string is returned.\nFILE FORMAT: DAED_ORB_EVTXY_TYP_PERYYY_LATZZ_SRTQQHz_XSC.csv\n  Parameter OrbitFilename:\n    If it is not empty then the rest arguments are ignored. \n    If it contains slashes then is is assumed it is a full path name and it is returned as it is.\n    If it does not contain slashes then the orbits path is added.\n    \n  Parameter EvtXY values:\n    EVTXS    X Event, Single Orbit\n    EVTXA    X Event, All Orbit\n    EVT1Y    1st Event: St Patrick’s day event [17 Mar 2015 – 20 Mar 2015]\n    EVT2Y    2nd Event ...  ...\n \n  Parameter TYP values:\n    LLA    Time,Latitude Longitude Altitude \n    VEL    Time,VmagVxVyVz\n    PTG    Time, X_GSE, Y_GSE, Z_GSE, RamX_GSE, RamY_GSE, RamZ_GSE\n  Parameter PerYYY values:\n    PER120, PER150    Perigee Altitude at 120km or 150km\n  Parameter LatZZ values:\n    LAT00, LAT40, LAT80    Perigee Latitude at 0°, 40° or 80°\n \n  Parameter SRXXHZ values:\n    SRT16Hz, SRT01Hz    Sampling rate 16Hz or 1HZ\n \n  Parameter SC values:\n    MSC    Mother Spacecraft\n    SSC    Sub-Spacecraft\n",
	)
	OrbitSelector_Btn.layout.min_width = '200px'
	OrbitSelector_Btn.style.button_color = 'YellowGreen'
	IRI16_Panel.children += (OrbitSelector_Btn,)
	# Create widgets for module's outputs
	OutputsPanel = widgets.VBox()
	OutputsPanel.layout.min_width = '300px'
	WIDGET_OrbitSelector_OrbitCSVfilename = widgets.Label(value='  --> OrbitCSVfilename  ')
	WIDGET_OrbitSelector_OrbitCSVfilename.layout.border = '1px dashed green'
	WIDGET_OrbitSelector_OrbitCSVfilename.layout.margin ='0px 40px 0px 0px' 
	WIDGET_OrbitSelector_OrbitCSVfilename.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_OrbitSelector_OrbitCSVfilename,)
	WIDGET_OrbitSelector_OrbitNetCDFfilename = widgets.Label(value='  --> OrbitNetCDFfilename  ')
	WIDGET_OrbitSelector_OrbitNetCDFfilename.layout.border = '1px dashed green'
	WIDGET_OrbitSelector_OrbitNetCDFfilename.layout.margin ='0px 40px 0px 0px' 
	WIDGET_OrbitSelector_OrbitNetCDFfilename.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_OrbitSelector_OrbitNetCDFfilename,)
	IRI16_Panel.children += (OutputsPanel,)
	########
	## GUI code for module 'Selector'
	########
	# Create widgets for module's inputs
	InputsPanel = widgets.VBox()
	InputsPanel.layout.min_width = '330px'
	global WIDGET_Selector_MeshgridCSVfilename
	WIDGET_Selector_MeshgridCSVfilename = widgets.Text(value="")
	WIDGET_Selector_MeshgridCSVfilename.description = 'MeshgridCSVfilename'
	WIDGET_Selector_MeshgridCSVfilename.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_Selector_MeshgridCSVfilename,)
	global WIDGET_Selector_MeshgridDatetime
	WIDGET_Selector_MeshgridDatetime = widgets.Text(value="Mar 17 2015 00:50:00.000000000")
	WIDGET_Selector_MeshgridDatetime.description = 'MeshgridDatetime'
	WIDGET_Selector_MeshgridDatetime.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_Selector_MeshgridDatetime,)
	global WIDGET_Selector_MeshgridAltitude
	WIDGET_Selector_MeshgridAltitude = widgets.FloatText(value=120)
	WIDGET_Selector_MeshgridAltitude.description = 'MeshgridAltitude'
	WIDGET_Selector_MeshgridAltitude.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_Selector_MeshgridAltitude,)
	global WIDGET_Selector_MeshgridLatitudeStep
	WIDGET_Selector_MeshgridLatitudeStep = widgets.FloatText(value=5)
	WIDGET_Selector_MeshgridLatitudeStep.description = 'MeshgridLatitudeStep'
	WIDGET_Selector_MeshgridLatitudeStep.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_Selector_MeshgridLatitudeStep,)
	global WIDGET_Selector_MeshgridLongitudeStep
	WIDGET_Selector_MeshgridLongitudeStep = widgets.FloatText(value=5)
	WIDGET_Selector_MeshgridLongitudeStep.description = 'MeshgridLongitudeStep'
	WIDGET_Selector_MeshgridLongitudeStep.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_Selector_MeshgridLongitudeStep,)
	global WIDGET_Selector_NoiseWx
	WIDGET_Selector_NoiseWx = widgets.FloatText(value=0)
	WIDGET_Selector_NoiseWx.description = 'NoiseWx'
	WIDGET_Selector_NoiseWx.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_Selector_NoiseWx,)
	global WIDGET_Selector_NoiseWy
	WIDGET_Selector_NoiseWy = widgets.FloatText(value=0)
	WIDGET_Selector_NoiseWy.description = 'NoiseWy'
	WIDGET_Selector_NoiseWy.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_Selector_NoiseWy,)
	global WIDGET_Selector_NoiseAx
	WIDGET_Selector_NoiseAx = widgets.FloatText(value=0)
	WIDGET_Selector_NoiseAx.description = 'NoiseAx'
	WIDGET_Selector_NoiseAx.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_Selector_NoiseAx,)
	global WIDGET_Selector_NoiseAy
	WIDGET_Selector_NoiseAy = widgets.FloatText(value=0)
	WIDGET_Selector_NoiseAy.description = 'NoiseAy'
	WIDGET_Selector_NoiseAy.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_Selector_NoiseAy,)
	global WIDGET_Selector_IRI16_variableForMeshgrid
	WIDGET_Selector_IRI16_variableForMeshgrid = widgets.Dropdown( options=['', 'Op', 'O2p', 'Te', 'Ti', 'Ne'], value='', description='Variable')
	WIDGET_Selector_IRI16_variableForMeshgrid.description = 'IRI16_variableForMeshgrid'
	WIDGET_Selector_IRI16_variableForMeshgrid.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_Selector_IRI16_variableForMeshgrid,)
	global WIDGET_Selector_IRI16_variableForOrbit
	WIDGET_Selector_IRI16_variableForOrbit = widgets.Dropdown( options=['', 'Op', 'O2p', 'Te', 'Ti', 'Ne'], value='', description='Variable')
	WIDGET_Selector_IRI16_variableForOrbit.description = 'IRI16_variableForOrbit'
	WIDGET_Selector_IRI16_variableForOrbit.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_Selector_IRI16_variableForOrbit,)
	IRI16_Panel.children += (InputsPanel,)
	# Create widget for the moudle black-box body
	Selector_Btn = widgets.Button (
		description='Selector',
		tooltip="Functions as a single point which collects all user input and distributes it to the various modules of this group",
	)
	Selector_Btn.layout.min_width = '200px'
	Selector_Btn.style.button_color = 'gold'
	IRI16_Panel.children += (Selector_Btn,)
	# Create widgets for module's outputs
	OutputsPanel = widgets.VBox()
	OutputsPanel.layout.min_width = '300px'
	WIDGET_Selector_MeshgridCSVfilenameOUT = widgets.Label(value='  --> MeshgridCSVfilenameOUT  ')
	WIDGET_Selector_MeshgridCSVfilenameOUT.layout.border = '1px dashed green'
	WIDGET_Selector_MeshgridCSVfilenameOUT.layout.margin ='0px 40px 0px 0px' 
	WIDGET_Selector_MeshgridCSVfilenameOUT.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_Selector_MeshgridCSVfilenameOUT,)
	WIDGET_Selector_MeshgridDatetimeOUT = widgets.Label(value='  --> MeshgridDatetimeOUT  ')
	WIDGET_Selector_MeshgridDatetimeOUT.layout.border = '1px dashed green'
	WIDGET_Selector_MeshgridDatetimeOUT.layout.margin ='0px 40px 0px 0px' 
	WIDGET_Selector_MeshgridDatetimeOUT.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_Selector_MeshgridDatetimeOUT,)
	WIDGET_Selector_MeshgridAltitudeOUT = widgets.Label(value='  --> MeshgridAltitudeOUT  ')
	WIDGET_Selector_MeshgridAltitudeOUT.layout.border = '1px dashed green'
	WIDGET_Selector_MeshgridAltitudeOUT.layout.margin ='0px 40px 0px 0px' 
	WIDGET_Selector_MeshgridAltitudeOUT.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_Selector_MeshgridAltitudeOUT,)
	WIDGET_Selector_MeshgridLatitudeStepOUT = widgets.Label(value='  --> MeshgridLatitudeStepOUT  ')
	WIDGET_Selector_MeshgridLatitudeStepOUT.layout.border = '1px dashed green'
	WIDGET_Selector_MeshgridLatitudeStepOUT.layout.margin ='0px 40px 0px 0px' 
	WIDGET_Selector_MeshgridLatitudeStepOUT.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_Selector_MeshgridLatitudeStepOUT,)
	WIDGET_Selector_MeshgridLongitudeStepOUT = widgets.Label(value='  --> MeshgridLongitudeStepOUT  ')
	WIDGET_Selector_MeshgridLongitudeStepOUT.layout.border = '1px dashed green'
	WIDGET_Selector_MeshgridLongitudeStepOUT.layout.margin ='0px 40px 0px 0px' 
	WIDGET_Selector_MeshgridLongitudeStepOUT.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_Selector_MeshgridLongitudeStepOUT,)
	WIDGET_Selector_NoiseWxOUT = widgets.Label(value='  --> NoiseWxOUT  ')
	WIDGET_Selector_NoiseWxOUT.layout.border = '1px dashed green'
	WIDGET_Selector_NoiseWxOUT.layout.margin ='0px 40px 0px 0px' 
	WIDGET_Selector_NoiseWxOUT.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_Selector_NoiseWxOUT,)
	WIDGET_Selector_NoiseWyOUT = widgets.Label(value='  --> NoiseWyOUT  ')
	WIDGET_Selector_NoiseWyOUT.layout.border = '1px dashed green'
	WIDGET_Selector_NoiseWyOUT.layout.margin ='0px 40px 0px 0px' 
	WIDGET_Selector_NoiseWyOUT.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_Selector_NoiseWyOUT,)
	WIDGET_Selector_NoiseAxOUT = widgets.Label(value='  --> NoiseAxOUT  ')
	WIDGET_Selector_NoiseAxOUT.layout.border = '1px dashed green'
	WIDGET_Selector_NoiseAxOUT.layout.margin ='0px 40px 0px 0px' 
	WIDGET_Selector_NoiseAxOUT.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_Selector_NoiseAxOUT,)
	WIDGET_Selector_NoiseAyOUT = widgets.Label(value='  --> NoiseAyOUT  ')
	WIDGET_Selector_NoiseAyOUT.layout.border = '1px dashed green'
	WIDGET_Selector_NoiseAyOUT.layout.margin ='0px 40px 0px 0px' 
	WIDGET_Selector_NoiseAyOUT.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_Selector_NoiseAyOUT,)
	WIDGET_Selector_IRI16_variableForMeshgrid_OUT = widgets.Label(value='  --> IRI16_variableForMeshgrid_OUT  ')
	WIDGET_Selector_IRI16_variableForMeshgrid_OUT.layout.border = '1px dashed green'
	WIDGET_Selector_IRI16_variableForMeshgrid_OUT.layout.margin ='0px 40px 0px 0px' 
	WIDGET_Selector_IRI16_variableForMeshgrid_OUT.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_Selector_IRI16_variableForMeshgrid_OUT,)
	WIDGET_Selector_IRI16_variableForOrbit_OUT = widgets.Label(value='  --> IRI16_variableForOrbit_OUT  ')
	WIDGET_Selector_IRI16_variableForOrbit_OUT.layout.border = '1px dashed green'
	WIDGET_Selector_IRI16_variableForOrbit_OUT.layout.margin ='0px 40px 0px 0px' 
	WIDGET_Selector_IRI16_variableForOrbit_OUT.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_Selector_IRI16_variableForOrbit_OUT,)
	IRI16_Panel.children += (OutputsPanel,)
	########
	## GUI code for module 'CreateCSV_Sphere'
	########
	# Create widgets for module's inputs
	InputsPanel = widgets.VBox()
	InputsPanel.layout.min_width = '330px'
	WIDGET_CreateCSV_Sphere_CSVfilename = widgets.Label(value='  Selector.MeshgridCSVfilenameOUT --> CSVfilename  ')
	WIDGET_CreateCSV_Sphere_CSVfilename.layout.border = '1px dashed blue'
	WIDGET_CreateCSV_Sphere_CSVfilename.layout.padding = '0px 10px 0px 10px'
	InputsPanel.children += (WIDGET_CreateCSV_Sphere_CSVfilename,)
	WIDGET_CreateCSV_Sphere_fixedDatetimeString = widgets.Label(value='  Selector.MeshgridDatetimeOUT --> fixedDatetimeString  ')
	WIDGET_CreateCSV_Sphere_fixedDatetimeString.layout.border = '1px dashed blue'
	WIDGET_CreateCSV_Sphere_fixedDatetimeString.layout.padding = '0px 10px 0px 10px'
	InputsPanel.children += (WIDGET_CreateCSV_Sphere_fixedDatetimeString,)
	WIDGET_CreateCSV_Sphere_fixedAltitude = widgets.Label(value='  Selector.MeshgridAltitudeOUT --> fixedAltitude  ')
	WIDGET_CreateCSV_Sphere_fixedAltitude.layout.border = '1px dashed blue'
	WIDGET_CreateCSV_Sphere_fixedAltitude.layout.padding = '0px 10px 0px 10px'
	InputsPanel.children += (WIDGET_CreateCSV_Sphere_fixedAltitude,)
	WIDGET_CreateCSV_Sphere_LatitudeStep = widgets.Label(value='  Selector.MeshgridLatitudeStepOUT --> LatitudeStep  ')
	WIDGET_CreateCSV_Sphere_LatitudeStep.layout.border = '1px dashed blue'
	WIDGET_CreateCSV_Sphere_LatitudeStep.layout.padding = '0px 10px 0px 10px'
	InputsPanel.children += (WIDGET_CreateCSV_Sphere_LatitudeStep,)
	WIDGET_CreateCSV_Sphere_LongitudeStep = widgets.Label(value='  Selector.MeshgridLongitudeStepOUT --> LongitudeStep  ')
	WIDGET_CreateCSV_Sphere_LongitudeStep.layout.border = '1px dashed blue'
	WIDGET_CreateCSV_Sphere_LongitudeStep.layout.padding = '0px 10px 0px 10px'
	InputsPanel.children += (WIDGET_CreateCSV_Sphere_LongitudeStep,)
	IRI16_Panel.children += (InputsPanel,)
	# Create widget for the moudle black-box body
	CreateCSV_Sphere_Btn = widgets.Button (
		description='CreateCSV_Sphere',
		tooltip="Creates a CSV file with the format [Time, Lat, Lon, Alt]\nTime and Alt are fixed values given as arguments and Lat and Lon are produced for the whole globe accodring to steps specified.\nCSVfilename: string, the csv filename to be written (overwrite mode)\nfixedDatetimeString: string, the date-time fixed text with format as example: Mar 17 2015 00:12:00.000. Fixed value for every csv entry.\nfixedAltitude: float positive, kilometers above earth. Fixed value for every csv entry.\nLatitudeStep: float positive, Latiude values range between 87.5 and -88.5. The step is how Lat should be incremented at each iteration.\nLongitudeStep: float positive, Longitude values range between -180 and 180.  The step is how Lon should be incremented at each iteration.\nRETURNS: the filename it has written (same as the argument)\n         the altitude (same as the argument)\n         the number of lines written",
	)
	CreateCSV_Sphere_Btn.layout.min_width = '200px'
	CreateCSV_Sphere_Btn.style.button_color = 'gold'
	IRI16_Panel.children += (CreateCSV_Sphere_Btn,)
	# Create widgets for module's outputs
	OutputsPanel = widgets.VBox()
	OutputsPanel.layout.min_width = '300px'
	WIDGET_CreateCSV_Sphere_theCSVfilename = widgets.Label(value='  --> theCSVfilename  ')
	WIDGET_CreateCSV_Sphere_theCSVfilename.layout.border = '1px dashed green'
	WIDGET_CreateCSV_Sphere_theCSVfilename.layout.margin ='0px 40px 0px 0px' 
	WIDGET_CreateCSV_Sphere_theCSVfilename.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_CreateCSV_Sphere_theCSVfilename,)
	WIDGET_CreateCSV_Sphere_theAltitude = widgets.Label(value='  --> theAltitude  ')
	WIDGET_CreateCSV_Sphere_theAltitude.layout.border = '1px dashed green'
	WIDGET_CreateCSV_Sphere_theAltitude.layout.margin ='0px 40px 0px 0px' 
	WIDGET_CreateCSV_Sphere_theAltitude.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_CreateCSV_Sphere_theAltitude,)
	WIDGET_CreateCSV_Sphere_NumberOfLinesWritten = widgets.Label(value='  --> NumberOfLinesWritten  ')
	WIDGET_CreateCSV_Sphere_NumberOfLinesWritten.layout.border = '1px dashed green'
	WIDGET_CreateCSV_Sphere_NumberOfLinesWritten.layout.margin ='0px 40px 0px 0px' 
	WIDGET_CreateCSV_Sphere_NumberOfLinesWritten.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_CreateCSV_Sphere_NumberOfLinesWritten,)
	IRI16_Panel.children += (OutputsPanel,)
	########
	## GUI code for module 'IRI16'
	########
	# Create widgets for module's inputs
	InputsPanel = widgets.VBox()
	InputsPanel.layout.min_width = '330px'
	WIDGET_IRI16_FilenameCSV_forMeshgrid = widgets.Label(value='  CreateCSV_Sphere.theCSVfilename --> FilenameCSV_forMeshgrid  ')
	WIDGET_IRI16_FilenameCSV_forMeshgrid.layout.border = '1px dashed blue'
	WIDGET_IRI16_FilenameCSV_forMeshgrid.layout.padding = '0px 10px 0px 10px'
	InputsPanel.children += (WIDGET_IRI16_FilenameCSV_forMeshgrid,)
	WIDGET_IRI16_Variable_forMeshgrid = widgets.Label(value='  Selector.IRI16_variableForMeshgrid_OUT --> Variable_forMeshgrid  ')
	WIDGET_IRI16_Variable_forMeshgrid.layout.border = '1px dashed blue'
	WIDGET_IRI16_Variable_forMeshgrid.layout.padding = '0px 10px 0px 10px'
	InputsPanel.children += (WIDGET_IRI16_Variable_forMeshgrid,)
	IRI16_Panel.children += (InputsPanel,)
	# Create widget for the moudle black-box body
	IRI16_Btn = widgets.Button (
		description='IRI16',
		tooltip="International Reference Ionosphere empirical model (IRI)\nfilename: the CSV file to read. Format: Time,Alt,Lon,Lat. The file defines the coordinates for which the IRI calculation will take place.\nParameter: selects which values to calculate. Possible values are:\n    ''   for all\n    'Op'  for Oxygen+\n    'O2p' for Oxygen2+\n    'Te' for Electron Temperature\n    'Ti' for Ion Temperature\nReturns: the name of the CSV file which contains the calculated values. Format: Time,Alt,Lon,Lat,Op,O2p,Te,Ti",
	)
	IRI16_Btn.layout.min_width = '200px'
	IRI16_Btn.style.button_color = 'salmon'
	IRI16_Panel.children += (IRI16_Btn,)
	# Create widgets for module's outputs
	OutputsPanel = widgets.VBox()
	OutputsPanel.layout.min_width = '300px'
	WIDGET_IRI16_MeshgridResultCSV = widgets.Label(value='  --> MeshgridResultCSV  ')
	WIDGET_IRI16_MeshgridResultCSV.layout.border = '1px dashed green'
	WIDGET_IRI16_MeshgridResultCSV.layout.margin ='0px 40px 0px 0px' 
	WIDGET_IRI16_MeshgridResultCSV.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_IRI16_MeshgridResultCSV,)
	IRI16_Panel.children += (OutputsPanel,)
	########
	## GUI code for module 'IRI16'
	########
	# Create widgets for module's inputs
	InputsPanel = widgets.VBox()
	InputsPanel.layout.min_width = '330px'
	WIDGET_IRI16_CSVfilenameOrbit = widgets.Label(value='  OrbitSelector.OrbitCSVfilename --> CSVfilenameOrbit  ')
	WIDGET_IRI16_CSVfilenameOrbit.layout.border = '1px dashed blue'
	WIDGET_IRI16_CSVfilenameOrbit.layout.padding = '0px 10px 0px 10px'
	InputsPanel.children += (WIDGET_IRI16_CSVfilenameOrbit,)
	WIDGET_IRI16_VariableOrbit = widgets.Label(value='  Selector.IRI16_variableForOrbit_OUT --> VariableOrbit  ')
	WIDGET_IRI16_VariableOrbit.layout.border = '1px dashed blue'
	WIDGET_IRI16_VariableOrbit.layout.padding = '0px 10px 0px 10px'
	InputsPanel.children += (WIDGET_IRI16_VariableOrbit,)
	IRI16_Panel.children += (InputsPanel,)
	# Create widget for the moudle black-box body
	IRI16_Btn = widgets.Button (
		description='IRI16',
		tooltip="International Reference Ionosphere empirical model (IRI)\nfilename: the CSV file to read. Format: Time,Alt,Lon,Lat. The file defines the coordinates for which the IRI calculation will take place.\nParameter: selects which values to calculate. Possible values are:\n    ''   for all\n    'Op'  for Oxygen+\n    'O2p' for Oxygen2+\n    'Te' for Electron Temperature\n    'Ti' for Ion Temperature\nReturns: the name of the CSV file which contains the calculated values. Format: Time,Alt,Lon,Lat,Op,O2p,Te,Ti",
	)
	IRI16_Btn.layout.min_width = '200px'
	IRI16_Btn.style.button_color = 'salmon'
	IRI16_Panel.children += (IRI16_Btn,)
	# Create widgets for module's outputs
	OutputsPanel = widgets.VBox()
	OutputsPanel.layout.min_width = '300px'
	WIDGET_IRI16_orbitCSVfilename = widgets.Label(value='  --> orbitCSVfilename  ')
	WIDGET_IRI16_orbitCSVfilename.layout.border = '1px dashed green'
	WIDGET_IRI16_orbitCSVfilename.layout.margin ='0px 40px 0px 0px' 
	WIDGET_IRI16_orbitCSVfilename.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_IRI16_orbitCSVfilename,)
	IRI16_Panel.children += (OutputsPanel,)
	########
	## GUI code for module 'SubGridVariability'
	########
	# Create widgets for module's inputs
	InputsPanel = widgets.VBox()
	InputsPanel.layout.min_width = '330px'
	WIDGET_SubGridVariability_NoiseWx = widgets.Label(value='  Selector.NoiseWxOUT --> NoiseWx  ')
	WIDGET_SubGridVariability_NoiseWx.layout.border = '1px dashed blue'
	WIDGET_SubGridVariability_NoiseWx.layout.padding = '0px 10px 0px 10px'
	InputsPanel.children += (WIDGET_SubGridVariability_NoiseWx,)
	WIDGET_SubGridVariability_NoiseWy = widgets.Label(value='  Selector.NoiseWyOUT --> NoiseWy  ')
	WIDGET_SubGridVariability_NoiseWy.layout.border = '1px dashed blue'
	WIDGET_SubGridVariability_NoiseWy.layout.padding = '0px 10px 0px 10px'
	InputsPanel.children += (WIDGET_SubGridVariability_NoiseWy,)
	WIDGET_SubGridVariability_NoiseAx = widgets.Label(value='  Selector.NoiseAxOUT --> NoiseAx  ')
	WIDGET_SubGridVariability_NoiseAx.layout.border = '1px dashed blue'
	WIDGET_SubGridVariability_NoiseAx.layout.padding = '0px 10px 0px 10px'
	InputsPanel.children += (WIDGET_SubGridVariability_NoiseAx,)
	WIDGET_SubGridVariability_NoiseAy = widgets.Label(value='  Selector.NoiseAyOUT --> NoiseAy  ')
	WIDGET_SubGridVariability_NoiseAy.layout.border = '1px dashed blue'
	WIDGET_SubGridVariability_NoiseAy.layout.padding = '0px 10px 0px 10px'
	InputsPanel.children += (WIDGET_SubGridVariability_NoiseAy,)
	WIDGET_SubGridVariability_filename = widgets.Label(value='  IRI16.orbitCSVfilename --> filename  ')
	WIDGET_SubGridVariability_filename.layout.border = '1px dashed blue'
	WIDGET_SubGridVariability_filename.layout.padding = '0px 10px 0px 10px'
	InputsPanel.children += (WIDGET_SubGridVariability_filename,)
	global WIDGET_SubGridVariability_ValueName
	WIDGET_SubGridVariability_ValueName = widgets.Text()
	WIDGET_SubGridVariability_ValueName.description = 'ValueName'
	WIDGET_SubGridVariability_ValueName.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_SubGridVariability_ValueName,)
	IRI16_Panel.children += (InputsPanel,)
	# Create widget for the moudle black-box body
	SubGridVariability_Btn = widgets.Button (
		description='SubGridVariability',
		tooltip="Inputs\n NoiseAx,NoiseAy: real:: Amplitude of noise Kind (values 0-100)\n NoiseWx,NoiseWy: real:: Wavenumber of noise Kind (values 0-1000)\n Filename: Character:: field to read orbit from\n ValueName: Character:: variable from model data to use and add SGV. Corresponds to CSV Column\nOutputs\n Plots..\nReturns\n a csv filename of format Time, Lat, Lon, Alt, value containing the altered-with-noise values",
	)
	SubGridVariability_Btn.layout.min_width = '200px'
	SubGridVariability_Btn.style.button_color = 'gold'
	IRI16_Panel.children += (SubGridVariability_Btn,)
	# Create widgets for module's outputs
	OutputsPanel = widgets.VBox()
	OutputsPanel.layout.min_width = '300px'
	WIDGET_SubGridVariability_VariabilityCSVfilename = widgets.Label(value='  --> VariabilityCSVfilename  ')
	WIDGET_SubGridVariability_VariabilityCSVfilename.layout.border = '1px dashed green'
	WIDGET_SubGridVariability_VariabilityCSVfilename.layout.margin ='0px 40px 0px 0px' 
	WIDGET_SubGridVariability_VariabilityCSVfilename.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_SubGridVariability_VariabilityCSVfilename,)
	IRI16_Panel.children += (OutputsPanel,)
	########
	## GUI code for module 'PlotGlobe'
	########
	# Create widgets for module's inputs
	InputsPanel = widgets.VBox()
	InputsPanel.layout.min_width = '330px'
	WIDGET_PlotGlobe_SurfaceFilename = widgets.Label(value='  IRI16.MeshgridResultCSV --> SurfaceFilename  ')
	WIDGET_PlotGlobe_SurfaceFilename.layout.border = '1px dashed blue'
	WIDGET_PlotGlobe_SurfaceFilename.layout.padding = '0px 10px 0px 10px'
	InputsPanel.children += (WIDGET_PlotGlobe_SurfaceFilename,)
	global WIDGET_PlotGlobe_SurfaceVariableToPlot
	WIDGET_PlotGlobe_SurfaceVariableToPlot = widgets.Text(value='ALFA')
	WIDGET_PlotGlobe_SurfaceVariableToPlot.description = 'SurfaceVariableToPlot'
	WIDGET_PlotGlobe_SurfaceVariableToPlot.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_PlotGlobe_SurfaceVariableToPlot,)
	global WIDGET_PlotGlobe_SurfaceColorbarTitle
	WIDGET_PlotGlobe_SurfaceColorbarTitle = widgets.Text()
	WIDGET_PlotGlobe_SurfaceColorbarTitle.description = 'SurfaceColorbarTitle'
	WIDGET_PlotGlobe_SurfaceColorbarTitle.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_PlotGlobe_SurfaceColorbarTitle,)
	global WIDGET_PlotGlobe_SurfaceColorscaleName
	WIDGET_PlotGlobe_SurfaceColorscaleName = widgets.Dropdown( options=['', 'None', 'Blackbody', 'Bluered', 'Blues', 'Earth', 'Electric', 'Greens', 'Greys', 'Hot', 'Jet', 'Picnic', 'Portland', 'Rainbow', 'RdBu', 'Reds', 'Viridis', 'YlGnBu', 'YlOrRd'], value='Jet', description='Variable')
	WIDGET_PlotGlobe_SurfaceColorscaleName.description = 'SurfaceColorscaleName'
	WIDGET_PlotGlobe_SurfaceColorscaleName.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_PlotGlobe_SurfaceColorscaleName,)
	WIDGET_PlotGlobe_OrbitFilename = widgets.Label(value='  SubGridVariability.VariabilityCSVfilename --> OrbitFilename  ')
	WIDGET_PlotGlobe_OrbitFilename.layout.border = '1px dashed blue'
	WIDGET_PlotGlobe_OrbitFilename.layout.padding = '0px 10px 0px 10px'
	InputsPanel.children += (WIDGET_PlotGlobe_OrbitFilename,)
	global WIDGET_PlotGlobe_OrbitVariableToPlot
	WIDGET_PlotGlobe_OrbitVariableToPlot = widgets.Text(value='ALFA')
	WIDGET_PlotGlobe_OrbitVariableToPlot.description = 'OrbitVariableToPlot'
	WIDGET_PlotGlobe_OrbitVariableToPlot.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_PlotGlobe_OrbitVariableToPlot,)
	global WIDGET_PlotGlobe_OrbitColorbarTitle
	WIDGET_PlotGlobe_OrbitColorbarTitle = widgets.Text()
	WIDGET_PlotGlobe_OrbitColorbarTitle.description = 'OrbitColorbarTitle'
	WIDGET_PlotGlobe_OrbitColorbarTitle.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_PlotGlobe_OrbitColorbarTitle,)
	global WIDGET_PlotGlobe_OrbitColorscaleName
	WIDGET_PlotGlobe_OrbitColorscaleName = widgets.Dropdown( options=['', 'None', 'Blackbody', 'Bluered', 'Blues', 'Earth', 'Electric', 'Greens', 'Greys', 'Hot', 'Jet', 'Picnic', 'Portland', 'Rainbow', 'RdBu', 'Reds', 'Viridis', 'YlGnBu', 'YlOrRd'], value='Jet', description='Variable')
	WIDGET_PlotGlobe_OrbitColorscaleName.description = 'OrbitColorscaleName'
	WIDGET_PlotGlobe_OrbitColorscaleName.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_PlotGlobe_OrbitColorscaleName,)
	global WIDGET_PlotGlobe_PlotTitle
	WIDGET_PlotGlobe_PlotTitle = widgets.Text()
	WIDGET_PlotGlobe_PlotTitle.description = 'PlotTitle'
	WIDGET_PlotGlobe_PlotTitle.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_PlotGlobe_PlotTitle,)
	IRI16_Panel.children += (InputsPanel,)
	# Create widget for the moudle black-box body
	PlotGlobe_Btn = widgets.Button (
		description='PlotGlobe',
		tooltip="Creates a 3D plot of an earth globe. Can plot a sphere surface and/or a satellite orbit. \nThe surface and the orbit are colored according to the selected variable values in the data files. \nThe data files can be of CSV o NetCDF4 format.\nThe same or different variables can be selected for plotting at the surface and the orbit.\nThe same or different colorscales can be selected for the surface and the orbit.\nThe same or different value range can be selected for the surface and the orbit by assigning the same or different colorbar titles.\nValid Colorscale Names: ‘Blackbody’, ‘Bluered’, ‘Blues’, ‘Earth’, ‘Electric’, ‘Greens’, ‘Greys’, ‘Hot’, ‘Jet’, ‘Picnic’, ‘Portland’, ‘Rainbow’, ‘RdBu’, ‘Reds’, ‘Viridis’, ‘YlGnBu’, ‘YlOrRd’\nARGUMENTS:\n  SurfaceFilename: The file which contains data for a surface. If empty string then no surface will be plotted. \n                   CSV format: Time,Lat,Lon,Alt,value. Contains the data for the sphere surface. \n                   NetCDF format: see Daedalus User Manual\n  SurfaceVariableToPlot: The name of the variable to read from the surface data file and plot upon the globe.\n  SurfaceColorbarTitle: A title to display above the colorbar which refers to the surface data\n  SurfaceColorscaleName: The name of the Colorscale to use for the surface data. \n                         In case an empty string is passed then a default HeatMap colorscale will be applied.\n                         In case None is passed then all points will be black irrespective of value.\n  OrbitFilename: counterpart of the corresponding argument for the Surface data\n  OrbitVariableToPlot: counterpart of the corresponding argument for the Surface data \n  OrbitColorbarTitle: counterpart of the corresponding argument for the Surface data \n  OrbitColorscaleName: counterpart of the corresponding argument for the Surface data\n  PlotTitle: A title which is displayed on the top of the plot.\nRETURNS: a string containing information about the Data\n",
	)
	PlotGlobe_Btn.layout.min_width = '200px'
	PlotGlobe_Btn.style.button_color = 'deeppink'
	IRI16_Panel.children += (PlotGlobe_Btn,)
	# Create widgets for module's outputs
	OutputsPanel = widgets.VBox()
	OutputsPanel.layout.min_width = '300px'
	IRI16_Panel.children += (OutputsPanel,)
	return IRI16_Panel


def Construct_TopLevel():
	#Create Containers
	TopLevel_Panel = widgets.Box()
	TopLevel_Panel.layout.overflow_x = 'scroll'
	TopLevel_Panel.layout.flex_flow = 'row'
	TopLevel_Panel.layout.display = 'flex'
	########
	## GUI code for module 'OrbitSelector'
	########
	# Create widgets for module's inputs
	InputsPanel = widgets.VBox()
	InputsPanel.layout.min_width = '330px'
	global WIDGET_OrbitSelector_Filename
	WIDGET_OrbitSelector_Filename = widgets.Text("DAED_ORB_Evt0_LLA_Per150_Lat80_Srt01Hz_Msc.nc")
	WIDGET_OrbitSelector_Filename.description = 'Filename'
	WIDGET_OrbitSelector_Filename.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_OrbitSelector_Filename,)
	global WIDGET_OrbitSelector_EvtXY
	WIDGET_OrbitSelector_EvtXY = widgets.Dropdown( options=['Evt0', 'Evt1', 'Evt2', 'Evt3', 'Evt4', 'Evt5', 'Evt6', 'Evt7', 'Evt8', 'Evt9'], description='Variable')
	WIDGET_OrbitSelector_EvtXY.description = 'EvtXY'
	WIDGET_OrbitSelector_EvtXY.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_OrbitSelector_EvtXY,)
	global WIDGET_OrbitSelector_TYP
	WIDGET_OrbitSelector_TYP = widgets.Dropdown( options=['LLA', 'VEL', 'PTG'], description='Variable')
	WIDGET_OrbitSelector_TYP.description = 'TYP'
	WIDGET_OrbitSelector_TYP.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_OrbitSelector_TYP,)
	global WIDGET_OrbitSelector_PerYYY
	WIDGET_OrbitSelector_PerYYY = widgets.Dropdown( options=['Per120', 'Per150'], description='Variable')
	WIDGET_OrbitSelector_PerYYY.description = 'PerYYY'
	WIDGET_OrbitSelector_PerYYY.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_OrbitSelector_PerYYY,)
	global WIDGET_OrbitSelector_LatZZ
	WIDGET_OrbitSelector_LatZZ = widgets.Dropdown( options=['Lat00', 'Lat40', 'Lat80'],  description='Variable')
	WIDGET_OrbitSelector_LatZZ.description = 'LatZZ'
	WIDGET_OrbitSelector_LatZZ.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_OrbitSelector_LatZZ,)
	global WIDGET_OrbitSelector_SRXXHZ
	WIDGET_OrbitSelector_SRXXHZ = widgets.Dropdown( options=['Srt16Hz', 'Srt01Hz'],  description='Variable')
	WIDGET_OrbitSelector_SRXXHZ.description = 'SRXXHZ'
	WIDGET_OrbitSelector_SRXXHZ.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_OrbitSelector_SRXXHZ,)
	global WIDGET_OrbitSelector_SC
	WIDGET_OrbitSelector_SC = widgets.Dropdown( options=['Msc', 'Ssc'],  description='Variable')
	WIDGET_OrbitSelector_SC.description = 'SC'
	WIDGET_OrbitSelector_SC.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_OrbitSelector_SC,)
	TopLevel_Panel.children += (InputsPanel,)
	# Create widget for the moudle black-box body
	OrbitSelector_Btn = widgets.Button (
		description='OrbitSelector',
		tooltip="Constructs and returns an orbit full-path filename from the orbit properties.\nThe filename includes all parameters required to select the appropriate orbit.\nAll filenames  have the same length and same number of parameters. \nIn case OrbitFilename and EvtXY are empty then an empty string is returned.\nFILE FORMAT: DAED_ORB_EVTXY_TYP_PERYYY_LATZZ_SRTQQHz_XSC.csv\n  Parameter OrbitFilename:\n    If it is not empty then the rest arguments are ignored. \n    If it contains slashes then is is assumed it is a full path name and it is returned as it is.\n    If it does not contain slashes then the orbits path is added.\n    \n  Parameter EvtXY values:\n    EVTXS    X Event, Single Orbit\n    EVTXA    X Event, All Orbit\n    EVT1Y    1st Event: St Patrick’s day event [17 Mar 2015 – 20 Mar 2015]\n    EVT2Y    2nd Event ...  ...\n \n  Parameter TYP values:\n    LLA    Time,Latitude Longitude Altitude \n    VEL    Time,VmagVxVyVz\n    PTG    Time, X_GSE, Y_GSE, Z_GSE, RamX_GSE, RamY_GSE, RamZ_GSE\n  Parameter PerYYY values:\n    PER120, PER150    Perigee Altitude at 120km or 150km\n  Parameter LatZZ values:\n    LAT00, LAT40, LAT80    Perigee Latitude at 0°, 40° or 80°\n \n  Parameter SRXXHZ values:\n    SRT16Hz, SRT01Hz    Sampling rate 16Hz or 1HZ\n \n  Parameter SC values:\n    MSC    Mother Spacecraft\n    SSC    Sub-Spacecraft\n",
	)
	OrbitSelector_Btn.layout.min_width = '200px'
	OrbitSelector_Btn.style.button_color = 'YellowGreen'
	TopLevel_Panel.children += (OrbitSelector_Btn,)
	# Create widgets for module's outputs
	OutputsPanel = widgets.VBox()
	OutputsPanel.layout.min_width = '300px'
	WIDGET_OrbitSelector_OrbitCSVfilename = widgets.Label(value='  --> OrbitCSVfilename  ')
	WIDGET_OrbitSelector_OrbitCSVfilename.layout.border = '1px dashed green'
	WIDGET_OrbitSelector_OrbitCSVfilename.layout.margin ='0px 40px 0px 0px' 
	WIDGET_OrbitSelector_OrbitCSVfilename.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_OrbitSelector_OrbitCSVfilename,)
	WIDGET_OrbitSelector_OrbitNetCDFfilename = widgets.Label(value='  --> OrbitNetCDFfilename  ')
	WIDGET_OrbitSelector_OrbitNetCDFfilename.layout.border = '1px dashed green'
	WIDGET_OrbitSelector_OrbitNetCDFfilename.layout.margin ='0px 40px 0px 0px' 
	WIDGET_OrbitSelector_OrbitNetCDFfilename.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_OrbitSelector_OrbitNetCDFfilename,)
	TopLevel_Panel.children += (OutputsPanel,)
	########
	## GUI code for module 'Selector'
	########
	# Create widgets for module's inputs
	InputsPanel = widgets.VBox()
	InputsPanel.layout.min_width = '330px'
	global WIDGET_Selector_EnableSubGridVariability
	WIDGET_Selector_EnableSubGridVariability = widgets.Checkbox( value=True)
	WIDGET_Selector_EnableSubGridVariability.description = 'EnableSubGridVariability'
	WIDGET_Selector_EnableSubGridVariability.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_Selector_EnableSubGridVariability,)
	TopLevel_Panel.children += (InputsPanel,)
	# Create widget for the moudle black-box body
	Selector_Btn = widgets.Button (
		description='Selector',
		tooltip="",
	)
	Selector_Btn.layout.min_width = '200px'
	Selector_Btn.style.button_color = 'YellowGreen'
	TopLevel_Panel.children += (Selector_Btn,)
	# Create widgets for module's outputs
	OutputsPanel = widgets.VBox()
	OutputsPanel.layout.min_width = '300px'
	WIDGET_Selector_EnableSubGridVariabilityOUT = widgets.Label(value='  --> EnableSubGridVariabilityOUT  ')
	WIDGET_Selector_EnableSubGridVariabilityOUT.layout.border = '1px dashed green'
	WIDGET_Selector_EnableSubGridVariabilityOUT.layout.margin ='0px 40px 0px 0px' 
	WIDGET_Selector_EnableSubGridVariabilityOUT.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_Selector_EnableSubGridVariabilityOUT,)
	TopLevel_Panel.children += (OutputsPanel,)
	########
	## GUI code for module 'Interpolator'
	########
	# Create widgets for module's inputs
	InputsPanel = widgets.VBox()
	InputsPanel.layout.min_width = '330px'
	global WIDGET_Interpolator_ModelFile
	WIDGET_Interpolator_ModelFile = widgets.Text(value="/home/NAS/TIEGCM_DATA/TIEGCM_EVT2_2017_Sept_S/tiegcm_s_24900.nc")
	WIDGET_Interpolator_ModelFile.description = 'ModelFile'
	WIDGET_Interpolator_ModelFile.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_Interpolator_ModelFile,)
	WIDGET_Interpolator_OrbitFile = widgets.Label(value='  OrbitSelector.OrbitNetCDFfilename --> OrbitFile  ')
	WIDGET_Interpolator_OrbitFile.layout.border = '1px dashed blue'
	WIDGET_Interpolator_OrbitFile.layout.padding = '0px 10px 0px 10px'
	InputsPanel.children += (WIDGET_Interpolator_OrbitFile,)
	TopLevel_Panel.children += (InputsPanel,)
	# Create widget for the moudle black-box body
	Interpolator_Btn = widgets.Button (
		description='Interpolator',
		tooltip="model-->string:: name of model eg TIEGCM\n model_data_file--> string:: model data files stored on DaedalusNAS\n orbit_file-->string :: orbit filename in Time-Lat-Lon-Alt format stored on DaedalusNAS\n save--> Logical:: if true saves interpolated values to directory\n VAR--> string:: variable to inteprolate, must be  one of the variables included in the model data\nOutputs:: Plots+ 1 csv file stored on DaedalusNAS in ModelOutputs/Interpolator",
	)
	Interpolator_Btn.layout.min_width = '200px'
	Interpolator_Btn.style.button_color = 'gold'
	TopLevel_Panel.children += (Interpolator_Btn,)
	# Create widgets for module's outputs
	OutputsPanel = widgets.VBox()
	OutputsPanel.layout.min_width = '300px'
	WIDGET_Interpolator_InterpolatorResultFilename = widgets.Label(value='  --> InterpolatorResultFilename  ')
	WIDGET_Interpolator_InterpolatorResultFilename.layout.border = '1px dashed green'
	WIDGET_Interpolator_InterpolatorResultFilename.layout.margin ='0px 40px 0px 0px' 
	WIDGET_Interpolator_InterpolatorResultFilename.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_Interpolator_InterpolatorResultFilename,)
	TopLevel_Panel.children += (OutputsPanel,)
	########
	## GUI code for module 'WakeEffectsCalc'
	########
	# Create widgets for module's inputs
	InputsPanel = widgets.VBox()
	InputsPanel.layout.min_width = '330px'
	WIDGET_WakeEffectsCalc_InputFilename = widgets.Label(value='  SGV_E.ResultFilename --> InputFilename  ')
	WIDGET_WakeEffectsCalc_InputFilename.layout.border = '1px dashed blue'
	WIDGET_WakeEffectsCalc_InputFilename.layout.padding = '0px 10px 0px 10px'
	InputsPanel.children += (WIDGET_WakeEffectsCalc_InputFilename,)
	TopLevel_Panel.children += (InputsPanel,)
	# Create widget for the moudle black-box body
	WakeEffectsCalc_Btn = widgets.Button (
		description='WakeEffectsCalc',
		tooltip="",
	)
	WakeEffectsCalc_Btn.layout.min_width = '200px'
	WakeEffectsCalc_Btn.style.button_color = 'DeepSkyBlue'
	TopLevel_Panel.children += (WakeEffectsCalc_Btn,)
	# Create widgets for module's outputs
	OutputsPanel = widgets.VBox()
	OutputsPanel.layout.min_width = '300px'
	WIDGET_WakeEffectsCalc_ResultFilename = widgets.Label(value='  --> ResultFilename  ')
	WIDGET_WakeEffectsCalc_ResultFilename.layout.border = '1px dashed green'
	WIDGET_WakeEffectsCalc_ResultFilename.layout.margin ='0px 40px 0px 0px' 
	WIDGET_WakeEffectsCalc_ResultFilename.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_WakeEffectsCalc_ResultFilename,)
	TopLevel_Panel.children += (OutputsPanel,)
	########
	## GUI code for module 'ChargingEffectsCalc'
	########
	# Create widgets for module's inputs
	InputsPanel = widgets.VBox()
	InputsPanel.layout.min_width = '330px'
	WIDGET_ChargingEffectsCalc_InputFilename = widgets.Label(value='  SGV_E.ResultFilename --> InputFilename  ')
	WIDGET_ChargingEffectsCalc_InputFilename.layout.border = '1px dashed blue'
	WIDGET_ChargingEffectsCalc_InputFilename.layout.padding = '0px 10px 0px 10px'
	InputsPanel.children += (WIDGET_ChargingEffectsCalc_InputFilename,)
	TopLevel_Panel.children += (InputsPanel,)
	# Create widget for the moudle black-box body
	ChargingEffectsCalc_Btn = widgets.Button (
		description='ChargingEffectsCalc',
		tooltip="",
	)
	ChargingEffectsCalc_Btn.layout.min_width = '200px'
	ChargingEffectsCalc_Btn.style.button_color = 'DeepSkyBlue'
	TopLevel_Panel.children += (ChargingEffectsCalc_Btn,)
	# Create widgets for module's outputs
	OutputsPanel = widgets.VBox()
	OutputsPanel.layout.min_width = '300px'
	WIDGET_ChargingEffectsCalc_ResultFilename = widgets.Label(value='  --> ResultFilename  ')
	WIDGET_ChargingEffectsCalc_ResultFilename.layout.border = '1px dashed green'
	WIDGET_ChargingEffectsCalc_ResultFilename.layout.margin ='0px 40px 0px 0px' 
	WIDGET_ChargingEffectsCalc_ResultFilename.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_ChargingEffectsCalc_ResultFilename,)
	TopLevel_Panel.children += (OutputsPanel,)
	########
	## GUI code for module 'DerivedProductsCalc'
	########
	# Create widgets for module's inputs
	InputsPanel = widgets.VBox()
	InputsPanel.layout.min_width = '330px'
	WIDGET_DerivedProductsCalc_InputFilename = widgets.Label(value='  InstrumentGNSS.ResultFilename --> InputFilename  ')
	WIDGET_DerivedProductsCalc_InputFilename.layout.border = '1px dashed blue'
	WIDGET_DerivedProductsCalc_InputFilename.layout.padding = '0px 10px 0px 10px'
	InputsPanel.children += (WIDGET_DerivedProductsCalc_InputFilename,)
	TopLevel_Panel.children += (InputsPanel,)
	# Create widget for the moudle black-box body
	DerivedProductsCalc_Btn = widgets.Button (
		description='DerivedProductsCalc',
		tooltip="",
	)
	DerivedProductsCalc_Btn.layout.min_width = '200px'
	DerivedProductsCalc_Btn.style.button_color = 'Tomato'
	TopLevel_Panel.children += (DerivedProductsCalc_Btn,)
	# Create widgets for module's outputs
	OutputsPanel = widgets.VBox()
	OutputsPanel.layout.min_width = '300px'
	WIDGET_DerivedProductsCalc_ResultFilename = widgets.Label(value='  --> ResultFilename  ')
	WIDGET_DerivedProductsCalc_ResultFilename.layout.border = '1px dashed green'
	WIDGET_DerivedProductsCalc_ResultFilename.layout.margin ='0px 40px 0px 0px' 
	WIDGET_DerivedProductsCalc_ResultFilename.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_DerivedProductsCalc_ResultFilename,)
	TopLevel_Panel.children += (OutputsPanel,)
	return TopLevel_Panel


def Construct_Interpolator():
	#Create Containers
	Interpolator_Panel = widgets.Box()
	Interpolator_Panel.layout.overflow_x = 'scroll'
	Interpolator_Panel.layout.flex_flow = 'row'
	Interpolator_Panel.layout.display = 'flex'
	########
	## GUI code for module 'OrbitSelector'
	########
	# Create widgets for module's inputs
	InputsPanel = widgets.VBox()
	InputsPanel.layout.min_width = '330px'
	global WIDGET_OrbitSelector_Filename
	WIDGET_OrbitSelector_Filename = widgets.Text("DAED_ORB_Evt0_LLA_Per150_Lat80_Srt01Hz_Msc.csv")
	WIDGET_OrbitSelector_Filename.description = 'Filename'
	WIDGET_OrbitSelector_Filename.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_OrbitSelector_Filename,)
	global WIDGET_OrbitSelector_EvtXY
	WIDGET_OrbitSelector_EvtXY = widgets.Dropdown( options=['Evt0', 'Evt1', 'Evt2', 'Evt3', 'Evt4', 'Evt5', 'Evt6', 'Evt7', 'Evt8', 'Evt9'], description='Variable')
	WIDGET_OrbitSelector_EvtXY.description = 'EvtXY'
	WIDGET_OrbitSelector_EvtXY.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_OrbitSelector_EvtXY,)
	global WIDGET_OrbitSelector_TYP
	WIDGET_OrbitSelector_TYP = widgets.Dropdown( options=['LLA', 'VEL', 'PTG'], description='Variable')
	WIDGET_OrbitSelector_TYP.description = 'TYP'
	WIDGET_OrbitSelector_TYP.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_OrbitSelector_TYP,)
	global WIDGET_OrbitSelector_PerYYY
	WIDGET_OrbitSelector_PerYYY = widgets.Dropdown( options=['Per120', 'Per150'], description='Variable')
	WIDGET_OrbitSelector_PerYYY.description = 'PerYYY'
	WIDGET_OrbitSelector_PerYYY.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_OrbitSelector_PerYYY,)
	global WIDGET_OrbitSelector_LatZZ
	WIDGET_OrbitSelector_LatZZ = widgets.Dropdown( options=['Lat00', 'Lat40', 'Lat80'],  description='Variable')
	WIDGET_OrbitSelector_LatZZ.description = 'LatZZ'
	WIDGET_OrbitSelector_LatZZ.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_OrbitSelector_LatZZ,)
	global WIDGET_OrbitSelector_SRXXHZ
	WIDGET_OrbitSelector_SRXXHZ = widgets.Dropdown( options=['Srt16Hz', 'Srt01Hz'],  description='Variable')
	WIDGET_OrbitSelector_SRXXHZ.description = 'SRXXHZ'
	WIDGET_OrbitSelector_SRXXHZ.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_OrbitSelector_SRXXHZ,)
	global WIDGET_OrbitSelector_SC
	WIDGET_OrbitSelector_SC = widgets.Dropdown( options=['Msc', 'Ssc'],  description='Variable')
	WIDGET_OrbitSelector_SC.description = 'SC'
	WIDGET_OrbitSelector_SC.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_OrbitSelector_SC,)
	Interpolator_Panel.children += (InputsPanel,)
	# Create widget for the moudle black-box body
	OrbitSelector_Btn = widgets.Button (
		description='OrbitSelector',
		tooltip="Constructs and returns an orbit full-path filename from the orbit properties.\nThe filename includes all parameters required to select the appropriate orbit.\nAll filenames  have the same length and same number of parameters. \nIn case OrbitFilename and EvtXY are empty then an empty string is returned.\nFILE FORMAT: DAED_ORB_EVTXY_TYP_PERYYY_LATZZ_SRTQQHz_XSC.csv\n  Parameter OrbitFilename:\n    If it is not empty then the rest arguments are ignored. \n    If it contains slashes then is is assumed it is a full path name and it is returned as it is.\n    If it does not contain slashes then the orbits path is added.\n    \n  Parameter EvtXY values:\n    EVTXS    X Event, Single Orbit\n    EVTXA    X Event, All Orbit\n    EVT1Y    1st Event: St Patrick’s day event [17 Mar 2015 – 20 Mar 2015]\n    EVT2Y    2nd Event ...  ...\n \n  Parameter TYP values:\n    LLA    Time,Latitude Longitude Altitude \n    VEL    Time,VmagVxVyVz\n    PTG    Time, X_GSE, Y_GSE, Z_GSE, RamX_GSE, RamY_GSE, RamZ_GSE\n  Parameter PerYYY values:\n    PER120, PER150    Perigee Altitude at 120km or 150km\n  Parameter LatZZ values:\n    LAT00, LAT40, LAT80    Perigee Latitude at 0°, 40° or 80°\n \n  Parameter SRXXHZ values:\n    SRT16Hz, SRT01Hz    Sampling rate 16Hz or 1HZ\n \n  Parameter SC values:\n    MSC    Mother Spacecraft\n    SSC    Sub-Spacecraft\n",
	)
	OrbitSelector_Btn.layout.min_width = '200px'
	OrbitSelector_Btn.style.button_color = 'YellowGreen'
	Interpolator_Panel.children += (OrbitSelector_Btn,)
	# Create widgets for module's outputs
	OutputsPanel = widgets.VBox()
	OutputsPanel.layout.min_width = '300px'
	WIDGET_OrbitSelector_OrbitCSVfilename = widgets.Label(value='  --> OrbitCSVfilename  ')
	WIDGET_OrbitSelector_OrbitCSVfilename.layout.border = '1px dashed green'
	WIDGET_OrbitSelector_OrbitCSVfilename.layout.margin ='0px 40px 0px 0px' 
	WIDGET_OrbitSelector_OrbitCSVfilename.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_OrbitSelector_OrbitCSVfilename,)
	WIDGET_OrbitSelector_OrbitNetCDFfilename = widgets.Label(value='  --> OrbitNetCDFfilename  ')
	WIDGET_OrbitSelector_OrbitNetCDFfilename.layout.border = '1px dashed green'
	WIDGET_OrbitSelector_OrbitNetCDFfilename.layout.margin ='0px 40px 0px 0px' 
	WIDGET_OrbitSelector_OrbitNetCDFfilename.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_OrbitSelector_OrbitNetCDFfilename,)
	Interpolator_Panel.children += (OutputsPanel,)
	########
	## GUI code for module 'Interpolator'
	########
	# Create widgets for module's inputs
	InputsPanel = widgets.VBox()
	InputsPanel.layout.min_width = '330px'
	global WIDGET_Interpolator_model
	WIDGET_Interpolator_model = widgets.Text(value="TIEGCM")
	WIDGET_Interpolator_model.description = 'model'
	WIDGET_Interpolator_model.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_Interpolator_model,)
	global WIDGET_Interpolator_model_data_file
	WIDGET_Interpolator_model_data_file = widgets.Text(value="tiegcm_s_24900.nc")
	WIDGET_Interpolator_model_data_file.description = 'model_data_file'
	WIDGET_Interpolator_model_data_file.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_Interpolator_model_data_file,)
	WIDGET_Interpolator_orbit_file = widgets.Label(value='  OrbitSelector.OrbitCSVfilename --> orbit_file  ')
	WIDGET_Interpolator_orbit_file.layout.border = '1px dashed blue'
	WIDGET_Interpolator_orbit_file.layout.padding = '0px 10px 0px 10px'
	InputsPanel.children += (WIDGET_Interpolator_orbit_file,)
	global WIDGET_Interpolator_save
	WIDGET_Interpolator_save = widgets.Checkbox( value=True, description='Save')
	WIDGET_Interpolator_save.description = 'save'
	WIDGET_Interpolator_save.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_Interpolator_save,)
	global WIDGET_Interpolator_VAR
	WIDGET_Interpolator_VAR = widgets.Text(value="O2")
	WIDGET_Interpolator_VAR.description = 'VAR'
	WIDGET_Interpolator_VAR.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_Interpolator_VAR,)
	Interpolator_Panel.children += (InputsPanel,)
	# Create widget for the moudle black-box body
	Interpolator_Btn = widgets.Button (
		description='Interpolator',
		tooltip="model-->string:: name of model eg TIEGCM\n model_data_file--> string:: model data files stored on DaedalusNAS\n orbit_file-->string :: orbit filename in Time-Lat-Lon-Alt format stored on DaedalusNAS\n save--> Logical:: if true saves interpolated values to directory\n VAR--> string:: variable to inteprolate, must be  one of the variables included in the model data\nOutputs:: Plots+ 1 csv file stored on DaedalusNAS in ModelOutputs/Interpolator",
	)
	Interpolator_Btn.layout.min_width = '200px'
	Interpolator_Btn.style.button_color = 'gold'
	Interpolator_Panel.children += (Interpolator_Btn,)
	# Create widgets for module's outputs
	OutputsPanel = widgets.VBox()
	OutputsPanel.layout.min_width = '300px'
	WIDGET_Interpolator_InterpolatedOrbitCSV = widgets.Label(value='  --> InterpolatedOrbitCSV  ')
	WIDGET_Interpolator_InterpolatedOrbitCSV.layout.border = '1px dashed green'
	WIDGET_Interpolator_InterpolatedOrbitCSV.layout.margin ='0px 40px 0px 0px' 
	WIDGET_Interpolator_InterpolatedOrbitCSV.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_Interpolator_InterpolatedOrbitCSV,)
	Interpolator_Panel.children += (OutputsPanel,)
	########
	## GUI code for module 'PlotGlobe'
	########
	# Create widgets for module's inputs
	InputsPanel = widgets.VBox()
	InputsPanel.layout.min_width = '330px'
	global WIDGET_PlotGlobe_SurfaceFilename
	WIDGET_PlotGlobe_SurfaceFilename = widgets.Text()
	WIDGET_PlotGlobe_SurfaceFilename.description = 'SurfaceFilename'
	WIDGET_PlotGlobe_SurfaceFilename.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_PlotGlobe_SurfaceFilename,)
	global WIDGET_PlotGlobe_SurfaceVariableToPlot
	WIDGET_PlotGlobe_SurfaceVariableToPlot = widgets.Text(value='ALFA')
	WIDGET_PlotGlobe_SurfaceVariableToPlot.description = 'SurfaceVariableToPlot'
	WIDGET_PlotGlobe_SurfaceVariableToPlot.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_PlotGlobe_SurfaceVariableToPlot,)
	global WIDGET_PlotGlobe_SurfaceColorbarTitle
	WIDGET_PlotGlobe_SurfaceColorbarTitle = widgets.Text()
	WIDGET_PlotGlobe_SurfaceColorbarTitle.description = 'SurfaceColorbarTitle'
	WIDGET_PlotGlobe_SurfaceColorbarTitle.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_PlotGlobe_SurfaceColorbarTitle,)
	global WIDGET_PlotGlobe_SurfaceColorscaleName
	WIDGET_PlotGlobe_SurfaceColorscaleName = widgets.Dropdown( options=['', 'None', 'Blackbody', 'Bluered', 'Blues', 'Earth', 'Electric', 'Greens', 'Greys', 'Hot', 'Jet', 'Picnic', 'Portland', 'Rainbow', 'RdBu', 'Reds', 'Viridis', 'YlGnBu', 'YlOrRd'], value='Jet', description='Variable')
	WIDGET_PlotGlobe_SurfaceColorscaleName.description = 'SurfaceColorscaleName'
	WIDGET_PlotGlobe_SurfaceColorscaleName.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_PlotGlobe_SurfaceColorscaleName,)
	WIDGET_PlotGlobe_OrbitFilename = widgets.Label(value='  Interpolator.InterpolatedOrbitCSV --> OrbitFilename  ')
	WIDGET_PlotGlobe_OrbitFilename.layout.border = '1px dashed blue'
	WIDGET_PlotGlobe_OrbitFilename.layout.padding = '0px 10px 0px 10px'
	InputsPanel.children += (WIDGET_PlotGlobe_OrbitFilename,)
	global WIDGET_PlotGlobe_OrbitVariableToPlot
	WIDGET_PlotGlobe_OrbitVariableToPlot = widgets.Text(value='ALFA')
	WIDGET_PlotGlobe_OrbitVariableToPlot.description = 'OrbitVariableToPlot'
	WIDGET_PlotGlobe_OrbitVariableToPlot.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_PlotGlobe_OrbitVariableToPlot,)
	global WIDGET_PlotGlobe_OrbitColorbarTitle
	WIDGET_PlotGlobe_OrbitColorbarTitle = widgets.Text()
	WIDGET_PlotGlobe_OrbitColorbarTitle.description = 'OrbitColorbarTitle'
	WIDGET_PlotGlobe_OrbitColorbarTitle.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_PlotGlobe_OrbitColorbarTitle,)
	global WIDGET_PlotGlobe_OrbitColorscaleName
	WIDGET_PlotGlobe_OrbitColorscaleName = widgets.Dropdown( options=['', 'None', 'Blackbody', 'Bluered', 'Blues', 'Earth', 'Electric', 'Greens', 'Greys', 'Hot', 'Jet', 'Picnic', 'Portland', 'Rainbow', 'RdBu', 'Reds', 'Viridis', 'YlGnBu', 'YlOrRd'], value='Jet', description='Variable')
	WIDGET_PlotGlobe_OrbitColorscaleName.description = 'OrbitColorscaleName'
	WIDGET_PlotGlobe_OrbitColorscaleName.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_PlotGlobe_OrbitColorscaleName,)
	global WIDGET_PlotGlobe_PlotTitle
	WIDGET_PlotGlobe_PlotTitle = widgets.Text()
	WIDGET_PlotGlobe_PlotTitle.description = 'PlotTitle'
	WIDGET_PlotGlobe_PlotTitle.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_PlotGlobe_PlotTitle,)
	Interpolator_Panel.children += (InputsPanel,)
	# Create widget for the moudle black-box body
	PlotGlobe_Btn = widgets.Button (
		description='PlotGlobe',
		tooltip="Creates a 3D plot of an earth globe. Can plot a sphere surface and/or a satellite orbit. \nThe surface and the orbit are colored according to the selected variable values in the data files. \nThe data files can be of CSV o NetCDF4 format.\nThe same or different variables can be selected for plotting at the surface and the orbit.\nThe same or different colorscales can be selected for the surface and the orbit.\nThe same or different value range can be selected for the surface and the orbit by assigning the same or different colorbar titles.\nValid Colorscale Names: ‘Blackbody’, ‘Bluered’, ‘Blues’, ‘Earth’, ‘Electric’, ‘Greens’, ‘Greys’, ‘Hot’, ‘Jet’, ‘Picnic’, ‘Portland’, ‘Rainbow’, ‘RdBu’, ‘Reds’, ‘Viridis’, ‘YlGnBu’, ‘YlOrRd’\nARGUMENTS:\n  SurfaceFilename: The file which contains data for a surface. If empty string then no surface will be plotted. \n                   CSV format: Time,Lat,Lon,Alt,value. Contains the data for the sphere surface. \n                   NetCDF format: see Daedalus User Manual\n  SurfaceVariableToPlot: The name of the variable to read from the surface data file and plot upon the globe.\n  SurfaceColorbarTitle: A title to display above the colorbar which refers to the surface data\n  SurfaceColorscaleName: The name of the Colorscale to use for the surface data. \n                         In case an empty string is passed then a default HeatMap colorscale will be applied.\n                         In case None is passed then all points will be black irrespective of value.\n  OrbitFilename: counterpart of the corresponding argument for the Surface data\n  OrbitVariableToPlot: counterpart of the corresponding argument for the Surface data \n  OrbitColorbarTitle: counterpart of the corresponding argument for the Surface data \n  OrbitColorscaleName: counterpart of the corresponding argument for the Surface data\n  PlotTitle: A title which is displayed on the top of the plot.\nRETURNS: a string containing information about the Data\n",
	)
	PlotGlobe_Btn.layout.min_width = '200px'
	PlotGlobe_Btn.style.button_color = 'deeppink'
	Interpolator_Panel.children += (PlotGlobe_Btn,)
	# Create widgets for module's outputs
	OutputsPanel = widgets.VBox()
	OutputsPanel.layout.min_width = '300px'
	Interpolator_Panel.children += (OutputsPanel,)
	return Interpolator_Panel


def Construct_MSISE00():
	#Create Containers
	MSISE00_Panel = widgets.Box()
	MSISE00_Panel.layout.overflow_x = 'scroll'
	MSISE00_Panel.layout.flex_flow = 'row'
	MSISE00_Panel.layout.display = 'flex'
	########
	## GUI code for module 'OrbitSelector'
	########
	# Create widgets for module's inputs
	InputsPanel = widgets.VBox()
	InputsPanel.layout.min_width = '330px'
	global WIDGET_OrbitSelector_Filename
	WIDGET_OrbitSelector_Filename = widgets.Text()
	WIDGET_OrbitSelector_Filename.description = 'Filename'
	WIDGET_OrbitSelector_Filename.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_OrbitSelector_Filename,)
	global WIDGET_OrbitSelector_EvtXY
	WIDGET_OrbitSelector_EvtXY = widgets.Dropdown( options=['', 'Evt0', 'Evt1', 'Evt2', 'Evt3', 'Evt4', 'Evt5', 'Evt6', 'Evt7', 'Evt8', 'Evt9'], description='Variable')
	WIDGET_OrbitSelector_EvtXY.description = 'EvtXY'
	WIDGET_OrbitSelector_EvtXY.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_OrbitSelector_EvtXY,)
	global WIDGET_OrbitSelector_TYP
	WIDGET_OrbitSelector_TYP = widgets.Dropdown( options=['LLA', 'VEL', 'PTG'], description='Variable')
	WIDGET_OrbitSelector_TYP.description = 'TYP'
	WIDGET_OrbitSelector_TYP.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_OrbitSelector_TYP,)
	global WIDGET_OrbitSelector_PerYYY
	WIDGET_OrbitSelector_PerYYY = widgets.Dropdown( options=['Per120', 'Per150'], description='Variable')
	WIDGET_OrbitSelector_PerYYY.description = 'PerYYY'
	WIDGET_OrbitSelector_PerYYY.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_OrbitSelector_PerYYY,)
	global WIDGET_OrbitSelector_LatZZ
	WIDGET_OrbitSelector_LatZZ = widgets.Dropdown( options=['Lat00', 'Lat40', 'Lat80'],  description='Variable')
	WIDGET_OrbitSelector_LatZZ.description = 'LatZZ'
	WIDGET_OrbitSelector_LatZZ.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_OrbitSelector_LatZZ,)
	global WIDGET_OrbitSelector_SRXXHZ
	WIDGET_OrbitSelector_SRXXHZ = widgets.Dropdown( options=['Srt16Hz', 'Srt01Hz'],  description='Variable')
	WIDGET_OrbitSelector_SRXXHZ.description = 'SRXXHZ'
	WIDGET_OrbitSelector_SRXXHZ.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_OrbitSelector_SRXXHZ,)
	global WIDGET_OrbitSelector_SC
	WIDGET_OrbitSelector_SC = widgets.Dropdown( options=['Msc', 'Ssc'],  description='Variable')
	WIDGET_OrbitSelector_SC.description = 'SC'
	WIDGET_OrbitSelector_SC.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_OrbitSelector_SC,)
	MSISE00_Panel.children += (InputsPanel,)
	# Create widget for the moudle black-box body
	OrbitSelector_Btn = widgets.Button (
		description='OrbitSelector',
		tooltip="Constructs and returns an orbit full-path filename from the orbit properties.\nThe filename includes all parameters required to select the appropriate orbit.\nAll filenames  have the same length and same number of parameters. \nIn case OrbitFilename and EvtXY are empty then an empty string is returned.\nFILE FORMAT: DAED_ORB_EVTXY_TYP_PERYYY_LATZZ_SRTQQHz_XSC.csv\n  Parameter OrbitFilename:\n    If it is not empty then the rest arguments are ignored. \n    If it contains slashes then is is assumed it is a full path name and it is returned as it is.\n    If it does not contain slashes then the orbits path is added.\n    \n  Parameter EvtXY values:\n    EVTXS    X Event, Single Orbit\n    EVTXA    X Event, All Orbit\n    EVT1Y    1st Event: St Patrick’s day event [17 Mar 2015 – 20 Mar 2015]\n    EVT2Y    2nd Event ...  ...\n \n  Parameter TYP values:\n    LLA    Time,Latitude Longitude Altitude \n    VEL    Time,VmagVxVyVz\n    PTG    Time, X_GSE, Y_GSE, Z_GSE, RamX_GSE, RamY_GSE, RamZ_GSE\n  Parameter PerYYY values:\n    PER120, PER150    Perigee Altitude at 120km or 150km\n  Parameter LatZZ values:\n    LAT00, LAT40, LAT80    Perigee Latitude at 0°, 40° or 80°\n \n  Parameter SRXXHZ values:\n    SRT16Hz, SRT01Hz    Sampling rate 16Hz or 1HZ\n \n  Parameter SC values:\n    MSC    Mother Spacecraft\n    SSC    Sub-Spacecraft\n",
	)
	OrbitSelector_Btn.layout.min_width = '200px'
	OrbitSelector_Btn.style.button_color = 'YellowGreen'
	MSISE00_Panel.children += (OrbitSelector_Btn,)
	# Create widgets for module's outputs
	OutputsPanel = widgets.VBox()
	OutputsPanel.layout.min_width = '300px'
	WIDGET_OrbitSelector_OrbitCSVfilename = widgets.Label(value='  --> OrbitCSVfilename  ')
	WIDGET_OrbitSelector_OrbitCSVfilename.layout.border = '1px dashed green'
	WIDGET_OrbitSelector_OrbitCSVfilename.layout.margin ='0px 40px 0px 0px' 
	WIDGET_OrbitSelector_OrbitCSVfilename.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_OrbitSelector_OrbitCSVfilename,)
	WIDGET_OrbitSelector_OrbitNetCDFfilename = widgets.Label(value='  --> OrbitNetCDFfilename  ')
	WIDGET_OrbitSelector_OrbitNetCDFfilename.layout.border = '1px dashed green'
	WIDGET_OrbitSelector_OrbitNetCDFfilename.layout.margin ='0px 40px 0px 0px' 
	WIDGET_OrbitSelector_OrbitNetCDFfilename.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_OrbitSelector_OrbitNetCDFfilename,)
	MSISE00_Panel.children += (OutputsPanel,)
	########
	## GUI code for module 'MSISE00'
	########
	# Create widgets for module's inputs
	InputsPanel = widgets.VBox()
	InputsPanel.layout.min_width = '330px'
	WIDGET_MSISE00_CSVfilename_forOrbit = widgets.Label(value='  OrbitSelector.OrbitCSVfilename --> CSVfilename_forOrbit  ')
	WIDGET_MSISE00_CSVfilename_forOrbit.layout.border = '1px dashed blue'
	WIDGET_MSISE00_CSVfilename_forOrbit.layout.padding = '0px 10px 0px 10px'
	InputsPanel.children += (WIDGET_MSISE00_CSVfilename_forOrbit,)
	global WIDGET_MSISE00_Variable_forOrbit
	WIDGET_MSISE00_Variable_forOrbit = widgets.Dropdown( options=['', 'Tn', 'O','O2','N2', 'rho'], value='', description='Variable')
	WIDGET_MSISE00_Variable_forOrbit.description = 'Variable_forOrbit'
	WIDGET_MSISE00_Variable_forOrbit.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_MSISE00_Variable_forOrbit,)
	MSISE00_Panel.children += (InputsPanel,)
	# Create widget for the moudle black-box body
	MSISE00_Btn = widgets.Button (
		description='MSISE00',
		tooltip="MSISE-90 (Mass Spectrometer - Incoherent Scatter) empirical model (MSISE)\nfilename: the CSV file to read. Format: Time,Alt,Lon,Lat. The file defines the coordinates for which the MSISE calculation will take place.\nParameter: selects which values to calculate. Possible values are:\n    ''   for all\n    'Tn'  for Neutral Temperature\n    'O' for Oxygen atoms\n    'O2' for Oxygen molecules\n    'N2' for Nitrogen molecules\n    'rho' for Mass Densities\nReturns: the name of the CSV file which contains the calculated values. Format: Time,Alt,Lon,Lat,Tn,O,O2,N2,rho",
	)
	MSISE00_Btn.layout.min_width = '200px'
	MSISE00_Btn.style.button_color = 'gold'
	MSISE00_Panel.children += (MSISE00_Btn,)
	# Create widgets for module's outputs
	OutputsPanel = widgets.VBox()
	OutputsPanel.layout.min_width = '300px'
	WIDGET_MSISE00_ObitResultCSV = widgets.Label(value='  --> ObitResultCSV  ')
	WIDGET_MSISE00_ObitResultCSV.layout.border = '1px dashed green'
	WIDGET_MSISE00_ObitResultCSV.layout.margin ='0px 40px 0px 0px' 
	WIDGET_MSISE00_ObitResultCSV.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_MSISE00_ObitResultCSV,)
	MSISE00_Panel.children += (OutputsPanel,)
	########
	## GUI code for module 'CreateCSV_Sphere'
	########
	# Create widgets for module's inputs
	InputsPanel = widgets.VBox()
	InputsPanel.layout.min_width = '330px'
	global WIDGET_CreateCSV_Sphere_CSVfilename
	WIDGET_CreateCSV_Sphere_CSVfilename = widgets.Text(value="")
	WIDGET_CreateCSV_Sphere_CSVfilename.description = 'CSVfilename'
	WIDGET_CreateCSV_Sphere_CSVfilename.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_CreateCSV_Sphere_CSVfilename,)
	global WIDGET_CreateCSV_Sphere_fixedDatetimeString
	WIDGET_CreateCSV_Sphere_fixedDatetimeString = widgets.Text(value="Mar 17 2015 00:50:00.000000000")
	WIDGET_CreateCSV_Sphere_fixedDatetimeString.description = 'fixedDatetimeString'
	WIDGET_CreateCSV_Sphere_fixedDatetimeString.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_CreateCSV_Sphere_fixedDatetimeString,)
	global WIDGET_CreateCSV_Sphere_fixedAltitude
	WIDGET_CreateCSV_Sphere_fixedAltitude = widgets.FloatText(value=120)
	WIDGET_CreateCSV_Sphere_fixedAltitude.description = 'fixedAltitude'
	WIDGET_CreateCSV_Sphere_fixedAltitude.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_CreateCSV_Sphere_fixedAltitude,)
	global WIDGET_CreateCSV_Sphere_LatitudeStep
	WIDGET_CreateCSV_Sphere_LatitudeStep = widgets.FloatText(value=5)
	WIDGET_CreateCSV_Sphere_LatitudeStep.description = 'LatitudeStep'
	WIDGET_CreateCSV_Sphere_LatitudeStep.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_CreateCSV_Sphere_LatitudeStep,)
	global WIDGET_CreateCSV_Sphere_LongitudeStep
	WIDGET_CreateCSV_Sphere_LongitudeStep = widgets.FloatText(value=5)
	WIDGET_CreateCSV_Sphere_LongitudeStep.description = 'LongitudeStep'
	WIDGET_CreateCSV_Sphere_LongitudeStep.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_CreateCSV_Sphere_LongitudeStep,)
	MSISE00_Panel.children += (InputsPanel,)
	# Create widget for the moudle black-box body
	CreateCSV_Sphere_Btn = widgets.Button (
		description='CreateCSV_Sphere',
		tooltip="Creates a CSV file with the format [Time, Lat, Lon, Alt]\nTime and Alt are fixed values given as arguments and Lat and Lon are produced for the whole globe accodring to steps specified.\nCSVfilename: string, the csv filename to be written (overwrite mode)\nfixedDatetimeString: string, the date-time fixed text with format as example: Mar 17 2015 00:12:00.000. Fixed value for every csv entry.\nfixedAltitude: float positive, kilometers above earth. Fixed value for every csv entry.\nLatitudeStep: float positive, Latiude values range between 87.5 and -88.5. The step is how Lat should be incremented at each iteration.\nLongitudeStep: float positive, Longitude values range between -180 and 180.  The step is how Lon should be incremented at each iteration.\nRETURNS: the filename it has written (same as the argument)\n         the altitude (same as the argument)\n         the number of lines written",
	)
	CreateCSV_Sphere_Btn.layout.min_width = '200px'
	CreateCSV_Sphere_Btn.style.button_color = 'gold'
	MSISE00_Panel.children += (CreateCSV_Sphere_Btn,)
	# Create widgets for module's outputs
	OutputsPanel = widgets.VBox()
	OutputsPanel.layout.min_width = '300px'
	WIDGET_CreateCSV_Sphere_theCSVfilename = widgets.Label(value='  --> theCSVfilename  ')
	WIDGET_CreateCSV_Sphere_theCSVfilename.layout.border = '1px dashed green'
	WIDGET_CreateCSV_Sphere_theCSVfilename.layout.margin ='0px 40px 0px 0px' 
	WIDGET_CreateCSV_Sphere_theCSVfilename.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_CreateCSV_Sphere_theCSVfilename,)
	WIDGET_CreateCSV_Sphere_theAltitude = widgets.Label(value='  --> theAltitude  ')
	WIDGET_CreateCSV_Sphere_theAltitude.layout.border = '1px dashed green'
	WIDGET_CreateCSV_Sphere_theAltitude.layout.margin ='0px 40px 0px 0px' 
	WIDGET_CreateCSV_Sphere_theAltitude.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_CreateCSV_Sphere_theAltitude,)
	WIDGET_CreateCSV_Sphere_NumberOfLinesWritten = widgets.Label(value='  --> NumberOfLinesWritten  ')
	WIDGET_CreateCSV_Sphere_NumberOfLinesWritten.layout.border = '1px dashed green'
	WIDGET_CreateCSV_Sphere_NumberOfLinesWritten.layout.margin ='0px 40px 0px 0px' 
	WIDGET_CreateCSV_Sphere_NumberOfLinesWritten.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_CreateCSV_Sphere_NumberOfLinesWritten,)
	MSISE00_Panel.children += (OutputsPanel,)
	########
	## GUI code for module 'MSISE00'
	########
	# Create widgets for module's inputs
	InputsPanel = widgets.VBox()
	InputsPanel.layout.min_width = '330px'
	WIDGET_MSISE00_CSVfilename_forMeshgrid = widgets.Label(value='  CreateCSV_Sphere.theCSVfilename --> CSVfilename_forMeshgrid  ')
	WIDGET_MSISE00_CSVfilename_forMeshgrid.layout.border = '1px dashed blue'
	WIDGET_MSISE00_CSVfilename_forMeshgrid.layout.padding = '0px 10px 0px 10px'
	InputsPanel.children += (WIDGET_MSISE00_CSVfilename_forMeshgrid,)
	global WIDGET_MSISE00_Variable_forMeshgrid
	WIDGET_MSISE00_Variable_forMeshgrid = widgets.Dropdown( options=['', 'Tn', 'O','O2','N2', 'rho'], value='', description='Variable')
	WIDGET_MSISE00_Variable_forMeshgrid.description = 'Variable_forMeshgrid'
	WIDGET_MSISE00_Variable_forMeshgrid.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_MSISE00_Variable_forMeshgrid,)
	MSISE00_Panel.children += (InputsPanel,)
	# Create widget for the moudle black-box body
	MSISE00_Btn = widgets.Button (
		description='MSISE00',
		tooltip="MSISE-90 (Mass Spectrometer - Incoherent Scatter) empirical model (MSISE)\nfilename: the CSV file to read. Format: Time,Alt,Lon,Lat. The file defines the coordinates for which the MSISE calculation will take place.\nParameter: selects which values to calculate. Possible values are:\n    ''   for all\n    'Tn'  for Neutral Temperature\n    'O' for Oxygen atoms\n    'O2' for Oxygen molecules\n    'N2' for Nitrogen molecules\n    'rho' for Mass Densities\nReturns: the name of the CSV file which contains the calculated values. Format: Time,Alt,Lon,Lat,Tn,O,O2,N2,rho",
	)
	MSISE00_Btn.layout.min_width = '200px'
	MSISE00_Btn.style.button_color = 'gold'
	MSISE00_Panel.children += (MSISE00_Btn,)
	# Create widgets for module's outputs
	OutputsPanel = widgets.VBox()
	OutputsPanel.layout.min_width = '300px'
	WIDGET_MSISE00_MeshgridResultCSV = widgets.Label(value='  --> MeshgridResultCSV  ')
	WIDGET_MSISE00_MeshgridResultCSV.layout.border = '1px dashed green'
	WIDGET_MSISE00_MeshgridResultCSV.layout.margin ='0px 40px 0px 0px' 
	WIDGET_MSISE00_MeshgridResultCSV.layout.padding ='0px 10px 0px 10px' 
	OutputsPanel.children += (WIDGET_MSISE00_MeshgridResultCSV,)
	MSISE00_Panel.children += (OutputsPanel,)
	########
	## GUI code for module 'PlotGlobe'
	########
	# Create widgets for module's inputs
	InputsPanel = widgets.VBox()
	InputsPanel.layout.min_width = '330px'
	WIDGET_PlotGlobe_SurfaceFilename = widgets.Label(value='  MSISE00.MeshgridResultCSV --> SurfaceFilename  ')
	WIDGET_PlotGlobe_SurfaceFilename.layout.border = '1px dashed blue'
	WIDGET_PlotGlobe_SurfaceFilename.layout.padding = '0px 10px 0px 10px'
	InputsPanel.children += (WIDGET_PlotGlobe_SurfaceFilename,)
	global WIDGET_PlotGlobe_SurfaceVariableToPlot
	WIDGET_PlotGlobe_SurfaceVariableToPlot = widgets.Text(value='ALFA')
	WIDGET_PlotGlobe_SurfaceVariableToPlot.description = 'SurfaceVariableToPlot'
	WIDGET_PlotGlobe_SurfaceVariableToPlot.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_PlotGlobe_SurfaceVariableToPlot,)
	global WIDGET_PlotGlobe_SurfaceColorbarTitle
	WIDGET_PlotGlobe_SurfaceColorbarTitle = widgets.Text()
	WIDGET_PlotGlobe_SurfaceColorbarTitle.description = 'SurfaceColorbarTitle'
	WIDGET_PlotGlobe_SurfaceColorbarTitle.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_PlotGlobe_SurfaceColorbarTitle,)
	global WIDGET_PlotGlobe_SurfaceColorscaleName
	WIDGET_PlotGlobe_SurfaceColorscaleName = widgets.Dropdown( options=['', 'None', 'Blackbody', 'Bluered', 'Blues', 'Earth', 'Electric', 'Greens', 'Greys', 'Hot', 'Jet', 'Picnic', 'Portland', 'Rainbow', 'RdBu', 'Reds', 'Viridis', 'YlGnBu', 'YlOrRd'], value='Jet', description='Variable')
	WIDGET_PlotGlobe_SurfaceColorscaleName.description = 'SurfaceColorscaleName'
	WIDGET_PlotGlobe_SurfaceColorscaleName.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_PlotGlobe_SurfaceColorscaleName,)
	WIDGET_PlotGlobe_OrbitFilename = widgets.Label(value='  MSISE00.ObitResultCSV --> OrbitFilename  ')
	WIDGET_PlotGlobe_OrbitFilename.layout.border = '1px dashed blue'
	WIDGET_PlotGlobe_OrbitFilename.layout.padding = '0px 10px 0px 10px'
	InputsPanel.children += (WIDGET_PlotGlobe_OrbitFilename,)
	global WIDGET_PlotGlobe_OrbitVariableToPlot
	WIDGET_PlotGlobe_OrbitVariableToPlot = widgets.Text(value='ALFA')
	WIDGET_PlotGlobe_OrbitVariableToPlot.description = 'OrbitVariableToPlot'
	WIDGET_PlotGlobe_OrbitVariableToPlot.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_PlotGlobe_OrbitVariableToPlot,)
	global WIDGET_PlotGlobe_OrbitColorbarTitle
	WIDGET_PlotGlobe_OrbitColorbarTitle = widgets.Text()
	WIDGET_PlotGlobe_OrbitColorbarTitle.description = 'OrbitColorbarTitle'
	WIDGET_PlotGlobe_OrbitColorbarTitle.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_PlotGlobe_OrbitColorbarTitle,)
	global WIDGET_PlotGlobe_OrbitColorscaleName
	WIDGET_PlotGlobe_OrbitColorscaleName = widgets.Dropdown( options=['', 'None', 'Blackbody', 'Bluered', 'Blues', 'Earth', 'Electric', 'Greens', 'Greys', 'Hot', 'Jet', 'Picnic', 'Portland', 'Rainbow', 'RdBu', 'Reds', 'Viridis', 'YlGnBu', 'YlOrRd'], value='Jet', description='Variable')
	WIDGET_PlotGlobe_OrbitColorscaleName.description = 'OrbitColorscaleName'
	WIDGET_PlotGlobe_OrbitColorscaleName.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_PlotGlobe_OrbitColorscaleName,)
	global WIDGET_PlotGlobe_PlotTitle
	WIDGET_PlotGlobe_PlotTitle = widgets.Text()
	WIDGET_PlotGlobe_PlotTitle.description = 'PlotTitle'
	WIDGET_PlotGlobe_PlotTitle.layout.border = '1px dashed blue'
	InputsPanel.children += (WIDGET_PlotGlobe_PlotTitle,)
	MSISE00_Panel.children += (InputsPanel,)
	# Create widget for the moudle black-box body
	PlotGlobe_Btn = widgets.Button (
		description='PlotGlobe',
		tooltip="Creates a 3D plot of an earth globe. Can plot a sphere surface and/or a satellite orbit. \nThe surface and the orbit are colored according to the selected variable values in the data files. \nThe data files can be of CSV o NetCDF4 format.\nThe same or different variables can be selected for plotting at the surface and the orbit.\nThe same or different colorscales can be selected for the surface and the orbit.\nThe same or different value range can be selected for the surface and the orbit by assigning the same or different colorbar titles.\nValid Colorscale Names: ‘Blackbody’, ‘Bluered’, ‘Blues’, ‘Earth’, ‘Electric’, ‘Greens’, ‘Greys’, ‘Hot’, ‘Jet’, ‘Picnic’, ‘Portland’, ‘Rainbow’, ‘RdBu’, ‘Reds’, ‘Viridis’, ‘YlGnBu’, ‘YlOrRd’\nARGUMENTS:\n  SurfaceFilename: The file which contains data for a surface. If empty string then no surface will be plotted. \n                   CSV format: Time,Lat,Lon,Alt,value. Contains the data for the sphere surface. \n                   NetCDF format: see Daedalus User Manual\n  SurfaceVariableToPlot: The name of the variable to read from the surface data file and plot upon the globe.\n  SurfaceColorbarTitle: A title to display above the colorbar which refers to the surface data\n  SurfaceColorscaleName: The name of the Colorscale to use for the surface data. \n                         In case an empty string is passed then a default HeatMap colorscale will be applied.\n                         In case None is passed then all points will be black irrespective of value.\n  OrbitFilename: counterpart of the corresponding argument for the Surface data\n  OrbitVariableToPlot: counterpart of the corresponding argument for the Surface data \n  OrbitColorbarTitle: counterpart of the corresponding argument for the Surface data \n  OrbitColorscaleName: counterpart of the corresponding argument for the Surface data\n  PlotTitle: A title which is displayed on the top of the plot.\nRETURNS: a string containing information about the Data\n",
	)
	PlotGlobe_Btn.layout.min_width = '200px'
	PlotGlobe_Btn.style.button_color = 'deeppink'
	MSISE00_Panel.children += (PlotGlobe_Btn,)
	# Create widgets for module's outputs
	OutputsPanel = widgets.VBox()
	OutputsPanel.layout.min_width = '300px'
	MSISE00_Panel.children += (OutputsPanel,)
	return MSISE00_Panel


