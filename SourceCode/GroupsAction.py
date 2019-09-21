from SourceCode.ModulesSourceCode.OrbitSelector.OrbitSelector import OrbitSelector
from SourceCode.ModulesSourceCode.Interpolator.Interpolator import Interpolator
from SourceCode.ModulesSourceCode.PlotGlobe.PlotGlobe import PlotGlobe
from SourceCode.ModulesSourceCode.OrbitSelector.OrbitSelector import OrbitSelector
from SourceCode.ModulesSourceCode.HWM14.HWM14 import HWM14
from SourceCode.ModulesSourceCode.CreateCSV_Sphere.CreateCSV_Sphere import CreateCSV_Sphere
from SourceCode.ModulesSourceCode.HWM14.HWM14 import HWM14
from SourceCode.ModulesSourceCode.PlotGlobe.PlotGlobe import PlotGlobe
from SourceCode.ModulesSourceCode.OrbitSelector.OrbitSelector import OrbitSelector
from SourceCode.ModulesSourceCode.IGRF.IGRF import IGRF
from SourceCode.ModulesSourceCode.CreateCSV_Sphere.CreateCSV_Sphere import CreateCSV_Sphere
from SourceCode.ModulesSourceCode.IGRF.IGRF import IGRF
from SourceCode.ModulesSourceCode.PlotGlobe.PlotGlobe import PlotGlobe
from SourceCode.ModulesSourceCode.OrbitSelector.OrbitSelector import OrbitSelector
from SourceCode.ModulesSourceCode.Selector.Selector import Selector
from SourceCode.ModulesSourceCode.CreateCSV_Sphere.CreateCSV_Sphere import CreateCSV_Sphere
from SourceCode.ModulesSourceCode.IRI16.IRI16 import IRI16
from SourceCode.ModulesSourceCode.IRI16.IRI16 import IRI16
from SourceCode.ModulesSourceCode.SubGridVariability.SubGridVariability import SubGridVariability
from SourceCode.ModulesSourceCode.PlotGlobe.PlotGlobe import PlotGlobe
from SourceCode.ModulesSourceCode.OrbitSelector.OrbitSelector import OrbitSelector
from SourceCode.ModulesSourceCode.MSISE00.MSISE00 import MSISE00
from SourceCode.ModulesSourceCode.CreateCSV_Sphere.CreateCSV_Sphere import CreateCSV_Sphere
from SourceCode.ModulesSourceCode.MSISE00.MSISE00 import MSISE00
from SourceCode.ModulesSourceCode.PlotGlobe.PlotGlobe import PlotGlobe
import SourceCode.PanelDisplayer as PanelDisplayer

def Display( *arg ): # in order to avoid circular dependent imports
	import SourceCode.DaedalusMazeGUI as DaedalusMazeGUI
	DaedalusMazeGUI.Display( *arg )

"""
Executes these Simulation Modules:
  OrbitSelector  Interpolator  PlotGlobe  
"""
def Execute_Interpolator():
	Display('Executing Interpolator:')
	OrbitSelector_OrbitFilename = OrbitSelector( PanelDisplayer.WIDGET_OrbitSelector_Filename.value, PanelDisplayer.WIDGET_OrbitSelector_EvtXY.value, PanelDisplayer.WIDGET_OrbitSelector_TYP.value, PanelDisplayer.WIDGET_OrbitSelector_PerYYY.value, PanelDisplayer.WIDGET_OrbitSelector_LatZZ.value, PanelDisplayer.WIDGET_OrbitSelector_SRXXHZ.value, PanelDisplayer.WIDGET_OrbitSelector_SC.value )
	Display( '    OrbitSelector( ' , PanelDisplayer.WIDGET_OrbitSelector_Filename.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_EvtXY.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_TYP.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_PerYYY.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_LatZZ.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_SRXXHZ.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_SC.value , ' ) results:')
	Display( '          OrbitFilename: ' + str(OrbitSelector_OrbitFilename) )
	Interpolator_InterpolatedOrbitCSV = Interpolator( PanelDisplayer.WIDGET_Interpolator_model.value, PanelDisplayer.WIDGET_Interpolator_model_data_file.value, OrbitSelector_OrbitFilename, PanelDisplayer.WIDGET_Interpolator_save.value, PanelDisplayer.WIDGET_Interpolator_VAR.value )
	Display( '    Interpolator( ' , PanelDisplayer.WIDGET_Interpolator_model.value,', ' , PanelDisplayer.WIDGET_Interpolator_model_data_file.value,', ' , OrbitSelector_OrbitFilename,', ' , PanelDisplayer.WIDGET_Interpolator_save.value,', ' , PanelDisplayer.WIDGET_Interpolator_VAR.value , ' ) results:')
	Display( '          InterpolatedOrbitCSV: ' + str(Interpolator_InterpolatedOrbitCSV) )
	PlotGlobe( PanelDisplayer.WIDGET_PlotGlobe_DataCSVfilename.value, Interpolator_InterpolatedOrbitCSV, PanelDisplayer.WIDGET_PlotGlobe_GraphTitle.value, PanelDisplayer.WIDGET_PlotGlobe_ColorbarTitle.value, PanelDisplayer.WIDGET_PlotGlobe_ColorscaleName.value )
	Display( '    PlotGlobe( ' , PanelDisplayer.WIDGET_PlotGlobe_DataCSVfilename.value,', ' , Interpolator_InterpolatedOrbitCSV,', ' , PanelDisplayer.WIDGET_PlotGlobe_GraphTitle.value,', ' , PanelDisplayer.WIDGET_PlotGlobe_ColorbarTitle.value,', ' , PanelDisplayer.WIDGET_PlotGlobe_ColorscaleName.value , ' ) results:')
	Display( '          -void-' )
	Display('')
"""
Executes these Simulation Modules:
  OrbitSelector  HWM14  CreateCSV_Sphere  HWM14  PlotGlobe  
"""
def Execute_HWM14():
	Display('Executing HWM14:')
	OrbitSelector_OrbitFilename = OrbitSelector( PanelDisplayer.WIDGET_OrbitSelector_Filename.value, PanelDisplayer.WIDGET_OrbitSelector_EvtXY.value, PanelDisplayer.WIDGET_OrbitSelector_TYP.value, PanelDisplayer.WIDGET_OrbitSelector_PerYYY.value, PanelDisplayer.WIDGET_OrbitSelector_LatZZ.value, PanelDisplayer.WIDGET_OrbitSelector_SRXXHZ.value, PanelDisplayer.WIDGET_OrbitSelector_SC.value )
	Display( '    OrbitSelector( ' , PanelDisplayer.WIDGET_OrbitSelector_Filename.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_EvtXY.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_TYP.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_PerYYY.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_LatZZ.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_SRXXHZ.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_SC.value , ' ) results:')
	Display( '          OrbitFilename: ' + str(OrbitSelector_OrbitFilename) )
	HWM14_OrbitResultCSV = HWM14( OrbitSelector_OrbitFilename, PanelDisplayer.WIDGET_HWM14_Variable_forOrbit.value )
	Display( '    HWM14( ' , OrbitSelector_OrbitFilename,', ' , PanelDisplayer.WIDGET_HWM14_Variable_forOrbit.value , ' ) results:')
	Display( '          OrbitResultCSV: ' + str(HWM14_OrbitResultCSV) )
	CreateCSV_Sphere_theCSVfilename, CreateCSV_Sphere_theAltitude, CreateCSV_Sphere_NumberOfLinesWritten = CreateCSV_Sphere( PanelDisplayer.WIDGET_CreateCSV_Sphere_CSVfilename.value, PanelDisplayer.WIDGET_CreateCSV_Sphere_fixedDatetimeString.value, PanelDisplayer.WIDGET_CreateCSV_Sphere_fixedAltitude.value, PanelDisplayer.WIDGET_CreateCSV_Sphere_LatitudeStep.value, PanelDisplayer.WIDGET_CreateCSV_Sphere_LongitudeStep.value )
	Display( '    CreateCSV_Sphere( ' , PanelDisplayer.WIDGET_CreateCSV_Sphere_CSVfilename.value,', ' , PanelDisplayer.WIDGET_CreateCSV_Sphere_fixedDatetimeString.value,', ' , PanelDisplayer.WIDGET_CreateCSV_Sphere_fixedAltitude.value,', ' , PanelDisplayer.WIDGET_CreateCSV_Sphere_LatitudeStep.value,', ' , PanelDisplayer.WIDGET_CreateCSV_Sphere_LongitudeStep.value , ' ) results:')
	Display( '          theCSVfilename: ' + str(CreateCSV_Sphere_theCSVfilename) )
	Display( '          theAltitude: ' + str(CreateCSV_Sphere_theAltitude) )
	Display( '          NumberOfLinesWritten: ' + str(CreateCSV_Sphere_NumberOfLinesWritten) )
	HWM14_MeshgridResultCSV = HWM14( CreateCSV_Sphere_theCSVfilename, PanelDisplayer.WIDGET_HWM14_Variable_forMeshgrid.value )
	Display( '    HWM14( ' , CreateCSV_Sphere_theCSVfilename,', ' , PanelDisplayer.WIDGET_HWM14_Variable_forMeshgrid.value , ' ) results:')
	Display( '          MeshgridResultCSV: ' + str(HWM14_MeshgridResultCSV) )
	PlotGlobe( HWM14_MeshgridResultCSV, HWM14_OrbitResultCSV, PanelDisplayer.WIDGET_PlotGlobe_PlotTitle.value, PanelDisplayer.WIDGET_PlotGlobe_ColorbarTitle.value, PanelDisplayer.WIDGET_PlotGlobe_ColorscaleName.value )
	Display( '    PlotGlobe( ' , HWM14_MeshgridResultCSV,', ' , HWM14_OrbitResultCSV,', ' , PanelDisplayer.WIDGET_PlotGlobe_PlotTitle.value,', ' , PanelDisplayer.WIDGET_PlotGlobe_ColorbarTitle.value,', ' , PanelDisplayer.WIDGET_PlotGlobe_ColorscaleName.value , ' ) results:')
	Display( '          -void-' )
	Display('')
"""
Executes these Simulation Modules:
  OrbitSelector  MSISE00  CreateCSV_Sphere  MSISE00  PlotGlobe  
"""
def Execute_MSISE00():
	Display('Executing MSISE00:')
	OrbitSelector_OrbitFilename = OrbitSelector( PanelDisplayer.WIDGET_OrbitSelector_Filename.value, PanelDisplayer.WIDGET_OrbitSelector_EvtXY.value, PanelDisplayer.WIDGET_OrbitSelector_TYP.value, PanelDisplayer.WIDGET_OrbitSelector_PerYYY.value, PanelDisplayer.WIDGET_OrbitSelector_LatZZ.value, PanelDisplayer.WIDGET_OrbitSelector_SRXXHZ.value, PanelDisplayer.WIDGET_OrbitSelector_SC.value )
	Display( '    OrbitSelector( ' , PanelDisplayer.WIDGET_OrbitSelector_Filename.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_EvtXY.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_TYP.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_PerYYY.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_LatZZ.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_SRXXHZ.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_SC.value , ' ) results:')
	Display( '          OrbitFilename: ' + str(OrbitSelector_OrbitFilename) )
	MSISE00_ObitResultCSV = MSISE00( OrbitSelector_OrbitFilename, PanelDisplayer.WIDGET_MSISE00_Variable_forOrbit.value )
	Display( '    MSISE00( ' , OrbitSelector_OrbitFilename,', ' , PanelDisplayer.WIDGET_MSISE00_Variable_forOrbit.value , ' ) results:')
	Display( '          ObitResultCSV: ' + str(MSISE00_ObitResultCSV) )
	CreateCSV_Sphere_theCSVfilename, CreateCSV_Sphere_theAltitude, CreateCSV_Sphere_NumberOfLinesWritten = CreateCSV_Sphere( PanelDisplayer.WIDGET_CreateCSV_Sphere_CSVfilename.value, PanelDisplayer.WIDGET_CreateCSV_Sphere_fixedDatetimeString.value, PanelDisplayer.WIDGET_CreateCSV_Sphere_fixedAltitude.value, PanelDisplayer.WIDGET_CreateCSV_Sphere_LatitudeStep.value, PanelDisplayer.WIDGET_CreateCSV_Sphere_LongitudeStep.value )
	Display( '    CreateCSV_Sphere( ' , PanelDisplayer.WIDGET_CreateCSV_Sphere_CSVfilename.value,', ' , PanelDisplayer.WIDGET_CreateCSV_Sphere_fixedDatetimeString.value,', ' , PanelDisplayer.WIDGET_CreateCSV_Sphere_fixedAltitude.value,', ' , PanelDisplayer.WIDGET_CreateCSV_Sphere_LatitudeStep.value,', ' , PanelDisplayer.WIDGET_CreateCSV_Sphere_LongitudeStep.value , ' ) results:')
	Display( '          theCSVfilename: ' + str(CreateCSV_Sphere_theCSVfilename) )
	Display( '          theAltitude: ' + str(CreateCSV_Sphere_theAltitude) )
	Display( '          NumberOfLinesWritten: ' + str(CreateCSV_Sphere_NumberOfLinesWritten) )
	MSISE00_MeshgridResultCSV = MSISE00( CreateCSV_Sphere_theCSVfilename, PanelDisplayer.WIDGET_MSISE00_Variable_forMeshgrid.value )
	Display( '    MSISE00( ' , CreateCSV_Sphere_theCSVfilename,', ' , PanelDisplayer.WIDGET_MSISE00_Variable_forMeshgrid.value , ' ) results:')
	Display( '          MeshgridResultCSV: ' + str(MSISE00_MeshgridResultCSV) )
	PlotGlobe( MSISE00_MeshgridResultCSV, MSISE00_ObitResultCSV, PanelDisplayer.WIDGET_PlotGlobe_PlotTitle.value, PanelDisplayer.WIDGET_PlotGlobe_ColorbarTitle.value, PanelDisplayer.WIDGET_PlotGlobe_ColorscaleName.value )
	Display( '    PlotGlobe( ' , MSISE00_MeshgridResultCSV,', ' , MSISE00_ObitResultCSV,', ' , PanelDisplayer.WIDGET_PlotGlobe_PlotTitle.value,', ' , PanelDisplayer.WIDGET_PlotGlobe_ColorbarTitle.value,', ' , PanelDisplayer.WIDGET_PlotGlobe_ColorscaleName.value , ' ) results:')
	Display( '          -void-' )
	Display('')
"""
Executes these Simulation Modules:
  OrbitSelector  Selector  CreateCSV_Sphere  IRI16  IRI16  SubGridVariability  PlotGlobe  
"""
def Execute_IRI16():
	Display('Executing IRI16:')
	OrbitSelector_OrbitFilename = OrbitSelector( PanelDisplayer.WIDGET_OrbitSelector_Filename.value, PanelDisplayer.WIDGET_OrbitSelector_EvtXY.value, PanelDisplayer.WIDGET_OrbitSelector_TYP.value, PanelDisplayer.WIDGET_OrbitSelector_PerYYY.value, PanelDisplayer.WIDGET_OrbitSelector_LatZZ.value, PanelDisplayer.WIDGET_OrbitSelector_SRXXHZ.value, PanelDisplayer.WIDGET_OrbitSelector_SC.value )
	Display( '    OrbitSelector( ' , PanelDisplayer.WIDGET_OrbitSelector_Filename.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_EvtXY.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_TYP.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_PerYYY.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_LatZZ.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_SRXXHZ.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_SC.value , ' ) results:')
	Display( '          OrbitFilename: ' + str(OrbitSelector_OrbitFilename) )
	Selector_MeshgridCSVfilenameOUT, Selector_MeshgridDatetimeOUT, Selector_MeshgridAltitudeOUT, Selector_MeshgridLatitudeStepOUT, Selector_MeshgridLongitudeStepOUT, Selector_NoiseWxOUT, Selector_NoiseWyOUT, Selector_NoiseAxOUT, Selector_NoiseAyOUT, Selector_IRI16_variableForMeshgrid_OUT, Selector_IRI16_variableForOrbit_OUT = Selector( PanelDisplayer.WIDGET_Selector_MeshgridCSVfilename.value, PanelDisplayer.WIDGET_Selector_MeshgridDatetime.value, PanelDisplayer.WIDGET_Selector_MeshgridAltitude.value, PanelDisplayer.WIDGET_Selector_MeshgridLatitudeStep.value, PanelDisplayer.WIDGET_Selector_MeshgridLongitudeStep.value, PanelDisplayer.WIDGET_Selector_NoiseWx.value, PanelDisplayer.WIDGET_Selector_NoiseWy.value, PanelDisplayer.WIDGET_Selector_NoiseAx.value, PanelDisplayer.WIDGET_Selector_NoiseAy.value, PanelDisplayer.WIDGET_Selector_IRI16_variableForMeshgrid.value, PanelDisplayer.WIDGET_Selector_IRI16_variableForOrbit.value )
	Display( '    Selector( ' , PanelDisplayer.WIDGET_Selector_MeshgridCSVfilename.value,', ' , PanelDisplayer.WIDGET_Selector_MeshgridDatetime.value,', ' , PanelDisplayer.WIDGET_Selector_MeshgridAltitude.value,', ' , PanelDisplayer.WIDGET_Selector_MeshgridLatitudeStep.value,', ' , PanelDisplayer.WIDGET_Selector_MeshgridLongitudeStep.value,', ' , PanelDisplayer.WIDGET_Selector_NoiseWx.value,', ' , PanelDisplayer.WIDGET_Selector_NoiseWy.value,', ' , PanelDisplayer.WIDGET_Selector_NoiseAx.value,', ' , PanelDisplayer.WIDGET_Selector_NoiseAy.value,', ' , PanelDisplayer.WIDGET_Selector_IRI16_variableForMeshgrid.value,', ' , PanelDisplayer.WIDGET_Selector_IRI16_variableForOrbit.value , ' ) results:')
	Display( '          MeshgridCSVfilenameOUT: ' + str(Selector_MeshgridCSVfilenameOUT) )
	Display( '          MeshgridDatetimeOUT: ' + str(Selector_MeshgridDatetimeOUT) )
	Display( '          MeshgridAltitudeOUT: ' + str(Selector_MeshgridAltitudeOUT) )
	Display( '          MeshgridLatitudeStepOUT: ' + str(Selector_MeshgridLatitudeStepOUT) )
	Display( '          MeshgridLongitudeStepOUT: ' + str(Selector_MeshgridLongitudeStepOUT) )
	Display( '          NoiseWxOUT: ' + str(Selector_NoiseWxOUT) )
	Display( '          NoiseWyOUT: ' + str(Selector_NoiseWyOUT) )
	Display( '          NoiseAxOUT: ' + str(Selector_NoiseAxOUT) )
	Display( '          NoiseAyOUT: ' + str(Selector_NoiseAyOUT) )
	Display( '          IRI16_variableForMeshgrid_OUT: ' + str(Selector_IRI16_variableForMeshgrid_OUT) )
	Display( '          IRI16_variableForOrbit_OUT: ' + str(Selector_IRI16_variableForOrbit_OUT) )
	CreateCSV_Sphere_theCSVfilename, CreateCSV_Sphere_theAltitude, CreateCSV_Sphere_NumberOfLinesWritten = CreateCSV_Sphere( Selector_MeshgridCSVfilenameOUT, Selector_MeshgridDatetimeOUT, Selector_MeshgridAltitudeOUT, Selector_MeshgridLatitudeStepOUT, Selector_MeshgridLongitudeStepOUT )
	Display( '    CreateCSV_Sphere( ' , Selector_MeshgridCSVfilenameOUT,', ' , Selector_MeshgridDatetimeOUT,', ' , Selector_MeshgridAltitudeOUT,', ' , Selector_MeshgridLatitudeStepOUT,', ' , Selector_MeshgridLongitudeStepOUT , ' ) results:')
	Display( '          theCSVfilename: ' + str(CreateCSV_Sphere_theCSVfilename) )
	Display( '          theAltitude: ' + str(CreateCSV_Sphere_theAltitude) )
	Display( '          NumberOfLinesWritten: ' + str(CreateCSV_Sphere_NumberOfLinesWritten) )
	IRI16_MeshgridResultCSV = IRI16( CreateCSV_Sphere_theCSVfilename, Selector_IRI16_variableForMeshgrid_OUT )
	Display( '    IRI16( ' , CreateCSV_Sphere_theCSVfilename,', ' , Selector_IRI16_variableForMeshgrid_OUT , ' ) results:')
	Display( '          MeshgridResultCSV: ' + str(IRI16_MeshgridResultCSV) )
	IRI16_orbitCSVfilename = IRI16( OrbitSelector_OrbitFilename, Selector_IRI16_variableForOrbit_OUT )
	Display( '    IRI16( ' , OrbitSelector_OrbitFilename,', ' , Selector_IRI16_variableForOrbit_OUT , ' ) results:')
	Display( '          orbitCSVfilename: ' + str(IRI16_orbitCSVfilename) )
	SubGridVariability_VariabilityCSVfilename = SubGridVariability( Selector_NoiseWxOUT, Selector_NoiseWyOUT, Selector_NoiseAxOUT, Selector_NoiseAyOUT, IRI16_orbitCSVfilename, PanelDisplayer.WIDGET_SubGridVariability_ValueName.value )
	Display( '    SubGridVariability( ' , Selector_NoiseWxOUT,', ' , Selector_NoiseWyOUT,', ' , Selector_NoiseAxOUT,', ' , Selector_NoiseAyOUT,', ' , IRI16_orbitCSVfilename,', ' , PanelDisplayer.WIDGET_SubGridVariability_ValueName.value , ' ) results:')
	Display( '          VariabilityCSVfilename: ' + str(SubGridVariability_VariabilityCSVfilename) )
	PlotGlobe( IRI16_MeshgridResultCSV, SubGridVariability_VariabilityCSVfilename, PanelDisplayer.WIDGET_PlotGlobe_PlotTitle.value, PanelDisplayer.WIDGET_PlotGlobe_ColorbarTitle.value, PanelDisplayer.WIDGET_PlotGlobe_ColorscaleName.value )
	Display( '    PlotGlobe( ' , IRI16_MeshgridResultCSV,', ' , SubGridVariability_VariabilityCSVfilename,', ' , PanelDisplayer.WIDGET_PlotGlobe_PlotTitle.value,', ' , PanelDisplayer.WIDGET_PlotGlobe_ColorbarTitle.value,', ' , PanelDisplayer.WIDGET_PlotGlobe_ColorscaleName.value , ' ) results:')
	Display( '          -void-' )
	Display('')
"""
Executes these Simulation Modules:
  OrbitSelector  IGRF  CreateCSV_Sphere  IGRF  PlotGlobe  
"""
def Execute_IGRF():
	Display('Executing IGRF:')
	OrbitSelector_OrbitFilename = OrbitSelector( PanelDisplayer.WIDGET_OrbitSelector_Filename.value, PanelDisplayer.WIDGET_OrbitSelector_EvtXY.value, PanelDisplayer.WIDGET_OrbitSelector_TYP.value, PanelDisplayer.WIDGET_OrbitSelector_PerYYY.value, PanelDisplayer.WIDGET_OrbitSelector_LatZZ.value, PanelDisplayer.WIDGET_OrbitSelector_SRXXHZ.value, PanelDisplayer.WIDGET_OrbitSelector_SC.value )
	Display( '    OrbitSelector( ' , PanelDisplayer.WIDGET_OrbitSelector_Filename.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_EvtXY.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_TYP.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_PerYYY.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_LatZZ.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_SRXXHZ.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_SC.value , ' ) results:')
	Display( '          OrbitFilename: ' + str(OrbitSelector_OrbitFilename) )
	IGRF_OrbitResultCSV = IGRF( OrbitSelector_OrbitFilename, PanelDisplayer.WIDGET_IGRF_Variable_forOrbit.value )
	Display( '    IGRF( ' , OrbitSelector_OrbitFilename,', ' , PanelDisplayer.WIDGET_IGRF_Variable_forOrbit.value , ' ) results:')
	Display( '          OrbitResultCSV: ' + str(IGRF_OrbitResultCSV) )
	CreateCSV_Sphere_theCSVfilename, CreateCSV_Sphere_theAltitude, CreateCSV_Sphere_NumberOfLinesWritten = CreateCSV_Sphere( PanelDisplayer.WIDGET_CreateCSV_Sphere_CSVfilename.value, PanelDisplayer.WIDGET_CreateCSV_Sphere_fixedDatetimeString.value, PanelDisplayer.WIDGET_CreateCSV_Sphere_fixedAltitude.value, PanelDisplayer.WIDGET_CreateCSV_Sphere_LatitudeStep.value, PanelDisplayer.WIDGET_CreateCSV_Sphere_LongitudeStep.value )
	Display( '    CreateCSV_Sphere( ' , PanelDisplayer.WIDGET_CreateCSV_Sphere_CSVfilename.value,', ' , PanelDisplayer.WIDGET_CreateCSV_Sphere_fixedDatetimeString.value,', ' , PanelDisplayer.WIDGET_CreateCSV_Sphere_fixedAltitude.value,', ' , PanelDisplayer.WIDGET_CreateCSV_Sphere_LatitudeStep.value,', ' , PanelDisplayer.WIDGET_CreateCSV_Sphere_LongitudeStep.value , ' ) results:')
	Display( '          theCSVfilename: ' + str(CreateCSV_Sphere_theCSVfilename) )
	Display( '          theAltitude: ' + str(CreateCSV_Sphere_theAltitude) )
	Display( '          NumberOfLinesWritten: ' + str(CreateCSV_Sphere_NumberOfLinesWritten) )
	IGRF_MeshgridResultCSV = IGRF( CreateCSV_Sphere_theCSVfilename, PanelDisplayer.WIDGET_IGRF_Variable_forMeshgrid.value )
	Display( '    IGRF( ' , CreateCSV_Sphere_theCSVfilename,', ' , PanelDisplayer.WIDGET_IGRF_Variable_forMeshgrid.value , ' ) results:')
	Display( '          MeshgridResultCSV: ' + str(IGRF_MeshgridResultCSV) )
	PlotGlobe( IGRF_MeshgridResultCSV, IGRF_OrbitResultCSV, PanelDisplayer.WIDGET_PlotGlobe_PlotTitle.value, PanelDisplayer.WIDGET_PlotGlobe_ColorbarTitle.value, PanelDisplayer.WIDGET_PlotGlobe_ColorscaleName.value )
	Display( '    PlotGlobe( ' , IGRF_MeshgridResultCSV,', ' , IGRF_OrbitResultCSV,', ' , PanelDisplayer.WIDGET_PlotGlobe_PlotTitle.value,', ' , PanelDisplayer.WIDGET_PlotGlobe_ColorbarTitle.value,', ' , PanelDisplayer.WIDGET_PlotGlobe_ColorscaleName.value , ' ) results:')
	Display( '          -void-' )
	Display('')
