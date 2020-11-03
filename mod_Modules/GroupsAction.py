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
from SourceCode.ModulesSourceCode.OrbitSelector.OrbitSelector import OrbitSelector
from SourceCode.ModulesSourceCode.Interpolator.Interpolator import Interpolator
from SourceCode.ModulesSourceCode.CreateNETCDFsphere.CreateNETCDFsphere import CreateNETCDFsphere
from SourceCode.ModulesSourceCode.Interpolator.Interpolator import Interpolator
from SourceCode.ModulesSourceCode.PlotGlobe.PlotGlobe import PlotGlobe
from SourceCode.ModulesSourceCode.OrbitSelector.OrbitSelector import OrbitSelector
from SourceCode.ModulesSourceCode.Selector.Selector import Selector
from SourceCode.ModulesSourceCode.Interpolator.Interpolator import Interpolator
from SourceCode.ModulesSourceCode.SGV_vi.SGV_vi import SGV_vi
from SourceCode.ModulesSourceCode.SGV_Ti.SGV_Ti import SGV_Ti
from SourceCode.ModulesSourceCode.SGV_Te.SGV_Te import SGV_Te
from SourceCode.ModulesSourceCode.SGV_Ni.SGV_Ni import SGV_Ni
from SourceCode.ModulesSourceCode.SGV_nix.SGV_nix import SGV_nix
from SourceCode.ModulesSourceCode.SGV_Ne.SGV_Ne import SGV_Ne
from SourceCode.ModulesSourceCode.SGV_TEC.SGV_TEC import SGV_TEC
from SourceCode.ModulesSourceCode.SGV_un.SGV_un import SGV_un
from SourceCode.ModulesSourceCode.SGV_Tn.SGV_Tn import SGV_Tn
from SourceCode.ModulesSourceCode.SGV_Nn.SGV_Nn import SGV_Nn
from SourceCode.ModulesSourceCode.SGV_nnx.SGV_nnx import SGV_nnx
from SourceCode.ModulesSourceCode.SGV_JEPP.SGV_JEPP import SGV_JEPP
from SourceCode.ModulesSourceCode.SGV_B.SGV_B import SGV_B
from SourceCode.ModulesSourceCode.SGV_E.SGV_E import SGV_E
from SourceCode.ModulesSourceCode.WakeEffectsCalc.WakeEffectsCalc import WakeEffectsCalc
from SourceCode.ModulesSourceCode.ChargingEffectsCalc.ChargingEffectsCalc import ChargingEffectsCalc
from SourceCode.ModulesSourceCode.InstrumentIDMRPA.InstrumentIDMRPA import InstrumentIDMRPA
from SourceCode.ModulesSourceCode.InstrumentTII.InstrumentTII import InstrumentTII
from SourceCode.ModulesSourceCode.InstrumentRWSCWS.InstrumentRWSCWS import InstrumentRWSCWS
from SourceCode.ModulesSourceCode.InstrumentIMS.InstrumentIMS import InstrumentIMS
from SourceCode.ModulesSourceCode.InstrumentNMS.InstrumentNMS import InstrumentNMS
from SourceCode.ModulesSourceCode.InstrumentACC.InstrumentACC import InstrumentACC
from SourceCode.ModulesSourceCode.InstrumentEPDS3.InstrumentEPDS3 import InstrumentEPDS3
from SourceCode.ModulesSourceCode.InstrumentEFI.InstrumentEFI import InstrumentEFI
from SourceCode.ModulesSourceCode.InstrumentMIP.InstrumentMIP import InstrumentMIP
from SourceCode.ModulesSourceCode.InstrumentMAG.InstrumentMAG import InstrumentMAG
from SourceCode.ModulesSourceCode.InstrumentGNSS.InstrumentGNSS import InstrumentGNSS
from SourceCode.ModulesSourceCode.DerivedProductsCalc.DerivedProductsCalc import DerivedProductsCalc
from SourceCode.ModulesSourceCode.Test1.Test1 import Test1
import SourceCode.PanelDisplayer as PanelDisplayer

def Display( *arg ): # in order to avoid circular dependent imports
	import SourceCode.DaedalusMazeGUI as DaedalusMazeGUI
	DaedalusMazeGUI.Display( *arg )

"""
Executes these Simulation Modules:
  OrbitSelector  IGRF  CreateCSV_Sphere  IGRF  PlotGlobe  
"""
def Execute_IGRF():
	Display('Executing IGRF:')
	OrbitSelector_OrbitCSVfilename, OrbitSelector_OrbitNetCDFfilename = OrbitSelector( PanelDisplayer.WIDGET_OrbitSelector_Filename.value, PanelDisplayer.WIDGET_OrbitSelector_EvtXY.value, PanelDisplayer.WIDGET_OrbitSelector_TYP.value, PanelDisplayer.WIDGET_OrbitSelector_PerYYY.value, PanelDisplayer.WIDGET_OrbitSelector_LatZZ.value, PanelDisplayer.WIDGET_OrbitSelector_SRXXHZ.value, PanelDisplayer.WIDGET_OrbitSelector_SC.value )
	Display( '    OrbitSelector( ' , PanelDisplayer.WIDGET_OrbitSelector_Filename.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_EvtXY.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_TYP.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_PerYYY.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_LatZZ.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_SRXXHZ.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_SC.value , ' ) results:')
	Display( '          OrbitCSVfilename: ' + str(OrbitSelector_OrbitCSVfilename) )
	Display( '          OrbitNetCDFfilename: ' + str(OrbitSelector_OrbitNetCDFfilename) )
	IGRF_OrbitResultCSV = IGRF( OrbitSelector_OrbitCSVfilename, PanelDisplayer.WIDGET_IGRF_Variable_forOrbit.value )
	Display( '    IGRF( ' , OrbitSelector_OrbitCSVfilename,', ' , PanelDisplayer.WIDGET_IGRF_Variable_forOrbit.value , ' ) results:')
	Display( '          OrbitResultCSV: ' + str(IGRF_OrbitResultCSV) )
	CreateCSV_Sphere_theCSVfilename, CreateCSV_Sphere_theAltitude, CreateCSV_Sphere_NumberOfLinesWritten = CreateCSV_Sphere( PanelDisplayer.WIDGET_CreateCSV_Sphere_CSVfilename.value, PanelDisplayer.WIDGET_CreateCSV_Sphere_fixedDatetimeString.value, PanelDisplayer.WIDGET_CreateCSV_Sphere_fixedAltitude.value, PanelDisplayer.WIDGET_CreateCSV_Sphere_LatitudeStep.value, PanelDisplayer.WIDGET_CreateCSV_Sphere_LongitudeStep.value )
	Display( '    CreateCSV_Sphere( ' , PanelDisplayer.WIDGET_CreateCSV_Sphere_CSVfilename.value,', ' , PanelDisplayer.WIDGET_CreateCSV_Sphere_fixedDatetimeString.value,', ' , PanelDisplayer.WIDGET_CreateCSV_Sphere_fixedAltitude.value,', ' , PanelDisplayer.WIDGET_CreateCSV_Sphere_LatitudeStep.value,', ' , PanelDisplayer.WIDGET_CreateCSV_Sphere_LongitudeStep.value , ' ) results:')
	Display( '          theCSVfilename: ' + str(CreateCSV_Sphere_theCSVfilename) )
	Display( '          theAltitude: ' + str(CreateCSV_Sphere_theAltitude) )
	Display( '          NumberOfLinesWritten: ' + str(CreateCSV_Sphere_NumberOfLinesWritten) )
	IGRF_MeshgridResultCSV = IGRF( CreateCSV_Sphere_theCSVfilename, PanelDisplayer.WIDGET_IGRF_Variable_forMeshgrid.value )
	Display( '    IGRF( ' , CreateCSV_Sphere_theCSVfilename,', ' , PanelDisplayer.WIDGET_IGRF_Variable_forMeshgrid.value , ' ) results:')
	Display( '          MeshgridResultCSV: ' + str(IGRF_MeshgridResultCSV) )
	PlotGlobe( IGRF_MeshgridResultCSV, PanelDisplayer.WIDGET_PlotGlobe_SurfaceVariableToPlot.value, PanelDisplayer.WIDGET_PlotGlobe_SurfaceColorbarTitle.value, PanelDisplayer.WIDGET_PlotGlobe_SurfaceColorscaleName.value, IGRF_OrbitResultCSV, PanelDisplayer.WIDGET_PlotGlobe_OrbitVariableToPlot.value, PanelDisplayer.WIDGET_PlotGlobe_OrbitColorbarTitle.value, PanelDisplayer.WIDGET_PlotGlobe_OrbitColorscaleName.value, PanelDisplayer.WIDGET_PlotGlobe_PlotTitle.value )
	Display( '    PlotGlobe( ' , IGRF_MeshgridResultCSV,', ' , PanelDisplayer.WIDGET_PlotGlobe_SurfaceVariableToPlot.value,', ' , PanelDisplayer.WIDGET_PlotGlobe_SurfaceColorbarTitle.value,', ' , PanelDisplayer.WIDGET_PlotGlobe_SurfaceColorscaleName.value,', ' , IGRF_OrbitResultCSV,', ' , PanelDisplayer.WIDGET_PlotGlobe_OrbitVariableToPlot.value,', ' , PanelDisplayer.WIDGET_PlotGlobe_OrbitColorbarTitle.value,', ' , PanelDisplayer.WIDGET_PlotGlobe_OrbitColorscaleName.value,', ' , PanelDisplayer.WIDGET_PlotGlobe_PlotTitle.value , ' ) results:')
	Display( '          -void-' )
	Display('')
"""
Executes these Simulation Modules:
  OrbitSelector  Interpolator  CreateNETCDFsphere  Interpolator  PlotGlobe  
"""
def Execute_Interpolator_NetCDF():
	Display('Executing Interpolator_NetCDF:')
	OrbitSelector_OrbitCSVfilename, OrbitSelector_OrbitNetCDFfilename = OrbitSelector( PanelDisplayer.WIDGET_OrbitSelector_Filename.value, PanelDisplayer.WIDGET_OrbitSelector_EvtXY.value, PanelDisplayer.WIDGET_OrbitSelector_TYP.value, PanelDisplayer.WIDGET_OrbitSelector_PerYYY.value, PanelDisplayer.WIDGET_OrbitSelector_LatZZ.value, PanelDisplayer.WIDGET_OrbitSelector_SRXXHZ.value, PanelDisplayer.WIDGET_OrbitSelector_SC.value )
	Display( '    OrbitSelector( ' , PanelDisplayer.WIDGET_OrbitSelector_Filename.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_EvtXY.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_TYP.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_PerYYY.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_LatZZ.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_SRXXHZ.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_SC.value , ' ) results:')
	Display( '          OrbitCSVfilename: ' + str(OrbitSelector_OrbitCSVfilename) )
	Display( '          OrbitNetCDFfilename: ' + str(OrbitSelector_OrbitNetCDFfilename) )
	Interpolator_InterpolatedOrbit = Interpolator( PanelDisplayer.WIDGET_Interpolator_ModelFile_forOrbit.value, OrbitSelector_OrbitNetCDFfilename )
	Display( '    Interpolator( ' , PanelDisplayer.WIDGET_Interpolator_ModelFile_forOrbit.value,', ' , OrbitSelector_OrbitNetCDFfilename , ' ) results:')
	Display( '          InterpolatedOrbit: ' + str(Interpolator_InterpolatedOrbit) )
	CreateNETCDFsphere_theSurfaceFilename = CreateNETCDFsphere( PanelDisplayer.WIDGET_CreateNETCDFsphere_fixedDatetimeString.value, PanelDisplayer.WIDGET_CreateNETCDFsphere_fixedAltitude.value, PanelDisplayer.WIDGET_CreateNETCDFsphere_LatitudeStep.value, PanelDisplayer.WIDGET_CreateNETCDFsphere_LongitudeStep.value )
	Display( '    CreateNETCDFsphere( ' , PanelDisplayer.WIDGET_CreateNETCDFsphere_fixedDatetimeString.value,', ' , PanelDisplayer.WIDGET_CreateNETCDFsphere_fixedAltitude.value,', ' , PanelDisplayer.WIDGET_CreateNETCDFsphere_LatitudeStep.value,', ' , PanelDisplayer.WIDGET_CreateNETCDFsphere_LongitudeStep.value , ' ) results:')
	Display( '          theSurfaceFilename: ' + str(CreateNETCDFsphere_theSurfaceFilename) )
	Interpolator_InterpolatedSurface = Interpolator( PanelDisplayer.WIDGET_Interpolator_ModelFile_forSurface.value, CreateNETCDFsphere_theSurfaceFilename )
	Display( '    Interpolator( ' , PanelDisplayer.WIDGET_Interpolator_ModelFile_forSurface.value,', ' , CreateNETCDFsphere_theSurfaceFilename , ' ) results:')
	Display( '          InterpolatedSurface: ' + str(Interpolator_InterpolatedSurface) )
	PlotGlobe( Interpolator_InterpolatedSurface, PanelDisplayer.WIDGET_PlotGlobe_SurfaceVariableToPlot.value, PanelDisplayer.WIDGET_PlotGlobe_SurfaceColorbarTitle.value, PanelDisplayer.WIDGET_PlotGlobe_SurfaceColorscaleName.value, Interpolator_InterpolatedOrbit, PanelDisplayer.WIDGET_PlotGlobe_OrbitVariableToPlot.value, PanelDisplayer.WIDGET_PlotGlobe_OrbitColorbarTitle.value, PanelDisplayer.WIDGET_PlotGlobe_OrbitColorscaleName.value, PanelDisplayer.WIDGET_PlotGlobe_PlotTitle.value )
	Display( '    PlotGlobe( ' , Interpolator_InterpolatedSurface,', ' , PanelDisplayer.WIDGET_PlotGlobe_SurfaceVariableToPlot.value,', ' , PanelDisplayer.WIDGET_PlotGlobe_SurfaceColorbarTitle.value,', ' , PanelDisplayer.WIDGET_PlotGlobe_SurfaceColorscaleName.value,', ' , Interpolator_InterpolatedOrbit,', ' , PanelDisplayer.WIDGET_PlotGlobe_OrbitVariableToPlot.value,', ' , PanelDisplayer.WIDGET_PlotGlobe_OrbitColorbarTitle.value,', ' , PanelDisplayer.WIDGET_PlotGlobe_OrbitColorscaleName.value,', ' , PanelDisplayer.WIDGET_PlotGlobe_PlotTitle.value , ' ) results:')
	Display( '          -void-' )
	Display('')
"""
Executes these Simulation Modules:
  Test1  
"""
def Execute_ExecuteTest1():
	Display('Executing ExecuteTest1:')
	Test1_txt = Test1(  )
	Display( '    Test1() results:')
	Display( '          txt: ' + str(Test1_txt) )
	Display('')
"""
Executes these Simulation Modules:
  OrbitSelector  HWM14  CreateCSV_Sphere  HWM14  PlotGlobe  
"""
def Execute_HWM14():
	Display('Executing HWM14:')
	OrbitSelector_OrbitCSVfilename, OrbitSelector_OrbitNetCDFfilename = OrbitSelector( PanelDisplayer.WIDGET_OrbitSelector_Filename.value, PanelDisplayer.WIDGET_OrbitSelector_EvtXY.value, PanelDisplayer.WIDGET_OrbitSelector_TYP.value, PanelDisplayer.WIDGET_OrbitSelector_PerYYY.value, PanelDisplayer.WIDGET_OrbitSelector_LatZZ.value, PanelDisplayer.WIDGET_OrbitSelector_SRXXHZ.value, PanelDisplayer.WIDGET_OrbitSelector_SC.value )
	Display( '    OrbitSelector( ' , PanelDisplayer.WIDGET_OrbitSelector_Filename.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_EvtXY.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_TYP.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_PerYYY.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_LatZZ.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_SRXXHZ.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_SC.value , ' ) results:')
	Display( '          OrbitCSVfilename: ' + str(OrbitSelector_OrbitCSVfilename) )
	Display( '          OrbitNetCDFfilename: ' + str(OrbitSelector_OrbitNetCDFfilename) )
	HWM14_OrbitResultCSV = HWM14( OrbitSelector_OrbitCSVfilename, PanelDisplayer.WIDGET_HWM14_Variable_forOrbit.value )
	Display( '    HWM14( ' , OrbitSelector_OrbitCSVfilename,', ' , PanelDisplayer.WIDGET_HWM14_Variable_forOrbit.value , ' ) results:')
	Display( '          OrbitResultCSV: ' + str(HWM14_OrbitResultCSV) )
	CreateCSV_Sphere_theCSVfilename, CreateCSV_Sphere_theAltitude, CreateCSV_Sphere_NumberOfLinesWritten = CreateCSV_Sphere( PanelDisplayer.WIDGET_CreateCSV_Sphere_CSVfilename.value, PanelDisplayer.WIDGET_CreateCSV_Sphere_fixedDatetimeString.value, PanelDisplayer.WIDGET_CreateCSV_Sphere_fixedAltitude.value, PanelDisplayer.WIDGET_CreateCSV_Sphere_LatitudeStep.value, PanelDisplayer.WIDGET_CreateCSV_Sphere_LongitudeStep.value )
	Display( '    CreateCSV_Sphere( ' , PanelDisplayer.WIDGET_CreateCSV_Sphere_CSVfilename.value,', ' , PanelDisplayer.WIDGET_CreateCSV_Sphere_fixedDatetimeString.value,', ' , PanelDisplayer.WIDGET_CreateCSV_Sphere_fixedAltitude.value,', ' , PanelDisplayer.WIDGET_CreateCSV_Sphere_LatitudeStep.value,', ' , PanelDisplayer.WIDGET_CreateCSV_Sphere_LongitudeStep.value , ' ) results:')
	Display( '          theCSVfilename: ' + str(CreateCSV_Sphere_theCSVfilename) )
	Display( '          theAltitude: ' + str(CreateCSV_Sphere_theAltitude) )
	Display( '          NumberOfLinesWritten: ' + str(CreateCSV_Sphere_NumberOfLinesWritten) )
	HWM14_MeshgridResultCSV = HWM14( CreateCSV_Sphere_theCSVfilename, PanelDisplayer.WIDGET_HWM14_Variable_forMeshgrid.value )
	Display( '    HWM14( ' , CreateCSV_Sphere_theCSVfilename,', ' , PanelDisplayer.WIDGET_HWM14_Variable_forMeshgrid.value , ' ) results:')
	Display( '          MeshgridResultCSV: ' + str(HWM14_MeshgridResultCSV) )
	PlotGlobe( HWM14_MeshgridResultCSV, PanelDisplayer.WIDGET_PlotGlobe_SurfaceVariableToPlot.value, PanelDisplayer.WIDGET_PlotGlobe_SurfaceColorbarTitle.value, PanelDisplayer.WIDGET_PlotGlobe_SurfaceColorscaleName.value, HWM14_OrbitResultCSV, PanelDisplayer.WIDGET_PlotGlobe_OrbitVariableToPlot.value, PanelDisplayer.WIDGET_PlotGlobe_OrbitColorbarTitle.value, PanelDisplayer.WIDGET_PlotGlobe_OrbitColorscaleName.value, PanelDisplayer.WIDGET_PlotGlobe_PlotTitle.value )
	Display( '    PlotGlobe( ' , HWM14_MeshgridResultCSV,', ' , PanelDisplayer.WIDGET_PlotGlobe_SurfaceVariableToPlot.value,', ' , PanelDisplayer.WIDGET_PlotGlobe_SurfaceColorbarTitle.value,', ' , PanelDisplayer.WIDGET_PlotGlobe_SurfaceColorscaleName.value,', ' , HWM14_OrbitResultCSV,', ' , PanelDisplayer.WIDGET_PlotGlobe_OrbitVariableToPlot.value,', ' , PanelDisplayer.WIDGET_PlotGlobe_OrbitColorbarTitle.value,', ' , PanelDisplayer.WIDGET_PlotGlobe_OrbitColorscaleName.value,', ' , PanelDisplayer.WIDGET_PlotGlobe_PlotTitle.value , ' ) results:')
	Display( '          -void-' )
	Display('')
"""
Executes these Simulation Modules:
  OrbitSelector  Selector  CreateCSV_Sphere  IRI16  IRI16  SubGridVariability  PlotGlobe  
"""
def Execute_IRI16():
	Display('Executing IRI16:')
	OrbitSelector_OrbitCSVfilename, OrbitSelector_OrbitNetCDFfilename = OrbitSelector( PanelDisplayer.WIDGET_OrbitSelector_Filename.value, PanelDisplayer.WIDGET_OrbitSelector_EvtXY.value, PanelDisplayer.WIDGET_OrbitSelector_TYP.value, PanelDisplayer.WIDGET_OrbitSelector_PerYYY.value, PanelDisplayer.WIDGET_OrbitSelector_LatZZ.value, PanelDisplayer.WIDGET_OrbitSelector_SRXXHZ.value, PanelDisplayer.WIDGET_OrbitSelector_SC.value )
	Display( '    OrbitSelector( ' , PanelDisplayer.WIDGET_OrbitSelector_Filename.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_EvtXY.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_TYP.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_PerYYY.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_LatZZ.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_SRXXHZ.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_SC.value , ' ) results:')
	Display( '          OrbitCSVfilename: ' + str(OrbitSelector_OrbitCSVfilename) )
	Display( '          OrbitNetCDFfilename: ' + str(OrbitSelector_OrbitNetCDFfilename) )
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
	IRI16_orbitCSVfilename = IRI16( OrbitSelector_OrbitCSVfilename, Selector_IRI16_variableForOrbit_OUT )
	Display( '    IRI16( ' , OrbitSelector_OrbitCSVfilename,', ' , Selector_IRI16_variableForOrbit_OUT , ' ) results:')
	Display( '          orbitCSVfilename: ' + str(IRI16_orbitCSVfilename) )
	SubGridVariability_VariabilityCSVfilename = SubGridVariability( Selector_NoiseWxOUT, Selector_NoiseWyOUT, Selector_NoiseAxOUT, Selector_NoiseAyOUT, IRI16_orbitCSVfilename, PanelDisplayer.WIDGET_SubGridVariability_ValueName.value )
	Display( '    SubGridVariability( ' , Selector_NoiseWxOUT,', ' , Selector_NoiseWyOUT,', ' , Selector_NoiseAxOUT,', ' , Selector_NoiseAyOUT,', ' , IRI16_orbitCSVfilename,', ' , PanelDisplayer.WIDGET_SubGridVariability_ValueName.value , ' ) results:')
	Display( '          VariabilityCSVfilename: ' + str(SubGridVariability_VariabilityCSVfilename) )
	PlotGlobe( IRI16_MeshgridResultCSV, PanelDisplayer.WIDGET_PlotGlobe_SurfaceVariableToPlot.value, PanelDisplayer.WIDGET_PlotGlobe_SurfaceColorbarTitle.value, PanelDisplayer.WIDGET_PlotGlobe_SurfaceColorscaleName.value, SubGridVariability_VariabilityCSVfilename, PanelDisplayer.WIDGET_PlotGlobe_OrbitVariableToPlot.value, PanelDisplayer.WIDGET_PlotGlobe_OrbitColorbarTitle.value, PanelDisplayer.WIDGET_PlotGlobe_OrbitColorscaleName.value, PanelDisplayer.WIDGET_PlotGlobe_PlotTitle.value )
	Display( '    PlotGlobe( ' , IRI16_MeshgridResultCSV,', ' , PanelDisplayer.WIDGET_PlotGlobe_SurfaceVariableToPlot.value,', ' , PanelDisplayer.WIDGET_PlotGlobe_SurfaceColorbarTitle.value,', ' , PanelDisplayer.WIDGET_PlotGlobe_SurfaceColorscaleName.value,', ' , SubGridVariability_VariabilityCSVfilename,', ' , PanelDisplayer.WIDGET_PlotGlobe_OrbitVariableToPlot.value,', ' , PanelDisplayer.WIDGET_PlotGlobe_OrbitColorbarTitle.value,', ' , PanelDisplayer.WIDGET_PlotGlobe_OrbitColorscaleName.value,', ' , PanelDisplayer.WIDGET_PlotGlobe_PlotTitle.value , ' ) results:')
	Display( '          -void-' )
	Display('')
"""
Executes these Simulation Modules:
  OrbitSelector  Selector  Interpolator  SGV_vi  SGV_Ti  SGV_Te  SGV_Ni  SGV_nix  SGV_Ne  SGV_TEC  SGV_un  SGV_Tn  SGV_Nn  SGV_nnx  SGV_JEPP  SGV_B  SGV_E  WakeEffectsCalc  ChargingEffectsCalc  InstrumentIDMRPA  InstrumentTII  InstrumentRWSCWS  InstrumentIMS  InstrumentNMS  InstrumentACC  InstrumentEPDS3  InstrumentEFI  InstrumentMIP  InstrumentMAG  InstrumentGNSS  DerivedProductsCalc  
"""
def Execute_TopLevel():
	Display('Executing TopLevel:')
	OrbitSelector_OrbitCSVfilename, OrbitSelector_OrbitNetCDFfilename = OrbitSelector( PanelDisplayer.WIDGET_OrbitSelector_Filename.value, PanelDisplayer.WIDGET_OrbitSelector_EvtXY.value, PanelDisplayer.WIDGET_OrbitSelector_TYP.value, PanelDisplayer.WIDGET_OrbitSelector_PerYYY.value, PanelDisplayer.WIDGET_OrbitSelector_LatZZ.value, PanelDisplayer.WIDGET_OrbitSelector_SRXXHZ.value, PanelDisplayer.WIDGET_OrbitSelector_SC.value )
	Display( '    OrbitSelector( ' , PanelDisplayer.WIDGET_OrbitSelector_Filename.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_EvtXY.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_TYP.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_PerYYY.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_LatZZ.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_SRXXHZ.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_SC.value , ' ) results:')
	Display( '          OrbitCSVfilename: ' + str(OrbitSelector_OrbitCSVfilename) )
	Display( '          OrbitNetCDFfilename: ' + str(OrbitSelector_OrbitNetCDFfilename) )
	Selector_EnableSubGridVariabilityOUT = Selector( PanelDisplayer.WIDGET_Selector_EnableSubGridVariability.value )
	Display( '    Selector( ' , PanelDisplayer.WIDGET_Selector_EnableSubGridVariability.value , ' ) results:')
	Display( '          EnableSubGridVariabilityOUT: ' + str(Selector_EnableSubGridVariabilityOUT) )
	Interpolator_InterpolatorResultFilename = Interpolator( PanelDisplayer.WIDGET_Interpolator_ModelFile.value, OrbitSelector_OrbitNetCDFfilename )
	Display( '    Interpolator( ' , PanelDisplayer.WIDGET_Interpolator_ModelFile.value,', ' , OrbitSelector_OrbitNetCDFfilename , ' ) results:')
	Display( '          InterpolatorResultFilename: ' + str(Interpolator_InterpolatorResultFilename) )
	SGV_vi_ResultFilename = SGV_vi( Interpolator_InterpolatorResultFilename, Selector_EnableSubGridVariabilityOUT )
	Display( '    SGV_vi( ' , Interpolator_InterpolatorResultFilename,', ' , Selector_EnableSubGridVariabilityOUT , ' ) results:')
	Display( '          ResultFilename: ' + str(SGV_vi_ResultFilename) )
	SGV_Ti_ResultFilename = SGV_Ti( Interpolator_InterpolatorResultFilename, Selector_EnableSubGridVariabilityOUT )
	Display( '    SGV_Ti( ' , Interpolator_InterpolatorResultFilename,', ' , Selector_EnableSubGridVariabilityOUT , ' ) results:')
	Display( '          ResultFilename: ' + str(SGV_Ti_ResultFilename) )
	SGV_Te_ResultFilename = SGV_Te( Interpolator_InterpolatorResultFilename, Selector_EnableSubGridVariabilityOUT )
	Display( '    SGV_Te( ' , Interpolator_InterpolatorResultFilename,', ' , Selector_EnableSubGridVariabilityOUT , ' ) results:')
	Display( '          ResultFilename: ' + str(SGV_Te_ResultFilename) )
	SGV_Ni_ResultFilename = SGV_Ni( Interpolator_InterpolatorResultFilename, Selector_EnableSubGridVariabilityOUT )
	Display( '    SGV_Ni( ' , Interpolator_InterpolatorResultFilename,', ' , Selector_EnableSubGridVariabilityOUT , ' ) results:')
	Display( '          ResultFilename: ' + str(SGV_Ni_ResultFilename) )
	SGV_nix_ResultFilename = SGV_nix( Interpolator_InterpolatorResultFilename, Selector_EnableSubGridVariabilityOUT )
	Display( '    SGV_nix( ' , Interpolator_InterpolatorResultFilename,', ' , Selector_EnableSubGridVariabilityOUT , ' ) results:')
	Display( '          ResultFilename: ' + str(SGV_nix_ResultFilename) )
	SGV_Ne_ResultFilename = SGV_Ne( Interpolator_InterpolatorResultFilename, Selector_EnableSubGridVariabilityOUT )
	Display( '    SGV_Ne( ' , Interpolator_InterpolatorResultFilename,', ' , Selector_EnableSubGridVariabilityOUT , ' ) results:')
	Display( '          ResultFilename: ' + str(SGV_Ne_ResultFilename) )
	SGV_TEC_ResultFilename = SGV_TEC( Interpolator_InterpolatorResultFilename, Selector_EnableSubGridVariabilityOUT )
	Display( '    SGV_TEC( ' , Interpolator_InterpolatorResultFilename,', ' , Selector_EnableSubGridVariabilityOUT , ' ) results:')
	Display( '          ResultFilename: ' + str(SGV_TEC_ResultFilename) )
	SGV_un_ResultFilename = SGV_un( Interpolator_InterpolatorResultFilename, Selector_EnableSubGridVariabilityOUT )
	Display( '    SGV_un( ' , Interpolator_InterpolatorResultFilename,', ' , Selector_EnableSubGridVariabilityOUT , ' ) results:')
	Display( '          ResultFilename: ' + str(SGV_un_ResultFilename) )
	SGV_Tn_ResultFilename = SGV_Tn( Interpolator_InterpolatorResultFilename, Selector_EnableSubGridVariabilityOUT )
	Display( '    SGV_Tn( ' , Interpolator_InterpolatorResultFilename,', ' , Selector_EnableSubGridVariabilityOUT , ' ) results:')
	Display( '          ResultFilename: ' + str(SGV_Tn_ResultFilename) )
	SGV_Nn_ResultFilename = SGV_Nn( Interpolator_InterpolatorResultFilename, Selector_EnableSubGridVariabilityOUT )
	Display( '    SGV_Nn( ' , Interpolator_InterpolatorResultFilename,', ' , Selector_EnableSubGridVariabilityOUT , ' ) results:')
	Display( '          ResultFilename: ' + str(SGV_Nn_ResultFilename) )
	SGV_nnx_ResultFilename = SGV_nnx( Interpolator_InterpolatorResultFilename, Selector_EnableSubGridVariabilityOUT )
	Display( '    SGV_nnx( ' , Interpolator_InterpolatorResultFilename,', ' , Selector_EnableSubGridVariabilityOUT , ' ) results:')
	Display( '          ResultFilename: ' + str(SGV_nnx_ResultFilename) )
	SGV_JEPP_ResultFilename = SGV_JEPP( Interpolator_InterpolatorResultFilename, Selector_EnableSubGridVariabilityOUT )
	Display( '    SGV_JEPP( ' , Interpolator_InterpolatorResultFilename,', ' , Selector_EnableSubGridVariabilityOUT , ' ) results:')
	Display( '          ResultFilename: ' + str(SGV_JEPP_ResultFilename) )
	SGV_B_ResultFilename = SGV_B( Interpolator_InterpolatorResultFilename, Selector_EnableSubGridVariabilityOUT )
	Display( '    SGV_B( ' , Interpolator_InterpolatorResultFilename,', ' , Selector_EnableSubGridVariabilityOUT , ' ) results:')
	Display( '          ResultFilename: ' + str(SGV_B_ResultFilename) )
	SGV_E_ResultFilename = SGV_E( Interpolator_InterpolatorResultFilename, Selector_EnableSubGridVariabilityOUT )
	Display( '    SGV_E( ' , Interpolator_InterpolatorResultFilename,', ' , Selector_EnableSubGridVariabilityOUT , ' ) results:')
	Display( '          ResultFilename: ' + str(SGV_E_ResultFilename) )
	WakeEffectsCalc_ResultFilename = WakeEffectsCalc( SGV_E_ResultFilename )
	Display( '    WakeEffectsCalc( ' , SGV_E_ResultFilename , ' ) results:')
	Display( '          ResultFilename: ' + str(WakeEffectsCalc_ResultFilename) )
	ChargingEffectsCalc_ResultFilename = ChargingEffectsCalc( SGV_E_ResultFilename )
	Display( '    ChargingEffectsCalc( ' , SGV_E_ResultFilename , ' ) results:')
	Display( '          ResultFilename: ' + str(ChargingEffectsCalc_ResultFilename) )
	InstrumentIDMRPA_ResultFilename = InstrumentIDMRPA( ChargingEffectsCalc_ResultFilename )
	Display( '    InstrumentIDMRPA( ' , ChargingEffectsCalc_ResultFilename , ' ) results:')
	Display( '          ResultFilename: ' + str(InstrumentIDMRPA_ResultFilename) )
	InstrumentTII_ResultFilename = InstrumentTII( ChargingEffectsCalc_ResultFilename )
	Display( '    InstrumentTII( ' , ChargingEffectsCalc_ResultFilename , ' ) results:')
	Display( '          ResultFilename: ' + str(InstrumentTII_ResultFilename) )
	InstrumentRWSCWS_ResultFilename = InstrumentRWSCWS( ChargingEffectsCalc_ResultFilename )
	Display( '    InstrumentRWSCWS( ' , ChargingEffectsCalc_ResultFilename , ' ) results:')
	Display( '          ResultFilename: ' + str(InstrumentRWSCWS_ResultFilename) )
	InstrumentIMS_ResultFilename = InstrumentIMS( ChargingEffectsCalc_ResultFilename )
	Display( '    InstrumentIMS( ' , ChargingEffectsCalc_ResultFilename , ' ) results:')
	Display( '          ResultFilename: ' + str(InstrumentIMS_ResultFilename) )
	InstrumentNMS_ResultFilename = InstrumentNMS( ChargingEffectsCalc_ResultFilename )
	Display( '    InstrumentNMS( ' , ChargingEffectsCalc_ResultFilename , ' ) results:')
	Display( '          ResultFilename: ' + str(InstrumentNMS_ResultFilename) )
	InstrumentACC_ResultFilename = InstrumentACC( ChargingEffectsCalc_ResultFilename )
	Display( '    InstrumentACC( ' , ChargingEffectsCalc_ResultFilename , ' ) results:')
	Display( '          ResultFilename: ' + str(InstrumentACC_ResultFilename) )
	InstrumentEPDS3_ResultFilename = InstrumentEPDS3( ChargingEffectsCalc_ResultFilename )
	Display( '    InstrumentEPDS3( ' , ChargingEffectsCalc_ResultFilename , ' ) results:')
	Display( '          ResultFilename: ' + str(InstrumentEPDS3_ResultFilename) )
	InstrumentEFI_ResultFilename = InstrumentEFI( ChargingEffectsCalc_ResultFilename )
	Display( '    InstrumentEFI( ' , ChargingEffectsCalc_ResultFilename , ' ) results:')
	Display( '          ResultFilename: ' + str(InstrumentEFI_ResultFilename) )
	InstrumentMIP_ResultFilename = InstrumentMIP( ChargingEffectsCalc_ResultFilename )
	Display( '    InstrumentMIP( ' , ChargingEffectsCalc_ResultFilename , ' ) results:')
	Display( '          ResultFilename: ' + str(InstrumentMIP_ResultFilename) )
	InstrumentMAG_ResultFilename = InstrumentMAG( ChargingEffectsCalc_ResultFilename )
	Display( '    InstrumentMAG( ' , ChargingEffectsCalc_ResultFilename , ' ) results:')
	Display( '          ResultFilename: ' + str(InstrumentMAG_ResultFilename) )
	InstrumentGNSS_ResultFilename = InstrumentGNSS( ChargingEffectsCalc_ResultFilename )
	Display( '    InstrumentGNSS( ' , ChargingEffectsCalc_ResultFilename , ' ) results:')
	Display( '          ResultFilename: ' + str(InstrumentGNSS_ResultFilename) )
	DerivedProductsCalc_ResultFilename = DerivedProductsCalc( InstrumentGNSS_ResultFilename )
	Display( '    DerivedProductsCalc( ' , InstrumentGNSS_ResultFilename , ' ) results:')
	Display( '          ResultFilename: ' + str(DerivedProductsCalc_ResultFilename) )
	Display('')
"""
Executes these Simulation Modules:
  OrbitSelector  Interpolator  PlotGlobe  
"""
def Execute_Interpolator():
	Display('Executing Interpolator:')
	OrbitSelector_OrbitCSVfilename, OrbitSelector_OrbitNetCDFfilename = OrbitSelector( PanelDisplayer.WIDGET_OrbitSelector_Filename.value, PanelDisplayer.WIDGET_OrbitSelector_EvtXY.value, PanelDisplayer.WIDGET_OrbitSelector_TYP.value, PanelDisplayer.WIDGET_OrbitSelector_PerYYY.value, PanelDisplayer.WIDGET_OrbitSelector_LatZZ.value, PanelDisplayer.WIDGET_OrbitSelector_SRXXHZ.value, PanelDisplayer.WIDGET_OrbitSelector_SC.value )
	Display( '    OrbitSelector( ' , PanelDisplayer.WIDGET_OrbitSelector_Filename.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_EvtXY.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_TYP.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_PerYYY.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_LatZZ.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_SRXXHZ.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_SC.value , ' ) results:')
	Display( '          OrbitCSVfilename: ' + str(OrbitSelector_OrbitCSVfilename) )
	Display( '          OrbitNetCDFfilename: ' + str(OrbitSelector_OrbitNetCDFfilename) )
	Interpolator_InterpolatedOrbitCSV = Interpolator( PanelDisplayer.WIDGET_Interpolator_model.value, PanelDisplayer.WIDGET_Interpolator_model_data_file.value, OrbitSelector_OrbitCSVfilename, PanelDisplayer.WIDGET_Interpolator_save.value, PanelDisplayer.WIDGET_Interpolator_VAR.value )
	Display( '    Interpolator( ' , PanelDisplayer.WIDGET_Interpolator_model.value,', ' , PanelDisplayer.WIDGET_Interpolator_model_data_file.value,', ' , OrbitSelector_OrbitCSVfilename,', ' , PanelDisplayer.WIDGET_Interpolator_save.value,', ' , PanelDisplayer.WIDGET_Interpolator_VAR.value , ' ) results:')
	Display( '          InterpolatedOrbitCSV: ' + str(Interpolator_InterpolatedOrbitCSV) )
	PlotGlobe( PanelDisplayer.WIDGET_PlotGlobe_SurfaceFilename.value, PanelDisplayer.WIDGET_PlotGlobe_SurfaceVariableToPlot.value, PanelDisplayer.WIDGET_PlotGlobe_SurfaceColorbarTitle.value, PanelDisplayer.WIDGET_PlotGlobe_SurfaceColorscaleName.value, Interpolator_InterpolatedOrbitCSV, PanelDisplayer.WIDGET_PlotGlobe_OrbitVariableToPlot.value, PanelDisplayer.WIDGET_PlotGlobe_OrbitColorbarTitle.value, PanelDisplayer.WIDGET_PlotGlobe_OrbitColorscaleName.value, PanelDisplayer.WIDGET_PlotGlobe_PlotTitle.value )
	Display( '    PlotGlobe( ' , PanelDisplayer.WIDGET_PlotGlobe_SurfaceFilename.value,', ' , PanelDisplayer.WIDGET_PlotGlobe_SurfaceVariableToPlot.value,', ' , PanelDisplayer.WIDGET_PlotGlobe_SurfaceColorbarTitle.value,', ' , PanelDisplayer.WIDGET_PlotGlobe_SurfaceColorscaleName.value,', ' , Interpolator_InterpolatedOrbitCSV,', ' , PanelDisplayer.WIDGET_PlotGlobe_OrbitVariableToPlot.value,', ' , PanelDisplayer.WIDGET_PlotGlobe_OrbitColorbarTitle.value,', ' , PanelDisplayer.WIDGET_PlotGlobe_OrbitColorscaleName.value,', ' , PanelDisplayer.WIDGET_PlotGlobe_PlotTitle.value , ' ) results:')
	Display( '          -void-' )
	Display('')
"""
Executes these Simulation Modules:
  OrbitSelector  MSISE00  CreateCSV_Sphere  MSISE00  PlotGlobe  
"""
def Execute_MSISE00():
	Display('Executing MSISE00:')
	OrbitSelector_OrbitCSVfilename, OrbitSelector_OrbitNetCDFfilename = OrbitSelector( PanelDisplayer.WIDGET_OrbitSelector_Filename.value, PanelDisplayer.WIDGET_OrbitSelector_EvtXY.value, PanelDisplayer.WIDGET_OrbitSelector_TYP.value, PanelDisplayer.WIDGET_OrbitSelector_PerYYY.value, PanelDisplayer.WIDGET_OrbitSelector_LatZZ.value, PanelDisplayer.WIDGET_OrbitSelector_SRXXHZ.value, PanelDisplayer.WIDGET_OrbitSelector_SC.value )
	Display( '    OrbitSelector( ' , PanelDisplayer.WIDGET_OrbitSelector_Filename.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_EvtXY.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_TYP.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_PerYYY.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_LatZZ.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_SRXXHZ.value,', ' , PanelDisplayer.WIDGET_OrbitSelector_SC.value , ' ) results:')
	Display( '          OrbitCSVfilename: ' + str(OrbitSelector_OrbitCSVfilename) )
	Display( '          OrbitNetCDFfilename: ' + str(OrbitSelector_OrbitNetCDFfilename) )
	MSISE00_ObitResultCSV = MSISE00( OrbitSelector_OrbitCSVfilename, PanelDisplayer.WIDGET_MSISE00_Variable_forOrbit.value )
	Display( '    MSISE00( ' , OrbitSelector_OrbitCSVfilename,', ' , PanelDisplayer.WIDGET_MSISE00_Variable_forOrbit.value , ' ) results:')
	Display( '          ObitResultCSV: ' + str(MSISE00_ObitResultCSV) )
	CreateCSV_Sphere_theCSVfilename, CreateCSV_Sphere_theAltitude, CreateCSV_Sphere_NumberOfLinesWritten = CreateCSV_Sphere( PanelDisplayer.WIDGET_CreateCSV_Sphere_CSVfilename.value, PanelDisplayer.WIDGET_CreateCSV_Sphere_fixedDatetimeString.value, PanelDisplayer.WIDGET_CreateCSV_Sphere_fixedAltitude.value, PanelDisplayer.WIDGET_CreateCSV_Sphere_LatitudeStep.value, PanelDisplayer.WIDGET_CreateCSV_Sphere_LongitudeStep.value )
	Display( '    CreateCSV_Sphere( ' , PanelDisplayer.WIDGET_CreateCSV_Sphere_CSVfilename.value,', ' , PanelDisplayer.WIDGET_CreateCSV_Sphere_fixedDatetimeString.value,', ' , PanelDisplayer.WIDGET_CreateCSV_Sphere_fixedAltitude.value,', ' , PanelDisplayer.WIDGET_CreateCSV_Sphere_LatitudeStep.value,', ' , PanelDisplayer.WIDGET_CreateCSV_Sphere_LongitudeStep.value , ' ) results:')
	Display( '          theCSVfilename: ' + str(CreateCSV_Sphere_theCSVfilename) )
	Display( '          theAltitude: ' + str(CreateCSV_Sphere_theAltitude) )
	Display( '          NumberOfLinesWritten: ' + str(CreateCSV_Sphere_NumberOfLinesWritten) )
	MSISE00_MeshgridResultCSV = MSISE00( CreateCSV_Sphere_theCSVfilename, PanelDisplayer.WIDGET_MSISE00_Variable_forMeshgrid.value )
	Display( '    MSISE00( ' , CreateCSV_Sphere_theCSVfilename,', ' , PanelDisplayer.WIDGET_MSISE00_Variable_forMeshgrid.value , ' ) results:')
	Display( '          MeshgridResultCSV: ' + str(MSISE00_MeshgridResultCSV) )
	PlotGlobe( MSISE00_MeshgridResultCSV, PanelDisplayer.WIDGET_PlotGlobe_SurfaceVariableToPlot.value, PanelDisplayer.WIDGET_PlotGlobe_SurfaceColorbarTitle.value, PanelDisplayer.WIDGET_PlotGlobe_SurfaceColorscaleName.value, MSISE00_ObitResultCSV, PanelDisplayer.WIDGET_PlotGlobe_OrbitVariableToPlot.value, PanelDisplayer.WIDGET_PlotGlobe_OrbitColorbarTitle.value, PanelDisplayer.WIDGET_PlotGlobe_OrbitColorscaleName.value, PanelDisplayer.WIDGET_PlotGlobe_PlotTitle.value )
	Display( '    PlotGlobe( ' , MSISE00_MeshgridResultCSV,', ' , PanelDisplayer.WIDGET_PlotGlobe_SurfaceVariableToPlot.value,', ' , PanelDisplayer.WIDGET_PlotGlobe_SurfaceColorbarTitle.value,', ' , PanelDisplayer.WIDGET_PlotGlobe_SurfaceColorscaleName.value,', ' , MSISE00_ObitResultCSV,', ' , PanelDisplayer.WIDGET_PlotGlobe_OrbitVariableToPlot.value,', ' , PanelDisplayer.WIDGET_PlotGlobe_OrbitColorbarTitle.value,', ' , PanelDisplayer.WIDGET_PlotGlobe_OrbitColorscaleName.value,', ' , PanelDisplayer.WIDGET_PlotGlobe_PlotTitle.value , ' ) results:')
	Display( '          -void-' )
	Display('')
