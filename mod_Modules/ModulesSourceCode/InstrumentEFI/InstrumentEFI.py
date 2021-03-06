import SourceCode.DaedalusGlobals as DaedalusGlobals
from datetime import datetime
from netCDF4 import Dataset
from shutil import copyfile
import os

"""
Please write here a description of the module's inputs, outputs and function
"""
def InstrumentEFI( InputFilename ):
    # copy input file to your dedicated folder (if not exists already) in order to process it
    simple_InputFilename = InputFilename[ InputFilename.rfind("/")+1 : InputFilename.rfind("_") ]
    ResultFilename = DaedalusGlobals.SyntheticScience_Files_Path + simple_InputFilename + "_SyntheticScience" + ".nc"
    if os.path.isfile(ResultFilename) == False :
        copyfile(InputFilename, ResultFilename)

    # open the file and store your calculations
    resultCDF=Dataset(ResultFilename, "a")
    '''
    Here you can write source code for storing calculations into the result NetCDF file
    '''
    resultCDF.EditTime   = str(datetime.now())
    resultCDF.close()    
        
    ##
    return ResultFilename
