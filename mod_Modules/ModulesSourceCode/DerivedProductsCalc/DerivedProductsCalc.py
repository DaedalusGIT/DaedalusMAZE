import SourceCode.DaedalusGlobals as DaedalusGlobals
from datetime import datetime
from netCDF4 import Dataset
from shutil import copyfile
import os

"""
Please write here a description of the module's inputs, outputs and function
"""
def DerivedProductsCalc( InputFilename ):
    ResultFilename = InputFilename
    if os.path.isfile(ResultFilename):
        # open the file and store your calculations
        resultCDF=Dataset(ResultFilename, "a")
        '''
        Here you can write source code for storing calculations into the result NetCDF file
        '''
        resultCDF.EditTime   = str(datetime.now())
        resultCDF.close()    
    ##
    return ResultFilename
