import numpy as np
import pandas as pd
import netCDF4
from netCDF4 import Dataset 
import datetime
import glob
import os
import random

MLT_duration_of_a_bucket = 1
MLAT_degrees_of_a_bucket = 2.5
ALT_distance_of_a_bucket   = 5
MLT_min = 0
MLT_max = 24
MLAT_min = -90
MLAT_max = +90
ALT_min = 100
ALT_max = 500
num_of_KP_bins = 2
KPsequence     = [ 0, 3 ] 
ALTsequence    = list( range( ALT_min,    ALT_max,    ALT_distance_of_a_bucket  ) )
MLATsequence = list( np.arange( MLAT_min, MLAT_max, MLAT_degrees_of_a_bucket) )
MLTsequence    = list( range( MLT_min,    MLT_max,    MLT_duration_of_a_bucket ) )

ResultFilename = ""
TypeOfCalculation = ""
Progress = 0  # integer 0 to 100
DistributionNumOfSlots = 0

def setDataParams( _MLT_min, _MLT_max, _MLT_duration_of_a_bucket, _MLAT_min, _MLAT_max, _MLAT_degrees_of_a_bucket, _ALT_min, _ALT_max, _ALT_distance_of_a_bucket, _num_of_KP_bins,_TypeOfCalculation, _DistributionNumOfSlots ):
    global MLT_min, MLT_max, MLT_duration_of_a_bucket, ALT_min, ALT_max, ALT_distance_of_a_bucket, MLAT_min, MLAT_max, MLAT_degrees_of_a_bucket, num_of_KP_bins,  TypeOfCalculation, DistributionNumOfSlots, KPsequence, ALTsequence, MLATsequence, MLTsequence, ResultFilename
    ####
    MLT_duration_of_a_bucket   = _MLT_duration_of_a_bucket
    MLAT_degrees_of_a_bucket = _MLAT_degrees_of_a_bucket
    ALT_distance_of_a_bucket   = _ALT_distance_of_a_bucket
    MLT_min = _MLT_min
    MLT_max = _MLT_max
    MLAT_min = _MLAT_min
    MLAT_max = _MLAT_max
    ALT_min = _ALT_min
    ALT_max = _ALT_max
    num_of_KP_bins = _num_of_KP_bins
    MLTsequence    = list( np.arange( MLT_min,    MLT_max,    MLT_duration_of_a_bucket ) )
    MLATsequence = list( np.arange( MLAT_min, MLAT_max, MLAT_degrees_of_a_bucket) )
    ALTsequence    = list( np.arange( ALT_min,    ALT_max,    ALT_distance_of_a_bucket  ) )
    if num_of_KP_bins == 1:
        KPsequence     = [ 0 ] 
    elif num_of_KP_bins == 2:
        KPsequence     = [ 0, 3 ] 
    elif num_of_KP_bins == 3:    
        KPsequence     = [ 0, 2, 4 ] 
    #
    TypeOfCalculation = _TypeOfCalculation
    DistributionNumOfSlots = _DistributionNumOfSlots
    # construct the results filename and check if exists so that you do not overwrite it
    ResultFilename = "results/" + TypeOfCalculation + "__"
    ResultFilename += "MLT" + "_" + str(MLT_min) + "_" + str(MLT_max) + "_" + str(MLT_duration_of_a_bucket) + "_"
    ResultFilename += "MLAT" + "_" + str(MLAT_min) + "_" + str(MLAT_max) + "_" + str(MLAT_degrees_of_a_bucket) + "_"
    ResultFilename += "ALT" + "_" + str(ALT_min) + "_" + str(ALT_max) + "_" + str(ALT_distance_of_a_bucket) + "_"
    ResultFilename += "Kp" + str(num_of_KP_bins) + "Bins"
    ResultFilename += ".nc"
    ResultFilename = ResultFilename.replace(".0", "")
    
    
                    
def init_ResultDataStructure():
    Buckets = dict()
    for aKP in KPsequence:
        for anALT in ALTsequence:
            for aMagLat in MLATsequence:
                for aMLT in MLTsequence:
                    Buckets[(aKP, anALT, aMagLat, aMLT, "Sum")] = 0
                    Buckets[(aKP, anALT, aMagLat, aMLT, "Len")] = 0
                    Buckets[(aKP, anALT, aMagLat, aMLT, "Vals")] = list()
                    Buckets[(aKP, anALT, aMagLat, aMLT, "Percentile10")] = 0
                    Buckets[(aKP, anALT, aMagLat, aMLT, "Percentile25")] = 0
                    Buckets[(aKP, anALT, aMagLat, aMLT, "Percentile50")] = 0
                    Buckets[(aKP, anALT, aMagLat, aMLT, "Percentile75")] = 0
                    Buckets[(aKP, anALT, aMagLat, aMLT, "Percentile90")] = 0
                    Buckets[(aKP, anALT, aMagLat, aMLT, "Variance")] = 0
                    Buckets[(aKP, anALT, aMagLat, aMLT, "Minimum")] = 0
                    Buckets[(aKP, anALT, aMagLat, aMLT, "Maximum")] = 0
                    Buckets[(aKP, anALT, aMagLat, aMLT, "Distribution")] = [0] * DistributionNumOfSlots
                    
    return Buckets



def LocatePositionInBuckets( aKp, anALT, aMagLat, aMLT ):
    kp_to_fall = alt_to_fall = maglat_to_fall = mlt_to_fall = None
    # find correct Alt
    for tmp in ALTsequence:
        if anALT>=tmp and anALT<tmp+ALT_distance_of_a_bucket:
            alt_to_fall=tmp
            break
    if alt_to_fall is None and anALT==ALTsequence[-1]+ALT_distance_of_a_bucket: alt_to_fall=ALTsequence[-1]
    # find correct kp
    if num_of_KP_bins == 1:
        kp_to_fall = 0
    elif num_of_KP_bins == 2:
        if aKp < 3: 
            kp_to_fall = 0
        else:
            kp_to_fall = 3
    elif num_of_KP_bins == 3:
        if aKp < 2: 
            kp_to_fall = 0
        elif aKp < 4: 
            kp_to_fall = 2
        else:
            kp_to_fall = 4
    # find correct MLT
    if MLTsequence[-1] < 24:
        for tmp in MLTsequence:
            if aMLT>=tmp and aMLT<tmp+MLT_duration_of_a_bucket: 
                mlt_to_fall=tmp
                break
        if mlt_to_fall is None and aMLT==MLTsequence[-1]+MLT_duration_of_a_bucket: mlt_to_fall=MLTsequence[-1] # for last position
    else:
        MLT_to_check = aMLT
        if MLT_to_check < MLTsequence[0]: MLT_to_check+=24
        for tmp in MLTsequence:
            if MLT_to_check>=tmp and MLT_to_check<tmp+MLT_duration_of_a_bucket: 
                mlt_to_fall=tmp
                break
        if mlt_to_fall is None and MLT_to_check==MLTsequence[-1]+MLT_duration_of_a_bucket: mlt_to_fall=MLTsequence[-1] # for last position
    # find correct MagLat
    for tmp in MLATsequence:
        if aMagLat>=tmp and aMagLat<tmp+MLAT_degrees_of_a_bucket:
            maglat_to_fall=tmp
            break
    if maglat_to_fall is None and aMagLat==MLATsequence[-1]+MLAT_degrees_of_a_bucket: maglat_to_fall=MLATsequence[-1]
    #
    return  kp_to_fall, alt_to_fall, maglat_to_fall, mlt_to_fall



def WriteResultsToCDF(ResultBuckets, ResultFilename, VariableName, Units):
    resultsCDF = Dataset( ResultFilename, 'w' )
    # add attributes defining the results file - TODO: add more attributes
    resultsCDF.DateOfCreation = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    resultsCDF.VariableName = VariableName
    resultsCDF.Units = Units
    resultsCDF.TypeOfCalculation = TypeOfCalculation
    # Create dimensions
    resultsCDF.createDimension( "dim_MLT", len(MLTsequence) )
    resultsCDF.createDimension( "dim_MLAT", len(MLATsequence) )
    resultsCDF.createDimension( "dim_KP", len(KPsequence) )
    resultsCDF.createDimension( "dim_ALT", len(ALTsequence) )
    resultsCDF.createDimension( "dim_distribution", DistributionNumOfSlots )
    # Create variables
    VAR_MLT = resultsCDF.createVariable( "MLT", "f4", "dim_MLT" )
    VAR_MLT[:] = MLTsequence
    VAR_MLAT = resultsCDF.createVariable( "MLAT", "f4", "dim_MLAT" )
    VAR_MLAT[:] = MLATsequence
    VAR_KP = resultsCDF.createVariable( "KP", "f4", "dim_KP" )
    VAR_KP[:] = KPsequence
    VAR_ALT = resultsCDF.createVariable( "ALT", "f4", "dim_ALT" )
    VAR_ALT[:] = ALTsequence
    #
    VAR_BinSums = resultsCDF.createVariable( "BinSums", "f4", ('dim_KP', 'dim_ALT', 'dim_MLAT', 'dim_MLT') )
    VAR_BinLens = resultsCDF.createVariable( "BinLens", "f4", ('dim_KP', 'dim_ALT', 'dim_MLAT', 'dim_MLT') )
    VAR_Percentile10=resultsCDF.createVariable( "Percentile10", "f4", ('dim_KP','dim_ALT','dim_MLAT','dim_MLT'))
    VAR_Percentile25=resultsCDF.createVariable( "Percentile25", "f4", ('dim_KP','dim_ALT','dim_MLAT','dim_MLT'))
    VAR_Percentile50=resultsCDF.createVariable( "Percentile50", "f4", ('dim_KP','dim_ALT','dim_MLAT','dim_MLT'))
    VAR_Percentile75=resultsCDF.createVariable( "Percentile75", "f4", ('dim_KP','dim_ALT','dim_MLAT','dim_MLT'))
    VAR_Percentile90=resultsCDF.createVariable( "Percentile90", "f4", ('dim_KP','dim_ALT','dim_MLAT','dim_MLT'))
    VAR_Variance=resultsCDF.createVariable( "Variance", "f4", ('dim_KP','dim_ALT','dim_MLAT','dim_MLT'))
    VAR_Minimum=resultsCDF.createVariable( "Minimum", "f4", ('dim_KP','dim_ALT','dim_MLAT','dim_MLT'))
    VAR_Maximum=resultsCDF.createVariable( "Maximum", "f4", ('dim_KP','dim_ALT','dim_MLAT','dim_MLT'))
    if DistributionNumOfSlots > 0:
        VAR_Distribution=resultsCDF.createVariable( "Distribution", "i4",('dim_KP','dim_ALT','dim_MLAT','dim_MLT','dim_distribution'))
    for aKP in KPsequence:
        for anALT in ALTsequence:
            for aMagLat in MLATsequence:
                for aMLT in MLTsequence:
                    vector = (KPsequence.index(aKP), ALTsequence.index(anALT), MLATsequence.index(aMagLat), MLTsequence.index(aMLT))
                    VAR_BinSums[vector] = ResultBuckets[(aKP, anALT, aMagLat, aMLT, "Sum")]
                    VAR_BinLens[vector] = ResultBuckets[(aKP, anALT, aMagLat, aMLT, "Len")]
                    VAR_Percentile10[vector] = ResultBuckets[(aKP, anALT, aMagLat, aMLT, "Percentile10")]
                    VAR_Percentile25[vector] = ResultBuckets[(aKP, anALT, aMagLat, aMLT, "Percentile25")]
                    VAR_Percentile50[vector] = ResultBuckets[(aKP, anALT, aMagLat, aMLT, "Percentile50")]
                    VAR_Percentile75[vector] = ResultBuckets[(aKP, anALT, aMagLat, aMLT, "Percentile75")]
                    VAR_Percentile90[vector] = ResultBuckets[(aKP, anALT, aMagLat, aMLT, "Percentile90")]
                    VAR_Variance[vector] = ResultBuckets[(aKP, anALT, aMagLat, aMLT, "Variance")]
                    VAR_Minimum[vector] = ResultBuckets[(aKP, anALT, aMagLat, aMLT, "Minimum")]
                    VAR_Maximum[vector] = ResultBuckets[(aKP, anALT, aMagLat, aMLT, "Maximum")]
                    if DistributionNumOfSlots > 0:
                        for i in range(0, DistributionNumOfSlots):
                            VAR_Distribution[vector+(i,)] = ResultBuckets[(aKP, anALT, aMagLat, aMLT, "Distribution")][i]
                    
    resultsCDF.close()
    

    
def LoadResultsCDF( _ResultFilename ):
    global MLT_min, MLT_max, MLT_duration_of_a_bucket
    global ALT_min, ALT_max, ALT_distance_of_a_bucket 
    global MLAT_min, MLAT_max, MLAT_degrees_of_a_bucket
    global num_of_KP_bins, ResultFilename, TypeOfCalculation
    global KPsequence, ALTsequence, MLATsequence, MLTsequence
    
    ResultFilename = _ResultFilename
    
    # open result file
    try:
        CDFroot = Dataset( ResultFilename, 'r' )
    except:
        print ( " !!!!!!!! WRONG FORMAT:", ResultFilename )

    # read atributes
    DateOfCreation = CDFroot.DateOfCreation
    VariableName = CDFroot.VariableName
    Units = CDFroot.Units
    TypeOfCalculation = CDFroot.TypeOfCalculation
    
    # read results
    MLTs    = CDFroot.variables['MLT'][:]         
    MLATs   = CDFroot.variables['MLAT'][:] 
    ALTs    = CDFroot.variables['ALT'][:]
    KPs     = CDFroot.variables['KP'][:]        
    BinSums = CDFroot.variables['BinSums'][:, :, :, :] 
    BinLens = CDFroot.variables['BinLens'][:, :, :, :]
    P10 = CDFroot.variables['Percentile10'][:, :, :, :]
    P25 = CDFroot.variables['Percentile25'][:, :, :, :]
    P50 = CDFroot.variables['Percentile50'][:, :, :, :]
    P75 = CDFroot.variables['Percentile75'][:, :, :, :]
    P90 = CDFroot.variables['Percentile90'][:, :, :, :]
    try:
        Variances = CDFroot.variables['Variance'][:, :, :, :]
        Minimums = CDFroot.variables['Minimum'][:, :, :, :]
        Maximums = CDFroot.variables['Maximum'][:, :, :, :]
        Distribution = CDFroot.variables['Distribution'][:, :, :, :, :]    
    except:
        pass
    
    # apply loaded info to data structures
    MLT_duration_of_a_bucket = MLTs[1] - MLTs[0]
    MLT_min = MLTs[0]
    MLT_max = MLTs[-1] + MLT_duration_of_a_bucket
    #
    if len(MLATs)==1:
        MLAT_degrees_of_a_bucket = 180
    else:
        MLAT_degrees_of_a_bucket = MLATs[1] - MLATs[0]
    MLAT_min = MLATs[0]
    MLAT_max = MLATs[-1] + MLAT_degrees_of_a_bucket
    #
    ALT_distance_of_a_bucket = ALTs[1] - ALTs[0]
    ALT_min = ALTs[0]
    ALT_max = ALTs[-1] + ALT_distance_of_a_bucket
    #
    num_of_KP_bins = len(KPs)

    # reconstruct sequences
    MLTsequence  = list( np.arange( MLT_min,  MLT_max,  MLT_duration_of_a_bucket ) )
    MLATsequence = list( np.arange( MLAT_min, MLAT_max, MLAT_degrees_of_a_bucket) )
    ALTsequence  = list( np.arange( ALT_min,  ALT_max,  ALT_distance_of_a_bucket  ) )
    
    if num_of_KP_bins == 1:
        KPsequence     = [ 0 ] 
    elif num_of_KP_bins == 2:
        KPsequence     = [ 0, 3 ] 
    elif num_of_KP_bins == 3:    
        KPsequence     = [ 0, 2, 4 ] 
        
    # assign data to the buckets
    ResultBuckets = init_ResultDataStructure()
    for idx_kp in range(0, len(BinSums)):
        for idx_alt in range(0, len(BinSums[0])):
            for idx_mlat in range(0, len(BinSums[0,0])):
                for idx_mlt in range(0, len(BinSums[0,0,0])):
                    aKP     = KPs[idx_kp]
                    anALT   = ALTs[idx_alt]
                    aMagLat = MLATs[idx_mlat]
                    aMLT    = MLTs[idx_mlt]
                    ResultBuckets[(aKP, anALT, aMagLat, aMLT, "Sum")] = BinSums[idx_kp, idx_alt, idx_mlat, idx_mlt]
                    ResultBuckets[(aKP, anALT, aMagLat, aMLT, "Len")] = BinLens[idx_kp, idx_alt, idx_mlat, idx_mlt]
                    ResultBuckets[(aKP, anALT, aMagLat, aMLT, "Percentile10")] = P10[idx_kp, idx_alt, idx_mlat, idx_mlt]
                    ResultBuckets[(aKP, anALT, aMagLat, aMLT, "Percentile25")] = P25[idx_kp, idx_alt, idx_mlat, idx_mlt]
                    ResultBuckets[(aKP, anALT, aMagLat, aMLT, "Percentile50")] = P50[idx_kp, idx_alt, idx_mlat, idx_mlt]
                    ResultBuckets[(aKP, anALT, aMagLat, aMLT, "Percentile75")] = P75[idx_kp, idx_alt, idx_mlat, idx_mlt]
                    ResultBuckets[(aKP, anALT, aMagLat, aMLT, "Percentile90")] = P90[idx_kp, idx_alt, idx_mlat, idx_mlt]
                    try:
                        ResultBuckets[(aKP, anALT, aMagLat, aMLT, "Variance")] = Variances[idx_kp, idx_alt, idx_mlat, idx_mlt]
                        ResultBuckets[(aKP, anALT, aMagLat, aMLT, "Minimum")] = Minimums[idx_kp, idx_alt, idx_mlat, idx_mlt]
                        ResultBuckets[(aKP, anALT, aMagLat, aMLT, "Maximum")] = Maximums[idx_kp, idx_alt, idx_mlat, idx_mlt]
                        
                        # for distribution
                        ResultBuckets[(aKP, anALT, aMagLat, aMLT, "Distribution")] = [0] * len(Distribution[0,0,0,0,:])
                        for i in range(0, len(Distribution[0,0,0,0,:])):
                            ResultBuckets[(aKP, anALT, aMagLat, aMLT, "Distribution")][i] = Distribution[idx_kp, idx_alt, idx_mlat, idx_mlt, i]
                    except:
                        pass
                    
    # verbose
    print( "Opened ", ResultFilename, ":"  )
    print("   DateOfCreation:", DateOfCreation, " TypeOfCalculation:", TypeOfCalculation )
    print("   Variable:", VariableName, " Units:", Units, " filled bins:", np.count_nonzero(BinSums), "/", len(BinSums)*len(BinSums[0])*len(BinSums[0][0])*len(BinSums[0][0][0]), " Num of Kp bins:", num_of_KP_bins)
    print("   MLT:", MLT_min, "-", MLT_max, "step", MLT_duration_of_a_bucket, "  Alt.:", ALT_min, "-", ALT_max, "step", ALT_distance_of_a_bucket, "  Mag.Lat.:", MLAT_min, "-", MLAT_max, "step", MLAT_degrees_of_a_bucket)
    print("\n")
    
    # clean up
    CDFroot.close()
    return ResultBuckets, BinSums, BinLens, VariableName, Units
    

'''
merges several cdf result files into a single one.
arguments:
    CDF_filenames: string with wildcards, describing the files to be merged
    mergedFilename: string with the filename of the final merged file
'''
def mergeResultCDFs( CDF_filenames, mergedFilename ):
    AllCDFfilenames = sorted( glob.glob( CDF_filenames ) )
    MergedBuckets = init_ResultDataStructure()
    for CDF_filename in AllCDFfilenames:
        # read a file
        ResultBuckets, BinSums, BinLens, VariableName, Units = LoadResultsCDF( CDF_filename )
        # merge it 
        for aKP in KPsequence:
            for anALT in ALTsequence:
                for aMagLat in MLATsequence:
                    for aMLT in MLTsequence:
                        MergedBuckets[(aKP, anALT, aMagLat, aMLT, "Sum")] += ResultBuckets[(aKP, anALT, aMagLat, aMLT, "Sum")]
                        MergedBuckets[(aKP, anALT, aMagLat, aMLT, "Len")] += ResultBuckets[(aKP, anALT, aMagLat, aMLT, "Len")]
        # delete file
        os.remove( CDF_filename )
    # save merged data
    WriteResultsToCDF(MergedBuckets, mergedFilename, VariableName, Units)
        
        
        
    
def TT():
        df = pd.DataFrame()
        df["Alts"] = [1,2,3,1,2,3,5,7,3]
        df["MLTs"] = [1,1,1,2,2,2,1,1,1]
        df["JH"]   = [0.11,0.2,0.3,0.1,0.2,0.3,0.5,0.7,0.3]
        print(df)
        print("-- - - - - \n")

        #df["bins"] = df.groupby('JH').Alts.apply(pd.cut, [0,1,2,3], include_lowest=True)  #labels=["One", "Two", "Three"]
        #print( df["bins"] )
        #print("-- - - - - \n")
        #print(df.groupby(['bins', 'JH']).size().unstack() )
        #df.groupby(['bins_X', 'S']).size().unstack()
        #print( type(df["bins"]), type(df["bins"][0]) )
        #print( df.groupby(df["bins"])['JH'].agg(['count', 'sum']) )

        df = df.assign(
            Alt_cut = df.groupby('JH').Alts.apply(pd.cut, [0,1,2,3], include_lowest=True),
            MLT_cut = df.groupby('JH').MLTs.apply(pd.cut, [0,1,2], include_lowest=True)
        )
        #print(d1)
        #d2 = d1.assign(cartesian=pd.Categorical(d1.filter(regex='_cut').apply(tuple, 1)))
        #print(d2)
        print( df.groupby(["Alt_cut","MLT_cut"])['JH'].agg(['count', 'sum']) )
        print(df)
    
#def appendTo_AllValuesBuckets( kp_to_fall, alt_to_fall, maglat_to_fall, mlt_to_fall, aValue ):
#    global AllValuesBuckets
#    AllValuesBuckets[ kp_to_fall, alt_to_fall, maglat_to_fall, mlt_to_fall, "Vals" ].append( aValue )


