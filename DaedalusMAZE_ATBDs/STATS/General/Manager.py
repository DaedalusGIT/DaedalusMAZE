'''
TODO: 
  read orbits
  store tmp as bytes not as strings
  rm tmp/* after finish
'''

# local imports
import Data

# system imports
import netCDF4
from netCDF4 import Dataset 
import os
import datetime
import time
import glob
import shutil
import math
import numpy as np
import multiprocessing
from pathlib import Path
import random
import pyglow
from array import array

Test = False
NUM_OF_PROCESSORS = 14
TMP_FOLDER = "/media/balukid/STATStmp/" #"results/tmp/"

def StartCalculating( NetCDF_files_path, ResultFilename, TypeOfCalculation ):
    startSecs = time.time()
    print( "START", datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S") )
    
    Allprocesses = list()
    AllCDFfiles = sorted( glob.glob( NetCDF_files_path, recursive=True ) )
    print( datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S") )
    print( "I will calculate '" + TypeOfCalculation + "' on", len(AllCDFfiles), "files:\n    ", NetCDF_files_path, "\n" )
    print( "Results will be stored in '" + ResultFilename + "'\n" )
    
    # DO NOT DELL PARTIAL FILES - THEY CAN BE USED TO CONTINUE CALCULATING AFTER AN INTERMEDIATE HALT
    # del older partial txt files - there is one file for each bucket containing all values in it
    #try:
    #    shutil.rmtree( TMP_FOLDER )
    #except:
    #    pass
    
    n = 0
    for CDF_file in AllCDFfiles:
        n += 1
        Data.Progress = int( 100 * n/221)
        if Test and "004" in CDF_file: break
        
        # spawn new process
        P = multiprocessing.Process(target=PROC_StatsCalculator, args=(n,CDF_file,TypeOfCalculation))
        Allprocesses.append(P)
        P.start()
        
        pause_spawning = True
        while pause_spawning:
            Num_of_alive_processes = 0        
            for P in Allprocesses:
                if P.is_alive():
                    Num_of_alive_processes += 1            
            if Num_of_alive_processes >= NUM_OF_PROCESSORS:
                pause_spawning = True
                time.sleep(12)
            else:
                pause_spawning = False
        
           
    # wait for all processes to terminate
    for T in Allprocesses: T.join()
        
    # every process creates a partial file, merge all of them into one
    #print( "Merging partial CDF files...",  datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S"))
    #Data.mergeResultCDFs( ResultFilename+".part*", ResultFilename )
    
    print( "Merging partial data files and calculating result values...",  datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S"))
    ResultBuckets = Data.init_ResultDataStructure()
    NumOfBins = len(Data.KPsequence) * len(Data.ALTsequence) * len(Data.MLATsequence) * len(Data.MLTsequence)
    CurrBinNum = 0
    for aKP in Data.KPsequence:
        for anALT in Data.ALTsequence:
            for aMagLat in Data.MLATsequence:
                for aMLT in Data.MLTsequence:
                    CurrBinNum += 1
                    Data.Progress = int( 100 * CurrBinNum/NumOfBins )
                    AllBinValues = list()
                    for i in range(1,222): # read all partial files for this bin 
                        partialDataFolder = TMP_FOLDER+"proc"+ f"{i:03}" +"/"
                        if os.path.isdir(partialDataFolder)==False:
                            #print( "There are no partial data files for process", i )
                            continue
                        partialTextFilename = partialDataFolder + str(aKP)+"_"+str(anALT)+"_"+str(aMagLat)+"_"+str(aMLT)+".txt"
                        if os.path.exists(partialTextFilename) == False: # no hits for this bin from this process
                            continue
                            
                        # FOR READING VALUES STORED AS TEXT:
                        #    f = open(partialTextFilename, "r")
                        #    line = f.read()
                        #    f.close()
                        #    AllBinValues += [float(s) for s in line.split(' ')]
                        f = open(partialTextFilename, "rb")
                        float_array = array('d')
                        float_array.frombytes(f.read())
                        AllBinValues = float_array.tolist()
                        f.close()
                        
                    if len(AllBinValues) > 0:
                        ResultBuckets[aKP, anALT, aMagLat, aMLT, "Sum"] = np.sum(AllBinValues)
                        ResultBuckets[aKP, anALT, aMagLat, aMLT, "Len"] = len(AllBinValues)
                        ResultBuckets[aKP, anALT, aMagLat, aMLT, "Percentile10"] = np.percentile(AllBinValues, 10)
                        ResultBuckets[aKP, anALT, aMagLat, aMLT, "Percentile25"] = np.percentile(AllBinValues, 25)
                        ResultBuckets[aKP, anALT, aMagLat, aMLT, "Percentile50"] = np.percentile(AllBinValues, 50)
                        ResultBuckets[aKP, anALT, aMagLat, aMLT, "Percentile75"] = np.percentile(AllBinValues, 75)
                        ResultBuckets[aKP, anALT, aMagLat, aMLT, "Percentile90"] = np.percentile(AllBinValues, 90)
                        ResultBuckets[aKP, anALT, aMagLat, aMLT, "Variance"] = np.var(AllBinValues)
                        ResultBuckets[aKP, anALT, aMagLat, aMLT, "Minimum"] = np.nanmin(AllBinValues)
                        ResultBuckets[aKP, anALT, aMagLat, aMLT, "Maximum"] = np.nanmax(AllBinValues)
                        
                        # calculate distribution
                        if Data.DistributionNumOfSlots > 0:
                            histo_values, histo_ranges = np.histogram(AllBinValues, bins=Data.DistributionNumOfSlots)
                            for i in range(0, Data.DistributionNumOfSlots):
                                ResultBuckets[aKP, anALT, aMagLat, aMLT, "Distribution"][i] = histo_values[i]
        
    Data.WriteResultsToCDF(ResultBuckets, ResultFilename, "Joule Heating", "W/m3")

    # DO NOT DELL PARTIAL FILES - THEY CAN BE USED TO CONTINUE CALCULATING AFTER AN INTERMEDIATE HALT
    # delete temporary files, which contain all values for each bin
    #try:
    #    shutil.rmtree( TMP_FOLDER )
    #except:
    #    pass
    
    # 
    finishSecs = time.time()
    print( "FINISH",  datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S"), " (", finishSecs-startSecs, "sec )")


    
    

    
                  
                  
    

'''
Reads a NetCDF file and saves all the values of the variable in files.
The variable is chosen by the <TypeOfCalculation> argument
The process saves several files in its own folder with name: TMP_FOLDER+"proc"+<ProcessNum>+"/"
The folder contains one binary file for each bin. The file contains all values of the variable which fall in the bin
'''
def PROC_StatsCalculator(ProcessNum, CDF_filename, TypeOfCalculation):
    # partialResultFilename IS OBSOLETE, TMP_FOLDER is used to save all values
    #partialResultFilename = Data.ResultFilename+".part"+str(ProcessNum)
    #if os.path.exists(partialResultFilename):
    #    print( "Process",ProcessNum, ":", partialResultFilename, "exists. I am not needed." )
    #    return # <<<<
    
    # check if the data of this process have already been calculated
    procfolder = TMP_FOLDER+"proc"+ f"{ProcessNum:03}" +"/"
    if os.path.isdir(procfolder):
        print( "Data for file", ProcessNum, "already calculated.", "Process", ProcessNum, "finished." )
        return # <<<<
    else:
        if os.path.isdir(TMP_FOLDER)==False: os.mkdir( TMP_FOLDER )
        os.mkdir( procfolder )    
    
    # open netCDF file 
    ######## while( os.path.exists("ReadingFile.flag") ): # wait until no other process is reading from disk/NFS
    ########    time.sleep(random.randint(8,20))
    print("Process",ProcessNum,"reading ",CDF_filename[CDF_filename.rfind('/')+1:], datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S"))
    #Path("ReadingFile.flag").touch() # raise a flag that this process is now reading a file, so that other processes wait
    try:
        CDFroot = Dataset( CDF_filename, 'r' )
    except:
        print ( " !!!!!!!! WRONG FORMAT:", CDF_filename )
        #os.remove("ReadingFile.flag") # lower the reading-file flag
        return
        
    # read the data from the netCDF file
    #TIMEs  = CDFroot.variables['time'][:] 
    if "JH" in TypeOfCalculation:        Ohmics = CDFroot.variables['Ohmic'][:, :, :, :]  # m/s
    if "PedCond" in TypeOfCalculation:   PEDs   = CDFroot.variables['SIGMA_PED'][:, :, :, :]
    if "HallCond" in TypeOfCalculation:  HALs   = CDFroot.variables['SIGMA_HAL'][:, :, :, :]
    if "EEX_si" in TypeOfCalculation:    EEXs   = CDFroot.variables['EEX_si'][:, :, :, :]
    if "EEY_si" in TypeOfCalculation:    EEYs   = CDFroot.variables['EEY_si'][:, :, :, :]
    if "ConvHeat" in TypeOfCalculation:  ConvH  = CDFroot.variables['Convection_heating'][:, :, :, :]
    if "WindHeat" in TypeOfCalculation:  WindH  = CDFroot.variables['Wind_heating'][:, :, :, :]
    if "JHminusWindHeat" in TypeOfCalculation:  WindH  = CDFroot.variables['Wind_heating'][:, :, :, :]        
        
    if "JHnoWindsEISCAT" in TypeOfCalculation: 
        UI = CDFroot.variables['UI_ExB'][:, :, :, :]
        VI = CDFroot.variables['VI_ExB'][:, :, :, :]
        WI = CDFroot.variables['WI_ExB'][:, :, :, :]
        TIMEs = CDFroot.variables['time'][:] 
        PEDs  = CDFroot.variables['SIGMA_PED'][:, :, :, :]
        LONs  = CDFroot.variables['lon'][:] 
        ZGs   = CDFroot.variables['ZG'][:, :, :, :] / 100000 # Geometric height stored in cm, converted to km
    #
    if "EISCAT" in TypeOfCalculation:    
        LATs   = CDFroot.variables['lat'][:] 
    MLATs   = CDFroot.variables['mlat_qdf'][:, :, :, :] 
    MLTs    = CDFroot.variables['mlt_qdf'][:, :, :, :]         
    ALTs    = CDFroot.variables['ZGMID'][:, :, :, :] / 100000 # Geometric height stored in cm, converted to km
    KPs     = CDFroot.variables['Kp'][:]
    
    #try:
    #    os.remove("ReadingFile.flag") # lower the reading-file flag
    #except:
    #    pass

    ResultBuckets = Data.init_ResultDataStructure().copy()
    num_of_unbinned_items = 0
    step = 1
    for idx_time in range(0, len(ALTs), step):
        if idx_time%20==0: print( idx_time, "of", len(ALTs), CDF_filename[CDF_filename.rfind("sech_")+4:CDF_filename.rfind("_JH")+1]+"...nc", datetime.datetime.now().strftime("%H:%M:%S") )
        for idx_lev in range(0, len(ALTs[0]), step):
            for idx_lat in range(0, len(ALTs[0,0]), step):
                for idx_lon in range(0, len(ALTs[0,0,0]), step):
                    
                    curr_alt_km = ALTs[idx_time, idx_lev, idx_lat, idx_lon] 
                    
                    # ignore values for out-of-range positions 
                    if curr_alt_km<Data.ALT_min or curr_alt_km>Data.ALT_max:
                        continue
                    if "EISCAT" in TypeOfCalculation:    
                        if LATs[idx_lat]<71 or LATs[idx_lat]>78.5:  #if LATs[idx_lat] < 60: 
                            continue
                        
                    curr_kp     = KPs[idx_time]
                    curr_mlt    = MLTs[idx_time, idx_lev, idx_lat, idx_lon]
                    curr_maglat = MLATs[idx_time, idx_lev, idx_lat, idx_lon]
                        
                    kp_to_fall,alt_to_fall,maglat_to_fall,mlt_to_fall = Data.LocatePositionInBuckets(curr_kp,curr_alt_km,curr_maglat, curr_mlt)
    
                    if kp_to_fall is None or alt_to_fall is None or maglat_to_fall is None or mlt_to_fall is None:
                        num_of_unbinned_items += 1
                        #if num_of_unbinned_items < 4:
                        #    print( "WARNING: [", curr_kp, curr_alt_km, curr_maglat, curr_mlt, "] -> [", kp_to_fall, alt_to_fall, maglat_to_fall, mlt_to_fall, "] not binned" )
                        break
                    else:
                        if TypeOfCalculation=="JHminusWindHeat" or TypeOfCalculation=="JHminusWindHeatEISCAT":
                            aValue = Ohmics[idx_time, idx_lev, idx_lat, idx_lon] - WindH[idx_time, idx_lev, idx_lat, idx_lon]
                            if aValue > 100: continue # ignore faulty large values
                        elif TypeOfCalculation=="JHEISCAT" or TypeOfCalculation=="JH":
                            aValue = Ohmics[idx_time, idx_lev, idx_lat, idx_lon]
                            if aValue > 100: continue # ignore faulty large values
                        elif "PedCond" in TypeOfCalculation:
                            aValue = PEDs[idx_time, idx_lev, idx_lat, idx_lon]
                        elif "HallCond" in TypeOfCalculation:
                            aValue = HALs[idx_time, idx_lev, idx_lat, idx_lon]
                        elif "EEX_si" in TypeOfCalculation:
                            aValue = EEXs[idx_time, idx_lev, idx_lat, idx_lon]
                        elif "EEY_si" in TypeOfCalculation:
                            aValue = EEYs[idx_time, idx_lev, idx_lat, idx_lon]
                        elif "ConvHeat" in TypeOfCalculation:
                            aValue = ConvH[idx_time, idx_lev, idx_lat, idx_lon]
                            if aValue > 100: continue # ignore faulty large values
                        elif "WindHeat" in TypeOfCalculation:
                            aValue = WindH[idx_time, idx_lev, idx_lat, idx_lon]
                        elif TypeOfCalculation=="JHnoWindsEISCAT": 
                            I = list()
                            B = list()
                            time_p = datetime.datetime(2015, 3, 15, 0, 0, 0) + datetime.timedelta(minutes=TIMEs[idx_time])
                            lat_p = LATs[idx_lat]
                            lon_p = LONs[idx_lon]
                            alt_p = ZGs[idx_time, idx_lev, idx_lat, idx_lon]
                            pt = pyglow.Point(time_p, lat_p, lon_p, alt_p) # pyglow igrf
                            pt.run_igrf()
                            B.append( pt.Bx )  # Be, Tesla  (si)
                            B.append( pt.By )  # Bn, Tesla  (si)
                            B.append( pt.Bz )  # Bu, Tesla  (si)
                            I.append( UI[idx_time,idx_lev,idx_lat,idx_lon] )
                            I.append( VI[idx_time,idx_lev,idx_lat,idx_lon] )
                            I.append( WI[idx_time,idx_lev,idx_lat,idx_lon] )
                            E = -1 * np.cross(I, B)
                            aValue = PEDs[idx_time,idx_lev,idx_lat,idx_lon] * np.dot(E, E)
                        else:
                            print("ERROR: UNRECOGNISED TypeOfCalculation '" + TypeOfCalculation + "'")
                            CDFroot.close()
                            return
                        
                        ResultBuckets[ kp_to_fall, alt_to_fall, maglat_to_fall, mlt_to_fall, "Vals" ].append( aValue )

    # close cdf
    print( "Process", ProcessNum, "num_of_unbinned_items =", num_of_unbinned_items )
    CDFroot.close()
    
    # save results
    #XML.WriteResultsToXML( ResultBuckets, "Ohmic", ResultFilename + ".part" + str(ProcessNum), CDF_filename )
    #Data.WriteResultsToCDF(ResultBuckets, partialResultFilename, "Joule Heating", "W/m3")
    
    # ---- save values of each bin in a txt file
    print("A", datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S") )
    for aKP in Data.KPsequence:
        for anALT in Data.ALTsequence:
            for aMagLat in Data.MLATsequence:
                for aMLT in Data.MLTsequence:    
                    if len( ResultBuckets[ aKP, anALT, aMagLat, aMLT, "Vals" ] ) > 0:
                        fname = str(aKP) + "_" + str(anALT) + "_" + str(aMagLat) + "_" + str(aMLT) + ".txt"
                        # FOR WRITING VALUES AS TEXT:
                        #    f = open( procfolder + fname, "w" )
                        #    f.write( ' '.join([str(i) for i in ResultBuckets[aKP, anALT, aMagLat, aMLT, "Vals"]]) )
                        f = open( procfolder + fname, "wb" )
                        float_array = array('d', ResultBuckets[aKP, anALT, aMagLat, aMLT, "Vals"])
                        float_array.tofile(f)
                        f.close()

    # ----
    print("Process",ProcessNum,"concluded.")
    
    