from __future__ import print_function

#import sys
#sys.path.insert(0, "../SUPPORT_FUNCTIONS")
import os
from SourceCode.ModulesSourceCode.IRI16.run_IRI16 import *
from SourceCode.ModulesSourceCode.IRI16.SUPPORT_FUNCTIONS.export_data import export_iri16
import SourceCode.ModulesSourceCode.IRI16.SUPPORT_FUNCTIONS.variable_classes as variable_classes
from SourceCode.ModulesSourceCode.IRI16.SUPPORT_FUNCTIONS.variable_classes import input_parameters_past
from SourceCode.ModulesSourceCode.IRI16.SUPPORT_FUNCTIONS.load_data import *
from shutil import copyfile
import pandas as pd

# start load datafile
#file to read as input from input folder, and parameter could be
# "Op", "O2p","Te","Ti", "Ne", "NOp" or "" to Output selection
def IRI16(file_to_read_full_path, parameter):
    print (file_to_read_full_path)
    if file_to_read_full_path is None or len(file_to_read_full_path)==0: 
        return "" 
    
    # copy input file so that can be used by IRI
    OnlyFilename =  file_to_read_full_path[ file_to_read_full_path.rfind('/')+1 : -4 ]
    #copyfile(file_to_read_full_path, "../../NAS/Data_Files/OrbitData/" + OnlyFilename + ".csv")


    #variable_classes.save_name = file_to_read_full_path
    #read_name = "SourceCode/ModulesSourceCode/input/" + variable_classes.save_name + ".csv"
    read_name=file_to_read_full_path
    if parameter=="":
        parameter="all"
    input_path=os.path.dirname(file_to_read_full_path)
    #export_name=OnlyFilename+"_IRI16_"+parameter+".csv"
    #export_name=os.path.dirname(file_to_read_full_path)+"/"+os.path.splitext(os.path.basename(file_to_read_full_path))[0]+"_IRI16_"+parameter+".csv"
    export_name="../../../home/NAS/Data_Files/ModelOutputs/IRI/"+os.path.splitext(os.path.basename(file_to_read_full_path))[0]+"_IRI16_"+parameter+".csv"
    #export_name_Checkdir="SourceCode/ModulesSourceCode/Outputs/"+export_name
    export_name_Checkdir=export_name
    
    # delete old export file
    if os.path.exists(export_name_Checkdir)==True:
        os.remove( export_name_Checkdir )
    
    if os.path.exists(export_name_Checkdir)==True:
        print("parameter for this input file already exported")
        return export_name
    else:
        input_data_file = pd.read_csv(read_name)

        no_columns = len(input_data_file.columns)
        no_rows = input_data_file.shape[0]

        # past run
        Input = [input_parameters_past() for _ in range(no_rows)]
        load_data_func_past(input_data_file, Input, no_rows)
        # end of loading data

        # start HWM run
        run_IRI16_func(Input, no_rows, no_columns)
        # end models run

        print('pass outputs to output dataframe')
        
        df_results = pd.DataFrame({'ne_iri16_cm-3': variable_classes.Outputs["ne_iri16"],
                                   'Te_iri16_K': variable_classes.Outputs["Te_iri16"],
                                   'Ti_iri16_K': variable_classes.Outputs["Ti_iri16"],
                                   'O+_iri16_cm-3': variable_classes.Outputs["Oplus_iri16"],
                                   'O2+_iri16_cm-3': variable_classes.Outputs["O2plus_iri16"],
                                   'NO+_iri16_cm-3': variable_classes.Outputs["NOplus_iri16"]
                                   },
                                  index=input_data_file.index)

        variable_classes.IRI16_df = pd.concat([input_data_file, df_results], axis=1,
                                              sort=False)  # sindesi tou pinaka eisagvgis dedomwnvn me ton pinaka eksodou
        # result_df.to_csv('test.csv', index=None)

        # save to csv
        print("start exporting csv")
        export_iri16(export_name,parameter)
        print("end exporting csv")

        variable_classes.Outputs = {
            "ne_iri16": [],  # 'ne_iri16_cm-3'
            "Te_iri16": [],  # 'Te_iri16_K'
            "Ti_iri16": [],  # 'Ti_iri16_K'
            "Oplus_iri16": [],  # 'O+_iri16_cm-3'
            "O2plus_iri16": [],  # 'O2+_iri16_cm-3'
            "NOplus_iri16": [],  # 'NO+_iri16_cm-3'
        }
        variable_classes.IRI16_df = pd.DataFrame()
        df_results= pd.DataFrame()
        return export_name_Checkdir


if __name__ == "__main__":
    file = str(sys.argv[1])
    parameter = str(sys.argv[2])
    call_IRI16(file, parameter)
