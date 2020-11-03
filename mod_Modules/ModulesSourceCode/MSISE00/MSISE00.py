from __future__ import print_function

#import sys
#sys.path.insert(0, "../SUPPORT_FUNCTIONS")
import os
from SourceCode.ModulesSourceCode.MSISE00.run_MSISE00 import *
from SourceCode.ModulesSourceCode.MSISE00.SUPPORT_FUNCTIONS.export_data import export_msise00
import SourceCode.ModulesSourceCode.MSISE00.SUPPORT_FUNCTIONS.variable_classes as variable_classes
from SourceCode.ModulesSourceCode.MSISE00.SUPPORT_FUNCTIONS.variable_classes import input_parameters_past
from SourceCode.ModulesSourceCode.MSISE00.SUPPORT_FUNCTIONS.load_data import *
from shutil import copyfile
import pandas as pd

# start load datafile
#file to read as input from input path given as first function argument and parameter could be
# "Tn", "O","O2","N2", "rho" or "" to Output selection
def MSISE00(file_to_read_full_path, parameter):
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
    #export_name=OnlyFilename+"_MSISE00_"+parameter+".csv"
    #export_name=os.path.dirname(file_to_read_full_path)+"/"+os.path.splitext(os.path.basename(file_to_read_full_path))[0]+"_MSISE00_"+parameter+".csv"
    export_name="../../../home/NAS/Data_Files/ModelOutputs/MSISE/"+os.path.splitext(os.path.basename(file_to_read_full_path))[0]+"_MSISE00_"+parameter+".csv"
    #export_name_Checkdir="SourceCode/ModulesSourceCode/Outputs/"+export_name
    export_name_Checkdir=export_name

    # delete old export file (WHY???)
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

        # start IRI run
        run_MSISE00_func(Input, no_rows, no_columns)
        # end models run

        print('pass outputs to output dataframe')
        df_results = pd.DataFrame({'Tn_msise00_K': variable_classes.Outputs["Tn_msise00"],
                                   'O_msise00_cm-3': variable_classes.Outputs["O_msise00"],
                                   'O2_msise00_cm-3': variable_classes.Outputs["O2_msise00"],
                                   'N2_msise00_cm-3': variable_classes.Outputs["N2_msise00"],
                                   'rho_msise00_g*cm-3': variable_classes.Outputs["rho_msise00"]
                                   },
                                  index=input_data_file.index)


        variable_classes.MSISE00_df = pd.concat([input_data_file, df_results], axis=1,
                                              sort=False)  # sindesi tou pinaka eisagvgis dedomwnvn me ton pinaka eksodou


        # save to csv
        print("start exporting csv")
        export_msise00(export_name,parameter)
        print("end exporting csv")

        variable_classes.Outputs = {
            "Tn_msise00": [],  # 'Tn_msise00_K'
            "O_msise00": [],  # 'O_msise00_cm-3'
            "O2_msise00": [],  # O2_msise00_cm-3'
            "N2_msise00": [],  # O2_msise00_cm-3'
            "rho_msise00": [],  # 'rho_msise00_g*cm-3'
        }
        variable_classes.MSISE00_df = pd.DataFrame()
        df_results= pd.DataFrame()
        return export_name_Checkdir


if __name__ == "__main__":
    file = str(sys.argv[1])
    parameter = str(sys.argv[2])
    MSISE00(file_full_path, parameter)
