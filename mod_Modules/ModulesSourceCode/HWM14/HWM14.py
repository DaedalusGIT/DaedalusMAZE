from __future__ import print_function

#import sys
#sys.path.insert(0, "../SUPPORT_FUNCTIONS")
import os
from SourceCode.ModulesSourceCode.HWM14.run_HWM14 import *
from SourceCode.ModulesSourceCode.HWM14.SUPPORT_FUNCTIONS.export_data import export_hwm14
import SourceCode.ModulesSourceCode.HWM14.SUPPORT_FUNCTIONS.variable_classes as variable_classes
from SourceCode.ModulesSourceCode.HWM14.SUPPORT_FUNCTIONS.variable_classes import input_parameters_past
from SourceCode.ModulesSourceCode.HWM14.SUPPORT_FUNCTIONS.load_data import *
from shutil import copyfile
import pandas as pd

# start load datafile
#file to read as input from input path given as first function argument and parameter could be
# "v", "u" or "" to Output selection
# v meridional wind (northward) to output table
# u zonal wind (eastward) to output 
def HWM14(file_to_read_full_path, parameter):
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
    #export_name=OnlyFilename+"_HWM14_"+parameter+".csv"
    #export_name=os.path.dirname(file_to_read_full_path)+"/"+os.path.splitext(os.path.basename(file_to_read_full_path))[0]+"_HWM14_"+parameter+".csv"
    export_name="../../../home/NAS/Data_Files/ModelOutputs/HWM/"+os.path.splitext(os.path.basename(file_to_read_full_path))[0]+"_HWM14_"+parameter+".csv"
    
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

        # start HWM41 run
        run_HWM14_func(Input, no_rows, no_columns)
        # end models run

        print('pass outputs to output dataframe')
        df_results = pd.DataFrame({'v_hwm14_m/s': variable_classes.Outputs["v_hwm14"],
                                   'u_hwm14_m/s': variable_classes.Outputs["u_hwm14"]
                                   },
                                  index=input_data_file.index)



        variable_classes.HWM14_df = pd.concat([input_data_file, df_results], axis=1,
                                              sort=False)  # sindesi tou pinaka eisagvgis dedomwnvn me ton pinaka eksodou


        # save to csv
        print("start exporting csv")
        export_hwm14(export_name,parameter)
        print("end exporting csv")

        variable_classes.Outputs = {
            "v_hwm14": [],  # 'v_hwm14_m/s' meridional wind (northward)
            "dynamic_range_v_hwm14": [],  # [0]->min,[1]->min_index, [2]->max, [3]->max_index
            "u_hwm14": [],  # 'u_hwm14_m/s' zonal wind (eastward)
            "dynamic_range_u_hwm14": [],  # [0]->min,[1]->min_index, [2]->max, [3]->max_index
        }
        variable_classes.HWM14_df = pd.DataFrame()
        df_results= pd.DataFrame()
        return export_name_Checkdir


if __name__ == "__main__":
    file = str(sys.argv[1])
    parameter = str(sys.argv[2])
    MSISE00(file_full_path, parameter)
