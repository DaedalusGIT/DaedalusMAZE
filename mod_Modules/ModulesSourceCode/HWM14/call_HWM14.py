from __future__ import print_function
from load_data import *
import sys
sys.path.insert(0, "../SUPPORT_FUNCTIONS")
import os
from export_data import *
from run_HWM14 import *
from variable_classes import *

# start load datafile
#file to read as input from input folder, and parameter could be
# "v", "u" or "" to Output selection
def call_HWM14(file_to_read,parameter):
    variable_classes.save_name = file_to_read
    read_name = "../input/" + variable_classes.save_name + ".csv"
    if parameter=="":
        parameter="all"
    export_name=file_to_read+"_HWM_"+parameter+".csv"
    export_name_Checkdir="../Outputs/"+export_name
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
        run_HWM14_func(Input, no_rows, no_columns)
        # end models run

        print('pass outputs to output dataframe')

        df_results = pd.DataFrame({'v_hwm14_m/s': variable_classes.Outputs["v_hwm14"],
                                   'u_hwm14_m/s': variable_classes.Outputs["u_hwm14"]
                                   },
                                  index=input_data_file.index)

        variable_classes.HWM14_df = pd.concat([input_data_file, df_results], axis=1,
                              sort=False)  # sindesi tou pinaka eisagvgis dedomwnvn me ton pinaka eksodou
        # result_df.to_csv('test.csv', index=None)


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
        return export_name

if __name__ == "__main__":
    file = str(sys.argv[1])
    parameter = str(sys.argv[2])
    call_HWM14(file,parameter)
