from __future__ import print_function
from load_data import *
import sys
sys.path.insert(0, "../SUPPORT_FUNCTIONS")
import os
from export_data import *
from run_MSISE00 import *
from variable_classes import *


# start load datafile
#file to read as input from input folder, and parameter could be
# "Tn", "O","O2","N2", "rho" or "" to Output selection
def call_MSISE00(file_to_read,parameter):
    variable_classes.save_name = file_to_read
    read_name = "../input/" + variable_classes.save_name + ".csv"
    if parameter=="":
        parameter="all"
    export_name=file_to_read+"_MSISE_"+parameter+".csv"
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
        run_MSISE00_func(Input, no_rows, no_columns)
        # end models run

        print('pass outputs to output dataframe')

        df_results = pd.DataFrame({'Tn_msise00_K': Outputs["Tn_msise00"],
                                   'O_msise00_cm-3': Outputs["O_msise00"],
                                   'O2_msise00_cm-3': Outputs["O2_msise00"],
                                   'N2_msise00_cm-3': Outputs["N2_msise00"],
                                   'rho_msise00_g*cm-3': Outputs["rho_msise00"]
                                   },
                                  index=input_data_file.index)

        variable_classes.MSISE00_df = pd.concat([input_data_file, df_results], axis=1,
                                                sort=False)  # sindesi tou pinaka eisagvgis dedomwnvn me ton pinaka eksodou
        # result_df.to_csv('test.csv', index=None)

        # save to csv
        print("start exporting csv")
        export_msise14(export_name,parameter)
        print("end exporting csv")

        variable_classes.Outputs = {
            "Tn_msise00": [],  # 'Tn_msise00_K'
            "O_msise00": [],  # 'O_msise00_cm-3'
            "O2_msise00": [],  # O2_msise00_cm-3'
            "N2_msise00": [],  # O2_msise00_cm-3'
            "rho_msise00": [],  # 'rho_msise00_g*cm-3'
        }
        variable_classes.MSISE00_df = pd.DataFrame()
        return export_name


if __name__ == "__main__":
    file = str(sys.argv[1])
    parameter = str(sys.argv[2])
    call_MSISE00(file, parameter)
