from __future__ import print_function
import pyglow
# from SourceCode.ModulesSourceCode.IGRF.SUPPORT_FUNCTIONS.variable_classes import *
import SourceCode.ModulesSourceCode.MSISE00.SUPPORT_FUNCTIONS.variable_classes as variable_classes

# pt = pyglow.Point(dn, lat, lon, alt) #for default indices
# pt = pyglow.Point(dn, lat, lon, alt, user_ind=True) #for user defined indices
min_MSISE00_alt = 0
max_MSISE00_alt = 1000


def run_MSISE00_func(dataframe, no_rows, no_columns):
    print('models started...')
    for i in range(no_rows):
        #print("run line  ", i, " of lines ", no_rows)

        pt = pyglow.Point(dataframe[i].datetime, float(dataframe[i].g_lat), float(dataframe[i].g_long),
                          float(dataframe[i].alt))
        # run iri2016...
        pt.run_msis()

        if (dataframe[i].alt <= max_MSISE00_alt) & (dataframe[i].alt >= min_MSISE00_alt):
            variable_classes.Outputs["Tn_msise00"].append(pt.Tn_msis)  # msise00_output[0]=Tn
            variable_classes.Outputs["O_msise00"].append(pt.nn['O'])  # msise00_output[1]=O
            variable_classes.Outputs["O2_msise00"].append(pt.nn['O2'])  # msise00_output[2]=O2
            variable_classes.Outputs["N2_msise00"].append(pt.nn['N2'])  # msise00_output[3]=O2
            variable_classes.Outputs["rho_msise00"].append(
                pt.rho)  # msise00_output[4]=Dn (Density) total mass density rho [grams/cm^3]
        else:
            variable_classes.Outputs["Tn_msise00"].append(None)  # msise00_output[0]=Tn
            variable_classes.Outputs["O_msise00"].append(None)  # msise00_output[1]=O
            variable_classes.Outputs["O2_msise00"].append(None)  # msise00_output[2]=O2
            variable_classes.Outputs["N2_msise00"].append(None)  # msise00_output[3]=O2
            variable_classes.Outputs["rho_msise00"].append(
               None)  # msise00_output[4]=Dn (Density) total mass density rho [grams/cm^3]

    
    print('MSISE finished')
    return
