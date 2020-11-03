from __future__ import print_function
import pyglow
# from SourceCode.ModulesSourceCode.IGRF.SUPPORT_FUNCTIONS.variable_classes import *
import SourceCode.ModulesSourceCode.IGRF.SUPPORT_FUNCTIONS.variable_classes as variable_classes
# pt = pyglow.Point(dn, lat, lon, alt) #for default indices
# pt = pyglow.Point(dn, lat, lon, alt, user_ind=True) #for user defined indices
min_IGRF_alt = 0
max_IGRF_alt = 30000


def run_IGRF_func(dataframe, no_rows, no_columns):
    print('IGRF started...')
    for i in range(no_rows):
        #print("run line  ", i, " of lines ", no_rows)

        pt = pyglow.Point(dataframe[i].datetime, float(dataframe[i].g_lat), float(dataframe[i].g_long),
                          float(dataframe[i].alt))
        # run iri2016...
        pt.run_igrf()

        if (dataframe[i].alt <= max_IGRF_alt) & (dataframe[i].alt >= min_IGRF_alt):
            variable_classes.Outputs["B"].append(pt.B)
            variable_classes.Outputs["Bx"].append(pt.Bx)
            variable_classes.Outputs["By"].append(pt.By)
            variable_classes.Outputs["Bz"].append(pt.Bz)
        else:
            variable_classes.Outputs["B"].append(None)
            variable_classes.Outputs["Bx"].append(None)
            variable_classes.Outputs["By"].append(None)
            variable_classes.Outputs["Bz"].append(None)  
    print('IGRF finished')
    return
