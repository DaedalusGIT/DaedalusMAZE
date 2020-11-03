from __future__ import print_function
import pyglow
#import sys
#sys.path.insert(0, "./SUPPORT_FUNCTIONS")
import SourceCode.ModulesSourceCode.HWM14.SUPPORT_FUNCTIONS.variable_classes as variable_classes
from SourceCode.ModulesSourceCode.HWM14.SUPPORT_FUNCTIONS.variable_classes import *

# pt = pyglow.Point(dn, lat, lon, alt) #for default indices
# pt = pyglow.Point(dn, lat, lon, alt, user_ind=True) #for user defined indices
min_HWM14_alt = 0
max_HWM14_alt = 450


def run_HWM14_func(dataframe, no_rows, no_columns):
    print('models started...')
    for i in range(no_rows):
        #print("run line  ", i, " of lines ", no_rows)

        pt = pyglow.Point(dataframe[i].datetime, float(dataframe[i].g_lat), float(dataframe[i].g_long),
                          float(dataframe[i].alt))
        pt.run_hwm()
        if (dataframe[i].alt <= max_HWM14_alt) & (dataframe[i].alt >= min_HWM14_alt):
            variable_classes.Outputs["v_hwm14"].append(pt.v)  # meridional wind (northward) to output table
            variable_classes.Outputs["u_hwm14"].append(pt.u)  # zonal wind (eastward) to output table
        else:
            variable_classes.Outputs["v_hwm14"].append(None)
            variable_classes.Outputs["u_hwm14"].append(None)

    print('HWM14 finished')
    return
