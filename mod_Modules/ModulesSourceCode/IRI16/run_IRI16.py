from __future__ import print_function
import pyglow
#import sys
#sys.path.insert(0, "../SUPPORT_FUNCTIONS")
import SourceCode.ModulesSourceCode.IRI16.SUPPORT_FUNCTIONS.variable_classes as variable_classes

# pt = pyglow.Point(dn, lat, lon, alt) #for default indices
# pt = pyglow.Point(dn, lat, lon, alt, user_ind=True) #for user defined indices
min_IRI16_alt = 60
max_IRI16_alt = 2000
min_IRI116_NOplus_O2plus = 300


def run_IRI16_func(dataframe, no_rows, no_columns):
    print('models started...')
    for i in range(no_rows):
#         print("run line  ", i, " of lines ", no_rows)

        pt = pyglow.Point(dataframe[i].datetime, float(dataframe[i].g_lat), float(dataframe[i].g_long),
                          float(dataframe[i].alt))
        # run iri2016...
        pt.run_iri()

        if (dataframe[i].alt <= max_IRI16_alt) & (dataframe[i].alt >= min_IRI16_alt):
            variable_classes.Outputs["ne_iri16"].append(pt.ne)  # iri2016_output[0]=Ne
            variable_classes.Outputs["Te_iri16"].append(pt.Te)  # iri2016_output[1]=Te
            variable_classes.Outputs["Ti_iri16"].append(pt.Ti)  # iri2016_output[2]=Ti

            # .nn list contains 'NO+', 'H+', 'O2+', 'HE+', 'O+' from iri model
            variable_classes.Outputs["Oplus_iri16"].append(pt.ni['O+'])  # iri2016_output[3]=Oplus
        else:
            variable_classes.Outputs["ne_iri16"].append(None)  # iri2016_output[0]=Ne
            variable_classes.Outputs["Te_iri16"].append(None)  # iri2016_output[1]=Te
            variable_classes.Outputs["Ti_iri16"].append(None)  # iri2016_output[2]=Ti

            # .nn list contains 'NO+', 'H+', 'O2+', 'HE+', 'O+' from iri model
            variable_classes.Outputs["Oplus_iri16"].append(None)  # iri2016_output[3]=Oplus
        if (dataframe[i].alt <= min_IRI116_NOplus_O2plus) & (dataframe[i].alt >= min_IRI16_alt):
            variable_classes.Outputs["O2plus_iri16"].append(pt.ni['O2+'])  # iri2016_output[4]=O2_plus
            variable_classes.Outputs["NOplus_iri16"].append(pt.ni['NO+'])  # iri2016_output[5]=NO_plus
        else:
            variable_classes.Outputs["O2plus_iri16"].append(None)  # iri2016_output[4]=O2_plus
            variable_classes.Outputs["NOplus_iri16"].append(None)  # iri2016_output[5]=NO_plus
    print('IRI16 finished')
    return
