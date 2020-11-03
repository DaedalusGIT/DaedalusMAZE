
from datetime import datetime, timedelta
import os.path
#import sys
#sys.path.insert(0, "../SUPPORT_FUNCTIONS")
# import variable_classes
import SourceCode.ModulesSourceCode.IGRF.SUPPORT_FUNCTIONS.variable_classes as variable_classes

def format_conv(x):  # round exponential to 2 decimal
    return "{:2,.2e}".format(float(x))


def export_hwm14(export_name,parameter):
    #create or check existense of Output directory
    # foldername = 'Outputs/'
    # try:
    #     os.mkdir(foldername)
    #     folderpath = os.path.abspath('ModulesSourceCode/'+'Outputs')
    # except:
    #     folderpath = os.path.abspath('ModulesSourceCode/'+'Outputs')
    #check for parameter selection and edit dataframe for export
    #in case of null parameter selection all parameters are exported
    if parameter=='v':
        variable_classes.HWM14_df.drop(
        ['u_hwm14_m/s'], axis=1,inplace=True)
    if parameter=='u':
        variable_classes.HWM14_df.drop(
        ['v_hwm14_m/s'], axis=1,inplace=True)


    # output_name = folderpath + '/' + export_name
    variable_classes.HWM14_df.to_csv(export_name, index=None)
    return


def export_iri16(export_name,parameter):
    #create or check existense of Output directory
    # foldername = 'Outputs/'
    # try:
    #     os.mkdir(foldername)
    #     folderpath = os.path.abspath('ModulesSourceCode/'+'Outputs')
    # except:
    #     folderpath = os.path.abspath('ModulesSourceCode/'+'Outputs')

    #check for parameter selection and edit dataframe for export
    #in case of null parameter selection all parameters are exported
    if parameter=='Op':
        variable_classes.IRI16_df.drop(
        ['NO+_iri16_cm-3', 'O2+_iri16_cm-3', 'Te_iri16_K', 'Ti_iri16_K',
        'ne_iri16_cm-3'], axis=1,inplace=True)
    if parameter=='O2p':
        variable_classes.IRI16_df.drop(
        ['NO+_iri16_cm-3', 'O+_iri16_cm-3', 'Te_iri16_K', 'Ti_iri16_K',
        'ne_iri16_cm-3'], axis=1,inplace=True)
    if parameter=='Te':
        variable_classes.IRI16_df.drop(
        ['NO+_iri16_cm-3', 'O+_iri16_cm-3', 'O2+_iri16_cm-3', 'Ti_iri16_K',
        'ne_iri16_cm-3'], axis=1,inplace=True)
    if parameter=='Ti':
        variable_classes.IRI16_df.drop(
        ['NO+_iri16_cm-3', 'O+_iri16_cm-3', 'Te_iri16_K', 'O2+_iri16_cm-3',
        'ne_iri16_cm-3'], axis=1,inplace=True)
    if parameter=='Ne':
        variable_classes.IRI16_df.drop(
        ['NO+_iri16_cm-3', 'O+_iri16_cm-3', 'Te_iri16_K', 'O2+_iri16_cm-3',
        'O2+_iri16_cm-3'], axis=1,inplace=True)

    # output_name = folderpath + '/' + export_name
    variable_classes.IRI16_df.to_csv(export_name, index=None)
    return


def export_msise00(export_name,parameter):
    # foldername = 'Outputs/'
    # try:
    #     os.mkdir(foldername)
    #     folderpath = os.path.abspath('ModulesSourceCode/'+'Outputs')
    # except:
    #     folderpath = os.path.abspath('ModulesSourceCode/'+'Outputs')

    #check for parameter selection and edit dataframe for export
    #in case of null parameter selection all parameters are exported
    if parameter=='Tn':
        variable_classes.MSISE00_df.drop(
        ['O_msise00_cm-3', 'O2_msise00_cm-3', 'N2_msise00_cm-3', 'rho_msise00_g*cm-3'],
         axis=1,inplace=True)
    if parameter=='O':
        variable_classes.MSISE00_df.drop(
        ['Tn_msise00_K', 'O2_msise00_cm-3', 'N2_msise00_cm-3', 'rho_msise00_g*cm-3'],
         axis=1,inplace=True)
    if parameter=='O2':
        variable_classes.MSISE00_df.drop(
        ['Tn_msise00_K', 'O_msise00_cm-3', 'N2_msise00_cm-3', 'rho_msise00_g*cm-3'],
         axis=1,inplace=True)
    if parameter=='N2':
        variable_classes.MSISE00_df.drop(
        ['Tn_msise00_K', 'O_msise00_cm-3', 'O2_msise00_cm-3', 'rho_msise00_g*cm-3'],
         axis=1,inplace=True)
    if parameter=='rho':
        variable_classes.MSISE00_df.drop(
        ['Tn_msise00_K', 'O_msise00_cm-3', 'O2_msise00_cm-3', 'N2_msise00_cm-3'],
         axis=1,inplace=True)

    # output_name = folderpath + '/' +export_name
    variable_classes.MSISE00_df.to_csv(export_name, index=None)
    return

def export_igrf(export_name,parameter):
    # foldername = 'Outputs/'
    # try:
    #     os.mkdir(foldername)
    #     folderpath = os.path.abspath('ModulesSourceCode/'+'Outputs')
    # except:
    #     folderpath = os.path.abspath('ModulesSourceCode/'+'Outputs')
    #check for parameter selection and edit dataframe for export
    #in case of null parameter selection all parameters are exported

    if parameter=='B':
        variable_classes.IGRF_df.drop(
        ['Bx_Tesla', 'By_Tesla', 'Bz_Tesla'],
         axis=1,inplace=True)
    if parameter=='Bx':
        variable_classes.IGRF_df.drop(
        ['B_Tesla', 'By_Tesla', 'Bz_Tesla'],
         axis=1,inplace=True)
    if parameter=='By':
        variable_classes.IGRF_df.drop(
        ['Bx_Tesla', 'B_Tesla', 'Bz_Tesla'],
         axis=1,inplace=True)
    if parameter=='Bz':
        variable_classes.IGRF_df.drop(
        ['Bx_Tesla', 'By_Tesla', 'B_Tesla'],
         axis=1,inplace=True)


    # output_name = folderpath + '/' +export_name
    variable_classes.IGRF_df.to_csv(export_name, index=None)
    return
