import numpy as np
import pandas as pd
import sympy as sy
# import math
from mpmath import *

mp.dps = 15
mp.pretty = True


class input_parameters_past:

    # /*
    #  *   NOTES ON INPUT VARIABLES:
    #  * input variables for pyglow_input
    #
    #  */

    def __init__(self, date=0, date_Julian=0, g_lat=0.0, g_long=0.0, alt=0.0,
                 datetime=0):
        self.date = date  # /* date */
        self.date_Julian = date_Julian  # /* date */
        self.alt = alt  # /* altitude in kilometers */
        self.g_lat = g_lat  # /* geodetic latitude */
        self.g_long = g_long  # /* geodetic longitude */
        self.datetime = datetime  # datetime calculated from above inputs


class input_parameters_future:

    # /*
    #  *   NOTES ON INPUT VARIABLES:
    #  * input variables for pyglow_input
    #
    #  */

    def __init__(self, date=0, date_Julian=0, g_lat=0.0, g_long=0.0, alt=0.0,
                 datetime=0, f107=0, f107a=0, f107p=0, ap_daily=0, ap=0, kp_daily=0, kp=0, dst=0, ae=0):
        self.date = date  # /* date */
        self.date_Julian = date_Julian  # /* date */
        self.alt = alt  # /* altitude in kilometers */
        self.g_lat = g_lat  # /* geodetic latitude */
        self.g_long = g_long  # /* geodetic longitude */
        self.datetime = datetime  # datetime calculated from above inputs

        self.f107 = f107  # /* daily F10.7 flux */
        self.f107a = f107a  # /* 81 day average of F10.7 flux (centered on doy) */
        self.f107p = f107p  # /* previous day F10.7 flux */
        self.ap_daily = ap_daily  # /* magnetic index(daily) */
        self.ap = ap  # /* magnetic index */
        self.kp_daily = kp_daily  # /* magnetic index(daily) */
        self.kp = kp  # /* magnetic index */
        self.dst = dst  # /* Disturbance storm time index index */
        self.ae = ae  # /* auroral electrojet index */


# Collision Frequency denoted with V
Outputs = {
    "v_hwm14": [],  # 'v_hwm14_m/s' meridional wind (northward)
    "dynamic_range_v_hwm14": [],  # [0]->min,[1]->min_index, [2]->max, [3]->max_index
    "u_hwm14": [],  # 'u_hwm14_m/s' zonal wind (eastward)
    "dynamic_range_u_hwm14": [],  # [0]->min,[1]->min_index, [2]->max, [3]->max_index
    "ne_iri16": [],  # 'ne_iri16_cm-3'
    "dynamic_range_ne_iri16": [],  # [0]->min,[1]->min_index, [2]->max, [3]->max_index
    "Te_iri16": [],  # 'Te_iri16_K'
    "dynamic_range_Te_iri16": [],  # [0]->min,[1]->min_index, [2]->max, [3]->max_index
    "Ti_iri16": [],  # 'Ti_iri16_K'
    "dynamic_range_Ti_iri16": [],  # [0]->min,[1]->min_index, [2]->max, [3]->max_index
    "Oplus_iri16": [],  # 'O+_iri16_cm-3'
    "dynamic_range_Oplus_iri16": [],  # [0]->min,[1]->min_index, [2]->max, [3]->max_index
    "O2plus_iri16": [],  # 'O2+_iri16_cm-3'
    "dynamic_range_O2plus_iri16": [],  # [0]->min,[1]->min_index, [2]->max, [3]->max_index
    "NOplus_iri16": [],  # 'NO+_iri16_cm-3'
    "dynamic_range_NOplus_iri16": [],  # [0]->min,[1]->min_index, [2]->max, [3]->max_index
    "Tn_msise00": [],  # 'Tn_msise00_K'
    "dynamic_range_Tn_msise00": [],  # [0]->min,[1]->min_index, [2]->max, [3]->max_index
    "O_msise00": [],  # 'O_msise00_cm-3'
    "dynamic_range_O_msise00": [],  # [0]->min,[1]->min_index, [2]->max, [3]->max_index
    "O2_msise00": [],  # O2_msise00_cm-3'
    "dynamic_range_O2_msise00": [],  # [0]->min,[1]->min_index, [2]->max, [3]->max_index
    "N2_msise00": [],  # O2_msise00_cm-3'
    "dynamic_range_N2_msise00": [],  # [0]->min,[1]->min_index, [2]->max, [3]->max_index
    "rho_msise00": [],  # 'rho_msise00_g*cm-3'
    "dynamic_range_rho_msise00": [],  # [0]->min,[1]->min_index, [2]->max, [3]->max_index
    "V_e": [],  # 'CollisionFrequency_e
    "dynamic_range_V_e": [],  # [0]->min,[1]->min_index, [2]->max, [3]->max_index
    "V_Oplus": [],  # 'CollisionFrequency_Oplus
    "dynamic_range_V_Oplus": [],  # [0]->min,[1]->min_index, [2]->max, [3]->max_index
    "V_NOplus": [],  # 'CollisionFrequency_NOplus
    "dynamic_range_V_NOplus": [],  # [0]->min,[1]->min_index, [2]->max, [3]->max_index
    "V_O2plus": [],  # 'CollisionFrequency_Oplus_
    "dynamic_range_V_O2plus": [],  # [0]->min,[1]->min_index, [2]->max, [3]->max_index
    "Omega_e": [],
    "Omega_Oplus": [],
    "Omega_O2plus": [],
    "Omega_NOplus": [],
    "r_e": [],
    "r_Oplus": [],
    "r_O2plus": [],
    "r_NOplus": [],
    "sigmap": [],
    "dynamic_range_sigmap": [],  # [0]->min,[1]->min_index, [2]->max, [3]->max_index
    "sigmap_error": [],
    "sigmap_error_percentage": [],
    "dynamic_range_sigmap_error": [],  # [0]->min,[1]->min_index, [2]->max, [3]->max_index
    "B": []  # magnetic field
}


class error_parameters:
    def __init__(self, irine_error=0.000, iriTe_error=0.000, iriTi_error=0.0, iriOplus_error=0.000,
                 iriO2plus_error=0.000,
                 iriNOplus_error=0.000, msiseO_error=0.000, msiseO2_error=0.000, msiseN2_error=0.000, B_error=0.000):
        self.irine_error = irine_error  #
        self.iriTe_error = iriTe_error  #
        self.iriTi_error = iriTi_error  #
        self.iriOplus_error = iriOplus_error  #
        self.iriO2plus_error = iriO2plus_error  #
        self.iriNOplus_error = iriNOplus_error  #

        self.msiseO_error = msiseO_error  #
        self.msiseO2_error = msiseO2_error  #
        self.msiseN2_error = msiseN2_error  #

        self.B_error = B_error  # magnetic field


HWM14_df = pd.DataFrame()
IRI16_df = pd.DataFrame()
MSISE00_df = pd.DataFrame()
CollisionFrequency_df = pd.DataFrame()
omega_df = pd.DataFrame()
r_df = pd.DataFrame()
sigma_p_df = pd.DataFrame()
IGRF_df = pd.DataFrame()

# constants for Collision Frequency Calculation
c1 = 2.33e-11
c2 = 1.21e-4
c3 = 1.82e-10
c4 = 3.6e-2
c5 = 8.9e-11
c6 = 5.7e-4
c7 = 4.27e-10
c8 = 2.44e-10
c9 = 2.44e-10
c10 = 6.64e-10
c11 = 3.67e-11
c12 = 6.82e-10
c13 = 2.59e-11
c14 = 2.31e-10
c15 = 4.13e-10
c16 = 0.75
c17 = 0.064
fb = 1.7  # Burnside Factor, the factor that connect theoretical and practical estimations of Oplus [from 1.7-0.3 to 1.7+0.7]

# apo edv prwtoi oroi
ArO = 15.9994  # atomic mass O=6
ArN = 14.0067  # atomic mass N=14
NAvog = 6.02214076e23

E = 1  # Electric field
E_error = 0.01 * E
q_e = 1.60217657e-17  # Coulomb
q_e_error = 0
m_NOplus = (ArO + ArN) / NAvog * 1000
m_NOplus_error = 0.01 * m_NOplus
m_O2plus = 2 * ArO / (NAvog * 1000)
m_O2plus_error = 0.01 * m_O2plus
m_Oplus = 1 * ArO / (NAvog * 1000)  # kg
m_Oplus_error = 0.01 * m_Oplus  # kg
m_e = 9.11e-37  # kg
m_e_error = 0.01 * m_e
CubicCm2CubicM = 1e-6

save_name=""
