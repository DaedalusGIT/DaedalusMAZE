import apexpy as ap
import aacgmv2 as avm
from aacgmv2 import _aacgmv2 as avm2
import numpy as np
from datetime import datetime


    
'''
Calculates Magnetic Coordinates and Magnetic Local Time.
The values are calculated for modified apex at 90 km - which is used by TIEGCM as well.
Note: in order to convert from GeocentricLat to GeodeticLat you can use geo_lat2geod_lat.
'''
def getMagneticProperties( Time, GeodeticLat, GeodeticLon ):
    apObj = ap.Apex( Time )    
    if GeodeticLon > 180: GeodeticLon -= 360
    MagneticLatitude, MagneticLongitude = apObj.geo2apex( GeodeticLat, GeodeticLon, 90 )
    MagneticLocalTime = apObj.mlon2mlt( MagneticLongitude, Time )
    return MagneticLatitude, MagneticLongitude, MagneticLocalTime


# GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
#                  GEOCENRIC - GEODETIC CONVERSIONS (start)
# GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG

def geod_lat2geo_lat(phi):
    # calculate geocentric latitude from geodetic latitude according to WGS 84
    a = 6378137  # meter semi major axis of earth
    f = 1 / 298.257  # flattening
    b = a - f * a  # semi minor axis
    e = ((a ** 2 - b ** 2) ** (1 / 2)) / a
    phi_rad = np.deg2rad(phi)
    geo_lat = np.arctan((1 - e ** 2) * np.tan(phi_rad))
    geo_lat = np.rad2deg(geo_lat)
    return geo_lat  # in degrees


def geo_lat2geod_lat(phi):
    # calculate geocentric latitude from geodetic latitude according to WGS 84
    a = 6378137  # meter semi major axis of earth
    f = 1 / 298.257  # flattening
    b = a - f * a  # semi minor axis
    e = ((a ** 2 - b ** 2) ** (1 / 2)) / a
    phi_rad = np.deg2rad(phi)
    geod_lat = np.arctan(np.tan(phi_rad) / (1 - e ** 2))
    geod_lat = np.rad2deg(geod_lat)
    return geod_lat  # in degrees

'''
Geographic to Geodetic transformation
lat_geo_phi: the geographic latitude in degrees
lon_geo_lmd: the geographic longitude in degrees
alt_geo    : the altitude from Earth's surface
'''
def geo2geod(lat_geo_phi, lon_geo_lmd, alt_geo):
    lat_geod_phi = geo_lat2geod_lat(lat_geo_phi)  # degrees
    lon_geod_lmd = lon_geo_lmd  # degrees
    alt_geod = alt_geo  # km
    return lat_geod_phi, lon_geod_lmd, alt_geod

#Geodetic to Geographic transformation
def geod2geo(lat_geod_phi, lon_geod_lmd, alt_geod):
    lat_geo_phi = geod_lat2geo_lat(lat_geod_phi)  # degrees
    lon_geo_lmd = lon_geod_lmd  # degrees
    alt_geo = alt_geod  # km
    return lat_geod_phi, lon_geod_lmd, alt_geod

'''
EXAMPLE USAGE:
lat_geod_phi = 50.8  # degrees
lon_geod_lmd = 4.36  # degrees
alt_geod = 6365.5100  # km

lat_geo_phi,lon_geo_lmd,alt_geo=geod2geo(lat_geod_phi,lon_geod_lmd,alt_geod)
print ("GEOD lat",lat_geod_phi,"GEOD lon",lon_geod_lmd, "GEOD alt",alt_geod)
print ("--->")
print ("GEO lat",lat_geo_phi,"GEO lon",lon_geo_lmd, "GEO alt",alt_geo)

lat_geod_phi,lon_geod_lmd,alt_geod=geo2geod(lat_geo_phi,lon_geo_lmd,alt_geo)
print ("GEO lat",lat_geo_phi,"GEO lon",lon_geo_lmd, "GEO alt",alt_geo)
print ("--->")
print ("GEOD lat",lat_geod_phi,"GEOD lon",lon_geod_lmd, "GEOD alt",alt_geod)
'''
# GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
#                  GEOCENRIC - GEODETIC CONVERSIONS (finish)
# GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG










