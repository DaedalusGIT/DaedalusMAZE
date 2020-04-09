import apexpy as ap
import numpy as np
from datetime import datetime

def GeographicLat_to_GeodeticLat( GeographicLat ):
    GeodeticLat = GeographicLat
    return GeodeticLat


def get_ApexPy_Properties( Time, GeographicLat, GeographicLon, Altitude ):
    # init
    apObj = ap.Apex( Time )
    GeodeticLat = GeographicLat_to_GeodeticLat( GeographicLat )
    if GeographicLon > 180: GeographicLon -= 360
    # calc
    MagneticLatitude, MagneticLongitude = apObj.geo2apex( GeodeticLat, GeographicLon, Altitude )
    MeanLocalTime = apObj.mlon2mlt( MagneticLongitude, Time )
    # 
    return MagneticLatitude, MagneticLongitude, MeanLocalTime