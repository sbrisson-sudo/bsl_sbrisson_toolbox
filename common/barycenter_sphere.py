#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 2022
@author: Sylvain Brisson sylvain.brisson@ens.fr
"""

import numpy as np
from numpy import cos,sin
import matplotlib.pyplot as plt


from obspy.core.event import read_events


from .my_io import read_receivers


import argparse

def compute_mean_point_stations_event(event, stations):
    """Compute the mean points to set the projection
    - event : obspy event
    - stations : dataframe with a lat and lon column"""
    
    # Compute the barycenter of receivers
    lat_mean_stn, lon_mean_stn = barycenter_on_sphere(stations["lat"], stations["lon"])
    
    origin = event.preferred_origin()
    
    # Barycenter of event and stations barycenter
    lat_mean,lon_mean = barycenter_on_sphere([lat_mean_stn,origin.latitude],[lon_mean_stn, origin.longitude])
    
    return lat_mean,lon_mean
        
def barycenter_on_sphere(lats,lons):
    """Compute the mean point of the clouds points lats,lons, all inputs and outputs in degrees"""
    
    lats,lons = np.asanyarray(lats)*np.pi/180., np.asanyarray(lons)*np.pi/180.
    
     # Convert lat/lon (must be in radians) to Cartesian coordinates for each location.
    X = cos(lats) * cos(lons)
    Y = cos(lats) * sin(lons)
    Z = sin(lats)

    # Compute average x, y and z coordinates.
    x_mean = X.sum()/len(X)
    y_mean = Y.sum()/len(Y)
    z_mean = Z.sum()/len(Z)

    # Convert average x, y, z coordinate to latitude and longitude.
    lon_mean = np.arctan2(y_mean, x_mean)
    Hyp = np.sqrt(x_mean**2 + y_mean**2 + z_mean**2)
    lat_mean = np.arcsin(z_mean/Hyp)
    
    
    return lat_mean*180./np.pi,lon_mean*180./np.pi


if __name__ == "__main__":
    
    
    parser = argparse.ArgumentParser()

    parser.add_argument("-r", dest="receivers_file", help="receivers.dat file", type=str, required=True)
    parser.add_argument("--event", help="quakeML event file", type=str, required=True)

    args = parser.parse_args()
    
    
    event = read_events(args.event)[0]
    stations = read_receivers(args.receivers_file)
    
    print(compute_mean_point_stations_event(event, stations))