#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 2022
@author: Sylvain Brisson sylvain.brisson@ens.fr
"""

import re
import numpy as np

import cartopy.crs as ccrs
import cartopy.feature as cfeature  
import matplotlib.pyplot as plt

plt.style.use("seaborn-paper")

import os
import sys

pathdir = os.path.dirname(__file__)
pathbase = os.path.abspath(os.path.join(pathdir, os.pardir))
sys.path.append(pathbase)

from common.my_io import read_receivers
from common.barycenter_sphere import barycenter_on_sphere, compute_mean_point_stations_event

import argparse

from obspy.core.event import read_events

if __name__ == "__main__":
    
    # setup command-line args
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--stations', type=str,help='receivers.dat file', required=True)
    parser.add_argument('--event', type=str,help='quakeML event file', required=True)
    parser.add_argument('--lonlat1', type=float,help='first profile point',nargs=2,required=True)
    parser.add_argument('--lonlat2', type=float,help='first profile point',nargs=2,required=True)
    parser.add_argument('-o', dest="outfile", type=str,help='out figure name')

    
    args = parser.parse_args()
    
    # getting station positions
    stations = read_receivers(args.stations)

    # getting event information
    event = read_events(args.event)[0]
    
    # setting/computing limits points of profil
    
    lonlat1 = args.lonlat1
    lonlat2 = args.lonlat2

    # compute mean point of stations
    lat_mean_stn,lon_mean_stn = barycenter_on_sphere(stations["lat"],stations["lon"])
    
    origin = event.preferred_origin()
        

    fig = plt.figure()
    
    lat0,lon0 = compute_mean_point_stations_event(event,stations)
    
    proj = ccrs.Robinson(
        central_longitude = lon0,
    )   
    
    ax = fig.add_subplot(1, 1, 1, projection=proj)

    ax.coastlines()
    ax.add_feature(cfeature.LAND)
    ax.gridlines(linestyle=":", color="k")
    
    lon1,lat1 = lonlat1 
    lon2,lat2 = lonlat2
    
    print([lon1,lon2,origin.longitude,lon_mean_stn],[lat1,lat2,origin.latitude,lat_mean_stn])
    
    ax.scatter(stations["lon"], stations["lat"], marker=".", color="darkgreen", s = 10, transform = ccrs.PlateCarree())
    
    ax.scatter([lon1,lon2,origin.longitude,lon_mean_stn],[lat1,lat2,origin.latitude,lat_mean_stn], transform = ccrs.PlateCarree())
    
    ax.plot((lon1,lon2), (lat1,lat2), "r",transform=ccrs.Geodetic())
    # ax.plot((lon1,lon_mean_stn), (lat1,lat_mean_stn), "b",transform=ccrs.Geodetic())
    

    
    ax.set_global()
    
    # showing / saving
    if args.outfile:
        if os.path.exists(args.outfile):
            if input(f"File {args.outfile} exists, overwrite ? [y/n] ") != "y":
                print(f"No file saved.")
                exit()
                
        print(f"Saving {args.outfile}")
        plt.savefig(args.outfile, dpi=500, bbox_inches='tight')

    else:
        plt.show()
    