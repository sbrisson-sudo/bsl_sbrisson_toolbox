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
from obspy.geodetics.base import gps2dist_azimuth




if __name__ == "__main__":
    
    # setup command-line args
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--receivers', type=str,help='receivers.dat file', required=True)
    parser.add_argument('--event', type=str,help='quakeML event file', required=True)
    parser.add_argument('--az', type=float,help='min max azimuth', nargs=2, default=[-180.,180.])

    args = parser.parse_args()

    # getting event information
    event = read_events(args.event)[0]

    origin = event.preferred_origin()
    
    # getting station positions
    stations = read_receivers(args.receivers)

    # computing stations azimuths
    def get_az(row):
        az = gps2dist_azimuth(origin.latitude, origin.longitude, row.lat, row.lon)[1]
        if az > 180. : return az - 360.
        return az

    stations["az"] = stations.apply(
        get_az,
        axis = 1
        )

    # filter on azimuth
    azmin,azmax = args.az
    stations = stations[(stations.az >= azmin) & (stations.az <= azmax)]

    azmin = stations["az"].min()
    azmax = stations["az"].max()


    # setting/computing limits points of profil
    
    lonlat1 = (0., 0.)
    lonlat2 = (100., 0.)

    # compute mean point of stations
    lat_mean_stn,lon_mean_stn = barycenter_on_sphere(stations["lat"],stations["lon"])
    

    fig = plt.figure()
    
    lat0,lon0 = compute_mean_point_stations_event(event,stations)
    
    proj = ccrs.Robinson(
        central_longitude = lon0,
    )   
    
    ax = fig.add_subplot(1, 1, 1, projection=proj)

    ax.set_title(f"azrange = [{azmin:.1f},{azmax:.1f}]")


    # connect event
    def onclick(event):

        global lonlat1,lonlat2, line_gc, point1, point2

        xy_data = (event.xdata, event.ydata)
        # convert from data to cartesian coordinates

        lon,lat = ccrs.PlateCarree().transform_point(*xy_data, src_crs=proj)

        if event.button == 1:
            # left click
            lonlat1 = (lon,lat)
            print(f"lonlat1 = {lon:.1f} {lat:.1f}")
        if event.button == 3:
            # right click
            lonlat2 = (lon,lat)
            print(f"lonlat2 = {lon:.1f} {lat:.1f}")

        lon1,lat1 = lonlat1 
        lon2,lat2 = lonlat2 

        # update plot
        line_gc.set_data((lon1,lon2), (lat1,lat2))
        point1.set_offsets(lonlat1)
        point2.set_offsets(lonlat2)

        fig.canvas.draw()
        
    cid = fig.canvas.mpl_connect('button_press_event', onclick)

    ax.coastlines()
    ax.add_feature(cfeature.LAND)
    ax.gridlines(linestyle=":", color="k")
    
    lon1,lat1 = lonlat1 
    lon2,lat2 = lonlat2
        
    ax.scatter(stations["lon"], stations["lat"], marker=".", color="darkgreen", s = 10, transform = ccrs.PlateCarree())


    ax.scatter([origin.longitude],[origin.latitude], transform = ccrs.PlateCarree(), color="r",ec="k")
    ax.scatter([lon_mean_stn],[lat_mean_stn], transform = ccrs.PlateCarree(), color="g",ec="k")
    point1 = ax.scatter(lon1,lat1, transform = ccrs.PlateCarree(), color="b")
    point2 = ax.scatter(lon2,lat2, transform = ccrs.PlateCarree(), color="b")
    
    # ax.plot((lon1,lon_mean_stn), (lat1,lat_mean_stn), "b",transform=ccrs.Geodetic())
    
    line_gc, = ax.plot((lon1,lon2), (lat1,lat2), "r",transform=ccrs.Geodetic())

    ax.set_global()
    
    # plt.ion()
    plt.show()
    