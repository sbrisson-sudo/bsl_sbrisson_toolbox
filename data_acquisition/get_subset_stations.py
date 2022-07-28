#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Thu Apr 07 2022
@author: Sylvain Brisson sylvain.brisson@ens.fr

Read a station list and return a new file containing the stations
verifying distance and/or azimuth criteria with respect to the source

"""

import argparse 

import sys,os

pathdir = os.path.dirname(__file__)
pathbase = os.path.abspath(os.path.join(pathdir, os.pardir))
sys.path.append(pathbase)

from common.my_io import read_receivers, write_receivers

from obspy import read_events
from obspy.geodetics.base import locations2degrees,gps2dist_azimuth

parser = argparse.ArgumentParser()

parser.add_argument("in_file", help="stations list file", type=str)
parser.add_argument("-e",dest="ev_file", help="event file", type=str)
parser.add_argument("-d",dest="dist_range", nargs=2, type=float)
parser.add_argument("-az",dest="az_range", nargs=2, type=float)
parser.add_argument("-o",dest="out_file", type=str, default="receivers.dat.subset")

args = parser.parse_args()

# read input data
stations = read_receivers(args.in_file)
ev = read_events(args.ev_file)[0]

ev_orig = ev.preferred_origin()
evla,evlon = ev_orig.latitude, ev_orig.longitude


# filter on distance
if args.dist_range:

    dmin,dmax = args.dist_range

    get_dist = lambda row : locations2degrees(evla,evlon,row.lat,row.lon)

    stations["dist"] = stations.apply(get_dist, axis=1)
    stations = stations[(stations["dist"] >= dmin) & (stations["dist"] <= dmax)]

    if stations.empty:
        print("Warning : no stations respecting the distance criteria")

# filter on azimuth
if args.az_range:

    azmin,azmax = args.az_range

    def get_az(row):
        az = gps2dist_azimuth(evla,evlon,row.lat,row.lon)[1]
        return az if az < 180 else az - 360

    stations["az"] = stations.apply(get_az, axis=1)

    stations = stations[(stations["az"] >= azmin) & (stations["az"] <= azmax)]

    if stations.empty:
        print("Warning : no stations respecting the azimuth criteria")

# save remaining stations




print(f"Saving remaining stations in {args.out_file}")
write_receivers(stations, args.out_file)





