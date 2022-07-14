#!/usr/bin/env python3

"""
Read the traces ASCII files from 2 runs using exactly the same stations, time information
compute the difference and save it to an obspy pickle file"""

import numpy as np
import pandas as pd
from numpy import pi, cos, sin, arccos

from obspy import Stream, Trace
from obspy.core import UTCDateTime
from obspy.geodetics import locations2degrees


import os
import sys
import re
import json

import argparse

from parse_macromesh import parse_macromesh

# parsing arguments
parser = argparse.ArgumentParser()

parser.add_argument("trace_dirs", help="2 traces files directories", type=str, nargs=2)
parser.add_argument("-o", dest="out_file", help="output file prefix", type=str, required = True)
parser.add_argument("-r", dest="receivers_file", help="path to the receivers.dat file", type=str, default="receivers.dat")
parser.add_argument("-m", dest="macromesh_file", help="path to the macromesh.dat file", type=str, default="macromesh.dat")
parser.add_argument("-rinfo", dest="recepteurs_info_file", help="path to the recepteurs.info file", type=str)
parser.add_argument("-d", dest="dist_range", help="distance bounds in degrees", nargs=2, type=float, default=[0.0, 180.0])

args = parser.parse_args()

trace_re = r"^U(?P<channel>[LTR])_(?P<network>[A-Z\d]{2})_(?P<station>[_A-Z\d]{3,4})"

# from weird cSEM convention to usual letters
components = {
    "L" : "R",
    "R" : "Z",
    "T" : "T"
}

# def angularDistance(θ1, λ1, θ2, λ2):
#     θ1, θ2, λ1, λ2 = map(lambda x : x*pi/180, [θ1, θ2, λ1, λ2])
#     return arccos( sin(θ1)*sin(θ2) + cos(θ1)*cos(θ2)*cos(λ2-λ1))* 180/pi

# reading receivers.dat

data = pd.read_csv(args.receivers_file, header = 2, sep = "\s+")
stations_coord = dict(zip(data["stn"], zip(data["lat"], data["lon:"])))
print(f"{args.receivers_file} read.")

# reading macromesh.dat or recepteurs.info

# reading macromesh.dat or recepteurs.info
Rt_km = 6371.
depth_km,lat_s,lon_s,t0_s = 0.0,0.0,0.0,0.0

config = parse_macromesh(args.macromesh_file)
coord  = config["source"]["coordinates"]

depth_km = Rt_km - coord["radius"]/1000.
lat_s = 90. - coord["colatitude"]
lon_s = coord["longitude"]

t0_s = config["source"]["origin_time"]

print(f"Event at ({lat_s:.1f}°,{lon_s:.1f}°), {depth_km:.0f}km, t0 = {t0_s}s.")

# reading traces data
# computing difference
# saving differene into an obspy stream object
traces = []
starttime = UTCDateTime()

trace_dir1, trace_dir2 = args.trace_dirs

files1 = os.listdir(trace_dir1)
files2 = os.listdir(trace_dir2)

d1,d2 = args.dist_range

for i,f in enumerate(files1):
    
    re_search = re.search(trace_re, f)

    if re_search:
        
        print(f"Converting {f} ({i/len(files1)*100:.0f}%).", end="\r")
        
        chan = components[re_search.group(1)]
        netw = re_search.group(2)
        stat = re_search.group(3)
        
        data1 = np.loadtxt(os.path.join(trace_dir1, f))
        data2 = np.loadtxt(os.path.join(trace_dir2, f))
        
        sr = 1/(data1[1,0]-data1[0,0]) # sampling rate in hertz
        N = data1.shape[0]
        
        lat_r,lon_r = stations_coord[stat]
        Δ = locations2degrees(lat_r,lon_r,lat_s,lon_s)
        
        if d1 <= Δ and d2 >= Δ:
        
            # Fill header attributes
            stats = {
                'network': netw, 
                'station': stat,
                'coordinates': {'latitude':lat_r, 'longitude':lon_r},
                'evla': lat_s,
                'evlo': lon_s,
                'evde': depth_km,
                'distance': Δ,
                'component': chan, 
                'npts': N, 
                'sampling_rate': sr,
                'starttime' : starttime,
                'event_origin_time' : starttime + t0_s
                }

            traces.append(
                # Trace(data = data[:,1], header = stats).slice(UTCDateTime()+t0_s)                
                Trace(data = data2[:,1] - data1[:,1], header = stats)
            )
            
            
                        
if len(traces) == 0:
    print("Error : no traces found.")
    exit()
        
stream = Stream(traces = traces)

out_file = args.out_file + ".pickle"

if os.path.exists(out_file):
    if input(f"File {out_file} exists, overwrite ? [y/n] ") != "y": exit()
    
print(f"Writting {out_file}")
stream.write(out_file, format='PICKLE')