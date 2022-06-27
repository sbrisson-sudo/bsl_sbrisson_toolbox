#!/usr/bin/env python3

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

parser.add_argument("--trace-dir", dest="trace_dir", help="traces files directory", type=str, default="traces")
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

t0 = config["source"]["origin_time"]



# def read_macromesh_dat(macromesh_file):
#     global depth_km, lat_s, lon_s, t0_s
#     with open(macromesh_file, "r") as io:
#         lines = io.read().splitlines()
#         i = 84
#         depth_km = Rt - float(lines[i-3])/1000
#         lat_s = 90 - float(lines[i-2]) # colatitude
#         lon_s = float(lines[i-1])
#         i = 117
#         t0_s = float(lines[i+1])
#     print("macromesh.dat read.")

# def read_recepteurs_info(recepteurs_info_file):
#     global depth_km, lat_s, lon_s, t0_s
#     with open(recepteurs_info_file, "r") as io:
#         lines = io.read().splitlines()
#         depth_km,lat_s,lon_s = map(float, lines[5].split())
#         lat_s = 90.0 - lat_s
#         depth_km = Rt - depth_km/1000
#         t0_s = 450.0 # pas d'info sur t0
#     print("recepteurs.info read.")

# def read_macromesh_json(macromesh_file):
#     global depth_km, lat_s, lon_s, t0_s, fmax
#     with open(macromesh_file, "r") as f:
#         data = json.load(f)
#         coords = data["source"]["coordinates"]
#         depth_km = Rt - coords["radius"]/1000.0
#         lat_s = 90. - coords["colatitude"]
#         lon_s = coords["longitude"]
#         t0_s =  data["source"]["origin_time"]
#     print("macromesh.json read.")
        
        
# if args.macromesh_file:
#     if os.path.splitext(args.macromesh_file)[-1] == ".dat":
#         read_macromesh_dat(args.macromesh_file)
#     else:
#         read_macromesh_json(args.macromesh_file)

# elif args.recepteurs_info_file:
#     read_recepteurs_info(args.recepteurs_info_file)
# elif os.path.exists("macromesh.json"):
#     read_macromesh_json("macromesh.json")
# elif os.path.exists("macromesh.dat"):
#     read_macromesh_dat("macromesh.dat")
# elif os.path.exists("recepteurs.info"):
#     read_recepteurs_info("recepteurs.info")
# else:
#     print("No file to read source coordinate found, exiting.")
#     exit()

print(f"Event at ({lat_s:.1f}°,{lon_s:.1f}°), {depth_km:.0f}km, t0 = {t0_s}s.")

# reading traces data inside Trace obspy type
traces = []
starttime = UTCDateTime()
files = os.listdir(args.trace_dir)

d1,d2 = args.dist_range

for i,f in enumerate(files):
    
    re_search = re.search(trace_re, f)
    if re_search:
        
        print(f"Converting {f} ({i/len(files)*100:.0f}%).", end="\r")
        
        chan = components[re_search.group(1)]
        netw = re_search.group(2)
        stat = re_search.group(3)
        
        data = np.loadtxt(os.path.join(args.trace_dir, f))
        
        sr = 1/(data[1,0]-data[0,0]) # sampling rate in hertz
        N = data.shape[0]
        
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
                Trace(data = data[:,1], header = stats)
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