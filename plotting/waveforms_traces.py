#!/usr/bin/env python3

import numpy as np
from numpy import pi, cos, sin, arccos
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys, os
import pandas as pd
import re
import json

from obspy.taup import TauPyModel
from obspy.taup.utils import get_phase_names
from obspy.geodetics import locations2degrees

pathdir = os.path.dirname(__file__)
pathbase = os.path.abspath(os.path.join(pathdir, os.pardir))
sys.path.append(pathbase)

from csem_tools.parse_macromesh import parse_macromesh


mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=['636EFA', 'EF553B', '00CC96', 'AB63FA', 'FFA15A', '19D3F3', 'FF6692', 'B6E880', 'FF97FF', 'FECB52'])
# mpl.use('tkagg')

# """
# Angular distance from (lat1, lat2, lon1, lon2)
# """
# def angularDistance(θ1, λ1, θ2, λ2):
#     θ1, θ2, λ1, λ2 = map(lambda x : x*pi/180, [θ1, θ2, λ1, λ2])
#     return arccos(sin(θ1)*sin(θ2) + cos(θ1)*cos(θ2)*cos(λ1-λ2))* 180/pi


import argparse

# parsing arguments
parser = argparse.ArgumentParser()

parser.add_argument("trace_files", help="traces files to plot", nargs="+", type=str)
parser.add_argument("-o", dest="out_file", help="output file prefix", type=str)
parser.add_argument("-r", dest="receivers_file", help="path to the receivers.dat file", type=str, default="receivers.dat")
parser.add_argument("-m", dest="macromesh_file", help="path to the macromesh.dat file", type=str, default="macromesh.dat")
parser.add_argument("-rinfo", dest="recepteurs_info_file", help="path to the recepteurs.info file", type=str)
parser.add_argument("-t", dest="time_range", help="time bounds", nargs=2, type=float, default=[None,None])
parser.add_argument("--phases", dest="phase_list", help="list of phases arrival to compute with taup and plot", nargs="+", type=str)
parser.add_argument("--title", dest="title", help="title to add to the figure (default : component)", type=str, default="")
parser.add_argument("--legend", dest="labels", help="labels to add to the traces (first label for the 3 first traces etc.", type=str, nargs="+")
parser.add_argument("--diff", help="For benchmarking : plot the differences between 2 set of 3 traces, they must have the same time axis.", action="store_true")
parser.add_argument("--legfmt", help="Legacy macromesh.dat format (no ulvz flags)", action="store_true")

args = parser.parse_args()

Rt = 6371.0

# tools to parse trace file names
stationNameRe = r"_[a-zA-Z0-9]*$"    
stationDirRe = r"U[LTR]"

dirLabels = {
    "L" : "Longitudinal",
    "R" : "Vertical",
    "T" : "Transverse"
}

dirAxes = {
    "L" : 0,
    "R" : 1,
    "T" : 2
}

firstTrace = {
    "L" : True,
    "R" : True,
    "T" : True,
}

# *** STATION INFORMATION ***
# reading receivers.dat

data_sta = pd.read_csv(args.receivers_file, header = 2, sep = "\s+")
stations_coord = dict(zip(data_sta["stn"], zip(data_sta["lat"], data_sta["lon:"])))
print(f"{args.receivers_file} read.")
stations = []

# *** SYNTHETICS INFORMATION ***
# reading macromesh.dat or recepteurs.info
Rt = 6371.
dep_s,lat_s,lon_s,t0 = 0.0,0.0,0.0,0.0
fmax = 0.0 # fréquence max de la source (heaviside)


def read_recepteurs_info(recepteurs_info_file):
    global dep_s, lat_s, lon_s, t0
    with open(recepteurs_info_file, "r") as f:
        lines = f.read().splitlines()
        dep_s,lat_s,lon_s = map(float, lines[5].split())
        lat_s = 90.0 - lat_s
        dep_s = Rt - dep_s/1000
        t0 = 450.0 # pas d'info sur t0
    print("recepteurs.info read.")
    
legacy_fmt = args.legfmt
        
if args.recepteurs_info_file:
    # reading recepteurs.info
    read_recepteurs_info(args.recepteurs_info_file)
    
else:
    # reading macromesh.dat
        
    try:
        config = parse_macromesh(args.macromesh_file, legacy_fmt)
        coord_ev = config["source"]["coordinates"]
        dep_s = Rt - coord_ev["radius"]/1000.
        lat_s = 90. - coord_ev["colatitude"]
        lon_s = coord_ev["longitude"]
        t0 = config["source"]["origin_time"]
        print(f"{args.macromesh_file} read.")
    except FileNotFoundError:
        raise Exception("macromesh.dat not found, use the -m or the -rinfo options.")


# reading time bounds
if args.time_range[0] or args.time_range[1]:
    t1,t2 = args.time_range
else:
    t1 = 0
    t2 = None

# *** PLOTTING ***

fig, axes = plt.subplots(3, 1, figsize = (15,5))

# labels for each 3 traces
if args.labels and not(len(args.labels) == len(args.trace_files)//3): 
        print("Error : wrong number of labels for the number of traces.")
        args.labels = None
      
if args.diff and not(len(args.trace_files)==6):
    print("Error : wrong number of traces for doing diff (should be 6).")
    args.diff = False
    
data = dict(zip(["L","R","T"], [[],[],[]]))
        
for (i,file) in enumerate(args.trace_files):

    trace = np.loadtxt(file)
    
    filename = os.path.basename(file)

    name = re.findall(stationNameRe, filename)[0][1:]
    stations.append(name if len(name)==4 else "_"+name)

    dir = re.findall(stationDirRe, filename)[0][1]
    time = trace[:,0]-t0
    
    imin = np.argmin(np.abs(time - t1))
    imax = np.argmin(np.abs(time - t2)) if t2 else len(trace[:,0])-1

    data[dir].append(trace[imin:imax,1])
    
for dir,data_dir in data.items():
        
    for i,trace in enumerate(data_dir):
    
        if dirAxes[dir] == 0 and args.labels:
            axes[dirAxes[dir]].plot(time[imin:imax],trace, lw=0.9, label=args.labels[i])
        else:
            axes[dirAxes[dir]].plot(time[imin:imax],trace, lw=0.9)
        
    if args.diff:
        axes[dirAxes[dir]].plot(time[imin:imax],data_dir[1]-data_dir[0], lw=0.9, label="difference" if dir=="L" else None)

for ax in axes[:2]: ax.xaxis.set_ticklabels([]) 
axes[2].set_xlabel("Time (s)")
for ax,dir in zip(axes, dirLabels.keys()): 
    axes[dirAxes[dir]].set_ylabel(dirLabels[dir])
    # if t2 < 1e9: ax.set_xlim([t1, t2])
    ax.grid()

# writting station infos

stations = list(set(stations)) # élimination doublons
    
name = stations[0]
lat_r,lon_r = stations_coord[name]
dist = locations2degrees(lat_s,lon_s,lat_r,lon_r)
prefix = args.title+" - " if args.title else ""
fig.suptitle(f"{prefix}Station {name} - Δ = {dist:.1f}°")

# writting legend info

if args.labels or args.diff:
    axes[0].legend(loc=2)

# compute arrival times

ms = 10

if args.phase_list and len(stations) == 1: 
    
    obspy_phase_lists = ["ttp", "tts", "ttbasic", "tts+", "ttp+", "ttall"]
    if len(args.phase_list) == 1 and args.phase_list[0] in obspy_phase_lists:
        phase_list = get_phase_names(args.phase_list[0])
    else:
        phase_list = args.phase_list
        
    model_name="prem"
    model = TauPyModel(model=model_name)
    
    print(f"Computing phases time arrivals for model {model_name} and the following phases : {phase_list}")
    
    
    arrivals = model.get_travel_times(
            source_depth_in_km = dep_s,
            distance_in_degree = dist,
            phase_list = phase_list,
            )
    
    print(len(arrivals))
    
    colors=['#636EFA', '#EF553B', '#00CC96', '#AB63FA', '#FFA15A', '#19D3F3', '#FF6692', '#B6E880', '#FF97FF', '#FECB52']
    
    for c,arrival in zip(colors,arrivals):
        
        for ax in axes[:2]:
            ax.plot(arrival.time, 0.8*ax.get_ylim()[1], "v", ms=ms, markeredgecolor="k")
        
        # legend only on one ax
        ax = axes[2]
        ax.plot(arrival.time, 0.8*ax.get_ylim()[1], "v", ms=ms, label = arrival.phase.name, markeredgecolor="k")
        ax.legend(title = f"Taup model : {model_name}", loc=2)

#-------------
# saving figure

# showing / saving
if args.out_file:
    if os.path.exists(args.out_file):
        if input(f"File {args.out_file} exists, overwrite ? [y/n] ") != "y":
            print(f"No file saved.")
            exit()
            
    print(f"Saving {args.out_file}")
    plt.savefig(args.out_file, dpi=500, bbox_inches='tight')

else:
    plt.show()