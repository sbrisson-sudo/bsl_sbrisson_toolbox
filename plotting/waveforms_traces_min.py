#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import pandas as pd
import re


mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=['636EFA', 'EF553B', '00CC96', 'AB63FA', 'FFA15A', '19D3F3', 'FF6692', 'B6E880', 'FF97FF', 'FECB52'])

import argparse

# parsing arguments
parser = argparse.ArgumentParser()

parser.add_argument("trace_files", help="traces files to plot", nargs="+", type=str)
parser.add_argument("-o", dest="out_file", help="output file prefix", type=str)
parser.add_argument("-t", dest="time_range", help="time bounds", nargs=2, type=float, default=[None,None])
parser.add_argument("--title", dest="title", help="title to add to the figure (default : component)", type=str, default="")
parser.add_argument("--labels", help="labels to add to the traces (first label for the 3 first traces etc.", type=str, nargs="+")

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
        
for (i,file) in enumerate(args.trace_files):

    trace = np.loadtxt(file)
    
    filename = os.path.basename(file)

    name = re.findall(stationNameRe, filename)[0][1:]

    dir = re.findall(stationDirRe, filename)[0][1]
    time = trace[:,0]
    data = trace[:,1]
    
    imin = np.argmin(np.abs(time - t1))
    imax = np.argmin(np.abs(time - t2)) if t2 else len(trace[:,0])-1

    if dirAxes[dir] == 0 and args.labels:
        axes[dirAxes[dir]].plot(time[imin:imax],data[imin:imax], lw=0.9, label=args.labels[i//3])
    else:
        axes[dirAxes[dir]].plot(time[imin:imax],data[imin:imax], lw=0.9)



for ax in axes[:2]: ax.xaxis.set_ticklabels([]) 
axes[2].set_xlabel("Time (s)")
for ax,dir in zip(axes, dirLabels.keys()): 
    axes[dirAxes[dir]].set_ylabel(dirLabels[dir])
    ax.grid()

# writting legend info

if args.labels:
    axes[0].legend(loc=2)

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