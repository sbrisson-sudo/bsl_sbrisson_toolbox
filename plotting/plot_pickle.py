#!/usr/bin/env python3

import numpy as np
from numpy import pi, sin, cos, arccos
import sys, os

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.transforms import blended_transform_factory
import obspy


from obspy.core import read
from obspy.taup import TauPyModel
from obspy.taup.utils import get_phase_names

import argparse

# parsing arguments
parser = argparse.ArgumentParser()

parser.add_argument("in_file", help="pickle traces file", type=str)
parser.add_argument("-t", dest="time_range", help="time bounds", nargs=2, type=float, default=[None,None])
parser.add_argument("-d", dest="dist_range", help="distance bounds", nargs=2, type=float, default=[None,None])
parser.add_argument("-c", dest="component", help="component to plot (R|T|Z)", type=str, default="Z")
parser.add_argument("-s", dest="scale", help="scaling of the traces", type=float, default=1.0)
parser.add_argument("-o", dest="out_file", help="output figure name", type=str, default="")
parser.add_argument("--norm", help="normalisation method", type=str, default="trace")
parser.add_argument("--title", dest="title", help="title to add to the figure (default : component)", type=str, default="")

parser.add_argument("--phases", dest="phase_list", help="list of phases arrival to compute with taup and plot", nargs="+", type=str)
parser.add_argument("--station-names", dest="plot_station_name", help="display station names", action="store_true")
parser.add_argument("--fill", help="fill the traces (blue and red)", action="store_true")

args = parser.parse_args()

# reading stream file

stream = read(args.in_file)
nb_traces = stream.count()

# get event location
trace = stream.traces[0]
lat_s = trace.stats.evla
lon_s = trace.stats.evlo
dep_s = trace.stats.evde
t0 = trace.stats.starttime

# select on component
stream = stream.select(component = args.component)
print(f"Select : {nb_traces  - stream.count()} traces removed based on component criteria.")
nb_traces = stream.count()

def select_distance(stream, d1, d2):
    """Select traces in stream based on distance to event."""
    to_remove = []
    for trace in stream.traces:
        if not(d1 <= trace.stats.distance and d2 >= trace.stats.distance):
            to_remove.append(trace)
    for trace in to_remove:
        stream.remove(trace)
    return stream


if args.dist_range[0] or args.dist_range[0]:
    d1, d2 = args.dist_range
    stream = select_distance(stream, d1, d2)
    print(f"Select : {nb_traces  - stream.count()} traces removed based on distance criteria.")

if args.time_range[0] or args.time_range[1]:
    t1,t2 = args.time_range
    start_time = trace.stats.starttime
    stream.trim(start_time + t1, start_time + t2)
else:
    t1 = 500.0


fig = plt.figure(tight_layout = True, figsize = (8,8))

stream.plot(
    type='section', 
    dist_degree=True, 
    ev_coord = (lat_s,lon_s),
    norm_method = args.norm,
    scale = args.scale,
    show=False, 
    reftime = t0,
    fig=fig,
    fillcolors = ("r","b") if args.fill else (None,None)
    )

ax = fig.axes[0]
ax.grid(True)


if args.plot_station_name:
    # ax.yaxis.set_label_position("right")
    # ax.yaxis.tick_right()
    # ax.legend(False)
    transform = blended_transform_factory(ax.transData, ax.transAxes)
    for tr in stream:
        ax.text(tr.stats.distance, 1.0,tr.stats.station, transform=transform, zorder=10, va="bottom", ha="center", fontfamily = "monospace", rotation = -45.)

# compute arrival times

if args.phase_list:
    
    obspy_phase_lists = ["ttp", "tts", "ttbasic", "tts+", "ttp+", "ttall"]
    if len(args.phase_list) == 1 and args.phase_list[0] in obspy_phase_lists:
        phase_list = get_phase_names(args.phase_list[0])
    else:
        phase_list = args.phase_list
        
    model_name="prem"
    model = TauPyModel(model=model_name)
    
    print(f"Computing phases time arrivals for model {model_name} and the following phases : {phase_list}")
    
    d1, d2 = ax.get_xlim()

    dist = np.linspace(d1, d2, 50)
    
    phase_arrival_time = dict(zip(phase_list,[[] for i in range(len(phase_list))]))
    phase_arrival_dist = dict(zip(phase_list,[[] for i in range(len(phase_list))]))

    arrivals_dist = [
        model.get_travel_times(
            source_depth_in_km=dep_s,
            distance_in_degree=d,
            phase_list = phase_list,
            ) for d in dist
    ]

    for d,arrivals in zip(dist,arrivals_dist):
        for arrival in arrivals:
            if arrival.time > 0.0: # la phase arrive bien à cette distance
                phase_arrival_time[arrival.phase.name].append(arrival.time) 
                phase_arrival_dist[arrival.phase.name].append(d) 

    colors=['#636EFA', '#EF553B', '#00CC96', '#AB63FA', '#FFA15A', '#19D3F3', '#FF6692', '#B6E880', '#FF97FF', '#FECB52']
    
    for c,phase in zip(colors,phase_list):
        ax.plot(phase_arrival_dist[phase], np.array(phase_arrival_time[phase]), c="w", lw=3)
        ax.plot(phase_arrival_dist[phase], np.array(phase_arrival_time[phase]), label=phase, c=c, lw=2)

    ax.legend(title = f"Taup model : {model_name}")
    
# title
component_names = {
    "Z" : "vertical",
    "R" : "radial",
    "T" : "transverse",
}    

prefix = args.title+" - " if args.title else ""
fig.suptitle(f"{prefix}{component_names[args.component]}")

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
