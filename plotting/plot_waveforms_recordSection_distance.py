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
parser.add_argument("-c", dest="component", help="component to plot (R|T|Z)", type=str, default="T")
parser.add_argument("-s", dest="scale", help="scaling of the traces", type=float, default=1.0)
parser.add_argument("-o", dest="out_file", help="output figure name", type=str, default="")
parser.add_argument("--norm", help="normalisation method", type=str, default="trace")
parser.add_argument("--title", dest="title", help="title to add to the figure (default : component)", type=str, default="")
parser.add_argument("--phases", dest="phase_list", help="list of phases arrival to compute with taup and plot", nargs="+", type=str)
parser.add_argument("--phase-ref", dest="phase_ref", help="phase to use for time reference", type=str, nargs="+")
parser.add_argument("--station-names", dest="plot_station_name", help="display station names", action="store_true")
parser.add_argument("--fill", help="fill the traces (blue and red)", action="store_true")

args = parser.parse_args()

# reading stream file
stream = read(args.in_file)
nb_traces = stream.count()

# get event location
trace = stream.traces[0]
lat_event = trace.stats.evla
lon_event = trace.stats.evlo
depth_event = trace.stats.evde
origin_time_event = trace.stats.event_origin_time

# select on component
stream = stream.select(component = args.component)

print(f"Select : {nb_traces  - stream.count()} traces removed based on component criteria.")
nb_traces = stream.count()

# select on distance
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
    
    
# ++ align traces on arrival time ++
if args.phase_ref:
    
    model_name = "prem"
    # model_name="iasp91"

    model = TauPyModel(model=model_name)
    
    print(f"Computing phase time arrivals for model {model_name} and the following phase : {args.phase_ref}")
    
    for tr in stream:
                                
        t_arr = model.get_travel_times(
            source_depth_in_km=depth_event,
            distance_in_degree=tr.stats.distance,
            phase_list = args.phase_ref,
            )[0].time
        
        # print(f"Δ={tr.stats.distance:.1f}° => tS/Sdiff={t_arr:.1f}s")
        
        tr.stats.starttime -= t_arr

if args.time_range[0] or args.time_range[1]:
    t1,t2 = args.time_range
    start_time = origin_time_event
    
    stream.trim(start_time + t1, start_time + t2)

fig = plt.figure(tight_layout = True, figsize = (8,8))

stream.plot(
    type='section', 
    dist_degree=True, 
    ev_coord = (lat_event,lon_event),
    norm_method = args.norm,
    scale = args.scale,
    show=False, 
    reftime = origin_time_event,
    fig=fig,
    fillcolors = ("r","b") if args.fill else (None,None)
    )

ax = fig.axes[0]
ax.grid(True)

if args.plot_station_name:

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
        
    # model_name="prem"
    model_name="iasp91"
    model = TauPyModel(model=model_name)
    
    print(f"Computing phases time arrivals for model {model_name} and the following phases : {phase_list}")
    
    d1, d2 = ax.get_xlim()

    dist = np.linspace(d1, d2, 50)
    
    if args.phase_ref:
        ref_phase_arr_time = [
            model.get_travel_times(
            source_depth_in_km=depth_event,
            distance_in_degree=d,
            phase_list = args.phase_ref)[0].time for d in dist
        ]
    
    phase_arrival_time = dict(zip(phase_list,[[] for i in range(len(phase_list))]))
    phase_arrival_dist = dict(zip(phase_list,[[] for i in range(len(phase_list))]))

    arrivals_dist = [
        model.get_travel_times(
            source_depth_in_km=depth_event,
            distance_in_degree=d,
            phase_list = phase_list,
            ) for d in dist
    ]

    for d,arrivals in zip(dist,arrivals_dist):
        for arrival in arrivals:
            if arrival.time > 0.0: # la phase arrive bien à cette distance
                phase_arrival_time[arrival.phase.name].append(arrival.time) 
                phase_arrival_dist[arrival.phase.name].append(d) 
         
    if args.phase_ref:
        for phase in phase_list:
            for i,d in enumerate(phase_arrival_dist[phase]):
                phase_arrival_time[phase][i] -= ref_phase_arr_time[np.where(dist == d)[0][0]]


    colors=['#636EFA', '#EF553B', '#00CC96', '#AB63FA', '#FFA15A', '#19D3F3', '#FF6692', '#B6E880', '#FF97FF', '#FECB52']
    
    for c,phase in zip(colors,phase_list):
        ax.plot(phase_arrival_dist[phase], np.array(phase_arrival_time[phase]), c="w", lw=3)
        ax.plot(phase_arrival_dist[phase], np.array(phase_arrival_time[phase]), label=phase, c=c, lw=2)

    ax.legend(title = f"Taup model : {model_name}")
    
# title
component_names = {
    "Z" : "Vertical",
    "R" : "Longitudinal",
    "T" : "Transverse",
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
