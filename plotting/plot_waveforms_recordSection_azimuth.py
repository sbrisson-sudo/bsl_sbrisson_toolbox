#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Thu Apr 07 2022
@author: Sylvain Brisson sylvain.brisson@ens.fr

Routine to plot traces store in a pickle file in a record section over the azimuth.
The azimuth is by defualt the azimuth of the station as seen from the source and in reference to the North Pole (default) or an ULVZ (in --ulvz option used).

"""

# Importations

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
from obspy.geodetics.base import gulvz_station2dist_azimuth,locations2degrees

import argparse

# Command line options

parser = argparse.ArgumentParser()

parser.add_argument("in_file", help="pickle traces file", type=str)
parser.add_argument("-t", dest="time_range", help="time bounds", nargs=2, type=float, default=[None,None])
parser.add_argument("-d", dest="dist_range", help="distance bounds", nargs=2, type=float, default=[None,None])

parser.add_argument("-az", dest="az_range", help="azimuth bounds", nargs=2, type=float, default=[None,None])
parser.add_argument("--ulvz", dest="ulvz", help="[latitude,longitude] of the ULVZ.", nargs=2, type=float, default=[90.0,0.0])

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

# ++ reading stream file ++
stream = read(args.in_file)
nb_traces = stream.count()
print(f"{stream.count()} traces loaded")


# ++ get event information ++
trace = stream.traces[0]
lat_s = trace.stats.evla
lon_s = trace.stats.evlo
dep_s = trace.stats.evde
origin_time_event = trace.stats.event_origin_time

# ++ select on component ++
stream = stream.select(component = args.component)
print(f"Select : {nb_traces  - stream.count()} traces removed based on component criteria.")
nb_traces = stream.count()

# ++ select on distance ++
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
    
    model_name="prem"
    model = TauPyModel(model=model_name)
    
    print(f"Computing phase time arrivals for model {model_name} and the following phase : {args.phase_ref}")
    
    for tr in stream:
        
        dist_b = locations2degrees(lat_s,lon_s,*list(tr.stats.coordinates.values()))
                        
        t_arr = model.get_travel_times(
            source_depth_in_km=dep_s,
            distance_in_degree=tr.stats.distance,
            phase_list = args.phase_ref,
            )[0].time
                # print(f"Δ={tr.stats.distance:.1f}° => tS/Sdiff={t_arr:.1f}s")
        
        tr.stats.starttime -= t_arr


# ++ trim on time ++
if args.time_range[0] or args.time_range[1]:
    t1,t2 = args.time_range
    start_time = origin_time_event
    stream.trim(start_time + t1, start_time + t2)


# ++ compute azimuth  = angle(station,source,ulvz) ++
def get_angle(source, ulvz, station):
    """Compute the angle between the ulvz and the station with respect to the source.

    Args:
        source  : (lat,lon) source, in degrees
        ulvz    : (lat,lon) source, in degrees
        station : (lat,lon) source, in degrees
    """
    
    # compute distances between event - pivot - station
    Δsource_station = locations2degrees(*source, *station)
    Δsource_ulvz = locations2degrees(*source, *ulvz)
    Δulvz_station = locations2degrees(*ulvz, *station)
    
    # compute azimuth
    az = arccos((cos(Δulvz_station*pi/180.) - cos(Δsource_ulvz*pi/180.)*cos(Δsource_station*pi/180.)) / (sin(Δsource_ulvz*pi/180.)*sin(Δsource_station*pi/180.)))*180./pi
            
    return az
                
for tr in stream:
    
    tr.stats.azimuth = get_angle(
        source=[lat_s,lon_s], 
        ulvz=args.ulvz, 
        station=list(tr.stats.coordinates.values())
        )*1000 # because the record section plotter plot in km for data in m

# ++ select on azimuth ++
def select_azimuth(stream, a1, a2):
    """Select traces in stream based on azimuth."""
    to_remove = []
    for trace in stream.traces:
        if not(a1 <= trace.stats.azimuth and a2 >= trace.stats.azimuth):
            to_remove.append(trace)
    for trace in to_remove:
        stream.remove(trace)
    return stream

if args.az_range[0] or args.az_range[0]:
    az1,az2 = args.az_range
    stream = select_azimuth(stream, az1, az2)
    print(f"Select : {nb_traces  - stream.count()} traces removed based on azimuth criteria.")

# ++ initiate plot ++
fig = plt.figure(tight_layout = True, figsize = (8,8))

def switch_az_dist(st):
    """Switch azimuth and distance fields of each traces in a stream object
    Meant to be able to use the stream.plot(type='record') function along azimuth"""
    for tr in stream:
        tr.stats.distance, tr.stats.azimuth = tr.stats.azimuth,tr.stats.distance
        
switch_az_dist(stream)

# ++ plot record section ++
stream.plot(
    type='section', 
    # dist_degree=True, 
    # ev_coord = (lat_s,lon_s),
    norm_method = args.norm,
    scale = args.scale,
    show=False, 
    reftime = origin_time_event,
    fig=fig,
    fillcolors = ("r","b") if args.fill else (None,None),
    orientation = "horizontal"
    )

switch_az_dist(stream)

ax = fig.axes[0]
ax.set_ylabel("Azimuth [°]")
ax.grid(True)

# ++ identify phases by computing arrival time for the higher azimuth ++
if args.phase_list:
    
    # find max azimuth trace
    tr_max = stream[0]
    az_max = tr_max.stats.azimuth
    for tr in stream[1:]:
        if tr.stats.azimuth > az_max:
            tr_max = tr
            az_max = tr_max.stats.azimuth
            
    # get arrival times
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
            distance_in_degree = tr_max.stats.distance,
            phase_list = phase_list,
            )
        
    colors=['#636EFA', '#EF553B', '#00CC96', '#AB63FA', '#FFA15A', '#19D3F3', '#FF6692', '#B6E880', '#FF97FF', '#FECB52']
    
    for c,arrival in zip(colors,arrivals):

        ax.plot(arrival.time, 0.95*ax.get_ylim()[1], "v", ms=10, label = arrival.phase.name, markeredgecolor="k", color=c)
        ax.legend(title = f"Taup model : {model_name}", loc=2)


# ++ write station codes ++
if args.plot_station_name:
    transform = blended_transform_factory(ax.transData, ax.transAxes)
    for tr in stream:
        # print(tr.stats.station, tr.stats.distance)
        ax.text(tr.stats.azimuth/1000, 1.0, tr.stats.station, transform=transform, zorder=10, va="bottom", ha="center", fontfamily = "monospace", rotation = -45.)

# ++ title ++
component_names = {
    "Z" : "Vertical",
    "R" : "Radial",
    "T" : "Transverse",
}    
prefix = args.title+" - " if args.title else ""
fig.suptitle(f"{prefix}{component_names[args.component]}")

# ++ showing / saving ++
if args.out_file:
    if os.path.exists(args.out_file):
        if input(f"File {args.out_file} exists, overwrite ? [y/n] ") != "y":
            print(f"No file saved.")
            exit()
            
    print(f"Saving {args.out_file}")
    plt.savefig(args.out_file, dpi=500, bbox_inches='tight')

else:
    plt.show()
