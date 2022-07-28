#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Thu Apr 07 2022
@author: Sylvain Brisson sylvain.brisson@ens.fr

Routine to plot traces store in a pickle file in a record section over the azimuth.
The azimuth is by defualt the azimuth of the station as seen from the source and in reference to the North Pole (default) or an ULVZ (in --ulvz option used).

"""

# Importations

import json 
from numpy import pi, sin, cos, arccos

import  os

import matplotlib.pyplot as plt
from matplotlib.transforms import blended_transform_factory

plt.rcParams.update({'font.size': 20})


from obspy.core import read
from obspy.taup import TauPyModel
from obspy.taup.utils import get_phase_names
from obspy.geodetics.base import locations2degrees, gps2dist_azimuth

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

parser.add_argument("--norm", help="normalisation method", type=str, default="stream")
parser.add_argument("--title", dest="title", help="title to add to the figure (default : component)", type=str, default="")

parser.add_argument("--phases", dest="phase_list", help="list of phases arrival to compute with taup and plot", nargs="+", type=str)
parser.add_argument("--phase-ref", dest="phase_ref", help="phase to use for time reference", type=str, nargs="+")

parser.add_argument("--station-names", dest="plot_station_name", help="display station names", action="store_true")
parser.add_argument("--event-metadata", dest="event_metadata", help="add text with event metadata", action="store_true")
parser.add_argument("--fill", help="fill the traces (blue and red)", action="store_true")
parser.add_argument("--colorbar", help="add azimuth colorbar", action="store_true")

parser.add_argument("--stn-list", dest="station_list", help="station list file", type=str)

parser.add_argument("--ulvz-metadata", dest="ulvz_metadata", help="ULVZ radius heigh and dlnvs", nargs=3, type=float)


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

# ++ select base on station name if asked ++
if args.station_list:
    stations_to_plot = open(args.station_list, "r").readlines()
    stations_to_plot = [s.strip() for s in stations_to_plot]
    to_remove = []
    for tr in stream:
        if tr.stats.station not in stations_to_plot:
            to_remove.append(tr)
    for trace in to_remove:
        stream.remove(trace)
    print(f"Select : {nb_traces  - stream.count()} traces removed based on station code")
    nb_traces = stream.count()


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

# get min_max distances
dmin,dmax = 180., 0. 
for tr in stream:
    if tr.stats.distance < dmin: dmin = tr.stats.distance
    if tr.stats.distance > dmax: dmax = tr.stats.distance

# ++ align traces on arrival time ++
if args.phase_ref:
    
    model_name="prem"
    model = TauPyModel(model=model_name)
    
    for tr in stream:
                                
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
        
    station_latlon = tr.stats.coordinates["latitude"],tr.stats.coordinates["longitude"]
    
    # tr.stats.azimuth = - get_angle(
    #     source=[lat_s,lon_s], 
    #     ulvz=args.ulvz, 
    #     station=station_latlon
    #     )*1000 # because the record section plotter plot in km for data in m
    
    # if tr.stats.coordinates["longitude"] < 0:
    #     tr.stats.azimuth *= -1

    # print(lat_s,lon_s,tr.stats.coordinates["latitude"],tr.stats.coordinates["longitude"])


    tr.stats.azimuth = gps2dist_azimuth(lat_s,lon_s,tr.stats.coordinates["latitude"],tr.stats.coordinates["longitude"])[1]


    if tr.stats.azimuth > 180.:
        tr.stats.azimuth -= 360.

    
    tr.stats.azimuth *= 1000. # routine de plot en km
    
    # print(f"{station_latlon} : az = {tr.stats.azimuth}")
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

if args.az_range[0] or args.az_range[1]:
    az1,az2 = args.az_range
    assert(az1 < az2)
    stream = select_azimuth(stream, az1*1000, az2*1000)
    print(f"Select : {nb_traces  - stream.count()} traces removed based on azimuth criteria.")

# generate colormap
azimuths = [tr.stats.azimuth/1000. for tr in stream]
az_min = min(azimuths)
az_max = max(azimuths)

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

from custom_obspy_imaging_waveform import WaveformPlotting
import matplotlib
import numpy as np

cNorm  = matplotlib.colors.Normalize(vmin=az_min, vmax=az_max)
scalarMap = matplotlib.cm.ScalarMappable(norm=cNorm, cmap="cool")

def get_max_abs_stream(stream):
    """Return the absolute maximum value of all the traces data"""
    norm_factor = 0.
    for tr in stream:
        tr_max = np.abs(tr.data.max())
        if tr_max > norm_factor:
            norm_factor = tr_max
    return norm_factor

norm_factor1 = get_max_abs_stream(stream)

t1,t2 = args.time_range
t3 = t1 + 0.1*(t2-t1)
start_time = origin_time_event
stream.trim(start_time + t1, start_time + t3)

norm_factor2 = get_max_abs_stream(stream)

scale2 = args.scale * norm_factor2 / norm_factor1

waveform = WaveformPlotting(
    stream=stream, 
    type='section', 
    norm_method = args.norm,
    scale = scale2,
    show=False, 
    reftime = origin_time_event,
    fig=fig,
    fillcolors = ("r","b") if args.fill else (None,None),
    orientation = "horizontal",
    color="distance",
    cmap=scalarMap.to_rgba,
    linewidth = 2
    )

waveform.plot_waveform()

if args.colorbar: plt.colorbar(scalarMap, label="Azimuth [°]", shrink=0.5)

switch_az_dist(stream)


# add event info

ax = fig.axes[0]

if args.event_metadata:

    from matplotlib.offsetbox import AnchoredText

    time_str = origin_time_event.date.strftime('%Y/%m/%d')

    event_metadata_str = f"Event metadata :\n- latitude {lat_s:.1f}°\n- longitude {lon_s:.1f}°\n- depth {dep_s:.1f}km\n- time {time_str}"

    props = dict(boxstyle='round', facecolor="white", edgecolor="#D7D7D7", alpha=0.75) # these are the text box parameters
    anchored_text = AnchoredText(event_metadata_str, loc=4, frameon=False, prop=dict(bbox=props))
    ax.add_artist(anchored_text)


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
    
    # compute ref phase tim arrival
    t_ref = 0.0
    if args.phase_ref:
        
        arrival_ref = model.get_travel_times(
            source_depth_in_km = dep_s,
            distance_in_degree = tr_max.stats.distance,
            phase_list = args.phase_ref,
            )[0]
        
        t_ref = arrival_ref.time
    
    for c,arrival in zip(colors,arrivals):
        
        # ax.plot(arrival.time-t_ref, 0.95*ax.get_ylim()[1], "v", ms=10, label = arrival.phase.name, markeredgecolor="k", color=c)
        ax.axvline(x=arrival.time-t_ref, color=c, label=arrival.phase.name)
        # ax.legend(title = f"Taup model : {model_name}", loc="upper right")


# ++ write station codes ++
if args.plot_station_name:
    transform = blended_transform_factory(ax.transData, ax.transAxes)
    for tr in stream:
        # print(tr.stats.station, tr.stats.distance)
        ax.text(tr.stats.distance/1000, 1.0, tr.stats.station, transform=transform, zorder=10, va="bottom", ha="center", fontfamily = "monospace", rotation = -45.)

# ++ title ++
component_names = {
    "Z" : "Vertical",
    "R" : "Radial",
    "T" : "Transverse",
}    
prefix = args.title+" - " if args.title else ""
title = f"{prefix}{component_names[args.component]} component"

# add distance information 

title += f" - $\Delta\in$[{dmin:.0f}°,{dmax:.0f}°]"

# fig.suptitle(title)


# outputting station list
with open("plotted_stations.txt", "w") as f:
    for tr in stream:
        f.write(tr.stats["station"] + "\n")



# ++ showing / saving ++
if args.out_file:
    if os.path.exists(args.out_file):
        if input(f"File {args.out_file} exists, overwrite ? [y/n] ") != "y":
            print(f"No file saved.")
            exit()
            
    print(f"Saving {args.out_file}")
    plt.savefig(args.out_file, dpi=500, bbox_inches='tight')
    # # also saving configuration
    # config_file = f"{args.out_file}.config"
    # print(f"Configuration options saved in {config_file}")
    # config_json = json.dumps(vars(args), indent = 3)
    # with open(config_file, "w") as out:
    #     out.write(config_json)


else:
    plt.show()
