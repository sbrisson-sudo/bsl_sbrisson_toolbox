#!/usr/bin/env python

# program to plot source and stations

# importations

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

import sys,os

from obspy.geodetics import locations2degrees
from obspy.core.event import read_events

pathdir = os.path.dirname(__file__)
pathbase = os.path.abspath(os.path.join(pathdir, os.pardir))

sys.path.append(pathbase)

from common.my_io import read_receivers
from common.barycenter_sphere import compute_mean_point_stations_event

import pandas as pd

import argparse


# config

parser = argparse.ArgumentParser()

parser.add_argument("receivers_file", help="receivers.dat file", type=str)
parser.add_argument("--latlon-proj",dest="lat0lon0", help="[lat,lon] of the source", nargs=2, type=float, default=[None,None])
parser.add_argument("--names", help="annotate station codes", action="store_true", default=False)
parser.add_argument("--event", help="quakeML event file", type=str)
parser.add_argument("--ulvz", help="[lat,lon] of an ulvz", nargs=2, type=float)

parser.add_argument("--dist", help="distance bounds (require source)", nargs=2, type=float, default=[None,None])
parser.add_argument("-o", dest="out_file", help="output figure name", type=str, default="")


args = parser.parse_args()

figsize=(8,8)

# to have better defined great circles
# https://stackoverflow.com/questions/40270990/cartopy-higher-resolution-for-great-circle-distance-line
class LowerThresholdOthographic(ccrs.Orthographic):
    @property
    def threshold(self):
        return 1e3
    
lat0,lon0 = args.lat0lon0 # center of the projection

    
# loading stations data
df = read_receivers(args.receivers_file)

# loading event data
if args.event:
    event = read_events(args.event)[0]
    if not(lat0):
        lat0,lon0 = compute_mean_point_stations_event(event,df)

# initiate figure

fig = plt.figure(figsize=figsize)
    
proj = LowerThresholdOthographic(
    central_latitude  = lat0,
    central_longitude = lon0,
)   

ax = fig.add_subplot(1, 1, 1, projection=proj)

ax.coastlines()
ax.add_feature(cfeature.LAND)
ax.gridlines(linestyle=":", color="k")

# plotting source (just position)
if args.event:

    origin = event.preferred_origin()
    evla, evlon = origin.latitude, origin.longitude
    
    ax.scatter(evlon,evla, marker="*", color="r", s = 100, transform = ccrs.PlateCarree(), label="source")

# plotting ulvz position
if args.ulvz:
    ulvzla, ulvzlon = args.ulvz
    ax.scatter(ulvzlon,ulvzla, marker="o", color="orange", s = 100, transform = ccrs.PlateCarree(), label="ulvz")



# ++ filter on distance ++
if args.event and args.dist[0]:
    dmin,dmax = args.dist
    df["dist"] = [locations2degrees(evla,evlon,??,??) for ??,?? in zip(df["lat"],df["lon"])]
    df = df.loc[(df['dist'] >= dmin) & (df['dist'] <= dmax)]

    
names = list(df["stn"])
lon = df["lon"].to_numpy()
lat = df["lat"].to_numpy()

ax.scatter(lon,lat, marker="^", color="g", s = 100, transform = ccrs.PlateCarree(), label = "stations", ec="k")

transform = ccrs.PlateCarree()._as_mpl_transform(ax)

if args.names:
    for name,(lon_s,lat_s) in zip(names, zip(lon,lat)):
        ax.annotate(name, xy=(lon_s,lat_s), xycoords=transform,ha='right', va='top')
        

# if len(names) == 1:
#     ax.plot((??s,??u),(??s,??u), "gray",transform=ccrs.Geodetic(), alpha=0.4)
#     ax.plot((??s,0),(??s,90.0), "gray",transform=ccrs.Geodetic(), alpha=0.4)
#     ax.plot((0,??u),(90.,??u),  "gray",transform=ccrs.Geodetic(), alpha=0.4)
#     ax.plot((??s,lon[0]),(??s,lat[0]), "gray",transform=ccrs.Geodetic(), alpha=0.4)
#     ax.plot( (lon[0],??u),(lat[0],??u), "gray",transform=ccrs.Geodetic(), alpha=0.4)
#     ax.plot( (lon[0],0),(lat[0],90.), "gray",transform=ccrs.Geodetic(), alpha=0.4)

# ++ adding azimuth info ++
# def get_absolute_angle(source, ulvz, station):
    # """Compute the angle between the ulvz and the station with respect to the source.

    # Args:
    #     source  : (lat,lon) source, in degrees
    #     ulvz    : (lat,lon) source, in degrees
    #     station : (lat,lon) source, in degrees
    # """
    
    # ??s,??s = (source[0]+90.)*pi/180., source[1]*pi/180.
    # ??u,??u = (ulvz[0]+90.)*pi/180., ulvz[1]*pi/180.
    # ??r,??r = (station[0]+90.)*pi/180., station[1]*pi/180.
    
    
    # a = arccos(cos(??s)*cos(??u) + sin(??s)*sin(??u)*cos(??s-??u))
    # b = arccos(cos(??r)*cos(??u) + sin(??r)*sin(??u)*cos(??r-??u))
    # c = arccos(cos(??s)*cos(??r) + sin(??s)*sin(??r)*cos(??s-??r))
        
    # d = (cos(b)-cos(a)*cos(c))/(sin(a)*sin(c))
    
    # if -1 < d and d < 1:
    #     return arccos(d)*180./pi
    # return 0.0

# def get_angle(source, ulvz, station):
#     """Compute the angle between the ulvz and the station with respect to the source.

#     Args:
#         source  : (lat,lon) source, in degrees
#         ulvz    : (lat,lon) source, in degrees
#         station : (lat,lon) source, in degrees
#     """
    
#     ??s,??s = source
#     ??u,??u = ulvz
#     ??r,??r = station
    
        
#     _,az_su,_ = gps2dist_azimuth(??s, ??s, ??u,  ??u)
#     _,az_sr,_ = gps2dist_azimuth(??s, ??s, ??r,  ??r)
        
#     return az_su - az_sr

# for lon_s,lat_s in zip(lon,lat):
        
#     az = get_angle(source=[??s,??s],ulvz=[??u,??u],station=[lat_s,lon_s])
    
#     ax.annotate(f"{az:.1f}", xy=(lon_s,lat_s), xycoords=transform,ha='left', va='bottom')




# showing figure 
ax.legend()
ax.set_global()

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