#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Thu Apr 07 2022
@author: Sylvain Brisson sylvain.brisson@ens.fr

Plot stations event (and ulvz) and associated great circle paths

"""
# importations

import argparse
import sys,os

pathdir = os.path.dirname(__file__)
pathbase = os.path.abspath(os.path.join(pathdir, os.pardir))
sys.path.append(pathbase)

from common.my_io import read_receivers
from common.barycenter_sphere import compute_mean_point_stations_event

import json

import numpy as np

import matplotlib.pyplot as plt
import matplotlib
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature

from obspy.core.event import read_events
from obspy.imaging.beachball import beach
from obspy.geodetics.base import locations2degrees, gps2dist_azimuth

import geopy
from geopy.distance import geodesic

import warnings
warnings.filterwarnings('ignore')

# constantes

R_EARTH_KM = 6371.0

# subroutines 

def plot_great_circles_obspy(event, stations, ax, color_az = True, cbar=True):
    """Plot Great circle
    - event : obspy event
    - stations : pandas dataframe with lat and lon entries
    """
    
    origin = event.preferred_origin() or event.origins[0]
    lon_s,lat_s = origin.longitude, origin.latitude

    print(color_az)

    if color_az:
        plot_great_circles_color_az(lat_s, lon_s, stations, ax, cbar=cbar)
    else:
        plot_great_circles(lat_s,lon_s,stations, ax)

def plot_great_circles(lat_s, lon_s, stations, ax):
    """Plot Great circles
    """
        
    alpha = np.exp(-len(stations.index)/500)
    # alpha = 0.2
        
    first = True
    for (lat_r,lon_r) in zip(stations["lat"], stations["lon"]):

        if first:
            ax.plot((lon_s,lon_r), (lat_s,lat_r), "k",transform=ccrs.Geodetic(), lw=1, alpha=alpha, label="great circle paths")
            first = False

        ax.plot((lon_s,lon_r), (lat_s,lat_r), "k",transform=ccrs.Geodetic(), lw=1, alpha=alpha)


def plot_great_circles_color_az(lat_s, lon_s, stations, ax, cbar=True):
    """Plot Great circles, idem with color dependant based on azimuth
    """
        
    alpha = np.exp(-len(stations.index)/500)
    
    # compute azimuth
    def get_az(row):
        az = gps2dist_azimuth(lat_s,lon_s,row.lat,row.lon)[1]
        if az > 180.:
            return az - 360.
        return az
    stations["az"] = stations.apply(get_az, axis=1)
        
    # generate colormap
    cNorm  = matplotlib.colors.Normalize(vmin=stations["az"].min(), vmax=stations["az"].max())
    scalarMap = matplotlib.cm.ScalarMappable(norm=cNorm, cmap="cool")

    first = True
    for (lat_r,lon_r,az) in zip(stations["lat"], stations["lon"], stations["az"]):

        ax.plot((lon_s,lon_r), (lat_s,lat_r),transform=ccrs.Geodetic(), lw=1, alpha=alpha, color=scalarMap.to_rgba(az))

    if cbar : ax.get_figure().colorbar(scalarMap, label="Azimuth [°]", shrink = 0.5)


def plot_event(event, ax, proj, add_title=False):
    """Plot Event BeachBall"""
    
    # get event info
    origin = event.preferred_origin() or event.origins[0]
    focmec = event.preferred_focal_mechanism() or event.focal_mechanisms[0]
    tensor = focmec.moment_tensor.tensor
    moment_list = [
        tensor.m_rr, tensor.m_tt, tensor.m_pp,
        tensor.m_rt, tensor.m_rp, tensor.m_tp]
    
    lon,lat = origin.longitude, origin.latitude
    
    # Add beachballs for two events
    xy = proj.transform_point(lon, lat, src_crs=ccrs.Geodetic())
    
    b = beach(moment_list, xy=xy, width=5e5, linewidth=1, facecolor="k")
    b.set_zorder(5)
    ax.add_collection(b)
    
    # add event info as title
    if add_title:
        title = f"Event : {event.preferred_origin().time.date.strftime('%d %b %Y')} (Mw = {event.preferred_magnitude().mag})"
        plt.title(title)
        
    
def plot_stations(stations, ax):    
    ax.plot(
        stations["lon"], stations["lat"], 
        marker="^", markerfacecolor=(0,1,0.1), markeredgecolor='k', ms=5,ls="",
        transform = ccrs.PlateCarree(), zorder=10)

def plot_hotspots(ax, write_names = False):
    '''
    Plot hotspots from Steinberger (2000)..
    '''
    
    hotspots_file = open(os.path.join(os.path.dirname(__file__), "data/hotspots.json"),"r")
    hotspots = json.load(hotspots_file)
    
    X = [xy[0] for xy in hotspots.values()]
    Y = [xy[1] for xy in hotspots.values()]
    
    if not(write_names):
    
        ax.plot(Y,X,ls='', marker="o", markerfacecolor=(0,1,0.1), markeredgecolor='k', ms=5, zorder=5, clip_on=False,transform = ccrs.PlateCarree(), label="hotspot")
        
    else:
        
        transform = ccrs.PlateCarree()._as_mpl_transform(ax)
        
        ax.plot(Y,X,ls='', marker="o", markerfacecolor=(0,1,0.1),
                markeredgecolor='k', ms=5,
                zorder=5, clip_on=False,transform = ccrs.PlateCarree(), label="hotspot")

        for name,xy in hotspots.items():
                        
            ax.annotate(name, xy[::-1], xycoords=transform)
        

def plot_plates(ax):
    '''
    Plot plate boundaries from the UTIG PLATES collection.
    '''

    for file in ["ridges.shp", "trenches.shp", "transforms.shp"]:
        
        fname = os.path.join(os.path.dirname(__file__), f"data/{file}")
        
        plateBD_feature = ShapelyFeature(Reader(fname).geometries(),ccrs.PlateCarree())
        ax.add_feature(plateBD_feature, facecolor=(0,0,0,0), edgecolor='black', lw=2.5)

        plateBD_feature = ShapelyFeature(Reader(fname).geometries(),ccrs.PlateCarree())
        ax.add_feature(plateBD_feature, facecolor=(0,0,0,0), edgecolor='red', lw=1)


# to have better defined great circles
# https://stackoverflow.com/questions/40270990/cartopy-higher-resolution-for-great-circle-distance-line
class LowerThresholdOthographic(ccrs.Orthographic):
    @property
    def threshold(self):
        return 1e3

    
if __name__ == "__main__":
    
    # setup command-line args
    parser = argparse.ArgumentParser() 
        
    parser.add_argument('--receivers', type=str,help='receivers.dat file')
    parser.add_argument('--event', type=str,help='quakeML event file')
    parser.add_argument("--ulvz", help="[lat,lon,radius] of an ulvz", nargs=3, type=float)

    parser.add_argument('--hotspots', action='store_true',help='plot main hotspots.')
    parser.add_argument('--hotspots_names', action='store_true',help='plot main hotspots names.')
    parser.add_argument('--plates', action='store_true',help='plot main plates boundaries.')
    
    parser.add_argument("-d", dest="dist_range", help="distance bounds", nargs=2, type=float, default=[0.0,180.0])
    
    parser.add_argument('-o', dest="outfile", type=str,help='out figure name')

    parser.add_argument('--no-title',dest="no_title", action='store_true')
    parser.add_argument('--no-colorbar',dest="no_cbar", action='store_true')
    parser.add_argument('--color-az',dest="color_az", action='store_true')

    args = parser.parse_args()

    # plotting
        
    fig = plt.figure(figsize=(8,8))    
    
    if args.receivers: stations = read_receivers(args.receivers)
    
    if args.event: event = read_events(args.event)[0]
    
    lat0,lon0 = compute_mean_point_stations_event(event,stations)
    
    proj = LowerThresholdOthographic(
        central_latitude=lat0,
        central_longitude=lon0
    )   
    
    ax = fig.add_subplot(1, 1, 1, projection=proj)
    ax.gridlines(linestyle=":", color="k")
    
    if args.plates: plot_plates(ax)
    
    if args.hotspots: plot_hotspots(ax, write_names=args.hotspots_names)
    
    if args.receivers and args.event: 
        
        # filter on distance
        dmin,dmax = args.dist_range
        if dmin > 0.0 or dmax < 180.0:
            ev_orig = event.preferred_origin()
            evla,evlon = ev_orig.latitude, ev_orig.longitude
            stations = stations[(locations2degrees(evla,evlon,stations["lat"],stations["lon"]) >= dmin) & (locations2degrees(evla,evlon,stations["lat"],stations["lon"]) <= dmax)]

        plot_great_circles_obspy(event,stations,ax,color_az=args.color_az,cbar=not(args.no_cbar))
        
    if args.event: plot_event(event,ax,proj)
    
    if args.receivers: plot_stations(stations, ax)
    
    if args.ulvz:

        lat,lon,radius = args.ulvz
        ulvz_coord = geopy.Point(lat,lon)

        lats_ulvz,lons_ulvz = [],[]

        for az in np.linspace(0.0, 360.0, 100):

            point = geodesic(kilometers=radius).destination(ulvz_coord, bearing=az)
            lats_ulvz.append(point.latitude)
            lons_ulvz.append(point.longitude)

        ax.plot(lons_ulvz,lats_ulvz,color="r",lw=2,transform=ccrs.Geodetic())
    
    ax.set_global()
    ax.coastlines()
    ax.add_feature(cfeature.LAND)

    # showing / saving
    if args.outfile:
        if os.path.exists(args.outfile):
            if input(f"File {args.outfile} exists, overwrite ? [y/n] ") != "y":
                print(f"No file saved.")
                exit()
                
        print(f"Saving {args.outfile}")
        plt.savefig(args.outfile, dpi=500, bbox_inches='tight')

    else:
        plt.show()