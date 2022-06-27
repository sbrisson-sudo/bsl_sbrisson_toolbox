#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Thu Apr 07 2022
@author: Sylvain Brisson sylvain.brisson@ens.fr

Plot stations event (and ulvz) and associated great circle paths with A3D background model

"""
# importations

import argparse
import pickle
import sys,os

pathdir = os.path.dirname(__file__)
pathbase = os.path.abspath(os.path.join(pathdir, os.pardir))
sys.path.append(pathbase)

from common.my_io import read_receivers
from common.setup import models_path_default
from common.barycenter_sphere import compute_mean_point_stations_event


import json

import numpy as np

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature

from obspy.core.event import read_events
from obspy.imaging.beachball import beach
from obspy.geodetics.base import locations2degrees

import pandas as pd



from Model1D import Model1D
from ModelA3d import ModelA3d
import pyspl

from UCBColorMaps import cmapSved, cmapXved

import warnings
warnings.filterwarnings('ignore')




# from FigTools import plot_plates, plot_hotspots, plot_gc_clean

# constantes

r_earth = 6371.0

# type of mean used with mean-removal turned on
# mean_type = 'harmonic'
mean_type = 'arithmetic'


# default parameters for plot models
model_name_d = "2.6S6X4"
model_file_d = os.path.join(models_path_default, "Model-2.6S6X4-ZeroMean.A3d")
parameter_d = "S" 
grid_file_d = os.path.join(models_path_default, "grid.6") 
rmean_d = False 
depth_and_level_d = [2880.0,3.0,7.265]
ref_model_file_d = os.path.join(models_path_default, "Model-2.6S6X4-ZeroMean_1D")


def get_interp(fname, dx = 1.0):
    """
    Build or load the spherical-spline interpolant for the specified
    knot-grid file.
    """
    cache = {}
    
    # try to load it
    sspl_cache_path = os.path.join(models_path_default, 'spline_interpolants.pkl')
    
    if os.path.exists(sspl_cache_path):
        with open(sspl_cache_path, 'rb') as f:
            print(f"spline_interpolants.pkl deserialized" )
            cache = pickle.load(f)
            
    if fname in cache:
        if cache[fname]['dx'] == dx:
            return cache[fname]['interp']
        
    print("Unable to deserialize spline interpolants, computing them...")
        
    # rebuild it
        
    grid = np.loadtxt(fname, skiprows=1)
    sspl = pyspl.SphericalSplines(grid[:,0], grid[:,1], grid[:,2])
    
    lons, lats = np.meshgrid(
        -180.0 + np.arange(0, 360.0 + 1e-5 * dx, dx),
        -90.0  + np.arange(0, 180.0 + 1e-5 * dx, dx))
    
    H = sspl.evaluate(lons.ravel(), lats.ravel())
    
    cache[fname] = {'dx': dx, 'interp': (lons, lats, H)}
    
    with open(sspl_cache_path, 'wb') as f:
        pickle.dump(cache, f)
        
    return lons, lats, H

def plot_model(ax,model_name=model_name_d, model_file=model_file_d, parameter=parameter_d, grid_file=grid_file_d, rmean=rmean_d, depth_and_level=depth_and_level_d, ref_model_file=ref_model_file_d):
    """Plot a map of the A3d file

    Args:
        model_file (_type_): _description_
        parameter (_type_): _description_
        grid_file (_type_): _description_
        rmean (_type_): _description_
        depth (_type_): _description_
        ref_model_file (_type_): _description_
    """
    
    # load the A3d model
    model = ModelA3d(model_file)
    model.load_from_file()
    
    # fetch the spherical spline interpolant
    lons, lats, H = get_interp(grid_file)
    
    # extract model coefficients
    coefs = model.get_parameter_by_name(parameter).get_values()
    
    # build the radial b-spline interpolant
    bspl = pyspl.CubicBSplines(model.get_bspl_knots())
    r = np.asarray([r_earth - depth_and_level[0]])
    V = bspl.evaluate(r)
    
    # sample the model (relative perturbations)
    x = H * (V * coefs).T
    
    # load the reference model, derive reference parameters
    ref = Model1D(ref_model_file)
    ref.load_from_file()
    vsv = ref.get_values(1000 * r, parameter='vsv')
    vsh = ref.get_values(1000 * r, parameter='vsh')
    vs0 = 0.001 * np.sqrt((2 * vsv ** 2 + vsh ** 2) / 3)
    xi0 = vsh ** 2 / vsv ** 2
    
    # xi is plotted with respect to isotropy
    if parameter == 'X': x = xi0 * (1.0 + x) - 1.0
    
    # rms strength at each depth
    weights = np.cos(np.radians(lats))
    weights = weights.reshape((weights.size,1))
    
    # mean removal (only valid for Vs)
    if rmean:
        assert parameter == 'S', 'Mean removal is only supported for Vs'
        vs = (1.0 + x) * vs0
        if mean_type == 'arithmetic':
            mean = (weights * vs).sum(axis = 0) / weights.sum()
        elif mean_type == 'harmonic':
            mean = weights.sum() / (weights / vs).sum(axis = 0)
        else:
            raise ValueError(f'Unrecognized mean type "{mean_type}"')
        
    
    # has a reference value been specified?
    if len(depth_and_level) == 3:
        depth, level, x1 = depth_and_level
        if parameter == 'S':
            vs = (x[:,0] + 1.0) * vs0[0]
            x[:,0] = (vs - x1) / x1
        elif parameter == 'X':
            # reminder: xi shifted to isotropy
            # back to absolute
            xi = x[:,0] + 1.0
            x[:,0] = (xi - x1) / x1
        else:
            raise ValueError(f'Unrecognized parameter "{parameter}"')
    else:
        depth, level = depth_and_level
        # only remove mean if we have not used an alternative reference
        if rmean:
            vs = (x[:,0] + 1.0) * vs0[0]
            x[:,0] = (vs - mean[0]) / mean[0]
    # to percent
    x[:,0] *= 100
    
    vmax = depth_and_level[1]
        
    # ++ plotting ++
    
        # + color map +
    # set up colormap (xi always recentered)
    if parameter == 'X':
        cmap, vmin = cmapXved(41, vmax)
        old_center = 1.0 - vmax / (vmax - vmin)
        # cmap = recenter_cmap(cmap, 0.5, old_center=old_center)
    elif parameter == 'S':
        cmap = cmapSved(41)    
    
    x_rs = x.reshape((lats.shape))
    cf = ax.contourf(lons,lats,x_rs, cmap=cmap,transform=ccrs.PlateCarree(),levels=np.linspace(x_rs.min(),x_rs.max(),100), zorder=1)
    
    print("Setting zorder")
    
    # add a colorbar
    
    cbar = plt.colorbar(cf, shrink=0.5, orientation="horizontal", format='%3.0f', label=f"{'$V_s$' if parameter=='S' else 'Xi'} - {depth:.0f}km ({model_name})", pad=0.05)
    
    
def plot_great_circles_obspy(event, stations, ax):
    """Plot Great circle
    - event : obspy event
    - stations : pandas dataframe with lat and lon entries
    """
    
    origin = event.preferred_origin() or event.origins[0]
    lon_s,lat_s = origin.longitude, origin.latitude
    
    plot_great_circles(lat_s,lon_s,stations, ax)

def plot_great_circles(lat_s, lon_s, stations, ax):
    """Plot Great circles
    """
        
    alpha = np.exp(-len(stations.index)/500)
    # alpha = 0.2
        
    first = True
    for (lat_r,lon_r) in zip(stations["lat"], stations["lon"]):
        if first:
            ax.plot((lon_s,lon_r), (lat_s,lat_r), "g",transform=ccrs.Geodetic(), lw=1, alpha=alpha, label="great circle paths")
            first = False
        ax.plot((lon_s,lon_r), (lat_s,lat_r), "g",transform=ccrs.Geodetic(), lw=1, alpha=alpha)


def plot_event(event, ax, proj):
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
    title = f"Event : {event.preferred_origin().time.date.strftime('%d %b %Y')} (Mw = {event.preferred_magnitude().mag})"
    plt.title(title)
    
    
def plot_stations(stations, ax):    
    ax.scatter(stations["lon"], stations["lat"], marker=".", color="darkgreen", s = 10, transform = ccrs.PlateCarree(), label = "stations", zorder=10)

def plot_hotspots(ax, write_names = False):
    '''Plot hotspots from Steinberger (2000)..
    '''
    
    hotspots_file = open(os.path.join(os.path.dirname(__file__), "data/hotspots.json"),"r")
    hotspots = json.load(hotspots_file)
    
    X = [xy[0] for xy in hotspots.values()]
    Y = [xy[1] for xy in hotspots.values()]
        
    # ax.scatter(X,Y,marker="o",color="g",s=1000,transform = ccrs.PlateCarree(),label="hotspots",zorder=10)
    
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
    '''Plot plate boundaries from the UTIG PLATES collection.
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
    parser.add_argument("--ulvz", help="[lat,lon] of an ulvz", nargs=2, type=float)

    parser.add_argument('--hotspots', action='store_true',help='plot main hotspots.')
    parser.add_argument('--hotspots_names', action='store_true',help='plot main hotspots names.')
    parser.add_argument('--plates', action='store_true',help='plot main plates boundaries.')
    
    parser.add_argument("--latlon", help="central latitude and longitude in degrees", nargs=2, type=float, default=[None,None])
    parser.add_argument("-d", dest="dist_range", help="distance bounds", nargs=2, type=float, default=[0.0,180.0])
    
    parser.add_argument('-o', dest="outfile", type=str,help='out figure name')

    args = parser.parse_args()
    
    # plotting
        
    fig = plt.figure(figsize=(8,8))
    
    lat0,lon0 = args.latlon
    
    
    if args.receivers: stations = read_receivers(args.receivers)
    
    if args.event: event = read_events(args.event)[0]
    
    if not(lat0) and args.receivers and args.event:
        lat0,lon0 = compute_mean_point_stations_event(event,stations)
    
    
    
    proj = LowerThresholdOthographic(
        central_latitude=lat0,
        central_longitude=lon0
    )   
    
    ax = fig.add_subplot(1, 1, 1, projection=proj)
    
    ax.gridlines(linestyle=":", color="k")

    plot_model(ax)
    
    if args.plates: plot_plates(ax)
    
    if args.hotspots: plot_hotspots(ax, write_names=args.hotspots_names)
    

        
    if args.receivers and args.event: 
        
        # filter on distance
        dmin,dmax = args.dist_range
        if dmin > 0.0 or dmax < 180.0:
            ev_orig = event.preferred_origin()
            evla,evlon = ev_orig.latitude, ev_orig.longitude
            stations = stations[(locations2degrees(evla,evlon,stations["lat"],stations["lon"]) >= dmin) & (locations2degrees(evla,evlon,stations["lat"],stations["lon"]) <= dmax)]
        plot_great_circles_obspy(event,stations,ax)
        
    if args.event: plot_event(event,ax,proj)
    
    if args.receivers: plot_stations(stations, ax)
    
    if args.ulvz:
        lat,lon = args.ulvz
        ax.scatter(lon, lat, marker="o", color="r", ec="k",lw=1.5, s = 100, transform = ccrs.PlateCarree(), label = "ULVZ", zorder=4)
    
    
    ax.set_global()
    ax.coastlines()
    ax.legend(loc='lower left')
    
    # if args.event and args.receivers:
        
    #     time = event.preferred_origin().time
    #     fig.suptitle(f"{event.origins[0].time.date.strftime('%d-%b-%Y')} - Mw={event.magnitudes[0].mag:1.1f}, $\Delta \in [{dmin:.0f}°,{dmax:.0f}°]$")

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
        
    # saving options
    optionFile = f".{os.path.basename(__file__)}_lastUsageOptions.txt"
    print(f"Saving options used in {optionFile}")
    with open(optionFile, "w") as f:
        f.write(" ".join(sys.argv[1:]))
    