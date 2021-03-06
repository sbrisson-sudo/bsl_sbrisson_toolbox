#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 2022
@author: Sylvain Brisson sylvain.brisson@ens.fr
"""

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature  

# plt.style.use("seaborn-paper")
plt.rcParams.update({'font.size': 13})

import os
import sys


pathdir = os.path.dirname(__file__)
pathbase = os.path.abspath(os.path.join(pathdir, os.pardir))
sys.path.append(pathbase)

from common.setup import models_path_default
from common.my_io import read_receivers, to_obspy_inventory
from common.barycenter_sphere import barycenter_on_sphere

from Model1D import Model1D
from ModelA3d import ModelA3d
import pyspl
from UCBColorMaps import cmapSved, cmapXved

from Sphere import delaz, gc_minor

import argparse

from obspy.core.event import read_events
from obspy.taup.ray_paths import get_ray_paths
from obspy.core.event import Catalog

from obspy.geodetics.base import gps2dist_azimuth, locations2degrees

import pickle


rEarth = 6371.0

def plot_model(lonlat1, lonlat2, modelConfig, interpConfig, vmax=3.0, param='S'):
    
    # distance, azimuth
    delta, az = delaz(lonlat1, lonlat2)
    
    # step size
    dx = delta / (interpConfig["nx"] - 1)
    
    # compute great-circle path
    lons, lats = gc_minor(lonlat1, lonlat2, dx)
    
    # setup depth sampling
    dr = (interpConfig["dmax"] - interpConfig["dmin"]) / (interpConfig["nr"] - 1)
    r = rEarth - interpConfig["dmax"] + dr * np.arange(interpConfig["nr"])
    
     # load the A3d model
    model = ModelA3d(modelConfig['modelfile'])
    model.load_from_file()
    coefs = model.get_parameter_by_name(param).get_values()
    
    # ++ interpolation phase 1 (model)  ++

    # sample the model (polar coords in great-circle plane)
    # load the grid
    grid = np.loadtxt(modelConfig['gridfiles'][param], skiprows=1)
    # compute the sspl interpolant
    sspl = pyspl.SphericalSplines(grid[:,0], grid[:,1], grid[:,2])
    H = sspl.evaluate(lons.ravel(), lats.ravel())
    # compute the bspl interpolant
    bspl = pyspl.CubicBSplines(model.get_bspl_knots())
    V = bspl.evaluate(r)
    # sample
    x = (H * (V * coefs).T).T
    
    # to percent
    x *= 100
    
    # define circular coord grid in the great-circle plane (polar coords) over
    t = dx * np.arange(interpConfig["nx"])
    
    # + color map +
    # set up colormap (xi always recentered)
    if param == 'X':
        cmap, vmin = cmapXved(41, vmax)
        old_center = 1.0 - vmax / (vmax - vmin)
        # cmap = recenter_cmap(cmap, 0.5, old_center=old_center)
    elif param == 'S':
        cmap = cmapSved(41)    
        
    # ++ plotting ++
    
    # setup plot (figure and axes)
    fig, ax = plt.subplots(figsize=(10,10),subplot_kw={'projection': 'polar'})
    
    ax.set_theta_direction(-1)
    ax.set_theta_offset(np.pi / 2.0 + delta*np.pi/360)

    nlev = 200
    pp = plt.contourf(t*np.pi/180.0, r, x, cmap=cmap, levels=nlev, vmin=-vmax, vmax=vmax)
    
    plt.setp(ax, rorigin=0, rmin=r.min(), rmax=r.max())
    ax.set_xlim((t.min()*np.pi/180., t.max()*np.pi/180.))

    ax.tick_params(labelleft=False, labelright=True,
                    labeltop=True, labelbottom=False)
    ax.yaxis.set_major_formatter('{x:.0f}km')
    
    # + add lonlat +
        
    lon,lat = lonlat1
    # ax.annotate(f"({lat:.0f}??N,{lon:.0f}??W)", xy=(t.min()*np.pi/180.,r.max()), annotation_clip=False, horizontalalignment='right', verticalalignment='bottom',) 
    # lon,lat = lonlat2
    # ax.annotate(f"({lat:.0f}??N,{lon:.0f}??W)", xy=(t.max()*np.pi/180.,r.max()), annotation_clip=False, horizontalalignment='left', verticalalignment='bottom',) 
        
    ax.grid(False)
    
    # + add colorbar +
    
    # m = plt.cm.ScalarMappable(cmap=cmap)
    # m.set_array(x)
    # m.set_clim(-vmax, vmax)
    # ticks = np.arange(-int(vmax-0.1), int(vmax))
    # plt.colorbar(m, orientation="horizontal", shrink=0.3, pad=-0.33, aspect=10, boundaries=np.linspace(-vmax, vmax, nlev+1), ticks=ticks, label=f"$V_S$ anomaly (%) - model {modelConfig['name']}")
    
    
    return ax
    

def plot_paths(event, stations, ax, lonlat1, phases=["S", "Sdiff"]):
    
    # model = TauPyModel(model="prem")
    
    origin = event.preferred_origin() or event.origins[0]
    evlon,evlat = origin.longitude, origin.latitude
    evdepth = origin.depth
    evDist,_ =delaz((evlon,evlat),lonlat1)
    
    # filter over great circle distance
        
    ax.scatter(evDist*np.pi/180, rEarth-evdepth/1000, zorder=10, clip_on=False, marker="*", color="r", ec="k", s=200, lw=1.5)

    stations["dist"] = stations.apply(
        lambda row: delaz((row["lon"],row["lat"]),lonlat1)[0],
        axis=1
    )
    
    stations["depth"] = np.ones(len(stations.index))*rEarth
    
    ax.scatter(stations["dist"]*np.pi/180, stations["depth"], color="g", marker=".", clip_on=False, zorder=10, ec="k")
    
    # plot paths
        
    paths = get_paths(stations, event, phases, lonlat1)
    
    alpha = np.exp(-len(stations.index)/200)
    
    colors=['#636EFA', '#EF553B', '#00CC96', '#AB63FA', '#FFA15A', '#19D3F3', '#FF6692', '#B6E880', '#FF97FF', '#FECB52']
    
    for c,(phase_name,phase_paths) in zip(colors, paths.items()):
        plt.plot([], [], color=c, label=phase_name)
        for path in phase_paths:
            plt.plot(path["dist"], path["radius"], color=c, alpha=alpha)
    

def get_paths(stations, event, phases, lonlat1):
    
    """Return the path in dist (radians) and radius (km)"""
    
    path_serial_file = "paths.pickle"
    
    if os.path.exists(path_serial_file):
        print(f">> Reading paths from {path_serial_file}")
        with open(path_serial_file, 'rb') as f:
            print(f"paths deserialized" )
            paths_cache = pickle.load(f)
            
        # check that the cache contains all the asked phases
        if all([p in paths_cache.keys() for p in phases]): return paths_cache
        print(">> One on more phases missing in the cache, computation needed")
    
    print(">> Computing ray paths")
    paths = get_ray_paths(inventory=to_obspy_inventory(stations), catalog=Catalog([event]), phase_list=phases, taup_model="prem", coordinate_system="RTP")
    
    paths_cache = dict(zip(phases, [[] for _ in range(len(phases))]))
    
    for gcircle, phase_name,_,_,_,_,_ in paths:
        
        r = gcircle[0,:]*rEarth
        lat_r = 90 - gcircle[1,:]*180/np.pi 
        lon_r = gcircle[2,:]*180/np.pi 
    
        dist = np.array([delaz(lonlat1, (lon,lat))[0] for (lat,lon) in zip(lat_r,lon_r)])*np.pi/180
        
        paths_cache[phase_name].append({"dist":dist, "radius":r})    
        
    print(f">> Writting paths in {path_serial_file}")
    with open(path_serial_file, 'wb') as f:
        pickle.dump(paths_cache, f)
    
    return paths_cache
    
    
def plot_ulvz_simple(lonlat_ulvz, lonlat1, ax):
    
    dist = delaz(lonlat_ulvz, lonlat1)[0]   
    ax.scatter(dist*np.pi/180., 3480., marker="o", color="r", ec='k', s=200, lw=1.5, zorder=10)

if __name__ == "__main__":
    
    # parameter to plot (A3d descriptor name)
    param = 'S'
    
    # model configuration
    modelConfig = {
        # simply for naming of output files
        'name': '2.6S6X4',
        # A3d file
        'modelfile': os.path.join(models_path_default,'Model-2.6S6X4-ZeroMean.A3d'),
        # 1D reference model (needed for plotting xi sections)
        'refmodelfile': os.path.join(models_path_default,'Model-2.6S6X4-ZeroMean_1D'),
        # spherical spline grid files
        'gridfiles':  {
            'S': os.path.join(models_path_default,'grid.6'),
            'X': os.path.join(models_path_default,'grid.4')}
        }
        
    interpConfig = {
        "dmin" : 0.0,
        "dmax" : 2890.0,
        "nx" : 100,
        "nr" : 200,
    }
    
    # setup command-line args
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--receivers', type=str,help='receivers.dat file')
    parser.add_argument('--event', type=str,help='quakeML event file')
    parser.add_argument('--lonlat1', type=float,help='first profile point',nargs=2,required=True)
    parser.add_argument('--lonlat2', type=float,help='first profile point',nargs=2,required=True)
    parser.add_argument('--az', type=float,help='min max azimuth',nargs=2, default=[0.,360.])
    parser.add_argument('--ulvz-lonlat', dest='ulvz', type=float,help='ulvz position',nargs=2)
    
    # parser.add_argument("-dmax", help="maximum station-section deistance to plot stations", type=float, default=180.0)
    
    parser.add_argument("--phases", dest="phase_list", help="list of phases arrival to compute with taup and plot", nargs="+", type=str, default=["S", "Sdiff"])
    parser.add_argument("-d", dest="dist_range", help="distance bounds", nargs=2, type=float, default=[0.0,180.0])

    parser.add_argument("--ulvz", help="[lat,lon,radius,heigh] of an ulvz", nargs=4, type=float)

    
    parser.add_argument('-o', dest="outfile", type=str,help='out figure name')
    
    args = parser.parse_args()
    
    # getting station positions
    if args.receivers: stations = read_receivers(args.receivers)

    # getting event information
    if args.event: event = read_events(args.event)[0]

    if args.receivers and args.event:
        # filter on distance
        dmin,dmax = args.dist_range
        if dmin > 0.0 or dmax < 180.0:

            ev_orig = event.preferred_origin()
            evla,evlon = ev_orig.latitude, ev_orig.longitude

            get_dist = lambda row : locations2degrees(evla,evlon,row.lat,row.lon)

            stations["dist"] = stations.apply(get_dist, axis=1)

            stations = stations[(stations["dist"] >= dmin) & (stations["dist"] <= dmax)]

            if stations.empty:
                print("Warning : no stations respecting the distance criteria")


    origin = event.preferred_origin()

    stations["az"] = stations.apply(
        lambda stn: gps2dist_azimuth(stn.lat, stn.lon, origin.latitude, origin.longitude)[2],
        axis = 1
        )

    # filter on azimuth

    azmin,azmax = args.az

    if azmin < 0.: azmin += 360.
    if azmax < 0.: azmax += 360.
    stations = stations[(stations.az >= azmin) & (stations.az <= azmax)]
    

    # setting/computing limits points of profil
    
    lonlat1 = args.lonlat1
    lonlat2 = args.lonlat2

    
    # bellow is an attempt to automatically choose the end points of the section

    # if (not(any(lonlat1)) or not(any(lonlat2))) and args.receivers and args.event :
        
    #     # compute mean point of stations
    #     lat_mean_stn,lon_mean_stn = barycenter_on_sphere(stations["lat"],stations["lon"])
        
    #     origin = event.preferred_origin()
    #     _,az,b_az = gps2dist_azimuth(origin.latitude, lat_mean_stn, origin.longitude, lon_mean_stn)
                
    #     # compute average event-station distance
    #     dist_mean = locations2degrees(origin.latitude, lat_mean_stn, origin.longitude, lon_mean_stn)
        
    #     # points to be taken for section limits
    #     dist_event_side = 0.1*dist_mean
    #     dist_stn_side = 0.2*dist_mean
        
    #     dist_event_side_km = degrees2kilometers(dist_event_side)
    #     dist_stn_side_km = degrees2kilometers(dist_stn_side)

        
    #     event_position = geopy.Point(origin.latitude, origin.longitude)
    #     mean_stn_position = geopy.Point(lat_mean_stn,lon_mean_stn)
        
    #     point1 = geodesic(kilometers=dist_event_side_km).destination(event_position, bearing=b_az)
    #     point2 = geodesic(kilometers=dist_stn_side_km).destination(mean_stn_position, bearing=az)
        
    #     lonlat1 = (point1.longitude, point1.latitude)
    #     lonlat2 = (point2.longitude, point2.latitude)
        
    #     # test everything ok
        
    #     import cartopy.crs as ccrs
    #     import cartopy.feature as cfeature  
        
    #     fig = plt.figure()
        
    #     lat0,lon0 = compute_mean_point_stations_event(event,stations)
        
    #     proj = ccrs.Robinson(
    #         central_longitude = lon0,
    #     )   
        
    #     ax = fig.add_subplot(1, 1, 1, projection=proj)

    #     ax.coastlines()
    #     ax.add_feature(cfeature.LAND)
    #     ax.gridlines(linestyle=":", color="k")
        
    #     lon1,lat1 = lonlat1 
    #     lon2,lat2 = lonlat2
        
    #     print([lon1,lon2,origin.longitude,lon_mean_stn],[lat1,lat2,origin.latitude,lat_mean_stn])
        
    #     ax.scatter([lon1,lon2,origin.longitude,lon_mean_stn],[lat1,lat2,origin.latitude,lat_mean_stn], transform = ccrs.PlateCarree())
        
    #     ax.plot((lon1,lon2), (lat1,lat2), "g",transform=ccrs.Geodetic())
    #     ax.plot((lon1,lon_mean_stn), (lat1,lat_mean_stn), "b",transform=ccrs.Geodetic())
        
    #     ax.set_global()
        
    #     plt.show()
        
    #     exit()
        
        
    phases = args.phase_list

    ax = plot_model(args.lonlat1, args.lonlat2, modelConfig, interpConfig, vmax=4)
    
    if args.event and args.receivers:
        plot_paths(event, stations, ax, args.lonlat1, phases=phases)
        
    if args.ulvz : 
        # plot_ulvz_simple(args.ulvz, args.lonlat1, ax)
        lat_ulvz,lon_ulvz,radius,heigh = args.ulvz

        exag_heigh = 1.
        h = heigh * exag_heigh

        lonref,latref = lonlat1
        dist = locations2degrees(lat_ulvz,lon_ulvz,latref,lonref)

        ddist = radius * 360./(2*np.pi*3480.)

        x_list = np.linspace(dist-ddist, dist+ddist,10) * np.pi/180
        x_list = np.concatenate([x_list,x_list[::-1],[(dist-ddist)*np.pi/180]])

        y_list = np.concatenate([np.ones(10)*3480,np.ones(10)*(3480+h),[3480]])

        plt.plot(x_list,y_list,"r", lw=1,clip_on=False)

        
    plt.setp(ax, rorigin=0, rmin=3480., rmax=rEarth)
    # plt.legend(loc=10, bbox_to_anchor=(0., 0.4, 0.05, 0.5), title="Phases :")


    # append subax with profil

    lat0,lon0 = barycenter_on_sphere([args.lonlat1[1],args.lonlat2[1]],[args.lonlat1[0],args.lonlat2[0]])

    fig = plt.gcf()

    # add box with orthographic projection
    sub_ax = fig.add_axes([0.65, 0.55, 0.2, 0.2],projection=ccrs.Orthographic(central_latitude=lat0, central_longitude=lon0))
    sub_ax.add_feature(cfeature.LAND)
    sub_ax.add_feature(cfeature.OCEAN)

    sub_ax.scatter(stations["lon"], stations["lat"], marker=".", color="darkgreen", s = 10, transform = ccrs.PlateCarree())
    sub_ax.scatter(origin.longitude, origin.latitude, marker="*", color="red",s=100, transform = ccrs.PlateCarree(), ec="k", zorder=10)
    
    lon1,lat1 = args.lonlat1
    lon2,lat2 = args.lonlat2

    # sub_ax.scatter([lon1,lon2],[lat1,lat2], transform = ccrs.PlateCarree())
    
    sub_ax.plot((lon1,lon2), (lat1,lat2), "r-|",transform=ccrs.Geodetic())
        
    sub_ax.set_global()

    # set figure title 
    # dt_str = origin.time.datetime.strftime("%d %b %Y")
    # sub_ax.set_title(f"Event : {dt_str}")
    

    if args.outfile:
        print(f"Saving {args.outfile}")
        plt.savefig(args.outfile, dpi=500, bbox_inches='tight')
    else:
        plt.show()
    
    # saving options
    optionFile = f".{os.path.basename(__file__)}_lastUsageOptions.txt"
    print(f"Saving options used in {optionFile}")
    with open(optionFile, "w") as f:
        f.write(" ".join(sys.argv[1:]))