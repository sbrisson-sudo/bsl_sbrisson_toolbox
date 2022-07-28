#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Thu Apr 07 2022
@author: Sylvain Brisson sylvain.brisson@ens.fr
"""
# importations

import argparse
import sys,os

pathdir = os.path.dirname(__file__)
pathbase = os.path.abspath(os.path.join(pathdir, os.pardir))
sys.path.append(pathbase)

from common.barycenter_sphere import barycenter_on_sphere



import matplotlib.pyplot as plt
import matplotlib
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature

import mpl_toolkits.axes_grid1 as axGrid


from obspy.core.event import read_events,Catalog

from obspy.imaging.beachball import beach

import pandas as pd

import warnings
warnings.filterwarnings('ignore')

# constantes

r_earth = 6371.0

def plot_event(event, ax, proj, s=5e5, facecolor="k"):
    """Plot Event BeachBall"""
    
    # get event info
    origin = event.preferred_origin() or event.origins[0]
    focmec = event.preferred_focal_mechanism() or event.focal_mechanisms[0]
    tensor = focmec.moment_tensor.tensor
    moment_list = [
        tensor.m_rr, tensor.m_tt, tensor.m_pp,
        tensor.m_rt, tensor.m_rp, tensor.m_tp]
    
    lon,lat = origin.longitude, origin.latitude
    
    xy = proj.transform_point(lon, lat, src_crs=ccrs.Geodetic())
    
    b = beach(moment_list, xy=xy, width=s, linewidth=1, facecolor=facecolor)
    b.set_zorder(5)
    ax.add_collection(b)

def plot_plates(ax):
    '''Plot plate boundaries from the UTIG PLATES collection.
    '''

    for file in ["ridges.shp", "trenches.shp", "transforms.shp"]:
        
        fname = os.path.join(os.path.dirname(__file__), f"data/{file}")
        
        plateBD_feature = ShapelyFeature(Reader(fname).geometries(),ccrs.PlateCarree())
        ax.add_feature(plateBD_feature, facecolor=(0,0,0,0), edgecolor='black', lw=2.5)

        plateBD_feature = ShapelyFeature(Reader(fname).geometries(),ccrs.PlateCarree())
        ax.add_feature(plateBD_feature, facecolor=(0,0,0,0), edgecolor='red', lw=1)


        
if __name__ == "__main__":
    
    # setup command-line args
    parser = argparse.ArgumentParser() 
        
    parser.add_argument('events', type=str,nargs="+",help='quakeML event files')

    parser.add_argument('-o', dest="outfile", type=str,help='out figure name')

    args = parser.parse_args()
    
    # plotting
        
    fig = plt.figure(figsize=(8,8))
            
    events = Catalog()
    for f in args.events:
        events.append(read_events(f)[0])

    # compute mangnitude size scale
    ev_mw, ev_depth = [],[]
    ev_lons,ev_lats = [],[]
    for ev in events:
        origin = ev.preferred_origin()
        ev_depth.append(origin.depth/1000)
        ev_mw.append(ev.preferred_magnitude().mag)
        ev_lats.append(origin.latitude)
        ev_lons.append(origin.longitude)

    # for Orthographic
    size_min = 1e5
    size_max = 5e5

    # for PlateCarre
    size_min = 1
    size_max = 3

    Mw_min,Mw_max = min(ev_mw), max(ev_mw)

    size_mw = lambda mw : (mw - Mw_min)/(Mw_max - Mw_min) * (size_max-size_min) + size_min

    depth_min,depth_max = min(ev_depth), max(ev_depth)

    cNorm  = matplotlib.colors.Normalize(vmin=depth_min, vmax=depth_max)
    depth_scalarMap = matplotlib.cm.ScalarMappable(norm=cNorm, cmap="viridis_r")
    
    lat0,lon0 = barycenter_on_sphere(ev_lats,ev_lons)
    
    # proj = LowerThresholdOthographic(
    #     central_latitude=lat0,
    #     central_longitude=lon0
    # )   

    # proj = ccrs.Robinson(
    #     central_longitude=lon0
    # )   

    # proj = ccrs.Mercator(
    #     central_longitude=lon0
    # )   

    proj = ccrs.PlateCarree(central_longitude=lon0)
    
    ax = fig.add_subplot(1, 1, 1, projection=proj)
    
    ax.coastlines()
    ax.add_feature(cfeature.LAND)
    ax.gridlines(linestyle=":", color="k")
    
    for event in events:
        plot_event(event,ax,proj,
        s=size_mw(event.preferred_magnitude().mag),
        facecolor=depth_scalarMap.to_rgba(event.preferred_origin().depth/1000))

    dax = axGrid.make_axes_locatable(ax)
    cax = dax.new_horizontal(size="4%", pad=0.05, axes_class=plt.Axes)
    fig.add_axes(cax)

    cbar = ax.get_figure().colorbar(depth_scalarMap, format='%.0fkm', label=f"Event depth", pad=0.05, cax=cax)
    
    # add plates
    plot_plates(ax)

    # add legend (not functional)

    # handles = []
    # labels = ["Mw 6.0", "Mw 7.0"]

    # for event,size in zip([events[0],events[1]],[size_min,size_max]):

    #     focmec = event.preferred_focal_mechanism() or event.focal_mechanisms[0]
    #     tensor = focmec.moment_tensor.tensor
    #     moment_list = [
    #         tensor.m_rr, tensor.m_tt, tensor.m_pp,
    #         tensor.m_rt, tensor.m_rp, tensor.m_tp]
        
    #     b = beach(moment_list, width=size, linewidth=1, facecolor="r")

    #     handles.append(beach)

    # fig.legend(handles=handles, labels=labels)

    # ax.set_global()
    # ax.set_extent([-20, 170, -40, 40])
    ax.set_xlim([-40, 30])
    ax.set_ylim([-20, 10])

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
