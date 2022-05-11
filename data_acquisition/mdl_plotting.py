#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 05 2022
@author: Sylvain Brisson sylvain.brisson@ens.fr
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from obspy.imaging.beachball import beach
import cartopy.crs as ccrs
import cartopy.feature as cfeature

import obspy

import os

def plot_event(event):
    """Plot a Event object.

    Args:
        event: obspy Event object
    """
    
    origin = event.origins[0]
    
    
    fig = plt.figure()
    
    proj = ccrs.Orthographic(
    central_latitude = origin.latitude,
    central_longitude=origin.longitude
    )   
    
    ax = fig.add_subplot(1, 1, 1, projection=proj)
    
    ax.gridlines(linestyle=":", color="k")
    ax.set_global()
    ax.coastlines()
    ax.stock_img()

    
    focmec = event.focal_mechanisms[0]
    tensor = focmec.moment_tensor.tensor
    moment_list = [tensor.m_rr, tensor.m_tt, tensor.m_pp,
               tensor.m_rt, tensor.m_rp, tensor.m_tp]


    xy = proj.transform_point(origin.longitude, origin.latitude, src_crs=ccrs.Geodetic())
    
    b = beach(moment_list, xy=xy, width=1e6, linewidth=1, facecolor="k")
    b.set_zorder(5)
    ax.add_collection(b)
    
    
def plot_stations(stations, networks=None):
    """Plot stations in an interactive plotly window

    Args:
        stations: obspy stations catalogue
        nw: list of networks to plot (default : all of them)
    """
    
    # loading data into panda dataframe
    
    catalog = {
        "nw_code"       : [],
        "station_code"  : [],
        "latitude"      : [],
        "longitude"     : [],
        "start_time"     : [],
        "end_time"       : [],
    }
    
    for sta in stations:
        catalog["nw_code"].append(sta._code)
        data_sta = sta.stations[0].__dict__
        catalog["station_code"].append(data_sta["_code"])
        catalog["latitude"].append(data_sta["_latitude"])
        catalog["longitude"].append(data_sta["_longitude"])
        catalog["start_time"].append(data_sta["start_date"])
        catalog["end_time"].append(data_sta["end_date"])
        
    df = pd.DataFrame(catalog)
    df.drop_duplicates("station_code", inplace=True)
        
    # plotting it
    
    plt.figure(figsize=(10,7))

    ax = plt.axes(projection=ccrs.PlateCarree())

    ax.add_feature(cfeature.LAND)
    ax.add_feature(cfeature.OCEAN)
    ax.add_feature(cfeature.COASTLINE)
    
    if not(networks): networks = list(set(list(df["nw_code"])))
        
    for nw, df_nw in df.groupby('nw_code'):
        if nw in networks:
            ax.scatter(df_nw["longitude"], df_nw["latitude"], transform = ccrs.PlateCarree(), marker="v", ec="k", label = nw, s = 100)
    
    ax.legend(title="Networks")
    
    # extents
    lat_min = df["latitude"].min()
    lat_max = df["latitude"].max()
    lon_min = df["longitude"].min()
    lon_max = df["longitude"].max()
    ax.set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())

if __name__ == "__main__":
    
    stations = obspy.read_inventory("../../../data/reproducting_yuan_2017/6.4_23-Dec-2018/stations/*")
    
    plot_stations(stations, networks=["CX"])
    
    plt.show()
