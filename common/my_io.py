#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 15 2022
@author: Sylvain Brisson sylvain.brisson@ens.fr

Receivers class
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from obspy.core.inventory.inventory import Inventory
from  obspy.core.inventory.network import Network
from obspy.core.inventory.station import Station

def read_receivers(file):

    """
    Load the receivers information from a receivers.dat tp a pandas dataframe
    """

    data = pd.read_csv(file, header = 2, sep = "\s+")
    data.rename(columns = {'lon:':'lon'}, inplace = True)
    
    return data

def write_receivers(stations, out_file = "receivers.dat"):
    """
    Write the dataframe to a receivers.dat file
    """

    with open(out_file, 'w') as f:

        f.write(f"Number of stations:\n{len(stations.index)}\nnw stn lat lon:\n")

        for index, row in stations.iterrows():

            f.write(f"{row['nw']:<2} {row['stn'][:4]:<4} {row['lat']:8.4f}  {row['lon']:8.4f}\n")



def to_obspy_inventory(stations):
    """Convert from the stations dataframe load from receivers.dat to a obspy station inventory"""

    stations_obspy = []
    # networks_obspy = []

    for _,row in stations.iterrows():

        stations_obspy.append(Station(row["stn"], row["lat"], row["lon"], 0.0))

        # networks_obspy.append(Network(row["nw"]))

    return Inventory(networks=[Network("SY", stations_obspy)])

if __name__ == "__main__":
    
    stations = read_receivers("receivers.dat")
    
    print(stations)
    
    stations_obspy = to_obspy_inventory(stations)
    
    stations_obspy.plot()
    plt.show()
    