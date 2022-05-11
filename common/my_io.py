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

    """Load the receivers information from a receivers.dat or a receivers.csv"""

    data = pd.read_csv(file, header = 2, sep = "\s+")
    data.rename(columns = {'lon:':'lon'}, inplace = True)
    
    return data

def to_obspy_inventory(stations):
    """Convert from the stations dataframe load from receivers.dat to a obspy station inventory"""
    stations_obspy = []
    for idx,row in stations.iterrows():
        stations_obspy.append(Station(row["stn"], row["lat"], row["lon"], 0.0))
    return Inventory(networks=[Network(row["nw"], stations_obspy)])

if __name__ == "__main__":
    
    stations = read_receivers("receivers.dat")
    
    print(stations)
    
    stations_obspy = to_obspy_inventory(stations)
    
    stations_obspy.plot()
    plt.show()
    