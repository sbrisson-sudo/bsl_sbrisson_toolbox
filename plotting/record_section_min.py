#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 11 2022
@author: Sylvain Brisson sylvain.brisson@ens.fr
"""

import numpy as np
import matplotlib.pyplot as plt

import os
import sys

from obspy.core import read

from obspy.taup import plot_travel_times

def get_event_information(stream):
    """get event informatoin from the metadata of the stream first trace
    """
    trace = stream.traces[0]
    lat_event = trace.stats.evla
    lon_event = trace.stats.evlo
    depth_event = trace.stats.evde
    origin_time_event = trace.stats.event_origin_time
    
    return lat_event, lon_event, depth_event, origin_time_event

if __name__ == "__main__":
    
    
    in_file = sys.argv[1]
    
    stream = read(in_file)
    stream = stream.select(component = "T")
    
    evla, evlon, evde, evtime = get_event_information(stream)
    
    fig = plt.figure(tight_layout = True, figsize = (8,8))
        
    stream.plot(
        type='section', 
        dist_degree=True, 
        ev_coord = (evla, evlon),
        show=False, 
        fig=fig,
        reftime = evtime,
        )
    
    ax = fig.axes[0]

    d1, d2 = ax.get_xlim()
    plot_travel_times(
        min_degrees=d1,
        max_degrees=d2,
        source_depth=evde, 
        phase_list=["S", "Sdiff"], 
        ax = ax,
        show=False)
    
    
    plt.show()
