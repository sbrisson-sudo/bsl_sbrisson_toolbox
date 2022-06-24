#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 05 2022
@author: Sylvain Brisson sylvain.brisson@ens.fr
"""

import numpy as np
import matplotlib.pyplot as plt

import os
import argparse

from obspy.core.event import read_events
from obspy.imaging.source import plot_radiation_pattern, _setup_figure_and_axes

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('event_file', type=str, help='obspy compatible event file')
    parser.add_argument('-o', dest="out_file", type=str,help='out figure name')
    args = parser.parse_args()
    
    events = read_events(args.event_file)
    event = events[0]
    
    # -- obspy default plotting routine
    event.plot(show=False)    
    # event.plot(kind = ['beachball', 'p_sphere', 's_sphere'])
    
    # -- custom plotting
    # fm = event.preferred_focal_mechanism() or event.focal_mechanisms[0]
    # mtensor = fm.moment_tensor.tensor
    # mt = [mtensor.m_rr, mtensor.m_tt, mtensor.m_pp,
    #       mtensor.m_rt, mtensor.m_rp, mtensor.m_tp]
    
    # kind = ['p_sphere', 's_sphere']
    
    # # fig, axes, kind_ = _setup_figure_and_axes(kind)
    # # plot_radiation_pattern(mt, kind=kind, coordinate_system='RTP', fig=fig, show=False)
    
    # plot_radiation_pattern(mt, kind=kind, coordinate_system='RTP', show=False)
    
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
