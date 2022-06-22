#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 2022
@author: Sylvain Brisson sylvain.brisson@ens.fr
"""

import numpy as np
import matplotlib.pyplot as plt

import os

from obspy.core.event import read_events

import argparse


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("event_file", type=str)
    parser.add_argument("-o", dest="out_file", type=str)
    args = parser.parse_args()

    event = read_events(args.event_file)[0]
    print(event.short_str())
    origin = event.preferred_origin()
    # print(f"depth = {origin.depth/1000} {f'+- {origin.depth_errors.uncertainty/1000}' km")
    print(f"depth = {origin.depth/1000} km")
    
    event.plot(show=False)

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
