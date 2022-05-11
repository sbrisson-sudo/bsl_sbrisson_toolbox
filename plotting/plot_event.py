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
    args = parser.parse_args()

    event = read_events(args.event_file)[0]
    
    event.plot()
    plt.show()
