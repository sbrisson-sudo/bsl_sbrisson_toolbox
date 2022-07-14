#!/usr/bin/env python

"""
Filter a obspy stream object
"""

import numpy as np
from numpy.fft import fft

import matplotlib.pyplot as plt
# plt.style.use("myStyle")

import pandas as pd

import re

import sys, os

import argparse 

from obspy import read

# CONFIGURATION

# bounds frequencies (mHz)
f1 = 0.01
f2 = 0.5


def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("in_file", help="pickle traces file", type=str)
    parser.add_argument("periods", help="T1 and T2, periods (in s) to filter pass band with",type=float, nargs=2)
    args = parser.parse_args()
    
    # ++ reading stream file ++
    stream = read(args.in_file)
    nb_traces = stream.count()
    print(f">> Stream object loaded ({stream.count()} traces)")
    
    # ++ filtering ++
    T1,T2 = args.periods
    print(f">> Filtering between {T1}s and {T2}s...")
    stream.filter('bandpass', freqmin=1/T2, freqmax=1/T1)
    print("Done.")
        
    out_file = os.path.splitext(args.in_file)[0]+f"_{T1:.0f}s_{T2:.0f}s.pickle"
    print(f">> Saving {out_file}")
    
    if os.path.exists(out_file):
        if input(f"File {out_file} exists, overwrite ? [y/n] ") != "y": exit()
    
    stream.write(out_file, format='PICKLE')

        
if __name__ == "__main__":
    
    main()