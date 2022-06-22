#!/usr/bin/env python

import numpy as np
from numpy.fft import fft

import matplotlib.pyplot as plt
plt.style.use("myStyle")

import pandas as pd

import re

import sys, os

import argparse 

# CONFIGURATION

# bounds frequencies (mHz)
f1 = 0.01
f2 = 0.5


def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("trace_file", help="ascii trace file",type=str)
    parser.add_argument("-o", dest="out_file", help="output figure name", type=str, default="")
    args = parser.parse_args()


    # reading trace
    data = np.loadtxt(args.trace_file)
    time = data[:,0]
    signal = data[:,1]
    N = len(signal)

    # sampling rate (Hz)
    sr = 1/(time[1] - time[0])


    # tapering (hanning)
    window = np.hanning(N)
    signalTappered = window*signal

    # fft
    signalFFT = fft(signalTappered)[:N//2]
    signalFFTmag = np.abs(signalFFT)
    freq = np.linspace(0.0, 0.5*sr, len(signalFFT))

    # plotting
    fig, ax = plt.subplots(tight_layout=True, figsize = (8, 5))
    fig.suptitle(f"{os.path.basename(args.trace_file)} - {int(time[-1]//3600+1)}h - Hanning taper")


    idx_f1 = np.argmin(np.abs(freq - f1))
    idx_f2 = np.argmin(np.abs(freq - f2))


    # plotting spectrum
    ax.plot(freq[idx_f1:idx_f2], signalFFTmag[idx_f1:idx_f2], ".-")


    ax.set_ylabel("Amplitude")
    ax.set_xlabel("Frequency (Hz)")

    # ++ showing / saving ++
    if args.out_file:
        if os.path.exists(args.out_file):
            if input(f"File {args.out_file} exists, overwrite ? [y/n] ") != "y":
                print(f"No file saved.")
                exit()
                
        print(f"Saving {args.out_file}")
        plt.savefig(args.out_file, dpi=500, bbox_inches='tight')

    else:
        plt.show()
        
if __name__ == "__main__":
    
    main()