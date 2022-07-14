#!/usr/bin/env python


"""
Compute the mean frequency specttul of all the ASCII traces files found"""

import numpy as np
from numpy.fft import fft

import matplotlib.pyplot as plt
# plt.style.use("myStyle")

import pandas as pd

import re

import sys, os

import argparse 

from tqdm import tqdm

# CONFIGURATION

re_trace_file = r"U(?P<comp>[LRT])_(?P<nw>[a-zA-Z0-9_]{2})_(?P<stn>[a-zA-Z0-9_]{3,4})"

# bounds frequencies (mHz)
f1 = 0.01
f2 = 0.5


def main():
    
    parser = argparse.ArgumentParser()
    # parser.add_argument("trace_file", help="ascii trace file",type=str)
    parser.add_argument("-o", dest="out_file", help="output figure name", type=str, default="")
    args = parser.parse_args()

    trace_files = []
    nb_file_max = 500
    i = 0

    for f in os.listdir():

        match = re.search(re_trace_file, f)

        if match:
            trace_files.append(f)

            i += 1
            if i > nb_file_max:
                break

            # print(match["comp"])
            # print(match["nw"])
            # print(match["stn"])

    # reading first trace
    data = np.loadtxt(trace_files[0])
    time = data[:,0]
    signal = data[:,1]

    N = len(signal)
    N = 18000

    signalFFT = fft(signal)[:N//2]

    Nfft = len(signalFFT)

    # sampling rate (Hz)
    sr = 1/(time[1] - time[0])

    freq = np.linspace(0.0, 0.5*sr, len(signalFFT))


    # tapering window (hanning)
    window = np.hanning(N)

    # loading data and computing FFT
    ffts_signals = np.zeros((len(trace_files),Nfft))
    print(">> Loading data and computing FFTs")
    for i,f in tqdm(enumerate(trace_files), total=len(trace_files)):

        data = np.loadtxt(f)
        signal = window*data[:N,1]

        ffts_signals[i,:] = np.abs(fft(signal)[:N//2])


    # computing mean
    print(">> Computing mean spectrum")
    mean_fft = np.mean(ffts_signals, axis=0)

    print(freq.shape)
    print(mean_fft.shape)

    # saving data
    print(">> Saving mean_spectrum.asc")
    stacked = np.dstack((freq, mean_fft)).reshape((freq.shape[0],2))
    np.savetxt("mean_spectrum.asc", stacked)


    # plotting
    fig, ax = plt.subplots(tight_layout=True, figsize = (8, 5))
    # fig.suptitle(f"{os.path.basename(args.trace_file)} - {int(time[-1]//3600+1)}h - Hanning taper")

    idx_f1 = np.argmin(np.abs(freq - f1))
    idx_f2 = np.argmin(np.abs(freq - f2))


    # plotting spectrum
    # ax.plot(freq[idx_f1:idx_f2], mean_fft[idx_f1:idx_f2], ".-")
    ax.plot(freq[idx_f1:idx_f2], mean_fft[idx_f1:idx_f2])
    # ax.plot(freq, mean_fft)


    ax.set_ylabel("Amplitude")
    ax.set_xlabel("Frequency (Hz)")

    ax.set_xscale('log')

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