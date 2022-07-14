#!/usr/bin/env python

"""
Convert multiple traces contain in a stream serialized file into ascii files
"""

from obspy import read

import numpy as np

import sys 

import matplotlib.pyplot as plt

import argparse

from tqdm import tqdm

def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("pkl_file", help="obspy pickle stream file",type=str)
    # parser.add_argument("synthetics", help="is it from csem ?",action="store_true")
    
    # parser.add_argument("--filter", help="periods bounds to pass band tapper", nargs=2, type=float, default=[None,None])
    args = parser.parse_args()
    
    print(">> Reading stream serialized file")
    st = read(args.pkl_file)

    print("Done.")

    for tr in tqdm(st, total=len(st)):

        metadata = tr.stats
        
        out_file = f"U{metadata.channel[-1]}_{metadata.network}_{metadata.station}"
        
        # if args.filter[0]:
        #     T1,T2 = args.filter
        #     st.filter('bandpass', freqmin=1/T2, freqmax=1/T1)
        #     out_file += f"_{T1:.0f}s_{T2:.0f}s"
        
        data = tr.data
        
        dt = 1/metadata.sampling_rate
        
        time = np.array([i*dt for i in range(len(data))])

        # print(dt*len(data))
        # exit()
        
        trace = np.array([time, data]).T
        
        # print(f"Saving {out_file}")
        
        np.savetxt(out_file, trace)
    
    # plt.plot(time, data)
    # plt.show()

if __name__ == "__main__":
    
    main()