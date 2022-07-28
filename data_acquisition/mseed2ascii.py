#!/usr/bin/env python

"""
Convert a trace from mseed format to ascii 2 column file, one for each component
"""

from obspy import read

import numpy as np

import os

import argparse

from tqdm import tqdm

def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("seed_files", help="mseed trace file",type=str, nargs="+")
    parser.add_argument("--filter", help="periods bounds to pass band tapper", nargs=2, type=float, default=[None,None])
    parser.add_argument("--out-dir", dest="out_dir",help="output_directory", type=str, default="traces_ascii")
    args = parser.parse_args()

    if not(os.path.exists(args.out_dir)):
        print(f"Creating directory {args.out_dir}")
        os.mkdir(args.out_dir)
    
    for seed_file in tqdm(args.seed_files):

        st = read(seed_file)
        
        metadata = st[0].stats

        print(metadata)
        exit()
        
        out_file = f"U{metadata.channel[2]}_{metadata.network}_{metadata.station}"
        
        if args.filter[0]:
            T1,T2 = args.filter
            st.filter('bandpass', freqmin=1/T2, freqmax=1/T1)
            out_file += f"_{T1:.0f}s_{T2:.0f}s"
        
        data = st[0].data
        
        dt = metadata.sampling_rate
        
        time = np.array([i*dt for i in range(len(data))])
        
        trace = np.array([time, data]).T
        
        # print(f"Saving {out_file}")
    
        np.savetxt(os.path.join(args.out_dir,out_file), trace)
    
    # plt.plot(time, data)
    # plt.show()

if __name__ == "__main__":
    
    main()