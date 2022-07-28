#!/usr/bin/env python3

from obspy import Stream,read
from obspy.geodetics import locations2degrees

import os

import argparse

# parsing arguments
parser = argparse.ArgumentParser()

parser.add_argument("dir_sac_files", help="directory containing the traces files to convert", type=str)
parser.add_argument("-o", dest="out_file", help="output file prefix", type=str, required = True)
parser.add_argument("-e", dest="ev_file", help="event file", type=str, required = True)

args = parser.parse_args()

# read event file 



st = Stream()

in_files = os.listdir(args.dir_sac_files)

for i,f in enumerate(in_files):

    print(f"Reading {f} ({i/len(in_files)*100:.0f}%).", end="\r")

    st2 = read(os.path.join(args.dir_sac_files,f),debug_headers=True)

    sac_metadata = st2[0].stats.sac

    evla,evlo,evde = sac_metadata["evla"],sac_metadata["evlo"],-sac_metadata["evel"]
    stla,stlo = sac_metadata["stla"],sac_metadata["stlo"]

    dist = locations2degrees(evla,evlo,stla,stlo)

    # add needed metadata

    stats = st2[0].stats
    stats.evla = evla 
    stats.evlo = evlo 
    stats.evde = evde 
    stats.distance = dist

    # merge streams

    st += st2
    


print()
            
out_file = args.out_file + ".pkl"

if os.path.exists(out_file):
    if input(f"File {out_file} exists, overwrite ? [y/n] ") != "y": exit()
    
print(f"Writting {out_file}")
st.write(out_file, format='PICKLE')