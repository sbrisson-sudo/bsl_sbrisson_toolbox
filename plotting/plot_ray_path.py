#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

from obspy.taup import TauPyModel
from obspy.taup.utils import get_phase_names

import argparse
import os

# parsing arguments
parser = argparse.ArgumentParser()

parser.add_argument("phase_list", help="list of phases arrival to plot path", nargs="+", type=str)
parser.add_argument("-d", dest="dist", help="distance in degree",nargs="+", type=float, default=[120.0])
parser.add_argument("--depth", dest="depth", help="depth in km", type=float, default=20.0)
parser.add_argument("--ulvz-dist", dest="ulvz_dist", help="distance source to ulvz", type=float, default=20.0)
parser.add_argument('-o', dest="out_file", type=str,help='out figure name')


args = parser.parse_args()


model = TauPyModel(model="iasp91")

obspy_phase_lists = ["ttp", "tts", "ttbasic", "tts+", "ttp+", "ttall"]
if len(args.phase_list) == 1 and args.phase_list[0] in obspy_phase_lists:
    phase_list = get_phase_names(args.phase_list[0])
else:
    phase_list = args.phase_list
    

arrivals = model.get_ray_paths(
    source_depth_in_km=args.depth, 
    distance_in_degree=args.dist[0], 
    phase_list=phase_list,
    )

ax = arrivals.plot_rays(show = False)
ax.annotate(f"{args.dist[0]:.0f}°", xy=(args.dist[0]*np.pi/180.0, 6371.0))


for d in args.dist[1:]:
    
    arrivals = model.get_ray_paths(
        source_depth_in_km=args.depth, 
        distance_in_degree=d, 
        phase_list=phase_list,
    )

    arrivals.plot_rays(show = False, ax = ax)
    
    ax.annotate(f"{d:.0f}°", xy=(d*np.pi/180.0, 6371.0))
    

depth = 3480.0
if args.ulvz_dist: 
    ax.scatter(args.ulvz_dist*np.pi/180.0, depth, color="r", zorder=5, ec="k", s=50)
    ax.annotate(f"{args.ulvz_dist:.0f}°", xy=(args.ulvz_dist*np.pi/180.0, depth))



ax.set_ylim([0.0, 6371.0])

if len(phase_list) > 10:    
    ax.legend(ncol=len(phase_list)//10, loc = "lower left")
    
# plt.title(f"Δ = {args.dist:.0f}°")    
    
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

