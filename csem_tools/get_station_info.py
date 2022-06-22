#!/usr/bin/env python3

import pandas as pd
import sys, os
from numpy import cos,sin,arccos,arcsin,pi

helpMessage = f"Usage : {os.path.basename(__file__)} <station file> <station names>\nOptions:\n-d lat lon : source\n-h : display this help message "
if len(sys.argv) < 3 : print(helpMessage) ; exit()

file = sys.argv[1]
argv = sys.argv[2:]

"""
Angular distance from (lat1, lat2, lon1, lon2)
"""
def angularDistance(θ1, λ1, θ2, λ2):
    θ1, θ2, λ1, λ2 = map(lambda x : x*pi/180, [θ1, θ2, λ1, λ2])
    return arccos(sin(θ1)*sin(θ2) + cos(θ1)*cos(θ2)*cos(λ1-λ2))* 180/pi

# parsing options
stations = []
srcLat, srcLon = None, None
for (i,arg) in enumerate(argv):
    if arg == "-h": 
        print(helpMessage)
        exit()
    if arg == "-d":
        srcLat = float(argv.pop(i+1))
        srcLon = float(argv.pop(i+1))

    else : stations.append(arg)

# reading receivers.dat
data = pd.read_csv(file, header = 2, sep = "\s+")
networks = dict(zip(data["stn"], data["nw"]))
positions = dict(zip(data["stn"], zip(data["lat"], data["lon:"])))

# giving station informations
for stn in stations:

    lat, lon = positions[stn]
    nw = networks[stn]

    if srcLat:
        print(f"{stn} : ({lat:.2f}°,{lon:.2f}°), Δ = {angularDistance(lat, lon, srcLat, srcLon):.0f}°, network : {nw}")
    else:
        print(f"{stn} : ({lat:.2f}°,{lon:.2f}°), network : {nw}")
