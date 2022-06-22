#!/usr/bin/env python3

# -------------------
# program to generate receivers at certain distances and azimuth from the source (relatively to the Norh pole)
# the source is at the equator (lon = 0.0, colat = 90.)
# the azimuth are the angle (source_northpole northpole_receiver)
# -------------------


import numpy as np
from numpy import pi,cos,sin,arccos

import warnings
warnings.filterwarnings("error")

# -- configuration

dist_deg = np.linspace(95,110,15)    # distance in degrees
az_deg = np.linspace(-10, 10,25)      # azimuth in degrees

# --


def generate_station_names(N):
    a,b,c,d = 65,65,65,65
    names = []
    for i in range(N):
        names.append(chr(a)+chr(b)+chr(c)+chr(d))
        d += 1
        if d==91:
            d = 65
            c += 1
        if c==91:
            c = 65
            b += 1
    return names
        
dist_deg = np.linspace(90,110,20)    # distance in degrees
az_deg = np.linspace(-10,10,25)      # azimuth in degrees
N_dist,N_az = len(dist_deg),len(az_deg)

dist_rad = dist_deg*pi/180.
az_rad = az_deg*pi/180.

# source coordinates
colat_source_deg = 90.0
lon_source_deg = 0.0

coφs = colat_source_deg*pi/180.
λs = lon_source_deg*pi/180.

# stations coodinates
colat_rad = np.zeros((N_dist,N_az))
lon_rad = np.zeros((N_dist,N_az))

for i,Δ in enumerate(dist_rad):
    for j,α in enumerate(az_rad):
        
        colat_rad[i,j] = arccos(cos(Δ)*cos(coφs) + sin(Δ)*sin(coφs)*cos(α))
        
        try:
            dλ = arccos((cos(Δ)-cos(colat_rad[i,j])*cos(coφs))/(sin(colat_rad[i,j])*sin(coφs)))
        except RuntimeWarning:
            dλ = np.pi

        if α < 0. : dλ *= -1
        
        lon_rad[i,j] = λs + dλ
        
lat_deg = (90. - colat_rad*180./pi).reshape(N_dist*N_az) 
lon_deg = (lon_rad*180./pi).reshape(N_dist*N_az)
stations = dict(zip(generate_station_names(N_dist*N_az), zip(lat_deg,lon_deg)))

# saving it into receivers.dat file

nw_code = "SY" # for synthetics

with open("receivers.dat", 'w') as out:
    header = ["Number of stations is:",len(lon_deg),"nw stn lat lon:"]
    out.writelines(f"{l}\n" for l in header)
    for name,(lat,lon) in stations.items():
        out.write(f"{nw_code:<2} {name[:4]:<4} {lat:8.4f}  {lon:8.4f}\n")
