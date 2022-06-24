#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 05 2022
@author: Sylvain Brisson sylvain.brisson@ens.fr

3d SH radiation pattern

"""

from turtle import width
import numpy as np
import matplotlib.pyplot as plt

plt.style.use("publication")

from numpy import cos, sin, pi
import pandas as pd

import os,sys
import argparse

# reading event from quakeML or CMTSOLUTION file
from obspy.core.event import read_events

# computing azimuth and gc distance
from obspy.geodetics import gps2dist_azimuth, locations2degrees

# computing ray paths
from obspy.taup import TauPyModel
model = TauPyModel(model='prem')


# reading receivers.dat
pathdir = os.path.dirname(__file__)
pathbase = os.path.abspath(os.path.join(pathdir, os.pardir))
sys.path.append(pathbase)

from common.my_io import read_receivers, to_obspy_inventory

# going from moment tensor to strike,dip,slip
from pyrocko import moment_tensor as pmt


# CONFIGURATION

kind = "SH"

#-----------------

def radiation_pattern_P(φf,δ,λ,φs,i):
    """Return the P radiation pattern, all angles in radians

    Args:
        phif (float): strike
        delta (float): dip
        lambda (float): slip angle
        phis (float): azimuth
        i (float): takeoff angle
    """
    
    dφ = φs - φf
        
    return cos(λ)*sin(δ)*sin(i)**2*sin(2*dφ) \
            - cos(λ)*cos(δ)*sin(2*i)*cos(dφ) \
            + sin(λ)*sin(2*δ)*(cos(i)**2 - sin(i)**2*sin(dφ)**2) \
            + sin(λ)*cos(2*δ)*sin(2*i)*sin(dφ)
            
def radiation_pattern_SH(φf,δ,λ,φs,i):
    """Return the SH radiation pattern, all angles in radians

    Args:
        phif (float): strike
        delta (float): dip
        lambda (float): slip angle
        phis (float): azimuth
        i (float): takeoff angle
    """
    
    dφ = φs - φf
        
    return cos(λ)*cos(δ)*cos(i)*sin(dφ) \
            + cos(λ)*sin(δ)*sin(i)*cos(2*dφ) \
            + sin(λ)*cos(2*δ)*cos(i)*cos(dφ) \
            - 0.5*sin(λ)*sin(2*δ)*sin(i)*sin(2*dφ)
            
            
def radiation_pattern_SV(φf,δ,λ,φs,i):
    """Return the SV radiation pattern, all angles in radians

    Args:
        phif (float): strike
        delta (float): dip
        lambda (float): slip angle
        phis (float): azimuth
        i (float): takeoff angle
    """
    
    dφ = φs - φf
        
    return sin(λ)*cos(2*δ)*cos(2*i)*sin(dφ) \
        - cos(λ)*cos(δ)*cos(2*i)*cos(dφ) \
        + 0.5*cos(λ)*sin(δ)*sin(2*i)*sin(2*dφ) \
        - 0.5*sin(λ)*sin(2*δ)*sin(2*i)*(1+sin(dφ)**2)


def norm(x): 
    """linear transofmation from [a,b] to [0,1]"""
    
    return (x - x.min())/(x.max()-x.min())

rp_function = {
    "P" : radiation_pattern_P,
    "SH" : radiation_pattern_SH,
    "SV" : radiation_pattern_SV
}

def aztkr2xyz(az,tk,r):
    """
    Convert (lat,lon) position into (x,y,z) on the unit sphere
    """
    theta = - az
    phi = tk
        
    x = r*sin(theta)*cos(phi)
    y = r*sin(theta)*sin(phi)
    z = - r*cos(theta)
    
    return x,y,z

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('event_file', type=str, help='obspy compatible event file')
    parser.add_argument('-o', dest="out_file", type=str,help='out figure name')
    parser.add_argument('-c', dest="component", type=str,help='P SH or SV', default="SH")
    parser.add_argument('-r', dest="rec_file", type=str,help='receivers.dat file path')

    args = parser.parse_args()
    
    events = read_events(args.event_file)
    event = events[0]
    
    # getting moment tensor 
    fm = event.preferred_focal_mechanism() or event.focal_mechanisms[0]
    mtensor = fm.moment_tensor.tensor
    mt_rtp = [mtensor.m_rr, mtensor.m_tt, mtensor.m_pp,
            mtensor.m_rt, mtensor.m_rp, mtensor.m_tp]
    
    # going from rtp to ned
    signs = [1, 1, 1, -1, 1, -1]
    indices = [1, 2, 0, 5, 3, 4]
    mt_ned = [sign * mt_rtp[ind] for sign, ind in zip(signs, indices)]
    
    # init pyrocko moment tensor
    m = pmt.MomentTensor(
        mnn = mt_ned[0], 
        mee = mt_ned[1],
        mdd = mt_ned[2], 
        mne = mt_ned[3], 
        mnd = mt_ned[4],    
        med = mt_ned[5],
    )
    
    print(m)
    
    # the 2 nodal planes, given in (strike,dip,slip)
    np1, np2 = m.both_strike_dip_rake()
        
    φ, δ, λ = np.array(np2)*pi/180.   
    
    # computing horizontal and vertical radiation pattern
    
    az = np.linspace(0.0, 2*pi, 100) # azimuth range
    tk = np.linspace(0.0, 2*pi, 100) # takeoff angle range
    
    az,tk = np.meshgrid(az, tk)
    
    rp =  norm(np.exp(np.exp(rp_function[kind](φ, δ, λ, az, tk))))
    # rp =  rp_function[kind](φ, δ, λ, az, tk)
    
    # plt.imshow(rp)
    # plt.show()
    # exit()
    
    # plotting
    
    fig, ax = plt.subplots(1, 1, figsize = (8,8),subplot_kw={'projection':'3d'})
    event_id = f"{event.magnitudes[0].mag:1.1f}_{event.origins[0].time.date.strftime('%d-%b-%Y')}"
    fig.suptitle(f"{kind} radiation pattern - {event_id}")
    
    x,y,z = aztkr2xyz(tk, az, rp)
    
    ax.plot_surface(x,y,z)
        
    ax.set_xticks([-1, 1],["S","N"])
    ax.set_yticks([-1, 1],["W","E"])
    ax.set_zticks([-1, 1],["Down","Up"])
    
    if args.rec_file:
        
        
        origin = event.preferred_origin() or event.origin[0]
        
        stations = read_receivers(args.rec_file)
        
        # plotting azimuth disctribution in horizontal plots
        stations["az"] = [
            gps2dist_azimuth(
            origin["latitude"],
            origin["longitude"],
            lat, lon)[1] * pi/180 for lat,lon in zip(stations["lat"],stations["lon"])
            ]
                
        # computin gc distance
        stations["dist"] = [
            locations2degrees(
                origin["latitude"],
                origin["longitude"],
                lat,lon) for lat,lon in  zip(stations["lat"],stations["lon"])
                ]
        
        # computing takeoff angle
        # obspy get_ray_path return a list of Arrival object 
        print("Computing takeoff-angles")
        stations["tk"] =  [
            model.get_ray_paths(0, d,phase_list=["S","Sdiff"])[0].takeoff_angle*pi/180 for d in stations["dist"]
            ]
        
        x_st,y_st,z_st = aztkr2xyz(stations["tk"], stations["az"], 1.0)
        
        for x_,y_,z_ in zip(x_st,y_st,z_st):
            plt.plot([0,x_],[0,y_],[0,z_], color = "r", zorder=10)
        
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
