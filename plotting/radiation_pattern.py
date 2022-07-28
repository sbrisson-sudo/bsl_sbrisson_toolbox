#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 05 2022
@author: Sylvain Brisson sylvain.brisson@ens.fr
"""

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 20})


# plt.style.use("publication")

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

nb_bins_az = 5
size_of_bins_az = 3.0

nb_bins_tk = 1
size_of_bins_tk = 1.0

center_az_np = True

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


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('event_file', type=str, help='obspy compatible event file')
    parser.add_argument('--only', type=str, help='SH | SV | P')
    parser.add_argument('-o', dest="out_file", type=str,help='out figure name')
    parser.add_argument('-r', dest="rec_file", type=str,help='receivers.dat file path')
    parser.add_argument('-az', dest="azimuth", type=float,help='azimuth for vertical plots')
    parser.add_argument('-baz', dest="back_azimuth", type=float,help='back azimuth for vertical plots')
    args = parser.parse_args()
    
    events = read_events(args.event_file)
    event = events[0]


    # getting if plotting SV SH P or only one
    if args.only:
        assert(args.only in ["SH", "SV", "P"])
    
    # getting moment tensor 
    fm = event.preferred_focal_mechanism() or event.focal_mechanisms[0]
    mtensor = fm.moment_tensor.tensor
    mt_rtp = [mtensor.m_rr, mtensor.m_tt, mtensor.m_pp,
            mtensor.m_rt, mtensor.m_rp, mtensor.m_tp]
    
    # going from rtp to ned
    signs = [1, 1, 1, -1, 1, -1]
    indices = [1, 2, 0, 5, 3, 4]
    mt_ned = [sign * mt_rtp[ind] for sign, ind in zip(signs, indices)]

    DYNCM2NM = 1e-7
    
    # init pyrocko moment tensor
    m = pmt.MomentTensor(
        mnn = mt_ned[0]*DYNCM2NM, 
        mee = mt_ned[1]*DYNCM2NM,
        mdd = mt_ned[2]*DYNCM2NM, 
        mne = mt_ned[3]*DYNCM2NM, 
        mnd = mt_ned[4]*DYNCM2NM,    
        med = mt_ned[5]*DYNCM2NM,
    )

    print(m)

    # the 2 nodal planes, given in (strike,dip,slip)
    np1, np2 = m.both_strike_dip_rake()
        
    φ, δ, λ = np.array(np2)*pi/180. 
    
    az_ref = 0.0 # direction for vertical plot  
    
    # loading receivers data
    
    if args.rec_file:
        
        origin = event.preferred_origin() or event.origin[0]
        
        stations = read_receivers(args.rec_file)
        
        manage_angle_np = (lambda a : a if a <= 180 else a - 360) if center_az_np else lambda a : a
        
        # plotting azimuth disctribution in horizontal plots
        stations["az"] = [
            manage_angle_np(gps2dist_azimuth(
            origin["latitude"],
            origin["longitude"],
            lat, lon)[1]) * pi/180 for lat,lon in zip(stations["lat"],stations["lon"])
            ]    
        
        az_ref = stations["az"].mean()    
            
        # plotting takeoff distribution in vertical plots
        
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
            (180. - model.get_ray_paths(0, d,phase_list=["S","Sdiff"])[0].takeoff_angle)*pi/180 for d in stations["dist"]
            ]

        # print(stations)
        
        # print(f"mean takeoff angle = {stations['tk'].mean()}")
    
    # computing horizontal and vertical radiation pattern
    
    if args.back_azimuth:
        az_ref = args.back_azimuth - 180.
    elif args.azimuth:
        az_ref = args.azimuth
    
    az = np.linspace(0.0, 2*pi, 100) # azimuth range
    tk = np.linspace(0.0, 2*pi, 100) # takeoff angle range

    if not(args.only):
    
        rp_p_az  = np.abs(radiation_pattern_P(φ, δ, λ, az, np.pi/2))
        rp_p_tk  = np.abs(radiation_pattern_P(φ, δ, λ, az_ref, tk))
        rp_sh_az = np.abs(radiation_pattern_SH(φ, δ, λ, az, np.pi/2))
        rp_sh_tk = np.abs(radiation_pattern_SH(φ, δ, λ, az_ref, tk))
        rp_sv_az = np.abs(radiation_pattern_SV(φ, δ, λ, az, np.pi/2))
        rp_sv_tk = np.abs(radiation_pattern_SV(φ, δ, λ, az_ref, tk))
        
        rp_max = max([rp.max() for rp in (rp_p_az, rp_p_tk, rp_sh_az, rp_sh_tk, rp_sv_az, rp_sv_tk)])
        
        rp_p_az  /= rp_max    
        rp_p_tk  /= rp_max    
        rp_sh_az /= rp_max        
        rp_sh_tk /= rp_max        
        rp_sv_az /= rp_max        
        rp_sv_tk /= rp_max        
            
        # plotting
        
        fig, axes = plt.subplots(2, 3,tight_layout=True, subplot_kw={'projection':'polar'}, figsize=(14,8))
        
        origin = event.preferred_origin() or event.origins[0]
        # fig.suptitle(f"Radiation patterns - Event : {event.origins[0].time.date.strftime('%d %b %Y')}, Mw = {event.magnitudes[0].mag:1.1f}, ({origin.latitude:0>2.1f}°,{origin.longitude:0>3.1f}°)")
        
        # colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
        # c1, c2 = colors[:2]

        c1,c2 = '#636EFA', '#EF553B'
        
        for ax in axes[0,:]: 
            ax.set_yticks([]) 
            ax.set_yticklabels([]) 
            ax.set_theta_zero_location("N")
            ax.set_theta_direction(-1)
            ax.set_xticks([0, pi/2, pi, 3*pi/2]) 
            ax.set_xticklabels(["N","E","S","W"])
            ax.set_ylim([0., 1.1])
            
        for ax in axes[1,:]: 
            ax.set_yticks([])
            ax.set_yticklabels([])
            ax.set_theta_zero_location("N")
            ax.set_theta_direction(-1)
            ax.set_xticks([0, pi/2, pi, 3*pi/2])
            ax.set_xticklabels(["Up",f"az={az_ref*180/np.pi:.1f}°","Down",""])
            ax.set_ylim([0., 1.1])

        

        axes[0,0].set_title("P")
        axes[0,0].text(-0.2,0.5, "Horizontal", rotation = 90, va = "center", transform=axes[0,0].transAxes)
        axes[1,0].text(-0.2, 0.5, "Vertical", rotation = 90, va = "center", transform=axes[1,0].transAxes)
        axes[0,0].plot(az, rp_p_az, c=c1)
        axes[1,0].plot(tk, rp_p_tk, c=c2)

        axes[0,1].set_title("SV")
        axes[0,1].plot(az, rp_sv_az, c=c1)
        axes[1,1].plot(tk, rp_sv_tk, c=c2)
        
        axes[0,2].set_title("SH")
        axes[0,2].plot(az, rp_sh_az, c=c1)
        axes[1,2].plot(tk, rp_sh_tk, c=c2)
        
        # add the receivers directions
        
        if args.rec_file:
            
            # not using pd.hist because of normalisation issues
            counts = stations["az"].value_counts(bins=nb_bins_az,normalize=True)
            bins = counts.index
            bins_center = [(b.left+b.right)/2 for b in bins]
            width = bins[0].right - bins[0].left
            
            size_of_bins_az = 0.8 / counts.values.max()
            
            for i in range(3):
                # stations.hist(column = "az", ax = axes[0,i])
                axes[0,i].bar(bins_center, counts.values*size_of_bins_az, alpha = 0.5, width=width, ec="k", color=c1)
                
            # plotting tk angles distribution

            counts = stations["tk"].value_counts(bins = nb_bins_tk,normalize=True)
            bins = counts.index
            bins_center = [(b.left+b.right)/2 for b in bins]
            width = bins[0].right - bins[0].left
            
            size_of_bins_tk = 0.8 / counts.values.max()
            
            for i in range(3):
                axes[1,i].bar(bins_center, counts.values*size_of_bins_tk, alpha = 0.5, width=width, ec="k", color=c2)

    else : 

        to_plot = args.only

        rp_fun = {
            "P" : radiation_pattern_P,
            "SV" : radiation_pattern_SV,
            "SH" : radiation_pattern_SH
        }[to_plot]

        rp_az  = np.abs(rp_fun(φ, δ, λ, az, np.pi/2))
        rp_tk  = np.abs(rp_fun(φ, δ, λ, az_ref, tk))
        
        rp_max = max([rp.max() for rp in (rp_az, rp_tk)])
        
        rp_az  /= rp_max    
        rp_tk  /= rp_max    

        # plotting
        
        fig, (ax1,ax2) = plt.subplots(1, 2,tight_layout=True, subplot_kw={'projection':'polar'}, figsize=(10,5))
        
        origin = event.preferred_origin() or event.origins[0]
        # fig.suptitle(f"Radiation pattern {args.only} - Event : {event.origins[0].time.date.strftime('%d %b %Y')}, Mw = {event.magnitudes[0].mag:1.1f}, ({origin.latitude:0>2.1f}°,{origin.longitude:0>3.1f}°)")
        
        # colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
        # c1, c2 = colors[:2]

        c1,c2 = '#636EFA', '#EF553B'
        
        
        ax1.set_yticks([]) 
        ax1.set_yticklabels([]) 
        ax1.set_theta_zero_location("N")
        ax1.set_theta_direction(-1)
        ax1.set_xticks([0, pi/2, pi, 3*pi/2]) 
        ax1.set_xticklabels(["N","E","S","W"])
        ax1.set_ylim([0., 1.1])
            
        ax2.set_yticks([])
        ax2.set_yticklabels([])
        ax2.set_theta_zero_location("N")
        ax2.set_theta_direction(-1)
        ax2.set_xticks([0, pi/2, pi, 3*pi/2])
        ax2.set_xticklabels(["Up",f"az={az_ref*180/np.pi:.1f}°","Down",""])
        ax2.set_ylim([0., 1.1])

        ax1.plot(az, rp_az, c=c1)
        ax2.plot(tk, rp_tk, c=c2)
        
        # add the receivers directions
        
        if args.rec_file:
            
            # not using pd.hist because of normalisation issues
            counts = stations["az"].value_counts(bins=nb_bins_az,normalize=True)
            bins = counts.index
            bins_center = [(b.left+b.right)/2 for b in bins]
            width = bins[0].right - bins[0].left
            
            size_of_bins_az = 0.8 / counts.values.max()

            ax1.bar(bins_center, counts.values*size_of_bins_az, alpha = 0.5, width=width, ec="k", color=c1)
                
            # plotting tk angles distribution

            counts = stations["tk"].value_counts(bins = nb_bins_tk,normalize=True)
            bins = counts.index
            bins_center = [(b.left+b.right)/2 for b in bins]
            width = bins[0].right - bins[0].left
            
            size_of_bins_tk = 0.8 / counts.values.max()
            
            ax2.bar(bins_center, counts.values*size_of_bins_tk, alpha = 0.5, width=width, ec="k", color=c2)
        
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
