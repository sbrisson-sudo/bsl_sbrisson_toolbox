#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 2022
@author: Sylvain Brisson sylvain.brisson@ens.fr
"""

import numpy as np
import matplotlib.pyplot as plt

from obspy import read_events

import argparse

import os

def write_macromesh_source_info(event):
    
    origin = event.preferred_origin()
    ev_lat,ev_lon,ev_dep = origin.latitude, origin.longitude, origin.depth
    
    tensor = event.preferred_focal_mechanism().moment_tensor.tensor
    
    Rt_m = 6371000.
    sourceDatfile = "source.dat"
    
    ev_colat = 90. - ev_lat 
    ev_radius = Rt_m - ev_dep
    
    factor_mt = 1e-15

    with open(sourceDatfile, 'w') as out:

        out.write(
        f"""# point source coordinates r (radius, m), theta (colatitude, °), phi (longitude, °)
{ev_radius}
{ev_colat:8.4f}
{ev_lon:8.4f}
# Centroid Moment Tensor ({factor_mt:.0e}N.m) (Mrr,Mtt,Mpp,Mrt,Mrp,Mtp)
{tensor.m_rr*factor_mt:.4e}
{tensor.m_tt*factor_mt:.4e}
{tensor.m_pp*factor_mt:.4e}
{tensor.m_rt*factor_mt:.4e}
{tensor.m_rp*factor_mt:.4e}
{tensor.m_tp*factor_mt:.4e}"""
        )
        
    print(f"Source information written in {sourceDatfile}")


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('event_file', type=str, help='obsoy compatible event file')
    args = parser.parse_args()
    
    event = read_events(args.event_file)[0]
    
    write_macromesh_source_info(event)
    
