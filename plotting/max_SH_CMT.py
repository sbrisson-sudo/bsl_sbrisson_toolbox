#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 01 2022
@author: Sylvain Brisson sylvain.brisson@ens.fr
"""

import numpy as np
from numpy import cos,sin,pi
import matplotlib.pyplot as plt

import os,sys

from scipy.optimize import minimize

from pyrocko import moment_tensor as pmt

pathdir = os.path.dirname(__file__)
pathbase = os.path.abspath(os.path.join(pathdir, os.pardir))
sys.path.append(pathbase)

from csem_tools.generateQuakeML import build_quakeML


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
            

if __name__ == "__main__":
    
    coord = {
        "colatitude" : 90.,
        "longitude" : 0.,
        "radius" : 300.0
    }
    
    tk = 165.
    az = 0.
    
    def fun(x):
        φ,δ,λ = x
        return -radiation_pattern_SH(φ,δ,λ,az,tk)
        
    res = minimize(
        fun,
        [0.,0.,0.]
    )
    
    print(res)
        
    strike,dip,slip = np.array(res.x)*180/pi 
    
    print(strike,dip,slip)
    
    # convert from strike-dip_slip to moment tensor
    
    magnitude = 6.  # Magnitude of the earthquake
    m0 = pmt.magnitude_to_moment(magnitude)  # convert the mag to moment


    m = pmt.MomentTensor(
       strike = strike, 
       dip = dip,
       rake = slip,
       scalar_moment=m0
    )
    
    Mrtp = m.m6_up_south_east()
    
    
    
    build_quakeML(coord, Mrtp)
    
