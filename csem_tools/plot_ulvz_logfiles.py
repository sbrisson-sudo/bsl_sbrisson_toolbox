#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 31 2022
@author: Sylvain Brisson sylvain.brisson@ens.fr
"""

import numpy as np
import matplotlib.pyplot as plt

import os
import sys

def loncolatr2xyz(lon,colat,radius):
    """
    Convert (lat,lon) position into (x,y,z) on the unit sphere
    """
    theta = colat*np.pi/180
    phi = (90 + lon)*np.pi/180
        
    x = radius*np.sin(theta)*np.sin(phi)
    y = radius*np.cos(theta)
    z = radius*np.sin(theta)*np.cos(phi)
    
    return x,y,z

if __name__ == "__main__":
    
    if len(sys.argv) == 1:
        print("Usage : plot_ulvz_logfile <ulvz_log_file>")
        
    logfile = sys.argv[1]
        
    data = np.loadtxt(logfile)
    
    print(f"{data.shape()} points loaded")
    
    # go from colat, lon, radius to XYZ
    x,y,z = loncolatr2xyz(data[:,0], data[:,1], data[:,2])
    
    fig, ax = plt.subplots(1, 1, figsize = (8,8),subplot_kw={'projection':'3d'})
    plt.axis('off')

    
    ax.scatter(x,y,z,c = data[:,3], cmap = "Reds")
    
    ax.set_title("ULVZ anomaly : geometry")
    plt.show()
    
    