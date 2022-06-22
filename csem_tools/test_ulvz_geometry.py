#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 31 2022
@author: Sylvain Brisson sylvain.brisson@ens.fr
"""

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature


import os


# Position of the anomaly                   
lon_ulvz = 0.           # longitude  phi_0  (deg) 
colat_ulvz = 30.          # colatitude theta_0 (deg)                 
# Size of the anomaly                       
radius_ulvz = 500.         # radius of the cylindrical anomaly (km)        
heigh_ulvz = 20.        # radius of the top of the cynlindrical anomaly 

CMB = 3480. # radius of the CMB
Re = 6371. # radius of the EArth

def in_cylinder(phi, theta, r):
    
    X_km = (phi-lon_ulvz)*(r*2*np.pi/360)*np.cos((theta-90.0)*np.pi/180)
    Y_km = (theta-colat_ulvz)*(r*2*np.pi/360)
    dist = np.sqrt(X_km**2+Y_km**2)
    
    print(X_km, Y_km, dist)
    
    if dist < radius_ulvz and CMB < r and r-CMB < heigh_ulvz:
        return 1 
    return 0 

def loncolatr2xyz(lon,colat,radius):
    """
    Convert (lat,lon) position into (x,y,z) on the unit sphere
    """
    theta = colat*np.pi/180
    phi = lon*np.pi/180
        
    x = radius * np.sin(theta)*np.sin(phi)
    y = radius * np.cos(theta)
    z = radius * np.sin(theta)*np.cos(phi)
    
    return x,y,z




v_in_cylinder = np.vectorize(in_cylinder)

if __name__ == "__main__":
    
    # lon = np.linspace(0.0, 360.0, 100)
    # global
    lon = np.linspace(-180.0, 180.0, 100)
    colat = np.linspace(0.0, 180.0, 100)
    
    #local around ulvz
    
    # lon = np.linspace(-20.0, 20.0, 100) + lon_ulvz
    # colat = np.linspace(-20.0, 20.0, 100) + colat_ulvz
    
    r = 3490.0
    
    lon,colat = np.meshgrid(lon,colat)
    
    r = np.ones(lon.shape) * r
    
    is_in_ulvz = v_in_cylinder(lon, colat, r)
    
    # plot it
    
    # x,y,z = loncolatr2xyz(lon, colat, r)
        
    # x = x.reshape(100*100)
    # y = y.reshape(100*100)
    # z = z.reshape(100*100)
    # is_in_ulvz  = is_in_ulvz.reshape(100*100)
    
    # fig, ax = plt.subplots(1, 1, figsize = (8,8),subplot_kw={'projection':'3d'})
    # # plt.axis('off')

    
    # points = ax.scatter(x,y,z,c = is_in_ulvz, cmap = "Reds")
    # plt.colorbar(points, shrink=0.25)
    # # ax.scatter(x,y,z, cmap = "Reds")
    
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Orthographic(central_latitude=90.-colat_ulvz, central_longitude=lon_ulvz))
    
    # ax = fig.add_subplot(1, 1, 1, projection=ccrs.Orthographic())
    
    ax.contourf(lon, 90-colat, is_in_ulvz, transform=ccrs.PlateCarree(), cmap = "Reds")
    
    ax.set_global()
    ax.coastlines()
    
    plt.show()
    
    

    
    ax.set_title("ULVZ anomaly : geometry")
    plt.show()
    
    
    pass