#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 2022
@author: Sylvain Brisson sylvain.brisson@ens.fr
"""

import numpy as np
import matplotlib.pyplot as plt

# plt.style.use("seaborn-paper")
# plt.rcParams.update({'font.size': 13})

import os
import sys

pathdir = os.path.dirname(__file__)
pathbase = os.path.abspath(os.path.join(pathdir, os.pardir))
sys.path.append(pathbase)

from common.setup import models_path_default

from Model1D import Model1D
from ModelA3d import ModelA3d
import pyspl
from UCBColorMaps import cmapSved, cmapXved

from Sphere import delaz, gc_minor

import pandas as pd

rEarth = 6371.0

def plot_model(modelConfig, interpConfig, vmax=3.0, param='S'):

    # config
    # rad_int1 = 3480
    # rad_int2 = 3880
    
    rad_int1 = 3480
    rad_int2 = 4500

    # define the grid to evaluate the model

    nh = 100
    lons = np.linspace(0.0, 360.0, nh)
    lats = np.zeros(nh)
    theta = np.linspace(0.0, 360, nh) # for plotting

    # get the SEMUCB-WM1 data

    nr_sem = 10
    r_sem = np.linspace(rad_int1, rad_int2, nr_sem) 

    # load the A3d model
    model = ModelA3d(modelConfig['modelfile'])
    model.load_from_file()
    coefs = model.get_parameter_by_name(param).get_values()
    
    # sample the model (polar coords in great-circle plane)
    # load the grid
    grid = np.loadtxt(modelConfig['gridfiles'][param], skiprows=1)
    # compute the sspl interpolant
    sspl = pyspl.SphericalSplines(grid[:,0], grid[:,1], grid[:,2])
    H = sspl.evaluate(lons.ravel(), lats.ravel())
    # compute the bspl interpolant
    bspl = pyspl.CubicBSplines(model.get_bspl_knots())
    V = bspl.evaluate(r_sem)
    # sample
    x_sem = (H * (V * coefs).T).T
    
    # to percent
    x_sem *= 100

        # + color map +
    # set up colormap (xi always recentered)
    if param == 'X':
        cmap, vmin = cmapXved(41, vmax)
        old_center = 1.0 - vmax / (vmax - vmin)
        # cmap = recenter_cmap(cmap, 0.5, old_center=old_center)
    elif param == 'S':
        cmap = cmapSved(41)  
    

    # get the prem data 

    # load prem 

    variables = [
        "radius", 
        "density", 
        "Vpv", 
        "Vsv", 
        "Qkappa", 
        "Qshear", 
        "Vph", 
        "Vsh", 
        "η",
    ]

    df = pd.read_csv("prem", sep = '\s+', skiprows=3, header=0, names = variables)

    nr_prem_up = 20
    nr_prem_bot = 20

    r_prem_up = np.linspace(rad_int2, 6371, nr_prem_up)
    r_prem_bot = np.linspace(0.0, rad_int1, nr_prem_bot)

    to_plot = "density"

    x_prem_up = np.tile(np.interp(r_prem_up, df["radius"]/1000, df[to_plot]), nh).reshape(nh,nr_prem_up).T
    x_prem_down = np.tile(np.interp(r_prem_bot, df["radius"]/1000, df[to_plot]), nh).reshape(nh,nr_prem_bot).T

    # plt.imshow(x_prem_up)
    # plt.show()
    # exit()

    # assemble data

    # r = np.linspace(0.0, 6371, nr_sem + nr_prem_up + nr_prem_bot)
    # x = np.vstack([x_prem_up, x_sem, x_prem_down])

  
    # ++ plotting ++
    
    # setup plot (figure and axes)
    fig, ax = plt.subplots(figsize=(10,10),subplot_kw={'projection': 'polar'})
    
    ax.set_theta_direction(-1)
    # ax.set_theta_offset(np.pi / 2.0 + delta*np.pi/360)

    nlev = 200

    # only one plot
    # pp = plt.contourf(theta*np.pi/180.0, r, x, cmap=cmap, levels=nlev, vmin=-vmax, vmax=vmax)
    
    # several contours
    plt.contourf(theta*np.pi/180.0, r_sem, x_sem, cmap=cmap, levels=nlev, vmin=-vmax, vmax=vmax)

    cmap = "binary"

# plt.contourf(theta*np.pi/180.0, r_prem_up, x_prem_up, cmap=cmap, levels=nlev,vmin=df[to_plot].min(),vmax=df[to_plot].max(), alpha=0.3)
# plt.contourf(theta*np.pi/180.0, r_prem_bot, x_prem_down, cmap=cmap, levels=nlev,vmin=df[to_plot].min(),vmax=df[to_plot].max(), alpha=0.3)

    plt.contourf(theta*np.pi/180.0, r_prem_up, x_prem_up, cmap=cmap, levels=nlev, alpha=0.2)
    plt.contourf(theta*np.pi/180.0, r_prem_bot, x_prem_down, cmap=cmap, levels=nlev, alpha=0.2)


    # add lines 

    def circle_line(r,lw=1):
        n_points = 60
        theta = np.linspace(0.0, 2*np.pi, n_points)
        plt.plot(theta, np.ones(n_points)*r, lw=lw, color="k")

    circle_line(rad_int1, lw=2)
    circle_line(rad_int2, lw=2)

    for r in np.linspace(rad_int1, rad_int2, 6):
        circle_line(r, lw=0.5)

    ax.vlines(np.linspace(0.0, 2*np.pi, 40), ymin=rad_int1, ymax=rad_int2, color="k", lw=0.5)



    ax.grid(False)

    ax.set_xticklabels([])
    ax.set_yticklabels([])
    
    return ax
    

if __name__ == "__main__":
    
    # parameter to plot (A3d descriptor name)
    param = 'S'
    
    # model configuration
    modelConfig = {
        # simply for naming of output files
        'name': '2.6S6X4',
        # A3d file
        'modelfile': os.path.join(models_path_default,'Model-2.6S6X4-ZeroMean.A3d'),
        # 1D reference model (needed for plotting xi sections)
        'refmodelfile': os.path.join(models_path_default,'Model-2.6S6X4-ZeroMean_1D'),
        # spherical spline grid files
        'gridfiles':  {
            'S': os.path.join(models_path_default,'grid.6'),
            'X': os.path.join(models_path_default,'grid.4')}
        }
        
    interpConfig = {
        "dmin" : 0.0,
        "dmax" : 2890.0,
        "nx" : 100,
        "nr" : 200,
    }

    ax = plot_model(modelConfig, interpConfig, vmax=4)
   
    plt.savefig("schema_csem_alpha0.2.png", dpi=600, bbox_inches='tight')

    plt.show()
    

    # if args.outfile:
    #     print(f"Saving {args.outfile}")
    #     
    # else:
    #     plt.show()
