#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 2022
@author: Sylvain Brisson sylvain.brisson@ens.fr
"""



# unfinished

import numpy as np
import matplotlib.pyplot as plt

import os,sys

pathdir = os.path.dirname(__file__)
pathbase = os.path.abspath(os.path.join(pathdir, os.pardir))
sys.path.append(pathbase)

from common.my_io import read_receivers
from common.setup import models_path_default

from A3Dmodel_map_greatCircles import get_interp

from Model1D import Model1D
from ModelA3d import ModelA3d
import pyspl


def extract_model_given_depth(depth = 2891.):
    
    # load the A3d model
    model = ModelA3d(model_file)
    model.load_from_file()
    
    # fetch the spherical spline interpolant
    lons, lats, H = get_interp(grid_file)
    

if __name__ == "__main__":
    pass