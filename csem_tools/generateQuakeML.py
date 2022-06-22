#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 2022
@author: Sylvain Brisson sylvain.brisson@ens.fr

Read source information in macromesh.json and build a quakeML event file.
"""

import numpy as np
import matplotlib.pyplot as plt

import os
import json

import argparse

from obspy import UTCDateTime
from obspy.core.event import (Catalog, Comment, Event, EventDescription,Origin, Magnitude, FocalMechanism, MomentTensor,Tensor, SourceTimeFunction)


from obspy.geodetics import FlinnEngdahl
_fe = FlinnEngdahl()


def _get_ressource_id(tag):
    return f"smi:local/sandwichCSEMsynthetic/BerkeleySeismologyLab/SylvainBrisson/{tag}"


def build_quakeML(coordinates, Mrtp):
    """Construct a quakeML file associated with the moment tensor Mrtp
    - coordinates : dictionary with latitude,longitude and radius entries"""
    
    origin_time = UTCDateTime()
    
    origin = Origin(
        resource_id=_get_ressource_id("origin"),
        time=UTCDateTime(),
        longitude=coordinates["longitude"],
        latitude=90-coordinates["colatitude"],
        # Depth is in meters.
        depth=(6371000.-coordinates["radius"]),
        )
    
    mag = Magnitude(
        resource_id=_get_ressource_id("magnitude"),
        mag=6.0, 
        # magnitude_type="Mb",
        # evaluation_status="preliminary",
        # origin_id=preliminary_origin.resource_id
        )
    
    Mrtp_idxs = {
        "m_rr":0,
        "m_pp":2,
        "m_tt":1,
        "m_rt":3,
        "m_rp":4,
        "m_tp":5,
    }
    
    tensor = Tensor(
        m_rr=Mrtp[Mrtp_idxs["m_rr"]],
        m_pp=Mrtp[Mrtp_idxs["m_pp"]],
        m_tt=Mrtp[Mrtp_idxs["m_tt"]],
        m_rt=Mrtp[Mrtp_idxs["m_rt"]],
        m_rp=Mrtp[Mrtp_idxs["m_rp"]],
        m_tp=Mrtp[Mrtp_idxs["m_tp"]],
    )
    
    mt = MomentTensor(
        tensor=tensor,
    )
    
    print(mt)
    
    foc_mec = FocalMechanism(
        resource_id=_get_ressource_id("focalMechanism"),
        moment_tensor = mt
    )
    
    ev = Event()
    
    ev.origins.append(origin)
    ev.focal_mechanisms.append(foc_mec)
    ev.magnitudes.append(mag)

    ev.preferred_origin_id = origin.resource_id.id
    ev.preferred_focal_mechanism_id = foc_mec.resource_id.id
    ev.preferred_magnitude_id = mag.resource_id.id
    
    print(">> Writting event.xml")
    ev.write(f"event.xml", format="quakeml") 

if __name__ == "__main__":
    
    # parsing arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", dest="macromesh_file", type=str, default="macromesh.json")
    args = parser.parse_args()
    
    # reading json file
    f = open(args.macromesh_file, "r")
    data = json.load(f)
    
    source = data["source"]
    coordinates = source["coordinates"]
    
    build_quakeML(Mrtp = source["Mrtp"])
    
    
    
    
    
    
