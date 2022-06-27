#!/usr/bin/env python3

import numpy as np


import json
import sys, os
import re
from functools import reduce

import argparse 

def bool_f2py(bool_f):
    
    if bool_f == ".true.":
        return True
    elif bool_f == ".false.":
        return False
    raise Exception(f"'{bool_f}' : Invalid entry, must be .true. or .false.")

def bool_py2f(bool_py):
    return "T" if bool_py else "F"

codeNatureInterfaces = {
        2 : "solide-solide",
        8 : "DtN",
    }

def number_of_nodes(nelh,degh,nelv,degv):
    """Number of nodes (if only conforming interfaces), no mortars"""
    n_nodes = 0
    for nv,dv in zip(nelv,degv):
        n_nodes += nelh**2*nv * degh*dv**2
    return n_nodes

def parse_macromesh(macromesh_file, legacy_fmt=False):

    f = open(macromesh_file, "r")
    config = {}

    # reading file
    # code horrible pour lire fichier horrible

    # header
    for _ in range(8): s = f.readline().strip()

    # database dir
    config["database_dir"] = f.readline().strip()

    # backup stuff
    for _ in range(15-10+1): s = f.readline().strip()

    # sem macromesh.dat header
    for _ in range(43-16+1): s = f.readline().strip()
    # reading sem options

    config["npro"] = int(f.readline().strip()[15:])
    config["square_npro"] = bool_f2py(f.readline().strip()[15:])
    config["nb_proc"] = config["npro"]**2*(1 if config["square_npro"] else 2)

    # SEM part
    config["SEM"] = {}
    config_sem = config["SEM"]
    config_sem["SEM_mesh"] = {}
    config_sem_mesh = config_sem["SEM_mesh"]

    nregion =  int(f.readline().strip()[15:])
    parametrique = bool_f2py(f.readline().strip()[15:])

    ncour = int(f.readline().strip()[15:])
    config_sem_mesh["nb_rings"] = ncour

    f.readline()

    depthInterfaces = []
    for i in range(ncour+1):
        depthInterfaces.append(float(f.readline().strip()))
    config_sem_mesh["interfaces_radius"] = depthInterfaces

    f.readline()
    natureInterfaces = []
    for i in range(ncour+1):
        natureInterfaces.append(int(f.readline().strip()))
    config_sem_mesh["interfaces_kind"] = natureInterfaces

    f.readline()
    geomInterfaces = []
    for i in range(ncour): 
        geomInterfaces.append(int(f.readline().strip()))
    config_sem_mesh["rings_geometry"] = geomInterfaces

    config_sem_mesh["nelh"] = int(f.readline().strip()[15:])

    f.readline()
    nelv = list(map(int, [f.readline().strip() for _ in range(ncour)]))
    config_sem_mesh["rings_nelv"] = nelv

    config_sem_mesh["degh"] = int(f.readline().strip()[15:])

    f.readline()
    degv = list(map(int, [f.readline().strip() for _ in range(ncour)]))
    config_sem_mesh["rings_degv"] = degv

    f.readline()
    anisoLayer = []
    for i in range(ncour): 
        anisoLayer.append(bool_f2py(f.readline().strip()))
    config_sem_mesh["rings_anisotropy"] = anisoLayer

    f.readline()

    config_sem["earth_model"] = f.readline().strip()[15:]
    config_sem["gravity"] = bool_f2py(f.readline().strip()[15:])
    config_sem["attenuation"] = bool_f2py(f.readline().strip()[15:])
    
    if not(legacy_fmt):
        f.readline()
        config_sem["ulvz"] = bool_f2py(f.readline().strip()[15:])
        config_sem["ulvz_method"] = f.readline().strip()[15:]

    for i in range(3): f.readline() # paraview

    f.readline()
    f.readline()
    f.readline()
    f.readline()

    # source
    config["source"] = {}
    config_source = config["source"]

    srcRadius = float(f.readline().strip())
    srcCoLat = float(f.readline().strip())
    srcLon = float(f.readline().strip())

    config_source["coordinates"] = {
        "longitude" : srcLon,
        "colatitude" : srcCoLat,
        "radius" : srcRadius,
    }

    f.readline()

    Mrtp = []
    for i in range(6): 
        Mrtp.append(float(f.readline().strip()))

    # macromesh order Mrr,Mtt,Mpp,Mrt,Mrp,Mtp
    Mrtp_dict = {
        "m_rr" : Mrtp[0],
        "m_pp" : Mrtp[2],
        "m_tt" : Mrtp[1],
        "m_rt" : Mrtp[3],
        "m_rp" : Mrtp[4],
        "m_tp" : Mrtp[5]
    }
    config_source["Mrtp"] = Mrtp_dict
    
    # compute moment and moment magnitude
    # normalisation unknown, not possible 
    
    # M0 = 1/np.sqrt(2) * np.sqrt(np.array(Mrtp)**2)
    # Mw = np.log10(M0)/1.5 - 10.73
    
    
    
    
        
    f.readline()
    config_source["geocentric_correction"] = bool_f2py(f.readline().strip())

    f.readline()
    f.readline()
    f.readline()

    config_source["time"] = {}
    config_time = config_source["time"]

    f.readline()
    config_time["dt"] = float(f.readline().strip())

    f.readline()
    config_time["nb_iter"] = int(f.readline().strip())

    f.readline()
    config_time["rate_output"] = int(f.readline().strip())

    f.readline()
    config_time["rate_save_disk"] = int(f.readline().strip())

    f.readline()
    config_time["rate_backup"] = int(f.readline().strip())

    f.readline()
    f.readline()

    f.readline()
    sourceName = f.readline().strip()

    if sourceName == "heavis":

        f.readline()
        f.readline()
        f.readline()

        sourceInfo = [
            float(f.readline().strip()) for i in range(4)
        ]

    else: 
        f.readline()
        sourceInfo = float(f.readline().strip())
        for i in range(5): f.readline()

    if sourceName == "heavis":  
        time_function = {
            "name" : "heavis",
            "corner_frequencies" : sourceInfo,
        }
    elif sourceName == "ricker":
        time_function = {
            "name" : "ricker",
            "central_frequency" : sourceInfo,
        }
    else:
        raise Exception(f"{sourceName} : unknown source name")

    config_source["time_function"] = time_function

    f.readline()
    config_source["origin_time"] = float(f.readline().strip())

    f.readline()
    f.readline()
    f.readline()

    config["modes"] = {}
    config_modes = config["modes"]

    config_modes["lmax"] = int(f.readline()[20:].strip())
    config_modes["nmax"] = int(f.readline()[20:].strip())

    nbDtn = int(f.readline().strip()[20:])

    modelesModes = []
    attenuationModes = []

    f.readline()

    for i in range(nbDtn):
        modelesModes.append(f.readline().strip()[21:])
        attenuationModes.append(bool_f2py(f.readline().strip()[20:]))

    assert(nbDtn == 2)
    config_modes["inner_sphere"] = {
        "file_name" : modelesModes[0],
        "attenuation" : attenuationModes[0],
    }
    config_modes["outer_shell"] = {
        "file_name" : modelesModes[1],
        "attenuation" : attenuationModes[1],
    }


    fmin = float(f.readline().strip()[20:])
    fmax_1_db = float(f.readline().strip()[20:])
    fmax_2_db = float(f.readline().strip()[20:])

    config_modes["database_freq"] = {
        "fmin" : fmin,
        "fmax1" : fmax_1_db,
        "fmax2" : fmax_2_db
    }

    # interp.dat
    f.readline()
    f.readline()
    f.readline()
    f.readline()
    f.readline()
    f.readline()

    config["interpolation"] = {}
    config_interp = config["interpolation"]

    config_interp["betaPerElmt"] = float(f.readline().split()[0])
    config_interp["Mcol"] = float(f.readline().split()[0])
    config_interp["GllPerParallel"] = float(f.readline().split()[0])
    config_interp["degPhi"] = int(f.readline().split()[0])
    config_interp["degTheta"] = int(f.readline().split()[0])

    fmax_1 = float(f.readline().split()[0])
    fmax_2 = float(f.readline().split()[0])

    config_modes["usage_freq"] = {
        "fmax1" : fmax_1,
        "fmax2" : fmax_2
    }

    # yannos.dat
    f.readline()

    config["yannos"] = {}

    f.readline()
    config["yannos"]["force_fmin"] = bool_f2py(f.readline().strip())

    f.readline()
    config["yannos"]["cancel_gravity"] = bool_f2py(f.readline().strip())

    f.readline()
    config["yannos"]["use_terf"] = bool_f2py(f.readline().strip())

    f.readline()
    config["yannos"]["check_modes"] = bool_f2py(f.readline().strip())

    f.readline()
    config["yannos"]["treshold_ray"] = float(f.readline().strip())

    f.readline()
    config["yannos"]["l_startlevel"] = int(f.readline().strip())

    # etc
    # yannos.dat
    f.readline()

    f.readline()
    t_start_coup = float(f.readline().strip())

    f.readline()
    t_end_coup = float(f.readline().strip())

    f.close()

    # Try to read modes models bounds

    for i,fmodel in enumerate(modelesModes):

        if os.path.exists(fmodel):

            with open(fmodel, "r") as f:
                lines = f.readlines()
                modelModeName = lines[0].strip()
                rmin = float(lines[6].split()[0])/1000
                rmax = float(lines[-1].split()[0])/1000
                modelModeBounds = (rmin,rmax)

            key = "inner_sphere" if i == 0 else "outer_shell"

            config_modes[key]["earth_model"] = modelModeName
            config_modes[key]["bounds"] = modelModeBounds

    return config

def pretty_print(config):
    """Print macromesh.dat configuration"""

    print(f">> Running over {config['nb_proc']} cpus")

    print("\n>> SEM parameters:")
    config_sem = config['SEM']
    print(f"- earth model: {config_sem['earth_model']}")
    if 'ulvz' in config_sem.keys():
        if config_sem['ulvz']:
            print(f"  + ULVZ (ulvz method = {config_sem['ulvz_method']})")
        else:
            print("- no ULVZ")
    print(f"- depth SEM shell : {config_sem['SEM_mesh']['interfaces_radius']}")
    print(f"- nb of elements : {config_sem['SEM_mesh']['nelh']**2*sum(config_sem['SEM_mesh']['rings_nelv'])}")
    print(f"""- nb of nodes : {number_of_nodes(
config_sem['SEM_mesh']['nelh'],
config_sem['SEM_mesh']['degh'],
config_sem['SEM_mesh']['rings_nelv'],
config_sem['SEM_mesh']['rings_degv'])}""")

    print("\n>> DtN operators database:")
    freq_db = config['modes']['database_freq']

    print(f"- fmin = {freq_db['fmin']:.2e} = {1/freq_db['fmin']:.1f}s")
    print(f"- fmax 1 (frequency domain) = {freq_db['fmax1']:.2e} = {1/freq_db['fmax1']:.1f}s")
    print(f"- fmax 2 (spectral domain) =  {freq_db['fmax2']:.2e} = {1/freq_db['fmax2']:.1f}s")

    for key in ["inner_sphere", "outer_shell"]:

        config_shell = config['modes'][key]

        print(f" - {key} : model file =  {config_shell['file_name']}")

        if 'earth_model' in config_shell.keys():
            print(f"    model name =  {config_shell['earth_model']}")
            print(f"    model bounds radius = {config_shell['bounds']}km")

    print(">> DtN operators usage:")
    usage_freq = config['modes']['usage_freq']
    print(f"  - fmax coupling (permanent regime) {usage_freq['fmax1']:.2e} = {1/usage_freq['fmax1']:.1f}s")
    print(f"  - fmax coupling (transitory regime) {usage_freq['fmax2']:.2e} = {1/usage_freq['fmax2']:.1f}s")

    source = config['source']
    coord = source['coordinates']
    print("\n>> Source parameters:")
    print("- position:")
    print(f"    - radius: {coord['radius']/1000}km (depth ~ {6370-coord['radius']/1000:.1f}km)")
    print(f"    - latitude: {90 - coord['colatitude']:.1f}°")
    print(f"    - longitude: {coord['longitude']:.1f}°")

    time_fct = source['time_function']
    print(f"- fonction source: {time_fct['name']}")
    if time_fct["name"] == "ricker":
        print(f"    - central frequency = {time_fct['central_frequency']}Hz = {1/time_fct['central_frequency']:.1f}s")
    else:
        for i in range(4):
            f = time_fct['corner_frequencies'][i]
            print(f"    - f{i} = {f:.2e}Hz = {1/f:.1f}s")
    print(f"- t0 = {source['origin_time']}s\n")


    # checks frequencies

    print(">> Checks (Yann recommandations)")

    ok = True

    assert(time_fct['name'] == "heavis")
    fmax_source = time_fct['corner_frequencies'][3]

    fmax_1,fmax_2,fmax_1_db,fmax_2_db = usage_freq['fmax1'],usage_freq['fmax2'],freq_db['fmax1'],freq_db['fmax2']

    if 1.25*fmax_source > fmax_1:
        print("- 1.25fmax_source < fmax_coupling not respected")
        ok = False

    if 2.*fmax_source < fmax_1:
        print("- 2fmax_source < fmax_coupling : you can decrease fmax_coupling")
        ok = False

    if fmax_1 > fmax_1_db:
        print("- fmax_coupling > fmax_database not respected, database not compatible")
        ok = False

    if fmax_2 > fmax_2_db:
        print("- fmax_coupling_trans > fmax_database_2 not respected, database not compatible")
        ok = False

    if 1.5*fmax_1 > fmax_2 or 2*fmax_1 < fmax_2:
        print("- 1.5fmax_coupling_trans < fmax_coupling_trans < 2fmax_coupling_trans not respected")
        ok = False

    if ok: print("Frequencies ok.")


def write_json(config, in_file):

    # writting json format
    if in_file == "macromesh.dat":
        out_file = "macromesh.json"
    else:
        outfile = os.path.join(os.basename(in_file), "macromesh.json")

    print(f"\n>> Writting {out_file}")
    config_json = json.dumps(config, indent = 3)
    with open(out_file, "w") as out:
        out.write(config_json)


def write_quakeML(config):
    
    # obspy imports
    from obspy import UTCDateTime
    from obspy.core.event import Catalog, Comment, Event, EventDescription,Origin, Magnitude, FocalMechanism, MomentTensor,Tensor, SourceTimeFunction

    resourceId = "synthetic_event"
    origin_time = UTCDateTime("2000/01/01")
    Rt = 6371000.0
    
    # Moment magnitude calculation in dyne * cm.

    source = config['source']
    Mrtp = source['Mrtp']
    coord = source['coordinates']

    m_0 = 1.0 / np.sqrt(2.0) * np.sqrt(
        Mrtp['m_rr'] ** 2 +
        Mrtp['m_tt'] ** 2 +
        Mrtp['m_pp'] ** 2 +
        2.0 * Mrtp['m_rt'] ** 2 +
        2.0 * Mrtp['m_rp'] ** 2 +
        2.0 * Mrtp['m_tp'] ** 2)
    m_w = 2.0 / 3.0 * (np.log10(m_0) - 16.1)
        
    origin = Origin(
        resource_id=resourceId + "_origin",
        time=source['origin_time'],
        longitude=coord['longitude'],
        latitude=90 - coord['colatitude'],
        depth=Rt - coord['radius'],
    )
    
    mag = Magnitude(
        resource_id=resourceId + "_magnitude",
        # Round to 2 digits.
        mag=round(m_w, 2),
        magnitude_type="mw",
        origin_id=origin.resource_id
    )
    
    foc_mec = FocalMechanism(
        resource_id=resourceId + "_focmec",
        triggering_origin_id=origin.resource_id
    )

    tensor = Tensor(
        m_rr = Mrtp["m_rr"],
        m_pp = Mrtp["m_pp"],
        m_tt = Mrtp["m_tt"],
        m_rt = Mrtp["m_rt"],
        m_rp = Mrtp["m_rp"],
        m_tp = Mrtp["m_tp"]
    )
    
    mt = MomentTensor(
        resource_id=resourceId + "_moment_tensor",
        derived_origin_id=origin.resource_id,
        moment_magnitude_id=mag.resource_id,
        # Convert to Nm.
        scalar_moment = m_0 / 1E7,
        tensor=tensor,
    )
    
    # Assemble everything.
    foc_mec.moment_tensor = mt

    ev = Event(
        resource_id=resourceId + "_event",
        event_type="earthquake"
        )

    ev.origins.append(origin)
    ev.magnitudes.append(mag)
    ev.focal_mechanisms.append(foc_mec)

    # Set the preferred items.
    ev.preferred_origin_id = origin.resource_id.id
    ev.preferred_magnitude_id = mag.resource_id.id
    ev.preferred_focal_mechanism_id = foc_mec.resource_id.id

    ev_file = "event.xml"
    print(f">> Writting event information in {ev_file}")
    ev.write(ev_file, format="QUAKEML")

def main():

    # parsing arguments
    parser = argparse.ArgumentParser()

    parser.add_argument("--file", dest="macromesh_file", type=str, default="macromesh.dat", required=False)
    parser.add_argument("--json", help="store parsong result in a json file", action="store_true", default=False)
    parser.add_argument("--event", help="store event info under obspy CMT solution format", action="store_true", default=False)
    parser.add_argument("--legacy", help="Use the old format (no ulvz specified)", action="store_true", default=False)

    args = parser.parse_args()

    config = parse_macromesh(args.macromesh_file, legacy_fmt = args.legacy)

    pretty_print(config)

    if args.json: write_json(config, args.macromesh_file)

    if args.event: write_quakeML(config)


if __name__ == "__main__":

    main()


