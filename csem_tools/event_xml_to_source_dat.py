#!/usr/bin/env python3

"""

Read a obspy compatible event file and write the 
information needed by csem in a source.dat file

"""

from obspy import read_events
import argparse
import numpy as np


R_EARTH_M  = 6371000.
NM2CSEM    = 1e-20
DYNCM2CSEM = 1e-27

def write_source_dat(event, out_file = "source.dat"):

    tensor = event.preferred_focal_mechanism().moment_tensor.tensor
    ev_lat = event.preferred_origin().latitude
    ev_lon = event.preferred_origin().longitude
    ev_dep = event.preferred_origin().depth
    origin_time = event.preferred_origin().time

    ev_colat = ev_lat + 90.
    ev_radius = R_EARTH_M - ev_dep

    # compute magnitude to check units
    m_0 = 1.0 / np.sqrt(2.0) * np.sqrt(
                (tensor.m_rr) ** 2 +
                (tensor.m_tt) ** 2 +
                (tensor.m_pp) ** 2 +
        2.0 *   (tensor.m_rt) ** 2 +
        2.0 *   (tensor.m_rp) ** 2 +
        2.0 *   (tensor.m_tp) ** 2)
    m_w = 2.0 / 3.0 * (np.log10(m_0) - 16.1)
    print(f"Mw = {m_w:.2f}")

    with open( out_file, 'w') as out:

        out.write(
        f"""# Position (radius(m), colatitude(°), longitude(°))
{ev_radius}
{ev_colat:8.4f}
{ev_lon:8.4f}
# Centroid Moment Tensor (10-20N.m) (m_rr,m_tt,m_pp,m_rt,m_rp,m_tp)
{tensor.m_rr * NM2CSEM:.5f}
{tensor.m_tt * NM2CSEM:.5f}
{tensor.m_pp * NM2CSEM:.5f}
{tensor.m_rt * NM2CSEM:.5f}
{tensor.m_rp * NM2CSEM:.5f}
{tensor.m_tp * NM2CSEM:.5f}
# Centroid Moment Tensor (10-27dyn.cm) (m_rr,m_tt,m_pp,m_rt,m_rp,m_tp)
{tensor.m_rr * DYNCM2CSEM:.5f}
{tensor.m_tt * DYNCM2CSEM:.5f}
{tensor.m_pp * DYNCM2CSEM:.5f}
{tensor.m_rt * DYNCM2CSEM:.5f}
{tensor.m_rp * DYNCM2CSEM:.5f}
{tensor.m_tp * DYNCM2CSEM:.5f}
"""
        )

    print(f"Event information for CSEM written in {out_file}")

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("event_file", type=str)
    parser.add_argument("-o", dest="out_file", help="output file name", type=str, default="source.dat")

    args = parser.parse_args()

    # read event
    event = read_events(args.event_file)[0]

    # write info
    write_source_dat(event, out_file=args.out_file)

