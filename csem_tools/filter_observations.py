#!/usr/bin/env python


"""

Read observations (as MSEED traces with associated XML station metadata) 
and filtered it with a CSEM source time function (defined by 4 frequencies
specified with -freq option)

Save the result in an obspy stream object.

IMPORTANT. It also apply a highpass filter at 20s (see in the code).

"""

import numpy as np
from scipy import fftpack

import os

import argparse 

from obspy import Stream, Trace, read, read_inventory
from obspy.geodetics import gps2dist_azimuth,locations2degrees
from obspy.core.event import read_events

from tqdm import tqdm

# CONFIGURATION

def heavis_filter(f1,f2,f3,f4,f):

    if f < f1 or f > f4 : return 0.

    if f < f2 : 
        return (-np.cos((f-f1)/(f2-f1)*np.pi) + 1.)/2
        
    if f > f3 : 
        return (np.cos((f-f3)/(f4-f3)*np.pi) + 1.)/2

    return 1.

def norm(x): return x/x.max()

def filter_stream(stream, f1,f2,f3,f4):

    """
    Apply the heaviside filter to a stream object
    """

    traces_filt = []

    for tr in tqdm(stream, total=len(stream)):

        # get data
        sig = tr.data
        N   = tr.stats.npts
        dt  = 1/tr.stats.sampling_rate # given in microsecond in mseed

        # tapering (hanning)
        window = np.hanning(N)
        sig_tappered = window*sig

        # fft (scipy)
        sig_fft = fftpack.fft(sig_tappered)
        freq = fftpack.fftfreq(sig_tappered.size, d=dt)

        # filtering 
        heavis_filt = np.array([heavis_filter(f1,f2,f3,f4,f) for f in freq])
        sig_fft_filt = sig_fft.copy()
        sig_fft_filt *= heavis_filt

        # ifft (scipy)
        sig_filt = np.real(fftpack.ifft(sig_fft_filt))

        # plt.plot(sig)
        # plt.plot(sig_filt)
        # plt.show()
        # exit()

        # writting back the trace 
        traces_filt.append(Trace(data = sig_filt, header = tr.stats))

    return Stream(traces = traces_filt)

def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("wf_st_dir", help="directories containing the waveforms then directrory containing the stations",type=str,nargs="+")

    parser.add_argument("--event", help="event_file",type=str,required=True)

    parser.add_argument("-freq", help="4 corners frequencies of the heaviside filter",type=float,nargs=4,required=True)
    parser.add_argument("-o", dest="out_file", help="output file prefix for serialized obspy file", type=str, required=True)
    parser.add_argument("-d", dest="dist_range", help="distance bounds", nargs=2, type=float, default=[None,None])

    
    args = parser.parse_args()

    # get directrories
    wf_dir,st_dir = args.wf_st_dir

    # get filter frequencies

    f1,f2,f3,f4 = args.freq
    # f1,f2,f3,f4 = 1.000e-03, 1.000e-02, 8.333e-02, 1.000e-01 

    print(">> Loading data")
    # load waveform data 
    st = Stream()
    for waveform in os.listdir(wf_dir):
        st += read(os.path.join(wf_dir, waveform))

    # loading station metadata
    stations = read_inventory(f"{st_dir}/*")

    # load event info
    event = read_events(args.event)[0]
    ev_lat = event.preferred_origin().latitude
    ev_lon = event.preferred_origin().longitude
    ev_dep = event.preferred_origin().depth
    origin_time = event.preferred_origin().time

    # computing aditionnal metadata
    for tr in st:
        
        st_coordinates = stations.get_coordinates(tr.id)
        st_lat = st_coordinates["latitude"]
        st_lon = st_coordinates["longitude"]
        
        tr.stats.coordinates = st_coordinates
        
        tr.stats.evla = ev_lat
        tr.stats.evlo = ev_lon
        tr.stats.evde = ev_dep/1000
        tr.stats.event_origin_time = origin_time
                
        _,b_az,_ = gps2dist_azimuth(st_lat, st_lon, ev_lat, ev_lon)
        dist = locations2degrees(st_lat, st_lon, ev_lat, ev_lon)
        
        tr.stats.back_azimuth = b_az # for rotation
        tr.stats.distance = dist

    # select on distance
    def select_distance(stream, d1, d2):
        """Select traces in stream based on distance to event."""
        to_remove = []
        for trace in stream.traces:
            if not(d1 <= trace.stats.distance and d2 >= trace.stats.distance):
                to_remove.append(trace)
        for trace in to_remove:
            stream.remove(trace)
        return stream

    nb_traces = len(st)

    if args.dist_range[0] or args.dist_range[0]:
        d1, d2 = args.dist_range
        st = select_distance(st, d1, d2)

    print(f"{nb_traces - len(st)} removed by distance filter")

    # interpolating
    print(">> Interpolating data")
    st.interpolate(sampling_rate=5.0)

    # filtering data       
    print(">> Applying heaviside filter")
    st = filter_stream(st, f1, f2, f3, f4)

    # second highpass filter
    print(">> Applying highpass filter 20s")
    st.filter('highpass', freq=1/20.)

    # rotating it
    stations = set([tr.stats.station for tr in st])
    st._trim_common_channels()

    for station in stations:                       
        try:
            st.select(station=station).rotate('NE->RT')
        except ValueError:
            print(f"Couldn't rotate:\n{st.select(station=station)}")

    # saving data
    out_file = args.out_file + ".pkl"
    if os.path.exists(out_file):
        if input(f"{out_file} exists, overiding ? [y/n] ") != "y":
            print("Exiting, not saved.")
            
    print(f">> Writting {out_file}")
    st.write(out_file, format='PICKLE')
    
        
            
if __name__ == "__main__":
    
    main()