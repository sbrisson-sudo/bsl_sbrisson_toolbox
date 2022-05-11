This directory contains plotting routines for plotting waveform data and UC Berkeley A3d models along with metadata (stations position, event information) and also plotting routines related to the sandwich cSEM codes.

For documentation on the use and the options of these routines use the `-h` option.

Inputs files :
- stations are given under the format used by cSEM (`receivers.dat`)
- events are given under the QuakeML format
- some metadata are sometimes read from the cSEM input file (`macromesh.dat`)
- waveform data are read from a obspy `stream` object (serialized into a file : a pickle file)
 

The data directory contains additionnal data for plots (hotspots and plates boundaries).

These are what they do:
- `plot_event.py` : plot and event (from a QuakeML file) : plot the beachball on an orthographic projection
- `plot_green.py` : plot the green functions produced by the `hsemm` executable
- `plot_path_map.py` : plot a A3d model as a map (with an orthographic projection). If provided the code also plot stations and event along with the great circles between them.
- `plot_path_section.py` : plot a A3d model as a mantle+crust section, possibly anlong with the source and receivers position and the path between them (calculated with taup in prem).
- `plot_pickle_azimuth.py` : plot a waveform record section along the azimuth
- `plot_pickle_distance.py` : plot a waveform record section along the distance
- `plot_pickle` ???
- `plot_ray_path.py` : plot path in prem from a source to a list of receivers for different phases
- `plot_station_source.py` : basic scatter plot of the source and the stations position
- `plot_traces.py` : plot waveforms (3 axes for R T and Z), possibly from differents stations/different methods