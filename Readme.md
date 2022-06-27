### SBRISSON BSL toolbox

For question you can contact me at sylvain.brisson@ens.fr

This repository contains scripts (mainly written in python) I developped during my internship at the Berkeley Seismology Lab working on ULVZ prospection and modelisation, with prof. Barbara Romanowicz.

The dependencies are:
- global python packages (available with conda or pip) and listed in the `environment.yml` file
- bsl python package `ucbpy` and bsl models files, to access these you will have to change the paths in `common/setup.py`

For documentation on the use and the options of these routines use the `-h` option.

Inputs files :
- stations are given under the format used by cSEM (`receivers.dat`)
- events are given under the QuakeML format
- some metadata are sometimes read from the cSEM input file (`macromesh.dat`)
- waveform data are read from a obspy `stream` object (serialized into a file : a pickle file)
 

These tools are distributed in different modules :
1. `csem_tools` : **tools to work with sandwich csem outputs and inputs**
   - `check_database_compatibility.py` : check that the DTN database is compatible with the frequencies written in the configuration file `macromesh.dat`
   - `convert_to_pickle.py` : convert the multiple ascii traces files returned by CSEM to a single serialized obspy stream object (to be loaded and used by obspy routines)
   - `generate_receivers.py` : generate a `receivers.dat` file with receivers sampling at regulat distance and azimuth position
   - `parse_macromesh.py` : parse the configuration file `macromesh.dat`, used by the other routines + can export the information under a json format, and the event information under a format compatible with an obspy-based workflow. 
   - `plot_trace.py` : plot one or several ascii traces files.
2. `data_acquisition` : **tools to get observation data**, the main workflow which used these scripts can be found on another repository at ???
   - `brkcmt` : convert from the bsl event id format to the cmt format
   - `distance_azimuth_domain.py` : define an obspy domain object delimited by distances bounds from the source and azimuth bounds (in reference to an arbitraty point).
   - `globalcmt_request.py` : get an event catalogue from the global CMT project from date,magnitude,depth bounds.
   - `get_event_globalCMT.py` : return an obspy-compatible event file from an CMT or BRK event id, calls `globalcmt_request.py`.
   - `mdl_plotting.py` : plot the domain defined in `distance_azimuth_domain.py`.
   - `mseed2ascii.py` : convert trace from mseed to 3 ascii files.

3. `plotting` : **plotting scripts** : for CSEM data but also BSL models
    - `A3Dmodel_map_greatCircles.py` :  plot a A3d model as a map (with an orthographic projection). If provided the code also plot stations and event along with the great circles between them.
    - `A3Dmodel_section_rayPaths.py` : plot a A3d model as a mantle+crust section, possibly anlong with the source and receivers position and the path between them (calculated with taup in prem).
    - `max_SH_CMT.py` : (not very usefull) for a source position and an azimuth and takeoff angle, find the moment tensor that maximize the SH radiation pattern in this direction
    - `plot_event.py` : obspy plotting event 
    - `plot_green.py` : plot the green functions (output by `hsemm`)
    - `radiation_pattern.py` : plot the horizontal and vertical radiation patterns for P,SV and SH
    - `radiation_pattern_3d.py` : plot the 3d radiation pattern
    - `radiation_pattern_3d_obspy.py` : obspy weird radiation plotting routine
    - `ray_path.py` : plot the path followed by the rays
    - `station_source.py` : plot the statoin and source position
    - `waveforms_recordSection_azimuth.py` : record section of traces along the azimuth (consider using the gui software)
    - `waveforms_recordSection_distance.py` : record section of traces along the distance (consider using the gui software)
    - `waveforms_traces.py` : plot ascii trace files

4. `postprocessing` : **miscellanous computation scripts** 
   - `filter.py` : band filter a obspy stream object
   - `spectrum.py` : compute and plot the frequencie content of traces

