## sandwich cSEM utilitaries

1. Plotting traces : `plotTraces.py`

    ```
    Usage : plotTraces.py <traces files>
    Options:
    -o <output figure file name>
    -t <t1> <t2> : time bounds
    -h : help
    Exemple : ./plotTraces.py U*_TA_W18A
    ```

2. Getting a summary of the macromesh.dat information : `parseMacromesh.py`

    ```
    Usage : parseMacromesh.py <macromesh.dat file>
    ```

3. Getting station position and distance to source : `getStationInfo.py`

    ```
    Usage : getStationInfo.py <station file> <station names>
    Options:
    -d lat lon : source
    -h : help 
    ```

