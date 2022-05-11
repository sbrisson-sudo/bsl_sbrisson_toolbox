#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Tue Apr 05 2022
@author: Sylvain Brisson sylvain.brisson@ens.fr

- Simple class to request CMT solution from globalCMT.org and load it into obspy Event object.
- Function to query an event from the jan76_dec20.ndk catalogue by id (brk or ndk)
"""


# HTTP request and HTML parser
from bs4 import BeautifulSoup
import requests

# obspy routine to load events from files
from obspy.core.event import read_events
import obspy

# usual imports
import os
import re
import warnings
import datetime as dt
import matplotlib.pyplot as plt
import linecache


class GlobalCMT_search:
    """Class to make a request to globalcmt.org adn download the data into a obspy Event object.
    """
    
    def __init__(self, **kwargs):
        
        self.date = kwargs.get("date", dt.datetime(1976, 1, 1))
        self.number_days = kwargs.get("nb_days", 1)
        
        self.Mw_min = kwargs.get("Mw_min", 0.0)
        self.Mw_max = kwargs.get("Mw_max", 10.0)
        
        self.lat_min,self.lat_max = kwargs.get("latitude", [-90, 90])
        self.lon_min,self.lon_max = kwargs.get("longitude", [-180, 180])
        self.depth_min = kwargs.get("depth_min", 0)
        self.depth_max = kwargs.get("depth_max", 1000)
        
        self._build_url()
        
    def _build_url(self):
        """Build the request url for the global cmt catalogue
        """
        
        self.url = f"https://www.globalcmt.org/cgi-bin/globalcmt-cgi-bin/CMT5/form?itype=ymd&yr={self.date.year}&mo={self.date.month}&day={self.date.day}&oyr=1976&omo=1&oday=1&jyr=1976&jday=1&ojyr=1976&ojday=1&otype=nd&nday={self.number_days}&lmw={self.Mw_min}&umw={self.Mw_max}&lms=0&ums=10&lmb=0&umb=10&llat={self.lat_min}&ulat={self.lat_max}&llon={self.lon_min}&ulon={self.lon_max}&lhd={self.depth_min}&uhd={self.depth_max}&lts=-9999&uts=9999&lpe1=0&upe1=90&lpe2=0&upe2=90&list=4"
        
    def get_cmt_solution(self):
        """Request the CMT solution and load it into obspy event object.
        """
        
        req = requests.get(self.url)
        
        if req.status_code != 200:
            msg = "Error Downloading file"
            raise Exception(msg)
        
        soup = BeautifulSoup(req.text, "html.parser")
        
        data = soup.find_all("pre")[1].text
        
        if data=="\n\n":
            msg = "Error : no event found"
            raise Exception(msg)
                
        # Create a temorary file to write the CMTSOLUTION data
        cmt_tmp_file = "cmt_solution.tmp"
        with open(cmt_tmp_file, "w") as f:
            f.write(data)
        
        # Read from the previous file into a Event Catalog object
        cat = read_events(cmt_tmp_file, format = "CMTSOLUTION")
        
        # Remove the temporary file
        os.remove(cmt_tmp_file)
        
        # Return only the biggest event found
        # if len(cat) > 1:
        #     warnings.warn(f"{len(cat)} events found, returning the biggest one.", ResourceWarning)
        return cat
    
    # def get_cmt_id(self):

# def getCMT_Id(event):
#     """Return the global cmt id of an obspy event (if in the catalogue)"""
    
#     # resource_id = event.resource_id.id
    
#     # # look for cmt id
#     # isCMT_id = r"/(?P<id>[A-Z]\d{6,12}[A-Z])/"
#     # id_search = re.search(isCMT_id, resource_id)
    
#     # if id_search:
#     #     return  id_search["id"]
    
#     # raise Exception("It seems that the event as not been obtained from globalcmt.")
    
#     # Test get from request
    
#     time = event.preferred_origin().time
    
#     dt = 
    
#     gcmt = GlobalCMT_search(
#         date = dt.datetime(time.year, time.month, time.day),
#         Mw_min = event.preferred_magnitude().mag - 0.1,
#     )
        
#     event = gcmt.get_cmt_solution()
    
def getEventById(id, format="ndk"):
    """Return a obspy event object from a global CMT id"""
        
    if format == "brk":
        stdout = os.popen(f'{os.path.join(os.path.dirname(__file__),"brk2cmt")} {id}')
        id = stdout.readline().strip()
    else:
        if format != "ndk":
            raise Exception("Unknown id format, only ndk and brk supported.")
        id = id.upper()
        
        
    # read identifiers and associated index in catalogue

    catalogue_file = os.path.join(os.path.dirname(__file__),"jan76_dec20.ndk")
    cmt_catalogue = {}

    with open(catalogue_file, "r") as f:

        for i,l in enumerate(f.readlines()[1::5]):
            id_cmt = l.split()[0]
            cmt_catalogue[id_cmt] = i 
            
    try: 
        idx = cmt_catalogue[id]
    except KeyError:
        raise Exception(f"Event {id} not found in globalCMT catalogue {catalogue_file}.")
        
    # read event
    event_ndk = "".join([linecache.getline(catalogue_file, idx*5+i) for i in range(1,6)])
    with open(f"{id_cmt}.ndk", "w") as f: f.write(event_ndk)
    event = read_events(f"{id_cmt}.ndk")[0]
    
    print(event)
    
    return event
        

if __name__ == "__main__":
    
    # Test get from request
    
    gcmt = GlobalCMT_search(
        date = dt.datetime(2014,8,18),
        Mw_min = 6.2,
    )
        
    event = gcmt.get_cmt_solution()
    
    print(event.short_str())
    
    # Test get from id
    
    id = "D14I2NRA"
    event = getEventById(id, format="brk")
    print(event.short_str())
    
    id = "C202012312312A"
    event = getEventById(id, format="ndk")
    print(event.short_str())
    
    


