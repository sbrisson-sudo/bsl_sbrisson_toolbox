#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Tue Apr 05 2022
@author: Sylvain Brisson sylvain.brisson@ens.fr

Get event from the global cmt catalogue (from id or from a day (+min Mw))
"""

# importations

import os
import argparse

import datetime as dt

from globalcmt_request import getEventById, GlobalCMT_search


parser = argparse.ArgumentParser()

# subparsers = parser.add_subparsers(help='sub-command help', required =True)

# parser_id = subparsers.add_parser('by_id', help='get by id')
# parser_id.add_argument('id', type=str, help='global CMT event ID')
# parser_id.add_argument('--brk', action="store_true", help='id is under brk format')
# parser_id.set_defaults(function=getById)


# parser_day = subparsers.add_parser('by_day', help='get by day and min magnitude')
# parser_day.add_argument("dt", help="YYYYMMDD", type=str)
# parser_day.add_argument("--Mw", help="minimal magnitude", type=float, default=5.5)
# parser_day.set_defaults(function=getBydt)


parser.add_argument('--id', type=str, help='global CMT event ID')
parser.add_argument('--brk', action="store_true", help='id is under brk format')


parser.add_argument("--dt", help="YYYYMMDD", type=str)
parser.add_argument("--Mw", help="minimal magnitude", type=float, default=5.5)

args = parser.parse_args()

if args.id:
    
    event = getEventById(args.id, format="brk" if args.brk else "ndk")

else:
    
    gcmt = GlobalCMT_search(
        date = dt.datetime.strptime(args.dt, "%Y%m%d"),
        Mw_min = 6.2,
    )
       
    events = gcmt.get_cmt_solution()
    
    if len(events) > 1:
        raise Exception("More than one event bigger than Mw={args.Mw} on {args.dt}.")

    event = events[0]
    

event_id = f"{event.magnitudes[0].mag:1.1f}_{event.origins[0].time.date.strftime('%d-%b-%Y')}"

print(event.short_str())

print(f">> Writting {event_id}.xml")
event.write(f"{event_id}.xml", format="quakeml") 
    

    


    
