#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Thu Apr 07 2022
@author: Sylvain Brisson sylvain.brisson@ens.fr

Read two stations list files (outputed by a record section plotting routine)
and return the common stations

"""

import argparse 


parser = argparse.ArgumentParser()

parser.add_argument("in_files", help="2 stations list files", nargs=2, type=str)

args = parser.parse_args()

in_file1,in_file2 = args.in_files

stations1 = open(in_file1, "r").readlines()
stations2 = open(in_file2, "r").readlines()

stations3 = list(set(stations1) & set(stations2))

print(f"{len(stations1)} stations in first file")
print(f"{len(stations2)} stations in second file")
print(f"{len(stations3)} stations in common")

out_file = "common_stations.txt"
open(out_file, "w").writelines(stations3)

print(f"Stations in common outputed in {out_file}")
