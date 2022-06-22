#!/usr/bin/env python3

from parse_macromesh import parse_macromesh
import os


def main():

    config_file = "macromesh.dat"

    if not(os.path.exists(config_file)):
        raise Exception(f"No config file ({config_file}) found, exiting..")
        

    config = parse_macromesh(config_file)
    freq = config["modes"]["database_freq"]
    

    db_dir = config["database_dir"]
    modes_db_dir = os.path.join(db_dir, "sol_ana")
    if not(os.path.exists(modes_db_dir)):
        raise Exception(f"Database directory ({modes_db_dir}) does not exist, exiting..")
    
    db_info_file = os.path.join(modes_db_dir,"info_1")
    if not(os.path.exists(db_info_file)):
        raise Exception("Database directory ({modes_db_dir}) does not contain the metadata file {os.path.basename(db_info_file)}")

    with open(os.path.join(db_info_file), 'r') as f:
        fmin,fmax1,fmax2 = map(float,f.readlines()[5].split())
    
    if fmin!=freq["fmin"] or fmax1!=freq["fmax1"] or fmax2!=freq["fmax2"]:

        print("The database is not compatible with the frequencies asked in the configuration file :")
        freq_db = {'fmin':fmin,'fmax1':fmax1,'fmax2':fmax2}
        print(f"- frequencies of the database : {freq_db}")
        print(f"- frequencies asked :           {freq}")

    else:
        print("The database is compatible with the frequencies asked in the configuration file")

if __name__ == "__main__":

    main()