#!/usr/bin/env python

"""
Plot different spectrums with associated labels, computed by the mean_spectrum routines
"""


import numpy as np 
import matplotlib.pyplot as plt
import sys


if __name__ == "__main__":

    in_files = sys.argv[1::2]
    labels = sys.argv[2::2]


    for in_file,label in zip(in_files, labels):

        data = np.loadtxt(in_file)

        w = data[:,0]
        h = data[:,1]

        # norm it
        h /= h.max()

        plt.plot(w, h, label=label)

    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Power")
    plt.xlim([0.0, 0.15])
    plt.grid()
    plt.legend()

    plt.show()
