#!/usr/bin/env python

from matplotlib import pyplot as plt
import numpy as np
import os

if len(os.sys.argv) != 2:
    print("Missing argument: Data file.")
    exit(1)

filename = os.sys.argv[1]
histname = filename.split('.')[0] + ".hist"

data = np.genfromtxt(filename).T
hist, bins = np.histogram(data[1], bins=25, weights=data[0])

with open(histname, "w") as f:
    f.write("# Bins Hist\n")
    for row in zip(bins, hist):
        f.write("{:7.3f} {:7.3f}\n".format(row[0], row[1]))
