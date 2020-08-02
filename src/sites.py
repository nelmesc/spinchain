#! /usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
import re
import sys
from pylab import *

# Set up the figure
plt.figure(figsize=(10,8))

# Create the axes
ax1=plt.subplot()

# Load the data from the sites scaling file
numSites = []
time1Excite = []
time2Excite = []
time3Excite = []
with open("scaling1Excite.dat", "r") as f:
    for line in f:
        split = line.split()
        numSites.append(int(split[0]))
        time1Excite.append(float(split[1]))
with open("scaling2Excite.dat", "r") as f:
    for line in f:
        split = line.split()
        time2Excite.append(float(split[1]))
with open("scaling3Excite.dat", "r") as f:
    for line in f:
        split = line.split()
        time3Excite.append(float(split[1]))

# Plot each line
plt.plot(numSites, time1Excite, "-g^", ms=10, lw=2, label="1 excitation")
plt.plot(numSites, time2Excite, "--b^", ms=10, lw=2, label="2 excitations")
plt.plot(numSites, time3Excite, "-.r^", ms=10, lw=2, label="3 excitations")

# Set the size of the axis ticks
ax1.tick_params(axis='y', pad=10, labelsize=20)
ax1.tick_params(axis='x', pad=10, labelsize=20)

# Set the names of the axis labels
plt.xlabel('number of sites', fontsize=20, color='black')
plt.ylabel('time taken / s', fontsize=20, color='black')
ax1.set_xlabel("number of sites", labelpad=10)
ax1.set_ylabel("time taken / s", labelpad=20)

# Add the legend
l=legend(loc='upper left', ncol=1, fontsize=20)

# Nicer spacing
plt.tight_layout()

# Save as a png
savefig('scalingSites.png', transparent=False)


