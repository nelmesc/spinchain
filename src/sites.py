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
times = []
with open("scaling1Excite.dat", "r") as f:
    for line in f:
        split = line.split()
        numSites.append(int(split[0]))
        times.append([])
        for val in split[1:]:
            times[-1].append(float(val))

# Calculate averages
avgs = []
errors = []
for vals in times:
    mean = sum(vals) / len(vals)
    sd = math.sqrt(sum([(x-mean)**2 for x in vals]))
    err = sd / sqrt(len(vals))
    avgs.append(mean)
    errors.append(err)

# Plot each line
plt.plot(numSites, avgs, "-g^", ms=10, lw=2, label="1 excitation (left end to right end)")
plt.errorbar(numSites, avgs, yerr=errors, fmt="none", color='grey', elinewidth=2,capthick=3,capsize=7)

# Line of best fit
# coef = np.polyfit(numSites,avgs,1,w=[1/a for a in errors])
# coef = np.polyfit(numSites,avgs,1)
coef = np.polyfit(numSites,avgs,1,w=[i for i, a in enumerate(errors)])
poly1d_fn = np.poly1d(coef) 
plt.plot(numSites[4:], poly1d_fn(numSites[4:]), '--r')

# Set the size of the axis ticks
ax1.tick_params(axis='y', pad=10, labelsize=20)
ax1.tick_params(axis='x', pad=10, labelsize=20)

# Set the names of the axis labels
plt.xlabel("", fontsize=20, color='black')
plt.ylabel("", fontsize=20, color='black')
plt.grid(True)
ax1.set_xlabel("spins in linear chain", labelpad=10)
ax1.set_ylabel("time to optimise for 99% fidelity (s)", labelpad=20)
plt.xlim(min(numSites)-1, max(numSites)+1)
ax1.xaxis.set_major_locator(MaxNLocator(integer=True))

# Nicer spacing
plt.tight_layout()

# Save as a png
savefig('scalingSites.png', transparent=False)


