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
ex1 = []
ex2 = []
ex3 = []
with open("scalingFixed.dat", "r") as f:
    for line in f:
        split = line.split()
        numSites.append(int(split[0]))
        ex1.append(float(split[1]))
        ex2.append(float(split[2]))
        ex3.append(float(split[3]))

# Plot each line
plt.plot(numSites, ex1, "-g^", ms=10, lw=2, label="1 excitation")
plt.plot(numSites, ex2, "-b^", ms=10, lw=2, label="2 excitations")
plt.plot(numSites, ex3, "-r^", ms=10, lw=2, label="3 excitations")

# Set the size of the axis ticks
ax1.tick_params(axis='y', pad=10, labelsize=20)
ax1.tick_params(axis='x', pad=10, labelsize=20)

# Set the names of the axis labels
plt.xlabel("", fontsize=20, color='black')
plt.ylabel("", fontsize=20, color='black')
plt.grid(True)
ax1.set_xlabel("spins in linear chain", labelpad=10)
ax1.set_ylabel("time for 200 generations [s]", labelpad=20)
plt.xlim(min(numSites)-1, max(numSites)+1)
ax1.xaxis.set_major_locator(MaxNLocator(integer=True))

plt.legend(prop={'size': 20})

# Nicer spacing
plt.tight_layout()

# Save as a png
savefig('scalingFixed.pdf', transparent=False)


