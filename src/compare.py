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
with open("compare.dat", "r") as f:
    for line in f:
        split = line.split()
        numSites.append(int(split[0]))
        ex1.append(float(split[1]))
        ex2.append(float(split[2]))

# Plot each line
plt.plot(numSites, ex1, "g^", ms=10, lw=2, label="genetic")
plt.plot(numSites, ex2, "bo", ms=10, lw=2, label="PST")

# Set the size of the axis ticks
ax1.tick_params(axis='y', pad=10, labelsize=20)
ax1.tick_params(axis='x', pad=10, labelsize=20)

# Force axes limits
ax1.set_xlim([0,7])
ax1.set_ylim([0.3,1.2])

# Set the names of the axis labels
plt.xlabel("", fontsize=20, color='black')
plt.ylabel("", fontsize=20, color='black')
plt.grid(True)
ax1.set_xlabel("site number, i", labelpad=10)
ax1.set_ylabel(r"$J_{i,i+1}~/~J_{max}$", labelpad=20)
plt.xlim(min(numSites)-1, max(numSites)+1)
ax1.xaxis.set_major_locator(MaxNLocator(integer=True))

plt.legend(prop={'size': 20}, loc="lower right")

# Nicer spacing
plt.tight_layout()

# Save as a png
savefig('compare.pdf', transparent=False)


