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

# Load the data from the parallel scaling file
cores = []
speedup = []
ideal = []
with open("scalingParallel.dat", "r") as f:
    for line in f:
        split = line.split()
        cores.append(int(split[0]))
        speedup.append(float(split[1]))
        ideal.append(int(split[0]))

# Plot each line
plt.plot(cores, speedup, "-g^", lw=2, ms=10, label="actual")
plt.plot(cores, ideal, "--b", lw=2, ms=10, label="ideal")

# Set the size of the axis ticks
ax1.tick_params(axis='y', labelsize=20)
ax1.tick_params(axis='x', labelsize=20)

# Set the names of the axis labels
plt.xlabel('cores used', fontsize=20, color='black')
plt.ylabel('speedup', fontsize=20, color='black')
ax1.set_xlabel("cores used", labelpad=10)
ax1.set_ylabel("speedup", labelpad=20)

# Add the legend
l=legend(loc="upper left", ncol=1, fontsize=20)

# Nicer spacing
plt.tight_layout()

# Save as a png
savefig('scalingParallel.png', transparent=False)

