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
dataSet = -1
cores = []
speedup = []
ideal = []
with open("scalingParallel.dat", "r") as f:
    for line in f:
        split = line.split()
        if split[0] == "1":
            dataSet += 1
            cores.append([])
            speedup.append([])
        cores[dataSet].append(int(split[0]))
        speedup[dataSet].append(float(split[1]))

for core in cores[dataSet]:
    ideal.append(int(core))

# Plot each line
for i in range(dataSet+1):

    customLabel = ""
    lineType = ":r^"
    if i == 0:
        lineType = ":r^"
        customLabel = "7 spins, 1 excitation"
    elif i == 1:
        lineType = "--go"
        customLabel = "7 spins, 3 excitations"
    elif i == 2:
        lineType = "-.bx"
        customLabel = "16 spins, 2 excitations"

    if len(customLabel) == 0:
        customLabel = "set " + str(i)

    plt.plot(cores[i], speedup[i], lineType, lw=2, ms=10, label=customLabel)

plt.plot(cores[0], ideal, "-k", lw=2, ms=10, label="ideal")

# Set the size of the axis ticks
ax1.tick_params(axis='y', labelsize=20)
ax1.tick_params(axis='x', labelsize=20)

# Set the names of the axis labels
plt.xlabel('cores used', fontsize=20, color='black')
plt.ylabel('speedup', fontsize=20, color='black')
plt.grid(True)
ax1.set_xlabel("cores used", labelpad=10)
ax1.set_ylabel("speedup", labelpad=20)

# Add the legend
l=legend(loc="upper left", ncol=1, fontsize=20)

# Nicer spacing
plt.tight_layout()

# Save as a png
savefig('scalingParallel.pdf', transparent=False)

