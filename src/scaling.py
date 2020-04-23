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

# Load the data from the scaling file
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
plt.plot(cores, speedup, color='green', lw=2, label="actual")
plt.plot(cores, ideal, color='blue', lw=2, label="ideal")

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

plt.tight_layout()

# Save as a png
savefig('scalingParallel.png', transparent=False)

# Clear to reuse for the site scaling
plt.clf()

# Set up the figure
plt.figure(figsize=(10,8))

# Create the axes
ax1=plt.subplot()

# Load the data from the scaling file
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
plt.plot(numSites, time1Excite, color='green', lw=2, label="1 excitation")
plt.plot(numSites, time2Excite, color='blue', lw=2, label="2 excitation")
plt.plot(numSites, time3Excite, color='red', lw=2, label="3 excitation")

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

plt.tight_layout()

# Save as a png
savefig('scalingSites.png', transparent=False)


