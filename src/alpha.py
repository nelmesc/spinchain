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
alphas = []
fidsBefore = []
fidsAfter = []
with open("eetest.dat", "r") as f:
    for line in f:
        split = line.split()
        alphas.append(float(split[0]))
        fidsBefore.append(float(split[1][:-1]) / 100.0)
        # fidsAfter.append(float(split[2][:-1]))

# Plot each line
plt.plot(alphas, fidsBefore, "-gx", ms=10, lw=2, label="before re-optimisation")
# plt.plot(alphas, fidsAfter, "--bx", ms=10, lw=2, label="after re-optimisation")

# Set the size of the axis ticks
ax1.tick_params(axis='y', pad=10, labelsize=20)
ax1.tick_params(axis='x', pad=10, labelsize=20)

# Set the names of the axis labels
plt.xlabel(r'$\alpha$', fontsize=20, color='black')
plt.ylabel(r"$F(t_f)$", fontsize=20, color='black')
ax1.set_xlabel(r"$\alpha$", labelpad=10)
ax1.set_ylabel(r"$F(t_f)$", labelpad=20)

# Add the legend
# l=legend(loc='lower right', ncol=1, fontsize=20)

# Nicer spacing
plt.tight_layout()

# Save as a png
savefig('alphaScaling.png', transparent=False)


