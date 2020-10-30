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

# Load the fitnesses from the genetic file
generation = []
worst = []
average = []
best = []
with open(sys.argv[1], "r") as f:
    for line in f:
        split = line.split()
        if len(split) > 5 and split[0] == "|" and split[2] == "|" and not split[1] == "Generation":
            generation.append(int(split[1]))
            worst.append(float(split[3]))
            average.append(float(split[5]))
            best.append(float(split[7]))

# Define the ranges
ax1.set_xlim([1,generation[-1]])
ax1.set_ylim([0,100])

# Plot each line
plt.plot(generation, best, "-", color='green', lw=2, label='best')
plt.plot(generation, average, "-", color='orange', lw=2, label='average')
plt.plot(generation, worst, "-", color='red', lw=2, label='worst')

# Set the size of the axis ticks
ax1.tick_params(axis='y', labelsize=20)
ax1.tick_params(axis='x', labelsize=20)

# Set the names of the axis labels
plt.xlabel('generation', fontsize=20, color='black')
plt.ylabel('fitness', fontsize=20, color='black')
plt.grid(True)
plt.tight_layout()

# Add the legend
l=legend(loc='lower right', ncol=1, fontsize=20)

# Save as a png
savefig('genetic.png', transparent=False)
