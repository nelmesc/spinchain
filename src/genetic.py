#! /usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
import re
import sys
from pylab import *

# Set up the figure
plt.figure(figsize=(10,7))

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

plt.plot(generation,worst,color='red',lw=2,label='worst')
plt.plot(generation,average,color='orange',lw=2,label='average')
plt.plot(generation,best,color='green',lw=2,label='best')

#SET SIZE OF AXIS TICKS
ax1.tick_params(axis='y', labelsize=20)
ax1.tick_params(axis='x', labelsize=20)

#SET LABEL NAMES
plt.xlabel('generation',fontsize=25,color='black')
plt.ylabel('fitness',fontsize=25,color='black')

#LEGEND
l=legend(bbox_to_anchor=(0, 1), loc='upper left', ncol=1,frameon=False,borderaxespad=0.,fontsize=20)

#SAVE AS PNG PICTURE
savefig('genetic.png',transparent=False)
plt.show()
