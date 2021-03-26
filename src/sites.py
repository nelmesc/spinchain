#! /usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
import re
import sys
from scipy.optimize import curve_fit
from matplotlib.ticker import MaxNLocator
import math

# Set up the figure
fig = plt.figure(figsize=(10,8))

# Create the axes
ax1 = plt.subplot()

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
    sd = math.sqrt(sum([(x-mean)**2 for x in vals]) / (len(vals)-1))
    err = sd / math.sqrt(len(vals))
    avgs.append(mean)
    errors.append(err)

# Plot each line
plt.plot(numSites, avgs, "g^", ms=10, lw=2, label="data from 100 repeats")
plt.errorbar(numSites, avgs, yerr=errors, fmt="none", color='grey', elinewidth=2,capthick=3,capsize=7)

# Linear regression
start = 0
coef = np.polyfit(numSites[start:],avgs[start:],1)
poly1d_fn = np.poly1d(coef) 
linear = ("y = "  
                + '{0:.2f}'.format(round(coef[0],2)) + r"$x$" + " "
                + '{0:+.2f}'.format(round(coef[1],2)))
# plt.plot(numSites[start:], poly1d_fn(numSites[start:]), '-r', label=linear)

# Quadratic regression
rounding = 3
def f(x, a, b):
    return a*(x**b)
coef, o = curve_fit(f, numSites, avgs)
rsquared = (o[0][1]*o[1][0]) / (o[0][0]*o[1][1])
quadratic = ("y = "  
                + '{0:.3f}'.format(round(coef[0],rounding)) + r"$x^{" + '{0:.3f}'.format(round(coef[1],rounding)) + "}$" + "    $R^2=" + str(round(rsquared, rounding)) + "$")
plt.plot(numSites, f(np.array(numSites), coef[0], coef[1]), '-b', label=quadratic)

# Exponential regression
coef = np.polyfit(numSites,np.log(avgs),1,w=numSites)
y = np.exp(coef[1]) * np.exp(coef[0]*np.array(numSites))
expo = ("y = "  
                + '{0:.2f}'.format(round(np.exp(coef[1]),2)) + " "
                + r"$e^{"+ '{0:.2f}'.format(round(coef[0],2)) + r"x}$")
plt.plot(numSites, y, '-k', label=expo)

# Add inset log plot 
ax2 = ax1.inset_axes([0.20, 0.4, 0.25, 0.25])
ax2.plot(numSites[start:], avgs[start:], "^g")
ax2.plot(numSites[start:], y[start:], "-k")
ax2.set_yscale("log")

# Set the size of the axis ticks
ax1.tick_params(axis='y', pad=10, labelsize=20)
ax1.tick_params(axis='x', pad=10, labelsize=20)

# Set the names of the axis labels
plt.xlabel("", fontsize=20, color='black')
plt.ylabel("", fontsize=20, color='black')
plt.grid(True)
ax1.set_xlabel("spins in linear chain", labelpad=10)
ax1.set_ylabel("time to optimise for 99% fidelity [s]", labelpad=20)
plt.xlim(min(numSites)-1, max(numSites)+1)
ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
plt.legend(prop={'size': 20})

# Nicer spacing
plt.tight_layout()

# Save as a png
plt.savefig('scalingSites.pdf', transparent=False)


