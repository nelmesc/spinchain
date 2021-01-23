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
    sd = math.sqrt(sum([(x-mean)**2 for x in vals]) / (len(vals)-1))
    err = sd / sqrt(len(vals))
    avgs.append(mean)
    errors.append(err)

# Plot each line
plt.plot(numSites, avgs, "g^", ms=10, lw=2, label="data from 100 repeats")
plt.errorbar(numSites, avgs, yerr=errors, fmt="none", color='grey', elinewidth=2,capthick=3,capsize=7)

rounding = 3

# Linear regression
start = 10
coef = np.polyfit(numSites[start:],avgs[start:],1)
poly1d_fn = np.poly1d(coef) 
linear = ("y = "  
                + '{0:.2f}'.format(round(coef[0],2)) + r"$x$" + " "
                + '{0:+.2f}'.format(round(coef[1],2)))
plt.plot(numSites[start:], poly1d_fn(numSites[start:]), '--r', label=linear)

# Quadratic regression
# from scipy.optimize import curve_fit
# def f(x, a, b):
    # return a*np.exp(b*x)
# coef, o = curve_fit(f, numSites[0:10], avgs[0:10])
# quadratic = ("y = "  
                # + '{0:.3f}'.format(round(coef[0],rounding)) + r"$x^2$" + " "
                # + '{0:+.3f}'.format(round(coef[1],rounding)))
# plt.plot(numSites[0:10], f(np.array(numSites[0:10]), coef[0], coef[1]), '--b', label=quadratic)
# def f(x, a, b):
    # return a*(x**2)+b*x
# coef, o = curve_fit(f, numSites, avgs)
# quadratic = ("y = "  
                # + '{0:.3f}'.format(round(coef[0],rounding)) + r"$x^2$" + " "
                # + '{0:+.3f}'.format(round(coef[1],rounding)) + r"$x$")
# plt.plot(numSites, f(np.array(numSites), coef[0], coef[1]), '--r', label=quadratic)
# def f(x, a):
    # return a*(x**2)
# coef, o = curve_fit(f, numSites, avgs)
# quadratic = ("y = "  
                # + '{0:.3f}'.format(round(coef[0],rounding)) + r"$x^2$")
# plt.plot(numSites, f(np.array(numSites), coef[0]), '--k', label=quadratic)

# Exponential regression
# coef = np.polyfit(numSites,np.log(avgs),1,w=numSites)
# y = np.exp(coef[1]) * np.exp(coef[0]*np.array(numSites))
# expo = ("y = "  
                # + '{0:.2f}'.format(round(np.exp(coef[1]),2)) + " "
                # + r"$e^{"+ '{0:.2f}'.format(round(coef[0],2)) + r"}$")
# plt.plot(numSites, y, '--k', label=expo)

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
savefig('scalingSites.pdf', transparent=False)


