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

rounding = 5

# Linear regression
# start = 10
# def f1(x, a, b):
    # return a*x+b
# coef, o = curve_fit(f1, numSites[start:], avgs[start:])
# rsquared = 1 - (np.sum((np.array(avgs[start:]) - f1(np.array(numSites[start:]), *coef))**2) / np.sum((avgs[start:] - np.mean(avgs[start:]))**2))
# linear = ("y = " + '{0:.2f}'.format(round(coef[0],2)) + r"$x$" + " " + '{0:+.2f}'.format(round(coef[1],2)) + "    $R^2=" + str(round(rsquared, rounding)) + "$")
# plt.plot(numSites[start:], f1(np.array(numSites[start:]), coef[0], coef[1]), '-r', label=linear)

# Quadratic regression 1
start = 0
end = 7
def f2(x, a, b):
    return a*(x**b)
coef, o = curve_fit(f2, numSites[start:end], avgs[start:end])
rsquared = 1 - (np.sum((np.array(avgs[start:end]) - f2(np.array(numSites[start:end]), *coef))**2) / np.sum((avgs[start:end] - np.mean(avgs[start:end]))**2))
quadratic = ("y = " + '{0:.5f}'.format(round(coef[0],rounding)) + r"$x^{" + '{0:.5f}'.format(round(coef[1],rounding)) + "}$" + "    $R^2=" + str(round(rsquared, rounding)) + "$")
plt.plot(numSites[start:end], f2(np.array(numSites[start:end]), coef[0], coef[1]), '-b', label=quadratic)

# Quadratic regression 2
start = 7
end = len(numSites)
def f2(x, a, b):
    return a*(x**b)
coef, o = curve_fit(f2, numSites[start:end], avgs[start:end])
rsquared = 1 - (np.sum((np.array(avgs[start:end]) - f2(np.array(numSites[start:end]), *coef))**2) / np.sum((avgs[start:end] - np.mean(avgs[start:end]))**2))
quadratic = ("y = " + '{0:.5f}'.format(round(coef[0],rounding)) + r"$x^{" + '{0:.5f}'.format(round(coef[1],rounding)) + "}$" + "    $R^2=" + str(round(rsquared, rounding)) + "$")
plt.plot(numSites[start:end], f2(np.array(numSites[start:end]), coef[0], coef[1]), '-r', label=quadratic)

# Exponential regression
# start = 7
# end = len(numSites)
# def f3(x, a, b):
    # return b*np.exp(a*x)
# coef, o = curve_fit(f3, numSites[start:end], avgs[start:end])
# rsquared = 1 - (np.sum((np.array(avgs[start:end]) - f3(np.array(numSites[start:end]), *coef))**2) / np.sum((avgs[start:end] - np.mean(avgs[start:end]))**2))
# expo = ("y = " + '{0:.2f}'.format(round(np.exp(coef[1]),2)) + " " + r"$e^{"+ '{0:.2f}'.format(round(coef[0],2)) + r"x}$" + "    $R^2=" + str(round(rsquared, rounding)) + "$")
# plt.plot(numSites[start:end], f3(np.array(numSites[start:end]), coef[0], coef[1]), '-g', label=expo)

# Add inset log plot 
# ax2 = ax1.inset_axes([0.20, 0.4, 0.25, 0.25])
# ax2.plot(numSites[start:], avgs[start:], "^g")
# ax2.plot(numSites[start:], y[start:], "-k")
# ax2.set_yscale("log")

# Set the size of the axis ticks
ax1.tick_params(axis='y', pad=10, labelsize=20)
ax1.tick_params(axis='x', pad=10, labelsize=20)

# Set the names of the axis labels
plt.xlabel("", fontsize=20, color='black')
plt.ylabel("", fontsize=20, color='black')
plt.grid(True)
ax1.set_xlabel("spins in linear chain", labelpad=10)
ax1.set_ylabel("time to optimise for 95% fidelity [s]", labelpad=20)
plt.xlim(min(numSites)-1, max(numSites)+1)
ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
plt.legend(prop={'size': 20})

# Nicer spacing
plt.tight_layout()

# Save as a png
plt.savefig('scalingSites.pdf', transparent=False)
ax1.set_xscale("log")
ax1.set_yscale("log")
plt.savefig('scalingSitesLog.pdf', transparent=False)


