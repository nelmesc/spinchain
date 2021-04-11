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
times90 = []
with open("scalingSites90.dat", "r") as f:
    for line in f:
        split = line.split()
        numSites.append(int(split[0]))
        times90.append([])
        for val in split[1:]:
            times90[-1].append(float(val))
times95 = []
with open("scalingSites95.dat", "r") as f:
    for line in f:
        split = line.split()
        times95.append([])
        for val in split[1:]:
            times95[-1].append(float(val))
times97 = []
with open("scalingSites97.dat", "r") as f:
    for line in f:
        split = line.split()
        times97.append([])
        for val in split[1:]:
            times97[-1].append(float(val))

# Calculate averages
avgs90 = []
errors90 = []
for vals in times90:
    mean = sum(vals) / len(vals)
    sd = math.sqrt(sum([(x-mean)**2 for x in vals]) / (len(vals)-1))
    err = sd / math.sqrt(len(vals))
    avgs90.append(mean)
    errors90.append(err)
avgs95 = []
errors95 = []
for vals in times95:
    mean = sum(vals) / len(vals)
    sd = math.sqrt(sum([(x-mean)**2 for x in vals]) / (len(vals)-1))
    err = sd / math.sqrt(len(vals))
    avgs95.append(mean)
    errors95.append(err)
avgs97 = []
errors97 = []
for vals in times97:
    mean = sum(vals) / len(vals)
    sd = math.sqrt(sum([(x-mean)**2 for x in vals]) / (len(vals)-1))
    err = sd / math.sqrt(len(vals))
    avgs97.append(mean)
    errors97.append(err)

# Params
start = 13
fittingStart = 13
rounding = 5

# Plot each line
plt.plot(numSites[start:], avgs90[start:], "r^", ms=10, lw=2, label="90%")
plt.errorbar(numSites[start:], avgs90[start:], yerr=errors90[start:], fmt="none", color='grey', elinewidth=2,capthick=3,capsize=7)
plt.plot(numSites[start:], avgs95[start:], "g^", ms=10, lw=2, label="95%")
plt.errorbar(numSites[start:], avgs95[start:], yerr=errors95[start:], fmt="none", color='grey', elinewidth=2,capthick=3,capsize=7)
plt.plot(numSites[start:], avgs97[start:], "b^", ms=10, lw=2, label="97%")
plt.errorbar(numSites[start:], avgs97[start:], yerr=errors97[start:], fmt="none", color='grey', elinewidth=2,capthick=3,capsize=7)

# Quadratic regression 2
end = len(numSites)
def f2(x, a, b):
    return a*(x**b)
coef, o = curve_fit(f2, numSites[fittingStart:end], avgs90[fittingStart:end])
rsquared = 1 - (np.sum((np.array(avgs90[fittingStart:end]) - f2(np.array(numSites[fittingStart:end]), *coef))**2) / np.sum((avgs90[fittingStart:end] - np.mean(avgs90[fittingStart:end]))**2))
quadratic = ("y = " + '{0:.2E}'.format(coef[0]) + " " + r"$x^{" + '{0:.5f}'.format(round(coef[1],rounding)) + "}$" + "    $R^2=" + str(round(rsquared, rounding)) + "$")
plt.plot(numSites[fittingStart:end], f2(np.array(numSites[fittingStart:end]), coef[0], coef[1]), '-r', label=quadratic)

# Quadratic regression 2
end = len(numSites)
def f2(x, a, b):
    return a*(x**b)
coef, o = curve_fit(f2, numSites[fittingStart:end], avgs95[fittingStart:end])
rsquared = 1 - (np.sum((np.array(avgs95[fittingStart:end]) - f2(np.array(numSites[fittingStart:end]), *coef))**2) / np.sum((avgs95[fittingStart:end] - np.mean(avgs95[fittingStart:end]))**2))
quadratic = ("y = " + '{0:.2E}'.format(coef[0]) + " " + r"$x^{" + '{0:.5f}'.format(round(coef[1],rounding)) + "}$" + "    $R^2=" + str(round(rsquared, rounding)) + "$")
plt.plot(numSites[fittingStart:end], f2(np.array(numSites[fittingStart:end]), coef[0], coef[1]), '-g', label=quadratic)

# Quadratic regression 2
end = len(numSites)
def f2(x, a, b):
    return a*(x**b)
coef, o = curve_fit(f2, numSites[fittingStart:end], avgs97[fittingStart:end])
rsquared = 1 - (np.sum((np.array(avgs97[fittingStart:end]) - f2(np.array(numSites[fittingStart:end]), *coef))**2) / np.sum((avgs97[fittingStart:end] - np.mean(avgs97[fittingStart:end]))**2))
quadratic = ("y = " + '{0:.2E}'.format(coef[0]) + " " + r"$x^{" + '{0:.5f}'.format(round(coef[1],rounding)) + "}$" + "    $R^2=" + str(round(rsquared, rounding)) + "$")
plt.plot(numSites[fittingStart:end], f2(np.array(numSites[fittingStart:end]), coef[0], coef[1]), '-b', label=quadratic)

# Exponential regression
end = len(numSites)
def f3(x, a, b):
    return b*np.exp(a*x)
coef, o = curve_fit(f3, numSites[fittingStart:end], avgs97[fittingStart:end])
rsquared = 1 - (np.sum((np.array(avgs97[fittingStart:end]) - f3(np.array(numSites[fittingStart:end]), *coef))**2) / np.sum((avgs97[fittingStart:end] - np.mean(avgs97[fittingStart:end]))**2))
expo = ("y = " + '{0:.2E}'.format(coef[1]) + " " + r"$e^{"+ '{0:.2f}'.format(coef[0]) + r"x}$" + "    $R^2=" + str(round(rsquared, rounding)) + "$")
plt.plot(numSites[fittingStart:end], f3(np.array(numSites[fittingStart:end]), coef[0], coef[1]), '-k', label=expo)

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
ax1.set_ylabel("time to optimise [s]", labelpad=20)
plt.xlim(min(numSites[start:])-1, max(numSites[start:])+1)
ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
plt.legend(prop={'size': 15})

# Nicer spacing
plt.tight_layout()

# Save as a png
plt.savefig('scalingSites.pdf', transparent=False)

# Repeat for the log plot
ax1.tick_params(axis='y', pad=10, labelsize=20)
ax1.tick_params(axis='x', pad=10, labelsize=20)
plt.xlabel("", fontsize=20, color='black')
plt.ylabel("", fontsize=20, color='black')
plt.grid(True)
ax1.set_xlabel("spins in linear chain", labelpad=10)
ax1.set_ylabel("time to optimise [s]", labelpad=20)
plt.xlim(min(numSites[start:])-1, max(numSites[start:])+1)
ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
plt.legend(prop={'size': 15})
ax1.set_xscale("log")
ax1.set_yscale("log")
plt.tight_layout()
plt.gcf().subplots_adjust(left=0.15)
plt.savefig('scalingSitesLog.pdf', transparent=False)


