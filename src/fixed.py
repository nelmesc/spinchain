#! /usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from matplotlib.ticker import MaxNLocator
import re
import sys

# Set up the figure
plt.figure(figsize=(10,8))

# Create the axes
ax1=plt.subplot()

# Load the data from the sites scaling file
numSites = []
ex1 = []
ex2 = []
ex3 = []
with open("scalingFixed.dat", "r") as f:
    for line in f:
        split = line.split()
        numSites.append(int(split[0]))
        ex1.append(float(split[1]))
        ex2.append(float(split[2]))
        ex3.append(float(split[3]))

# Plot each line
plt.plot(numSites, ex1, "g^", ms=10, lw=2, label="1 excitation")
plt.plot(numSites, ex2, "b^", ms=10, lw=2, label="2 excitations")
plt.plot(numSites, ex3, "r^", ms=10, lw=2, label="3 excitations")

# To give the fits more detail
bonusSites = np.arange(3,8,0.1)

# Linear fit
rounding = 3
def f1(x, a):
    return a*x
coef, o = curve_fit(f1, numSites, ex1)
rsquared = 1 - (np.sum((np.array(ex1) - f1(np.array(numSites), *coef))**2) / np.sum((ex1 - np.mean(ex1))**2))
linear = ("y = " + '{0:.2f}'.format(round(coef[0],2)) + r"$x$" + "    $R^2=" + str(round(rsquared, rounding)) + "$")
plt.plot(bonusSites, f1(np.array(bonusSites), *coef), '--g', label=linear)

# Quadratic fit
def f2(x, a, b):
    return a*(x**2) + b*x
coef, o = curve_fit(f2, numSites, ex2)
rsquared = 1 - (np.sum((np.array(ex2) - f2(np.array(numSites), *coef))**2) / np.sum((ex2 - np.mean(ex2))**2))
linear2 = ("y = " + '{0:.2f}'.format(round(coef[0],2)) + r"$x^2$" + " " + '{0:+.2f}'.format(round(coef[1],2)) + r"$x$" + "    $R^2=" + str(round(rsquared, rounding)) + "$")
plt.plot(bonusSites, f2(np.array(bonusSites), *coef), '--b', label=linear2)

# Polynomial fit
def f3(x, a, b, c):
    return a*(x**3) + b*(x**2) + c*x
coef, o = curve_fit(f3, numSites, ex3)
rsquared = 1 - (np.sum((np.array(ex3) - f3(np.array(numSites), *coef))**2) / np.sum((ex3 - np.mean(ex3))**2))
linear3 = ("y = " + '{0:.2f}'.format(round(coef[0],2)) + r"$x^3$" + '{0:+.2f}'.format(round(coef[1],2)) + r"$x^2$" + " " + '{0:+.2f}'.format(round(coef[2],2)) + r"$x$"  + "    $R^2=" + str(round(rsquared, rounding)) + "$")
plt.plot(bonusSites, f3(np.array(bonusSites), *coef), '--r', label=linear3)

# Set the size of the axis ticks
ax1.tick_params(axis='y', pad=10, labelsize=20)
ax1.tick_params(axis='x', pad=10, labelsize=20)

# Set the names of the axis labels
plt.xlabel("", fontsize=20, color='black')
plt.ylabel("", fontsize=20, color='black')
plt.grid(True)
ax1.set_xlabel("spins in linear chain", labelpad=10)
ax1.set_ylabel("time for 200 generations [s]", labelpad=20)
plt.xlim(min(numSites)-1, max(numSites)+1)
ax1.xaxis.set_major_locator(MaxNLocator(integer=True))

plt.legend(prop={'size': 20})

# Nicer spacing
plt.tight_layout()

# Save as a png
plt.savefig('scalingFixed.pdf', transparent=False)


