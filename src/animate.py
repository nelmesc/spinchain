#!/usr/bin/python3
import sys
from PIL import Image, ImageDraw
import numpy as np
import matplotlib.pyplot as plt
from io import StringIO

# Sizes in pixels
networkWidth = 800
dynamicsWidth = 800
height = 560
padding = 20

# TODO
startTime = 1
endTime = 10
steps = 100
skip = 5
gifSeconds = 5

# Load data from file
genomes = []
fidelities = []
with open(sys.argv[1], "r") as f:
    for line in f:
        split = line.split()
        genomes.append(split[0])
        fidelities.append([])
        for val in split[1:]:
            fidelities[-1].append(float(val))

networkImages = []
dynamicsImages = []
combinedImages = []

# For each instance
for index in range(0, len(genomes), skip):

    # DEBUG TODO
    print(index)

    # Create image for dynamics
    im1 = Image.new('RGB', (networkWidth, height), (255, 255, 255))
    draw1 = ImageDraw.Draw(im1)

    # Go through the genome TODO

        # Get the start letter

        # Get the end letter

        # Get the coupling val, inverting if needed

        # Get the vis info

    # Add the first node TODO

    # For each coupling TODO

        # Try to draw it if one node already placed

    # Generate dynamics graph as in dynamics.py but with slight formatting tweaks
    x = np.linspace(startTime, endTime, steps)
    y = np.array(fidelities[index])
    fig = plt.figure(figsize=(10,7), dpi=80)
    ax1 = plt.subplot()
    plt.title("generation " + str(index), fontsize=20)
    ax1.set_xlim([startTime, endTime])
    ax1.set_ylim([0.0, 1.0])
    plt.plot(x,y,color='black',lw=2)
    ax1.tick_params(axis='y', labelsize=20)
    ax1.tick_params(axis='x', labelsize=20)
    plt.xlabel('$\mathrm{time \cdot J_{max}}$',fontsize=25,color='black')
    plt.ylabel('${\cal{F}}(t)$',fontsize=25,color='black')
    plt.tight_layout()

    # Convert this matplotlib image to a PIL image
    fig.canvas.draw()
    im2 = Image.frombytes('RGB', fig.canvas.get_width_height(), fig.canvas.tostring_rgb())
    plt.close()

    # Create image for the combination
    im3 = Image.new('RGB', (dynamicsWidth+networkWidth+padding, height), (255, 255, 255))

    # Combine images
    im3.paste(im1, (0, 0))
    im3.paste(im2, (networkWidth+padding, 0))

    # Add all to arrays
    networkImages.append(im1)
    dynamicsImages.append(im2)
    combinedImages.append(im3)

    if index == 0:
        im1First = im1
        im2First = im2
        im3First = im3

# Save gif
perFrame = (gifSeconds * 1000 * skip) / len(genomes)
im1First.save('network.gif', save_all=True, append_images=networkImages, optimize=True, duration=perFrame, loop=0)
im2First.save('dynamics.gif', save_all=True, append_images=dynamicsImages, optimize=True, duration=perFrame, loop=0)
im3First.save('combined.gif', save_all=True, append_images=combinedImages, optimize=True, duration=perFrame, loop=0)
