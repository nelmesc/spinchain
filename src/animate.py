#!/usr/bin/python3
import sys
from PIL import Image, ImageDraw, ImageFont
import numpy as np
import matplotlib.pyplot as plt
from io import StringIO
import math

# Sizes in pixels
networkWidth = 800
dynamicsWidth = 800
height = 560
padding = 20
nodeSize = 40
betweenNodes = 100
maxSize = 5000
couplingWidth = 5

# TODO
startTime = 1
endTime = 10
steps = 100
skip = 2
gifSeconds = 20

# Heatmap function, val between 0 and 1 -> colour tuple
def heatmap(val):

    # blue, cyan, green, yellow, red
    cols = [(0,0,255),(0,255,255),(0,255,0),(255,255,0),(255,0,0)]
    vals = [0.000,0.250,0.500,0.750,1.000]

    indexLow = int(val * 4.0)
    indexHigh = indexLow + 1

    if indexHigh >= len(vals): return cols[-1]

    ratio = (val-vals[indexLow])/(vals[indexHigh]-vals[indexLow])
    col = [int(cols[indexLow][i] + (cols[indexHigh][i]-cols[indexLow][i])*ratio) for i in range(3)]

    return tuple(col)

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

# Init image arrays
networkImages = []
dynamicsImages = []
combinedImages = []

# For each instance
for index in range(0, len(genomes), skip):

    # Create image for dynamics
    im1 = Image.new('RGB', (maxSize, maxSize), (255, 255, 255))
    draw1 = ImageDraw.Draw(im1)

    # Create empty coupling array
    couplings = [] 
    lettersAdded = {}
    genome = genomes[index]
    maxCoupling = -99999

    # Determine how many digits are used per coupling
    digitsPer = 0
    i = 2
    while genome[i] in "0123456789":
        digitsPer += 1
        i += 1
    visStart = genome.find("#")

    # Go through the genome 
    for i in range(0, visStart, 2+digitsPer):

        newCoupling = {}

        # If in ascii order, coupling is positive
        if ord(genome[i]) < ord(genome[i+1]):

            # Get the letters
            newCoupling["a"] = genome[i]
            newCoupling["b"] = genome[i+1]

            # Get the coupling val
            newCoupling["val"] = int(genome[i+2:i+2+digitsPer])

        # If in reverse ascii order, coupling is negative
        else:

            # Get the letters
            newCoupling["a"] = genome[i+1]
            newCoupling["b"] = genome[i]

            # Get the coupling val
            newCoupling["val"] = -int(genome[i+2:i+3+digitsPer])

        # Get the vis info
        newCoupling["dir"] = int(genome[visStart+1+int(i/(2+digitsPer))], 16) * 22.5 * (math.pi/180.0)

        if abs(newCoupling["val"]) > maxCoupling: maxCoupling = abs(newCoupling["val"])

        # Add to the couplings array
        couplings.append(newCoupling)

    # Add the first node 
    x = maxSize/2
    y = maxSize/2
    draw1.ellipse([x-nodeSize/2,y-nodeSize/2,x+nodeSize/2,y+nodeSize/2],fill=(0,0,0))
    lettersAdded[couplings[0]["a"]] = [x, y]

    # Keep track of min/max x/y values
    minX = x
    maxX = x
    minY = y
    maxY = y

    # Repeat to make sure all nodes placed
    for iteration in range(10):

        # For each coupling
        for coupling in couplings:

            # If the first node in the coupling has already been placed
            if coupling["a"] in list(lettersAdded.keys()) and coupling["b"] not in list(lettersAdded.keys()):
                prevX = lettersAdded[coupling["a"]][0]
                prevY = lettersAdded[coupling["a"]][1]
                x = prevX + math.cos(coupling["dir"])*betweenNodes
                y = prevY + math.sin(coupling["dir"])*betweenNodes
                draw1.line([prevX,prevY,x,y], fill=heatmap(coupling["val"]/maxCoupling), width=couplingWidth)
                draw1.ellipse([x-nodeSize/2,y-nodeSize/2,x+nodeSize/2,y+nodeSize/2],fill=(0,0,0))
                draw1.ellipse([prevX-nodeSize/2,prevY-nodeSize/2,prevX+nodeSize/2,prevY+nodeSize/2],fill=(0,0,0))
                lettersAdded[coupling["b"]] = [x, y]
                if x < minX: minX = x
                if x > maxX: maxX = x
                if y < minY: minY = y
                if y > maxY: maxY = y

            # If the second node in the coupling has already been placed
            elif coupling["b"] in list(lettersAdded.keys()) and coupling["a"] not in list(lettersAdded.keys()):
                prevX = lettersAdded[coupling["b"]][0]
                prevY = lettersAdded[coupling["b"]][1]
                x = prevX + math.cos(coupling["dir"])*betweenNodes
                y = prevY + math.sin(coupling["dir"])*betweenNodes
                draw1.line([prevX,prevY,x,y], fill=heatmap(coupling["val"]/maxCoupling), width=couplingWidth)
                draw1.ellipse([x-nodeSize/2,y-nodeSize/2,x+nodeSize/2,y+nodeSize/2],fill=(0,0,0))
                draw1.ellipse([prevX-nodeSize/2,prevY-nodeSize/2,prevX+nodeSize/2,prevY+nodeSize/2],fill=(0,0,0))
                lettersAdded[coupling["a"]] = [x, y]
                if x < minX: minX = x
                if x > maxX: maxX = x
                if y < minY: minY = y
                if y > maxY: maxY = y

            # If both nodes are there, draw the coupling anyway
            elif coupling["a"] in list(lettersAdded.keys()) and coupling["b"] in list(lettersAdded.keys()):
                prevX = lettersAdded[coupling["a"]][0]
                prevY = lettersAdded[coupling["a"]][1]
                x = lettersAdded[coupling["b"]][0]
                y = lettersAdded[coupling["b"]][1]
                draw1.line([prevX,prevY,x,y], fill=heatmap(coupling["val"]/maxCoupling), width=couplingWidth)
                draw1.ellipse([x-nodeSize/2,y-nodeSize/2,x+nodeSize/2,y+nodeSize/2],fill=(0,0,0))
                draw1.ellipse([prevX-nodeSize/2,prevY-nodeSize/2,prevX+nodeSize/2,prevY+nodeSize/2],fill=(0,0,0))

    # Crop so everything is centered, leaving space for the title and the colour bar
    im1 = im1.crop((minX-nodeSize, minY-nodeSize-40, maxX+nodeSize, maxY+nodeSize+40))

    # Write the generation title above it
    draw1 = ImageDraw.Draw(im1)
    fnt = ImageFont.truetype("Pillow/Tests/fonts/FreeSans.ttf", 20)
    draw1.text((im1.size[0]/2-70,15), "generation " + str(index), font=fnt, fill=(0,0,0))

    # Draw a heatmap bar
    for x in range(0, im1.size[0]):
        col = heatmap(x/float(im1.size[0]))
        draw1.line([x, im1.size[1], x, im1.size[1]-30], fill=col, width=1)

    # Scale if too large
    if im1.size[0] > networkWidth:
        im1 = im1.resize((networkWidth, int(im1.size[1]*(networkWidth/im1.size[0]))))
    if im1.size[1] > height:
        im1 = im1.resize((int(networkWidth*(height/img1.size[1])), height))

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
    im3.paste(im1, (int(networkWidth/2-(maxX-minX)/2-nodeSize), int(height/2-(maxY-minY)/2-nodeSize)))
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
