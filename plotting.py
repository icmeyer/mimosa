import re
import numpy as np
from matplotlib import collections as mc
import matplotlib.pyplot as plt

from tools import get_trailing_numbers

def plotlines(lines='', circles = '', length = 1.26):
    fig, ax = plt.subplots()

    if lines:
        lc = mc.LineCollection(lines[0], colors = lines[1],  linewidths=1)
        # lc = mc.LineCollection(lines[0], lines[1])
        ax.add_collection(lc)
    if circles:
        for circle in circles:
            pltcircle = plt.Circle(circle[0], circle[1], color=circle[2],
                                   fill = False)
            ax.add_artist(pltcircle)
    # ax.autoscale()
    ax.set_xlim([-0.1,length+0.1])
    ax.set_ylim([-0.1,length+0.1])
    ax.margins(0.1)
    ax.set_aspect(1.0)
    plt.show()

def plot_from_rays(rays, regions, MATERIALS, length=1.26):
    fuelcolors = ['tab:orange', 'tab:green', 'tab:red', 'tab:purple',
                  'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']
    linesegs = []
    linecols = []
    for ray in rays:
        for segment in ray.segments:
            linesegs.append([segment.r0,segment.r1])
            mat = regions[segment.region].mat
            mod_counter = 0
            if mat == 'mod':
                mod_counter += 1
                linecols.append([0, 0, 1.0-mod_counter/9, 1])
            elif mat == 'fuel':
                linecols.append(fuelcolors[segment.region%9-1])
    lines = [linesegs, linecols]
    plotlines(lines=lines,length=length)



if __name__ == '__main__':
    segments = [[(0, 1), (1, 1)], [(2, 3), (-3, 2)], [(0, 2), (2, 3)]]
    colors = np.array([(1, 0, 0, 1), (0, 1, 0, 1), (0, 0, 1, 1)])
    lines = [segments, colors]
    circles = [[(0,0),1,'b'],[(1,1),2,'r']]
    plotlines(lines, circles = circles)
    # plotlines(lines, c)
    
