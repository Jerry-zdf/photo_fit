#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import sys
from mpl_toolkits.axisartist.axislines import SubplotZero
from matplotlib.ticker import MultipleLocator
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['mathtext.fontset'] = 'cm'


k= float(sys.argv[1])
l= int(sys.argv[2])

font_size = 12
line_width=1.0
marker_size=1.5


path_input = "//home/mateusz/workspace/photo_fit/input/"
path_fit = "/home/mateusz/workspace/photo_fit/output/"

def get_input(kk , ll):
    inpp = path_input + "z1_k" + "%.3f" % kk + "_l" + "%i" % ll + ".dat"
    input = np.loadtxt(inpp).transpose()
    return input[0], input[1] * (2*ll+1), input[2] * (2*ll+1)

def get_fit(kk , ll):
    fitp = path_fit + "fit_z1_k" + "%.3f" % kk + ".dat"
    fit = np.loadtxt(fitp, skiprows=(2 + ll*11), usecols=(1, 2, 3), max_rows=10)

    grid = np.linspace(start=0, stop=30, num=1000)

    fit_re = np.zeros_like(grid)
    fit_im = np.zeros_like(grid)

    for trip in fit:
        fit_re += grid**ll * trip[1] * np.exp(-trip[0] * grid**2) * trip[0]**(0.75 + 0.5*ll) * (2*ll+1)
        fit_im += grid**ll * trip[2] * np.exp(-trip[0] * grid**2) * trip[0]**(0.75 + 0.5*ll) * (2*ll+1)

    return grid, fit_re, fit_im

fg, fr, fi = get_fit(k, l)
ig, ir, ii = get_input(k, l)



fig, ax = plt.subplots(1, figsize=(5, 3.5))
plt.tight_layout()

ax.spines['left'].set_position('zero')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_position('zero')
ax.spines['top'].set_color('none')

ax.plot(fg, fr, color="red" , linestyle=(0, ())          , linewidth=line_width)
ax.plot(fg, fi, color="blue", linestyle=(0, (3, 1, 1, 1)), linewidth=line_width)

ax.plot(ig, ir, 'o', markersize=marker_size, color="darkred"  )
ax.plot(ig, ii, 's', markersize=marker_size, color="darkblue" )

ax.set_xlabel("$r$", fontsize=font_size)
#ax.set_ylabel("$k =$ "+"%.3f" % k + "     " + "$l=$ "+"%i" % l +, fontsize=font_size)

bottom, top = plt.ylim()
ypos = top - (top - bottom) / 7.

ax.text(24.0, ypos, "$k =$ "+"%.3f" % k + "\n$l=$ "+"%i" % l, fontsize=font_size)

ax.tick_params(axis='x', which="both",direction="inout")
ax.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)

plt.savefig("fit_k"+ "%.3f" % k + "_l" + "%i" % l +".png", dpi=300)