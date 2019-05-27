#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import sys
from mpl_toolkits.axisartist.axislines import SubplotZero
from matplotlib.ticker import MultipleLocator
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['mathtext.fontset'] = 'cm'


kklist= sys.argv[1:]

l= [0, 1, 3, 7]
linestyles = ['-', '--', '-.', ':']
markerstyles = ['o', 'v', 's', 'D']
colors = ['red', 'blue', 'green', 'orange']

font_size = 12
line_width=1.0
marker_size=1.8

grid_stride = 7


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

for kk in kklist:
    k = float(kk)

    # REAL PATR
    plt.clf()
    fig, ax = plt.subplots(1, figsize=(5, 3.5))
    #plt.tight_layout()

    ax.spines['left'].set_position('zero')
    ax.spines['right'].set_color('none')
    ax.spines['bottom'].set_position('zero')
    ax.spines['top'].set_color('none')

    for i in range(len(l)):
        fg, fr, fi = get_fit(k, l[i])
        ig, ir, ii = get_input(k, l[i])

        ax.plot(fg, fr, color=colors[i] , linestyle=linestyles[i], linewidth=line_width, label="$l = $" + str(l[i]))
        ax.plot(ig[::grid_stride], ir[::grid_stride], markerstyles[i], markersize=marker_size, color="dark" + colors[i]  )

    #ax.set_xlabel("$r$", fontsize=font_size)
    #ax.set_ylabel("$k =$ "+"%.3f" % k + "     " + "$l=$ "+"%i" % l +, fontsize=font_size)

    bottom, top = plt.ylim()
    ypos = - (top - bottom) / 15.

    ax.text(32.0, ypos, "$r$", fontsize=font_size)

    ax.tick_params(axis='x', which="both",direction="inout")
    ax.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)

    ax.legend(fontsize=font_size*0.8)

    plt.title("$k  = $"+"%.3f" % k + "\n(real part)", fontsize=font_size)

    plt.savefig("fit_k"+ "%.3f" % k +"_real.png", dpi=300)



    # IMAG PATR
    plt.clf()
    fig, ax = plt.subplots(1, figsize=(5, 3.5))
    #plt.tight_layout()

    ax.spines['left'].set_position('zero')
    ax.spines['right'].set_color('none')
    ax.spines['bottom'].set_position('zero')
    ax.spines['top'].set_color('none')

    for i in range(len(l)):
        fg, fr, fi = get_fit(k, l[i])
        ig, ir, ii = get_input(k, l[i])

        ax.plot(fg, fi, color=colors[i] , linestyle=linestyles[i], linewidth=line_width, label="$l = $" + str(l[i]))
        ax.plot(ig[::grid_stride], ii[::grid_stride], markerstyles[i], markersize=marker_size, color="dark" + colors[i]  )

    #ax.set_xlabel("$r$", fontsize=font_size)
    #ax.set_ylabel("$k =$ "+"%.3f" % k + "     " + "$l=$ "+"%i" % l +, fontsize=font_size)

    bottom, top = plt.ylim()
    ypos = - (top - bottom) / 15.

    ax.text(32.0, ypos, "$r$", fontsize=font_size)

    ax.tick_params(axis='x', which="both",direction="inout")
    ax.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)

    ax.legend(fontsize=font_size*0.8)

    plt.title("$k  = $"+"%.3f" % k + "\n(imaginary part)", fontsize=font_size)

    plt.savefig("fit_k"+ "%.3f" % k +"_imag.png", dpi=300)