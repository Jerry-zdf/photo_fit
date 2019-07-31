#!/usr/bin/python3

import mpmath
import os.path
import sys
import numpy as np

def o_function(l, k, r):
    res = mpmath.hyp1f1(-1j / k + l, 2 * l + 2, -2 * 1j * k * r)
    res *= (1./(2*l + 1)) * (k * r)**l
    res *= mpmath.gamma(1 + 1j / k)
    res *= mpmath.rf(-1j / k, l)
    res *= (-2j)**l / mpmath.fac(2 * l)
    res *= mpmath.exp(mpmath.pi/(2 * k))
    return res

l_max = 6
r_min = 0.01
r_step = 0.1
r_max = 100.

inputpath = "input2"

kvals = open(sys.argv[1], 'r').read().split()

for kval in kvals:
    k = float(kval)
    for l in range(l_max + 1):
        file_name = "input2/z1_k" + ("%.3f" % k) + "_l" + str(l) + ".dat"
        if os.path.isfile(file_name):
            continue

        output_file = open(file_name, 'w')
        r = r_min

        while r < r_max:
            z = o_function(l, k ,r)
            output_file.write("%10.6f %25.15e %25.15e\n" % (r, z.real, z.imag))    
            r += r_step

        output_file.close()
    
    os.system("build/Release/./photo_fit config.txt parameters.txt " + kval)
    

