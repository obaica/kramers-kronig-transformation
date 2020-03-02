#!/usr/bin/env python
# -*- coding: utf-8 -*

'''
This script is used to calculate the real part of the dielectric function from the imaginary part.
Usage: python3 kk.py IMAG.in output
Part of the core codes are copied from the script@utf (https://github.com/utf/kramers-kronig). 
Author: Gang Tang
Date: 03/03/2020
Email: gangtang@outlook.com
'''

import numpy as np
import sys
import math
import scipy.integrate

def symmetrise(a):
    """Turn a XX, YY, ZZ, XY, YZ, XZ array into a symmetrical 3x3 matrix"""
    return [[a[0], a[3], a[5]], [a[3], a[1], a[4]], [a[5], a[4], a[2]]]
            
IMAG = sys.argv[1]

file_imag = open(IMAG,'r')
line = file_imag.readlines()
file_imag.close()

eps_imag = []
energies = []

for i in range(len(line)):
    eps_imag.append([float(line[i].strip().split()[1]),float(line[i].strip().split()[2]),float(line[i].strip().split()[3]),float(line[i].strip().split()[4]),float(line[i].strip().split()[5]),float(line[i].strip().split()[6])])
    energies.append(float(line[i].strip().split()[0]))

eps_imag = np.array([symmetrise(a) for a in eps_imag])

def kkr(de, eps_imag, cshift=1e-6):
    eps_imag = np.array(eps_imag)
    nedos = eps_imag.shape[0]
    cshift = complex(0, cshift)
    w_i = np.arange(0, nedos*de, de, dtype=np.complex_)
    w_i = np.reshape(w_i, (nedos, 1, 1))
    
    def integration_element(w_r):
        factor = w_i / (w_i**2 - w_r**2 + cshift)
        total = np.sum(eps_imag * factor, axis=0)
        return total * (2/math.pi) * de + np.diag([1, 1, 1])

    return np.real([integration_element(w_r) for w_r in w_i[:,0,0]])

eps_real_calc = kkr(energies[1] - energies[0], eps_imag)

b=[]
b=eps_real_calc

output = sys.argv[2]
with open(output,'w') as myfile:
    for i in range(len(b)):
        myfile.write('%15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f\n' %(energies[i],b[i][0][0],b[i][1][1],b[i][2][2],b[i][0][1],b[i][1][2],b[i][0][2]))
    myfile.close()

