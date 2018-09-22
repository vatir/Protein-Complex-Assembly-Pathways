# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 10:01:06 2015

@author: Koan
"""

import sys
sys.path.append('C:\Users\Koan\Dropbox\Proteasome Project\VSProjects\Cluster')

import numpy as np
from scipy.integrate import odeint
from collections import OrderedDict
import itertools

from ODERun import ODEFunc

Func = ODEFunc()
Func.GetSDataFromFiles('CRNData\\species3')
Func.GetRDataFromFiles('CRNData\\react3')

def unique(seq): # Probably Fast might use later
    seen = set()
    seen_add = seen.add
    return [ x for x in seq if not (x in seen or seen_add(x))]

Reactions = Func.Reactions
Species = np.arange(len(Func.Species))
Complexes = unique(unique([(x,) for x in np.array(Reactions)[:,3]])+unique(map(tuple,map(np.sort,np.array(Reactions)[:,:2]))))

S = np.mat(np.zeros([len(Species),len(Complexes)]))
E = np.mat(np.zeros([len(Complexes),len(Reactions)]))

for i, Complex in enumerate(Complexes):
	if len(Complex) == 2:
		S[Complex[0], i] = 1.0
		S[Complex[1], i] = 1.0
	else:
		S[Complex[0], i] = 1.0

for i, Reaction in enumerate(Reactions):
	for j, Complex in enumerate(Complexes):
		if tuple(map(int,map(np.sort,Reaction[:2]))) == Complex:
			E[j,i] = 1
		if (Reaction[3],) == Complex:
			E[j,i] = -1

L = S*E