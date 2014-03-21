#!/usr/bin/env python
# -*- coding:Utf-8 -*-
from __future__ import division
from pylab import *
import scipy.linalg as LS
from gg_math import *
from solveS import TauGrid


def ComputeJGrid(edd, deltgrid, alpha):
	siz = len(deltgrid)
	A = zeros((siz, siz))
	A[0, 0] =  - (edd[0] / deltgrid[0] + 1)
	A[0, 1] =  (edd[1] / deltgrid[1])
	for i in range(1, siz - 1):
		A[i, i - 1] = edd[i - 1] / deltgrid[i - 1]
		A[i, i] = - (2 * edd[i] / deltgrid[i] + 1 - alpha)
		A[i, i + 1] =  edd[i + 1] / deltgrid[i + 1]
	A[-1, -2] =  - (edd[-2] / deltgrid[-2])
	A[-1, -1] =  (edd[-1] / deltgrid[-1]) # = 1 /3 dB/dTmax no -1 here!
	return A

def ComputeJGrid(edd, deltgrid, alpha2, h0):
	siz = len(deltgrid)
	A = zeros((siz, siz))
	A[0, 0] =  - (edd[0] / deltgrid[0] + h0)
	A[0, 1] =  (edd[1] / deltgrid[0])
	for i in range(1, siz - 1):
		alpha = 2 / (deltgrid[i - 1] ** 2 + deltgrid[i] * deltgrid[i - 1])
		A[i, i - 1] = alpha * edd[i - 1]
		A[i, i] = -alpha * edd[i] * (1 + deltgrid[i - 1] / deltgrid[i]) + alpha2 - 1
		A[i, i + 1] = alpha * edd[i + 1]* (deltgrid[i - 1] / deltgrid[i]) 
	A[-1, -2] =  - (edd[-2] / deltgrid[-1])
	A[-1, -1] =  (edd[-1] / deltgrid[-1]) # = 1 /3 dB/dTmax no -1 here!
	return A


def SolveU(grid, deltgrid, edding, alpha, beta, h0):
	""" Solves the equation 6-42 in Mihalas"""
	b = ones((len(grid))) * -beta
	Jeq = ComputeJGrid(edding, deltgrid, alpha, h0)
	b[-1] = 0
	b[0] = 0 #Always true... we have h0
	J = LS.solve(Jeq, b)
	return J

def ComputeS(grid, deltgrid, edding, alpha, beta, h0):
	""" Solve equation 6-42 in Mihalas and returns S"""
	J = SolveU(grid, deltgrid, edding, alpha, beta,h0)
	return alpha * J + beta

if __name__=="__main__":
	print "Hello world"
	grid, deltgrid, deltgrid2 = TauGrid()
	edding = array([-0.07159427, -0.04723267, -0.03318441, -0.02133561, -0.00999498, 0.00131193, 0.01289635, 0.02500979, 0.03792283, 0.0522074, 0.06912271, 0.09090694, 0.12061064, 0.16061824, 0.20938418, 0.25907735, 0.29873668, 0.32005077, 0.30677671, 0.51854665])

	epsilon = 0.1
	h0 = 0.1
	alpha = 1 - epsilon
	beta = epsilon
	print "Solution"
	print SolveU(grid, deltgrid2, edding, alpha, beta, h0)
