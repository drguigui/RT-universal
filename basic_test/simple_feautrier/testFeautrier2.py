#!/usr/bin/env python
# -*- coding:Utf-8 -*-
from __future__ import division
from pylab import *
import scipy.linalg as LS
from gg_math import *
from solveU import ComputeS
from solveS import SolveS, ComputeEddington, TauGrid, ComputeVfromU




def SolveFeautrier():
	nbang = 8
	grid, tmp, deltgrid = TauGrid(decnb = 20)
	anglerad, gmu, gwt = gaussangles(nbang)
	I0 = zeros((nbang))
	S = zeros((len(grid)))
	I0[-1] = 1
	#S += 1 #IlovePython

	epsilon = 0.00001
	alpha = 1 - epsilon
	beta = epsilon


	for i in range(5):
		urez = SolveS(nbang, grid, deltgrid, I0, S, anglerad, gmu, gwt)
		edding, h0 = ComputeEddington(urez, gmu, gwt)
		S = ComputeS(grid, deltgrid, edding, alpha, beta, h0)

		print i, urez
		print "================="
	
	#print urez

	print "We now compute v"
	v =  ComputeVfromU(urez, deltgrid, gmu)
	
	print "Now I+"
	print urez + v
	print "Now I -"
	print urez - v



if __name__=="__main__":
	SolveFeautrier()


