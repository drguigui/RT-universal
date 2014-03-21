#!/usr/bin/env python
# -*- coding:Utf-8 -*-
from __future__ import division
from pylab import *
import scipy.linalg as LS
from gg_math import *


def ComputeAGridOld(mu, deltgrid):
	siz = len(deltgrid)
	A = zeros((siz, siz))
	A[0, 0] =  - (mu / deltgrid[0] + 1)
	A[0, 1] =  (mu / deltgrid[1])
	for i in range(1, siz - 1):
		A[i, i - 1] = mu ** 2 / deltgrid[i - 1]
		A[i, i] = - (2 *mu ** 2 / deltgrid[i] + 1)
		A[i, i + 1] =  mu ** 2 / deltgrid[i + 1]
	A[-1, -2] =  - (mu / deltgrid[-2])
	#A[-1, -1] =  (mu / deltgrid[-1] - 1) # I guess the -1 should not be here. To check
	A[-1, -1] =  (mu / deltgrid[-1]) # I guess the -1 should not be here. To check
	return A


def ComputeAGrid(mu, deltgrid):
	siz = len(deltgrid)
	A = zeros((siz, siz))
	A[0, 0] =  - (mu / deltgrid[0] + 1)
	A[0, 1] =  (mu / deltgrid[0])
	for i in range(1, siz - 1):
		alpha = 2 * mu ** 2 / (deltgrid[i - 1] ** 2 + deltgrid[i] * deltgrid[i - 1])
		A[i, i - 1] = alpha
		A[i, i] = -alpha * (1 + deltgrid[i - 1] / deltgrid[i]) - 1
		A[i, i + 1] = alpha * (deltgrid[i - 1] / deltgrid[i]) 
	A[-1, -2] =  -(mu / deltgrid[-1])
	A[-1, -1] =  (mu / deltgrid[-1]) # I guess the -1 should not be here. To check
	return A

def SolveS(nbang, grid, deltgrid, I0, S, anglerad, gmu, gwt):
	""" Solve the equation 6.17 in Mihalas
	for only one frequency
	"""
	urez = zeros((nbang, len(grid)))
	for ang in range(nbang):
		mu = gmu[ang]
		# b = gwt[ang]
		Sr = -S.copy()
		Sr[0] = -I0[ang]
		Sr[-1] = 0
		A = ComputeAGrid(mu, deltgrid)
		u = LS.solve(A, Sr)
		urez[ang, :] = u
	return urez

def ComputeEddington(urez, gmu, gwt):
	""" Computes the Eddington factors for each altitude"""
	nbang, nbalt = urez.shape
	print urez.shape

	eddu = zeros((nbalt))
	eddo = zeros((nbalt))
	eddh = zeros((nbalt))
	for i in range(nbang):
		eddu += gwt[i] * gmu[i] ** 2 * urez[i, :]
		eddh += gwt[i] * gmu[i] * urez[i, :]
		eddo += gwt[i] * urez[i, :]

	edding = eddu / eddo
	return edding, (eddh / eddo)[0]

def ComputeVfromU(u, deltgrid, gmu):
	v = zeros((u).shape)
	for mu in range(len(gmu)):
		for i in range(len(v[0,:]) - 1):
			v[mu, i] = gmu[mu] / deltgrid[i] * (u[mu, i + 1] - u[mu, i])
		v[mu, -1] = gmu[mu] / deltgrid[-1] * (u[mu, -1] - u[mu, -2])
	return v



if __name__=="__main__":
	nbang = 8
	grid, deltgrid, deltgrid2 = TauGrid()
	anglerad, gmu, gwt = gaussangles(nbang)

	I0 = zeros((nbang))
	S = zeros((len(grid)))
	I0[-1] = 42
	S += 1 #IlovePython

	print "We solve the equation system"
	#urez = SolveS(nbang, grid, deltgrid, I0, S, anglerad, gmu, gwt)
	urez = SolveS(nbang, grid, deltgrid2, I0, S, anglerad, gmu, gwt)
	print urez

	print "Tau faible:", urez[0,:]
	eddfact, h0 = ComputeEddington(urez, gmu, gwt)
	print "Eddington factors"
	print eddfact

