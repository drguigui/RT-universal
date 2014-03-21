#!/usr/bin/env python
# -*- coding:Utf-8 -*-
from __future__ import division
from pylab import *


def ComputeWidth(vE):
	""" Computes the width of the grid vE"""
	resu = zeros((len(vE)))
	resu[0] = (vE[1] - vE[0]) / 2.
	for i in range(1, len(vE) - 1):
		resu[i] = (vE[i + 1] - vE[i - 1]) / 2.
	resu[-1] = (vE[-1] - vE[-2]) / 2.
	return resu


def ComputeWidth2(vE):
	""" Computes the width of the grid vE"""
	resu = zeros((len(vE)))
	for i in range(0, len(vE) - 1):
		resu[i] = (vE[i + 1] - vE[i]) 
	resu[-1] = resu[-2] 
	return resu



def TauGrid(tmin=1E-3, tmax=10, decnb = 5):
	""" Computes a grid in tau """
	leng = log(tmax / tmin) / log(10) * decnb
	lgrid = (linspace(log(tmin), log(tmax), leng ))
	grid = exp(lgrid.copy())
	return grid, ComputeWidth(grid), ComputeWidth2(grid)



def gaussangles(m):
	"""Calcul des angles gaussiens"""
	tol=1E-15
	if m<5:
		tol=1E-30 # truc pour la precision
	en=float(m)
	np1=m+1
	nnp1=m*np1
	cona=float(m-1)/float(8*m**3)
	lim= m // 2 + 1 ##### ATTENTION
	gmu=zeros(m,float)
	gwt=zeros(m,float)
	if m<1:
		gmu[0]=0.5
		gwt[0]=1
		return
	for k in range(1, lim):
		t=(4*k-1)*pi/float(4*m+2)
		x=cos(t+cona/tan(t))
		xi =x
		ggtop=-1
		tmp=0.
		pm2=1.
		while (abs(xi-x)>tol or ggtop==-1):
			x=xi
			ggtop=ggtop+1
			pm2=1.
			pm1=x
			p=0.
			for nn in range(2,m+1):
				p=((2*nn-1)*x*pm1-(nn-1)*pm2)/nn
				pm2=pm1
				pm1=p
			tmp=1./(1.-x**2)
			ppr=en*(pm2-x*p)*tmp
			p2pri=(2.*x*ppr-nnp1*p)*tmp
			xi=x-(p/ppr)*(1.+(p/ppr)*p2pri/(2.*ppr))
		gmu[k-1]=-x
		gwt[k-1]=2./(tmp*(en*pm2)**2)
		gmu[m-k]=-gmu[k-1]
		gwt[m-k]=gwt[k-1]
	if m%2!=0:
		gmu[lim]=0
		prod=1.
		for  k in range(3,m,2):
			prod*=float(k)/float(k-1)
		gwt[lim]=2./prod**2
	truc=gmu.copy()
	anglerad=gmu.copy()
	for i in range(0,len(gmu)):
		gmu[i]=0.5*gmu[i]+0.5
		gwt[i]=0.5*gwt[i]
		truc[i]=arccos(gmu[i])*180/pi
		anglerad[i]=arccos(gmu[i])
	#print gmu, gwt, truc, 1/gwt, anglerad
	print "Angles for m = ", m, " : ", truc
	print "Soit, en radians:", anglerad
	print "Soit, Mu = ", gmu
	#return dict(zip(anglerad,gwt))
	return anglerad, gmu, gwt

if __name__=="__main__":
	gaussangles(8)
	gaussangles(16)
	gaussangles(32)

