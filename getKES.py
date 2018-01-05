#! /usr/bin/env python
from netCDF4 import Dataset as NetCDFFile 
import numpy as np
import math
import scipy.special as sp
from sympy.functions.elementary.miscellaneous import sqrt
from pylab import *
###########################################################################################################################
#													FUNCTION 1
###########################################################################################################################
def KES ( array_BLM, N, R_sat, getField ):
# Pour faire des spectres matrices completes. BLM_m est equivalent a
# Blm_m seulement ce dernier est tronque en L et M. Le calcul de 
# l energie sont correspond a E(n)=1/2*Blm^2 pour m=0 (coeffs reel) et 
# E(n)=1/2*(Blm*conj(Blm) + Blm*conj(Blm)) pour m/=0 avec une somme sur
# les m=-l a m 
# si j ai une streamfunction psi projetee sur HS alors j applique un facteur
# CL = l(l+1)/(2*r^2) pcq (P. Augier et al. 2013):
# E 	= 1/2 (V.V) 
#	= 1/2  r^2/l(l+1) Real(rot(u)*rot(u) + div(u)*div(u))
# avec rot(u) = Laplacian psi and div(u) = Laplacian phi
# avec psi et phi la fonction courant et le portentel de vitesse et 
# * le complex conjuge. On a finalement:
# E	= l(l+1)/2r^2  (psi_lm* psi_lm + phi_lm* phi_lm)
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#===================================================================
#		    Kinetic Energy Spectra (KES)
#===================================================================
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	n = (np.arange(0.,N,1.))
	nn = np.transpose(range(0,N)+np.zeros([N-1, N]))
	Cl = n*(n+1.)/(2.*R_sat**2.) #see: P. Augier & E. Lindborg 2013 expression (12-13)
	CCl = nn*(nn+1.)/(2.*R_sat**2.)#see: P. Augier & E. Lindborg 2013 expression (12-13)

# ==> details of the energy variables
#Eno[l] ac l=0:L pour m=0  ac L=N
#En[l] ac l=0:L sum_m=1^M  ac M=N
#Emn[l,m] ac m=1:M ac M=N
# Energie-----------------
	if('sf' in getField) and (not 'vp' in getField):	# --> 'sf'
		print 'get sf Gauss Coeffs ####################'
		BLM = array_BLM['BLM0']
		#
		Eno = Cl*abs(BLM[:,0])**2						
		En = Cl*np.real(np.sum(BLM[:,1:]*conj(BLM[:,1:]) + BLM[:,1:]*conj(BLM[:,1:]),1))
		Emn = CCl*np.real(BLM[:,1:]*conj(BLM[:,1:]) + BLM[:,1:]*conj(BLM[:,1:]))

	elif ('vp' in getField) and (not 'sf' in getField):	# --> 'vp'
		print 'get vp Gauss Coeffs ####################'
		BLM = array_BLM['BLM0']
		#
		Eno = Cl*abs(BLM[:,0])**2						
		En = Cl*np.real(np.sum(BLM[:,1:]*conj(BLM[:,1:]) + BLM[:,1:]*conj(BLM[:,1:]),1))
		Emn = CCl*np.real(BLM[:,1:]*conj(BLM[:,1:]) + BLM[:,1:]*conj(BLM[:,1:]))

	elif ('sf' and 'vp' in getField):			# --> 'sf' & 'vp'
		print 'get sf & vp Gauss Coeffs ####################'
		if(getField[0]=='sf'):
			BLM_sf = array_BLM['BLM0']
			BLM_vp = array_BLM['BLM1']
		else:
			BLM_sf = array_BLM['BLM1']
			BLM_vp = array_BLM['BLM0']
		#
		Eno = Cl*(abs(BLM_sf[:,0])**2 + abs(BLM_vp[:,0])**2)			
		En = Cl*np.real(np.sum(BLM_sf[:,1:]*conj(BLM_sf[:,1:]) + BLM_sf[:,1:]*conj(BLM_sf[:,1:])  +  BLM_vp[:,1:]*conj(BLM_vp[:,1:]) + BLM_vp[:,1:]*conj(BLM_vp[:,1:]),1))
		Emn = CCl*np.real(BLM_sf[:,1:]*conj(BLM_sf[:,1:]) + BLM_sf[:,1:]*conj(BLM_sf[:,1:])  +  BLM_vp[:,1:]*conj(BLM_vp[:,1:]) + BLM_vp[:,1:]*conj(BLM_vp[:,1:]))

	else:							# --> 'T' or 'M' or 'Z'
		print 'get Gauss Coeffs ####################'
		BLM = array_BLM['BLM0']
		#
		Eno = (1./2.)*abs(BLM[:,0])**2			
		En = (1./2.)*np.real(np.sum(BLM[:,1:]*conj(BLM[:,1:]) + BLM[:,1:]*conj(BLM[:,1:]),1))
		Emn = (1./2.)*np.real(BLM[:,1:]*conj(BLM[:,1:]) + BLM[:,1:]*conj(BLM[:,1:]))#


	return ( Eno, En, Emn )
