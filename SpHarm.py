#! /usr/bin/env python
from netCDF4 import Dataset as NetCDFFile 
import numpy as np
import math
import scipy.special as sp
from math import factorial
from pylab import *
from decimal import *
###########################################################################################################################
#													FUNCTION 1
###########################################################################################################################
def getGaussCoeffs ( N, data_file):
# GaussCoeffs mais dans une matrice M_br et M_bi les coeffs de Gauss
# des fichiers xspec_'-field'_xx. Les coeffs reel et imaginaire des
# harmonique sont obtenus par la transformation: 
# Blm = C*(M_br-M_bi) + i* C*(M_br+M_bi), la partie m=0 est reel pure.
# seul un facteur 1/2 intervient. Cette transformation est faite a 
# partir des coeffs obtenus par SPHEREPACK.
# BLM est une matrice contenant tous les M=0:l et l=0:L.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#===================================================================
#			Load Data
#===================================================================
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	M_lm = np.zeros([N,N])
	M_br = np.zeros([N,N])
	M_bi = np.zeros([N,N])
	#=====>>  Load Here
	g = open(data_file, "r")
	# Je passe les en tete..
	header1 = g.readline()
	header2 = g.readline()
	header3 = g.readline()
	header4 = g.readline()
	field = header3[2:8]
	il=-1
	im=0
	# Boucle sur les lignes de g
	for line in g:
		columns = line.split()
		lm = (columns[0])
		br = (columns[1])
		bi = (columns[2])
		il=il+1
		if(lm != '#'):
			M_lm[il,im] = lm
			M_br[il,im] = br
			M_bi[il,im] = bi
			br=(float(br))
			bi=(float(bi))
		elif(lm == '#'):
			im=im+1
			il=-1
	g.close()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#===================================================================
# 			Arrange Coeffs
#===================================================================
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	M = np.arange(0,N,1.)# M[M]	
	CS_M = ((-1)**M[1:])# CS_M[M] avec M=1:L
	BLM_m = (1/(2*CS_M*math.sqrt(2)) * (M_br[:,1:] - M_bi[:,1:])) + 1j* (1/(2*CS_M*math.sqrt(2)) * (M_br[:,1:] + M_bi[:,1:])) #Blm[l,m]
	BLM_mo = 0.5 * M_br[:,0]
	BLM = np.c_[BLM_mo.reshape((-1, 1)),BLM_m]# BLM[L,M]

	return (BLM, field)

###########################################################################################################################
#													FUNCTION 2
###########################################################################################################################
def projectionSH ( BLM, L, M, pas ):
#===================================================================
#		Spherical Harmonics Projection
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Variables ........................................
# L and M maximum latitudinal and zonal degrees. Pas fixes the steps
# Phi and theta. x defines as cos(theta) and f is the real function
# once projected on spherical harmonics. Mm is a table giving the 
# value of zonal order m and eimphi is the azimuthal function sets
# as a table. Blm_m is a table of the transformation from initial real
# gauss coefficients obtained in spherepack with  Shsec computing coeff
# from a scalar function. Blm_mo, same than Blm_m for m=0. Plm are the
# Legender polynomials normalized in geodesy fashion by norm_geod. The
# spherical harmonic function is Ylm in its complexe form and transform
# into Ylm_p, Ylm_o, Ylm_m, in its real form for respectively the positive
# zero and negative zonal orders. f compute the final scalar field once 
# summed. 
# Here it starts.....................................
	#L=10 # nombre max de l
	#M=4 # nombre max de m
	#pas=0.005
	theta=(np.arange(0,(math.pi+pas),pas))
	phi=np.arange(0,(2*math.pi+pas),pas)
	x=np.cos(theta) # cos(-pi/2):cos(-pi/2)
	mm=len(x) # nombre de point que l on souhaite
	f_f = np.zeros([len(x),len(phi)]) # fonction total du champ de vitesse
	f_R = np.zeros([len(x),len(phi)]) # fonction des residus soit m/=0
# Zonal function.....................................
	Mm = range(0,M+1)+np.zeros([len(phi), M+1])#Mm[cst(len(phi)),m]
	eimphi=np.transpose(np.exp(1j * phi.reshape((-1, 1)) * Mm))#eimphi[m,phi]
# Blm gauss coefficients.............................
	Blm = BLM[0:L+1,0:M+1] # Ici on tronque BLM pour alleger la projection SH
	CS = ((-1)**Mm[0,1:]) # facteur (-1)^m servant plus loin.
# LEGENDRE POLYNOMIALS...............................
	Plm = np.zeros([len(x),L+1,M+1])
	norm_geod = np.zeros([L+1,M+1]) 
	for m in range(0, M+1):
		for l in range(m, L+1):
			norm_geod = (((2*np.float(l)+1)*factorial(np.float(l)-np.float(m)))/(2*factorial(np.float(l)+np.float(m))))**(0.5)
			Plm[:,l,m] = (norm_geod * sp.lpmv(m,np.float(l),x))#Plm2(mm,l,m)
# Spherical harmonics projection......................
#Iphi=100
	Ylm_p = np.zeros([len(x),len(phi)])
	Ylm_o = np.zeros([len(x)])
	Ylm_m = np.zeros([len(x),len(phi)])
	for itheta in range(0, len(theta)):
# complex final form 
# Here Blm are complex coefficients rearanged from spherepack coeff making the decomposed field on the form:
# f(theta,phi) = sum_l sum_m Blm * Plm(theta) * eimphi
		for iphi in range(0, len(phi)):
			Ylm = Blm[0:L+1,1:M+1] * Plm[itheta,:,1:] * eimphi[1:,iphi]

			Ylm_p[itheta,iphi] = np.real(sum((1/math.sqrt(2)) * (CS * conj(Ylm) + CS * Ylm)))
			Ylm_o[itheta] = np.real(sum(Blm[0:L+1,0] * Plm[itheta,:,0]))
			Ylm_m[itheta,iphi] = np.real(sum(((1j)/math.sqrt(2)) * (CS * conj(Ylm) - CS * Ylm)))

			f_f[itheta,iphi] =  Ylm_p[itheta,iphi] + Ylm_o[itheta] + Ylm_m[itheta,iphi]
			f_R[itheta,iphi] =  Ylm_p[itheta,iphi] + Ylm_m[itheta,iphi]

	return (f_f,f_R,Ylm_o,theta,phi,x)
