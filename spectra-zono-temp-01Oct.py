#! /usr/bin/env python
from netCDF4 import Dataset as NetCDFFile 
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import numpy as np
import math
#from mpl_toolkits.basemap import Basemap
import scipy.special as sp
import matplotlib.patches as mpatches
import matplotlib.font_manager as font_manager
from pylab import *
from decimal import *
# Simon's codes
import SpHarm
import getKES
####################################################################
# 		    LOAD filePTS ....
####################################################################
# file of spectral coefficient I upload for SpHarm.py
# getField = ['sf','vp,'T','M','Z']
# getField = ['sf','vp']
# ----------------------
f = open('filePTS.zono.temp', "r")
print f.readline()
print f.readline()
print f.readline()
f.readline()
# ----------------------
getField= (f.readline().split()[1:])
N=int(f.readline().split()[1])
# ----------------------
# Arrays to step in time for averaging. Each array has a given AST:
# "Average Step Time" which give the sleeping for averaging. One has
# to concatenate all the arrays in stept_file.
iAT11,iAT12 = f.readline().split()[1:3]
iAT21,iAT22 = f.readline().split()[1:3]
iAT31,iAT32 = f.readline().split()[1:3]
iAT41,iAT42 = f.readline().split()[1:3]
#
arrayT1 = np.arange(int(iAT11),int(iAT12),1) 
arrayT2 = np.arange(int(iAT21),int(iAT22),1) 
arrayT3 = np.arange(int(iAT31),int(iAT32),1) 
arrayT4 = np.arange(int(iAT41),int(iAT42),1) 
#
AST1 = int(f.readline().split()[1])
AST2 = int(f.readline().split()[1])
AST3 = int(f.readline().split()[1])
AST4 = int(f.readline().split()[1])
#
if(AST1!=0) and (AST2==0):
	stept_file = arrayT1
	AST4=AST1
elif(AST1!=0) and (AST2!=0) and (AST3==0):
	stept_file = np.concatenate((arrayT1, arrayT2),axis=0)
	AST4=AST2
elif(AST1!=0) and (AST2!=0) and (AST3!=0) and (AST4==0):
	stept_file = np.concatenate((arrayT1, arrayT2, arrayT3),axis=0)
	AST4=AST3
elif(AST1!=0) and (AST2!=0) and (AST3!=0) and (AST4!=0):
	stept_file = np.concatenate((arrayT1, arrayT2, arrayT3, arrayT4),axis=0)
else:
	print 'erreur dans filePTS.zono.temp'
print 'Loaded.....'
# --------------------------------------------------
print f.readline()
print f.readline()
print f.readline()
# Omega rotation rate
omega_sat = float(f.readline().split()[1])# 1.64*10.**(-6.) rotation rate of Saturn in s-1 from Gastin-2014-PEPI
# Planetary radius
R_sat = float(f.readline().split()[1])
# Spectra's constant
Ck = float(f.readline().split()[1])
epsilon = float(f.readline().split()[1]) # 1.5*10.**(-7.) m^2/s^3
Cz = float(f.readline().split()[1])
# Brunt Vaisala frequency [NBV] and atmosphere deepness [H]
H = float(f.readline().split()[1])
NBV = float(f.readline().split()[1])
print 'Loaded.....'
####################################################################
####################################################################



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#===================================================================
#			LOAD DATA
#===================================================================
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NT = len(stept_file) # nombre iteration en temps
Eno = np.zeros([N,NT])
En = np.zeros([N,NT])
Emn = np.zeros([N,N-1,NT]) # n inclue pas le m=0
#
EZ = np.zeros([NT])
ER = np.zeros([NT])
ET = np.zeros([NT])
#
for it in range(0, NT):
	print 'it file = ',stept_file[it]
	print 'itg = ',it
	# ----------> Set array that contain suitable varaibles.
	data_file = ['']*len(getField) # empty matrix for strings
	array_BLM = {}
	array_field = {}
	# ----------> Extract & Array data on field selected
	for iField in range(0, len(getField)):
		# Name the data file
		data_file[iField] = 'xspec_'+getField[iField]+'_'+str(stept_file[it])
		# Load Gauss Coefficients
		(BLM, field) = SpHarm.getGaussCoeffs(N, data_file[iField])
		# Varriables in an array
		array_BLM['BLM{0}'.format(iField)] = BLM
		array_field['field{0}'.format(iField)] = field	
		del BLM, field
	# ----------> Get Energy from array for the field
	# Kinetic Energy Spectra (KES)
	(Eno[:,it], En[:,it], Emn[:,:,it]) = getKES.KES ( array_BLM, N, R_sat, getField )

# ====================== ENERGY BUDGET ============================
# ===> Cas m=0
	print 'IcI------------------------- Field is ',getField,' -----------'
	print 'Spectral Coeffs'
	print 'E_mo =',np.sum(Eno[:,it])
# ===> Cas m/=0
	print 'E_R =',np.sum(En[:,it])
	print 'E_T =',np.sum(Eno[:,it]) + np.sum(En[:,it])
	print 'End'
	
	EZ[it] = np.sum(Eno[:,it])
	ER[it] = np.sum(En[:,it])
	ET[it] = np.sum(Eno[:,it]) + np.sum(En[:,it])
# ==================================================================

print '-------------ITERATIONS TERMINEES////////////////////////////'

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#===================================================================
#			SOME AVERAGE
#===================================================================
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
_En_ = np.zeros([N,NT])
_Eno_ = np.zeros([N,NT])
_Emn_ = np.zeros([N,N-1,NT])
#
try: arrayT2
except NameError: arrayT2 = [0]
try: arrayT3
except NameError: arrayT3 = [0]
try: arrayT4
except NameError: arrayT4 = [0]
#
try: AST2
except NameError: AST2 = 0
try: AST3
except NameError: AST3 = 0
try: AST4
except NameError: AST4 = 0

for itg in range(0, NT):
		# => first to arrayT1
	if(itg <= len(arrayT1)-AST1):
		_En_[:,itg] = np.mean(En[:,itg:itg+AST1],1)
		_Eno_[:,itg] = np.mean(Eno[:,itg:itg+AST1],1)
		_Emn_[:,:,itg] = np.mean(Emn[:,:,itg:itg+AST1],2)
		# => next to arrayT2
	if(itg > len(arrayT1)-1 and itg <= len(arrayT1)+len(arrayT2)-AST2):
		_En_[:,itg] = np.mean(En[:,itg:itg+AST2],1)
		_Eno_[:,itg] = np.mean(Eno[:,itg:itg+AST2],1)
		_Emn_[:,:,itg] = np.mean(Emn[:,:,itg:itg+AST2],2)
		# => next to arrayT3
	if(itg > len(arrayT1)+len(arrayT2)-1 and itg <= len(arrayT1)+len(arrayT2)+len(arrayT3)-AST3):
		_En_[:,itg] = np.mean(En[:,itg:itg+AST3],1)
		_Eno_[:,itg] = np.mean(Eno[:,itg:itg+AST3],1)
		_Emn_[:,:,itg] = np.mean(Emn[:,:,itg:itg+AST3],2)
		# => next to arrayT4
	if(itg > len(arrayT1)+len(arrayT2)+len(arrayT3)-1 and itg <= len(arrayT1)+len(arrayT2)+len(arrayT3)+len(arrayT4)-AST4):
		_En_[:,itg] = np.mean(En[:,itg:itg+AST4],1)
		_Eno_[:,itg] = np.mean(Eno[:,itg:itg+AST4],1)
		_Emn_[:,:,itg] = np.mean(Emn[:,:,itg:itg+AST4],2)
		#
	if(itg >= 3):
		_EZ_ = np.mean(EZ[itg-3:itg])
		_ER_ = np.mean(ER[itg-3:itg])
		_ET_ = np.mean(ET[itg-3:itg])

####################################################################
####################################################################

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#===================================================================
# 			PLOT SPECTRA
#===================================================================
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# To make quantities non-dimensional we use [omega] = 1/s as unity  
# for time and [R_sat] = m as length unity.
n = np.arange(0.,N,1.)
n_long = np.arange(0.,N+500,1.)
# Constant beta and coriolis parameter f varries in latitude. Here a
# typical value is choosen at mid-latitude pi/4.
angle = np.pi/4.
beta_sat = (2*omega_sat*np.cos(angle))/R_sat
f = 2*omega_sat*np.sin(angle)
URMS = 29.
# ---------------------------------------------------------------------
#	       	Numerical scaling results:
# Indexes n_R is the Rhines scale; n_beta is the transition scale in
# the zonostrophic spectra; n_z is the scale at which the flow beco-
# mes isotropic; n_D is the first Rossby deformation radius. R_beta
# is the zonostrophic index giving the extension on the zonostrophic
# inertial range.
# ---------------------------------------------------------------------
print '================================================================'
print '-------------------- SCALES ------------------------------------'
n_R = R_sat * (beta_sat/(2*URMS))**0.5 # rad/m
n_beta = ((Cz/Ck)**(3./10.))*((((beta_sat/beta_sat)**3.)/(epsilon/((R_sat**2.)*(omega_sat**3.))))**(1./5.))
n_z = n_beta*(2.*n_beta)**(3./7.)
L_D = (NBV*H)/(np.pi*f)
n_D = (np.pi*R_sat)/L_D
R_beta = n_beta/n_R
print 'URMS =',URMS
print 'n_beta =', n_beta 
print 'n_R = ', n_R
print 'n_z =', n_z
print 'R_beta = ', R_beta
print '----------------------'
print '----------------------'
gamma = (3./4.)*((10.*Cz)**(-5./6.))*R_beta**(10./3.)
print 'theory & numerics'
print 'zmf =',gamma/(1.+gamma), '&', _EZ_/_ET_
print 'theory & numerics'
print 'nzmf =',1./(1.+gamma), '&', _ER_/_ET_
print '-----------------------------------------------------------------'
print '####################################################################'

####################################################################
#			Instantaneous spectra
####################################################################
#font = {'family' : 'helvetica',
#        'size'   : 12}
# -------------------numerical spectra
fig = plt.figure()
plt.loglog(n,Eno/((omega_sat*R_sat)**2.), '-', color='red', linewidth=1.5) # PSD de m=0, E = Eno/(omega.R)^2 = non-dim.
plt.loglog(n,En/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1.5) # PSD des residus, E = En/(omega.R)^2 = non-dim.
#plt.loglog(n,Emn[:,0:40:10]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=2, alpha=0.2)
#plt.loglog(n,n_Emn/((omega_sat*R_sat)**2.), '-', color='green', linewidth=1)
# -------------------Theoretical spectra
plt.loglog(n_long[1:],Ck*((epsilon/((R_sat**2.)*(omega_sat**3.)))**(2./3.))*(n_long[1:]**(-5./3.)), '--', color='black', linewidth=2.5)
plt.loglog(n_long[1:],Cz*((beta_sat/beta_sat)**2.)*(n_long[1:]**(-5.)), '--', color='darkred', linewidth=2.5)
plt.loglog(n_long[1:],0.5*Ck*((epsilon/((R_sat**2.)*(omega_sat**3.)))**(2./3.))*(n_long[1:]**(-8./3.)), '--', color='green', linewidth=2.5)
# -------------------Typpical scales
plt.axvline(n_beta, color='grey', linestyle='dashed', linewidth=2.5)
plt.axvline(n_R, color='grey', linestyle='dashed', linewidth=2.5)
plt.axvline(n_z, color='grey', linestyle='dashed', linewidth=2.5)
plt.axvline(n_D, color='black', linestyle='solid', linewidth=1)
plt.xlabel('n')
plt.ylabel('zonal velocity energy distribution')
#plt.ylim((10**(-16),10**(-5)))
plt.xlim((1,N+5))
plt.grid()
plt.show(fig)
####################################################################
#			Averaged spectra
####################################################################
figA, axA = plt.subplots()
trans = mtransforms.blended_transform_factory(axA.transData, axA.transAxes)
# -------------------Theoretical spectra
plt.loglog(n_long[1:],Ck*((epsilon/((R_sat**2.)*(omega_sat**3.)))**(2./3.))*(n_long[1:]**(-5./3.)), '--', color='black', linewidth=2.5)
#plt.loglog(n_long[1:],Cz*((beta_sat/beta_sat)**2.)*(n_long[1:]**(-5.)), '--', color='darkred', linewidth=2.5)
#plt.loglog(n_long[1:],0.5*Ck*((epsilon/((R_sat**2.)*(omega_sat**3.)))**(2./3.))*(n_long[1:]**(-8./3.)), '--', color='green', linewidth=2.0)
# -------------------Typpical scales
#plt.axvline(n_beta, color='grey', linestyle='dashed', linewidth=2.5,  alpha=0.7)
plt.axvline(n_R, color='grey', linestyle='dashed', linewidth=2.5,  alpha=0.7)
#plt.axvline(n_z, color='grey', linestyle='dashed', linewidth=2.5, alpha=0.7)
#plt.axvline(n_D, color='black', linestyle='solid', linewidth=1.5)
axA.fill_between(np.arange(n_D,n_D+100,1), 0, 1, facecolor='grey', alpha=0.4, transform=trans)
# -------------------numerical spectra
#plt.loglog(n,Eno[:,2]/((omega_sat*R_sat)**2.), '-', color='red', linewidth=1.5) # PSD de m=0, 
plt.loglog(n,En[:,2]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
#
#plt.loglog(n,_En_[:,0]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
#plt.loglog(n,_En_[:,1]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
plt.loglog(n,_En_[:,3]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
plt.loglog(n,_En_[:,10]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
#plt.loglog(n,_En_[:,20]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
#plt.loglog(n,_En_[:,30]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
#plt.loglog(n,_En_[:,40]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
#plt.loglog(n,_En_[:,70]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
#plt.loglog(n,_En_[:,100]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
#plt.loglog(n,_En_[:,200]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0

#plt.loglog(n,_En_[:,1]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1.5) # PSD des m/=0
#plt.loglog(n,_En_[:,3]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
#plt.loglog(n,_En_[:,7]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
#plt.loglog(n,_En_[:,14]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
#plt.loglog(n,_En_[:,30]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
#plt.loglog(n,_En_[:,60]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
#plt.loglog(n,_En_[:,80]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
#
plt.loglog(n,_En_[:,-AST4]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=2.5) # PSD des m/=0
plt.xlabel('n')
plt.ylabel('Ek')
plt.ylim((10**(-16),10**(-5)))
plt.xlim((1,N+5))
plt.grid()
plt.show(figA)
#################################################################################
#			   Final spectra
#################################################################################
figB, axB = plt.subplots()
trans = mtransforms.blended_transform_factory(axB.transData, axB.transAxes)
# -------------------Theoretical spectra
plt.loglog(n_long[1:],Ck*((epsilon/((R_sat**2.)*(omega_sat**3.)))**(2./3.))*(n_long[1:]**(-5./3.)), '--', color='black', linewidth=2.5)
plt.loglog(n_long[1:],Cz*((beta_sat/beta_sat)**2.)*(n_long[1:]**(-5.)), '--', color='darkred', linewidth=2.5)
plt.loglog(n_long[1:],0.5*Ck*((epsilon/((R_sat**2.)*(omega_sat**3.)))**(2./3.))*(n_long[1:]**(-8./3.)), '--', color='darkgreen', linewidth=2.5)
# -------------------Typpical scales
plt.axvline(n_beta, color='grey', linestyle='dashed', linewidth=2.5,  alpha=0.7)
plt.axvline(n_R, color='grey', linestyle='dashed', linewidth=2.5,  alpha=0.7)
plt.axvline(n_z, color='grey', linestyle='dashed', linewidth=2.5, alpha=0.7)
#plt.axvline(n_D, color='black', linestyle='solid', linewidth=1)
axB.fill_between(np.arange(n_D,n_D+100,1), 0, 1, facecolor='grey', alpha=0.4, transform=trans)
# -------------------numerical spectra
#plt.loglog(n,Eno[:,2]/((omega_sat*R_sat)**2.), '-', color='red', linewidth=1.5) # PSD de m=0, 
#plt.loglog(n,En[:,2]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1.5) # PSD des m/=0
#
plt.loglog(n,_Eno_[:,-AST4]/((omega_sat*R_sat)**2.), '-', color='darkred', linewidth=2) # PSD des m/=0
#
plt.loglog(n,_En_[:,-AST4]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=2) # PSD des m/=0
#plt.loglog(n,_n_Emn_[:,80]/((omega_sat*R_sat)**2.), '-', color='darkgreen', linewidth=2) # PSD des m/=0
plt.loglog(n,_Emn_[:,0,-AST4]/((omega_sat*R_sat)**2.), '-', color='darkgreen', linewidth=1) # PSD des m/=0
plt.loglog(n,_Emn_[:,1,-AST4]/((omega_sat*R_sat)**2.), '-', color='darkgreen', linewidth=1) # PSD des m/=0
plt.loglog(n,_Emn_[:,2,-AST4]/((omega_sat*R_sat)**2.), '-', color='darkgreen', linewidth=1) # PSD des m/=0
plt.loglog(n,_Emn_[:,10,-AST4]/((omega_sat*R_sat)**2.), '-', color='darkgreen', linewidth=1) # PSD des m/=0
plt.xlabel('n')
plt.ylabel('Ek')
plt.ylim((10**(-16),10**(-5)))
plt.xlim((1,N+5))
plt.grid()
plt.show(figB)
####################################################################
#			Temporal energy evolution
####################################################################
figC = plt.figure()
plt.plot(stept_file,EZ, '.', color='darkred', linewidth=2.5)
plt.plot(stept_file,ER, '.', color='black', linewidth=2.5)
plt.plot(stept_file,ET, '.', color='navy', linewidth=2.5)
plt.yscale('log')
plt.ylabel('Ek')
plt.xlabel('Time')
plt.show(figC)
