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
####################################################################
perturb = False
scales = True
####################################################################
# 		    LOAD filePTS ....
####################################################################
E=np.load('EE.npz')
const=np.load('const.npz')
_En_=E['_En_']
_Eno_=E['_Eno_']
_Emn_=E['_Emn_']
_EZ_=E['_EZ_']
_ET_=E['_ET_']
_ER_=E['_ER_']
R_sat=const['R_sat']
omega_sat=const['omega_sat']
Ck=const['Ck']
Cz=const['Cz']
NBV = const['NBV']
H = const['H']

#################################
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
print _En_.shape

####################################################################
#                   COMPUTATIONS...
####################################################################

######################################################################
# To make quantities non-dimensional we use [omega] = 1/s as unity  
# for time and [R_sat] = m as length unity.
n = np.arange(0.,N,1.)
n_long = np.arange(0.,N+500,1.)
#######################
## calculate epsilon ##
#######################
nnn = 150 # to align the line
value = _En_[nnn,-AST4]/((omega_sat*R_sat)**2.)
fac = Ck*(n_long[nnn]**(-5./3.))
norm = ((R_sat**2.)*(omega_sat**3.))**(2./3.)
fac = fac / norm
epsilon = value/fac
epsilon = epsilon**(3./2.)
print "EPSILON",epsilon
#######################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Constant beta and coriolis parameter f varries in latitude. Here a
# typical value is choosen at mid-latitude pi/4.
angle = np.pi/4.
beta_sat = (2*omega_sat*np.cos(angle))/R_sat
f = 2*omega_sat*np.sin(angle)
URMS = 29.
# ---------------------------------------------------------------------
#               Numerical scaling results:
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
if _EZ_ is not None:
  print 'theory & numerics'
  print 'zmf =',gamma/(1.+gamma), '&', _EZ_/_ET_
  print 'theory & numerics'
  print 'nzmf =',1./(1.+gamma), '&', _ER_/_ET_
  print '-----------------------------------------------------------------'
print '####################################################################'


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

#
#plt.loglog(n,_En_[:,0]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
plt.loglog(n,_En_[:,1]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
plt.loglog(n,_En_[:,3]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
plt.loglog(n,_En_[:,5]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
plt.loglog(n,_En_[:,8]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
plt.loglog(n,_En_[:,10]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
plt.loglog(n,_En_[:,12]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
plt.loglog(n,_En_[:,14]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
plt.loglog(n,_En_[:,19]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
plt.loglog(n,_En_[:,20]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
plt.loglog(n,_En_[:,21]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
plt.loglog(n,_En_[:,22]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
plt.loglog(n,_En_[:,25]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
plt.loglog(n,_En_[:,30]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
plt.loglog(n,_En_[:,40]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
plt.loglog(n,_En_[:,50]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
plt.loglog(n,_En_[:,60]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
plt.loglog(n,_En_[:,70]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
plt.loglog(n,_En_[:,80]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
plt.loglog(n,_En_[:,90]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
plt.loglog(n,_En_[:,100]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
plt.loglog(n,_En_[:,110]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
plt.loglog(n,_En_[:,130]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
plt.loglog(n,_En_[:,140]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
plt.loglog(n,_En_[:,150]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
plt.loglog(n,_En_[:,200]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
plt.loglog(n,_En_[:,210]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
plt.loglog(n,_En_[:,220]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
plt.loglog(n,_En_[:,240]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
plt.loglog(n,_En_[:,260]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
plt.loglog(n,_En_[:,280]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
plt.loglog(n,_En_[:,300]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
plt.loglog(n,_En_[:,320]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
plt.loglog(n,_En_[:,325]/((omega_sat*R_sat)**2.), '-', color='red', linewidth=1) # PSD des m/=0
plt.loglog(n,_En_[:,:]/((omega_sat*R_sat)**2.), '-', color='red', linewidth=1) # PSD des m/=0

plt.loglog(n,_En_[:,-AST4]/((omega_sat*R_sat)**2.), '-', color='blue', linewidth=1) # PSD des m/=0



#plt.loglog(n,_En_[:,1]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1.5) # PSD des m/=0
#plt.loglog(n,_En_[:,3]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
#plt.loglog(n,_En_[:,7]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
#plt.loglog(n,_En_[:,14]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
#plt.loglog(n,_En_[:,30]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
#plt.loglog(n,_En_[:,60]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
#plt.loglog(n,_En_[:,80]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1) # PSD des m/=0
#
plt.xlabel('n')
plt.ylabel('Ek')
plt.ylim((10**(-16),10**(-5)))

plt.grid()
plt.show(figA)
#################################################################################
#			   Final spectra
#################################################################################
figB, axB = plt.subplots()
trans = mtransforms.blended_transform_factory(axB.transData, axB.transAxes)
plt.loglog(n[:-1],_Eno_[:-1,-AST4]/((omega_sat*R_sat)**2.), '-', color='darkred', linewidth=2, label="zonal (m=0)") # PSD des m/=0
plt.loglog(n[:-1],_En_[:-1,-AST4]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=2, label="residual (m>0)") # PSD des m/=0
# -------------------Theoretical spectra
plt.loglog(n_long[1:],Cz*((beta_sat/beta_sat)**2.)*(n_long[1:]**(-5.)), '--', color='darkred', linewidth=2.5, label=r'$n^{-5}$')
plt.loglog(n_long[1:],Ck*((epsilon/((R_sat**2.)*(omega_sat**3.)))**(2./3.))*(n_long[1:]**(-5./3.)), '--', color='black', linewidth=2.5, label=r'$n^{-5/3}$')
nnn = 50 ; zefit = (n_long[nnn]**(-3.)) / (_En_[nnn,-AST4]/((omega_sat*R_sat)**2.))
plt.loglog(n_long[1:],n_long[1:]**(-3.)/zefit, '--', color='gray', linewidth=2.5, label=r'$n^{-3}$')
if perturb:
  plt.loglog(n_long[1:],0.5*Ck*((epsilon/((R_sat**2.)*(omega_sat**3.)))**(2./3.))*(n_long[1:]**(-8./3.)), '--', color='darkgreen', linewidth=2.5)
# -------------------Typpical scales
if scales:
  plt.axvline(n_beta, color='grey', linestyle='dotted', linewidth=2.5,  alpha=0.7)
  plt.axvline(n_R, color='grey', linestyle='dotted', linewidth=2.5,  alpha=0.7)
  plt.axvline(n_z, color='grey', linestyle='dotted', linewidth=2.5, alpha=0.7)
  ##plt.axvline(n_D, color='black', linestyle='solid', linewidth=1)
  #axB.fill_between(np.arange(n_D,n_D+100,1), 0, 1, facecolor='grey', alpha=0.4, transform=trans)
# -------------------numerical spectra
#plt.loglog(n,Eno[:,2]/((omega_sat*R_sat)**2.), '-', color='red', linewidth=1.5) # PSD de m=0, 
#plt.loglog(n,En[:,2]/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1.5) # PSD des m/=0
#

#####################"
if perturb:
  #plt.loglog(n,_n_Emn_[:,80]/((omega_sat*R_sat)**2.), '-', color='darkgreen', linewidth=2) # PSD des m/=0
  plt.loglog(n,_Emn_[:,0,-AST4]/((omega_sat*R_sat)**2.), '-', color='darkgreen', linewidth=1) # PSD des m/=0
  plt.loglog(n,_Emn_[:,1,-AST4]/((omega_sat*R_sat)**2.), '-', color='darkgreen', linewidth=1) # PSD des m/=0
  plt.loglog(n,_Emn_[:,2,-AST4]/((omega_sat*R_sat)**2.), '-', color='darkgreen', linewidth=1) # PSD des m/=0
  plt.loglog(n,_Emn_[:,10,-AST4]/((omega_sat*R_sat)**2.), '-', color='darkgreen', linewidth=1) # PSD des m/=0
#####################

plt.xlabel('n')
plt.ylabel('Ek')
plt.ylim((10**(-14),10**(-5)))
plt.xlim(1.,1000.)
plt.grid()
plt.legend(loc="best")
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
