#! /usr/bin/env python
from netCDF4 import Dataset as NetCDFFile 
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import numpy as np
import math
import scipy.special as sp
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches
from pylab import *
from decimal import *
# Simon's codes
import SpHarm
import getKES
import getSumEnergy
# Aymeric's codes
import ppplot
####################################################################
# 		    LOAD filePTS ....
####################################################################
#-For SpHarm.py, which field to get:
# getField = ['sf','vp,'T','M','Z']
# ----------------------
f = open('filePTS.zono', "r")
print f.readline()
print f.readline()
print f.readline()
f.readline()
# ----------------------
getField= (f.readline().split()[1:])
FieldNetCDEF = f.readline().split()[1]
Niteration = f.readline().split()[1]
N=int(f.readline().split()[1])
L=int(f.readline().split()[1])
M=int(f.readline().split()[1])
it=int(f.readline().split()[1])
iz=int(f.readline().split()[1])
iplt=int(f.readline().split()[1])
pas=float(f.readline().split()[1])
# file of velocity/streamfunction I upload to plot. Il doit contenir
# une seule profondeur iz=0 and une seule iteration en temps it=0.  
# Le generer avec sfvp_t.py --> time_range=[1200,1201]
#nc = NetCDFFile('sfvpData-1200-1200.nc')
nc = NetCDFFile(FieldNetCDEF)
print 'Loaded.....'
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
# ----------> Set array that contain suitable varaibles.
data_file = ['']*len(getField) # empty matrix for strings
array_BLM = {}
array_field = {}
array_f_f = {}
array_f_R = {}
array_Ylm_o = {}
# ----------> Extract & Array data on field selected
for iField in range(0, len(getField)):
	# Name the data file
	data_file[iField] = 'xspec_'+getField[iField]+'_'+Niteration
	# Load Gauss Coefficients
	(BLM, field) = SpHarm.getGaussCoeffs(N, data_file[iField])
	# Spherical Harmonics Projection
	(f_f,f_R,Ylm_o,theta,phi,x) = SpHarm.projectionSH( BLM, L, M, pas )
	# Varriables in an array
	array_BLM['BLM{0}'.format(iField)] = BLM
	array_field['field{0}'.format(iField)] = field	
	array_f_f['f_f{0}'.format(iField)] = f_f
	array_f_R['f_R{0}'.format(iField)] = f_R
	array_Ylm_o['Ylm_o{0}'.format(iField)] = Ylm_o
	del BLM, field, f_f, f_R, Ylm_o
# ----------> Get Energy from array for the field
# Kinetic Energy Spectra (KES)
(Eno, En, Emn) = getKES.KES ( array_BLM, N, R_sat, getField )
# Energetic Integrals
(E_mo, E_R, E_T, URMS, latSH, field)=getSumEnergy.SumEnergy(nc, array_Ylm_o, array_field, array_f_f, array_f_R, getField , Eno, En, pas, x, theta, phi, R_sat, it, iz, iplt)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#===================================================================
# 			PLOT SPECTRA
#===================================================================
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# To make quantities non-dimensional we use [omega] = 1/s as unity  
# for time and [R_sat] = m as length unity.
epsilon_dim = (R_sat**2.)*(omega_sat**3.)
#
n = np.arange(0.,N,1.)
n_long = np.arange(0.,N+500,1.)
# Constant beta and coriolis parameter f varries in latitude. Here a
# typical value is choosen at mid-latitude pi/4.
angle = np.pi/4.
beta_sat = (2*omega_sat*np.cos(angle))/R_sat
f = 2*omega_sat*np.sin(angle)
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
print 'zmf =',gamma/(1.+gamma), '&', np.sum(Eno)/(np.sum(Eno) + np.sum(En))
print 'theory & numerics'
print 'nzmf =',1./(1.+gamma), '&', np.sum(En)/(np.sum(Eno) + np.sum(En))
print '-----------------------------------------------------------------'
print '####################################################################'
# -------------------numerical spectra
#fig = plt.figure()
font = {'family' : 'Times New Roman',
        'size'   : 12}
fig, ax = plt.subplots()
trans = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
plt.loglog(n,Eno/((omega_sat*R_sat)**2.), '-', color='darkred', linewidth=1.5) # PSD de m=0, E = Eno/(omega.R)^2 = non-dim.
plt.loglog(n,En/((omega_sat*R_sat)**2.), '-', color='black', linewidth=1.5) # PSD des residus, E = En/(omega.R)^2 = non-dim
plt.loglog(n,Emn[:,0:3:1]/((omega_sat*R_sat)**2.), '-', color='darkgreen', linewidth=2, alpha=0.2)
# -------------------Theoretical spectra
plt.loglog(n_long[1:],Ck*((epsilon/((R_sat**2.)*(omega_sat**3.)))**(2./3.))*(n_long[1:]**(-5./3.)), '--', color='black', linewidth=2.0)
plt.loglog(n_long[1:],Cz*((beta_sat/beta_sat)**2.)*(n_long[1:]**(-5.)), '--', color='darkred', linewidth=2.0)
plt.loglog(n_long[1:],0.5*Ck*((epsilon/((R_sat**2.)*(omega_sat**3.)))**(2./3.))*(n_long[1:]**(-8./3.)), '--', color='darkgreen', linewidth=2.0)
# -------------------Typpical scales
plt.axvline(n_beta, color='grey', linestyle='dashed', linewidth=2)
plt.axvline(n_R, color='grey', linestyle='dashed', linewidth=2)
plt.axvline(n_z, color='grey', linestyle='dashed', linewidth=2)
ax.fill_between(np.arange(n_D,n_D+10,1), 0, 1, facecolor='grey', alpha=0.4, transform=trans)
plt.xlabel('n',**font)
plt.ylabel('E(n)',**font)
plt.title(field+' velocity spectra',**font)
#plt.ylim((10**(-16),10**(-5)))
#plt.xlim((1,N+5))
plt.grid()
plt.show(fig)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
