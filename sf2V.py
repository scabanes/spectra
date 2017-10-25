#! /usr/bin/env python
#
def Fsf2v ( lat, sflat, R ):

# Ici on fait une boucle sur un profil en latitude pour trouver sa derivee premiere, soit la 
# vitesse zonal u_sf issue de la streamfunction sflat; u_sf = 1/(r) * d(sf)/dtheta. Est 
# est le rayon spherique. la derivee partielle est obtenue a partir d'une difference finie centree.
# une version plus aboutie de ces calculs est dans sfvp.py
	import numpy as np
	import math
	import scipy.special as sp
	import matplotlib.pyplot as plt

	u_sf = np.zeros(len(lat))
	# zonal component u_sf = 1/(r) * d(sf)/dtheta
	for iLat in range(1, len(lat)-1):
		u_sf[iLat] = (1/R)*(sflat[iLat+1]-sflat[iLat-1])/(2.*(np.pi/180.)*(lat[iLat]-lat[iLat-1]))
	
	#plt.plot(lat,u_sf)
	#plt.show()	

	return (u_sf)
