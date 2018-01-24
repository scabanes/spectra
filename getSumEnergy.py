#! /usr/bin/env python
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import scipy.special as sp
from pylab import *
# Here is WindsPharm
from windspharm.standard import VectorWind #permet de travailler avec la streamfunction. Package a installer
# Simon's codes
import sf2V
###########################################################################################################################
#													FUNCTION 1
###########################################################################################################################
def SumEnergy (nc, array_Ylm_o, array_field, array_f_f, array_f_R, getField , Eno, En, pas, x, theta, phi, R_sat, it, iz, iplt):
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#===================================================================
#		    	  Get NetCDEF data
#===================================================================
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	lat = nc.variables['latitude'][:] #(128,)
	lon = nc.variables['longitude'][:]
	ulocity = nc.variables['u'][:,:,:,:] # [t?,Prof?,theta,phi] Zonal velocity
	vlocity = nc.variables['v'][:,:,:,:] # [t?,Prof?,theta,phi] meridional velocity
	sflocity = nc.variables['sf'][:,:,:,:] # [t?,Prof?,theta,phi] meridional velocity
	vplocity = nc.variables['vp'][:,:,:,:] # [t?,Prof?,theta,phi] meridional velocity
	# Varribles of the field analyzed (Afield) & the velocity to compare with (Avitesse)
	Avitesse = np.zeros([ulocity.shape[0],ulocity.shape[1],ulocity.shape[2],ulocity.shape[3]])
	Afield = np.zeros([ulocity.shape[0],ulocity.shape[1],ulocity.shape[2],ulocity.shape[3]])
	u = np.zeros([ulocity.shape[2],ulocity.shape[3]])
	v = np.zeros([ulocity.shape[2],ulocity.shape[3]])
	u_sf_mo = 0.
	# !!!! Attention Phi=0 est decale de pi.
	# Pour ploter les profil en latitude il faut determiner 
	# une valeur de phi a laquelle on veut se positionner.
	# de plus il y a un decallage entre phi=0 dans ma proj-
	# ection harmonique spherique et les donnees NetCDEF. 
	# que l on corrige pour les plots a suivre avec la val- 
	# eur de Ratio. Enfin, latSH est la latitude en degres
	# issu du pas SpHarm pour retomber sur le vecteur lat 
	# des donnees NetCDEF:
	Iphi=100
	Ratio = Iphi*float(len(lon))/float(len(phi))
	latSH = (180/math.pi)*(theta-math.pi/2)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#===================================================================
#		    Select field to analyze
#===================================================================
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	#------- if getfield
	if('sf' in getField) and (not 'vp' in getField):	# --> 'sf'
		Afield = sflocity
		Ylm_o = array_Ylm_o['Ylm_o0']
		f_f = array_f_f['f_f0']
		f_R = array_f_R['f_R0']
		field = array_field['field0']
		# Je cherche le profil de vitesse m=0 a partir
		# de sf est l equivalent de _u_sf_= np.mean(u_
		# sf[:,:],1):
		(u_sf_mo) = sf2V.Fsf2v(latSH, Ylm_o, R_sat)
		# Velocity to Compare with, ici la vitesse non-
		# divergente contenant toute l energie de sf:
		w = VectorWind(ulocity[it,iz,:,:], vlocity[it,iz,:,:],rsphere=R_sat)#windspharm
		u, v = w.nondivergentcomponent()#windspharm
		Avitesse[it,iz,:,:] = 0.

	elif ('vp' in getField) and (not 'sf' in getField):	# --> 'vp'
		Afield = vplocity
		Ylm_o = array_Ylm_o['Ylm_o0']
		f_f = array_f_f['f_f0']
		f_R = array_f_R['f_R0']
		field = array_field['field0']
		# Velocity to Compare with, ici l ecoulement
		# divergent contenant toute l energie du vp:
		w = VectorWind(ulocity[it,iz,:,:], vlocity[it,iz,:,:],rsphere=R_sat)#windspharm
		u, v = w.divergentcomponent()#windspharm
		Avitesse[it,iz,:,:] = 0.

	elif ('sf' and 'vp' in getField):			# --> 'sf' & 'vp'
		Afield = sflocity
		Ylm_o =   np.zeros([len(theta)])# 0.
		f_f =  np.zeros([len(theta),len(phi)])# 0.
		f_R =  np.zeros([len(theta),len(phi)])# 0.
		field = array_field['field0'] + array_field['field1']
		# Velocity to Compare with, ecouolement dive-
		# rgent et non-divergent soit vitesse totale:
		u[:,:] = 0.
		v[:,:] = 0.
		Avitesse[it,iz,:,:] =  np.sqrt((vlocity**2) + (ulocity**2))#vitesse to compare with

	elif ('T' in getField):					# --> 'T'
		Afield = np.sqrt((vlocity**2) + (ulocity**2))
		Ylm_o = array_Ylm_o['Ylm_o0']
		f_f = array_f_f['f_f0']
		f_R = array_f_R['f_R0']
		field = array_field['field0']
		u[:,:] = 0.
		v[:,:] = 0.
		Avitesse = Afield

	elif ('M' in getField):					# --> 'M'
		Afield = vlocity
		Ylm_o = array_Ylm_o['Ylm_o0']
		f_f = array_f_f['f_f0']
		f_R = array_f_R['f_R0']
		field = array_field['field0']
		u[:,:] = 0.
		v[:,:] = 0.
		Avitesse = Afield

	elif ('Z' in getField):					# --> 'Z'
		Afield = ulocity
		Ylm_o = array_Ylm_o['Ylm_o0']
		f_f = array_f_f['f_f0']
		f_R = array_f_R['f_R0']
		field = array_field['field0']
		u[:,:] = 0.
		v[:,:] = 0.
		Avitesse = Afield
	#-------END if getfield
	# Azimuthal average
	_Afield_= np.mean(Afield[it,iz,:,:],1)
	_Avitesse_= np.mean(Avitesse[it,iz,:,:],1)
	_u_ = np.mean(u,1)
	_v_ = np.mean(v,1)
	# Plot if iplt==1
	if(iplt==1):
		# --------------------------------------------
		#		PLOT Profils
		# --------------------------------------------
		# ===> Cas m=0
		fig0 = plt.figure()
		plt.plot(lat[:],_Afield_, '-', color='red')
		plt.plot(latSH,Ylm_o, '.', color='red')
		plt.xlabel('latitude')
		plt.ylabel('velocity en m/s')
		#---------------------------------------------
		# ===> Cas m/=0
		plt.plot(lat[:],Afield[it,iz,:,(int(Ratio+len(lon)/2))], '-', color='blue')
		plt.plot(latSH,f_f[:,Iphi], '.', color='blue')
		# 
		plt.plot(lat[:],Afield[it,iz,:,(int(Ratio+len(lon)/2))]-_Afield_, '-', color='black')
		plt.plot(latSH,f_R[:,Iphi], '.', color='black')
		plt.xlabel('latitude')
		plt.ylabel('velocity en m/s')
		plt.title(str(field)+' velocity')
		plt.show()
		# --------------------------------------------
		#		PLOT Maps
		# --------------------------------------------
		cso = plt.contourf((180/math.pi)*(phi), (180/math.pi)*(theta-math.pi/2),f_f,15,cmap=plt.cm.jet)
		colorbar(cso, orientation='vertical')
		plt.title(str(field)+' velocity')
		plt.show()
		#
		cs1 = plt.contourf((180/math.pi)*(phi), (180/math.pi)*(theta-math.pi/2),f_R,15,cmap=plt.cm.jet)
		colorbar(cs1, orientation='vertical')
		plt.title(str(field)+' velocity m/=0')
		plt.show()
		#
		if ('sf' and 'vp' in getField):
			figA = plt.figure()
			cA = plt.contour(lon, lat, Afield[it,iz,:,:],12,color,'k')
			fmt = ticker.LogFormatterMathtext()
			fmt.create_dummy_axis()
			plt.clabel(cA, fontsize=7, inline=1, fmt=fmt)
			plt.title(str(field)+' velocity')
			plt.xlabel('latitude')
			plt.ylabel('longitude')
			plt.show(figA)
		else:
			figA = plt.figure()
			v = np.linspace(-80, 80, 30) # for zonal flow
			#v = np.linspace(-50, 50, 30) # for merid flow
			cA = plt.contourf(lon, lat, Afield[it,iz,:,:],v,cmap=cm.seismic)
			#colorbar(cA, orientation='vertical')
			xxx = plt.colorbar(ticks=v)
			print xxx
			plt.title(str(field)+' velocity')
			plt.xlabel('latitude')
			plt.ylabel('longitude')
			plt.show(figA)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#===================================================================
#		    	     Sum Energy
#===================================================================
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# les polynomes de Legendre on la normalisation geodesique: N = sqrt[(2l+1)/2 * (l-m)!/(l+m)!]
# on retrouve bien l expression de la densite d energie suivante, 
# ----------dans le cas m=0; 
# ici N = sqrt[(2l+1)/2]
# [E] = Int_theta f^2 sin(theta) dtheta = Int_theta | Blm sqrt[(2l+1)/2] Plm |^2 sin(theta) dtheta = 1*|Blm|^2
# La normalisation des polynome de Legendre leur donne la propriete suivante 
# Int_theta sqrt((2l+1)/2) Plm sin(theta) dtheta = 1; pour m=0. avoir [E] = [1/2 * v^2] = m^2/s^2.
# ----------dans le cas m/=0;
# [E] = I_theta Int_phi f^2 sin(theta) dphi dtheta = 
#	I_theta Int_phi Blm*conj(Blm) Ylm*conj(Ylm) sin(theta) dphi dtheta =  
#	2*pi*Blm*conj(Blm)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#			Some quantities
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	deltaT = pas#np.pi/len(_Uzonal_)
	sinTheta = np.sin(np.arange(0,np.pi+deltaT,deltaT))
	deltaP = (2.*np.pi)/len(lon)
	deltaTT = np.pi/len(lat)
	sinTTheta = np.sin(np.arange(0,np.pi,deltaTT))
	ELat_RR= np.zeros([len(lat), len(lon)])
	ELat_TT= np.zeros([len(lat), len(lon)])
	ELat_R= np.zeros([len(x),len(phi)])
	ELat_T= np.zeros([len(x),len(phi)])
	#
	ELat_RR_u= np.zeros([len(lat), len(lon)])
	ELat_RR_v= np.zeros([len(lat), len(lon)])
	ELat_TT_u= np.zeros([len(lat), len(lon)])
	ELat_TT_v= np.zeros([len(lat), len(lon)])
	#
	for iLon in range(0, len(lon)):
		# ------------------------------- if ('sf') or ('vp') alone.
		# On travaille sur m=0 & m/=0 sur u et v independemment, pl-
		# utot que sur VT=sqrt(u^2 + v^2). Le resultat est different.
		if('sf' in getField) and (not 'vp' in getField) or ('vp' in getField) and (not 'sf' in getField):
			ELat_RR_u[:,iLon] = (abs(u[:,iLon]-_u_)**2)*sinTTheta*deltaTT*deltaP
			ELat_RR_v[:,iLon] = (abs(v[:,iLon]-_v_)**2)*sinTTheta*deltaTT*deltaP
			ELat_TT_u[:,iLon] = (abs(u[:,iLon])**2)*sinTTheta*deltaTT*deltaP
			ELat_TT_v[:,iLon] = (abs(v[:,iLon])**2)*sinTTheta*deltaTT*deltaP
			E_MMO_u = (1/2.)*np.sum(((abs(_u_)**2)*sinTTheta*deltaTT))
			E_MMO_v = (1/2.)*np.sum(((abs(_v_)**2)*sinTTheta*deltaTT))
			#
			ELat_RR[:,iLon] = ELat_RR_u[:,iLon]+ELat_RR_v[:,iLon]
			ELat_TT[:,iLon] = ELat_TT_u[:,iLon]+ELat_TT_v[:,iLon]
			E_MMO = E_MMO_u+E_MMO_v
		# ------------------------------- if ('T'), ('M'), ('Z') or ('sf' & 'vp')		
		else:
			ELat_RR[:,iLon] = (abs(Avitesse[it,iz,:,iLon]-_Avitesse_)**2)*sinTTheta*deltaTT*deltaP
			ELat_TT[:,iLon] = (abs(Avitesse[it,iz,:,iLon])**2)*sinTTheta*deltaTT*deltaP
			E_MMO = (1/2.)*np.sum(((abs(_Avitesse_)**2)*sinTTheta*deltaTT))
	
	for iphi in range(0, len(phi)):
		ELat_R[:,iphi] = (abs(f_R[:,iphi])**2)*sinTheta*deltaT*deltaT #ici deltaT=deltaP=pas
		ELat_T[:,iphi] = (abs(f_f[:,iphi])**2)*sinTheta*deltaT*deltaT #ici deltaT=deltaP=pas
		E_MO = (1/2.)*np.sum(((abs(Ylm_o)**2)*sinTheta*deltaT))
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#				 M=0
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ===> Cas m=0	
	print '================================================================'
	print '                 Field :',field,'    '	
	print '================================================================'
	print '-------------------------------------------------------------'
	print '-------------------ENERGY, m=0-------------------------------'
	print '-------------------------------------------------------------'
	print '-------------------------------------- PROFIL m=0 TO COMPARE:'
	print 'E_mo =',E_MMO
	print '----------------------- profil SpHarm: si (T-M or Z);(L >> 0)'# Ylm_o doit etre transf en vitesse si sf.
	print 'E_mo =',E_MO
	print '--------------------------------------------- SPECTRAL COEFFS'
	print 'E_mo =',np.sum(Eno)
	print '--------------------------- profil SpHarm: si ( sf );(L >> 0)'# m=0 of sf is diff of azim_mean(Vtotale).
	print 'E_mo =',(1/2.)*np.sum(((abs(u_sf_mo)**2)*sinTheta*deltaT)) # + Valeur un peu inf: derivee approximee ds sf2V.
	print ' '
	print ' '
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#				M/=0
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ===> Cas m/=0	
	print '-------------------------------------------------------------'
	print '-------------------ENERGY, m/=0-------------------------------'
	print '-------------------------------------------------------------'	
	print '------------------------------------ PROFIL NCDEF TO COMPARE:'
	print 'E_R =',(1/(2.*2.*np.pi))*np.sum(ELat_RR) #!!!si sf: residus Vitesse total uv_sf st diff des residus de sf.
	print 'E_T =',(1/(2.*2.*np.pi))*np.sum(ELat_TT),'(E_mo=',E_MMO,')'
	print '------------------- profil SpHarm: si (T-M or Z);(M & L >> 0)'
	print 'E_R =',(1/(2.*2.*np.pi))*np.sum(ELat_R)
	print 'E_T =',(1/(2.*2.*np.pi))*np.sum(ELat_T),'(E_mo=',E_MO,')'
	print '--------------------------------------------- SPECTRAL COEFFS'
	print 'E_R =',np.sum(En)
	print 'E_T =',np.sum(Eno) + np.sum(En),'(E_mo=',np.sum(Eno),')'
	print 'End'
	print '================================================================'
	E_mo = np.sum(Eno)
	E_R = np.sum(En)
	E_T = np.sum(Eno) + np.sum(En)

	URMS = np.sqrt(np.mean(np.square(Avitesse[it,iz,:,:])))
	#URMS = np.mean(Avitesse[it,iz,:,:])
	#URMS = np.sqrt(np.sum(ulocity[it,iz,150:200,:]**2)/(len(x)*len(phi)))
	print 'URMS = ',URMS
	return (E_mo, E_R, E_T, URMS, latSH, field)

