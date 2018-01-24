#! /usr/bin/env python
from netCDF4 import Dataset as NetCDFFile 
import numpy as np
import matplotlib.pyplot as plt
from windspharm.standard import VectorWind
from windspharm.tools import prep_data
from windspharm.tools import get_recovery
import time
import sys
# Cette routine permet d ajouter la stream-function et le potentiel
# de vitesse aux autres composantes dans un unique fichier appele
# sfvpData.nc. Windspharm permet d obtenir sf et vp. On selectionne
# un seul layer iz en profondeur.
####################################################################
# 		    PTS: "Parameters to set"
####################################################################
zefile = 'DRAG90days_DISSIP10000_year1-10_uv_z5_512_every200d.nc' # netcdef file to load
time_range=[1000,1010] # ex: if [0,10] selected [0,9] is finally implemented;
iz=0 # select a layer on the vertical index iz.
R_sat = 58.232*10.**(6.) # en m

zefile = '/scratch/cnt0026/lmd7357/aspigaheat/reference/saturn/simu_REF_10days/Xhistins_28.nc'
#zefile = '/scratch/cnt0026/lmd7357/aspigaheat/reference/saturn/simu_REF_10days_graddiv5000_order1/Xhistins_28.nc'
zefile = '/store/aspigaheat/dynamico-giant/saturn/simu_REF3_10years/Xhistins_80.nc'
zefile = '/scratch/cnt0026/lmd7357/aspigaheat/reference/saturn/Xhistins_80.nc'
zefile = '/scratch/cnt0026/lmd7357/aspigaheat/reference/saturn/simu_REF_10days_graddiv5000_order1/Xhistins_1.nc'
time_range=[0,50]
iz=5
####################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 			  Data
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##-----------------
nc = NetCDFFile(zefile) # load file

# get dimensions from netCDF file
zedim = nc.dimensions.keys()
xcoord = zedim[1]
ycoord = zedim[2]
zcoord = zedim[3]
tcoord = zedim[4]
print tcoord,zcoord,ycoord,xcoord

nbiz=1 # On est pour le moment limite a un seul iz.
# On forme la martice presniv_AS qui n a qu une seule valeure. 
# C est pour quelle est une longueure.
presnivs_AS = np.zeros([nbiz])
#
print '----------------------------------------------------------------------------------'
print '----- Initial file uploaded'
for varname in nc.variables.keys():
	var = nc.variables[varname]
	print varname, var.dtype, var.dimensions, var.shape
print '----------------------------------------------------------------------------------'
#
lat_AS = nc.variables[ycoord][:] #(128,)
lon_AS = nc.variables[xcoord][:]
presnivs_AS[nbiz-1] = nc.variables[zcoord][iz]
print nc.variables[zcoord][iz]
#
u_AS = nc.variables['u'][time_range[0]:time_range[1],iz,:,:] # [t,Prof,theta,phi] Zonal velocity
v_AS = nc.variables['v'][time_range[0]:time_range[1],iz,:,:] # [t,Prof,theta,phi] meridional velocity
#-----------------------------------------------------------------------
#			WRITE NETCDF
#-----------------------------------------------------------------------
# creer un fichier netcdef --> http://www.ceda.ac.uk/static/media/uploads/ncas-reading-2015/11_create_netcdf_python.pdf
dataset = NetCDFFile('sfvpData-'+str(time_range[0])+'-'+str(time_range[1]-1)+'.nc', 'w', format='NETCDF3_CLASSIC')
	#----Dimensions-----------------------------------------------------------------------------
presnivs = dataset.createDimension(zcoord, len(presnivs_AS))
presnivss = dataset.createVariable(zcoord, np.float64, (zcoord,))
presnivss[:] = presnivs_AS
	#
longitude = dataset.createDimension(xcoord, len(lon_AS))
longitudes = dataset.createVariable(xcoord, np.float64, (xcoord,))
longitudes[:] = lon_AS
	#
latitude = dataset.createDimension(ycoord, len(lat_AS))
latitudes = dataset.createVariable(ycoord, np.float64, (ycoord,))
latitudes[:] = lat_AS  
	#
time_counter = dataset.createDimension(tcoord, None)
	#
v = dataset.createVariable('v', np.float64, (tcoord,zcoord,ycoord,xcoord))
u = dataset.createVariable('u', np.float64, (tcoord,zcoord,ycoord,xcoord))
streamfunction = dataset.createVariable('sf', np.float64, (tcoord,zcoord,ycoord,xcoord))
velocitypotential = dataset.createVariable('vp', np.float64, (tcoord,zcoord,ycoord,xcoord))
	#-----Data-----------------------------------------------------------------------------------
	# Create a VectorWind instance to handle the computation of streamfunction and
	# velocity potential.
	# Compute le champ de vitesse compatible avec la fonction vectorfield.
u_prep, info = prep_data(u_AS, 'tyx')#'tzyx' prepare les donnees pour vectrofield
v_prep, info = prep_data(v_AS, 'tyx')#'tzyx' prepare les donnees pour vectrofield
del u_AS, v_AS
w = VectorWind(u_prep, v_prep,rsphere=R_sat)
sf, vp = w.sfvp()
recover = get_recovery(info) # resize to the initial shape
uo, vo, sfo, vpo = recover(u_prep, v_prep, sf, vp)
del u_prep, v_prep, sf, vp
print uo.shape, sfo.shape
v[:,nbiz-1,:,:] = vo
u[:,nbiz-1,:,:] = uo
streamfunction[:,nbiz-1,:,:] = sfo
velocitypotential[:,nbiz-1,:,:] = vpo

dataset.close()

print '///////////////////',dataset.file_format,'named: sfvpData.nc writen///////////////' 

print '----------------------------------------------------------------------------------'
print '------ Created file netcdf'
ncData = NetCDFFile('sfvpData-'+str(time_range[0])+'-'+str(time_range[1]-1)+'.nc')
for varname in ncData.variables.keys():
	var = ncData.variables[varname]
	print varname, var.dtype, var.dimensions, var.shape
print '----------------------------------------------------------------------------------'

####################################################################
# 			OUTPUTS
####################################################################
# sfvpData-xxx-xxx.nc : 
# A file netcdet that contain velocities, u,v and sf, vp the stream-
# function and the velocity potential. 
