use ```git clone --recursive``` to get submodules windspharm and pyspharm

-------------------------------------------------------------------------

to install

```
mkdir $HOME/python_setup 
cd pyspharm
python setup.py install --home=$HOME/python_setup
cd ..
cd windspharm
python setup.py install --home=$HOME/python_setup
cd ..
```

in your ```.bashrc``` (or else) environment file add
```
export PYTHONPATH=$PYTHONPATH:$HOME/python_setup/lib/python/
export PYTHONPATH=$PYTHONPATH:$HOME/python_setup/lib/python/windspharm-1.6.dev0-py2.7.egg
```

-------------------------------------------------------------------------




- Statistical analysis of a 2D winds layer. Compute gauss coefficients
from projections on spherical harmonic functions. -

------------------------------------------------------------------
------------------------------------------------------------------
# It contains:
	executable:		sfvp_t.py	
				spectra_analysis_scalar_01Oct
				spectra-zono-01Oct.py
				spectra-zono-temp-01Oct.py

	functions		getKES.py
				SpHarm.py
				getSumEnergy.py
				sf2V.py

	parameters file		filePTS.zono
				filePTS.zono.temp

# It requires to install the following library:
	windspharm:	http://ajdawson.github.io/windspharm/

------------------------------------------------------------------
------------------------------------------------------------------

Start: 

 --> "sfvp_t.py"

# Use "sfvp_t.py" to extract wind field from a netcdf file containing
"latitude", "longitude", "presnivs", "u", "v".
# Select in Parameters To Set: 1) the netcdf file nc=NetCDFFile('xxx.nc')
2) the time_range[] you want to compute, according the time iteration 
available in xxx.nc file. 3) select a layer on the vertical index iz.
iz=0 if one layer only. 4) select R_sat, radius in meter of your planet.
# You obtain a new netcdf file named "sfvpData-tr1-tr2.nc" with tr1
and tr2 your time range, containing the stream-function and velocity
potential.

 --> "spectra_analysis_scalar_01Oct"

# Use "spectra_analysis_scalar_01Oct" to extract Gauss coefficients
on the basis of spherical harmonic functions.
# This fortran exe has arguments detailed with "-h" argument. however
the following are commonly used: 
"spectra_analysis_scalar_01Oct -field x1 -istp x2 -tstp x3 -mt x4 FILE"
- x1 : select the field to analyse; x1= 1 (zonal flow) 2 (merid flow) 
3 (total velocity) 4 (stream-function) 5 (velocity potential)
- x2 : initial time to compute from "sfvpData--.nc" (default: 0)
- x3 : time step to compute from "sfvpData--.nc" (default: 1)
- x4 : final time to compute from "sfvpData--.nc" (default: 0)
# It results files "xspec_x1_xx" created at each time step xx and for 
the selected field x1.

 --> "spectra-zono-01Oct.py"

# From files "xspec_x1_xx" we plot spectra of the axisymmetric and
non-axisymmetric component of the flow and we compute integrated 
energetic value over the whole layer surface of these two flow 
components.
# To run "spectra-zono-01Oct.py" one has to set parameters in the
file "filePTS.zono".

 --> "spectra-zono-temp-01Oct.py"

# From files "xspec_x1_xx" we plot time evolution spectra of the 
axisymmetric and non-axisymmetric component of the flow and we  
compute integrated energetic value over the whole layer surface 
of these two flow  components.
# To run "spectra-zono-temp-01Oct.py" one has to set parameters in the
file "filePTS.zono.temp".
