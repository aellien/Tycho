#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Last modification: 07/2021
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# MODULES

import xspec as xs
import bxa.xspec as bxa
import numpy as np
import matplotlib.pyplot as plt
import os
import ray
from datetime import datetime
from def_plot_qq_post import plot_posterior_predictions, plot_qq

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

startTime = datetime.now()

# paths, lists & variables
#path_scripts = '/home/ellien/Tycho/scripts'
#path_spectra = '/home/ellien/Tycho/data/10095/spectra_pub'
#path_data    = '/home/ellien/Tycho/data/10095/repro'
#
#path_bxa     = '/home/ellien/Tycho/analysis_pub/acis_bxa_region_2_sqrtcutoffplpshock/'

path_scripts = '/home/ellien/Tycho/scripts'
path_spectra = '/home/ellien/Tycho/data'
path_data    = '/home/ellien/Tycho/data'

path_bxa     = '/n03data/ellien/Tycho/analysis_pub/acis_bxa_region_2_pownei/'

# xspec model
xs.AllModels.clear()
model_name = 'TBabs(powerlaw+nei)'
#_______________________________________________________________________________
# par  comp
#   1    1   TBabs      nH         10^22    1.00000      +/-  0.0
#   2    2   powerlaw   PhoIndex            1.00000      +/-  0.0
#   3    2   powerlaw   norm                1.00000      +/-  0.0
#   4    3   nei        kT         keV      1.00000      +/-  0.0
#   5    3   nei        Abundanc            1.00000      frozen
#   6    3   nei        Tau        s/cm^3   1.00000E+11  +/-  0.0
#   7    3   nei        Redshift            0.0          frozen
#   8    3   nei        norm                1.00000      +/-  0.0
#_______________________________________________________________________________

list_input_par  = [ [     0.7, 0.001,      0,      0,    2.0,   2.0 ],
                [           3,  0.01,      1,      1,      5,     5 ],
                [       1e-02,  0.01,  1e-05,  1e-05,  1e-01, 1e-01 ],
                [         1.0,  0.01,    0.1,    0.1,  3e+00, 3e+00 ],
                [           1,  0.01,    0.1,    0.1,  1e+00, 1e+00 ],
                [       1e+08, 1e+08,  1e+08,  1e+08,  5e+09, 5e+09 ],
                [           0,  0.01, -0.999, -0.999,     10,    10 ],
                [       1e-03,  0.01,  1e-04,  1e-04,  5e-03, 5e-03 ] ]

model = xs.Model( model_name )

for k in range( 1, model.nParameters + 1 ):
    model(k).values = list_input_par[ k - 1 ]

ncomp = len( model.componentNames )

# xspec xset parameters
xs.Xset.abund = "angr"
xs.Xset.cosmo = "70 0 0.73"
xs.Xset.xsect = "bcmc"

# xspec fitting parameters
xs.Fit.statMethod  = 'cstat'
xs.Fit.method      = "leven 1000 0.01"
#xs.Fit.query       = "no"

os.chdir( path_spectra ) # Necessary to bypass the abscence of absolute path in
                         # anscillary file names in primary file header.

xs.AllData -= "*" # Clear all previous spectra.

# xspec read spectrum
infilename = 'spec_region_2.pi'
infilepath = os.path.join( path_spectra, infilename )

spectrum = xs.Spectrum( infilepath )
os.chdir( path_scripts ) # go back to script directory.

spectrum.ignore( "**-0.5,7.-**" )
#SPECTRUM.ignore( "1-35,480-1024" )

# bxa priors
transformations = []
for k in range( 1, model.nParameters + 1 ):
    if k in [ 3, 5, 6, 8 ]:
        transformations.append( bxa.create_loguniform_prior_for( model, model(k) ) )
    elif k in [ 7 ]:
        pass
    else:
        transformations.append( bxa.create_uniform_prior_for( model, model(k) ) )

# bxa solver
solver = bxa.BXASolver( transformations = transformations, outputfiles_basename = path_bxa )
results = solver.run( resume = 'overwrite', log_dir = os.path.join( path_bxa, 'logs' ) )

plot_qq( solver, path_bxa )
plot_posterior_predictions( solver, path_bxa )

print(datetime.now() - startTime)
