#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Last modification: 10/2021
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# MODULES

import bxa.xspec as bxa
import xspec as xs
import matplotlib.pyplot as plt
import pdb
import numpy as np
import os
from def_model_dictionary import *

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def plot_qq( solver, path_bxa ):
    '''Quantile plot from BXA documentation & examples.
    '''

    plt.figure( figsize = ( 7 , 7 ) )

    with bxa.XSilence():
    	solver.set_best_fit()
    	bxa.qq.qq( prefix = path_bxa, markers=5, annotate = True)

    plt.savefig( path_bxa + 'plots/qq_model_deviations.pdf', bbox_inches = 'tight' )
    plt.close()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def plot_posterior_predictions( solver, path_bxa ):
    '''Posterior predictions plot from BXA documentation & examples.
    '''
    plt.figure()
    data = solver.posterior_predictions_convolved( nsamples = 100 )
    binned = bxa.binning( outputfiles_basename = path_bxa, \
                          bins = data[ 'bins' ], widths = data[ 'width' ], \
                          data = data[ 'data' ], models = data[ 'models' ])

    for point in binned[ 'marked_binned' ]:
    	plt.errorbar( marker = 'o', zorder =-1, **point )
    plt.xlim( binned[ 'xlim' ])
    plt.ylim( binned[ 'ylim' ][0], binned[ 'ylim' ][1] * 2)
    plt.gca().set_yscale( 'log' )

    if xs.Plot.xAxis == 'keV':
    	plt.xlabel( 'Energy [keV]' )
    elif xs.Plot.xAxis == 'channel':
    	plt.xlabel('Channel')

    plt.ylabel( 'Counts/s/cm$^2$' )

    plt.savefig( path_bxa + 'plots/convolved_posterior.pdf', bbox_inches = 'tight' )
    plt.close()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def plot_best_fit_model_components( results, model_name, path_plots, Emin = 0.5, Emax = 7 ):
    '''Plot best fit model and its components.
    '''


    with bxa.solver.XSilence():

        declare_sqrtcutoffpl()

        # xspec model
        xs.AllModels.clear()
        model = xs.Model( model_name )

        paramnames = []
        for k in range( 1, model.nParameters + 1 ):
            pn = model.__call__(k).name
            paramnames.append( pn )

        for i, par in enumerate(paramnames):
            i += 1
            for j, bxapar in enumerate(results['paramnames']):

                pval = results['posterior']['mean'][j]

                if bxapar[:3] == 'log':
                    bxapar = bxapar[4:-1]
                    pval = np.power(10, results['posterior']['mean'][j])

                if bxapar == par:
                    model(i).values = pval
                    results['paramnames'].pop(j)
                    results['posterior']['mean'] = np.delete( results['posterior']['mean'], j )
                    break

        # write BXA fit to qdp
        xs.Plot.commands   = ()
        xs.Plot.device     = "/null"
        xs.Plot.add        = False
        xs.Plot.xLog       = True
        xs.Plot.yLog       = True
        xs.Plot.xAxis      = 'keV'
        cwd = os.getcwd()
        os.chdir( path_plots )
        qdpf = 'modeldp.qdp'
        xs.Plot.addCommand('wd %s' %(qdpf) )
        xs.Plot('model')

        # read qdp file
        e     = []
        de    = []
        tot   = []
        ncomp = len(model.componentNames) - 1
        if ncomp > 1:
            comps = [ [] for _ in range(ncomp) ]

        with open( qdpf, 'r' ) as file:
            for i in range(3):
                file.readline() # Ignore header.

            for line in file:
                split = line.split()
                if split[0] == 'NO':
                    print('plot_spectrum: met <<NO>> I stop reading.')
                    break

                if ( np.float(split[0]) > Emin ) & ( np.float(split[0]) < Emax ):
                    e.append(np.float(split[0]))
                    de.append(np.float(split[1]))
                    tot.append(np.float(split[2]))

                    try:
                        if ncomp > 1:
                            for i, c in enumerate(comps):
                                c.append(np.float(split[3+i]))
                    except:
                        pdb.set_trace()

        os.remove(qdpf)
        os.chdir(cwd)

        # plot model
        fig, ax = plt.subplots(1)
        ax.plot( e, tot, color = 'black', label = '%s' %(model_name), linewidth = 1, alpha = 0.8 )
        if ncomp > 1:
            for i, c in enumerate(comps):
                ax.plot( e, c, label = '%s' %(model.componentNames[i+1]), linewidth = 1, alpha = 0.8 )

        ax.set_xlabel('Energy (keV)')
        ax.set_ylabel(r'counts/cm$^2$/s/keV')
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.grid( which = 'major', linewidth = 1, alpha = 0.2)
        ax.legend()
        ax.set_ylim(bottom = 1E-7)

        plt.savefig(os.path.join(path_plots, 'model.pdf'), format = 'pdf')
        plt.close()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if __name__ == '__main__':

    # modules
    import os
    import json
    import glob
    from def_model_dictionary import *

    # Paths, lists & variables
    path_scripts  = '/home/ellien/Tycho/scripts/test'
    path_spectra  = '/home/ellien/Tycho/data/10095/spectra_pub'
    path_data     = '/home/ellien/Tycho/data/10095/repro'
    path_analysis = '/home/ellien/Tycho/analysis_pub/out2'

    model_dir = model_dictionnary()

    # iterate over all models
    for path_bxa in glob.glob( os.path.join( path_analysis, 'acis*2*' )):

        print(path_bxa)

        path_plots = os.path.join( path_bxa, 'plots' )
        path_results = os.path.join( path_bxa, 'info/results.json' )

        if os.path.exists(path_results):

            # Read results
            file = open( path_results )
            results = json.load( file )

            # plot
            key = path_bxa.split('_')[-1]
            model_name = model_dir[key]
            plot_best_fit_model_components( results, model_name, path_plots )
