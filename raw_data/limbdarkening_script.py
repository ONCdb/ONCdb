#! /usr/bin/env python3

import os, sys
sys.path.append('/Users/gbruno/python/source/ExoCTK/')
from ExoCTK import core
from ExoCTK.ldc import ldcfit as lf
from ExoCTK.ldc import ldcplot as lp
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy as np
import pdb

# Path for the models
#fits_files = '/Users/gbruno/idl/projects/limb_test/phxinten/'
fits_files = '/Users/gbruno/python/projects/ExoCTK/kurucz/'
joe_fits_files = '/user/jfilippazzo/Models/ATLAS9/'
teff, logg, FeH = 5250, 4.0, 0.0

def sing_v_ldc():
    model_grid = core.ModelGrid(joe_fits_files)
    K_band = core.Filter('Kepler.K')

    # Plot LDC
    x = lf.ldc(teff, logg, FeH, model_grid, ['quadratic', '4-parameter'], mu_min=0.00, bandpass=K_band)
    
    # Plot Sing
    mu = x['mu']
    c1, c2, c3, c4 = 0.71364387, -0.68152910, 1.3952835, -0.62918712
    a, b = 0.49340705, 0.19758911
    sing4 = 1. - c1*(1. - mu**0.5) - c2*(1. - mu) - c3*(1. - mu**1.5) - c4*(1. - mu**2)
    sing2 = 1. - a*(1. - mu) - b*(1. - mu)**2
    plt.plot(mu, sing2, 'g-', label = 'Sing quadratic')
    plt.plot(mu, sing4, 'r-', label = 'Sing 4 coeffs')

def ld_bins():

    bins = np.arange(1.125, 1.650, 0.00453128*9)
    # Stellar parameters
    #teff, logg, FeH = 2500, 4.5, 0.0
    # Initialize grid
    model_grid = core.ModelGrid(fits_files)

    K_band = core.Filter('Kepler.K')
    plt.plot(*K_band.rsr)

    # Open text file to save results
    fname = 'test1_kurucz_ldc_fix_Kband.dat'
    ldfile = open(fname, 'w')
    ldfile.write('bin_i\tbin_f\tu1\tu2\n')
    ldfile.write('-----\t-----\t--\t--\n')
    for i in range(1, len(bins)):
        # Compute LD coefficients
        print(i)
        model_grid.customize(Teff_rng=(5000,5500), logg_rng=(2.0,4.5), FeH_rng=(-.5,.5), wave_rng=(bins[i - 1], bins[i]))
        print(model_grid.data)
        #model_grid.customize(Teff_rng=(2000,3000), logg_rng=(4.0,5.0), FeH_rng=(-0.5,0.5), wave_rng=(bins[i - 1], bins[i]))
        interpolation = lf.ldc(teff, logg, FeH, model_grid, '4-parameter', mu_min = 0.00, plot=True)
        plt.close('all')
        coeffs = interpolation['4-parameter']['coeffs']
        #coeffs, mu, radius = lf.ldc(teff, logg, FeH, model_grid, '4-parameter', plot=False)

        print(bins[i - 1], bins[i], coeffs[0], coeffs[1], coeffs[2], coeffs[3])

        # Write to file
        #ldfile.write('%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n' %(bins[i - 1], bins[i], coeffs[0], coeffs[1], coeffs[2], coeffs[3]))
    # Whole range
    model_grid.customize(Teff_rng=(5000,5500), logg_rng=(2.0, 4.5), FeH_rng=(-0.5,0.5), wave_rng=(0.3, 1.0))
    #coeffs, mu, radius = lf.ldc(teff, logg, FeH, model_grid, '4-parameter', plot=False, bandpass = K_band)
    style = ['--', '-']
    for ii, laws in enumerate(['quadratic', '4-parameter']):
        interpolation = lf.ldc(teff, logg, FeH, model_grid, laws, ls = style[ii], plot=True, mu_min = 0.00, bandpass = K_band)
        coeffs = interpolation[laws]['coeffs']
        print('0.3', '1.0', coeffs, interpolation[laws]['err'])#, coeffs[1], coeffs[2], coeffs[3])
        #print('0.3', '1.0', coeffs[0], coeffs[1], coeffs[2], coeffs[3])
        #ldfile.write('%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n' %(float('0.3'), float('1.0'), coeffs[0], coeffs[1], coeffs[2], coeffs[3]))
    ldfile.close()
    mu = interpolation['mu']
    c1, c2, c3, c4 = 0.71364387, -0.68152910, 1.3952835, -0.62918712
    a, b = 0.49340705, 0.19758911
    sing4 = 1. - c1*(1. - mu**0.5) - c2*(1. - mu) - c3*(1. - mu**1.5) - c4*(1. - mu**2)
    sing2 = 1. - a*(1. - mu) - b*(1. - mu)**2
    plt.plot(mu, sing2, 'g-', label = 'Sing quadratic')
    plt.plot(mu, sing4, 'r-', label = 'Sing 4 coeffs')
    plt.xlabel('$\mu$', fontsize = 18)
    plt.ylabel('$I(\mu)/I(\mu = 1)$', fontsize = 18)
    plt.title('Teff = ' + str(teff) + ', log g = ' + str(logg) + ', [Fe/H] = ' + str(FeH))
    plt.legend(loc = 'best')
    plt.savefig('LDCvsSing.png')

    return
