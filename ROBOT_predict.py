#!/usr/bin/python

#!/usr/bin/env python
""" ROBOT_predict.py -- A script for predicting the Real/Bogus score of a template/science/subtraction triplet with ROBOT.

Usage: ROBOT_predict.py [-h] [-v] [-q] [--debug] [-m MODEL] <template> <science> <subtraction> <ra> <dec

Arguments:
    template (string)
        path to a FITS or PNG template file

    science (string)
        path to a FITS or PNG science file

    subtraction (string)
        path to a FITS or PNG subtraction file
	
    ra (float)
    	Right Ascension of target (decimal format)
	
    dec (float)
    	Declination of target (decimal format)
    	
Options:
    -h, --help                              Show this screen
    -v, --verbose                           Show extra information [default: False]   
    -q, --quietmode                         Minimize print to screen. This is useful when this function is called in another function. [default: False]  
    --debug                                 Output more for debugging [default: False]
    -m MODEL, --model MODEL                 Preloaded Model (load_ROBOT_model.py) [default: None]

Examples:
    bash: ROBOT_predict.py -v temp.fits sci.fits sub.fits 100.5 42.0
"""



import docopt
import numpy as np
import sys, os
import matplotlib.pyplot as plt
from tensorflow import keras
from keras.models import Sequential, Model
from keras.layers import Activation, Dense, Dropout, Flatten, Input, Concatenate, Conv2D, MaxPooling2D, GlobalAveragePooling2D, AveragePooling2D
import keras.backend as K
import astropy
from astropy.wcs import WCS
from astropy.io import fits
from astropy.nddata.utils import Cutout2D

__author__	= "Simon Goode"
__license__	= "MIT"
__version__	= "0.1"
__date__	= "2022-01-11"
__maintainer__	= "Simon Goode"
__email__	= "sgoode@swin.edu.au"

#########################################
# ====== House Keeping Functions ====== #
#########################################
'''These functions standardise verbose, debug printing'''
def print_verbose_string(printme,verbose=False,underscores=False):
    if verbose:
        if underscores:
            print("@" * len(f"VERBOSE: {printme}"),file=sys.stdout)
            print(f"VERBOSE: {printme}",file=sys.stdout)
            print("@" * len(f"VERBOSE: {printme}"),file=sys.stdout)
        else:
            print(f"VERBOSE: {printme}",file=sys.stdout)

def print_debug_string(printme,debugmode=False,underscores=False):
    if debugmode:
        if underscores:
            print("@" * len(f"DEBUG  : {printme}"),file=sys.stdout)
            print(f"DEBUG  : {printme}",file=sys.stdout)
            print("@" * len(f"DEBUG  : {printme}"),file=sys.stdout)
        else:
            print(f"DEBUG  : {printme}",file=sys.stdout)

'''These functions help organise or delete files'''
def clearit(fname):
    if os.path.isfile(fname):
        os.remove(fname)
    return None


#########################################
# ====== Supplementary Functions ====== #
#########################################

def EdgeCheck(img):
	if len(np.where([-.01<i<.01 for i in img.flatten()])[0]) >= 62:
		# Raise a flag if img is close to a CCD edge
		return True
	else:
		return False

#########################################
# =========== Main Function =========== #
#########################################

def ROBOT_predict(model, template, science, subtraction, ra, dec, verbose=False, debugmode=False):
	if model=='None':
		from load_ROBOT_model import load_ROBOT_model
		m = '/fred/oz100/pipes/DWF_PIPE/ROBOT_pipe/CONFIG/exp_35_remade.h5'
		w = '/fred/oz100/pipes/DWF_PIPE/ROBOT_pipe/CONFIG/exp_35_remade_weights.h5'
		model = load_ROBOT_model(m, w, verbose=verbose, debugmode=debugmode)
	
	arr = np.empty((1,31,31,3))
	norm_vals = [25.759, 7.685, 25.762]
	for i,f in enumerate([science, subtraction, template]):
		if f.endswith('.fits'):
			# Do FITS analysis
			try:
				with fits.open(f) as hdu:
					w = WCS(hdu[0].header)
					head = hdu[0].header
					xlim=head['NAXIS1']
					ylim=head['NAXIS2']

					pixcrd_im = np.array([[xlim, ylim]], np.float_)
					world_im = w.wcs_pix2world(pixcrd_im, 1)
					pixx_im, pixy_im = world_im[0][0], world_im[0][1]				
				
					corners=w.calc_footprint()
					corner_1 = corners[0]
					corner_2 = corners[1]
					corner_3 = corners[2]
					corner_4 = corners[3] 
					differnce = corner_1 - corner_2 

					pixcrd = np.array([[ra, dec]], np.float_)
					worldpix = w.wcs_world2pix(pixcrd, 1)
					pixx, pixy = worldpix[0][0], worldpix[0][1]
			
					try:
						cutout = Cutout2D(hdu[0].data, (pixx, pixy), 31, wcs= w)
						if i==2 or i==0:
							edge_flag = EdgeCheck(cutout.data)
							if edge_flag:
								return 0.
							
						arr[:,:,:,i] = cutout.data/norm_vals[i]
			
					except astropy.nddata.utils.NoOverlapError:
						hdu.close()
						return 0.
				hdu.close()
			except OSError:
				print_verbose_string(f'OSError caught: Empty or Corrupt FITS file. Returning NaN.', verbose=verbose)
				return np.nan
		else:
			return 'Filetype not recognised {f}'
			
		
	return model.predict(arr)[0][0]
		

		

############################################################################
####################### BODY OF PROGRAM STARTS HERE ########################
############################################################################

if __name__ == "__main__":

    # Read in input arguments
    arguments       = docopt.docopt(__doc__)
    verbose         = arguments['--verbose']
    quietmode       = arguments['--quietmode']
    debugmode       = arguments['--debug']    
    model           = arguments['--model']

    template        = arguments['<template>']
    science         = arguments['<science>']
    subtraction     = arguments['<subtraction>']
    
    if debugmode:
        print(arguments)  

    _ = ROBOT_predict(model, template, science, subtraction, verbose=verbose, debugmode=debugmode)
