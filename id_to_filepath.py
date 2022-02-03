#!/usr/bin/python

#!/usr/bin/env python
""" id_to_filepath.py -- A script for finding a detection fits file with a detection index. 

Usage: id_to_filepath.py [-h] [-v] [-q] [--debug] <file> <id> <listpath> <impath>

Arguments:
    file (string)
    	Path to a light curve file (fmt: MJD, g_mag, g_mag_err, upper_lim_mag (space-separated))
    
    id (int)
    	Index of SMOTE in the light curve file
	
    listpath (string)
    	Path to an ordered image list

    impath (string)
	Path to ccds for the appropriate field / date of observations
	
Options:
    -h, --help                              Show this screen
    -v, --verbose                           Show extra information [default: False]   
    -q, --quietmode                         Minimize print to screen. This is useful when this function is called in another function. [default: False]  
    --debug                                 Output more for debugging [default: False]

Examples:
    bash: id_to_filepath.py /fred/oz100/NOAO_archive/archive_NOAO_data/data_outputs/2015/01/4hr/g_band/single/lightcurves/files/DWF040748.225-545713.829_150116 61 /fred/oz100/NOAO_archive/archive_NOAO_data/scripts/create_lc/image_mjd_lists/FINAL_2015_01_4hr_150116.ascii /fred/oz100/NOAO_archive/archive_NOAO_data/data_outputs/2015/01/4hr/g_band/single/
"""



import docopt
import os, sys
import pandas as pd
import numpy as np
import math
import pickle as p
from astropy.io import fits
from astropy.wcs import WCS

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
'''These functions are used within the Main function'''
def RAsex_to_RAdec(fRAsex):
	frah = float(fRAsex[0:2])
	fram = float(fRAsex[2:4])
	fras = float(fRAsex[4:])
	return ((1/1 * frah) + (1/60 *fram) + (1/3600 *fras))* (360/24)
    
def DEsex_to_DEdec(fDEsex):
	fded = float(fDEsex[0:3])
	fdem = float(fDEsex[3:5])
	fdes = float(fDEsex[5:])    
	fDEdec = (math.fabs(fded)*3600.0+fdem*60.0+fdes)/3600.0
	if fDEsex[0] == '-':
		fDEdec = fDEdec * -1
	return fDEdec

def SearchCCDs(impath, ordered_ims, det_id, RA, DEC):
	for n, im in enumerate(ordered_ims):
		if n != det_id:
			continue
		im_path = f'{impath}{im[:26]}/ccds/'
		for fitsim in os.listdir(im_path):
			with fits.open(im_path+fitsim) as hdu:
				w = WCS(hdu[0].header)
				corners=w.calc_footprint()
				corner_1 = corners[0]
				corner_2 = corners[1]
				corner_3 = corners[2]

				if  corner_1[0] <= RA <=corner_2[0] and corner_1[1] >= DEC >= corner_3[1]:
					return im_path+fitsim

#########################################
# =========== Main Function =========== #
#########################################

def id_to_filepath(filepath, det_id, listpath, impath, verbose=False, debugmode=False):
	fname = filepath.split('/')[-1]
	print_debug_string(f'Selecting templates for {fname}', debugmode=debugmode)
	lc = pd.read_csv(filepath, delimiter=' ', header=0, error_bad_lines=False)
	
	day = fname.split('_')[1]
	ordered_ims = np.loadtxt(listpath, skiprows = 1, usecols=[0], dtype= str)
	RA = RAsex_to_RAdec(fname[3:13])
	DEC = DEsex_to_DEdec(fname[13:24])
	fitspath = SearchCCDs(impath, ordered_ims, det_id, RA, DEC)
	return fitspath
	
############################################################################
####################### BODY OF PROGRAM STARTS HERE ########################
############################################################################

if __name__ == "__main__":

    # Read in input arguments
    arguments       = docopt.docopt(__doc__)
    verbose         = arguments['--verbose']
    quietmode       = arguments['--quietmode']
    debugmode       = arguments['--debug']
    
    filepath        = arguments['<file>']
    det_id          = arguments['<id>']
    listpath        = arguments['<listpath>']
    impath          = arguments['<impath>']

    if debugmode:
        print(arguments)  

    _ = id_to_filepath(filepath, int(det_id), listpath, impath)
    with open('out.p', 'wb') as f:
        p.dump(_, f)
