#!/usr/bin/python

#!/usr/bin/env python
""" isolate_SMOTE_detection.py -- A script for finding the detection index of a SMOTE. This script assumes that the SMOTE is the brightest detection in a light curve.

Usage: isolate_SMOTE_detection.py [-h] [-v] [-q] [--debug] <file>

Arguments:
    file (string)
    	Path to a light curve file (fmt: MJD, g_mag, g_mag_err, upper_lim_mag (space-separated))
	
Options:
    -h, --help                              Show this screen
    -v, --verbose                           Show extra information [default: False]   
    -q, --quietmode                         Minimize print to screen. This is useful when this function is called in another function. [default: False]  
    --debug                                 Output more for debugging [default: False]

Examples:
    bash: isolate_SMOTE_detection.py /fred/oz100/NOAO_archive/archive_NOAO_data/data_outputs/2015/01/4hr/g_band/single/lightcurves/files/DWF040748.225-545713.829_150116
"""



import docopt
import os, sys
import pandas as pd
import numpy as np

__author__	= "Simon Goode"
__license__	= "MIT"
__version__	= "1.0"
__date__	= "2022-01-07"
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

	
#########################################
# =========== Main Function =========== #
#########################################

def isolate_SMOTE_detection(filepath, verbose=False, debugmode=False):
	fname = filepath.split('/')[-1]
	print_debug_string(f'Isolating SMOTE detection of {fname}', debugmode=debugmode)
	lc = pd.read_csv(filepath, delimiter=' ', header=0, error_bad_lines=False)
	smote_id = lc.g_mag.idxmax()  # Assumes the brightest detection is the SMOTE
	row = lc.iloc[smote_id]
	print_verbose_string(f'SMOTE Info:\n{row}', verbose=verbose)
	print_debug_string(f'Returning ID of {smote_id}', debugmode=debugmode)
	return smote_id
		
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

    if debugmode:
        print(arguments)  

    _ = isolate_SMOTE_detection(filepath)
    print(_)
