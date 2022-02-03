#!/usr/bin/python

#!/usr/bin/env python
""" selection_criteria.py -- A script for deciding if a single light curve is a SMOTE candidate

Usage: selection_criteria.py [-h] [-v] [-q] [--debug] <file>

Arguments:
    file (string)
    	Path to a light curve file (fmt: MJD, g_mag, g_mag_err, upper_lim_mag (space-separated))
	
Options:
    -h, --help                              Show this screen
    -v, --verbose                           Show extra information [default: False]   
    -q, --quietmode                         Minimize print to screen. This is useful when this function is called in another function. [default: False]  
    --debug                                 Output more for debugging [default: False]

Examples:
    bash: python selection_criteria.py /fred/oz100/NOAO_archive/archive_NOAO_data/data_outputs/2015/01/4hr/g_band/single/lightcurves/files/DWF040748.225-545713.829_150116
"""



import docopt
import os, sys
import pandas as pd
import numpy as np
import pickle as p

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
def OutOfNowhereCriteria(mags, verbose=False, debugmode=False):
	print_debug_string('Executing OutOfNowhereCriteria', debugmode=debugmode)
	p = False  # True/False flag for passing/failing the criteria
	_ = [i for i in mags if i != 0.]
	if len(_) != 1:
		print_verbose_string('Out of nowhere criteria failed - more/less than 1 detection', verbose=verbose)
	elif mags[0] == _[0] or mags[-1] == _[0]:
		print_verbose_string('Out of nowhere criteria failed - detection is first/last in lightcurve', verbose=verbose)
	elif _[0] > 25:
		print_verbose_string('Out of nowhere criteria failed - detection is below 25th magnitude', verbose=verbose)
	else:
		print_verbose_string('SUCCESS - Out of nowhere criteria', verbose=verbose)
		p = True
	return p
	
#########################################
# =========== Main Function =========== #
#########################################

def selection_criteria(filepath, verbose=False, debugmode=False):
	lc = pd.read_csv(filepath, delimiter=' ', header=0, error_bad_lines=False)
	fname = filepath.split('/')[-1]
	print_verbose_string(f'Analysing light curve: {fname}', verbose=verbose)
	
	p = OutOfNowhereCriteria(lc['g_mag'].tolist(), verbose=verbose, debugmode=debugmode)
	if p:
		return 'Out of Nowhere Candidate'
	else:
		return False
		
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

    _ = selection_criteria(filepath, verbose=verbose, debugmode=debugmode)
    with open('out.p', 'wb') as f:
        p.dump(_, f)
