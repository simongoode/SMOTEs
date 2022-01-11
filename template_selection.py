#!/usr/bin/python

#!/usr/bin/env python
""" template_selection.py -- A script for deciding which files should make a SMOTE template.

Usage: template_selection.py [-h] [-v] [-q] [--debug] <file> <id> [-n NUM] [-s SEP]

Arguments:
    file (string)
    	Path to a light curve file (fmt: MJD, g_mag, g_mag_err, upper_lim_mag (space-separated))
    
    id (int)
    	Index of SMOTE in the light curve file
Options:
    -h, --help                              Show this screen
    -v, --verbose                           Show extra information [default: False]   
    -q, --quietmode                         Minimize print to screen. This is useful when this function is called in another function. [default: False]  
    --debug                                 Output more for debugging [default: False]
    -n NUM, --number NUM                    Number of images to stack in template [default: 3]
    -s SEP, --separator SEP                 Number of detections to ignore around the SMOTE [default: 3]

Examples:
    bash: template_selection.py -n 5 -s 2 /fred/oz100/NOAO_archive/archive_NOAO_data/data_outputs/2015/01/4hr/g_band/single/lightcurves/files/DWF040748.225-545713.829_150116 61
"""



import docopt
import os, sys
import pandas as pd
import numpy as np

__author__	= "Simon Goode"
__license__	= "MIT"
__version__	= "1.0"
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

	
#########################################
# =========== Main Function =========== #
#########################################

def template_selection(filepath, dnum, sep, stacknum):
	fname = filepath.split('/')[-1]
	print_debug_string(f'Selecting templates for {fname}', debugmode=debugmode)
	lc = pd.read_csv(filepath, delimiter=' ', header=0, error_bad_lines=False)
	maxd = len(lc)-1
	
	dets = np.linspace(0,maxd,num=maxd+1, dtype=int)  # Array of detection indices
	ineligible = [dnum]  # Initialise the 'ineligible' list - images with indices in this list are not considered for the template
	for i in range(sep):  # Add the points nearby the SMOTE, defined by sep
		if dnum+(i+1) <= maxd:
			ineligible.append(dnum+(i+1))
		if dnum-(i+1) >= 0:
			ineligible.append(dnum-(i+1))

	eligible = [x for x in dets if x < dnum and x not in ineligible]  # Of the remaining images, try to find stacknum worth of images before the SMOTE
	temps = eligible[-stacknum:]

	if len(eligible) < stacknum:  # If there aren't enough images before the SMOTE, try looking after the SMOTE
		print_debug_string('Not enough detections before the SMOTE.', debugmode=debugmode)
		eligible = [x for x in dets if x > dnum and x not in ineligible]
		temps = eligible[:stacknum]

	if len(eligible) < stacknum:  # If there aren't enough images before or after the SMOTE, start collecting image before and after the SMOTE, in order of closeness.
		print_debug_string('Not enough detections after the SMOTE.', debugmode=debugmode)
		temps = []                  # When the number of images collected is equal to stacknum, use those. Else, it will gather all available images, assuming stacknum couldn't be reached
		for i in range(maxd):
			if len(temps) == stacknum:
				break
			if dnum-(i+1) >= 0 and dnum-(i+1) not in ineligible and dnum-(i+1) not in temps:
				temps.append(dnum-(i+1))
			if len(temps) == stacknum:
				break
			if dnum+(i+1) <= maxd and dnum+(i+1) not in ineligible and dnum+(i+1) not in temps:
				temps.append(dnum+(i+1))

	return temps  # Produce a list of indices. Images with these indices should be used to make a stacked template.

		
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
    smote_id        = arguments['<id>']
    
    sep             = arguments['--separator']
    stacknum        = arguments['--number']

    if debugmode:
        print(arguments)  

    _ = template_selection(filepath, int(smote_id), int(sep), int(stacknum))
