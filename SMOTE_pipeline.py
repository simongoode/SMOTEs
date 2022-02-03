#!/usr/bin/python

#!/usr/bin/env python
""" SMOTE_pipeline.py -- An example of the wrapper script

Usage: SMOTE_pipeline.py [-h] [-v] [-q] [--debug] <file> <listpath> <impath> <saveloc>

Arguments:
    file (string)
    	Path to a light curve file (fmt: MJD, g_mag, g_mag_err, upper_lim_mag (space-separated))

    listpath (string)
    	Path to an ordered image list

    impath (string)
	Path to ccds for the appropriate field / date of observations
	
    saveloc (string)
    	Path to a directory where data products will be saved
		
Options:
    -h, --help                              Show this screen
    -v, --verbose                           Show extra information [default: False]   
    -q, --quietmode                         Minimize print to screen. This is useful when this function is called in another function. [default: False]  
    --debug                                 Output more for debugging [default: False]

Examples:
    bash: python SMOTE_pipeline.py /fred/oz100/NOAO_archive/archive_NOAO_data/data_outputs/2015/01/4hr/g_band/single/lightcurves/files/DWF040748.225-545713.829_150116 /fred/oz100/NOAO_archive/archive_NOAO_data/scripts/create_lc/image_mjd_lists/FINAL_2015_01_4hr_150116.ascii /fred/oz100/NOAO_archive/archive_NOAO_data/data_outputs/2015/01/4hr/g_band/single/
"""



import docopt
import os, sys, math
import pandas as pd
import numpy as np
import pickle as p

from selection_criteria import selection_criteria
from isolate_SMOTE_detection import isolate_SMOTE_detection
from template_selection import template_selection
from id_to_filepath import id_to_filepath
from load_ROBOT_model import load_ROBOT_model
from ROBOT_predict import ROBOT_predict

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
	
#########################################
# =========== Main Function =========== #
#########################################

def SMOTE_pipeline(filepath, listpath, impath, saveloc, verbose=False, debugmode=False):
	print_verbose_string(f'Pipeline: Running selection_criteria.py', verbose=verbose)
	criteria = selection_criteria(filepath, verbose=verbose, debugmode=debugmode)
	
	if criteria != False:
		filename = filepath.split('/')[-1]
		print_verbose_string(f'Pipeline: selection criteria passed: {criteria}', verbose=verbose)
		
		smote_det = isolate_SMOTE_detection(filepath, verbose=verbose, debugmode=debugmode)
		print_verbose_string(f'Pipeline: SMOTE ID: {smote_det}', verbose=verbose)
		
		template_dets = template_selection(filepath, smote_det, sep=3, stacknum=3, verbose=verbose, debugmode=debugmode)
		print_verbose_string(f'Pipeline: Template IDs: {template_dets}', verbose=verbose)
		
		temp_files = [id_to_filepath(filepath, x, listpath, impath, verbose=verbose, debugmode=debugmode) for x in template_dets]
		with open(f'{saveloc}{filename}/templates.ascii', 'w+') as f:
			for x in temp_files:
				f.write(f'{x}\n')
			f.close()
		print_verbose_string(f'Pipeline: Template files saved: {saveloc}{filename}/templates.ascii', verbose=verbose)


		smote_file = id_to_filepath(filepath, smote_det, listpath, impath, verbose=verbose, debugmode=debugmode)	
		print_verbose_string(f'Pipeline: SMOTE file: {smote_file}', verbose=verbose)
		
		os.system(f'swarp @{saveloc}{filename}/templates.ascii -VERBOSE_TYPE QUIET -IMAGEOUT_NAME {saveloc}{filename}/template.fits -WEIGHTOUT_NAME {saveloc}{filename}/temp_weights.fits -XML_NAME {saveloc}{filename}/xml_template.xml')
		os.system(f'swarp {smote_file} -VERBOSE_TYPE QUIET -IMAGEOUT_NAME {saveloc}{filename}/SMOTE.fits -WEIGHTOUT_NAME {saveloc}{filename}/SMOTE_weights.fits -XML_NAME {saveloc}{filename}/xml_SMOTE.xml')
		print_verbose_string(f'Pipeline: Template stacked: {saveloc}{filename}/template.fits', verbose=verbose)

		os.system(f'python /fred/oz100/sgoode/dataexplore/datavis/fits/align_image.py --swarp /apps/skylake/software/compiler/gcc/6.4.0/swarp/2.38.0/bin/swarp {saveloc}{filename}/template.fits {saveloc}{filename}/SMOTE.fits -o {saveloc}{filename}/ -q')
		print_verbose_string(f'Pipeline: Template/SMOTE aligned: {saveloc}{filename}/template.resamp.fits, SMOTE.resamp.fits', verbose=verbose)

		os.system(f'python /fred/oz100/sgoode/dataexplore/datavis/fits/subtract_image.py -s {saveloc}{filename}/subtraction.fits -o --sextractor /apps/skylake/software/mpi/gcc/6.4.0/openmpi/3.0.0/sextractor/2.19.5/bin/sex {saveloc}{filename}/template.resamp.fits {saveloc}{filename}/SMOTE.resamp.fits')
		print_verbose_string(f'Pipeline: Subtraction created: {saveloc}{filename}/subtraction.fits', verbose=verbose)

		model = load_ROBOT_model('/fred/oz100/pipes/DWF_PIPE/ROBOT_pipe/CONFIG/exp_35_remade.h5', '/fred/oz100/pipes/DWF_PIPE/ROBOT_pipe/CONFIG/exp_35_remade_weights.h5', verbose=verbose, debugmode=debugmode)
		print_verbose_string(f'Pipeline: ROBOT Loaded', verbose=verbose)
		
		temp = f'{saveloc}{filename}/template.resamp.fits'
		sci = f'{saveloc}{filename}/SMOTE.resamp.fits'
		sub = f'{saveloc}{filename}/subtraction.fits'
		ra = RAsex_to_RAdec(filename[3:13])
		dec = DEsex_to_DEdec(filename[13:24])
		pred = ROBOT_predict(model, temp, sci, sub, ra, dec, verbose=verbose, debugmode=debugmode)
		print_verbose_string(f'Pipeline: {filename} score: {pred}', verbose=verbose)
		
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
    listpath        = arguments['<listpath>']
    impath          = arguments['<impath>']
    saveloc         = arguments['<saveloc>']
    
    if debugmode:
        print(arguments)  

    _ = SMOTE_pipeline(filepath, listpath, impath, saveloc, verbose=verbose, debugmode=debugmode)
    
