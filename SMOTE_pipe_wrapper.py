#!/usr/bin/python

#!/usr/bin/env python
""" SMOTE_pipe_wrapper.py -- The SMOTE pipeline wrapper

Usage: SMOTE_pipe_wrapper.py [-h] [-v] [-q] [--debug] <file_list> <listpath> <impath> <saveloc>

Arguments:
    file_list (string)
    	Path to a list of light curve files (fmt: MJD, g_mag, g_mag_err, upper_lim_mag (space-separated))

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
    bash: python SMOTE_pipe_wrapper.py 2015_01_4hr_test.list /fred/oz100/NOAO_archive/archive_NOAO_data/scripts/create_lc/image_mjd_lists/FINAL_2015_01_4hr_150115.ascii /fred/oz100/NOAO_archive/archive_NOAO_data/data_outputs/2015/01/4hr/g_band/single/ /fred/oz100/sgoode/SMOTEs/Data/2015_01_4hr/
"""



import docopt
import os, sys, math
import pandas as pd
import numpy as np
from astropy.io import ascii

from selection_criteria import selection_criteria
from isolate_SMOTE_detection import isolate_SMOTE_detection
from template_selection import template_selection
from id_to_filepath import id_to_filepath
from load_ROBOT_model import load_ROBOT_model
from ROBOT_predict import ROBOT_predict
from PlotInfo import PlotInfo
from match_catalogs import find_closest_in_cat


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

def SMOTE_pipe_wrapper(file_list, listpath, impath, saveloc, model, verbose=False, debugmode=False):
	for f in open(file_list, 'r').readlines():
		filepath = f.split('\n')[0]		
		print_verbose_string(f'Pipeline: Running selection_criteria.py', verbose=verbose)
		criteria, smote_det = selection_criteria(filepath, verbose=verbose, debugmode=debugmode)
	
		if criteria != False:
		
			print_verbose_string(f'Pipeline: selection criteria passed: {criteria}', verbose=verbose)
			smote_file = id_to_filepath(filepath, smote_det, listpath, impath, verbose=verbose, debugmode=debugmode)	
			if smote_file == None:
				print_verbose_string(f"Pipeline: ...but the SMOTE file doesn't exist!", verbose=verbose)
				continue
				
			filename = filepath.split('/')[-1]
			
			#smote_det = isolate_SMOTE_detection(filepath, verbose=verbose, debugmode=debugmode)
			print_verbose_string(f'Pipeline: SMOTE ID: {smote_det}', verbose=verbose)
			print_verbose_string(f'Pipeline: SMOTE file: {smote_file}', verbose=verbose)
			
			template_dets = template_selection(filepath, smote_det, sep=3, stacknum=3, verbose=verbose, debugmode=debugmode)
			print_verbose_string(f'Pipeline: Template IDs: {template_dets}', verbose=verbose)
			
			temp_files = [id_to_filepath(filepath, x, listpath, impath, verbose=verbose, debugmode=debugmode) for x in template_dets]
			if temp_files == [None, None, None]:  # Hotfix for dithering
				if smote_det > 3:
					print_debug_string(f'temp_files: {temp_files} - trying ids [0,1,2]...', debugmode=debugmode)
					temp_files = [id_to_filepath(filepath, x, listpath, impath, verbose=verbose, debugmode=debugmode) for x in [0,1,2]]
				elif smote_det < 3:
					print_debug_string(f'temp_files: {temp_files} - trying ids [-1,-2,-3]...', debugmode=debugmode)
					temp_files = [id_to_filepath(filepath, x, listpath, impath, verbose=verbose, debugmode=debugmode) for x in [-1,-2,-3]]

			elif None in temp_files:  # If one of the template files is missing, try again with sep=6
				print_debug_string(f'temp_files: {temp_files} - again with sep=6...', debugmode=debugmode)
				template_dets = template_selection(filepath, smote_det, sep=6, stacknum=3, verbose=verbose, debugmode=debugmode)
				temp_files = [id_to_filepath(filepath, x, listpath, impath, verbose=verbose, debugmode=debugmode) for x in template_dets]
				
			if None in temp_files:  # After trying again, if still problematic, skip the file.
				print_verbose_string(f'Pipeline: ERROR: Template files not found ({temp_files})', verbose=verbose)
				continue
			
			if not os.path.exists(f'{saveloc}{filename}'):
				os.makedirs(f'{saveloc}{filename}')
				
			with open(f'{saveloc}{filename}/templates.ascii', 'w+') as f:
				for x in temp_files:
					f.write(f'{x}\n')
				f.close()
			print_verbose_string(f'Pipeline: Template files saved: {saveloc}{filename}/templates.ascii', verbose=verbose)
			
			os.system(f'swarp @{saveloc}{filename}/templates.ascii -VERBOSE_TYPE QUIET -IMAGEOUT_NAME {saveloc}{filename}/template.fits -WEIGHTOUT_NAME {saveloc}{filename}/temp_weights.fits -XML_NAME {saveloc}{filename}/xml_template.xml')
			os.system(f'swarp {smote_file} -VERBOSE_TYPE QUIET -IMAGEOUT_NAME {saveloc}{filename}/SMOTE.fits -WEIGHTOUT_NAME {saveloc}{filename}/SMOTE_weights.fits -XML_NAME {saveloc}{filename}/xml_SMOTE.xml')
			print_verbose_string(f'Pipeline: Template stacked: {saveloc}{filename}/template.fits', verbose=verbose)

			os.system(f'python /fred/oz100/sgoode/dataexplore/datavis/fits/align_image.py --swarp /apps/skylake/software/compiler/gcc/6.4.0/swarp/2.38.0/bin/swarp {saveloc}{filename}/template.fits {saveloc}{filename}/SMOTE.fits -o {saveloc}{filename}/ -q')
			print_verbose_string(f'Pipeline: Template/SMOTE aligned: {saveloc}{filename}/template.resamp.fits, SMOTE.resamp.fits', verbose=verbose)

			os.system(f'python /fred/oz100/sgoode/dataexplore/datavis/fits/subtract_image.py -s {saveloc}{filename}/subtraction.fits -o --sextractor /apps/skylake/software/mpi/gcc/6.4.0/openmpi/3.0.0/sextractor/2.19.5/bin/sex {saveloc}{filename}/template.resamp.fits {saveloc}{filename}/SMOTE.resamp.fits')
			print_verbose_string(f'Pipeline: Subtraction created: {saveloc}{filename}/subtraction.fits', verbose=verbose)
		
			temp = f'{saveloc}{filename}/template.resamp.fits'
			sci = f'{saveloc}{filename}/SMOTE.resamp.fits'
			sub = f'{saveloc}{filename}/subtraction.fits'
			ra = RAsex_to_RAdec(filename[3:13])
			dec = DEsex_to_DEdec(filename[13:24])
			try:
				pred = ROBOT_predict(model, temp, sci, sub, ra, dec, verbose=verbose, debugmode=debugmode)
				print_verbose_string(f'Pipeline: {filename} score: {pred}', verbose=verbose)
			except TypeError:
				print_verbose_string(f'Pipeline: ERROR: ROBOT_predict - buffer too small for requested array', verbose=verbose)
				continue
			try:
				os.system(f'python run_sourceextractor.py -v --debug --spreadmodel --sextractor /apps/skylake/software/mpi/gcc/6.4.0/openmpi/3.0.0/sextractor/2.19.5/bin/sex --psfex /apps/skylake/software/mpi/gcc/6.4.0/openmpi/3.0.0/psfex/3.21.1/bin/psfex --savecats {saveloc}{filename}/sourceextractor_catalogues {sci}')
			
				smote_cat = ascii.read(f'{saveloc}{filename}/sourceextractor_catalogues/SMOTE.resamp.cat')
				sub_cat = ascii.read(f'{saveloc}{filename}/sourceextractor_catalogues/subtraction.cat')
				smote_cat_df = pd.DataFrame(smote_cat.as_array())
				sub_cat_df = pd.DataFrame(sub_cat.as_array())
				smote_catline, smote_d2d = find_closest_in_cat(smote_cat_df, ra, dec)
				sub_catline, sub_d2d = find_closest_in_cat(sub_cat_df, ra, dec)
			
			
				f = open(f'{saveloc}{filename}/data.csv', 'w+')
				f.write('Filename, Type, ROBOT_Score, RA, Dec,'\
					'CLASS_STAR_SMOTE, ELLIPTICITY_SMOTE, SPREAD_MODEL_SMOTE, MAG_MODEL_SMOTE, MAGERR_MODEL_SMOTE,'\
					'CLASS_STAR_SUB, ELLIPTICITY_SUB, SPREAD_MODEL_SUB, MAG_MODEL_SUB, MAGERR_MODEL_SUB, Notes\n')
				
				f.write(f'{filename}, {criteria}, {pred}, {ra}, {dec},'\
					f'{smote_catline["CLASS_STAR"].values}, {smote_catline["ELLIPTICITY"].values}, {smote_catline["SPREAD_MODEL"].values},'\
					f'{smote_catline["MAG_MODEL"].values}, {smote_catline["MAGERR_MODEL"].values},'\
					f'{sub_catline["CLASS_STAR"].values}, {sub_catline["ELLIPTICITY"].values}, {sub_catline["SPREAD_MODEL"].values},'\
					f'{sub_catline["MAG_MODEL"].values}, {sub_catline["MAGERR_MODEL"].values}, -')
				f.close()
		
			except:
				print_verbose_string('Failed to Source Extract! Skipping this one...', verbose=verbose)
				
				f = open(f'{saveloc}{filename}/data.csv', 'w+')
				f.write('Filename, Type, ROBOT_Score, RA, Dec,'\
					'CLASS_STAR_SMOTE, ELLIPTICITY_SMOTE, SPREAD_MODEL_SMOTE, MAG_MODEL_SMOTE, MAGERR_MODEL_SMOTE,'\
					'CLASS_STAR_SUB, ELLIPTICITY_SUB, SPREAD_MODEL_SUB, MAG_MODEL_SUB, MAGERR_MODEL_SUB, Notes\n')
				
				f.write(f'{filename}, {criteria}, {pred}, {ra}, {dec},'\
					f'{np.nan}, {np.nan}, {np.nan},'\
					f'{np.nan}, {np.nan},'\
					f'{np.nan}, {np.nan}, {np.nan},'\
					f'{np.nan}, {np.nan}, Sextractor Failed')
				f.close()
				
			#print(f'{filename}, {criteria}, {pred}, {ra}, {dec},'\
				#f'{smote_catline["CLASS_STAR"].values}, {smote_catline["ELLIPTICITY"].values}, {smote_catline["SPREAD_MODEL"].values},'\
				#f'{smote_catline["MAG_MODEL"].values}, {smote_catline["MAGERR_MODEL"].values},'\
				#f'{sub_catline["CLASS_STAR"].values}, {sub_catline["ELLIPTICITY"].values}, {sub_catline["SPREAD_MODEL"].values},'\
				#f'{sub_catline["MAG_MODEL"].values}, {sub_catline["MAGERR_MODEL"].values}')

			#row = pd.DataFrame([filename, smote_det, ra, dec, pred])
			#row.to_csv(

############################################################################
####################### BODY OF PROGRAM STARTS HERE ########################
############################################################################

if __name__ == "__main__":

    # Read in input arguments
    arguments       = docopt.docopt(__doc__)
    verbose         = arguments['--verbose']
    quietmode       = arguments['--quietmode']
    debugmode       = arguments['--debug']
    
    file_list        = arguments['<file_list>']
    listpath        = arguments['<listpath>']
    impath          = arguments['<impath>']
    saveloc         = arguments['<saveloc>']
    
    if debugmode:
        print(arguments)  

    model = load_ROBOT_model('/fred/oz100/pipes/DWF_PIPE/ROBOT_pipe/CONFIG/exp_35_remade.h5', '/fred/oz100/pipes/DWF_PIPE/ROBOT_pipe/CONFIG/exp_35_remade_weights.h5', verbose=verbose, debugmode=debugmode)
    _ = SMOTE_pipe_wrapper(file_list, listpath, impath, saveloc, model, verbose=verbose, debugmode=debugmode)
    print('###############################################################')
    print('############ SMOTE pipeline successfully finished. ############')
    print('####### Exceptions below are a known bug in Tensorflow. #######')
    print('###############################################################')
