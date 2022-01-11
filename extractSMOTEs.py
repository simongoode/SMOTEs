#!/usr/bin/python

#!/usr/bin/env python
""" extractSMOTEs.py -- This script will iterate through a light curve masterlist and return a pd.DataFrame (and .csv copy) of a list of eligible SMOTE candidates. 

Usage: extractSMOTEs.py [-h] [-v] [-q] [--debug] [-o SAVELOC] [-l LISTLOC] <year> <month> <field>

Arguments:
    year (int)
    	xxxx
    month (string)
    	xxxx
    field (string)
    	xxxx
	
Options:
    -h, --help                              Show this screen
    -v, --verbose                           Show extra information [default: False]   
    -q, --quietmode                         Minimize print to screen. This is useful when this function is called in another function. [default: False]  
    --debug                                 Output more for debugging [default: False]
    -o SAVELOC, --out SAVELOC               Saved output as [default: ../Candidate_Dataframes/candidates.csv]
    -l LISTLOC, --list LISTLOC              Light curve list location [default: /fred/oz100/NOAO_archive/archive_NOAO_data/data_outputs/]

Examples:
    bash: python extractSMOTEs.py --out '/dir/save/here' --list '/dir/list/here' 2015 01 4hr
"""



import docopt
import os, sys
import pandas as pd
import numpy as np

__author__	= "Simon Goode"
__license__	= "MIT"
__version__	= "1.0"
__date__	= "2022-01-06"
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
def Extract_Info(f):  #Extract information about a SMOTE candidate
	nondet = []
	c = 0
	for line in open(f, 'r').readlines():
		if line.split(' ')[0] == 'MJD':
			continue
		elif float(line.split(' ')[1]) == 0.:
			nondet.append(float(line.split(' ')[3]))
		else:
			mjd = float(line.split(' ')[0])
			mag = float(line.split(' ')[1])
			magerr = float(line.split(' ')[2])
			index = c
		c += 1
			
	return index, mjd, mag, magerr, np.mean(nondet), np.std(nondet)/np.sqrt(len(nondet))

#########################################
# =========== Main Function =========== #
#########################################

def extract_SMOTEs(year, month, field, savedir, listdir, verbose=False, debugmode=False):
	print_verbose_string(f'Creating dataframe for {year}_{month}_{field}', underscores=False, verbose=verbose)
	lightcurve_files = listdir+f'{year}/{month}/{field}/g_band/single/lightcurves/files/'
	print_debug_string(f'Searching for light curves in {listdir}{year}/{month}/{field}/g_band/single/lightcurves/files/', underscores=False, debugmode=debugmode)
	
	#Initialise dataframe and arrays
	data = pd.DataFrame()
	files,ids,mjds,mags,magerrs,mean_non,mean_nonerr = [], [], [], [], [], [], []
	
	#End script early if in debug mode
	debug_counter = 0
	
	for filename in os.listdir(lightcurve_files):
		if debugmode and debug_counter > 10:
			print_debug_string(f'Breaking loop (debug_counter = {debug_counter})', debugmode=debugmode)
			break
			
		print_debug_string(f'Opening {filename}', debugmode=debugmode)
		if filename[:3] == 'DWF':
			try:
				print_debug_string(f'  Reading into dataframe', debugmode=debugmode)
				df = pd.read_csv(lightcurve_files+filename, delimiter=' ', header=0, error_bad_lines=False)
				det_num = len(df.g_mag)-len(df[df.g_mag == 0])
				print_debug_string(f'  {det_num} detection(s) found', debugmode=debugmode)
				
				### Selection Criteria ###
				if det_num == 1:
					print_debug_string(f'    Single detection criteria met', debugmode=debugmode)
					if df.g_mag.iloc[0] == 0 and df.g_mag.iloc[-1] == 0:
						a,b,c,d,e,f = Extract_Info(lightcurve_files+filename)
						print_debug_string(f'      Extracting Info: {a}, {b}, {c}, {d}, {e}, {f}', debugmode=debugmode)
						if c > 25:  #Eliminate Source Extractor errors where magnitude is ~100
							print_debug_string(f'    Source Extractor error (g_mag = {c})', debugmode=debugmode)
							continue
						else:
							files.append(filename)
							ids.append(a)
							mjds.append(b)
							mags.append(c)
							magerrs.append(d)
							mean_non.append(e)
							mean_nonerr.append(f)
							print_debug_string(f'      SMOTE candidate recorded', debugmode=debugmode)
							debug_counter += 1
					else:
						print_debug_string('    Detection either first or last in light curve', debugmode=debugmode)
				
				print_verbose_string(f'{filename} complete', verbose=verbose)
			except FileNotFoundError:
				print_verbose_string(f'Bad File: {filename}', verbose=verbose)
		elif not quietmode:
			print_debug_string(f'File does not begin with "DWF"', debugmode=debugmode)
			print_verbose_string(f'Not a recognised file: {filename}', verbose=verbose)
		
	data['Filename'] = files
	data['Detection Index'] = ids
	data['Detection MJD'] = mjds
	data['Detection Magnitude'] = mags
	data['Detection Magnitude Error'] = magerrs
	data['Mean Nondetection Magnitude'] = mean_non
	data['Mean Nondetection Magnitude Error'] = mean_nonerr
	
	data.to_csv(f'{savedir}')
	print_verbose_string(f'Dataframe Saved! ({savedir})', verbose=verbose)
	
	return data
	
############################################################################
####################### BODY OF PROGRAM STARTS HERE ########################
############################################################################

if __name__ == "__main__":

    # Read in input arguments
    arguments       = docopt.docopt(__doc__)
    verbose         = arguments['--verbose']
    quietmode       = arguments['--quietmode']
    debugmode       = arguments['--debug']
    
    year            = arguments['<year>']
    month           = arguments['<month>']
    field           = arguments['<field>']
    
    savedir         = arguments['--out']
    listdir         = arguments['--list']

    if debugmode:
        print(arguments)  

    _ = extract_SMOTEs(year, month, field, savedir, listdir, verbose=verbose, debugmode=debugmode)
