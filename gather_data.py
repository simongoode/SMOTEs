#!/usr/bin/python

#!/usr/bin/env python
""" gather_data.py -- Collect all data.csv files into a dataframe.
Usage: gather_data.py [-h] [-v] [-q] [--debug] <year> <month> <field> <date>

Arguments:
    year (int)
    	Year of processed data
    month (int)
    	Month of processed data
    field (string)
    	Field name of processed data
    date (string)	
    	Date of processed data. Input 'all' to combine all dates for the given year/month/field.
	
Options:
    -h, --help                              Show this screen
    -v, --verbose                           Show extra information [default: False]   
    -q, --quietmode                         Minimize print to screen. This is useful when this function is called in another function. [default: False]  
    --debug                                 Output more for debugging [default: False]

Examples:
    bash: python gather_data.py 2015 01 4hr 150115
"""



import docopt
import os, sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


__author__	= "Simon Goode"
__license__	= "MIT"
__version__	= "1.0"
__date__	= "2022-04-27"
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
# =========== Main Function =========== #
#########################################

def gather_data(year, month, field, date, verbose=False, debugmode=False):
	print_debug_string(f'Gathering data for {year}_{month}_{field}_{date}', debugmode=debugmode)
	
	dir_ext = f'/fred/oz100/sgoode/SMOTEs/Data/{year}_{month}_{field}'
	if date.lower() == 'all':
		try:
			#### Combine all dates
			out_file = f'{dir_ext}/{year}_{month}_{field}_all.df'
			
			rows_Filename =           	[]
			rows_Type =               	[]
			rows_ROBOT_Score =        	[]
			rows_RA =                 	[]
			rows_Dec =                	[]
			rows_CLASS_STAR_SMOTE =   	[]
			rows_ELLIPTICITY_SMOTE =  	[]
			rows_SPREAD_MODEL_SMOTE = 	[]
			rows_MAG_MODEL_SMOTE =    	[]
			rows_MAGERR_MODEL_SMOTE = 	[]
			rows_CLASS_STAR_SUB =     	[]
			rows_ELLIPTICITY_SUB =    	[]
			rows_SPREAD_MODEL_SUB =   	[]
			rows_MAG_MODEL_SUB =      	[]
			rows_MAGERR_MODEL_SUB =   	[]
			rows_Notes =              	[]
			
			for d in os.listdir(dir_ext):
				try:
					f = open(f'{dir_ext}/{d}/data.csv', 'r')
					dat = f.readlines()[1]
					f.close()
					
					rows_Filename.append(dat.split(',')[0].strip())
					rows_Type.append(dat.split(',')[1].strip())
					rows_ROBOT_Score.append(dat.split(',')[2].strip())
					rows_RA.append(dat.split(',')[3].strip())
					rows_Dec.append(dat.split(',')[4].strip())
					rows_CLASS_STAR_SMOTE.append(dat.split(',')[5].strip(' [] '))
					rows_ELLIPTICITY_SMOTE.append(dat.split(',')[6].strip(' [] '))
					rows_SPREAD_MODEL_SMOTE.append(dat.split(',')[7].strip(' [] '))
					rows_MAG_MODEL_SMOTE.append(dat.split(',')[8].strip(' [] '))
					rows_MAGERR_MODEL_SMOTE.append(dat.split(',')[9].strip(' [] '))
					rows_CLASS_STAR_SUB.append(dat.split(',')[10].strip(' [] '))
					rows_ELLIPTICITY_SUB.append(dat.split(',')[11].strip(' [] '))
					rows_SPREAD_MODEL_SUB.append(dat.split(',')[12].strip(' [] '))
					rows_MAG_MODEL_SUB.append(dat.split(',')[13].strip(' [] '))
					rows_MAGERR_MODEL_SUB.append(dat.split(',')[14].strip(' [] '))
					rows_Notes.append(dat.split(',')[15].strip())
					
				except Exception as e:
					print_debug_string(e, debugmode=debugmode)
					continue
			
			df = pd.DataFrame()
			df['Filename'] = rows_Filename
			df['Type'] = rows_Type
			df['ROBOT_Score'] = rows_ROBOT_Score
			df['RA'] = rows_RA
			df['Dec'] = rows_Dec
			df['CLASS_STAR_SMOTE'] = rows_CLASS_STAR_SMOTE
			df['ELLIPTICITY_SMOTE'] = rows_ELLIPTICITY_SMOTE
			df['SPREAD_MODEL_SMOTE'] = rows_SPREAD_MODEL_SMOTE
			df['MAG_MODEL_SMOTE'] = rows_MAG_MODEL_SMOTE
			df['MAGERR_MODEL_SMOTE'] = rows_MAGERR_MODEL_SMOTE
			df['CLASS_STAR_SUB'] = rows_CLASS_STAR_SUB
			df['ELLIPTICITY_SUB'] = rows_ELLIPTICITY_SUB
			df['SPREAD_MODEL_SUB'] = rows_SPREAD_MODEL_SUB
			df['MAG_MODEL_SUB'] = rows_MAG_MODEL_SUB
			df['MAGERR_MODEL_SUB'] = rows_MAGERR_MODEL_SUB
			df['Notes'] = rows_Notes
			
			print_verbose_string(df, verbose=verbose)
			if len(df) > 0:
				print_verbose_string(f'Saving dataframe: {out_file}', verbose=verbose)
				df.to_pickle(out_file)
			
		except Exception as e:
			#### What went wrong?
			print_debug_string(e, debugmode=debugmode)
			pass
	else:
		try:
			#### Combine all dates
			out_file = f'{dir_ext}/{year}_{month}_{field}_{date}.df'
			
			rows_Filename =           	[]
			rows_Type =               	[]
			rows_ROBOT_Score =        	[]
			rows_RA =                 	[]
			rows_Dec =                	[]
			rows_CLASS_STAR_SMOTE =   	[]
			rows_ELLIPTICITY_SMOTE =  	[]
			rows_SPREAD_MODEL_SMOTE = 	[]
			rows_MAG_MODEL_SMOTE =    	[]
			rows_MAGERR_MODEL_SMOTE = 	[]
			rows_CLASS_STAR_SUB =     	[]
			rows_ELLIPTICITY_SUB =    	[]
			rows_SPREAD_MODEL_SUB =   	[]
			rows_MAG_MODEL_SUB =      	[]
			rows_MAGERR_MODEL_SUB =   	[]
			rows_Notes =              	[]
			
			for d in os.listdir(dir_ext):
				if d.endswith(f'{date}'):
					try:
						f = open(f'{dir_ext}/{d}/data.csv', 'r')
						dat = f.readlines()[1]
						f.close()
						
						rows_Filename.append(dat.split(',')[0].strip())
						rows_Type.append(dat.split(',')[1].strip())
						rows_ROBOT_Score.append(dat.split(',')[2].strip())
						rows_RA.append(dat.split(',')[3].strip())
						rows_Dec.append(dat.split(',')[4].strip())
						rows_CLASS_STAR_SMOTE.append(dat.split(',')[5].strip(' [] '))
						rows_ELLIPTICITY_SMOTE.append(dat.split(',')[6].strip(' [] '))
						rows_SPREAD_MODEL_SMOTE.append(dat.split(',')[7].strip(' [] '))
						rows_MAG_MODEL_SMOTE.append(dat.split(',')[8].strip(' [] '))
						rows_MAGERR_MODEL_SMOTE.append(dat.split(',')[9].strip(' [] '))
						rows_CLASS_STAR_SUB.append(dat.split(',')[10].strip(' [] '))
						rows_ELLIPTICITY_SUB.append(dat.split(',')[11].strip(' [] '))
						rows_SPREAD_MODEL_SUB.append(dat.split(',')[12].strip(' [] '))
						rows_MAG_MODEL_SUB.append(dat.split(',')[13].strip(' [] '))
						rows_MAGERR_MODEL_SUB.append(dat.split(',')[14].strip(' [] '))
						rows_Notes.append(dat.split(',')[15].strip())
					
					except Exception as e:
						print_debug_string(e, debugmode=debugmode)
						continue
			
			df = pd.DataFrame()
			df['Filename'] = rows_Filename
			df['Type'] = rows_Type
			df['ROBOT_Score'] = rows_ROBOT_Score
			df['RA'] = rows_RA
			df['Dec'] = rows_Dec
			df['CLASS_STAR_SMOTE'] = rows_CLASS_STAR_SMOTE
			df['ELLIPTICITY_SMOTE'] = rows_ELLIPTICITY_SMOTE
			df['SPREAD_MODEL_SMOTE'] = rows_SPREAD_MODEL_SMOTE
			df['MAG_MODEL_SMOTE'] = rows_MAG_MODEL_SMOTE
			df['MAGERR_MODEL_SMOTE'] = rows_MAGERR_MODEL_SMOTE
			df['CLASS_STAR_SUB'] = rows_CLASS_STAR_SUB
			df['ELLIPTICITY_SUB'] = rows_ELLIPTICITY_SUB
			df['SPREAD_MODEL_SUB'] = rows_SPREAD_MODEL_SUB
			df['MAG_MODEL_SUB'] = rows_MAG_MODEL_SUB
			df['MAGERR_MODEL_SUB'] = rows_MAGERR_MODEL_SUB
			df['Notes'] = rows_Notes
			
			print_verbose_string(df, verbose=verbose)
			if len(df) > 0:
				print_verbose_string(f'Saving dataframe: {out_file}', verbose=verbose)
				df.to_pickle(out_file)
				
		except Exception as e:
			#### What went wrong?
			print_debug_string(e, debugmode=debugmode)
			pass

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
    date            = arguments['<date>']
            
    if debugmode:
        print(arguments)  

    _ = gather_data(year, month, field, date, verbose=verbose, debugmode=debugmode)
