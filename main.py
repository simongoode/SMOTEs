#Load ROBOT CNN
#Read light curves & apply selection criteria
	## FOR EACH LIGHT CURVE:
	#--> Find out-of-nowhere SMOTE candidates
	#--> Find quiescent-source SMOTE candidates
	
	## IF LIGHT CURVE IS A SMOTE CANDIDATE:
	#--> Append SMOTE candidates to a dataframe
	#New bits
	#--> Create SMOTE candidate directories / subdirectories
	#--> Add SMOTE ccd to dir (dir/full/smote.fits)
	#--> Create swarp-stacked templates (dir/full/template.fits)
	#--> Create aligned ccds (dir/full/aligned/.)
	#--> Create subtraction ccd (dir/full/subtraction.fits)
	#--> Create cutouts of SMOTE frame of temp, sci and sub ccds. (dir/cutouts/.)
	#--> Run image cutouts through ROBOT CNN
	#--> Append ROBOT score to dataframe

#Functions needed
	#def Initialize_dir():
		#Checks if candidate dir exists and creates one if not
		#Adds all subdirs ('full' and 'cutouts')
		#Adds the SMOTE-detection ccd to dir/full
	
	#def CreateTemplate():
		#Use swarp to create stacked template in dir/full
		
	#def AlignCCDs(im1, im2):
		#Use align_image.py to create aligned images in dir/full/aligned/.
	
	#def SubtractImages(im1, im2):
		#Use subtract_image.py to create subtraction in dir/full
		
	#def CreateCutouts([list of fits files]):
		#For each fits file in the list, create a cutout of the SMOTE in dir/cutouts
		
	#def ROBOT(temp, sci, sub):
		#Use the ROBOT CNN to produce a RB score of SMOTE
		
		
		
### Import packages ###
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd



### Function definitions ###
def Save_LC(f):
	x = []
	y = []
	yerr = []
	nondet = []
	for line in open(f, 'r').readlines():
		if line.split(' ')[0] == 'MJD':
			continue
		else:
			x.append(float(line.split(' ')[0]))
			y.append(float(line.split(' ')[1]))
			yerr.append(float(line.split(' ')[2]))
			nondet.append(float(line.split(' ')[3]))
	
	plt.figure(f, dpi=140)
	plt.errorbar(x, y, yerr, fmt=',', color='black', capsize=2.)
	plt.scatter(x, nondet, marker='^', color='blue', facecolor='none')
	plt.ylim(23, 10)
	#plt.save()

def Extract_Info(f):
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



### Arguments ###
yr = 2015
mo = '01'
field = '4hr'
overwrite = False

### Main Code ###
### Create Dataframe ###
if overwrite or not os.path.isfile('/fred/oz100/sgoode/SMOTEs/Candidate_Dataframes/{}_{}_{}.csv'.format(yr, mo, field)):
	print('Creating dataframe for {}_{}_{}'.format(yr, mo, field))
	wdir = '/fred/oz100/NOAO_archive/archive_NOAO_data/data_outputs/{}/{}/{}/g_band/single/lightcurves/files/'.format(yr, mo, field)
	data = pd.DataFrame()
	files,ids,mjds,mags,magerrs,mean_non,mean_nonerr = [], [], [], [], [], [], []
	for filename in os.listdir(wdir):
		if filename[:3] == 'DWF':	
			try:
				df = pd.read_csv(wdir+filename, delimiter=' ', header=0, error_bad_lines=False)
				det_num = len(df.g_mag)-len(df[df.g_mag == 0])
			
				### Selection Criteria ###
				if det_num == 1:
					if df.g_mag.iloc[0] == 0 and df.g_mag.iloc[-1] == 0:
						a,b,c,d,e,f = Extract_Info(wdir+filename)
						if c > 25:  #Eliminate Source Extractor errors where magnitude is ~100
							continue
						else:
							files.append(filename)
							ids.append(a)
							mjds.append(b)
							mags.append(c)
							magerrs.append(d)
							mean_non.append(e)
							mean_nonerr.append(f)
			except:
				print('Bad File: {}'.format(filename))
		else:
			print('Not a recognised file: {}'.format(filename))
		
	data['Filename'] = files
	data['Detection Index'] = ids
	data['Detection MJD'] = mjds
	data['Detection Magnitude'] = mags
	data['Detection Magnitude Error'] = magerrs
	data['Mean Nondetection Magnitude'] = mean_non
	data['Mean Nondetection Magnitude Error'] = mean_nonerr
	data.to_csv('/fred/oz100/sgoode/SMOTEs/Candidate_Dataframes/{}_{}_{}.csv'.format(yr, mo, field))

else:
	data = pd.read_csv('/fred/oz100/sgoode/SMOTEs/Candidate_Dataframes/{}_{}_{}.csv'.format(yr, mo, field))

print(data)
