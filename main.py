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
import math
from astropy.wcs import WCS
from astropy.io import fits


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

def Initialise_Dirs(datapath, wd, c_name):
	if datapath[-1] != '/':
		datapath = datapath+'/'
		
	if not os.path.exists(datapath):
		print('{} does not exist. Please create before continuing!'.format(datapath))
		return
		
	if not os.path.exists(wd+c_name):
		os.makedirs(wd+c_name)
		os.makedirs(wd+c_name+'/full/')
		os.makedirs(wd+c_name+'/cutouts/')
		os.makedirs(wd+c_name+'/full/aligned/')
		
	elif not os.path.exists(wd+c_name+'/full/'):
		os.makedirs(wd+c_name+'/full/')
	elif not os.path.exists(wd+c_name+'/cutouts/'):
		os.makedirs(wd+c_name+'/cutouts/')
	elif not os.path.exists(wd+c_name+'/full/aligned/'):
		os.makedirs(wd+c_name+'/full/aligned/')
	return

def RAdec_to_RAsex(fRAdec):
	fratotsec = (math.fabs(float(fRAdec))*3600.0)
	frah2 = (math.modf(fratotsec/3600.0)[1])
	fram2 = (math.modf((fratotsec-(frah2*3600.0))/60.0)[1])
	fras2 = (fratotsec-(frah2*3600.0)-(fram2*60.0))
	if round(fras2, 2) == 60.00:
		fram2 = fram2 + 1
		fras2 = 0
	if round(fram2, 2) == 60.00:
		frah2 = frah2 + 1
		fram2 = 0
	if round(fram2, 2) == 60.00:
		frah2 = frah2 + 1
		fram2 = 0
	if int(frah2) == 24 and (int(fram2) != 0 or int(fras2) != 0):
		frah2 = frah2 - 24
	return '%02i' % frah2 + ' ' + '%02i' % fram2 + ' ' + ('%.3f' % float(fras2)).zfill(6)


def DEdec_to_DEsex(fDEdec):
	fdetotsec = (math.fabs(float(fDEdec))*3600.0)
	fded2 = (math.modf(fdetotsec/3600.0)[1])
	fdem2 = (math.modf((fdetotsec-(fded2*3600.0))/60.0)[1])
	fdes2 = (fdetotsec-(fded2*3600.0)-(fdem2*60.0))
	if float(fDEdec) < 0:
		fded2sign = '-'
	else:
		fded2sign = '+'
	return fded2sign + '%02i' % fded2 + ' ' + '%02i' % fdem2 + ' ' + ('%.2f' % float(fdes2)).zfill(5)


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
    
def Template_Selection(dnum, sep, stacknum, maxd):
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
    eligible = [x for x in dets if x > dnum and x not in ineligible]
    temps = eligible[:stacknum]

  if len(eligible) < stacknum:  # If there aren't enough images before or after the SMOTE, start collecting image before and after the SMOTE, in order of closeness.
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

def SMOTE_Selection(dnum, n, maxd):
  dets = np.linspace(0,maxd,num=maxd+1, dtype=int)  # Array of detection indices
  eligible = [dnum]  # Initialise the 'eligible' list - images with indices in this list are saved
  for i in range(n):  # Add the points nearby the SMOTE, defined by n
    if dnum+(i+1) <= maxd:
      eligible.append(dnum+(i+1))
    if dnum-(i+1) >= 0:
      eligible.append(dnum-(i+1))
  eligible.sort()

  return eligible  # Produce a list of indices. Images with these indices should be saved.
  
  

### Arguments ###
yr = 2015
mo = '01'
field = '4hr'
overwrite = False
verbose = False
datadir = '/fred/oz100/sgoode/SMOTEs/Data/'
wd = datadir+'{}_{}_{}/'.format(yr, mo, field)

### Main Code ###
### Create/Load Dataframe ###
if overwrite or not os.path.isfile('/fred/oz100/sgoode/SMOTEs/Candidate_Dataframes/{}_{}_{}.csv'.format(yr, mo, field)):
	if verbose:
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
	if verbose:
		print('Dataframe saved! (/fred/oz100/sgoode/SMOTEs/Candidate_Dataframes/{}_{}_{}.csv)'.format(yr, mo, field))

else:
	data = pd.read_csv('/fred/oz100/sgoode/SMOTEs/Candidate_Dataframes/{}_{}_{}.csv'.format(yr, mo, field))
	if verbose:
		print('Dataframe successfully loaded: /fred/oz100/sgoode/SMOTEs/Candidate_Dataframes/{}_{}_{}.csv'.format(yr, mo, field))



### Process Candidates ###
for i,r in data.iterrows():
	### Initialise Directories ###
	Initialise_Dirs(datadir, wd, r['Filename'])
	
	
	### Prepare Candidate Information ###
	day = r['Filename'].split('_')[1]
	ordered_im_listpath  = '/fred/oz100/NOAO_archive/archive_NOAO_data/scripts/create_lc/image_mjd_lists/FINAL_'+str(yr)+'_'+mo+'_'+field+'_'+day+'.ascii'
	ordered_ims = np.loadtxt(ordered_im_listpath, skiprows = 1, usecols=[0], dtype= str)
	RA_test = RAsex_to_RAdec(r['Filename'][3:13])
	DEC_test = DEsex_to_DEdec(r['Filename'][13:24])
	path_cutout = '/fred/oz100/sgoode/SMOTEs/Data/{}_{}_{}/{}/cutouts/'.format(yr,mo,field,r['Filename'])
	last_ccd = 1
	found = False
	max_det = len(ordered_ims)
	det_num = r['Detection Index']
	temps = Template_Selection(det_num, 3, 3, max_det)
	nearby = SMOTE_Selection(det_num, 3, max_det)
	t_file = open(wd+'{}/template_files.ascii'.format(r['Filename']), 'w+')
	n_file = open(wd+'{}/nearby_files.ascii'.format(r['Filename']), 'w+')
	if verbose:
		print(r['Filename'])
	for n,im in enumerate(ordered_ims):
		if n not in temps and n not in nearby:
			continue
		im_path = '/fred/oz100/NOAO_archive/archive_NOAO_data/data_outputs/'+str(yr)+'/'+mo+'/'+ field+'/g_band/single/'+im[:26]+'/ccds/'
		if verbose:
			print('looking in last ccd ({})'.format(last_ccd))
		for fitsim in os.listdir(im_path):
			if fitsim.endswith('ext_{}.fits'.format(last_ccd)):
				with fits.open(im_path+fitsim) as hdu:
					w = WCS(hdu[0].header)
					corners=w.calc_footprint()
					corner_1 = corners[0]
					corner_2 = corners[1]
					corner_3 = corners[2]

					if  corner_1[0] <= RA_test <=corner_2[0] and corner_1[1] >= DEC_test >= corner_3[1]:
						if verbose:
							print('found it')
						if n in temps:
							t_file.write(im_path+fitsim+'\n')
						if n in nearby:
							n_file.write(im_path+fitsim+'\n')
						found = True
						break
		    
		if found == False:
			if verbose:		    
				print('looking in other ccds')
			for fitsim in os.listdir(im_path): 
				with fits.open(im_path+fitsim) as hdu:
					w = WCS(hdu[0].header)
					corners=w.calc_footprint()
					corner_1 = corners[0]
					corner_2 = corners[1]
					corner_3 = corners[2]

					if  corner_1[0] <= RA_test <=corner_2[0] and corner_1[1] >= DEC_test >= corner_3[1]:
						last_ccd = fitsim.split('ext_')[1].split('.')[0]
						if verbose:
							print('found in ccd {}'.format(last_ccd))
						if n in temps:
							t_file.write(im_path+fitsim+'\n')
						if n in nearby:
							n_file.write(im_path+fitsim+'\n')
						break

	t_file.close()
	n_file.close()
	break  # Break for testing
