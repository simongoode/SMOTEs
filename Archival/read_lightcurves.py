import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord

### Function Definitions ###
def Plot_LC(f):
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
	
	plt.figure(dpi=140)
	plt.errorbar(x, y, yerr, fmt=',', color='black', capsize=2.)
	plt.scatter(x, nondet, marker='^', color='blue', facecolor='none')
	plt.ylim(23, 10)
	plt.show()

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
	
### Main Code ###
yr = 2017
mo = '02'
field = 'Antlia'		
wdir = '/fred/oz100/NOAO_archive/archive_NOAO_data/data_outputs/{}/{}/{}/g_band/single/lightcurves/files/'.format(yr, mo, field)
data = pd.DataFrame()
files,ids,mjds,mags,magerrs,mean_non,mean_nonerr = [], [], [], [], [], [], []
for filename in os.listdir(wdir):
	if filename[:3] == 'DWF':	
		try:
			df = pd.read_csv(wdir+filename, delimiter=' ', header=0, error_bad_lines=False)
			det_num = len(df.g_mag)-len(df[df.g_mag == 0])
			if det_num == 1:
				if df.g_mag.iloc[0] == 0 and df.g_mag.iloc[-1] == 0:
					a,b,c,d,e,f = Extract_Info(wdir+filename)
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
print(data)
data.to_csv('/fred/oz100/sgoode/SMOTEs/Candidate_Dataframes/{}_{}_{}.csv'.format(yr, mo, field))
