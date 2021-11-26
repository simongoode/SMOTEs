import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.io import fits
import math
import os
from astropy.nddata.utils import Cutout2D
import pandas as pd

###################-------- Please Entre your Parametres FOR CANVIS ---------#######################
fname = '2015_01_4hr.csv'
df = pd.read_csv('/fred/oz100/sgoode/SMOTEs/Candidate_Dataframes/{}'.format(fname))
year = str(fname.split('_')[0])
month = str(fname.split('_')[1])
field = str(fname.split('_')[2].split('.')[0])
dates = []
for i,r in df.iterrows():
	dates.append(r['Filename'].split('_')[1])
dates = np.unique(dates)
object_ids = df.Filename.tolist()
size = 121
files_all = []

##############################------------------Funcations ---------------###############################

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


###################-------- Candidate Visualisation  ---------#######################   

counter = 0
for day in dates:
    ordered_im_listpath  = '/fred/oz100/NOAO_archive/archive_NOAO_data/scripts/create_lc/image_mjd_lists/FINAL_'+year+'_'+month+'_'+field+'_'+day+'.ascii'
    ordered_ims = np.loadtxt(ordered_im_listpath, skiprows = 1, usecols=[0], dtype= str)
    for i in object_ids:
	    if i[25:] != day:
	    	continue

	    RA_test = RAsex_to_RAdec(i[3:13])
	    DEC_test = DEsex_to_DEdec(i[13:24])
	    path_cutout = '/fred/oz100/sgoode/SMOTEs/Images/SMOTE_Candidates/{}_{}_{}/{}/'.format(year,month,field,i)
	    single_cutout = '/fred/oz100/sgoode/SMOTEs/Images/SMOTE_Singles/{}_{}_{}/{}/'.format(year,month,field,i)
	    if not os.path.exists(path_cutout):
	    	os.makedirs(path_cutout)
	    if not os.path.exists(single_cutout):
	    	os.makedirs(single_cutout)
	    last_ccd = 1
	    found = False
	    max_det = len(ordered_ims)
	    det_num = df[df['Filename'] == i]['Detection Index'].tolist()[0]
	    dets = np.linspace(det_num-3, det_num+3, num=7, dtype=int)
	    dets = [x for x in dets if x not in [-1,-2,max_det,max_det+1]]
	    files = []
	    
	    for n,im in enumerate(ordered_ims[dets]):
		    im_path = '/fred/oz100/NOAO_archive/archive_NOAO_data/data_outputs/'+year+'/'+month+'/'+ field+'/g_band/single/'+im[:26]+'/ccds/'
		    #print('looking in last ccd ({})'.format(last_ccd))
		    for fitsim in os.listdir(im_path):
			    if fitsim.endswith('ext_{}.fits'.format(last_ccd)):
				    with fits.open(im_path+fitsim) as hdu:
					    w = WCS(hdu[0].header)
					    head = hdu[0].header
					    xlim=head['NAXIS1']
					    ylim=head['NAXIS2']

					    pixcrd_im = np.array([[xlim, ylim]], np.float_)
					    world_im = w.wcs_pix2world(pixcrd_im, 1)
					    pixx_im, pixy_im = world_im[0][0], world_im[0][1]				

					    corners=w.calc_footprint()
					    corner_1 = corners[0]
					    corner_2 = corners[1]
					    corner_3 = corners[2]
					    corner_4 = corners[3] 
					    differnce = corner_1 - corner_2 

					    pixcrd = np.array([[RA_test, DEC_test]], np.float_)
					    worldpix = w.wcs_world2pix(pixcrd, 1)
					    pixx, pixy = worldpix[0][0], worldpix[0][1]

					    if  corner_1[0] <= RA_test <=corner_2[0] and corner_1[1] >= DEC_test >= corner_3[1]:
						    files.append(im_path+fitsim)
						    found = True
						    break
		    
		    if found == False:		    
			    #print('looking in other ccds')
			    for fitsim in os.listdir(im_path): 
				    with fits.open(im_path+fitsim) as hdu:
					    w = WCS(hdu[0].header)
					    head = hdu[0].header
					    xlim=head['NAXIS1']
					    ylim=head['NAXIS2']		
	
					    pixcrd_im = np.array([[xlim, ylim]], np.float_)
					    world_im = w.wcs_pix2world(pixcrd_im, 1)
					    pixx_im, pixy_im = world_im[0][0], world_im[0][1]				

					    corners=w.calc_footprint()
					    corner_1 = corners[0]
					    corner_2 = corners[1]
					    corner_3 = corners[2]
					    corner_4 = corners[3] 
					    differnce = corner_1 - corner_2 

					    pixcrd = np.array([[RA_test, DEC_test]], np.float_)
					    worldpix = w.wcs_world2pix(pixcrd, 1)
					    pixx, pixy = worldpix[0][0], worldpix[0][1]	

					    if  corner_1[0] <= RA_test <=corner_2[0] and corner_1[1] >= DEC_test >= corner_3[1]:
						    files.append(im_path+fitsim)
						    break
	    counter += 1
	    files_all.append(files)
	    print('Finished processing {}/{} candidates'.format(counter, len(object_ids)))
	    #print('cutouts for ' + str(i) + ' made, see them here: ' + str(path_cutout))
print(files_all)
df['Files'] = files_all
df.to_csv('/fred/oz100/sgoode/SMOTEs/Candidate_Dataframes/{}'.format(fname))
