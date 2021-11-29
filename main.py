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
