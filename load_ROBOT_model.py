#!/usr/bin/python

#!/usr/bin/env python
""" load_ROBOT_model.py -- A script for loading the ROBOT ML model (script does NOT provide predictions).

Usage: load_ROBOT_model.py [-h] [-v] [-q] [--debug] [-m MODEL] [-w WEIGHTS]

Arguments:

	
Options:
    -h, --help                              Show this screen
    -v, --verbose                           Show extra information [default: False]   
    -q, --quietmode                         Minimize print to screen. This is useful when this function is called in another function. [default: False]  
    --debug                                 Output more for debugging [default: False]
    -m MODEL, --model MODEL                 Keras model to load. [default: /fred/oz100/pipes/DWF_PIPE/ROBOT_pipe/CONFIG/exp_35_remade.h5]
    -w WEIGHTS, --weights WEIGHTS           Keras model weights to load. [default: /fred/oz100/pipes/DWF_PIPE/ROBOT_pipe/CONFIG/exp_35_remade_weights.h5]
    
Examples:
    bash: load_ROBOT_model.py -v -m trained_model.h5 -w model_weights.h5
"""



import docopt
import numpy as np
import sys, os
from tensorflow import keras
from keras.models import Sequential, Model
from keras.layers import Activation, Dense, Dropout, Flatten, Input, Concatenate, Conv2D, MaxPooling2D, GlobalAveragePooling2D, AveragePooling2D
import keras.backend as K

__author__	= "Simon Goode"
__license__	= "MIT"
__version__	= "0.1"
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
def MCC(y_true, y_pred):
    y_pred_pos = K.round(K.clip(y_pred, 0, 1))
    y_pred_neg = 1 - y_pred_pos

    y_pos = K.round(K.clip(y_true, 0, 1))
    y_neg = 1 - y_pos

    tp = K.sum(y_pos * y_pred_pos)
    tn = K.sum(y_neg * y_pred_neg)

    fp = K.sum(y_neg * y_pred_pos)
    fn = K.sum(y_pos * y_pred_neg)

    numerator = (tp * tn - fp * fn)
    denominator = K.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))

    return numerator / (denominator + K.epsilon())
    
    
def f1_score(y_true, y_pred): #taken from old keras source code
    true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
    possible_positives = K.sum(K.round(K.clip(y_true, 0, 1)))
    predicted_positives = K.sum(K.round(K.clip(y_pred, 0, 1)))
    precision = true_positives / (predicted_positives + K.epsilon())
    recall = true_positives / (possible_positives + K.epsilon())
    f1_val = 2*(precision*recall)/(precision+recall+K.epsilon())
    return f1_val

    
#########################################
# =========== Main Function =========== #
#########################################

def load_ROBOT_model(model, weights, verbose=False, debugmode=False):
	print_verbose_string(f'load_ROBOT_model: Loading model ({model})', verbose=verbose)
	loaded_model = keras.models.load_model(model, custom_objects={'f1_score':f1_score, 'MCC':MCC})
	print_verbose_string(f'load_ROBOT_model: Loading weights ({weights})', verbose=verbose)
	loaded_model.load_weights(weights)
	print_verbose_string(f'load_ROBOT_model: Model loaded successfully', verbose=verbose)
	if verbose:
		loaded_model.summary()
	
	return loaded_model

		

############################################################################
####################### BODY OF PROGRAM STARTS HERE ########################
############################################################################

if __name__ == "__main__":

    # Read in input arguments
    arguments       = docopt.docopt(__doc__)
    verbose         = arguments['--verbose']
    quietmode       = arguments['--quietmode']
    debugmode       = arguments['--debug']    
    model           = arguments['--model']
    weights         = arguments['--weights']

    if debugmode:
        print(arguments)  

    _ = load_ROBOT_model(str(model), str(weights), verbose=verbose, debugmode=debugmode)
    '''
    test_arr = np.zeros((1,31,31,3))
    mdl = load_ROBOT_model(model, weights, verbose=verbose, debugmode=debugmode)
    pred = mdl.predict(test_arr)[0][0]
    #print(pred)
    '''
