import warnings
warnings.filterwarnings('ignore')
import pandexo.engine.justdoit as jdi # THIS IS THE HOLY GRAIL OF PANDEXO
import numpy as np
import os
exo_dict = jdi.load_exo_dict()
exo_dict['observation']['sat_level'] = 80    #saturation level in percent of full well 
exo_dict['observation']['sat_unit'] = '%' 
exo_dict['observation']['noccultations'] = 2 #number of transits 
exo_dict['observation']['R'] = None          #fixed binning. I usually suggest ZERO binning.. you can always bin later 
                                             #without having to redo the calcualtion
exo_dict['observation']['baseline'] = 1.0    #fraction of time in transit versus out = in/out
exo_dict['observation']['baseline_unit'] = 'frac' 
exo_dict['observation']['noise_floor'] = 0   #this can be a fixed level or it can be a filepath 
exo_dict['star']['type'] = 'phoenix'        #phoenix or user (if you have your own)
exo_dict['star']['mag'] = 8.0               #magnitude of the system
exo_dict['star']['ref_wave'] = 1.25         #For J mag = 1.25, H = 1.6, K =2.22.. etc (all in micron)
exo_dict['star']['temp'] = 5500             #in K 
exo_dict['star']['metal'] = 0.0             # as log Fe/H
exo_dict['star']['logg'] = 4.0  
exo_dict['planet']['type'] = 'constant'
exo_dict['planet']['depth'] = 0.01                      #other options include "um","nm" ,"Angs", "secs" (for phase curves)
exo_dict['planet']['transit_duration'] = 2.0*60.0*60.0 

print('Starting TEST run')
jdi.run_pandexo(exo_dict, ['NIRSpec G140H'], save_file=False)
print('SUCCESS') 
