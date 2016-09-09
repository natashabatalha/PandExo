import os
import json
import numpy as np


class SetDefaultModes():
    """
    This class contains functionality for setting observing modes for exoplanet observations 
    This is NOT a complete set of observing modes. Instead, if offers a starting point for 
    choosing one instrument specification. There is one function for each instrument. 

    """

    def __init__(self, inst):
        self.instrument = inst[0:inst.find(' ')].lower()
        print self.instrument
        self.config = inst[inst.find(' ')+1:len(inst)].lower()   
        print self.config
    
    def pick(self): 
        getattr(self, self.instrument)()
        try: 
            return getattr(self, self.instrument)()
        except: 
            print "INVALID INSTRUMENT NAME"
                                   
    def niriss(self):
        #need to add in functionality for two orders, right now you only have option to do sub80
        order = int(self.config[7])
        
        subarray = {'2':'sub80', '1':'sub80'}
        with open(os.path.join(os.path.dirname(__file__), "reference",
                               "niriss_input.json")) as data_file:
            pandeia_data = json.load(data_file)
            pandeia_data["strategy"]["order"] = order
            pandeia_data["configuration"]["detector"]["subarray"] = subarray[str(order)]
        return pandeia_data
        
    def nirspec(self):
        filters = {'g140m':'f070lp','g140h':'f070lp',
                    'g235m':'f170lp','g235h':'f170lp',
                    'g395m':'f290lp','g395h':'f290lp',
                    'prism': 'clear'
                    }
        with open(os.path.join(os.path.dirname(__file__), "reference",
                               "nirspec_input.json")) as data_file:
            pandeia_data = json.load(data_file)
            pandeia_data["configuration"]["instrument"]["disperser"] = self.config
            pandeia_data["configuration"]["instrument"]["filter"] = filters[self.config]                                  
        return pandeia_data
        
    def nircam(self):
        with open(os.path.join(os.path.dirname(__file__), "reference",
                               "nircam_input.json")) as data_file:
            pandeia_data = json.load(data_file)
            pandeia_data["configuration"]["instrument"]["filter"]  = self.config           
        return pandeia_data

    def miri(self):
        #only lrs slit less     
        with open(os.path.join(os.path.dirname(__file__), "reference",
                               "miri_input.json")) as data_file:
                pandeia_data = json.load(data_file)
        return pandeia_data