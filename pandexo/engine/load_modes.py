import os
import json
import numpy as np


class SetDefaultModes():
    """Class to load default instrument mode dicts
    
    This class contains functionality for loading observing modes for exoplanet observations 
    This is NOT a complete set of ALL possibleobserving modes. Instead, it offers a starting point for 
    choosing one instrument specification. There is one function for each instrument. For example, 
    It sets slitless mode for MIRI LRS. If users are interested in other specific observation modes 
    they should load in the dictionary and then edit individual keys. 
    
    Included modes are: 
       - "MIRI LRS"
       - "NIRISS SOSS"
       - "NIRSpec G140M"
       - "NIRSpec G140H"
       - "NIRSpec G235M"
       - "NIRSpec G235H"
       - "NIRSpec G395M"
       - "NIRSpec G395H"
       - "NIRSpec Prism"
       - "NIRCam F322W2"
       - "NIRCam F444W"
       - "WFC3 G102"
       - "WFC3 G141"
       
    Parameters
    ----------
    inst : str
        Allowable strings listed above 
    """

    def __init__(self, inst):
        self.instrument = inst[0:inst.find(' ')].lower()
        self.config = inst[inst.find(' ')+1:len(inst)].lower()   
    
    def pick(self): 
        """Points to specific instrument based on key choice
        """
        try: 
            return getattr(self, self.instrument)()
        except: 
            print("INVALID INSTRUMENT NAME")
            return 
                                   
    def wfc3(self):
        """Handles WFC3 template
        """
        #wfc3_input
        with open(os.path.join(os.path.dirname(__file__), "reference",
                               "wfc3_input.json")) as data_file:
            pandeia_data = json.load(data_file)
            pandeia_data["configuration"]["instrument"]["disperser"] = self.config
        return pandeia_data
    
    def niriss(self):
        """Handles NIRISS template
        """
        #need to add in functionality for two orders, right now you only have option to do sub80
        
        with open(os.path.join(os.path.dirname(__file__), "reference",
                               "niriss_input.json")) as data_file:
            pandeia_data = json.load(data_file)
        return pandeia_data

    def nirspec(self):
        """Handles NIRSpec template
        """
        filters = {'g140m':'f100lp','g140h':'f100lp',
                    'g235m':'f170lp','g235h':'f170lp',
                    'g395m':'f290lp','g395h':'f290lp',
                    'prism': 'clear'
                    }
        with open(os.path.join(os.path.dirname(__file__), "reference",
                               "nirspec_input.json")) as data_file:
            pandeia_data = json.load(data_file)
            pandeia_data["configuration"]["instrument"]["disperser"] = self.config
            pandeia_data["configuration"]["instrument"]["filter"] = filters[self.config]    
            if self.config == 'prism':
                pandeia_data["configuration"]["detector"]["subarray"] = 'sub512'

        return pandeia_data

    def nircam(self):
        """Handles NIRCam template
        """
        with open(os.path.join(os.path.dirname(__file__), "reference",
                               "nircam_input.json")) as data_file:
            pandeia_data = json.load(data_file)
            pandeia_data["configuration"]["instrument"]["filter"]  = self.config           
        return pandeia_data

    def miri(self):
        """Handles MIRI template
        """
        #only lrs slit less     
        with open(os.path.join(os.path.dirname(__file__), "reference",
                               "miri_input.json")) as data_file:
                pandeia_data = json.load(data_file)
        return pandeia_data
