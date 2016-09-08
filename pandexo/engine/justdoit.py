import numpy as np
from pandeia.engine.instrument_factory import InstrumentFactory

ALL = {"MIRI LRS":False,
       "NIRISS SOSS":False,
       "NIRSpec G140M":False,
       "NIRSpec G140H":False,
       "NIRSpec G235M":False,
       "NIRSpec G235H":False,
       "NIRSpec G395M":False,
       "NIRSpec G395H":False,
       "NIRSpec Prism":False,
       "NIRCam F322W2":False,
       "NIRCam F444W":False}


def print_instruments():
    print "Choose from the following:"
    print ALL.keys()
    return

def load_exo_dict():
    pandexo_input = {
        "star":{
                "type" : "user or phoenix", 
                "starpath" : "file path",
                "w_unit": "Angs,cm,um,cm or Hz",
                "f_unit": "W/m2/um, FLAM, Jy, or erg/s/cm2/Hz",
                "mag": "magnitude",
                "ref_wave": "corresponding ref wave",
                "temp": "Only if phoenix, (in K)",
                "metal": "in log Fe/H",
                "logg": "cgs"
            },

    "planet" :{
	        "type": "user",
            "exopath" : "file path",
            "w_unit": "Angs,cm,um,cm or Hz",
            "f_unit": "rp/r* or fp/f*"
            },

    "observation": {
        "sat_level": "in % sat level",
        "transit_duration": "in seconds",
        "noccultations": "num transits",
        "wave_bin": "in micron",
        "fraction": "time in/out", 
        "noise_floor":"constant number or file name"
        }
    }
    print "Replace all inputs before feeding to run_modes:"
    print pandexo_input
    return pandexo_input   
    
def load_mode_dict()



def get_pce(instrument):

    obsmode = {
               'instrument': instrument,
               'mode': mode,
               'filter': config['filter'],
               'aperture': config['aperture'],
               'disperser': config['disperser']
               }
                             
    conf = {'instrument': obsmode}

    i = InstrumentFactory(config=conf)
    wr = i.get_wave_range()
    wave = np.linspace(wr['wmin'], wr['wmax'], num=500)
    pce = i.get_total_eff(wave)

    return wave,pce


def run_pandexo(modes_2_run, pandexo_input):
    for i in modes_2_run: 
        print "Starting PandExo Sim: "+ i 
        
        