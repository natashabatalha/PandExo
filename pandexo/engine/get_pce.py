import numpy as np
from pandeia.engine.instrument_factory import InstrumentFactory

def get_pce(instrument,mode,config):

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
