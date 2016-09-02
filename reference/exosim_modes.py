class SetDefaultModes():
    """
    This class contains functionality for setting observing modes for exoplanet observations 
    This is NOT a complete set of observing modes. Instead, if offers a starting point for 
    choosing one instrument specification. There is one function for each instrument. 

    Units are set in here, although can handle: cm, micron, ang, nm, and  W/m2/um, FLAM, Jy, or erg/s/cm2/Hz
    Star: Angs, FLAM 
    Planet: cm, depth
    """

    def __init__(self, filename_star, filename_planet, jmag):
		self.filename_star =  filename_star
		self.filename_planet = filename_planet
		self.jmag = jmag        

    def niriss(self):
        #planet setup 
        star = {
            'starpath' : self.filename_star,
                'w_unit': 'Angs',
                'f_unit': 'FLAM', 
                'logg' : 'g40',
                'mag': self.jmag,
                'ref_wave': 12500.0,
                'wmin': 0.5,
                'wmax': 30.0,
                'temp': 4000.0 #nothing is using this.. 
                }

        planet = {
            'exopath' : self.filename_planet,
            'w_unit':'nm',
            'f_unit': 'rp/r*', 
            }
        observation = {
        'ncoadd': 1, #number of integrations ot lump together 
        'sat_level': 45000.0, #e, based on Valenti email
        'transit_duraton': 2.9*60.0*60.0, #from ingress to egress 
        'noccultations': 1,
        'sv_on': False, #turn on stellar variability 
        'drift_on': False,
        'wave_bin': .00,
        'fraction': 1.0
        }
        exo_input = {
        'star': star,
        'planet': planet, 
        'observation': observation

        }

        # set up default source configuration
        source = {
            'id': 1,
            'target': True,
            'position': {
                'ang_unit': 'arcsec',
                'xoff': 0.0,
                'yoff': 0.0,
            },
            'shape': {
                'major': 0.0,
                'minor': 0.0
            },
            'spectrum': {
                'normalization': {
                     'type': 'none',
                     'norm_waveunit': 'microns',
                     'norm_fluxunit': 'mjy',
                     'norm_flux': 1.0,
                     'norm_wave': 0.7},
                'sed': {
                    'sed_type': 'input',
                    'wmin': 0.3,
                    'wmax': 6.0,
                    'spectrum': [],
                            },
                'lines': []
                    }
        }
        scene=[source]

        # set up the rest of the calculation
        background = 'medium'

        effects = {
                'background': True,
                'ipc': True,
                'saturation': True
                }
        noise =   {
                'crs': False,
                'darkcurrent': False,
                'ffnoise': False,
                'readnoise': False,
                'rn_correlation': False
                }

        obsmode = {
            'instrument': 'niriss',
            'mode': 'soss',
            'filter': None,
            'aperture': 'soss',
            'disperser': 'gr700xd_1'
        }
        exp_config = {
            'subarray': 'substrip80',
            'readmode': 'nisrapid',
            'ngroup': 10,
            'nint': 1,
            'nexp': 10
        }
        strategy = {
            'method': 'specapphot',
            'aperture_key': 'disperser',
            'aperture_size': 0.1,
            'sky_annulus': [0.3, 0.4],
            #'target': [0.0, 0.0]
#'dithers': [{'x': 0.0, 'y': 0.0}]
            'units': 'arcsec',
            'target_type': '',
            'target_source': 1,
            'target_xy': [0.0, 0.0]
        }


        calc_input = {
            'scene': scene,
            'background': background,
            'calculation': {'effects': effects,
                            'noise': noise},
            'configuration': {'instrument': obsmode,
                              'detector': exp_config},
            'strategy': strategy
        }
        return {'calc_input':calc_input, 'exo_input':exo_input}

    def nirspec(self):
        #planet setup 
		star = {
            'starpath' : self.filename_star,
                'w_unit': 'Angs',
                'f_unit': 'FLAM', 
                'logg' : 'g40',
                'mag': self.jmag,
                'ref_wave': 12500.0,
                'wmin': 0.5,
                'wmax': 30.0,
                'temp': 4000.0 #nothing is using this.. 
                }

		planet = {
            'exopath' : self.filename_planet,
            'w_unit':'nm',
            'f_unit': 'rp/r*',
            }
		observation = {
        'ncoadd': 1, #number of integrations ot lump together 
        'sat_level': 45000.0, #e, based on KLAUS email
        'transit_duraton': 2.9*60.0*60.0, #from ingress to egress 
        'noccultations': 1,
        'sv_on': False, #turn on stellar variability 
        'drift_on': False,
        'wave_bin': .0,
        'fraction': 1.0
		}
		exo_input = {
        'star': star,
        'planet': planet, 
        'observation': observation

        }

        # set up default source configuration
		source = {
            'id': 1,
            'target': True,
            'position': {
                'ang_unit': 'arcsec',
                'xoff': 00.0,
                'yoff': 00.0,
            },
            'shape': {
                'major': 0.0,
                'minor': 0.0
            },
            'spectrum': {
                'normalization': {
                    'type': 'none',
                        'norm_waveunit': 'microns',
                            'norm_fluxunit': 'mjy',
                                'norm_flux': 1.0,
                                    'norm_wave': 0.7},
                'sed': {
                    'sed_type': 'input',
                    'wmin': 0.3,
                    'wmax': 6.0,
                    'spectrum': [],
                        },
                'lines': []
                }
                        }
		scene=[source]

        # set up the rest of the calculation
		background = 'medium'

		effects = {
                'background': True,
                'ipc': True,
                'saturation': True
                }
		noise =   {
                'crs': False,
                'darkcurrent': False,
                'ffnoise': False,
                'readnoise': False,
                'rn_correlation': False
                }

		obsmode = {
            'instrument': 'nirspec',
            'mode': 'fixed_slit',
            'filter': 'f100lp',
            'aperture': 's1600a1',
            'disperser': 'g140m'
        }
		exp_config = {
            'subarray': 's1600a1',
            'readmode': 'nrsrapid',
            'ngroup': 10,
            'nint': 1,
            'nexp': 10
        }
		strategy = {
            'method': 'specapphot',
            #'aperture_key': 'disperser',
            'aperture_size': 0.1,
            'sky_annulus': [0.3, 0.4],
            #'target': [0.0, 0.0]
#'dithers': [{'x': 0.0, 'y': 0.0}]
            'units': 'arcsec',
            'target_type': '',
            'target_source': 1,
            'target_xy': [0.0, 0.0]
        }

		calc_input = {
            'scene': scene,
            'background': background,
            'calculation': {'effects': effects,
                            'noise': noise},
            'configuration': {'instrument': obsmode,
                              'detector': exp_config},
            'strategy': strategy
        }
		return {'calc_input':calc_input, 'exo_input':exo_input}
    

    def nircam(self):
        #planet setup 
        star = {
            'starpath' : self.filename_star,
                'w_unit': 'Angs',
                'f_unit': 'FLAM', 
                'logg' : 'g40',
                'mag': self.jmag,
                'ref_wave': 12500.0,
                'wmin': 0.5,
                'wmax': 30.0,
                'temp': 4000.0 #nothing is using this.. 
                }

        planet = {
            'exopath' : self.filename_planet,
            'w_unit':'nm',
            'f_unit': 'rp/r*', 
            }
        observation = {
        'ncoadd': 1, #number of integrations ot lump together 
        'sat_level': 45000.0, #e, based on Valenti email 
        'transit_duraton': 2.9*60.0*60.0, #from ingress to egress 
        'noccultations': 1,
        'sv_on': False, #turn on stellar variability 
        'drift_on': False,
        'wave_bin': .0,
        'fraction': 1.0
        }
        exo_input = {
        'star': star,
        'planet': planet, 
        'observation': observation

        }

        # set up default source configuration
        source = {
            'id': 1,
            'target': True,
            'position': {
                'ang_unit': 'arcsec',
                'xoff': 0.0,
                'yoff': 0.0,
            },
            'shape': {
                'pa': -45.0,
                'major': 0.6,
                'minor': 0.2
            },
            # flat continuum spectrum normalized to 0.1 mJy
            'spectrum': {
                'w_unit': 'micron',
                'f_unit': 'millijansky',
                'type': 'input',
                'spectrum': None,
                'wmin': 0.6,
                'wmax': 2.5
            }
        }
        scene=[source]

        # set up the rest of the calculation
        background = 'medium'

        effects = {
                'background': True,
                'ipc': True,
                'saturation': True
                }
        noise =   {
                'crs': False,
                'darkcurrent': False,
                'ffnoise': False,
                'readnoise': False,
                'rn_correlation': False
                }

        obsmode = {
            'instrument': 'nircam',
            'mode': 'grism',
            #'filter': NULL ,
            'aperture': 'sw',
            'disperser': 'grism'
        }
        exp_config = {
            'subarray': '256x256',
            'readmode': 'rapid',
            'ngroup': 10,
            'nint': 1,
            'nexp': 10
        }
        strategy = {
            'method': 'imagingapphot',
            'aperture_size': 0.08,
            'sky_annulus': [0.2, 0.4],
            'target': [0.0, 0.0],
            'dithers': [{'x': 0.0, 'y': 0.0}]
        }

        calc_input = {
            'scene': scene,
            'background': background,
            'calculation': {'effects': effects,
                            'noise': noise},
            'configuration': {'instrument': obsmode,
                              'detector': exp_config},
            'strategy': strategy
        }
        return {'calc_input':calc_input, 'exo_input':exo_input}



    def miri(self):
        #planet setup 
        star = {
            'starpath' : self.filename_star,
                'w_unit': 'nm',
                'f_unit': 'FLAM', 
                'logg' : 'g40',
                'mag': self.jmag,
                'ref_wave': 1.2500,
                'wmin': 5,
                'wmax': 20.0,
                'temp': 4000.0 #nothing is using this.. 
                }

        planet = {
            'exopath' : self.filename_planet,
            'w_unit':'um', 
            'f_unit': 'fp/f*', 
            }
        observation = {
        'ncoadd': 1, #number of integrations ot lump together 
        'sat_level': 45000.0, #e, based on Valenti email 
        'transit_duraton': 2.9*60.0*60.0, #from ingress to egress 
        'noccultations': 1,
        'sv_on': False, #turn on stellar variability 
        'drift_on': False,
        'wave_bin': .0,
        'fraction': 1.0
        }
        exo_input = {
        'star': star,
        'planet': planet, 
        'observation': observation

        }

        # set up default source configuration
        source = {
            'id': 1,
            'target': True,
            'position': {
                'ang_unit': 'arcsec',
                'xoff': 0.0,
                'yoff': 0.0,
            },
            'shape': {
                'major': 0.0,
                'minor': 0.0
            },
            'spectrum': {
                'normalization': {
                     'type': 'none',
                     'norm_waveunit': 'nm',
                     'norm_fluxunit': 'mjy',
                     'norm_flux': 1.0,
                     'norm_wave': 0.7},
                'sed': {
                    'sed_type': 'input',
                    'wmin': 1250,
                    'wmax': 20000.,
                    'spectrum': [],
                            },
                'lines': []
                    }
		}
		
        scene=[source]

        # set up the rest of the calculation
        background = 'medium'
        
        effects = {
                'background': True,
                'ipc': True,
                'saturation': True
                }

        noise =   {
                'crs': False,
                'darkcurrent': False,
                'ffnoise': False,
                'readnoise': False,
                'rn_correlation': False
                }
        obsmode = {
           'instrument': 'miri',
           'mode': 'lrsslitless',
           'filter': None,
           'aperture': 'imager',
           'disperser': 'p750l'
           }
        exp_config = {
              'subarray': 'slitlessprism',
              'readmode': 'fast',
              'ngroup': 10,
              'nint': 1,
              'nexp': 10
              }
        strategy = {
            'method': 'specapphot',
            'aperture_size': 0.42,
            'sky_annulus': [0.5,1.5],
            'target': [0.0, 0.0],
            'dithers': [{'x':0.0,'y':0.0}]
            }


        calc_input = {
            'scene': scene,
            'background': background,
            'calculation': {'effects': effects,
                            'noise': noise},
            'configuration': {'instrument': obsmode,
                              'detector': exp_config},
            'strategy': strategy
        }
        return {'calc_input':calc_input, 'exo_input':exo_input}
