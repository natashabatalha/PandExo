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
            'type': 'user',
            'starpath' : self.filename_star,
                'w_unit': 'Angs',
                'f_unit': 'FLAM', 
                'logg' : 0.0,
                'metal': 0.0, 
                'mag': self.jmag,
                'ref_wave': 1.25,#2.159,
                'wmin': 0.5,
                'wmax': 30.0,
                'temp': 4000.0 #nothing is using this.. 
                }

        planet = {
            'type': 'user',
            'exopath' : self.filename_planet,
            'w_unit':'nm',
            'f_unit': 'rp/r*', 
            }
        observation = {
        'ncoadd': 1, #number of integrations ot lump together 
        'sat_level': 100.0, #e, based on Valenti email
        'transit_duration': 2.0*60.0*60.0, #0, #from ingress to egress 
        'noccultations': 1,
        'wave_bin': .0,
        'fraction': 1.0, 
        'noise_floor': 00.0 #ppm 
        }
        exo_input = {
        'online': True,
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
                     'type': 'none',},
                'sed': {
                    'sed_type': 'input',
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
                'saturation': True,
                }
        noise =   {
                'crs': True,
                'darkcurrent': True,
                'ffnoise': False,
                'readnoise': True,
                'rn_correlation': True
                }

        obsmode = {
            'instrument': 'niriss',
            'mode': 'soss',
            'filter': None,
            'aperture': 'soss',
            'disperser': 'gr700xd'
        }
        exp_config = {
            'subarray': 'substrip80',
            'readmode': 'nisrapid',
            'ngroup': 10,
            'nint': 1,
            'nexp': 10
        }
        strategy = {
            'method': 'soss',
            'order': 1,
            'background_subtraction': False,
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
        return {'pandeia_input':calc_input, 'pandexo_input':exo_input}

    def nirspec(self):
        #planet setup 
		star = {
            'type': 'user',
            'starpath' : self.filename_star,
                'w_unit': 'Angs',
                'f_unit': 'FLAM', 
                'logg' : 0.0,
                'metal': 0.0, 
                'mag': self.jmag,
                'ref_wave': 1.25,#2.159,
                'wmin': 0.5,
                'wmax': 30.0,
                'temp': 4000.0 #nothing is using this.. 
                }

		planet = {
            'type': 'user',
            'exopath' : self.filename_planet,
            'w_unit':'nm',
            'f_unit': 'rp/r*',
            }
		observation = {
        'ncoadd': 1, #number of integrations ot lump together 
        'sat_level': 100.0, #e, based on KLAUS email
        'transit_duration': 2.0*60.0*60.0, #from ingress to egress 
        'noccultations': 1,
        'wave_bin': .0,
        'fraction': 1.0,
        'noise_floor': 00.0
		} 
        
        
		exo_input = {
        'online': True,
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
                    'type': 'none'},
                'sed': {
                    'sed_type': 'input',
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
                'crs': True,
                'darkcurrent': True,
                'ffnoise': False,
                'readnoise': True,
                'rn_correlation': True
                }

		obsmode = {
            'instrument': 'nirspec',
            'mode': 'fixed_slit',
            'filter': 'f070lp',
            'aperture': 's1600a1',
            'disperser': 'g140h'
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
            'background_subtraction':True,
            'aperture_size': 0.15,
            'sky_annulus': [0.3,0.5],
            'target_xy': [0.0, 0.0],
            'reference_wavelength':None,
            'units': 'arcsec'}



		calc_input = {
            'scene': scene,
            'background': background,
            'calculation': {'effects': effects,
                            'noise': noise},
            'configuration': {'instrument': obsmode,
                              'detector': exp_config},
            'strategy': strategy
        }
		return {'pandeia_input':calc_input, 'pandexo_input':exo_input}
    

    def nircam(self):
        #planet setup 
        star = {
            'type': 'user',
            'starpath' : self.filename_star,
                'w_unit': 'Angs',
                'f_unit': 'FLAM', 
                'logg' : 0.0,
                'metal': 0.0, 
                'mag': self.jmag,
                'ref_wave': 1.25,#2.159,
                'wmin': 0.5,
                'wmax': 30.0,
                'temp': 4000.0 #nothing is using this.. 
                }

        planet = {
            'type': 'user',
            'exopath' : self.filename_planet,
            'w_unit':'nm',
            'f_unit': 'rp/r*', 
            }
        observation = {
        'ncoadd': 1, #number of integrations ot lump together 
        'sat_level': 100.0, #% full well, based on Valenti email 
        'transit_duration': 2.0*60.0*60.0, #from ingress to egress 
        'noccultations': 1,
        'wave_bin': .0,
        'fraction': 1.0,
        'noise_floor': 0.0
        }
        exo_input = {
        'online':True,
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
            #'shape': {
            #    'major': 0.0,
            #    'minor': 0.0
            #},
            'spectrum': {
                'normalization': {
                    'type': 'none',},

            'sed': {
                    'sed_type': 'input',
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
                'crs': True,
                'darkcurrent': True,
                'ffnoise': False,
                'readnoise': True,
                'rn_correlation': True
                }

        obsmode = {
            'instrument': 'nircam',
            'mode': 'ssgrism',
            'filter': 'f322w' ,
            'aperture': 'lw',
            'disperser': 'grismr'
        }
        exp_config = {
            'subarray': 'sub64',
            'readmode': 'rapid',
            'ngroup': 10,
            'nint': 1,
            'nexp': 10
        }
        strategy = {
            'method': 'specapphot',
            'background_subtraction':True,
            'aperture_size': 0.15,
            'sky_annulus': [0.3,0.5],
            'target_xy': [0.0, 0.0],
            'reference_wavelength':None,
            'units': 'arcsec'}

        calc_input = {
            'scene': scene,
            'background': background,
            'calculation': {'effects': effects,
                            'noise': noise},
            'configuration': {'instrument': obsmode,
                              'detector': exp_config},
            'strategy': strategy
        }
        return {'pandeia_input':calc_input, 'pandexo_input':exo_input}



    def miri(self):
        #planet setup 
        star = {
            'type': 'user',
            'starpath' : self.filename_star,
                'w_unit': 'Angs',
                'f_unit': 'FLAM', 
                'logg' : 0.0,
                'metal': 0.0, 
                'mag': self.jmag,
                'ref_wave': 1.25,#2.159,
                'wmin': 5,
                'wmax': 20.0,
                'temp': 4000.0 #nothing is using this.. 
                }

        planet = {
            'type': 'user',
            'exopath' : self.filename_planet,
            'w_unit':'um', 
            'f_unit': 'rp/r*' 
            }
        observation = {
        'ncoadd': 1, #number of integrations ot lump together 
        'sat_level': 100.0, #e, based on Valenti email 
        'transit_duration': 2.0*60.0*60.0, #from ingress to egress 
        'noccultations': 1,
        'wave_bin': 0.0,
        'fraction': 1.0,
        'noise_floor': 0.0
        }
        
        exo_input = {
        'online':True,
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
            'spectrum': {
                'normalization': {
                     'type': 'none',
                     'norm_waveunit': 'nm',
                     'norm_fluxunit': 'mjy',
                     'norm_flux': 1.0,
                     'norm_wave': 0.7},
                'sed': {
                    'sed_type': 'input',
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
                'crs': True,
                'darkcurrent': True,
                'ffnoise': False,
                'readnoise': True,
                'rn_correlation': True
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
            'background_subtraction':True,
            'aperture_size': 0.42,
            'sky_annulus': [0.5,1.5],
            'target_xy': [0.0, 0.0],
            'reference_wavelength':None,
            'units': 'arcsec'}


        calc_input = {
            'scene': scene,
            'background': background,
            'calculation': {'effects': effects,
                            'noise': noise},
            'configuration': {'instrument': obsmode,
                              'detector': exp_config},
            'strategy': strategy
        }
        return {'pandeia_input':calc_input, 'pandexo_input':exo_input}
