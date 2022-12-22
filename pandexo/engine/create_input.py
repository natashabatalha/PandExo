import numpy as np 
import pickle
import pandas as pd
from sqlalchemy import *   
import astropy.units as u 
import astropy.constants as c
import os 
from astropy.modeling.models import BlackBody
import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    import pysynphot as psyn
    
def outTrans(input) :
    """Compute out of transit spectra
    
    Computes the out of transit spectra by normalizing flux to specified 
    magnitude and convert to specified Pandeia units of milliJy and microns.  
  
    Parameters 
    ----------
    input : dict
        stellar scene which includes parameters to extract phoenix database or a 
        filename which points to a stellar spectrum 
    
    Return
    ------
    dict 
        contains wave and flux_out_trans
    """ 

    ref_wave = float(input['ref_wave'])
    mag = float(input['mag'])

    ################# USER ####################################
    if input['type'] == 'user':
        if isinstance(input['starpath'], dict):
            star = input['starpath']
        else: #if isinstance(input['starpath'], str):
            star = np.genfromtxt(input['starpath'], dtype=(float, float),
                                 names='w, f')
        #get flux 
        flux = star['f'] #star.field(input['logg'])
        #get wavelength and reference wavelength for mag normalization
        wave = star['w'] #star.field('WAVELENGTH')
        
        #sort if not in ascending order 
        sort = np.array([wave,flux]).T
        sort= sort[sort[:,0].argsort()]
        wave = sort[:,0]
        flux = sort[:,1] 
        if input['w_unit'] == 'um':
            PANDEIA_WAVEUNITS = 'um'
            
        elif input['w_unit'] == 'nm':
            PANDEIA_WAVEUNITS = 'nm'
      
        elif input['w_unit'] == 'cm' :
            PANDEIA_WAVEUNITS = 'cm'
     
        elif input['w_unit'] == 'Angs' :
            PANDEIA_WAVEUNITS = 'angstrom'

        elif input['w_unit'] == 'Hz' :
            PANDEIA_WAVEUNITS = 'Hz'

        else: 
            raise Exception('Units are not correct. Pick um, nm, cm, hz, or Angs')        

        #convert to photons/s/nm/m^2 for flux normalization based on 
        #http://www.gemini.edu/sciops/instruments/integration-time-calculators/itc-help/source-definition
        if input['f_unit'] == 'Jy':
            PANDEIA_FLUXUNITS = 'jy' 
        elif input['f_unit'] == 'FLAM' :
            PANDEIA_FLUXUNITS = 'FLAM'
        elif input['f_unit'] == 'erg/cm2/s/Hz':
            flux = flux*1e23
            PANDEIA_FLUXUNITS = 'jy' 
        else: 
            raise Exception('Units are not correct. Pick FLAM or Jy or erg/cm2/s/Hz')

        sp = psyn.ArraySpectrum(wave, flux, waveunits=PANDEIA_WAVEUNITS, fluxunits=PANDEIA_FLUXUNITS)        #Convert evrything to nanometer for converstion based on gemini.edu  
        sp.convert("nm")
        sp.convert('jy')

    ############ PHOENIX ################################################
    elif input['type'] =='phoenix':
        #make sure metal is not out of bounds
        if input['metal'] > 0.5: input['metal'] = 0.5
        sp = psyn.Icat("phoenix", input['temp'], input['metal'], input['logg'])
        sp.convert("nm")
        sp.convert("jy")
        wave = sp.wave
        flux = sp.flux
        input['w_unit'] ='nm'
        input['f_unit'] = 'jy'
        
    else: 
        raise Exception('Wrong input type for stellar spectra')
    

    ############ NORMALIZATION ################################################
    refdata = os.environ.get("PYSYN_CDBS")

    all_bps = {"H": 'bessell_h_004_syn.fits',
                 "J":'bessell_j_003_syn.fits' ,
                 "K": 'bessell_k_003_syn.fits'}

    if (ref_wave <= 1.3) & (ref_wave >= 1.2):
        filt = 'J'
    elif (ref_wave <= 1.7) & (ref_wave >= 1.6):
        filt = 'H'
    elif (ref_wave <= 2.3) & (ref_wave >= 2.1):
        filt = 'K'
    else:
        raise Exception('Only J H and K zeropoints are included')

    bp_path = os.path.join(refdata, "comp", "nonhst", all_bps[filt])
    if not os.path.exists(bp_path): 
        raise Exception("Oops! PandExo 2.0 now requires users to download this file https://archive.stsci.edu/hlsps/reference-atlases/hlsp_reference-atlases_hst_multi_everything_multi_v11_sed.tar it will untar with the structure grp/redcat/trds. Please place the directories nonhst and comp into this folder: "+refdata)

    bp = psyn.FileBandpass(bp_path)

    sp.convert('angstroms')
    bp.convert('angstroms')

    rn_sp = sp.renorm(mag, 'vegamag', bp)


    rn_sp.convert("microns")
    rn_sp.convert("mjy")

    flux_out_trans = rn_sp.flux
    wave = rn_sp.wave
    return {'flux_out_trans': flux_out_trans, 'wave': wave,'phoenix':sp} 


def bothTrans(out_trans, planet,star=None) :
    """Calculates in transit flux 
    
    Takes output from `outTrans`, which is the normalized stellar flux, and 
    creates either a transit transmission spectrum, phase curve or emission spectrum. 
    Magnitude 
  
    Parameters 
    ----------
    out_trans: dict 
        includes dictionary from `outTrans` output. 
    planet: dict
        dictionary with direction to planet spectra, wavelength and flux units
    star: dict
        (Optional) dictionary within exo_input with stellar information. Only 
        used when scaling Fortney Grid spectra to get (rp/r*)^2

    Return
    ------
    dict 
        dictionary with out of transit flux, in transit flux, original model 
        and corresponding wavelengths
    """ 
    
    if planet['type'] =='user':
        if isinstance(planet['exopath'], dict):
            load_file = planet['exopath']
        else: #if isinstance(planet['exopath'], str):
            load_file = np.genfromtxt(planet['exopath'], dtype=(float, float),
                names='w, f')
        #get wavelength 
        wave_planet = load_file['w']
        #get planet flux 
        flux_planet = load_file['f']
        
        #sort if not in ascending order 
        sort = np.array([wave_planet,flux_planet]).T
        sort= sort[sort[:,0].argsort()]
        wave_planet = sort[:,0]
        flux_planet = sort[:,1]
    
    ############## IF USER SELECTS CONSTANT VALUE ##################   
    elif planet['type'] == 'constant':
        rplan = (planet['radius']*u.Unit(planet['r_unit'])).to(u.km)
        rstar = (star['radius']*u.Unit(star['r_unit'])).to(u.km)

        #constant transit depth
        if planet['f_unit'] == 'rp^2/r*^2':
            wave_planet = np.linspace(0.5,15,1000)
            planet['depth'] = float(rplan**2 / rstar**2)
            flux_planet = np.linspace(0.5,15,1000)*0 + planet['depth']
            planet['w_unit'] = 'um'

        #constant fp/f* (using out_trans from user)
        elif planet['f_unit'] == 'fp/f*':
            planet['w_unit'] = 'um'
            wave_planet = out_trans['wave'][(out_trans['wave']>0.5) & (out_trans['wave']<15)]
            flux_star = (out_trans['phoenix'].flux*(u.Jy)).to(u.mJy)[(out_trans['wave']>0.5) & (out_trans['wave']<15)]
            #MAKING SURE TO ADD IN SUPID PI FOR PER STERADIAN!!!!
            bb = BlackBody(temperature=planet['temp']*u.K) 
            flux_planet = (bb(wave_planet*u.micron)*np.pi*u.sr).to(u.mJy)
            # ( bb planet / pheonix sed ) * (rp/r*)^2
            flux_planet = np.array((flux_planet/flux_star) * (rplan/rstar)**2.0)

    ############## IF USER SELECTS TO PULL FROM GRID ##################
    elif planet['type'] =='grid':
        try:
            db = create_engine('sqlite:///'+os.environ.get('FORTGRID_DIR'))
            header= pd.read_sql_table('header',db)
        except:
            raise Exception('Fortney Grid File Path is incorrect, or not initialized')

        #radius of star
        try:
            rstar = (star['radius']*u.Unit(star['r_unit'])).to(u.km)
        except:
            raise Exception("Radius of Star not supplied for scaling. Check exo_input['star']['radius']")

        #radius of planet
        try:
            rplan = (planet['radius']*u.Unit(planet['r_unit'])).to(u.km)
        except: 
            planet['radius'] = (1.25*c.R_jup).to(u.km)
            rplan = planet['radius']
            print('Default Planet Radius of 1.25 Rj given')

        #clouds 
        if planet['cloud'].find('flat') != -1: 
            planet['flat'] = int(planet['cloud'][4:])
            planet['ray'] = 0 
        elif planet['cloud'].find('ray') != -1: 
            planet['ray'] = int(planet['cloud'][3:])
            planet['flat'] = 0 
        elif int(planet['cloud']) == 0: 
            planet['flat'] = 0 
            planet['ray'] = 0     
        else:
            planet['flat'] = 0 
            planet['ray'] = 0 
            print('No cloud parameter not specified, default no clouds added')
        
        #chemistry 
        if planet['chem'] == 'noTiO': 
            planet['noTiO'] = True
            planet['eqchem'] = True 
        if planet['chem'] == 'eqchem': 
            planet['noTiO'] = False
            planet['eqchem'] = True 
            #grid does not allow clouds for cases with TiO
            planet['flat'] = 0 
            planet['ray'] = 0 

        #we are only using gravity of 25 and scaling by mass from there 
        fort_grav = 25.0*u.m/u.s/u.s
        df = header.loc[(header.gravity==fort_grav) & (header.temp==planet['temp'])
                           & (header.noTiO==planet['noTiO']) & (header.ray==planet['ray']) &
                           (header.flat==planet['flat'])]
        wave_planet=np.array(pd.read_sql_table(df['name'].values[0],db)['wavelength'])[::-1]

        r_lambda=np.array(pd.read_sql_table(df['name'].values[0],db)['radius'])*u.km
        z_lambda = r_lambda- (1.25*u.R_jup).to(u.km) #all fortney models have fixed 1.25 radii

        #scale with planetary mass 
        try:
            mass = (planet['mass']*u.Unit(planet['m_unit'])).to(u.kg)
            gravity = c.G*(mass)/(rplan.to(u.m))**2.0 #convert radius to m for gravity units
            #scale lambbda (this technically ignores the fact that scaleheight is altitude dependent)
            #therefore, it will not be valide for very very low gravities
            z_lambda = z_lambda*fort_grav/gravity
        except: 
            #keep original z lambda 
            gravity=25.0
            z_lambda = z_lambda*fort_grav/fort_grav
            print('Default Planet Gravity of 25 m/s2 given')  
        
        #create new wavelength dependent R based on scaled ravity
        r_lambda = z_lambda + rplan


        #finally compute (rp/r*)^2
        flux_planet = np.array(r_lambda**2/rstar**2)[::-1]
        planet['w_unit'] = 'um'
        planet['f_unit'] = 'rp^2/r*^2'
    else: 
        raise Exception("Incorrect Planet Type. Options are 'user','constant','grid'") 

    #Convert wave to micron 
    if planet['w_unit'] == 'um':
        wave_planet = wave_planet
    elif planet['w_unit'] == 'nm':
        wave_planet = wave_planet*1e-3
    elif planet['w_unit'] == 'cm':
        wave_planet = wave_planet*1e4
    elif planet['w_unit'] == 'Angs' :
        wave_planet = wave_planet*1e-4
    elif planet['w_unit'] == 'Hz' :
        wave_planet = 3e17/wave_planet
    elif planet['w_unit'] == 'sec' :
        wave_planet = wave_planet
    else: 
        raise Exception('Units are not correct. Pick um, nm, cm, Angs or sec.')

    if planet['w_unit'] == 'sec' :
        #star flux to feed into pandeia
        time = wave_planet
        flux_star = out_trans['flux_out_trans']
        wave_star = out_trans['wave']
        if planet['f_unit'] == 'fp/f*' :
            flux_planet = flux_planet 
        else: 
            print("Seconds with rp^2/r*^2 units not an option. Switch to Fp/F*")
            return 
        
        return {'time':time, 'wave':wave_star,'flux_out_trans':flux_star, 'planet_phase':flux_planet,
                'model_wave':time, 'model_spec':flux_planet, 'frac':(1.+flux_planet)}    
        
    else:
        #star flux to calc transit depth
        flux_star = out_trans['flux_out_trans']
        wave_star = out_trans['wave']

    
        #give them same wave min and wave max 
        wavemin = max([min(wave_planet), min(wave_star),0.5])
        wavemax = min([max(wave_planet),max(wave_star),15])
    
        flux_planet = flux_planet[(wave_planet>wavemin) & (wave_planet<wavemax)]
        wave_planet = wave_planet[(wave_planet>wavemin) & (wave_planet<wavemax)]
    
        flux_out_trans = np.interp(wave_planet, wave_star, flux_star)

        #convert to 1-depth 
        if planet['f_unit'] == 'rp^2/r*^2' :
            depth_fraction = 1.-flux_planet 
            flux_in_trans = depth_fraction*flux_out_trans
        elif planet['f_unit'] == 'fp/f*':
            depth_fraction = (1.0 + flux_planet)
            flux_in_trans = flux_out_trans*(1.0 + flux_planet)        
        else: 
            raise Exception('Units are not correct. Pick rp^2/r*^2 or fp/f*')
    
        results= {'wave':wave_planet, 'flux_in_trans': flux_in_trans, 'flux_out_trans':flux_out_trans,
                    'model_wave':wave_planet, 'model_spec': flux_planet, 'frac':depth_fraction} 
    return results

def hst_spec(planet,star) :
    """Calculates in transit flux 
    
    Takes output from `outTrans`, which is the normalized stellar flux, and 
    creates either a transit transmission spectrum, phase curve or emission spectrum. 
    Magnitude 
  
    Parameters 
    ----------
    planet: dict
        dictionary with direction to planet spectra, wavelength and flux units
    star: dict
        dictionary within exo_input with stellar information. Only 
        used when scaling Fortney Grid spectra to get (rp/r*)^2

    Return
    ------
    dict 
        dictionary with out of transit flux, in transit flux, original model 
        and corresponding wavelengths
    """ 
    
    if planet['type'] =='user':
        load_file = np.genfromtxt(planet['exopath'], dtype=(float, float), names='w, f')   
        #get wavelength 
        wave_planet = load_file['w']
        #get planet flux 
        flux_planet = load_file['f']
        
        #sort if not in ascending order 
        sort = np.array([wave_planet,flux_planet]).T
        sort= sort[sort[:,0].argsort()]
        wave_planet = sort[:,0]
        flux_planet = sort[:,1]
    
    ############## IF USER SELECTS CONSTANT VALUE ##################   
    elif planet['type'] == 'constant':
        rplan = (planet['radius']*u.Unit(planet['r_unit'])).to(u.km)
        rstar = (star['radius']*u.Unit(star['r_unit'])).to(u.km)

        #constant transit depth
        if planet['f_unit'] == 'rp^2/r*^2':
            wave_planet = np.linspace(0.1,3,500)
            planet['depth'] = float(rplan**2 / rstar**2)
            flux_planet = np.linspace(0.1,3,500)*0 + planet['depth']
            planet['w_unit'] = 'um'

        #constant fp/f* (using out_trans from user)
        elif planet['f_unit'] == 'fp/f*':
            planet['w_unit'] = 'um'
            wave_planet = np.linspace(0.1,3,500)

            bb = BlackBody(temperature=star['temp']*u.K) 
            flux_star = (bb(wave_planet*u.micron)*np.pi*u.sr).to(u.mJy)
            #MAKING SURE TO ADD IN SUPID PI FOR PER STERADIAN!!!!
            bb = BlackBody(temperature=planet['temp']*u.K) 
            flux_planet = (bb(wave_planet*u.micron)*np.pi*u.sr).to(u.mJy)
            # ( bb planet / pheonix sed ) * (rp/r*)^2
            flux_planet = np.array((flux_planet/flux_star) * (rplan/rstar)**2.0)


    ############## IF USER SELECTS TO PULL FROM GRID ##################
    elif planet['type'] =='grid':
        try:
            db = create_engine('sqlite:///'+os.environ.get('FORTGRID_DIR'))
            header= pd.read_sql_table('header',db)
        except:
            raise Exception('Fortney Grid File Path is incorrect, or not initialized')

        #radius of star
        try:
            rstar = (star['radius']*u.Unit(star['r_unit'])).to(u.km)
        except:
            raise Exception("Radius of Star not supplied for scaling. Check exo_input['star']['radius']")

        #radius of planet
        try:
            rplan = (planet['radius']*u.Unit(planet['r_unit'])).to(u.km)
        except: 
            planet['radius'] = (1.25*c.R_jup).to(u.km)
            rplan = planet['radius']
            print('Default Planet Radius of 1.25 Rj given')

        #clouds 
        if planet['cloud'].find('flat') != -1: 
            planet['flat'] = int(planet['cloud'][4:])
            planet['ray'] = 0 
        elif planet['cloud'].find('ray') != -1: 
            planet['ray'] = int(planet['cloud'][3:])
            planet['flat'] = 0 
        elif int(planet['cloud']) == 0: 
            planet['flat'] = 0 
            planet['ray'] = 0     
        else:
            planet['flat'] = 0 
            planet['ray'] = 0 
            print('No cloud parameter not specified, default no clouds added')
        
        #chemistry 
        if planet['chem'] == 'noTiO': 
            planet['noTiO'] = True
            planet['eqchem'] = True 
        if planet['chem'] == 'eqchem': 
            planet['noTiO'] = False
            planet['eqchem'] = True 
            #grid does not allow clouds for cases with TiO
            planet['flat'] = 0 
            planet['ray'] = 0 

        #we are only using gravity of 25 and scaling by mass from there 
        fort_grav = 25.0*u.m/u.s/u.s
        df = header.loc[(header.gravity==fort_grav) & (header.temp==planet['temp'])
                           & (header.noTiO==planet['noTiO']) & (header.ray==planet['ray']) &
                           (header.flat==planet['flat'])]
        wave_planet=np.array(pd.read_sql_table(df['name'].values[0],db)['wavelength'])[::-1]

        r_lambda=np.array(pd.read_sql_table(df['name'].values[0],db)['radius'])*u.km
        z_lambda = r_lambda- (1.25*u.R_jup).to(u.km) #all fortney models have fixed 1.25 radii

        #scale with planetary mass 
        try:
            mass = (planet['mass']*u.Unit(planet['m_unit'])).to(u.kg)
            gravity = c.G*(mass)/(rplan.to(u.m))**2.0 #convert radius to m for gravity units
            #scale lambbda (this technically ignores the fact that scaleheight is altitude dependent)
            #therefore, it will not be valide for very very low gravities
            z_lambda = z_lambda*fort_grav/gravity
        except: 
            #keep original z lambda 
            gravity=25.0
            z_lambda = z_lambda*fort_grav/fort_grav
            print('Default Planet Gravity of 25 m/s2 given')  
        
        #create new wavelength dependent R based on scaled ravity
        r_lambda = z_lambda + rplan


        #finally compute (rp/r*)^2
        flux_planet = np.array(r_lambda**2/rstar**2)[::-1]
        planet['w_unit'] = 'um'
        planet['f_unit'] = 'rp^2/r*^2'
    else: 
        raise Exception("Incorrect Planet Type. Options are 'user','constant','grid'") 

    #Convert wave to micron 
    if planet['w_unit'] == 'um':
        wave_planet = wave_planet
    elif planet['w_unit'] == 'nm':
        wave_planet = wave_planet*1e-3
    elif planet['w_unit'] == 'cm':
        wave_planet = wave_planet*1e4
    elif planet['w_unit'] == 'Angs' :
        wave_planet = wave_planet*1e-4
    elif planet['w_unit'] == 'Hz' :
        wave_planet = 3e17/wave_planet
    else: 
        raise Exception('Units are not correct. Pick um, nm, cm, Angs or sec.')

    
    return wave_planet, flux_planet
