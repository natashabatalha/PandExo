import pysynphot as psyn
import numpy as np 
import pickle
import pandas as pd
from sqlalchemy import *   
import astropy.units as u 
import astropy.constants as c
import os 
from  astropy.modeling import blackbody as bb

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
    
    if input['type'] == 'user':
        star = np.genfromtxt(input['starpath'], dtype=(float, float), names='w, f') #pyfits.getdata(input['starpath'],1)
        #get flux 
        flux = star['f'] #star.field(input['logg'])
        #get wavelength and reference wavelength for mag normalization
        wave = star['w'] #star.field('WAVELENGTH')
        
        #sort if not in ascending order 
        sort = np.array([wave,flux]).T
        sort= sort[sort[:,0].argsort()]
        wave = sort[:,0]
        flux = sort[:,1] 
        sp=np.nan
        
    elif input['type'] =='phoenix':
        #make sure metal is not out of bounds
        if input['metal'] > 0.5: input['metal'] = 0.5
        sp = psyn.Icat("phoenix", input['temp'], input['metal'], input['logg'])
        sp.convert("microns")
        sp.convert("jy")
        wave = sp.wave
        flux = sp.flux
        input['w_unit'] ='um'
        input['f_unit'] = 'Jy'
        
    else: 
        raise Exception('Wrong input type for stellar spectra')
        
    ref_wave = float(input['ref_wave'])
    ref_wave = ref_wave*1e3
    
    
        #Convert evrything to nanometer for converstion based on gemini.edu  
    if input['w_unit'] == 'um':
        wave = wave*1e3
        
    elif input['w_unit'] == 'nm':
        wave = wave
  
    elif input['w_unit'] == 'cm' :
        wave = wave*1e7
 
    elif input['w_unit'] == 'Angs' :
        wave = wave*1e-1

    elif input['w_unit'] == 'Hz' :
        wave = 3e17/wave

    else: 
        raise Exception('Units are not correct. Pick um, nm, cm or Angs')

    #convert to photons/s/nm/m^2 for flux normalization based on 
    #http://www.gemini.edu/sciops/instruments/integration-time-calculators/itc-help/source-definition
    if input['f_unit'] == 'Jy':
        flux = flux*1.509e7/wave #eq. C
        
    elif input['f_unit'] == 'W/m2/um':
        flux = flux*wave/1.988e-13 #eq. D 
    elif input['f_unit'] == 'FLAM' :
        flux = flux*wave/1.988e-14 #eq. E
    elif input['f_unit'] == 'erg/s/cm2/Hz':
        flux = flux*1.509e30/wave #*4.0*np.pi
    else: 
        raise Exception('Units are not correct. Pick W/m2/um, FLAM, Jy, or erg/s/cm2/Hz')
    
    #normalize to specific mag 
    
    mag = float(input['mag'])
    norm_flux = np.interp(ref_wave, wave, flux) 

    #get zero point for J H and K 
    if (ref_wave <= 1.3e3) & (ref_wave >= 1.2e3):
        zeropoint = 1.97e7
    elif (ref_wave <= 1.7e3) & (ref_wave >= 1.6e3):
        zeropoint = 9.6e6
    elif (ref_wave <= 2.3e3) & (ref_wave >= 2.1e3):
        zeropoint = 4.5e6
    else:
        raise Exception('Only J H and K zeropoints are included')

    flux = flux/norm_flux*zeropoint*10**(-mag/2.5)
    
    #return to Pandeia units... milliJy and micron 
    flux_out_trans = flux*wave/1.509e7*1e3 #inverse of eq. C times 1e3 to get to milliJy instead of Jy 
    wave = wave*1e-3  #nm to micron

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
            flux_planet = (bb.blackbody_nu(wave_planet*u.micron, planet['temp']*u.K)*np.pi*u.sr).to(u.mJy)
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
            flux_star = (bb.blackbody_nu(wave_planet*u.micron, star['temp']*u.K)*np.pi*u.sr).to(u.mJy)
            #MAKING SURE TO ADD IN SUPID PI FOR PER STERADIAN!!!!
            flux_planet = (bb.blackbody_nu(wave_planet*u.micron, planet['temp']*u.K)*np.pi*u.sr).to(u.mJy)
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
