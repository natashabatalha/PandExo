import pysynphot as psyn
import numpy as np 
import pickle
       
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

    #convert to photons/s/nm/m^2 for flux normalization based on gemini.edu
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
    flux = flux/norm_flux*1.97e7*10**(-mag/2.5)
    
    #return to Pandeia units... milliJy and micron 
    flux_out_trans = flux*wave/1.509e7*1e3 #inverse of eq. C times 1e3 to get to milliJy instead of Jy 

    wave = wave*1e-3  #nm to micron

    return {'flux_out_trans': flux_out_trans, 'wave': wave} 


def bothTrans(out_trans, planet) :
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
        
    elif planet['type'] == 'database':
        raise Exception("Empty Database")
    else: 
        raise Exception("Incorrect Planet File") 

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


