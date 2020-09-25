import numpy as np
import pandas as pd
from .create_input import hst_spec
import batman
from .RECTE import RECTE
from scipy import optimize

def wfc3_GuessNOrbits(trdur):
    '''Predict number of HST orbits

    Predict number of HST orbits for transit observation when not provided by the user.

    Parameters
    ----------
    trdur : float
        transit duration in days

    Returns
    -------
    float
        number of requested orbits per transit (including discarded thermal-settling orbit)
    '''
    # Compute # of HST orbits during transit
    # ~96 minutes per HST orbit
    orbitsTr = trdur*24.*60/96
    if orbitsTr <= 1.5:
        norbits = 4.
    elif orbitsTr <= 2.0:
        norbits = 5.
    else:
        norbits = np.ceil(orbitsTr*2+1)

    return norbits


def wfc3_GuessParams(jmag, disperser, scanDirection, subarray, obsTime, maxScanHeight=180., maxExptime=150., targetFluence=30000., hmag=None):
    '''Predict nsamp and samp_seq when values not provided by the user.

    Parameters
    ----------
    jmag : float
        J-band magnitude
    disperser : str
        grism ('G141' or 'G102')
    scanDirection : str
        spatial scan direction ('Forward' or 'Round Trip')
    subarray : str
        Subarray aperture ('grism256' or 'grism512')
    obsTime : float
        Available observing time per HST orbit in seconds
    maxScanHeight : float
        (optional) maximum scan height in pixels
    maxExptime : float
        (Optional) default=150.0, maximum exposure time in seconds
    targetFluence : float
        (Optional) Desired fluence in electrons per pixel
    hmag : float
        (Optional) H-band magnitude

    Returns
    -------
    float
        nsamp--number of up-the-ramp samples (1..15)
    str
        samp_seq--time between non-destructive reads
    '''
    allnsamp = np.arange(1, 16)
    allsampseq = ['spars5', 'spars10', 'spars25']
    maxDutyCycle = 0
    for samp_seq in allsampseq:
        for nsamp in allnsamp:
            exptime, tottime, scanRate, scanHeight, fluence = wfc3_obs(jmag, disperser, scanDirection,
                                                                       subarray, nsamp, samp_seq, targetFluence, hmag)
            # Compute duty cycle and compare
            # Exposure time should be less than 2.5 minutes to achieve good time resolution
            ptsOrbit = np.floor(obsTime/tottime)
            dutyCycle = (exptime*(ptsOrbit-1))/50./60*100
            if (dutyCycle > maxDutyCycle) and (exptime < maxExptime) and (scanHeight < maxScanHeight):
                maxDutyCycle = dutyCycle
                bestsampseq = samp_seq
                bestnsamp = nsamp

    return bestnsamp, bestsampseq


def wfc3_obs(jmag, disperser, scanDirection, subarray, nsamp, samp_seq, targetFluence=30000., hmag=None):
    '''Determine the recommended exposure time, scan rate, scan height, and overheads.

    Parameters
    ----------
    jmag : float
        J-band magnitude
    disperser : str
        Grism ('G141' or 'G102')
    scanDirection : str
        spatial scan direction ('Forward' or 'Round Trip')
    subarray : str
        Subarray aperture ('grism256' or 'grism512')
    nsamp : float
        Number of up-the-ramp samples (1..15)
    samp_seq : str
        Time between non-destructive reads ('SPARS5', 'SPARS10', or 'SPARS25')
    targetFluence : float
        (Optional) Desired fluence in electrons per pixel
    hmag : float
        (Optional) H-band magnitude

    Returns
    -------
    float
        exptime--exposure time in seconds
    float
        tottime--total frame time including overheads in seconds
    float
        scanRate--recommented scan rate in arcsec/s
    float
        scanHeight--scan height in pixels
    float
        fluence--maximum pixel fluence in electrons
    '''
    # Estimate exposure time
    if subarray == 'grism512':
        # GRISM512
        if samp_seq == 'spars5':
            exptime = 0.853 + (nsamp-1)*2.9215  # SPARS5
        elif samp_seq == 'spars10':
            exptime = 0.853 + (nsamp-1)*7.9217  # SPARS10
        elif samp_seq == 'spars25':
            exptime = 0.853 + (nsamp-1)*22.9213  # SPARS25
        else:
            print(("****HALTED: Unknown SAMP_SEQ: %s" % samp_seq))
            return
    else:
        # GRISM256
        if samp_seq == 'spars5':
            exptime = 0.280 + (nsamp-1)*2.349  # SPARS5
        elif samp_seq == 'spars10':
            exptime = 0.278 + (nsamp-1)*7.3465  # SPARS10
        elif samp_seq == 'spars25':
            exptime = 0.278 + (nsamp-1)*22.346  # SPARS25
        else:
            print(("****HALTED: Unknown SAMP_SEQ: %s" % samp_seq))
            return

    # Recommended scan rate
    if hmag == None:
        hmag = jmag
    #scanRate = np.round(1.9*10**(-0.4*(hmag-5.9)), 3)  # arcsec/s
    #scanRate = np.round(2363./targetFluence*10**(-0.4*(jmag-9.75)), 3) # arcsec/s
    scanRate = (2491./targetFluence)*10**(-0.4*(jmag-9.75)) - (161/targetFluence)*10**(-0.4*(jmag-hmag))
    if disperser == 'g102':
        # G102/G141 flux ratio is ~0.8
        scanRate *= 0.8
    # Max fluence in electrons/pixel
    #fluence = (5.5/scanRate)*10**(-0.4*(hmag-15))*2.4  # electrons
    #fluence = (2363./scanRate)*10**(-0.4*(jmag-9.75))  # electrons
    fluence = (2491./scanRate)*10**(-0.4*(jmag-9.75)) - (161/scanRate)*10**(-0.4*(jmag-hmag))
    if disperser == 'g102':
        # WFC3_ISR_2012-08 states that the G102/G141 scale factor is 0.96 DN/electron
        fluence *= 0.96
        # G102/G141 flux ratio is ~0.8
        fluence *= 0.8
    # Scan height in pixels
    scanRatep = scanRate/0.121  # pixels/s
    scanHeight = scanRatep*exptime  # pixels
    '''
    #Quadratic correlation between scanRate and read overhead
    foo     = np.array([0.0,0.1,0.3,0.5,1.0,2.0,3.0,4.0])/0.121
    bar     = np.array([ 40, 40, 41, 42, 45, 50, 56, 61])
    c       = np.polyfit(foo,bar,2)
    model   = c[2] + c[1]*foo + c[0]*foo**2
    #c = [  6.12243227e-04,   6.31621064e-01,   3.96040946e+01]
    '''
    # Define instrument overheads (in seconds)
    c = [6.12243227e-04, 6.31621064e-01, 3.96040946e+01]
    read = c[2] + c[1]*scanRatep + c[0]*scanRatep**2
    # Correlation between scanHeight/scanRate and pointing overhead was determined elsewhere
    if scanDirection == 'Round Trip':
        # Round Trip scan direction doesn't have to return to starting point, therefore no overhead
        pointing = 0.
    elif scanDirection == 'Forward':
        c = [3.18485340e+01,   3.32968829e-02,   1.65687590e-02,
             7.65510038e-01,  -6.24504499e+01,   5.51452028e-03]
        pointing = c[0]*(1 - np.exp(-c[2]*(scanHeight-c[4]))) + \
            c[1]*scanHeight + c[3]*scanRatep + c[5]*scanRatep**2
    else:
        print(("****HALTED: Unknown scan direction: %s" % scanDirection))
        return
    # Estimate total frame time including overheads
    tottime = exptime+read+pointing  # seconds

    return exptime, tottime, scanRate, scanHeight, fluence


def wfc3_TExoNS(dictinput):
    '''Compute Transit depth uncertainty

    Compute the transit depth uncertainty for a defined system and number of spectrophotometric channels.
    Written by Kevin Stevenson      October 2016

    Parameters
    ----------
    dictinput : dict
        dictionary containing instrument parameters and exoplanet specific parameters. {"pandeia_input":dict1, "pandexo_input":dict1}

    Returns
    -------
    float
        deptherr--transit depth uncertainty per spectrophotometric channel
    float
        chanrms--light curve root mean squarerms
    float
        ptsOrbit--number of HST frames per orbit
    '''
    pandeia_input = dictinput['pandeia_input']
    pandexo_input = dictinput['pandexo_input']

    jmag = pandexo_input['star']['jmag']
    if np.ma.is_masked(jmag):
        print("Jmag not found.")
    try:
        hmag = pandexo_input['star']['hmag']
    except:
        hmag = jmag
        print("Hmag not found. Assuming no color dependence in the stellar type.")
    trdur = pandexo_input['planet']['transit_duration']
    numTr = pandexo_input['observation']['noccultations']

    schedulability = pandeia_input['strategy']['schedulability']
    scanDirection = pandeia_input['strategy']['scanDirection']
    nchan = pandeia_input['strategy']['nchan']
    norbits = pandeia_input['strategy']['norbits']
    useFirstOrbit = pandeia_input['strategy']['useFirstOrbit']
    try:
        targetFluence = pandeia_input['strategy']['targetFluence']
    except:
        targetFluence = 30000.
        print("Assuming a target fluence of 30,000 electrons.")
    disperser = pandeia_input['configuration']['instrument']['disperser'].lower(
    )
    subarray = pandeia_input['configuration']['detector']['subarray'].lower()
    nsamp = pandeia_input['configuration']['detector']['nsamp']
    samp_seq = pandeia_input['configuration']['detector']['samp_seq']

    try:
        samp_seq = samp_seq.lower()
    except:
        pass

    if disperser == 'g141':
        # Define reference Jmag, flux, variance, and exposure time for GJ1214
        refmag = 9.750
        refflux = 2.32e8
        refvar = 2.99e8
        refexptime = 88.436
    elif disperser == 'g102':
        # Define reference Jmag, flux, variance, and exposure time for WASP12
        refmag = 10.477
        refflux = 8.26e7
        refvar = 9.75e7
        refexptime = 103.129
    else:
        print(("****HALTED: Unknown disperser: %s" % disperser))
        return

    # Determine max recommended scan height
    if subarray == 'grism512':
        maxScanHeight = 430
    elif subarray == 'grism256':
        maxScanHeight = 180
    else:
        print(("****HALTED: Unknown subarray aperture: %s" % subarray))
        return

    # Define maximum frame time
    maxExptime = 150.

    # Define available observing time per HST orbit in seconds
    if str(schedulability) == '30':
        obsTime = 51.3*60
    elif str(schedulability) == '100':
        obsTime = 46.3*60
    else:
        print(("****HALTED: Unknown schedulability: %s" % schedulability))
        return

    # Compute recommended number of HST orbits and compare to user specified value
    guessorbits = wfc3_GuessNOrbits(trdur)
    if norbits == None:
        norbits = guessorbits
    elif norbits != guessorbits:
        print(("****WARNING: Number of specified HST orbits does not match number of recommended orbits: %0.0f" % guessorbits))

    if nsamp == 0 or nsamp == None or samp_seq == None or samp_seq == "none":
        # Estimate reasonable values
        nsamp, samp_seq = wfc3_GuessParams(jmag, disperser, scanDirection, subarray, obsTime, maxScanHeight, maxExptime, targetFluence, hmag)

    # Calculate observation parameters
    exptime, tottime, scanRate, scanHeight, fluence = wfc3_obs(jmag, disperser, scanDirection, subarray,
                                                               nsamp, samp_seq, targetFluence, hmag=hmag)
    if scanHeight > maxScanHeight:
        print(("****WARNING: Computed scan height exceeds maximum recommended height of %0.0f pixels." % maxScanHeight))
    if exptime > maxExptime:
        print(("****WARNING: Computed frame time (%0.0f seconds) exceeds maximum recommended duration of %0.0f seconds." % (exptime, maxExptime)))

    # Compute number of data points (frames) per orbit
    ptsOrbit = np.floor(obsTime/tottime)
    # First point (frame) is always low, ignore when computing duty cycle
    dutyCycle = (exptime*(ptsOrbit-1))/50./60*100

    # Compute number of non-destructive reads per orbit
    readsOrbit = ptsOrbit*(nsamp+1)

    # Look for mid-orbit buffer dumps
    if (subarray == 'grism256') and (readsOrbit >= 300) and (exptime <= 43):
        print(
            "****WARNING: Observing plan may incur mid-orbit buffer dumps.  Check with APT.")
    if (subarray == 'grism512') and (readsOrbit >= 120) and (exptime <= 100):
        print(
            "****WARNING: Observing plan may incur mid-orbit buffer dumps.  Check with APT.")

    # Compute number of HST orbits per transit
    # ~96 minutes per HST orbit
    orbitsTr = trdur*24.*60/96

    # Estimate number of good points during planet transit
    # First point in each HST orbit is flagged as bad; therefore, subtract from total

    if orbitsTr < 0.5:
        # Entire transit fits within one HST orbit
        ptsInTr = ptsOrbit * orbitsTr/0.5 - 1
    elif orbitsTr <= 1.5:
        # Assume one orbit centered on mid-transit time
        ptsInTr = ptsOrbit - 1
    elif orbitsTr < 2.:
        # Assume one orbit during full transit and one orbit during ingress/egress
        ptsInTr = ptsOrbit * \
            (np.floor(orbitsTr) +
             np.min((1, np.remainder(orbitsTr-np.floor(orbitsTr)-0.5, 1)/0.5))) - 2
    else:
        # Assume transit contains 2+ orbits timed to maximize # of data points.
        ptsInTr = ptsOrbit * (np.floor(orbitsTr) + np.min(
            (1, np.remainder(orbitsTr-np.floor(orbitsTr), 1)/0.5))) - np.ceil(orbitsTr)

    # Estimate number of good points outside of transit
    # Discard first HST orbit
    ptsOutTr = (ptsOrbit-1) * (norbits-1) - ptsInTr

    # Compute transit depth uncertainty per spectrophotometric channel
    ratio = 10**((refmag - jmag)/2.5)
    flux = ratio*refflux*exptime/refexptime
    fluxvar = ratio*refvar*exptime/refexptime
    chanflux = flux/nchan
    chanvar = fluxvar/nchan
    chanrms = np.sqrt(chanvar)/chanflux*1e6  # ppm
    inTrrms = chanrms/np.sqrt(ptsInTr*numTr)  # ppm
    outTrrms = chanrms/np.sqrt(ptsOutTr*numTr)  # ppm
    deptherr = np.sqrt(inTrrms**2 + outTrrms**2)  # ppm

    info = {"Number of HST orbits": norbits,
            "Use first orbit":   useFirstOrbit,
            "WFC3 parameters: NSAMP": nsamp,
            "WFC3 parameters: SAMP_SEQ": samp_seq.upper(),
            "Scan Direction": scanDirection,
            "Recommended scan rate (arcsec/s)": scanRate,
            "Scan height (pixels)": scanHeight,
            "Maximum pixel fluence (electrons)": fluence,
            "exposure time": exptime,
            "Estimated duty cycle (outside of Earth occultation)": dutyCycle,
            "Transit depth uncertainty(ppm)": deptherr,
            "Number of channels": nchan,
            "Number of Transits": numTr}

    return {"spec_error": deptherr/1e6,
            "light_curve_rms": chanrms/1e6,
            "nframes_per_orb": ptsOrbit,
            "info": info}


def calc_start_window(eventType, rms, ptsOrbit, numOrbits, depth, inc, aRs, period, windowSize, ecc=0, w=90., duration=None, offset=0., useFirstOrbit=False):
    '''Calculate earliest and latest start times

    Plot earliest and latest possible spectroscopic light curves for given start window size

    Parameters
    ----------
    eventType : str
        'transit' or 'eclipse'
    rms : float
        light curve root-mean-square
    ptsOrbit : int
        number of frames per HST orbit
    numOrbits :  int
        number of HST orbits per visit
    depth : float
        transit/eclipse depth
    inc : float
        orbital inclination in degrees
    aRs : float
        Semi-major axis in units of stellar radii (a/R*)
    period : float
        orbital period in days
    windowSize : float
        observation start window size in minutes
    ecc : float
        (Optional) eccentricity
    w : float
        (Optional) longitude of periastron in degrees
    duration : float
        (Optional) full transit/eclipse duration in days
    offset : float
        (Optional) manual offset in observation start time, in minutes
    useFirstOrbit : bool
        (Optional) whether to use first orbit

    Returns
    -------
    float
        minphase--earliest observation start phase
    float
        maxphase--latest observation start phase
    '''
    ptsOrbit  = int(ptsOrbit)
    numOrbits = int(numOrbits)

    hstperiod = 96./60/24                 # HST orbital period, days
    punc = windowSize/120./24/period  # Half start window size, in phase
    cosi = np.cos(inc*np.pi/180)     # Cosine of the inclination
    rprs = np.sqrt(depth)            # Planet-star radius ratio

    params = batman.TransitParams()
    if eventType == 'transit':
        midpt = period
        b = aRs*cosi*(1-ecc**2)/(1+ecc*np.sin(w*np.pi/180))  # Impact parameter
        # Account for planet speed on eccentric orbits
        sfactor = np.sqrt(1-ecc**2)/(1+ecc*np.sin(w*np.pi/180))
        # limb darkening coefficients
        params.u = [0.1, 0.1]
    elif eventType == 'eclipse':
        #midpt = period/2*(1+4*ecc*np.cos(w*np.pi/180)/np.pi)
        midpt = calculate_tsec(period, ecc, w*np.pi/180, inc*np.pi/180, t0=0, tperi=None, winn_approximation=False)
        b = aRs*cosi*(1-ecc**2)/(1-ecc*np.sin(w*np.pi/180))  # Impact parameter
        # Account for planet speed on eccentric orbits
        sfactor = np.sqrt(1-ecc**2)/(1-ecc*np.sin(w*np.pi/180))
        # limb darkening coefficients
        params.u = [0.0, 0.0]
    else:
        print(("****HALTED: Unknown event type: %s" % eventType))
        return
    params.t0 = midpt/period          # phase of transit/eclipse
    params.per = 1.                    # orbital period, units are orbital phase
    # planet radius (in units of stellar radii)
    params.rp = rprs
    # semi-major axis (in units of stellar radii)
    params.a = aRs
    params.inc = inc                   # orbital inclination (in degrees)
    params.ecc = ecc                   # eccentricity
    params.w = w                     # longitude of periastron (in degrees)
    params.limb_dark = "quadratic"           # limb darkening model

    # Transit/eclipse duration (in days)
    if duration == None:
        duration = period/np.pi * \
            np.arcsin(
                1./aRs*np.sqrt(((1+rprs)**2-(aRs*cosi)**2)/(1-cosi**2)))*sfactor
    phase1 = (midpt + duration/2. - hstperiod*(numOrbits-2 -
                                               useFirstOrbit) - hstperiod/2 + offset/24./60)/period
    phase2 = (midpt - duration/2. - hstperiod*2 + offset/24./60)/period
    minphase = (phase1+phase2)/2-punc
    maxphase = (phase1+phase2)/2+punc

    # Compute light curves at extremes of HST start window
    npts = int(4 * ptsOrbit * numOrbits)
    phdur = duration/period
    
    phase1 = np.linspace(minphase+(1-useFirstOrbit)*hstperiod/period,
                         minphase+hstperiod/period*(numOrbits-1)+hstperiod/period/2, int(npts))
    phase2 = np.linspace(maxphase+(1-useFirstOrbit)*hstperiod/period,
                         maxphase+hstperiod/period*(numOrbits-1)+hstperiod/period/2, int(npts))
    m = batman.TransitModel(params, phase1)
    trmodel1 = m.light_curve(params)
    m = batman.TransitModel(params, phase2)
    trmodel2 = m.light_curve(params)
    obsphase1 = []
    obsphase2 = []
    for i in range(numOrbits):
        obsphase1 = np.r_[obsphase1, np.linspace(
            minphase+hstperiod/period*i, minphase+hstperiod/period*i+hstperiod/period/2, int(ptsOrbit))]
        obsphase2 = np.r_[obsphase2, np.linspace(
            maxphase+hstperiod/period*i, maxphase+hstperiod/period*i+hstperiod/period/2, int(ptsOrbit))]
    m = batman.TransitModel(params, obsphase1)
    obstr1 = m.light_curve(params) + np.random.normal(0, rms, obsphase1.shape)
    m = batman.TransitModel(params, obsphase2)
    obstr2 = m.light_curve(params) + np.random.normal(0, rms, obsphase2.shape)

    return {'obsphase1': obsphase1, 'light_curve_rms': rms, 'obstr1': obstr1, 'obsphase2': obsphase2,
            'obstr2': obstr2, 'minphase': minphase, 'maxphase': maxphase, 'phase1': phase1,
            'phase2': phase2, 'trmodel1': trmodel1, 'trmodel2': trmodel2, 'eventType': eventType, 'planet period':period}


def planet_spec(planet, star, w_unit, disperser, deptherr, nchan, smooth=None):
    '''Plot exoplanet transmission/emission spectrum

    Parameters
    ----------
    planet: dict
        planet dictionary from exo_input
    star : dict
        star dictionary from exo_input
    w_unit : str
        wavelength unit (um or nm)
    disperser :
        grism (g102 or g141)
    deptherr : float
        simulated transit/eclipse depth uncertainty
    nchan : float
        number of spectrophotometric channels
    smooth : float
        (Optional)length of smoothing kernel

     Returns
     -------
     dict
        contains following keys {'model_wave','model_spec','binwave','binspec',
        'error','wmin','wmax'}
     '''
    # Load model wavelengths and spectrum
    mwave, mspec = hst_spec(planet, star)  # np.loadtxt(specfile, unpack=True)
    # Convert wavelength to microns
    # if w_unit == 'um':
    #    pass
    # elif w_unit == 'nm':
    #    mwave /= 1000.
    # else:
    #    print(("****HALTED: Unrecognized wavelength unit: '%s'" % w_unit))
    #    return

    # Smooth model spectrum (optional)
    if smooth != None:
        from .hst_smooth import smooth as sm
        mspec = sm(mspec, smooth)

    # Determine disperser wavelength boundaries
    if disperser == 'g141':
        wmin = 1.125
        wmax = 1.650
    elif disperser == 'g102':
        wmin = 0.84
        wmax = 1.13
    else:
        print(("****HALTED: Unrecognized disperser name: '%s'" % disperser))
        return

    # Determine wavelength bins
    binsize = (wmax - wmin)/nchan
    wave_low = np.round([i for i in np.linspace(wmin, wmax-binsize, nchan)], 3)
    wave_hi = np.round([i for i in np.linspace(wmin+binsize, wmax, nchan)], 3)
    binwave = (wave_low + wave_hi)/2.

    # Create simulated spectrum by binning model spectrum and addding uncertainty
    binspec = np.zeros(nchan)
    for i in range(nchan):
        ispec = np.where((mwave >= wave_low[i])*(mwave <= wave_hi[i]))
        binspec[i] = np.mean(mspec[ispec])
    binspec += np.random.normal(0, deptherr, nchan)

    return {'model_wave': mwave, 'model_spec': mspec, 'binwave': binwave, 'binspec': binspec,
            'error': deptherr, 'wmin': wmin, 'wmax': wmax}


def compute_sim_lightcurve(exposureDict, lightCurveDict, calRamp=False):
    """Compute simulated HST light curves

    Function to take the fluence and error estimates from wfc3_TExoNS
    and the model light curve from calc_start_window to simulate observed
    light curves. Ramp effect systemacts simulated by RECTE can be
    included by turning on the calRamp switch

    Parameters
    ----------
    exposureDict : dict
    Includes information for the observation. This dictionary is
    returned by function wfc3_TExoNS. Relevant keys in the dictionary
    are ['info']["Maximum pixel fluence (electrons)"] and
    ['info']['exposure time']

    lightCurveDict : dict
    Includes the model light curve. The light curves are estimed by
    calc_start_window

    calRamp : bool (default: False)
    Switch to turn on/off RECTE ramp calculation. Calculate realistic
    ramp effect systematics when the switch is turned on.

    Returns
    -------
    dict
    Resulting light curve (unit: e/pixel) for the earliest and latest time

    """
    fluence = exposureDict['info']["Maximum pixel fluence (electrons)"]
    exptime = exposureDict['info']['exposure time']
    obst1 = (lightCurveDict['obsphase1'] -
          lightCurveDict['obsphase1'][0]) *\
          lightCurveDict['planet period'] * 86400 # in seconds
    obst2 = (lightCurveDict['obsphase2'] -
          lightCurveDict['obsphase2'][0]) *\
          lightCurveDict['planet period'] * 86400 # in seconds
    counts1 = lightCurveDict['obstr1'] * fluence
    counts2 = lightCurveDict['obstr2'] * fluence
    # use RECTE to calculate the ramp if calRamp option is turned on
    if calRamp:
        counts1 = RECTE(counts1 / exptime,
                       obst1,
                       exptime)
        counts2 = RECTE(counts2 / exptime,
                       obst2,
                       exptime)
    count_noise = lightCurveDict['light_curve_rms'] * fluence
    resultDict = lightCurveDict.copy()
    resultDict['counts1'] = counts1
    resultDict['counts2'] = counts2
    resultDict['count_noise'] = count_noise
    resultDict['ramp_included'] = calRamp
    resultDict['model_counts1'] = resultDict['trmodel1'] * fluence
    resultDict['model_counts2'] = resultDict['trmodel2'] * fluence
    return resultDict


def compute_sim_hst(dictinput,verbose=False):
    """Sets up HST simulations

    Function to set up explanet observations for HST only and
    compute simulated spectrum.

    Parameters
    ----------
    dictinput : dict
        instrument and pandexo dictionaries in format {"pandeia_input":dict1, "pandexo_input":dict2}

    Returns
    -------
    dict
        All hst output info needed to plot simulated data, light curves timing info
    """
    pandexo_input = dictinput['pandexo_input']
    pandeia_input = dictinput['pandeia_input']

    disperser = pandeia_input['configuration']['instrument']['disperser'].lower()

    # add a switch for ramp calculation
    calRamp = pandeia_input['strategy']['calculateRamp']

    if not pandeia_input['strategy']['useFirstOrbit']:
        if verbose:print("Dropping first orbit designed by observation strategy")
        if verbose:print("Do not calculate ramp profile")
        calRamp = False
    numorbits = pandeia_input['strategy']['norbits']
    nchan = pandeia_input['strategy']['nchan']
    windowSize = pandeia_input['strategy']['windowSize']
    useFirstOrbit = pandeia_input['strategy']['useFirstOrbit']
    calc_type = pandexo_input['planet']['type']
    jmag = pandexo_input['star']['jmag']
    hmag = pandexo_input['star']['hmag']
    specfile = pandexo_input['planet']['exopath']
    w_unit = pandexo_input['planet']['w_unit']
    f_unit = pandexo_input['planet']['f_unit']
    depth = pandexo_input['planet']['depth']
    inc = pandexo_input['planet']['i']
    aRs = pandexo_input['planet']['ars']
    period = pandexo_input['planet']['period']
    ecc = pandexo_input['planet']['ecc']
    w = pandexo_input['planet']['w']
    try:
        offset = pandeia_input['strategy']['offset']
    except:
        offset = 0.
    
    # check to see if ecc or w was provided
    if (type(ecc) != float) and (type(ecc) != int):
        ecc = 0.0
    if (type(w) != float) and (type(w) != int):
        w = 90.0
    if (type(windowSize) != float) and (type(windowSize) != int):
        windowSize = 20.0
    if f_unit == "rp^2/r*^2":
        eventType = 'transit'
    elif f_unit == "fp/f*":
        eventType = 'eclipse'
    elif calc_type == 'grid':
        eventType = 'transit'
    else:
        raise Exception('Units are not correct. Pick rp^2/r*^2 or fp/f*')

    a = wfc3_TExoNS(dictinput)

    b = calc_start_window(eventType, a['light_curve_rms'], a['nframes_per_orb'], a['info']['Number of HST orbits'],
        depth, inc, aRs, period, windowSize, ecc, w, useFirstOrbit=useFirstOrbit, offset=offset)
    c = planet_spec(pandexo_input['planet'], pandexo_input['star'],
                    w_unit, disperser, a['spec_error'], nchan, smooth=20)
    info_div = create_out_div(a['info'], b['minphase'], b['maxphase'])
    simLightCurve = compute_sim_lightcurve(a, b, calRamp=calRamp)
    return {"wfc3_TExoNS": a,
            "calc_start_window": b,
            "planet_spec": c,
            "light_curve": simLightCurve,
            "info_div": info_div}


def create_out_div(input_dict, minphase, maxphase):
    """Function to render input dicts in html format for web front end

    Parameters
    ----------
    input_dict : dict
        any input dictionary

    Returns
    -------
    div
        html rendered table
    """
    input_dict['Start observations between orbital phases'] = str(
        minphase)+'-'+str(maxphase)
    input_div = pd.DataFrame.from_dict(input_dict, orient='index')
    input_div.columns = ['Value']
    input_div = input_div.to_html()
    input_div = '<table class="table table-striped"> \n' + \
        input_div[36:len(input_div)]
    input_div = input_div.encode()
    return input_div

def calculate_tsec(period, ecc, omega, inc, t0 = None, tperi = None, winn_approximation = False):
    ''' Function to calculate the time of secondary eclipse.

        This uses Halley's method (Newton-Raphson, but using second derivatives) to first find the true anomaly (f) at which secondary eclipse occurs,
        then uses this to get the eccentric anomaly (E) at secondary eclipse, which gives the mean anomaly (M) at secondary
        eclipse using Kepler's equation. This finally leads to the time of secondary eclipse using the definition of the mean
        anomaly (M = n*(t - tau) --- here tau is the time of pericenter passage, n = 2*pi/period the mean motion).
        Time inputs can be either the time of periastron passage directly or the time of transit center. If the latter, the
        true anomaly for primary transit will be calculated using Halley's method as well, and this will be used to get the
        time of periastron passage.

        Parameters
        ----------
        period : float
            The period of the transit in days.
        ecc : float
            Eccentricity of the orbit
        omega : float
            Argument of periastron passage (in radians)
        inc : string
            Inclination of the orbit (in radians)
        t0 : float
            The transit time in BJD or HJD (will be used to get time of periastron passage).
        tperi : float
            The time of periastron passage in BJD or HJD (needed if t0 is not supplied).
        winn_approximation : boolean
            If True, the approximation in Winn (2010) is used --- (only valid for not very eccentric and inclined orbits).
        Returns
        -------
        tsec : float
            The time of secondary eclipse '''
    # Check user is not playing trick on us:
    if period<0:
        raise Exception('Period cannot be a negative number.')
    if ecc<0 or ecc>1:
        raise Exception('Eccentricity (e) is out of bounds (0 < e < 1).')

    # Use true anomaly approximation given in Winn (2010) as starting point:
    f_occ_0 = (-0.5*np.pi) - omega
    if not winn_approximation:
        f_occ = optimize.newton(drsky, f_occ_0, fprime = drsky_prime, fprime2 = drsky_2prime, args = (ecc, omega, inc,))
    else:
        f_occ = f_occ_0
    # Define the mean motion, n:
    n = 2.*np.pi/period

    # If time of transit center is given, use it to calculate the time of periastron passage. If no time of periastron
    # or time-of-transit center given, raise error:
    if tperi is None:
        # For this, find true anomaly during transit. Use Winn (2010) as starting point:
        f_tra_0 = (np.pi/2.) - omega
        if not winn_approximation:
            f_tra = optimize.newton(drsky, f_tra_0, fprime = drsky_prime, fprime2 = drsky_2prime, args = (ecc, omega, inc,))
        else:
            f_tra = f_tra_0
        # Get eccentric anomaly during transit:
        E = getE(f_tra, ecc)

        # Get mean anomaly during transit:
        M = getM(E, ecc)

        # Get time of periastron passage from mean anomaly definition:
        tperi = t0 - (M/n)

    elif (tperi is None) and (t0 is None):
        raise ValueError('The time of periastron passage or time-of-transit center has to be supplied for the calculation to work.')

    # Get eccentric anomaly:
    E = getE(f_occ, ecc)

    # Get mean anomaly during secondary eclipse:
    M = getM(E, ecc)

    # Get the time of secondary eclipse using the definition of the mean anomaly:
    tsec = (M/n) + tperi

    # Note returned time-of-secondary eclipse is the closest to the time of periastron passage and/or time-of-transit center. Check that
    # the returned tsec is the *next* tsec to the time of periastron or t0 (i.e., the closest *future* tsec):
    if t0 is not None:
        tref = t0
    else:
        tref = tperi
    if tref > tsec:
        while True:
            tsec += period
            if tsec > tref:
                break
    return tsec

def drsky_2prime(x, ecc, omega, inc):
    ''' Second derivative of function drsky. This is the second derivative with respect to f of the drsky function.
    Parameters
    ----------
    x : float
      True anomaly
    ecc : float
      Eccentricity of the orbit
    omega : float
      Argument of periastron passage (in radians)
    inc : float
      Inclination of the orbit (in radians)

    Returns
    -------
    drsky_2prime : float
      Function evaluated at x, ecc, omega, inc'''

    sq_sini = np.sin(inc)**2
    sin_o_p_f = np.sin(x+omega)
    cos_o_p_f = np.cos(x+omega)
    ecosf = ecc*np.cos(x)
    esinf = ecc*np.sin(x)

    f1 = esinf - esinf*sq_sini*(sin_o_p_f**2)
    f2 = -sq_sini*(ecosf + 4.)*(sin_o_p_f*cos_o_p_f)

    return f1+f2

def drsky_prime(x, ecc, omega, inc):
    ''' Derivative of function drsky. This is the first derivative with respect to f of the drsky function.
    Parameters
    ----------
    x : float
      True anomaly
    ecc : float
      Eccentricity of the orbit
    omega : float
      Argument of periastron passage (in radians)
    inc : float
      Inclination of the orbit (in radians)

    Returns
    -------
    drsky_prime : float
      Function evaluated at x, ecc, omega, inc'''

    sq_sini = np.sin(inc)**2
    sin_o_p_f = np.sin(x+omega)
    cos_o_p_f = np.cos(x+omega)
    ecosf = ecc*np.cos(x)
    esinf = ecc*np.sin(x)

    f1 = (cos_o_p_f**2 - sin_o_p_f**2)*(sq_sini)*(1. + ecosf)
    f2 = -ecosf*(1 - (sin_o_p_f**2)*(sq_sini))
    f3 = esinf*sin_o_p_f*cos_o_p_f*sq_sini

    return f1+f2+f3

def drsky(x, ecc, omega, inc):
    ''' Function whose roots we wish to find to obtain time of secondary (and primary) eclipse(s)
    When one takes the derivative of equation (5) in Winn (2010; https://arxiv.org/abs/1001.2010v5), and equates that to zero (to find the
    minimum/maximum of said function), one gets to an equation of the form g(x) = 0. This function (drsky) is g(x), where x is the true
    anomaly.
    Parameters
    ----------
    x : float
      True anomaly
    ecc : float
      Eccentricity of the orbit
    omega : float
      Argument of periastron passage (in radians)
    inc : float
      Inclination of the orbit (in radians)

    Returns
    -------
    drsky : float
      Function evaluated at x, ecc, omega, inc '''

    sq_sini = np.sin(inc)**2
    sin_o_p_f = np.sin(x+omega)
    cos_o_p_f = np.cos(x+omega)

    f1 = sin_o_p_f*cos_o_p_f*sq_sini*(1. + ecc*np.cos(x))
    f2 = ecc*np.sin(x)*(1. - sin_o_p_f**2 * sq_sini)
    return f1 - f2

def getE(f,ecc):
    """ Function that returns the eccentric anomaly
    Note normally this is defined in terms of cosines (see, e.g., Section 2.4 in Murray and Dermott), but numerically
    this is troublesome because the arccosine doesn't handle negative numbers by definition (equation 2.43). That's why
    the arctan version is better as signs are preserved (derivation is also in the same section, equation 2.46).
    Parameters
    ----------
    f : float
      True anomaly
    ecc : float
      Eccentricity
    Returns
    -------
    E : float
      Eccentric anomaly """

    return 2. * np.arctan(np.sqrt((1.-ecc)/(1.+ecc))*np.tan(f/2.))

def getM(E, ecc):
    """ Function that returns the mean anomaly using Kepler's equation
    Parameters
    ----------
    E : float
      Eccentric anomaly
    ecc: float
      Eccentricity
    Returns
    -------
    M : float
      Mean anomaly """

    return E - ecc*np.sin(E)
