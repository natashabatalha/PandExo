
import numpy as np

"""
ToDo:
    write another function to predict (nsamp, samp_seq, norbits) if values aren't given
"""

def wfc3_GuessParams(hmag, disperser, scanDirection, subarray, obsTime, maxScanHeight=180., maxExptime=150.):
    '''
    Predict nsamp and samp_seq when values not provided by the user.
    
    PARAMETERS
    ----------
    hmag            : float, H-band magnitude
    maxScanHeight   : (optional) float, maximum scan height in pixels
    
    OUTPUT
    ------
    nsamp       : float, number of up-the-ramp samples (1..15)
    samp_seq    : string, time between non-destructive reads
    
    HISTORY
    -------
    Written by Kevin Stevenson      October 2016
    '''
    #Assume H-band
    '''
    # Recommended scan rate
    scanRate    = np.round(2*10**(-0.4*(hmag-5.9)),3)   #arcsec/s
    scanRatep   = scanRate/0.121                        #pixels/s
    # Calculate max exposure time
    # Default max scan height is 180 pixels (256 subarray)
    maxexptime  = maxScanHeight/scanRatep
    # Calculate number of samples in each mode
    # Assumes GRISM256, but GRISM512 is similar
    nsamp5      = np.min([15,np.floor((maxexptime - 0.280)/2.349)])     #SPARS5
    nsamp10     = np.min([15,np.floor((maxexptime - 0.278)/7.3465)])    #SPARS10
    nsamp25     = np.min([15,np.floor((maxexptime - 0.278)/22.346)])    #SPARS25
    #print(nsamp5,nsamp10,nsamp25)
    
    if nsamp25 >= 4:
        samp_seq    = 'spars25'
        nsamp       = nsamp25
    elif nsamp10 >= 2:
        samp_seq    = 'spars10'
        nsamp       = nsamp10
    else:
        samp_seq    = 'spars5'
        nsamp       = nsamp5
    '''
    allnsamp    = np.arange(1,16)
    allsampseq  = ['spars5', 'spars10', 'spars25']
    maxDutyCycle= 0
    for samp_seq in allsampseq:
        for nsamp in allnsamp:
            exptime, tottime, scanRate, scanHeight, fluence = wfc3_obs(hmag, disperser, scanDirection, 
                                                                       subarray, nsamp, samp_seq)
            # Compute duty cycle and compare
            # Exposure time should be less than 2.5 minutes to achieve good time resolution
            ptsOrbit    = np.floor(obsTime/tottime)
            dutyCycle   = (exptime*(ptsOrbit-1))/50./60*100
            if (dutyCycle > maxDutyCycle) and (exptime < maxExptime):
                maxDutyCycle    = dutyCycle
                bestsampseq     = samp_seq
                bestnsamp       = nsamp
    
    return bestnsamp, bestsampseq

def wfc3_GuessNOrbits(trdur):
    '''
    Predict number of HST orbits for transit observation when not provided by the user.
    
    PARAMETERS
    ----------
    trdur           : array, transit duration in seconds
    
    OUTPUT
    ------
    norbits     : float, number of requested orbits per transit (including discarded thermal-settling orbit)
    
    HISTORY
    -------
    Written by Kevin Stevenson      October 2016
    '''
    # Compute # of HST orbits during transit
    # ~96 minutes per HST orbit
    orbitsTr    = trdur/96./60.
    if orbitsTr <= 1.5:
        norbits = 4.
    elif orbitsTr <= 2.0:
        norbits = 5.
    else:
        norbits = np.ceil(orbitsTr*2+1)
    
    return norbits

def wfc3_obs(hmag, disperser, scanDirection, subarray, nsamp, samp_seq):
    '''
    Determine the recommended exposure time, scan rate, scan height, and overheads.
    '''
    
    # Estimate exposure time
    if subarray == 'grism512':
        if samp_seq == 'spars5':
            exptime     = 0.853 + (nsamp-1)*2.9215   #SPARS5
        elif samp_seq == 'spars10':
            exptime     = 0.853 + (nsamp-1)*7.9217   #SPARS10
        elif samp_seq == 'spars25':
            exptime     = 0.853 + (nsamp-1)*22.9213  #SPARS25
        else:
            print("Unknown SAMP_SEQ: %s" % samp_seq)
            return
    else:
        #GRISM256
        if samp_seq == 'spars5':
            exptime     = 0.280 + (nsamp-1)*2.349   #SPARS5
        elif samp_seq == 'spars10':
            exptime     = 0.278 + (nsamp-1)*7.3465  #SPARS10
        elif samp_seq == 'spars25':
            exptime     = 0.278 + (nsamp-1)*22.346  #SPARS25
        else:
            print("Unknown SAMP_SEQ: %s" % samp_seq)
            return
    
    # Recommended scan rate
    scanRate    = np.round(1.9*10**(-0.4*(hmag-5.9)),3)   #arcsec/s
    # Max fluence in electrons/pixel
    fluence     = (5.5/scanRate)*10**(-0.4*(hmag-15))*2.4   #electrons
    if disperser == 'g102':
        #WFC3_ISR_2012-08 indicates that G102/G141 scale factor is 0.96 DN/electron
        fluence *= 0.96
    # Scan height in pixels
    scanRatep   = scanRate/0.121                        #pixels/s
    scanHeight  = scanRatep*exptime
    '''
    #Quadratic correlation between scanRate and read overhead
    foo     = np.array([0.0,0.1,0.3,0.5,1.0,2.0,3.0,4.0])/0.121
    bar     = np.array([ 40, 40, 41, 42, 45, 50, 56, 61])
    c       = np.polyfit(foo,bar,2)
    model   = c[2] + c[1]*foo + c[0]*foo**2
    #c = [  6.12243227e-04,   6.31621064e-01,   3.96040946e+01]
    '''
    # Define instrument overheads (in seconds)
    read        = 39.60 + 6.316e-1*scanRatep + 6.122e-4*scanRatep**2
    #Correlation between scanHeight/scanRate and pointing overhead was determined elsewhere
    if scanDirection == 'Round Trip':
        #Round Trip scan direction doesn't have to return to starting point, therefore no overhead
        pointing    = 0.
    elif scanDirection == 'Forward':
        c0,c1,c2,c3,c4,c5= [  3.18485340e+01,   3.32968829e-02,   1.65687590e-02,   
                              7.65510038e-01,  -6.24504499e+01,   5.51452028e-03]
        pointing    = c0*(1 - np.exp(-c2*(scanHeight-c4))) + c1*scanHeight + c3*scanRatep +c5*scanRatep**2
    else:
        print("Unknown scan direction: %s" % scanDirection)
        return
    tottime     = exptime+read+pointing         #seconds
    #print("read",read)
    #print("pointing",pointing)
    
    return exptime, tottime, scanRate, scanHeight, fluence

def wfc3_TExoNS(hmag, trdur, numTr, nchan, disperser, scanDirection='Forward', subarray='grism256', nsamp=0, samp_seq=None, norbits=None, schedulability='100'):
    '''
    Compute the transit depth uncertainty for a defined system and number of spectrophotometric channels.
    
    PARAMETERS
    ----------
    hmag            : float, H-band magnitude
    trdur           : array, transit duration in seconds
    numTr           : array, number of observed transits/eclipses
    nchan           : float, number of spectrophotometric channels
    disperser       : string, grism (g102 or g141)
    scanDirection   : string, spatial scan direction (Forward or Round Trip)
    subarray        : string, subarray aperture (grism256 or grism512)
    nsamp           : float, number of up-the-ramp samples (1..15)
    samp_seq        : string, time between non-destructive reads
    norbits         : float, number of requested orbits per transit (including discarded thermal-settling orbit)
    schedulability  : string, % time HST can observe target (30 or 100)
    
    OUTPUT
    ------
    deptherr    : float, transit depth uncertainty per spectrophotometric channel
    
    HISTORY
    -------
    Written by Kevin Stevenson      October 2016
    '''
    #Assumptions: H-band
    disperser   = disperser.lower()
    subarray    = subarray.lower()
    if isinstance(samp_seq, str):
        samp_seq    = samp_seq.lower()
    
    if disperser == 'g141':
        # Define reference Hmag, flux, variance, and exposure time for GJ1214
        refmag      = 9.094
        refflux     = 2.32e8
        refvar      = 2.99e8
        refexptime  = 88.436
    elif disperser == 'g102':
        # Define reference Hmag, flux, variance, and exposure time for WASP12
        refmag      = 10.228
        refflux     = 8.26e7
        refvar      = 9.75e7
        refexptime  = 103.129
    else:
        print("Unknown disperser: %s" % disperser)
        return
    
    # Determine max recommended scan height
    if subarray == 'grism512':
        maxScanHeight = 430
    elif subarray == 'grism256':
        maxScanHeight = 180
    else:
        print("Unknown subarray aperture: %s" % subarray)
        return
    
    # Define maximum frame time
    maxExptime  = 150.
    
    # Estimate reasonable number of HST orbits
    guessorbits = wfc3_GuessNOrbits(trdur)
    if norbits == None:
        norbits = guessorbits
    elif norbits != guessorbits:
        print("****WARNING: Number of specified HST orbits does not match number of recommended orbits: %0.0f" % guessorbits)
    
    # Determine seconds of observing time per HST orbit
    if schedulability == '30':
        obsTime   = 51.3*60
    elif schedulability == '100':
        obsTime   = 46.3*60
    else:
        print("Unknown schedulability: %s" % schedulability)
        return
    
    if nsamp == 0 or samp_seq == None or samp_seq == "None":
        # Estimate reasonable values
        nsamp, samp_seq = wfc3_GuessParams(hmag, disperser, scanDirection, subarray, 
                                           obsTime, maxScanHeight, maxExptime)
    # Calculate observation parameters
    exptime, tottime, scanRate, scanHeight, fluence = wfc3_obs(hmag, disperser, scanDirection, 
                                                               subarray, nsamp, samp_seq)
    if scanHeight > maxScanHeight:
        print("****WARNING: Computed scan height exceeds maximum recommended height of %0.0f pixels." % maxScanHeight)
    if exptime > maxExptime:
        print("****WARNING: Computed frame time (%0.0f seconds) exceeds maximum recommended duration of %0.0f seconds." % (exptime,maxExptime))
    
    # Compute # of data points per orbit
    # First point is always bad, but take into account later
    ptsOrbit    = np.floor(obsTime/tottime)
    dutyCycle   = (exptime*(ptsOrbit-1))/50./60*100
    
    # Compute # of non-destructive reads per orbit
    readsOrbit  = ptsOrbit*(nsamp+1)
    
    # Look for mid-orbit buffer dumps
    if (subarray == 'grism256') and (readsOrbit >= 300) and (exptime <= 43):
        print("****WARNING: Observing plan may incur mid-orbit buffer dumps.  Check with APT.")
    if (subarray == 'grism512') and (readsOrbit >= 120) and (exptime <= 100):
        print("****WARNING: Observing plan may incur mid-orbit buffer dumps.  Check with APT.")
    
    # Compute # of HST orbits per transit
    # ~96 minutes per HST orbit
    orbitsTr    = trdur/96./60.
    
    # Estimate # of good points during planet transit
    # First point in each HST orbit is flagged as bad; therefore, subtract from total
    #foo = []
    #bar = np.arange(0.1,4,0.05)
    #for orbitsTr in bar:
    if orbitsTr < 0.5:
        # Entire transit fits within one HST orbit
        ptsInTr  = ptsOrbit * orbitsTr/0.5 - 1
    elif orbitsTr <= 1.5:
        # Assume one orbit centered on mid-transit time
        ptsInTr = ptsOrbit - 1
    elif orbitsTr < 2.:
        # Assume one orbit during full transit and one orbit during ingress/egress
        ptsInTr = ptsOrbit * (np.floor(orbitsTr) + np.min((1,np.remainder(orbitsTr-np.floor(orbitsTr)-0.5,1)/0.5))) - 2
    else:
        # Assume transit contains 2+ orbits timed to maximize # of data points.
        ptsInTr  = ptsOrbit * (np.floor(orbitsTr) + np.min((1,np.remainder(orbitsTr-np.floor(orbitsTr),1)/0.5))) - np.ceil(orbitsTr)
    #foo.append(ptsInTr)
    #plt.figure(1)
    #plt.clf()
    #plt.plot(bar,foo,'o-')
    
    # Estimate # of good points outside of transit
    # Discard first HST orbit
    ptsOutTr    = (ptsOrbit-1) * (norbits-1) - ptsInTr
    
    # Compute transit depth uncertainty per spectrophotometric channel
    ratio       = 10**((refmag - hmag)/2.5)
    flux        = ratio*refflux*exptime/refexptime
    fluxvar     = ratio*refvar*exptime/refexptime
    chanflux    = flux/nchan
    chanvar     = fluxvar/nchan
    chanrms     = np.sqrt(chanvar)/chanflux*1e6     #ppm
    #print("chanrms",chanrms)
    inTrrms     = chanrms/np.sqrt(ptsInTr*numTr)    #ppm
    outTrrms    = chanrms/np.sqrt(ptsOutTr*numTr)   #ppm
    #print("in/out",inTrrms,outTrrms)
    deptherr    = np.sqrt(inTrrms**2 + outTrrms**2) #ppm
    
    print("Number of HST orbits: %0.0f" % norbits)
    print("WFC3 parameters: NSAMP = %0.0f, SAMP_SEQ = %s" %(nsamp,samp_seq.upper()))
    print("Recommended scan rate: %0.3f arcsec/s" % scanRate)
    print("Scan height: %0.1f pixels" % scanHeight)
    print("Maximum pixel fluence: %0.0f electrons" % fluence)
    print("Estimated duty cycle (outside of Earth occultation): %0.1f%%" % dutyCycle)
    print("Transit depth uncertainty: %0.1f ppm for each of %0.0f channel(s)" % (deptherr, nchan))
    
    return deptherr/1e6, chanrms/1e6, ptsOrbit

def calc_StartWindow(rms, ptsOrbit, numorbits, depth, inc, aRs, period, tunc, duration=None, offset=0.):
    '''
    Plot earliest and latest possible spectroscopic light curves for given start window size
    
    offset in minutes
    '''
    import matplotlib.pyplot as plt
    import batman
    
    hstperiod   = 96./60/24             # HST orbital period, days
    gsacq       = 0.                    # Guide star acquisition, days
    punc        = tunc/60./24/period    # Start time uncertainty, orbital period
    cosi        = np.cos(inc*np.pi/180) # Cosine of the inclination
    rprs        = np.sqrt(depth)        # Planet-star radius ratio
    if duration == None:                # Transit duration
        duration    = period/np.pi*np.arcsin(1./aRs*np.sqrt(((1+rprs)**2-(aRs*cosi)**2)/(1-cosi**2)))

    params = batman.TransitParams()
    params.t0 = 1.                       #time of inferior conjunction
    params.per = 1.                      #orbital period
    params.rp = rprs                     #planet radius (in units of stellar radii)
    params.a = aRs                       #semi-major axis (in units of stellar radii)
    params.inc = inc                     #orbital inclination (in degrees)
    params.ecc = 0.                      #eccentricity
    params.w = 90.                       #longitude of periastron (in degrees)
    params.u = [0.0, 0.0]                #limb darkening coefficients
    params.limb_dark = "quadratic"       #limb darkening model

    phase1 = (period + duration/2. - hstperiod*(numorbits-2) - hstperiod/2 + offset/24./60)/period
    phase2 = (period - duration/2. - hstperiod*2 + offset/24./60)/period
    minphase    = (phase1+phase2)/2-punc
    maxphase    = (phase1+phase2)/2+punc
    #print(minphase,maxphase)

    #Plot extremes of possible HST observations
    npts    = 4 * ptsOrbit * numorbits
    phdur   = duration/period
    phase1  = np.linspace(minphase+hstperiod/period,minphase+hstperiod/period*(numorbits-1)+hstperiod/period/2,npts)
    phase2  = np.linspace(maxphase+hstperiod/period,maxphase+hstperiod/period*(numorbits-1)+hstperiod/period/2,npts)
    m = batman.TransitModel(params, phase1)
    trmodel1= m.light_curve(params)
    m = batman.TransitModel(params, phase2)
    trmodel2= m.light_curve(params)
    obsphase1   = []
    obsphase2   = []
    for i in range(numorbits):
        obsphase1 = np.r_[obsphase1, np.linspace(minphase+hstperiod/period*i,minphase+hstperiod/period*i+hstperiod/period/2,ptsOrbit)]
        obsphase2 = np.r_[obsphase2, np.linspace(maxphase+hstperiod/period*i,maxphase+hstperiod/period*i+hstperiod/period/2,ptsOrbit)]
    m = batman.TransitModel(params, obsphase1)
    obstr1= m.light_curve(params) + np.random.normal(0, rms, obsphase1.shape)
    m = batman.TransitModel(params, obsphase2)
    obstr2= m.light_curve(params) + np.random.normal(0, rms, obsphase2.shape)
    
    plt.figure(1, figsize=(12,4))
    plt.clf()
    plt.subplot(121)
    plt.title('Earliest Start Time')
    plt.errorbar(obsphase1, obstr1, rms, fmt='go')
    plt.plot(phase1, trmodel1, 'b-', lw=2)
    ylim1=plt.ylim()
    xlim1=plt.xlim()
    plt.ylabel("Flux")
    plt.xlabel("Orbital Phase")
    plt.subplot(122)
    plt.title('Latest Start Time')
    plt.errorbar(obsphase2, obstr2, rms, fmt='ro')
    plt.plot(phase2, trmodel2, 'b-', lw=2)
    ylim2=plt.ylim()
    xlim2=plt.xlim()
    plt.ylabel("Flux")
    plt.xlabel("Orbital Phase")
    #Put both subplots onto same x,y scale
    ylim=[np.min((ylim1[0],ylim2[0])), np.max((ylim1[1],ylim2[1]))]
    xlim=[np.min((xlim1[0],xlim2[0])), np.max((xlim1[1],xlim2[1]))]
    plt.ylim(ylim)
    plt.xlim(xlim)
    plt.subplot(121)
    plt.ylim(ylim)
    plt.xlim(xlim)
    
    #plt.savefig('PlanSpec.png')
    
    return minphase, maxphase

def plot_PlanSpec(specfile, w_unit, disperser, deptherr, nchan, smooth=None):
    '''
    Plot exoplanet transmission/emission spectrum
    '''
    import matplotlib.pyplot as plt
    import smooth as sm
    # Load model wavelengths and spectrum
    mwave, mspec = np.loadtxt(specfile, unpack=True)
    # Convert wavelength to microns
    if w_unit == 'um':
        pass
    elif w_unit == 'nm':
        mwave /= 1000.
    else:
        print("****WARNGING: Unrecognized wavelength unit '%s'" % w_unit)
        return
    
    # Smooth model spectrum (optional)
    if smooth != None:
        mspec = sm.smooth(mspec, smooth)
    
    # Determine disperser wavelength boundaries
    if disperser == 'g141':
        wmin = 1.125
        wmax = 1.650
    elif disperser == 'g102':
        wmin = 0.8
        wmax = 1.2
    else:
        print("WARNING: Unrecognized disperser name '%s'" % disperser)
    
    # Determine wavelength bins
    binsize     = (wmax - wmin)/nchan
    wave_low    = np.round([i for i in np.linspace(wmin, wmax-binsize, nchan)],3)
    wave_hi     = np.round([i for i in np.linspace(wmin+binsize, wmax, nchan)],3)
    binwave     = (wave_low + wave_hi)/2.
    
    # Create simulated spectrum by binning model spectrum and addding uncertainty
    binspec     = np.zeros(nchan)
    for i in range(nchan):
        ispec       = np.where((mwave >= wave_low[i])*(mwave <= wave_hi[i]))
        binspec[i]  = np.mean(mspec[ispec])
    binspec    += np.random.normal(0,deptherr,nchan)
    
    plt.figure(1, figsize=(12,4))
    plt.clf()
    plt.plot(mwave, mspec, '-k')
    plt.errorbar(binwave, binspec, deptherr, fmt='bo', ms=8, label='Simulated Spectrum')
    plt.xlim(wmin, wmax)
    plt.ylim(np.min(binspec)-2*deptherr, np.max(binspec)+2*deptherr)
    plt.legend(loc='upper left')
    plt.xlabel("Wavelength ($\mu m$)")
    plt.ylabel("Depth")
    
    return
