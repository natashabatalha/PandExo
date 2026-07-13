import os
import sys
import json
import warnings
import numpy as np
import pandas as pd
from copy import deepcopy 
from astropy.io import fits
from . import create_input as create
from .compute_noise import ExtractSpec
import astropy.units as u
import pickle

#constant parameters.. consider putting these into json file 
#max groups in integration
max_ngroup = {'nirspec':65535,
              'miri':65535,
              'niriss':65535,
              'nircam':100}
#minimum number of integrations
min_nint_trans = 3
DHS_F150W_MIN_WAVELENGTH = 0.96
MIRI_LRS_ALLOWED_SUBARRAYS = {
    "lrsslitless": ("slitlessprism", "slitlessprism_ip", "slitlessprism_ips"),
    "lrsslit": ("subslit", "full"),
}
DHS_READOUT_PATTERNS = (
    "rapid", "bright1", "dhs3", "dhs4", "dhs5", "dhs6", "dhs7"
)
DHS_READOUT_CADENCE_FRAMES = {
    "rapid": 1,
    "bright1": 2,
    "dhs3": 3,
    "dhs4": 4,
    "dhs5": 5,
    "dhs6": 6,
    "dhs7": 7,
}
DHS_DATA_EXCESS_RECOMMENDED_LIMIT_GB = 15.0
DHS_DATA_EXCESS_LOWER_THRESHOLD_GB = 5.0
DHS_SUSTAINABLE_DATA_RATE_GB_PER_HOUR = 3.132
DHS_NO_TA_FIXED_OVERHEAD_SECONDS = 694.0
DHS_INITIAL_SLEW_SECONDS = 2100.0
# Peak raw group rates calibrated from STScI's published RAPID, NGROUPS=2
# DHS data-excess rates. These account for each supported DHS subarray size.
DHS_RAW_GROUP_RATE_GB_PER_HOUR = {
    "sub41s1_2-spectra": 13.728,
    "sub82s2_4-spectra": 13.893,
    "sub164s4_8-spectra": 13.983,
    "sub260s4_8-spectra": 14.013,
}

#refdata directory
default_refdata_directory = os.environ.get("pandeia_refdata")


def _instrument_factory(config):
    from pandeia.engine.instrument_factory import InstrumentFactory

    return InstrumentFactory(config=config)


def _perform_calculation(*args, **kwargs):
    from pandeia.engine.perform_calculation import perform_calculation

    return perform_calculation(*args, **kwargs)


def _build_default_calc(*args, **kwargs):
    from pandeia.engine.calc_utils import build_default_calc

    return build_default_calc(*args, **kwargs)

def sort_by_wave_order(value, wave_order):
    if isinstance(value, np.ndarray) and value.shape[:1] == (len(wave_order),):
        return value[wave_order]
    if isinstance(value, list):
        return [sort_by_wave_order(item, wave_order) for item in value]
    if isinstance(value, tuple):
        return tuple(sort_by_wave_order(item, wave_order) for item in value)
    return value


def select_calculation(planet_wave_unit, nsuperstripe, is_dhs=False):
    """Choose the noise calculation from the detector readout mode."""
    use_slope = is_dhs or nsuperstripe > 1
    if planet_wave_unit == 'sec':
        return 'phase_spec_slope' if use_slope else 'phase_spec_fml'
    if use_slope:
        return 'slope method'
    return 'fml'


def nircam_dhs_no_ta_overhead(tframe):
    """Return the assumed no-TA scheduling and initial-slew overhead.

    The fixed 694-second scheduling component follows the documented NIRCam
    overhead model for one grism time-series exposure without target
    acquisition: visit scripts, guide-star acquisition, subarray and wheel
    configuration, OSS compilation, exposure setup/cleanup, fine-guide
    shutdown, and end-of-visit activities. Frame synchronization contributes
    another half frame. The standard 2,100-second initial slew is then added
    because it contributes to the nominal data allocation used by APT.

    https://jwst-docs.stsci.edu/jppom/visit-overheads-timing-model/instrument-specific-overheads/nircam-overheads#gsc.tab=0

    Parameters
    ----------
    tframe : float
        Detector frame time in seconds.

    Returns
    -------
    float
        Scheduling overhead plus initial slew in seconds.
    """
    if not np.isfinite(tframe) or tframe <= 0:
        raise ValueError("tframe must be a positive finite value")
    return (
        DHS_NO_TA_FIXED_OVERHEAD_SECONDS
        + 0.5 * tframe
        + DHS_INITIAL_SLEW_SECONDS
    )


def estimate_dhs_data_excess(
        subarray, readout_pattern, ngroup, exposure_hours,
        allocation_overhead_seconds=0.0):
    """Estimate NIRCam DHS data excess for one APT exposure.

    The estimate follows the data-rate relation underlying Table 2 of the
    STScI NIRCam short-wavelength grism time-series recommendations. It
    subtracts JWST's sustainable 3.132 GB/hour allocation from the generated
    data rate, then scales the excess rate by the full exposure duration.

    Parameters
    ----------
    subarray : str
        PandExo/Pandeia DHS subarray name, such as
        ``"sub260s4_8-spectra"``.
    readout_pattern : str
        One of ``RAPID``, ``BRIGHT1``, or ``DHS3`` through ``DHS7``. Matching
        is case-insensitive.
    ngroup : int
        Number of groups per integration.
    exposure_hours : float
        Elapsed duration of the APT exposure, including integration resets,
        in hours.
    allocation_overhead_seconds : float, optional
        Additional scheduling and slew duration over which APT accrues nominal
        data allocation. The default of zero reproduces the approximate rates
        in the STScI DHS recommendations.

    Returns
    -------
    tuple of float
        Estimated data-excess rate in GB/hour and total data excess in GB.
        A negative rate means the generated rate is below the sustainable
        allocation; the corresponding total data excess is reported as zero.

    Raises
    ------
    ValueError
        If the subarray, readout pattern, number of groups, or duration is not
        supported by the estimate.

    Notes
    -----
    This is an estimate for choosing a viable readout pattern. Users should
    still verify final data-volume constraints in APT.
    """
    subarray = str(subarray).lower()
    readout_pattern = str(readout_pattern).lower()
    if subarray not in DHS_RAW_GROUP_RATE_GB_PER_HOUR:
        raise ValueError(f"Unsupported NIRCam DHS subarray: {subarray}")
    if readout_pattern not in DHS_READOUT_CADENCE_FRAMES:
        raise ValueError(
            f"Unsupported NIRCam DHS readout pattern: {readout_pattern}"
        )
    if int(ngroup) != ngroup or ngroup < 1:
        raise ValueError("ngroup must be a positive integer")
    if exposure_hours < 0:
        raise ValueError("exposure_hours must be non-negative")
    if allocation_overhead_seconds < 0:
        raise ValueError("allocation_overhead_seconds must be non-negative")

    ngroup = int(ngroup)
    cadence_frames = DHS_READOUT_CADENCE_FRAMES[readout_pattern]
    clock_frames = 2 + (ngroup - 1) * cadence_frames
    generated_rate = (
        DHS_RAW_GROUP_RATE_GB_PER_HOUR[subarray]
        * ngroup
        / clock_frames
    )
    excess_rate = generated_rate - DHS_SUSTAINABLE_DATA_RATE_GB_PER_HOUR
    overhead_allocation = (
        DHS_SUSTAINABLE_DATA_RATE_GB_PER_HOUR
        * allocation_overhead_seconds
        / 3600.0
    )
    total_data_excess = excess_rate * exposure_hours - overhead_allocation
    return excess_rate, max(0.0, total_data_excess)


def validate_miri_lrs_subarray(conf):
    """Validate MIRI LRS mode/subarray combinations supported by Pandeia."""
    instrument = conf.get("instrument", {})
    detector = conf.get("detector", {})
    if str(instrument.get("instrument", "")).lower() != "miri":
        return

    mode = str(instrument.get("mode", "")).lower()
    if mode not in MIRI_LRS_ALLOWED_SUBARRAYS:
        return

    subarray = str(detector.get("subarray", "")).lower()
    allowed = MIRI_LRS_ALLOWED_SUBARRAYS[mode]
    if subarray not in allowed:
        allowed_display = ", ".join(item.upper() for item in allowed)
        raise ValueError(
            f"MIRI LRS {mode.upper()} supports only these subarrays: "
            f"{allowed_display}. Got {subarray.upper()}."
        )


def _pandeia_1d_values_at_wave(pand_dict, key, wave):
    """Return a Pandeia 1D diagnostic sampled on PandExo's wavelength grid."""
    pandeia_1d = pand_dict.get('1d') or {}
    if key not in pandeia_1d:
        return np.zeros(len(wave), dtype=float)

    diagnostic_wave = np.asarray(pandeia_1d[key][0], dtype=float)
    diagnostic_value = np.asarray(pandeia_1d[key][1], dtype=float)
    wave = np.asarray(wave, dtype=float)

    if (
        diagnostic_wave.shape == wave.shape
        and np.allclose(diagnostic_wave, wave, equal_nan=True)
    ):
        return diagnostic_value

    finite = np.isfinite(diagnostic_wave) & np.isfinite(diagnostic_value)
    if not np.any(finite):
        return np.zeros(len(wave), dtype=float)

    diagnostic_wave = diagnostic_wave[finite]
    diagnostic_value = diagnostic_value[finite]
    order = np.argsort(diagnostic_wave, kind='mergesort')
    return np.interp(
        wave,
        diagnostic_wave[order],
        diagnostic_value[order],
        left=0.0,
        right=0.0,
    )


def mask_fully_saturated_final_spectrum(finalspec, full_saturation_mask):
    """Set final-spectrum values to NaN where Pandeia reports full saturation."""
    full_saturation_mask = np.asarray(full_saturation_mask, dtype=bool)
    if not np.any(full_saturation_mask):
        finalspec['full_saturation_mask'] = full_saturation_mask
        return finalspec

    for key in ['spectrum', 'spectrum_w_rand', 'error_w_floor']:
        values = np.asarray(finalspec[key], dtype=float).copy()
        values[full_saturation_mask] = np.nan
        finalspec[key] = values

    finalspec['full_saturation_mask'] = full_saturation_mask
    return finalspec


def nirspec_valid_channel_mask(conf, extracted_noise, full_saturation=None):
    """Mask NIRSpec channels that Pandeia marks as unobserved."""
    instrument = conf.get('instrument', {})
    if (
        extracted_noise is None
        or str(instrument.get('instrument', '')).lower() != 'nirspec'
    ):
        return None

    extracted_noise = np.asarray(extracted_noise, dtype=float)
    valid_noise = np.isfinite(extracted_noise) & (extracted_noise > 0.0)
    if full_saturation is None:
        return valid_noise

    full_saturation = np.asarray(full_saturation, dtype=float)
    if full_saturation.shape != extracted_noise.shape:
        return valid_noise

    fully_saturated = (
        np.isfinite(full_saturation)
        & (full_saturation > 0.0)
    )
    return valid_noise | fully_saturated


def dhs_f150w_wavelength_mask(conf, wave):
    """Build a wavelength mask for the low-throughput edge of DHS F150W.

    NIRCam DHS calculations can include F150W/F150W2 wavelength samples below
    the useful bandpass where the throughput is effectively zero. Those samples
    make the raw and binned diagnostic plots difficult to read, so PandExo
    trims them before downstream binning and plot packaging.

    Parameters
    ----------
    conf : dict
        Pandeia ``configuration`` dictionary for the calculation. The helper
        inspects ``conf["instrument"]`` to determine whether the setup is a
        NIRCam DHS calculation with a filter name starting with ``"f150w"``.
    wave : array-like
        Wavelength samples in microns.

    Returns
    -------
    numpy.ndarray or None
        Boolean mask selecting wavelengths greater than or equal to
        ``DHS_F150W_MIN_WAVELENGTH``. Returns ``None`` when the calculation is
        not a NIRCam DHS F150W/F150W2 setup, so callers can leave other modes
        unchanged.
    """
    instrument = conf.get('instrument', {})
    if str(instrument.get('instrument', '')).lower() != 'nircam':
        return None

    aperture = str(instrument.get('aperture', '')).lower()
    mode = str(instrument.get('mode', '')).lower()
    filter_name = str(instrument.get('filter', '')).lower()
    if 'dhs' not in aperture and mode != 'dhs':
        return None
    if not filter_name.startswith('f150w'):
        return None

    return np.asarray(wave, dtype=float) >= DHS_F150W_MIN_WAVELENGTH


def _no_valid_spectral_channels_message(conf, scalar):
    """Describe an all-invalid spectral extraction without requiring Pandeia."""
    instrument = conf.get('instrument', {})
    detector = conf.get('detector', {})
    details = [
        'All spectral channels have non-positive or non-finite '
        'Pandeia extracted_noise values.',
        'The target may be fully saturated, or the requested setup may '
        'not place any valid spectral channels on the detector.',
        'Try a fainter target, fewer groups, or a different subarray.',
    ]
    for key, value in [
        ('disperser', instrument.get('disperser')),
        ('filter', instrument.get('filter')),
        ('subarray', detector.get('subarray')),
        ('readout_pattern', detector.get('readout_pattern')),
        ('fraction_saturation', scalar.get('fraction_saturation')),
        ('sat_ngroups', scalar.get('sat_ngroups')),
    ]:
        if value is not None:
            details.append(f'{key}={value}')
    return ' '.join(details)


def _saturation_warning_message(conf, scalar, pandeia_warning, saturation_kind):
    """Describe a Pandeia saturation warning with setup context."""
    instrument = conf.get('instrument', {})
    detector = conf.get('detector', {})
    details = [
        f'Pandeia reports {saturation_kind} saturation for this observation.',
        str(pandeia_warning),
    ]
    for key, value in [
        ('instrument', instrument.get('instrument')),
        ('mode', instrument.get('mode')),
        ('disperser', instrument.get('disperser')),
        ('filter', instrument.get('filter')),
        ('subarray', detector.get('subarray')),
        ('readout_pattern', detector.get('readout_pattern')),
        ('fraction_saturation', scalar.get('fraction_saturation')),
        ('sat_ngroups', scalar.get('sat_ngroups')),
    ]:
        if value is not None:
            details.append(f'{key}={value}')
    return ' '.join(details)


def _pandeia_warning_is_active(value):
    """Return whether a Pandeia warning field contains a real warning."""
    if value is None or value is False:
        return False
    if isinstance(value, str):
        return value != 'All good'
    return True


def validate_saturation_state(conf, pand_dict, extracted_noise):
    """Warn on partial saturation and fail when no spectral channel is usable."""
    pandeia_warnings = pand_dict.get('warnings') or {}
    partial_warning = pandeia_warnings.get('partial_saturated')
    full_warning = pandeia_warnings.get('full_saturated')
    scalar = pand_dict.get('scalar') or {}

    if extracted_noise is not None:
        extracted_noise = np.asarray(extracted_noise, dtype=float)
        valid_noise = np.isfinite(extracted_noise) & (extracted_noise > 0.0)
        if not np.any(valid_noise):
            raise ValueError(
                _no_valid_spectral_channels_message(conf, scalar)
            )

    if _pandeia_warning_is_active(partial_warning):
        warnings.warn(
            _saturation_warning_message(
                conf, scalar, partial_warning, 'partial'
            ),
            UserWarning,
            stacklevel=2,
        )

    if _pandeia_warning_is_active(full_warning):
        warnings.warn(
            _saturation_warning_message(
                conf, scalar, full_warning, 'full'
            ),
            UserWarning,
            stacklevel=2,
        )


def is_phase_spec(calculation):
    """Return whether a selected calculation is a phase-curve calculation."""
    return calculation.startswith('phase_spec')


def compute_full_sim(dictinput,verbose=False): 
    """Top level function to set up exoplanet obs. for JW
    
    Function to set up explanet observations for JWST only and 
    compute simulated spectrum. It uses STScI's Pandeia to compute 
    instrument throughputs and WebbPSF to compute PSFs. 
    
    Parameters
    ----------
    dictinput : dict
        dictionary containing instrument parameters and exoplanet specific 
        parameters. {"pandeia_input":dict1, "pandexo_input":dict1}
    verbose : bool 
        (Optional) prints out check points throughout code 
    
    Returns
    -------
    dict
        large dictionary with 1d, 2d simualtions, timing info, instrument info, warnings
    
    Examples
    --------  
      
    >>> from .pandexo.engine.jwst import compute_full_sim 
    >>> from .pandexo.engine.justplotit import jwst_1d_spec
    >>> a = compute_full_sim({"pandeia_input": pandeiadict, "pandexo_input":exodict})
    >>> jwst_1d_spec(a)
    .. image:: 1d_spec.png
    
    Notes
    -----
    It is much easier to run simulations through either **run_online** or **justdoit**. **justdoit** contains functions to create input dictionaries and **run_online** contains web forms to create input dictionaries.
    
    See Also
    -------- 
        pandexo.engine.justdoit.run_pandexo : Best function for running pandexo runs
        pandexo.engine.run_online : Allows running functions through online interface
    """
    pandeia_input = dictinput['pandeia_input']
    pandexo_input = dictinput['pandexo_input']    	

    def log(message):
        if isinstance(verbose, str):
            print(f"[{verbose}] {message}", flush=True)
        elif verbose:
            print(message, flush=True)
	
    #which instrument 
    instrument = pandeia_input['configuration']['instrument']['instrument']
    conf = pandeia_input['configuration']
    validate_miri_lrs_subarray(conf)

    #now fix DHS #of spectra depending on the subarray
    is_dhs = 'dhs' in conf['instrument']['aperture']
    requested_dhs_readout = str(
        conf['detector'].get('readout_pattern', '')
    ).lower()
    optimize_dhs_readout = is_dhs and requested_dhs_readout == 'optimize'
    if is_dhs:
        if optimize_dhs_readout:
            # Pandeia does not recognize PandExo's optimization sentinel.
            conf['detector']['readout_pattern'] = DHS_READOUT_PATTERNS[0]
        elif requested_dhs_readout not in DHS_READOUT_PATTERNS:
            allowed = ', '.join(pattern.upper() for pattern in DHS_READOUT_PATTERNS)
            raise ValueError(
                f"Unsupported NIRCam DHS readout pattern: "
                f"{requested_dhs_readout}. Choose optimize or one of {allowed}."
            )
    if is_dhs:
        subarray = pandeia_input['configuration']['detector']['subarray']
        nspectra = int(subarray.split('-spectra')[0][-1])#2*int(substripe[substripe.find('stripe')+6])
        pandeia_input['configuration']['instrument']['aperture'] = f'dhs0spec{nspectra}'

        #if long wave setup with dhs is asked for change to lw grism 
        if (('32' in pandeia_input['configuration']['instrument']['filter']) or 
            ('44' in pandeia_input['configuration']['instrument']['filter'])): 
            pandeia_input['configuration']['instrument']['mode']='lw_tsgrism'
            pandeia_input['configuration']['instrument']['aperture']='lw' 
            pandeia_input['configuration']['instrument']['disperser']='grismr' 
    #if optimize is in the ngroups section, this will throw an error 
    #so create temp conf with 2 groups 
    if 'optimize' in str(conf['detector']['ngroup']): 
        conf_temp = deepcopy(conf) 
        if 'dhs' in conf['instrument']['aperture']: 
            #for DHS also need to swap to bright mode to get only the highest throughput spectra 
            conf_temp['instrument']['aperture'] = 'dhs0bright'
        conf_temp['detector']['ngroup'] = 2
    else: 
        conf_temp = conf

    i = _instrument_factory(config=conf_temp)
    
    #detector parameters
    det_pars = i.read_detector_pars()
    fullwell = det_pars.get('fullwell', det_pars.get('saturation_fullwell'))
    if fullwell is None:
        raise KeyError(
            "Detector parameters do not include 'fullwell' or "
            "'saturation_fullwell'."
        )
    rn = det_pars.get('rn', det_pars.get('readnoise'))
    if rn is None:
        raise KeyError(
            "Detector parameters do not include 'rn' or 'readnoise'."
        )
    mingroups = det_pars['mingroups']
        
    #exposure parameters 
    exp_pars = i.the_detector.exposure_spec
    tframe = exp_pars.tframe
    nframe = exp_pars.nframe
    nskip = exp_pars.ndrop2
    nsuperstripe = int(getattr(exp_pars, "nsuperstripe", 1) or 1)
    # Multistripe readouts use Pandeia's slope-derived MULTIACCUM noise.
    # SOSS exposes this as nsuperstripe > 1; DHS is identified separately
    # because its multi-spectrum readouts report nsuperstripe == 1.
    # Legacy readouts retain PandExo's historical first-minus-last method.
    calculation = select_calculation(
        pandexo_input['planet']['w_unit'], nsuperstripe, is_dhs=is_dhs
    )
    exposure_time_per_int = getattr(exp_pars, "exposure_time", None)
    exposure_time_inputs = {
        "tfffr": getattr(exp_pars, "tfffr", None),
        "nreset1": getattr(exp_pars, "nreset1", 1),
        "ndrop1": getattr(exp_pars, "ndrop1", 0),
        "ndrop3": getattr(exp_pars, "ndrop3", 0),
    }

    sat_unit = pandexo_input['observation']['sat_unit']

    if sat_unit =='%':
        sat_level = pandexo_input['observation']['sat_level']/100.0*fullwell
    elif sat_unit =='e':
        sat_level = pandexo_input['observation']['sat_level']
    else: 
        raise Exception("Saturation Level Needs Units: % fullwell or Electrons ")
    

    
    #parameteres needed from exo_input
    mag = pandexo_input['star']['mag']
    
    
    noccultations = pandexo_input['observation']['noccultations']
    R = pandexo_input['observation']['R']

    noise_floor = pandexo_input['observation']['noise_floor']

    
    #get stellar spectrum and in transit spec
    star_spec = create.outTrans(pandexo_input['star'])
    #get rstar if user calling from grid 
    both_spec = create.bothTrans(star_spec, pandexo_input['planet'], star=pandexo_input['star'])
    out_spectrum = np.array([both_spec['wave'], both_spec['flux_out_trans']])
    
    #get transit duration from phase curve or from input 
    if is_phase_spec(calculation):
        transit_duration = max(both_spec['time']) - min(both_spec['time'])
    else: 
        #convert to seconds, then remove quantity and convert back to float 
        transit_duration = float((pandexo_input['planet']['transit_duration']*u.Unit(pandexo_input['planet']['td_unit'])).to(u.second)/u.second)

    #amount of exposure time out-of-occultation, as a fraction of in-occ time 
    try:
        expfact_out = pandexo_input['observation']['fraction'] 
        log("WARNING: key input fraction has been replaced with new 'baseline option'. See notebook example")
        pandexo_input['observation']['baseline'] = pandexo_input['observation']['fraction'] 
        pandexo_input['observation']['baseline_unit'] ='frac'
    except:
        if pandexo_input['observation']['baseline_unit'] =='frac':
            expfact_out = pandexo_input['observation']['baseline'] 
        elif pandexo_input['observation']['baseline_unit'] =='total':
            expfact_out = transit_duration/( pandexo_input['observation']['baseline'] - transit_duration)
        elif pandexo_input['observation']['baseline_unit'] =='total_hrs':
            expfact_out = transit_duration/( pandexo_input['observation']['baseline']*3600.0 - transit_duration)
        else: 
            raise Exception("Wrong units for baseine: either 'frac' or 'total' or 'total_hrs' accepted")

    #add to pandeia input 
    pandeia_input['scene'][0]['spectrum']['sed']['spectrum'] = out_spectrum
    
    if isinstance(pandeia_input["configuration"]["detector"]["ngroup"], (float,int)):
        m = {"ngroup":int(pandeia_input["configuration"]["detector"]["ngroup"]), "tframe":tframe,
            "nframe":nframe,"mingroups":mingroups,"nskip":nskip,
            "nsuperstripe":nsuperstripe,
            "exposure_time_per_int":exposure_time_per_int,
            "exposure_time_ngroup":int(pandeia_input["configuration"]["detector"]["ngroup"])}
        m.update(exposure_time_inputs)
    else:
        #run pandeia once to determine max exposure time per int and get exposure params
        log("Optimization Reqested: Computing Duty Cycle")
        m = {"maxexptime_per_int":compute_maxexptime_per_int(pandeia_input, sat_level) , 
            "tframe":tframe,"nframe":nframe,"mingroups":mingroups,"nskip":nskip,
            "nsuperstripe":nsuperstripe}
        m.update(exposure_time_inputs)
        log("Finished Duty Cycle Calc")

    #calculate all timing info
    max_ngroup_instrument = max_ngroup[pandeia_input["configuration"]["instrument"]["instrument"]]
    timing, flags = compute_timing(m,transit_duration,expfact_out,noccultations,max_ngroup_instrument)

    if optimize_dhs_readout:
        selected = None
        rejected_readouts = []
        for readout_pattern in DHS_READOUT_PATTERNS:
            conf['detector']['readout_pattern'] = readout_pattern
            candidate_conf = deepcopy(conf)
            if 'dhs' in candidate_conf['instrument']['aperture']:
                candidate_conf['instrument']['aperture'] = 'dhs0bright'
            candidate_conf['detector']['ngroup'] = 2
            candidate_instrument = _instrument_factory(config=candidate_conf)
            candidate_detector = candidate_instrument.read_detector_pars()
            candidate_exposure = candidate_instrument.the_detector.exposure_spec
            candidate_m = {
                "tframe": candidate_exposure.tframe,
                "nframe": candidate_exposure.nframe,
                "mingroups": candidate_detector['mingroups'],
                "nskip": candidate_exposure.ndrop2,
                "nsuperstripe": int(
                    getattr(candidate_exposure, "nsuperstripe", 1) or 1
                ),
                "tfffr": getattr(candidate_exposure, "tfffr", None),
                "nreset1": getattr(candidate_exposure, "nreset1", 1),
                "ndrop1": getattr(candidate_exposure, "ndrop1", 0),
                "ndrop3": getattr(candidate_exposure, "ndrop3", 0),
            }
            if "maxexptime_per_int" in m:
                candidate_m["maxexptime_per_int"] = m["maxexptime_per_int"]
            else:
                candidate_m["ngroup"] = m["ngroup"]
            candidate_timing, candidate_flags = compute_timing(
                candidate_m,
                transit_duration,
                expfact_out,
                noccultations,
                max_ngroup_instrument,
            )
            allocation_overhead = nircam_dhs_no_ta_overhead(
                candidate_exposure.tframe
            )
            excess_rate, data_excess = estimate_dhs_data_excess(
                conf['detector']['subarray'],
                readout_pattern,
                candidate_timing['APT: Num Groups per Integration'],
                candidate_timing['Transit+Baseline, no overhead (hrs)'],
                allocation_overhead_seconds=allocation_overhead,
            )
            candidate_timing[
                'Estimated DHS Data Excess Rate (GB/hr)'
            ] = excess_rate
            candidate_timing['Estimated DHS Data Excess (GB)'] = data_excess
            candidate_timing[
                'Assumed DHS Allocation Overhead (sec)'
            ] = allocation_overhead
            selected = (
                readout_pattern, candidate_instrument, candidate_timing,
                candidate_flags
            )

            saturation_limited = (
                candidate_flags['flag_default'].startswith(
                    'Optimized NGROUPS below minimum'
                )
            )
            if saturation_limited or (
                    data_excess <= DHS_DATA_EXCESS_RECOMMENDED_LIMIT_GB):
                break
            rejected_readouts.append(readout_pattern.upper())

        readout_pattern, i, timing, flags = selected
        conf['detector']['readout_pattern'] = readout_pattern
        exp_pars = i.the_detector.exposure_spec
        tframe = exp_pars.tframe
        nframe = exp_pars.nframe
        nskip = exp_pars.ndrop2
        nsuperstripe = int(getattr(exp_pars, "nsuperstripe", 1) or 1)
        flags['flag_dhs_readout'] = (
            f"Selected {readout_pattern.upper()}"
            + (
                f" after {', '.join(rejected_readouts)} exceeded the "
                f"{DHS_DATA_EXCESS_RECOMMENDED_LIMIT_GB:g} GB recommendation."
                if rejected_readouts else "."
            )
            + " Estimate assumes no target acquisition and a standard "
            "2,100-second initial slew."
        )
    elif is_dhs:
        allocation_overhead = nircam_dhs_no_ta_overhead(tframe)
        excess_rate, data_excess = estimate_dhs_data_excess(
            conf['detector']['subarray'],
            requested_dhs_readout,
            timing['APT: Num Groups per Integration'],
            timing['Transit+Baseline, no overhead (hrs)'],
            allocation_overhead_seconds=allocation_overhead,
        )
        timing['Estimated DHS Data Excess Rate (GB/hr)'] = excess_rate
        timing['Estimated DHS Data Excess (GB)'] = data_excess
        timing['Assumed DHS Allocation Overhead (sec)'] = allocation_overhead
        flags['flag_dhs_readout'] = (
            f"User selected {requested_dhs_readout.upper()}; readout pattern "
            "optimization was not performed. Estimate assumes no target "
            "acquisition and a standard 2,100-second initial slew."
        )
    
    #Simulate out trans and in transit
    log("Starting Out of Transit Simulation")
    out = perform_out(pandeia_input, pandexo_input,timing, both_spec)

    #extract extraction area before dict conversion
    extraction_area = out.extraction_area
    out = out.as_dict()
    out.pop('3d')
    update_timing_measurement_time(timing, out['scalar']['measurement_time'])
    log("End out of Transit")

    #Remove effects of Quantum Yield from shot noise 
    out = remove_QY(out, instrument)

    #this kind of redundant going to compute inn from out instead 
    #keep perform_in but change inputs to (out, timing, both_spec)
    log("Starting In Transit Simulation")
    inn = perform_in(pandeia_input, pandexo_input,timing, both_spec, out, calculation)
    log("End In Transit")
    


    #compute warning flags for timing info 
    warnings = add_warnings(out, timing, sat_level/fullwell, flags, instrument) 

    compNoise = ExtractSpec(inn, out, rn, extraction_area, timing)
    
    #slope method is pandeia's pure noise calculation (taken from SNR)
    #contains correlated noise, RN, dark current, sky, 
    #uses MULTIACCUM formula so we deviated from this. 
    #could eventually come back to this if Pandeia adopts First-Last formula
    if calculation == 'slope method': 
        #Extract relevant info from pandeia output (1d curves and wavelength) 
        #extracted flux in units of electron/s
        w = out['1d']['extracted_flux'][0]
        result = compNoise.run_slope_method()

    #derives noise from 2d postage stamps. Doing this results in a higher 
    #1d flux rate than the Pandeia gets from extracting its own. 
    #this should be used to benchmark Pandeia's 1d extraction  
    elif calculation == '2d extract':
        w = out['1d']['extracted_flux'][0]
        result = compNoise.run_2d_extract()
    
    #this is the historical noise calculation used for legacy readouts. It derives
    #its own calculation of readnoise and does not use MULTIACUMM 
    #noise formula  
    elif calculation == 'fml':
        w = out['1d']['extracted_flux'][0]
        result = compNoise.run_f_minus_l()
    
    elif calculation == 'phase_spec_fml':
        result = compNoise.run_phase_spec_fml()
        w = result['time']
    elif calculation == 'phase_spec_slope':
        result = compNoise.run_phase_spec_slope()
        w = result['time']
    else:
        result = None
        raise Exception('WARNING: Calculation method not found.')
        
    varin = result['var_in_1d']
    varout = result['var_out_1d']
    extracted_flux_out = result['photon_out_1d']
    extracted_flux_inn = result['photon_in_1d']
    extracted_flux_per_int_out = result.get('photon_out_1d_per_int')
    pandeia_extracted_noise = None
    pandeia_snr_int = None
    pandeia_full_saturation = None
    if not is_phase_spec(calculation):
        pandeia_extracted_noise = np.asarray(
            out['1d']['extracted_noise'][1], dtype=float
        )
        pandeia_full_saturation = _pandeia_1d_values_at_wave(
            out, 'n_full_saturated', w
        )
        pandeia_snr_int = [
            np.asarray(out['1d']['sn'][0], dtype=float),
            np.asarray(out['1d']['sn'][1], dtype=float),
        ]

    input_wave_order = np.argsort(w, kind='mergesort')
    if not np.array_equal(input_wave_order, np.arange(len(w))):
        w = w[input_wave_order]
        varin = varin[input_wave_order]
        varout = varout[input_wave_order]
        extracted_flux_out = extracted_flux_out[input_wave_order]
        extracted_flux_inn = extracted_flux_inn[input_wave_order]
        if extracted_flux_per_int_out is not None:
            extracted_flux_per_int_out = extracted_flux_per_int_out[input_wave_order]
        result['rn[out,in]'] = sort_by_wave_order(result['rn[out,in]'], input_wave_order)
        result['bkg[out,in]'] = sort_by_wave_order(result['bkg[out,in]'], input_wave_order)
        if pandeia_extracted_noise is not None:
            pandeia_extracted_noise = pandeia_extracted_noise[input_wave_order]
        if pandeia_full_saturation is not None:
            pandeia_full_saturation = pandeia_full_saturation[input_wave_order]
        if pandeia_snr_int is not None:
            pandeia_snr_int = sort_by_wave_order(pandeia_snr_int, input_wave_order)

    validate_saturation_state(conf, out, pandeia_extracted_noise)

    valid_channel = nirspec_valid_channel_mask(
        conf, pandeia_extracted_noise, pandeia_full_saturation
    )
    if valid_channel is not None:
        w = w[valid_channel]
        varin = varin[valid_channel]
        varout = varout[valid_channel]
        extracted_flux_out = extracted_flux_out[valid_channel]
        extracted_flux_inn = extracted_flux_inn[valid_channel]
        if extracted_flux_per_int_out is not None:
            extracted_flux_per_int_out = extracted_flux_per_int_out[valid_channel]
        pandeia_full_saturation = pandeia_full_saturation[valid_channel]
        result['rn[out,in]'] = sort_by_wave_order(result['rn[out,in]'], valid_channel)
        result['bkg[out,in]'] = sort_by_wave_order(result['bkg[out,in]'], valid_channel)
        pandeia_snr_int = sort_by_wave_order(pandeia_snr_int, valid_channel)

    dhs_wavelength_channel = None
    if not is_phase_spec(calculation):
        dhs_wavelength_channel = dhs_f150w_wavelength_mask(conf, w)
    if dhs_wavelength_channel is not None:
        w = w[dhs_wavelength_channel]
        varin = varin[dhs_wavelength_channel]
        varout = varout[dhs_wavelength_channel]
        extracted_flux_out = extracted_flux_out[dhs_wavelength_channel]
        extracted_flux_inn = extracted_flux_inn[dhs_wavelength_channel]
        if extracted_flux_per_int_out is not None:
            extracted_flux_per_int_out = extracted_flux_per_int_out[dhs_wavelength_channel]
        if pandeia_extracted_noise is not None:
            pandeia_extracted_noise = pandeia_extracted_noise[dhs_wavelength_channel]
        if pandeia_full_saturation is not None:
            pandeia_full_saturation = pandeia_full_saturation[dhs_wavelength_channel]
        result['rn[out,in]'] = sort_by_wave_order(result['rn[out,in]'], dhs_wavelength_channel)
        result['bkg[out,in]'] = sort_by_wave_order(result['bkg[out,in]'], dhs_wavelength_channel)
        pandeia_snr_int = sort_by_wave_order(pandeia_snr_int, dhs_wavelength_channel)

        
    #bin the data according to user input 
    if R != None: 
        wbin = bin_wave_to_R(w, R)
        saturated_channel = (
            np.zeros(len(w), dtype=float)
            if pandeia_full_saturation is None
            else (pandeia_full_saturation > 0.0).astype(float)
        )

        photon_out_bin = uniform_tophat_sum(wbin, w,extracted_flux_out)
        photon_in_bin = uniform_tophat_sum(wbin,w, extracted_flux_inn)
        if extracted_flux_per_int_out is None:
            electron_per_int_bin = photon_out_bin / result.get(
                'real_nint_out', result.get('nint_out', 1)
            )
        else:
            electron_per_int_bin = uniform_tophat_sum(
                wbin, w, extracted_flux_per_int_out
            )
        var_in_bin = uniform_tophat_sum(wbin, w,varin)
        var_out_bin = uniform_tophat_sum(wbin,w, varout)
        full_saturation_bin = (
            uniform_tophat_sum(wbin, w, saturated_channel) > 0.0
        )

        valid_photon = photon_out_bin > 0
        wbin = wbin[valid_photon]
        photon_in_bin = photon_in_bin[valid_photon]
        electron_per_int_bin = electron_per_int_bin[valid_photon]
        var_in_bin = var_in_bin[valid_photon]
        var_out_bin = var_out_bin[valid_photon]
        full_saturation_bin = full_saturation_bin[valid_photon]
        photon_out_bin = photon_out_bin[valid_photon]
    else: 
        wbin = w
        photon_out_bin = extracted_flux_out
        if extracted_flux_per_int_out is None:
            electron_per_int_bin = photon_out_bin / result.get(
                'real_nint_out', result.get('nint_out', 1)
            )
        else:
            electron_per_int_bin = extracted_flux_per_int_out
        full_saturation_bin = (
            np.zeros(len(wbin), dtype=bool)
            if pandeia_full_saturation is None
            else pandeia_full_saturation > 0.0
        )
        valid_photon = photon_out_bin > 0
        wbin = wbin[valid_photon]
        photon_in_bin = extracted_flux_inn
        photon_in_bin = photon_in_bin[valid_photon]
        electron_per_int_bin = electron_per_int_bin[valid_photon]
        var_in_bin = varin
        var_in_bin = var_in_bin[valid_photon]
        var_out_bin = varout
        var_out_bin = var_out_bin[valid_photon]
        full_saturation_bin = full_saturation_bin[valid_photon]
        photon_out_bin = photon_out_bin[valid_photon]
        
    
    if is_phase_spec(calculation):
        to = timing["Measurement Time per Integration (sec)"]
        ti = timing["Measurement Time per Integration (sec)"]
        nint_in = 1
        nint_out = 1
    else: 
        #otherwise error propagation and account for different 
        #times in transit and out 
        to = result['on_source_out']
        ti = result['on_source_in']
        nint_in = result['nint_in']
        nint_out = result['nint_out']
    var_tot = (to/ti/photon_out_bin)**2.0 * var_in_bin + (photon_in_bin*to/ti/photon_out_bin**2.0)**2.0 * var_out_bin
    error_spec = np.sqrt(var_tot)
        
    #factor in occultations to noise 
    nocc = timing['Number of Transits']
    error_spec = error_spec / np.sqrt(nocc)
        
    #Add in user specified noise floor 
    error_spec_nfloor = add_noise_floor(noise_floor, wbin, error_spec) 

    #add in random noise for the simulated spectrum 
    np.random.seed()
    rand_noise= error_spec_nfloor*(np.random.randn(len(wbin)))
    raw_spec = (photon_out_bin/to-photon_in_bin/ti)/(photon_out_bin/to)
    sim_spec = raw_spec + rand_noise 
    
    #if secondary tranist, multiply spectra by -1 
    if pandexo_input['planet']['f_unit'] == 'fp/f*':
        sim_spec = -1.0*sim_spec
        raw_spec = -1.0*raw_spec    

    wave_order = np.argsort(wbin, kind='mergesort')
    if not np.array_equal(wave_order, np.arange(len(wbin))):
        wbin = wbin[wave_order]
        photon_out_bin = photon_out_bin[wave_order]
        photon_in_bin = photon_in_bin[wave_order]
        electron_per_int_bin = electron_per_int_bin[wave_order]
        var_in_bin = var_in_bin[wave_order]
        var_out_bin = var_out_bin[wave_order]
        error_spec = error_spec[wave_order]
        error_spec_nfloor = error_spec_nfloor[wave_order]
        raw_spec = raw_spec[wave_order]
        sim_spec = sim_spec[wave_order]
        full_saturation_bin = full_saturation_bin[wave_order]
   
    #package processed data
    finalspec = {'wave':wbin,
              'spectrum': raw_spec,
              'spectrum_w_rand':sim_spec,
              'error_w_floor':error_spec_nfloor}
    finalspec = mask_fully_saturated_final_spectrum(
        finalspec, full_saturation_bin
    )
    
    rawstuff = {
                'electrons_out':photon_out_bin*nocc, 
                'electrons_in':photon_in_bin*nocc,
                'electron_per_int':electron_per_int_bin,
                'snr_int': (
                    pandeia_snr_int
                    if pandeia_snr_int is not None
                    else [out['1d']['sn'][0], out['1d']['sn'][1]]
                ),
                'var_in':var_in_bin*nocc, 
                'var_out':var_out_bin*nocc,
                'e_rate_out':photon_out_bin/to,
                'e_rate_in':photon_in_bin/ti,
                'wave':wbin,
                'error_no_floor':error_spec, 
                'rn[out,in]':sort_by_wave_order(result['rn[out,in]'], wave_order),
                'bkg[out,in]':sort_by_wave_order(result['bkg[out,in]'], wave_order)
                }
 
    result_dict = as_dict(out,both_spec ,finalspec, 
                timing, mag, sat_level, warnings,
                pandexo_input['planet']['f_unit'], rawstuff,calculation)

    return result_dict 
    
def compute_maxexptime_per_int(pandeia_input, sat_level):
    """Computes optimal maximum exposure time per integration
    
    Function to simulate 2d jwst image with 2 groups, 1 integration, 1 exposure 
    and return the maximum time 
    for one integration before saturation occurs. If saturation has 
    already occured, returns maxexptime_per_int as np.nan. This then 
    tells Pandexo to set min number of groups (ngroups =2). This avoids 
    error if saturation occurs. This routine assumes that min ngroups is 2. 
    
    Parameters
    ----------
    pandeia_input : dict 
        pandeia dictionary input 
    sat_level : int or float
        user defined saturation level in units of electrons
    
    Returns
    ------- 
    float 
        Maximum exposure time per integration before specified saturation level
    
    Examples
    --------
    
    >>> max_time = compute_maxexptime_per_int(pandeia_input, 50000.0)
    >>> print(max_time)
    12.0
    """
    
    #run once to get 2d rate image 
    pandeia_input['configuration']['detector']['ngroup'] = int(2 )
    pandeia_input['configuration']['detector']['nint'] = 1 
    pandeia_input['configuration']['detector']['nexp'] = 1
    
    report = _perform_calculation(pandeia_input, dict_report=False)
    report_dict = report.as_dict() 

    # count rate on the detector in e-/second/pixel
    #det = report_dict['2d']['detector']
    det = report.signal.rate_plus_bg_list[0]['fp_pix']

    timeinfo = report_dict['information']['exposure_specification']
    #totaltime = timeinfo['tgroup']*timeinfo['ngroup']*timeinfo['nint']
    
    maxdetvalue = np.max(det)


    #maximum time before saturation per integration 
    #based on user specified saturation level

    try:
        maxexptime_per_int = sat_level/maxdetvalue
    except: 
        maxexptime_per_int = np.nan
    
    return maxexptime_per_int
        
def compute_timing(m,transit_duration,expfact_out,noccultations,max_ngroup_instrument=65536.0):
    """Computes all timing info for observation
    
    Computes all JWST specific timing info for observation including. Some pertinent 
    JWST terminology:

        - frame: The result of sequentially clocking and digitizing all pixels in a rectangular area of an SCA. **Full-fame readout** means to digitize all pixels in an SCA, including reference pixels. **Frame** also applies to the result of clocking and digitizing a subarray on an SCA.
        - group: One or more consecutively read frames. There are no intervening resets. Frames may be averaged to form a group but for exoplanets the read out scheme is always 1 frame = 1 group
        - integration: The end result of resetting the detector and then non-destructively sampling it one or more times over a finite period of time before resetting the detector again. This is a unit of data for which signal is proportional to intensity, and it consists of one or more GROUPS.
        - exposure: The end result of one or more INTEGRATIONS over a finite period of time.  EXPOSURE defines the contents of a single FITS file.
    
    Parameters
    ---------
    m : dict 
        Dictionary output from **compute_maxexptime_per_int**
    transit_duration : float or int 
        transit duration in seconds 
    expfact_out : float or int 
        fraction of time spent in transit versus out of transit 
    noccultations : int 
        number of transits 
    
    Returns
    ------- 
    timing : dict
        All timing info
    warningflag : dict 
        Warning flags
    
    Examples
    --------
    >>> timing, flags = compute_timing(m, 2*60.0*60.0, 1.0, 1.0)
    >>> print((list(timing.keys())))
    ['Number of Transits', 'Num Integrations Out of Transit', 'Num Integrations In Transit', 
    'APT: Num Groups per Integration', 'Seconds per Frame', 'Observing Efficiency (%)', 'On Source Time(sec)', 
    'Exposure Time Per Integration (secs)', 'Reset time Plus 30 min TA time (hrs)',
    'APT: Num Integrations per Occultation', 'Transit Duration']
    """
    tframe = m['tframe']
    nframe = m['nframe']
    nskip = m['nskip']
    mingroups = m['mingroups']
    nsuperstripe = int(m.get("nsuperstripe", 1) or 1)
    if nsuperstripe < 1:
        nsuperstripe = 1
    def _clocktime_per_int(ngroups):
        if (m.get("exposure_time_per_int") is not None and
                ngroups == m.get("exposure_time_ngroup")):
            return m["exposure_time_per_int"]
        if nsuperstripe > 1 and m.get("tfffr") is not None:
            full_cycle = nsuperstripe * (
                m["tfffr"]
                + tframe * (
                    m.get("nreset1", 1)
                    + m.get("ndrop1", 0)
                    + (ngroups - 1.0)*(nframe + nskip)
                    + nframe
                    + m.get("ndrop3", 0)
                )
            )
            return full_cycle
        return (
            1.0 + nframe + (ngroups - 1.0) * (nframe + nskip)
        ) * tframe
    def _timing_values(ngroups):
        if ngroups == 1:
            frame_zero_dead = 0
        else:
            frame_zero_dead = -1

        science_time_per_int = (ngroups + frame_zero_dead)*tframe*(nframe+nskip)
        measurement_time_per_int = science_time_per_int*nsuperstripe
        clocktime_per_int = _clocktime_per_int(ngroups)
        eff = measurement_time_per_int/float(nsuperstripe)/clocktime_per_int
        nint_per_occultation = transit_duration/clocktime_per_int
        nint_in = np.ceil(nint_per_occultation)
        nint_out = np.ceil(nint_in/expfact_out)

        return {
            "frame_zero_dead": frame_zero_dead,
            "science_time_per_int": science_time_per_int,
            "measurement_time_per_int": measurement_time_per_int,
            "clocktime_per_int": clocktime_per_int,
            "eff": eff,
            "nint_per_occultation": nint_per_occultation,
            "nint_in": nint_in,
            "nint_out": nint_out,
            "exptime_per_int": ngroups*tframe,
        }

    optimized_ngroups = "maxexptime_per_int" in m
    if optimized_ngroups:
        #are we starting with a exposure time ?
        maxexptime_per_int = m['maxexptime_per_int']
    else:
        #or a pre defined number of groups specified by user
        ngroups_per_int = m['ngroup']
        
    flag_default = "All good"
    flag_high = "All good"
    flag_min_nint = "All good"
    if optimized_ngroups:
        #Frist, if maxexptime_per_int has been defined (from above), compute ngroups_per_int
        
        #number of frames in one integration is the maximum time beofre exposure 
        #divided by the time it takes for one frame. Note this does not include 
        #reset frames 

        nframes_per_int = np.floor(maxexptime_per_int/tframe)
    
        #for exoplanets nframe =1 an nskip always = 0 so ngroups_per_int 
        #and nframes_per_int area always the same 
        ngroups_per_int = np.floor(nframes_per_int/(nframe + nskip)) 
    
        #put restriction on number of groups 
        #there is a hard limit to the maximum number groups. 
        #if you exceed that limit, set it to the maximum value instead.
        #also set another check for saturation

        if ngroups_per_int > max_ngroup_instrument:
            ngroups_per_int = max_ngroup_instrument
            flag_high = f"Optimized NGROUPS above maximum ({max_ngroup_instrument}). SET TO NGROUPS={max_ngroup_instrument}"
 
        if (ngroups_per_int < mingroups) | np.isnan(ngroups_per_int):
            ngroups_per_int = mingroups
            nframes_per_int = mingroups
            flag_default = f"Optimized NGROUPS below minimum ({mingroups}). SET TO NGROUPS={mingroups}"

    elif 'ngroups_per_int' in locals(): 
        #if it maxexptime_per_int been defined then set nframes per int 
        nframes_per_int = ngroups_per_int*(nframe+nskip)
        
        #if that didn't work its because maxexptime_per_int is nan .. run calc with mingroups
    else:
        #if maxexptime_per_int is nan then just ngroups and nframe to 2 
        #for the sake of not returning error
        ngroups_per_int = mingroups
        nframes_per_int = mingroups
        flag_default = f"Something went wrong. SET TO NGROUPS={mingroups}"

    timing_values = _timing_values(ngroups_per_int)
    frame_zero_dead = timing_values["frame_zero_dead"]
    science_time_per_int = timing_values["science_time_per_int"]
    measurement_time_per_int = timing_values["measurement_time_per_int"]
    exptime_per_int = timing_values["exptime_per_int"]
    clocktime_per_int = timing_values["clocktime_per_int"]
    eff = timing_values["eff"]
    nint_per_occultation = timing_values["nint_per_occultation"]
    nint_in = timing_values["nint_in"]
    nint_out = timing_values["nint_out"]
    
    #you would never want a single integration in transit. 
    #here we assume that for very dim things, you would want at least 
    #3 integrations in transit 
    if nint_in < min_nint_trans:
        original_ngroups_per_int = ngroups_per_int
        original_nint_in = nint_in

        if optimized_ngroups:
            ngroups_per_int = np.max(
                [mingroups, np.floor(ngroups_per_int/min_nint_trans)]
            )
            timing_values = _timing_values(ngroups_per_int)
            while (
                    timing_values["nint_in"] < min_nint_trans and
                    ngroups_per_int > mingroups):
                ngroups_per_int -= 1
                timing_values = _timing_values(ngroups_per_int)

            frame_zero_dead = timing_values["frame_zero_dead"]
            science_time_per_int = timing_values["science_time_per_int"]
            measurement_time_per_int = timing_values["measurement_time_per_int"]
            exptime_per_int = timing_values["exptime_per_int"]
            clocktime_per_int = timing_values["clocktime_per_int"]
            eff = timing_values["eff"]
            nint_per_occultation = timing_values["nint_per_occultation"]
            nint_in = timing_values["nint_in"]
            nint_out = timing_values["nint_out"]
            flag_min_nint = (
                f"Optimized NGROUPS would produce {int(original_nint_in)} "
                f"in-transit integrations. Reduced NGROUPS from "
                f"{int(original_ngroups_per_int)} to {int(ngroups_per_int)} "
                f"to require at least {min_nint_trans} in-transit integrations."
            )
            if nint_in < min_nint_trans:
                flag_min_nint = (
                    f"Optimized NGROUPS would produce {int(original_nint_in)} "
                    f"in-transit integrations. Reduced NGROUPS from "
                    f"{int(original_ngroups_per_int)} to the minimum allowed "
                    f"value of {int(ngroups_per_int)}, but this still produces "
                    f"fewer than {min_nint_trans} in-transit integrations."
                )
        else:
            flag_min_nint = (
                f"User-specified NGROUPS produces {int(nint_in)} in-transit "
                f"integrations, below the recommended minimum of "
                f"{min_nint_trans}. NGROUPS was not changed."
            )
        
    if nint_out < min_nint_trans:
        nint_out = min_nint_trans

    effective_nint_in = nint_in / float(nsuperstripe)
    effective_nint_out = nint_out / float(nsuperstripe)
    # measurement_time_per_int is the Pandeia value for one nint=1
    # calculation. In SOSS multistripe modes it includes all superstripes,
    # so the per-wavelength on-source time must use effective_nint_*.
    on_source_in = effective_nint_in * measurement_time_per_int
    on_source_out = effective_nint_out * measurement_time_per_int
    effective_on_source_in = effective_nint_in * measurement_time_per_int
    effective_on_source_out = effective_nint_out * measurement_time_per_int
   
    timing = {
        "Transit Duration" : (transit_duration)/60.0/60.0,
        "Seconds per Frame" : tframe,
        "Time/Integration incl reset (sec)":clocktime_per_int,
        "Measurement Time per Integration (sec)":measurement_time_per_int,
        "APT: Num Groups per Integration" :int(ngroups_per_int), 
        "Num Integrations Out of Transit":int(nint_out),
        "Num Integrations In Transit":int(nint_in),
        "Num Superstripes":nsuperstripe,
        "Effective Integrations Out of Transit":effective_nint_out,
        "Effective Integrations In Transit":effective_nint_in,
        "On Source Time Out of Transit":on_source_out,
        "On Source Time In Transit":on_source_in,
        "Effective On Source Time Out of Transit":effective_on_source_out,
        "Effective On Source Time In Transit":effective_on_source_in,
        "APT: Num Integrations per Occultation":int(nint_out+nint_in),
        "Observing Efficiency (%)": eff*100.0,
        "Transit+Baseline, no overhead (hrs)": (nint_out+nint_in)*clocktime_per_int/60.0/60.0, 
        "Number of Transits": noccultations,
        "Zero Frame Efficiency Loss":frame_zero_dead
        }      

    return timing, {
        'flag_default':flag_default,
        'flag_high':flag_high,
        'flag_min_nint':flag_min_nint,
    }

def update_timing_measurement_time(timing, measurement_time_per_int):
    """Sync timing metadata to the measurement_time reported by Pandeia.

    PandExo runs Pandeia with nint=1 and scales integrations internally. For
    SOSS multistripe modes, Pandeia's one-integration measurement_time includes
    all superstripes, so per-wavelength time uses the effective integration
    counts in the timing dictionary.
    """
    nsuperstripe = float(timing.get("Num Superstripes", 1) or 1)
    nint_in = timing["Num Integrations In Transit"]
    nint_out = timing["Num Integrations Out of Transit"]
    effective_nint_in = nint_in / nsuperstripe
    effective_nint_out = nint_out / nsuperstripe

    timing["Measurement Time per Integration (sec)"] = measurement_time_per_int
    timing["Effective Integrations In Transit"] = effective_nint_in
    timing["Effective Integrations Out of Transit"] = effective_nint_out
    timing["On Source Time In Transit"] = effective_nint_in * measurement_time_per_int
    timing["On Source Time Out of Transit"] = effective_nint_out * measurement_time_per_int
    timing["Effective On Source Time In Transit"] = effective_nint_in * measurement_time_per_int
    timing["Effective On Source Time Out of Transit"] = effective_nint_out * measurement_time_per_int

def remove_QY(pandeia_dict, instrument):
    """Removes Quantum Yield from Pandeia Fluxes. Place Holder. 
    
    Parameters
    ----------
    pandeia_dict : dict 
        pandeia output dictionary
    instrument : str
        instrument running
    
    Returns
    -------
    dict 
        same exact dictionary with extracted_flux = extracted_flux/QY
    """
    if instrument == 'niriss':
        try:
            qy = fits.open(os.path.join(default_refdata_directory,'jwst', instrument,'qe' ,'jwst_niriss_h2rg_qe_20221003172003.fits'))
        except: 
            raise Exception('PANDEIA REFERENCE DATA NEEDS TO BE UPDATED')

        x_grid = pandeia_dict['1d']['extracted_flux'][0]
        qy_on_grid = np.interp(x_grid, qy[1].data['WAVELENGTH'], qy[1].data['CONVERSION'])
    elif instrument == 'nirspec':
        try:
            qy = fits.open(os.path.join(default_refdata_directory,'jwst', instrument,'qe' ,'jwst_nirspec_qe_20160902193401.fits'))
        except: 
            raise Exception('PANDEIA REFERENCE DATA NEEDS TO BE UPDATED')        
        x_grid = pandeia_dict['1d']['extracted_flux'][0]
        qy_on_grid = np.interp(x_grid, qy[1].data['WAVELENGTH'], qy[1].data['CONVERSION'])
    else:
        #nircam and miri currently have no qy effects
        qy_on_grid = 1.0

    pandeia_dict['1d']['extracted_flux'][1] = pandeia_dict['1d']['extracted_flux'][1]/qy_on_grid
    return pandeia_dict

def perform_out(pandeia_input, pandexo_input,timing, both_spec):
    """Runs pandeia for the out of transit data
    
    Parameters
    ----------
    pandeia_input : dict 
        pandeia specific input info 
    pandexo_input : dict 
        exoplanet specific observation info 
    timing : dict 
        timing dictionary from **compute_timing** 
    both_spec : dict 
        dictionary transit spectra computed from **createInput.bothTrans** 

    Returns
    -------
    dict 
        pandeia output dictionary for out of transit data 
    """
    #pandeia inputs, simulate one integration at a time 
    pandeia_input['configuration']['detector']['ngroup'] = int(timing['APT: Num Groups per Integration'])
    pandeia_input['configuration']['detector']['nint'] = 1#int(timing['Num Integrations Out of Transit'])
    pandeia_input['configuration']['detector']['nexp'] = 1 

    report_out = _perform_calculation(pandeia_input, dict_report=False)
    
    return report_out

    
def perform_in(pandeia_input, pandexo_input,timing, both_spec, out, calculation): 
    """Computes in transit data 
    
    Runs Pandeia for the in transit data or computes the in transit simulation 
    from the out of transit pandeia run 
    
    Parameters
    ----------
    pandeia_input : dict 
        pandeia specific input info 
    pandexo_input : dict 
        exoplanet specific observation info 
    timing : dict 
        timing dictionary from **compute_timing** 
    both_spec : dict 
        dictionary transit spectra computed from **createInput.bothTrans** 
    out : dict 
        out of transit dictionary from **perform_in**
    calculation : str
        key which speficies the kind of noise calcualtion 
        (2d extract, slope method, fml, phase_spec_fml, phase_spec_slope).
        Selected automatically from the detector readout mode.

    Returns
    -------
    dict
        pandeia output dictionary 
    """
    
    #function to run pandeia for in transit
    if is_phase_spec(calculation):
        #return the phase curve since it's all we need 
        report_in = {'time': both_spec['time'],'planet_phase': both_spec['planet_phase']}
    elif calculation == 'fml':
        #for FML method, we only use the flux rate calculated in pandeia so 
        #can compute in transit flux rate without running pandeia a third time     
        report_in = deepcopy(out)
        
        transit_depth = np.interp(report_in['1d']['extracted_flux'][0],
                                    both_spec['wave'], both_spec['frac'])
        report_in['1d']['extracted_flux'][1] = report_in['1d']['extracted_flux'][1]*transit_depth
    else: 
        #only run pandeia a third time if doing slope method and need accurate run for the 
        #nint and timing
        pandeia_input['configuration']['detector']['ngroup'] = int(timing['APT: Num Groups per Integration'])
        pandeia_input['configuration']['detector']['nint'] = 1#int(timing['Num Integrations In Transit'])
        pandeia_input['configuration']['detector']['nexp'] = 1
  
        in_transit_spec = np.array([both_spec['wave'], both_spec['flux_in_trans']])
    
        pandeia_input['scene'][0]['spectrum']['sed']['spectrum'] = in_transit_spec

        report_in = _perform_calculation(pandeia_input, dict_report=True)
        instrument = pandeia_input['configuration']['instrument']['instrument']
        #remove QY effects 
        report_in = remove_QY(report_in, instrument)
        report_in.pop('3d')
    
    return report_in
          
def add_warnings(pand_dict, timing, sat_level, flags,instrument): 
    """Add warnings for front end 
    
    Adds in necessary warning flags for a JWST observation usually associated with 
    too few or too many groups or saturation. Alerts user if saturation level is higher 
    than 80 percent and if the number of groups is less than 5. Or, if the full well is 
    greater than 80. These warnings are currently very arbitrary. Will be updated as 
    better JWST recommendations are made. 
    
    Parameters
    ----------
    pand_dict : 
        output from pandeia run 
    timing : dict 
        output from **compute_timing** 
    sat_level : int or float 
        user specified saturation level in fractional (00/100)
    flags : dict 
        warning flags taken from output of **compute_timing**
    instrument : str 
        Only allowable strings are: "nirspec", "niriss", "nircam", "miri"
    
    Returns
    -------
    dict 
        all warnings 

    Notes
    -----    
    These are warnings are just suggestions and are not yet required.   
    """

    ngroups_per_int = timing['APT: Num Groups per Integration']
  
    #check for saturation 
    try:  
        flag_nonl = pand_dict['warnings']['partial_saturated']
    except: 
        flag_nonl = "All good"    
    try: 
        flag_sat = pand_dict['warnings']['full_saturated']
    except: 
        flag_sat = "All good"
        
    #check for too small number of groups
    flag_low = "All good"
    flag_perc = "All good"

    if (sat_level > .80) & (ngroups_per_int <3):
        flag_low = "% full well>80% & only " + str(ngroups_per_int) + " groups"
    if ngroups_per_int==1:
        flag_low+='. Ngroups=1 is a new mode since Cycle 4 and has not been rigorously tested. Proceed with caution.'
    if (sat_level > .80): 
        flag_perc = "% full well>80%"

     
    warnings = {
            "Group Number Too Low?" : flag_low,
            "Group Number Too High?": flags["flag_high"],
            "Non linear?" : flag_nonl,
            "Saturated?" : flag_sat,
            "% full well high?": flag_perc, 
            "Num Groups Reset?": flags["flag_default"],
            "Minimum Integrations?": flags.get("flag_min_nint", "All good")
    }

    if 'Estimated DHS Data Excess (GB)' in timing:
        data_excess = timing['Estimated DHS Data Excess (GB)']
        if data_excess <= DHS_DATA_EXCESS_LOWER_THRESHOLD_GB:
            data_excess_warning = "All good"
        elif data_excess <= DHS_DATA_EXCESS_RECOMMENDED_LIMIT_GB:
            data_excess_warning = (
                f"Estimated data excess is {data_excess:.2f} GB, above the "
                f"{DHS_DATA_EXCESS_LOWER_THRESHOLD_GB:g} GB lower threshold. "
                "This is acceptable for DHS, but verify the setup in APT."
            )
        else:
            data_excess_warning = (
                f"Estimated data excess is {data_excess:.2f} GB, above the "
                f"{DHS_DATA_EXCESS_RECOMMENDED_LIMIT_GB:g} GB recommended "
                "limit. Verify and revise the setup in APT."
            )
        warnings['DHS Readout Optimization'] = flags.get(
            'flag_dhs_readout', 'User-specified readout pattern was not changed.'
        )
        warnings['DHS Data Excess?'] = data_excess_warning

    return warnings     
    
def add_noise_floor(noise_floor, wave_bin, error_spec):
    """Add in noise floor 
    
    This adds in a user speficied noise floor. Does not add the noise floor in quadrature 
    isntead it sets error[error<noise_floor] = noise_floor. If a wavelength dependent 
    noise floor is given and the wavelength ranges are off, it interpolates the out of 
    range noise floor. 
    
    Parameters
    ----------
    noise_floor : str or int 
        file with two column [wavelength, noise(ppm)] or single number with constant noise floor in ppm 
    wave_bin : array of float
        final binned wavelength grid from simulation 
    error_spec : array of float 
        final computed error on the planet spectrum in units of rp^2/r*^2 or fp/f*
    
    Returns
    -------
    array of float 
        error_spec-- new error
    
    Examples
    --------
    
    >>> import numpy as np
    >>> wave = np.linspace(1,2.7,10)
    >>> error = np.zeros(10)+1e-6
    >>> newerror = add_noise_floor(20, wave, error)
    >>> print(newerror)
    [  2.00000000e-05   2.00000000e-05   2.00000000e-05   2.00000000e-05
       2.00000000e-05   2.00000000e-05   2.00000000e-05   2.00000000e-05
       2.00000000e-05   2.00000000e-05]
    """
    #add user specified noise floor 
    if (type(noise_floor)==float) | (type(noise_floor) == int):
        error_spec[error_spec<noise_floor*1e-6] = noise_floor*1e-6
    elif (type(noise_floor)==str):
        read_noise = np.genfromtxt(noise_floor, dtype=(float, float), names='w, n')
        w_overlap = (wave_bin>=min(read_noise['w'])) & (wave_bin<=max(read_noise['w'])) 
        wnoise = wave_bin[w_overlap]
        noise = np.zeros(len(wave_bin))
        noise[w_overlap] = np.interp(wnoise , read_noise['w'], read_noise['n'])
        noise[(wave_bin>max(read_noise['w']))] = read_noise['n'][read_noise['w'] == max(read_noise['w'])]
        noise[(wave_bin<min(read_noise['w']))] = read_noise['n'][read_noise['w'] == min(read_noise['w'])]
        error_spec[error_spec<noise*1e-6] = noise[error_spec<noise*1e-6]*1e-6
    else: 
        raise ValueError('Noise Floor added was not integer or file')
    return error_spec

def bin_wave_to_R(w, R):
    """Creates new wavelength axis at specified resolution
    
    Parameters
    ----------
    w : list of float or numpy array of float
        Wavelength axis to be rebinned 
    R : float or int 
        Resolution to bin axis to 
    
    Returns
    -------
    list of float
        New wavelength axis at specified resolution
    
    Examples
    --------
    
    >>> newwave = bin_wave_to_R(np.linspace(1,2,1000), 10)
    >>> print((len(newwave)))
    11
    """
    wave = []
    tracker = min(w)
    i = 1 
    ind= 0
    firsttime = True
    while(tracker<max(w)):
        if i <len(w)-1:
            dlambda = w[i]-w[ind]
            newR = w[i]/dlambda
            if (newR < R) & (firsttime):
                tracker = w[ind]
                wave += [tracker]
                ind += 1
                i += 1 
                firsttime = True
            elif newR < R:
                tracker = w[ind]+dlambda/2.0
                wave +=[tracker]
                ind = (np.abs(w-tracker)).argmin()
                i = ind+1
                firsttime = True
            else:
                firsttime = False            
                i+=1    
        else:
            tracker = max(w)
            wave += [tracker]
    return np.array(wave)
    
def uniform_tophat_sum(newx,x, y):
    """Adapted from Mike R. Line to rebin spectra
    
    Sums groups of points in certain wave bin 
    
    Parameters
    ----------
    newx : list of float or numpy array of float
        New wavelength grid to rebin to 
    x : list of float or numpy array of float 
        Old wavelength grid to get rid of 
    y : list of float or numpy array of float 
        New rebinned y axis 
    
    Returns
    -------
    array of floats 
        new wavelength grid 
        
    Examples 
    --------
        
    >>> from .pandexo.engine.jwst import uniform_tophat_sum
    >>> oldgrid = np.linspace(1,3,100)
    >>> y = np.zeros(100)+10.0
    >>> newy = uniform_tophat_sum(np.linspace(2,3,3), oldgrid, y)
    >>> newy
    array([ 240.,  250.,  130.])
    """
    newx = np.array(newx)
    szmod=newx.shape[0]
    delta=np.zeros(szmod)
    ynew=np.zeros(szmod)
    delta[0:-1]=newx[1:]-newx[:-1]  
    delta[szmod-1]=delta[szmod-2] 
    #pdb.set_trace()
    for i in range(szmod-1):
        i=i+1
        loc=np.where((x >= newx[i]-0.5*delta[i-1]) & (x < newx[i]+0.5*delta[i]))
        ynew[i]=np.sum(y[loc])
    loc=np.where((x > newx[0]-0.5*delta[0]) & (x < newx[0]+0.5*delta[0]))
    ynew[0]=np.sum(y[loc])
    return ynew

def target_acq(instrument, both_spec, warning): 
    """Contains functionality to compute optimal TA strategy 

    Takes pandexo normalized flux from create_input and checks for saturation, or 
    if SNR is below the minimum requirement for each. Then adds warnings and 2d displays 
    and target acq info to final output dict 

    Parameters
    ----------
    instrument : str 
        possible options are niriss, nirspec, miri and nircam 
    both_spec : dict
        output dictionary from **create_input** 
    warning : dict 
        output dictionary from **add_warnings** 

    Retruns
    -------

    """
    out_spectrum = np.array([both_spec['wave'], both_spec['flux_out_trans']])

    #this automatically builds a default calculation 
    #I got reasonable answers for everything so all you should need to do here is swap out (instrument = 'niriss', 'nirspec','miri' or 'nircam')
    c = _build_default_calc(telescope='jwst', instrument=instrument, mode='target_acq', method='taphot')
    c['scene'][0]['spectrum']['sed'] = {'sed_type':'input','spectrum':out_spectrum}
    c['scene'][0]['spectrum']['normalization']['type'] = 'none'
    rphot = _perform_calculation(c, dict_report=True)

    #check warnings (pandeia doesn't return values for these warnings, so try will fail if all good)
    try: 
        warnings['TA Satruated?'] = rphot['warnings']['saturated']
    except:
        warnings['TA Satruated?'] = 'All good'

    try: 
        warnings['TA SNR Threshold'] = rphot['warnings']['ta_snr_threshold']
    except:
        warnings['TA SNR Threshold'] = 'All good'

    #build TA dict 
    ta = {'sn':rphot['scalar']['sn'],
            'ngroup': rphot['input']['configuration']['detector']['ngroup'], 
            'saturation':rphot['2d']['saturation']}


def _table_html(rows):
    def format_value(value):
        if isinstance(value, (float, np.floating)):
            if np.isfinite(value) and float(value).is_integer():
                return str(int(value))
            return f'{value:.6f}'
        return str(value)

    table = pd.DataFrame(rows, columns=['Parameter', 'Value'])
    table = table.set_index('Parameter')
    table.index.name = None
    table = table.to_html(formatters={'Value': format_value})
    return '<table class="table table-striped pandexo-summary-table"> \n' + table[36:len(table)]


def _jwst_instrument_name(instrument):
    names = {
        'miri': 'MIRI',
        'nircam': 'NIRCam',
        'niriss': 'NIRISS',
        'nirspec': 'NIRSpec',
    }
    return names.get(str(instrument).lower(), instrument)


def _jwst_template_name(instrument, mode, aperture, filter_name, paired_filter):
    instrument = str(instrument).lower()
    mode = str(mode).lower()
    aperture = str(aperture).lower()
    if instrument == 'miri':
        return 'MIRI Low Resolution Spectroscopy'
    if instrument == 'nirspec':
        return 'NIRSpec Bright Object Time Series'
    if instrument == 'niriss':
        return 'NIRISS Single Object Slitless Spectroscopy'
    if instrument == 'nircam':
        return 'NIRCam Grism Time Series'
    return mode


def _is_nircam_short_wave_filter(filter_name):
    return str(filter_name).lower().startswith(
        ('f070w', 'f090w', 'f115w', 'f150w', 'f150w2', 'f200w')
    )


def _nircam_channel_mode(mode, aperture, paired_filter):
    mode = str(mode).lower()
    aperture = str(aperture).lower()
    if mode == 'sw_tsgrism' or 'dhs' in aperture:
        return 'GRISM'
    if paired_filter is not None:
        return 'GRISM'
    return 'IMAGING'


def _display_subarray(subarray):
    if subarray is None:
        return None
    display_subarray = str(subarray).split(' (')[0].upper()
    if display_subarray.endswith('_PRM'):
        return f'{display_subarray[:-4]}_PRISM'
    return display_subarray


def _nircam_output_channels(subarray):
    subarray = str(subarray).lower()
    if 'noutputs=1' in subarray:
        return 1
    if subarray.startswith('subgrism'):
        return 4
    prefix = subarray.split('_', 1)[0]
    stripe_marker = prefix.rfind('s')
    if stripe_marker >= 0:
        try:
            return int(prefix[stripe_marker + 1:])
        except ValueError:
            pass
    return None


def _upper_or_none(value):
    if value is None:
        return 'None'
    return str(value).upper()


def _integer_display(value):
    """Return whole-number values as integers for display."""
    try:
        integer_value = int(value)
    except (TypeError, ValueError):
        return value
    return integer_value if integer_value == value else value


def _nircam_pupil_rows(filter_name, paired_filter):
    short_filter = None
    long_filter = None
    if filter_name is not None:
        if _is_nircam_short_wave_filter(filter_name):
            short_filter = filter_name
        else:
            long_filter = filter_name
    if paired_filter is not None:
        if _is_nircam_short_wave_filter(paired_filter):
            short_filter = paired_filter
        else:
            long_filter = paired_filter
    short_pupil_filter = 'None'
    long_pupil_filter = 'None'
    if short_filter is not None:
        short_pupil_filter = f'GDHS0+{_upper_or_none(short_filter)}'
    if long_filter is not None:
        long_pupil_filter = f'GRISMR+{_upper_or_none(long_filter)}'
    if short_filter is None and long_filter is not None:
        short_pupil_filter = 'CHOOSE THIS USING ETC'
    return [
        ('Short Pupil+Filter', short_pupil_filter),
        ('Long Pupil+Filter', long_pupil_filter),
    ]


def build_timing_display_div(out, timing):
    """Build browser-facing JWST APT-input and calculation-detail tables.

    Returns
    -------
    tuple of bytes
        HTML for the APT-input table followed by HTML for the calculation-detail
        table. The surrounding section headings are owned by ``view.html``.
    """
    configuration = out['input']['configuration']
    instrument_config = configuration['instrument']
    detector_config = configuration['detector']

    instrument = instrument_config.get('instrument')
    mode = instrument_config.get('mode')
    aperture = instrument_config.get('aperture')
    disperser = instrument_config.get('disperser')
    filter_name = instrument_config.get('filter')
    paired_filter = instrument_config.get('pandexofilterpair')
    readout_pattern = detector_config.get(
        'readout_pattern',
        detector_config.get('readmode')
    )
    nstripes = int(timing.get('Num Superstripes', 1) or 1)

    apt_rows = [
        ('Instrument', _jwst_instrument_name(instrument)),
        (
            'Template',
            _jwst_template_name(
                instrument, mode, aperture, filter_name, paired_filter
            )
        ),
    ]
    if str(instrument).lower() == 'nircam':
        apt_rows.append(
            ('SW Channel Mode', _nircam_channel_mode(mode, aperture, paired_filter))
        )
    apt_rows.append(('Subarray', _display_subarray(detector_config.get('subarray'))))
    if str(instrument).lower() == 'nircam':
        apt_rows.append(
            ('No. of Output Channels', _nircam_output_channels(detector_config.get('subarray')))
        )
    if str(instrument).lower() == 'nircam':
        apt_rows.extend(_nircam_pupil_rows(filter_name, paired_filter))
    elif str(instrument).lower() == 'nirspec':
        apt_rows.append(
            (
                'Grating/Filter',
                f'{_upper_or_none(disperser)}/{_upper_or_none(filter_name)}'
            )
        )
    elif str(instrument).lower() != 'niriss' and not (
        str(instrument).lower() == 'miri'
        and str(mode).lower() in ('lrsslitless', 'lrsslit')
    ):
        apt_rows.append(('Filter', filter_name))
    if (
        str(instrument).lower() == 'miri'
        and str(mode).lower() in ('lrsslitless', 'lrsslit')
    ):
        apt_rows.append(('Dither', 'None'))
    apt_rows.extend([
        ('Readout Pattern', _upper_or_none(readout_pattern)),
        (
            'Groups per Integration',
            timing['APT: Num Groups per Integration']
        ),
        (
            'Integrations per Occultation',
            timing['APT: Num Integrations per Occultation']
        ),
    ])

    calculation_rows = [
        ('Transit Duration (hr)', timing['Transit Duration']),
        ('Number of Transits', _integer_display(timing['Number of Transits'])),
        (
            'Transit + Baseline, No Overhead (hr)',
            timing['Transit+Baseline, no overhead (hrs)']
        ),
        ('Observing Efficiency (%)', timing['Observing Efficiency (%)']),
        ('Frame Time (sec)', timing['Seconds per Frame']),
        (
            'Integrations In Transit',
            _integer_display(timing['Num Integrations In Transit'])
        ),
        (
            'Integrations Out of Transit',
            _integer_display(timing['Num Integrations Out of Transit'])
        ),
    ]

    if 'Estimated DHS Data Excess (GB)' in timing:
        calculation_rows.extend([
            (
                'Estimated DHS Data Excess (GB)',
                timing['Estimated DHS Data Excess (GB)']
            ),
            (
                'Assumed No-TA Scheduling + Slew Overhead (sec)',
                timing['Assumed DHS Allocation Overhead (sec)']
            ),
        ])

    if nstripes > 1:
        calculation_rows.extend([
            ('Number of Stripes', nstripes),
            (
                'Elapsed Time per APT Integration incl. Reset (sec)',
                timing['Time/Integration incl reset (sec)']
            ),
            (
                'Science Time per Full Multistripe Cycle excl. Reset (sec)',
                timing['Measurement Time per Integration (sec)']
            ),
            (
                'Science Time per Stripe excl. Reset (sec)',
                timing['Measurement Time per Integration (sec)'] / float(nstripes)
            ),
        ])
    else:
        calculation_rows.extend([
            (
                'Elapsed Time per Integration incl. Reset (sec)',
                timing['Time/Integration incl reset (sec)']
            ),
            (
                'Science Time per Integration excl. Reset (sec)',
                timing['Measurement Time per Integration (sec)']
            ),
        ])

    return _table_html(apt_rows).encode(), _table_html(calculation_rows).encode()


def as_dict(out, both_spec ,binned, timing, mag, sat_level, warnings, punit, unbinned,calculation): 
    """Format dictionary for output data 
    
    Takes all output from jwst run and converts it to simple dictionary 
    
    Parameters
    ----------
    out : dict 
        output dictionary from **compute_out**
    both_spec : dict 
        output dictionary from **createInput.bothTrans**
    binned : dict 
        dictionary from **wrapper** 
    timing : dict 
        dictionary from **compute_timing**
    mag : dict 
        magnitude of system 
    sat_level : float or int
        saturation level in electrons 
    warnings : dict 
        warning dictionary from **add_warnings**
    punit : "fp/f*" or "rp^2/r*^2"
        unit of supplied spectra options are: only options are fp/f* or rp^2/r*^2
    unbinned : dict 
        unbinned raw data from **wrapper**
    calculation : str 
        noise calculation type
    
    Returns
    -------
    dict 
        compressed dictionary 
    """
    #for emission spectrum

    p=1.0
    if punit == 'fp/f*': p = -1.0
    apt_div, calculation_div = build_timing_display_div(out, timing)

    warnings_div = pd.DataFrame.from_dict(warnings, orient='index')
    warnings_div.columns = ['Value']
    warnings_div = warnings_div.to_html()
    warnings_div = (
        '<table class="table table-striped pandexo-summary-table"> \n'
        + warnings_div[36:len(warnings_div)]
    )
    warnings_div = warnings_div.encode()
    
    map_dhs_names = {'sub40stripe1_dhs':'SUB40S1_2-SPECTRA',
                     'sub80stripe2_dhs':'SUB80S2_4-SPECTRA',
                     'sub160stripe4_dhs':'SUB160S4_8-SPECTRA',
                     'sub256stripe4_dhs':'SUB256S4_8-SPECTRA'
            }
    
    subarray = out['input']['configuration']['detector']['subarray']
    for idhs in map_dhs_names.keys(): 
        subarray = subarray.replace(idhs, f'{idhs} (ETC Name)/ {map_dhs_names[idhs]} (APT Name)')

    input_dict = {
   	 "Target Mag": mag , 
   	 "Saturation Level (electons)": sat_level, 
   	 "Instrument": out['input']['configuration']['instrument']['instrument'], 
   	 "Mode": out['input']['configuration']['instrument']['mode'], 
   	 "Aperture": out['input']['configuration']['instrument']['aperture'], 
   	 "Disperser": out['input']['configuration']['instrument']['disperser'], 
   	 "Subarray": subarray, 
   	 "Readmode": out['input']['configuration']['detector']['readout_pattern'], 
 	 "Filter": out['input']['configuration']['instrument']['filter'],
 	 "Primary/Secondary": punit
    }
    
    input_div = pd.DataFrame.from_dict(input_dict, orient='index')
    input_div.columns = ['Value']
    input_div = input_div.to_html()
    input_div = (
        '<table class="table table-striped pandexo-summary-table"> \n'
        + input_div[36:len(input_div)]
    )
    input_div = input_div.encode()
    
    #add calc type to input dict (doing it here so it doesn't output on webpage
    input_dict["Calculation Type"]= calculation
    
    final_dict = {
    'OriginalInput': {'model_spec':both_spec['model_spec'],
                     'model_wave' : both_spec['model_wave'],
                     'star_spec': both_spec['flux_out_trans']},
    'RawData': unbinned,
    'FinalSpectrum': binned,
        
    #pic output 
    'PandeiaOutTrans': out, 

    #all timing info 
    'timing': timing,
    'warning':warnings,
    'input':input_dict,
    
    #divs for html rendering    
    'timing_div': apt_div + b'\n' + calculation_div,
    'apt_div': apt_div,
    'calculation_div': calculation_div,
    'input_div':input_div,
    'warnings_div':warnings_div,
    }
    return final_dict

    
