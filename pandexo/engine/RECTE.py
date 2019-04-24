# ! /usr/bin/env pythonp
"""ramp effect model

Version 1.0.0: published version

Version 0.1: Adapted original IDL code (From D. Apai) to python by Yifan Zhou

"""
import numpy as np
import itertools


def RECTE(
        cRates,
        tExp,
        exptime=180,
        trap_pop_s=0,
        trap_pop_f=0,
        dTrap_s=50,
        dTrap_f=10,
        dt0=0,
        lost=0,
        mode='scanning'
):
    """This function calculates HST/WFC3/IR ramp effect profile based on
the charge trapping explanation developed in Zhou et al. (2017).

    Parameters
    ----------
    cRates : numpy.array
		intrinsic count rate of each exposures, unit: e/s

    tExp : numpy.array
		time stamps for the exposures, unit: seconds

    exptime : numpy.array or float
		(default 180 seconds) exposure time

    trap_pop_s : float or numpy.array
		(default 0) number of occupied slow population
        charge traps before the very beginning of the observation

	trap_pop_f: float or numpy.array
		(default 0) number of occupied fast population
        charge traps before the very beginning of the observation

    dTrap_s: float or numpy.array
		(default [0]) number of additional charges trapped
        by slow population traps during earth occultation

    dTrap_f : float or numpy.array
		(default [0]) number of additional charges trapped
        by fast population traps during earth occultation

    dt0: float
		(default 0) exposure time before the very beginning
        of the observation. It could be due to guidence adjustment

    lost: float
		(default 0, no lost) fraction of trapped electrons that are
        not eventually detected

    mode : string
		(default scanning, scanning or staring, or others),
        for scanning mode observation , the pixel no longer receive
        photons during the overhead time, in staring mode,
        the pixel keps receiving elctrons

	Returns
	-------
	numpy.array
		observed counts


    Example
	-------

    see Examples and Cookbook

    """
    nTrap_s = 1525.38  # 1320.0
    eta_trap_s = 0.013318  # 0.01311
    tau_trap_s = 1.63e4
    nTrap_f = 162.38
    eta_trap_f = 0.008407
    tau_trap_f = 281.463
    try:
        dTrap_f = itertools.cycle(dTrap_f)
        dTrap_s = itertools.cycle(dTrap_s)
        dt0 = itertools.cycle(dt0)
    except TypeError:
        dTrap_f = itertools.cycle([dTrap_f])
        dTrap_s = itertools.cycle([dTrap_s])
        dt0 = itertools.cycle([dt0])
    obsCounts = np.zeros(len(tExp))
    trap_pop_s = min(trap_pop_s, nTrap_s)
    trap_pop_f = min(trap_pop_f, nTrap_f)
    for i in range(len(tExp)):
        try:
            dt = tExp[i+1] - tExp[i]
        except IndexError:
            dt = exptime
        f_i = cRates[i]
        c1_s = eta_trap_s * f_i / nTrap_s + 1 / tau_trap_s  # a key factor
        c1_f = eta_trap_f * f_i / nTrap_f + 1 / tau_trap_f
        # number of trapped electron during one exposure
        dE1_s = (eta_trap_s * f_i / c1_s - trap_pop_s) * \
            (1 - np.exp(-c1_s * exptime))
        dE1_f = (eta_trap_f * f_i / c1_f - trap_pop_f) * \
            (1 - np.exp(-c1_f * exptime))
        dE1_s = min(trap_pop_s + dE1_s, nTrap_s) - trap_pop_s
        dE1_f = min(trap_pop_f + dE1_f, nTrap_f) - trap_pop_f
        trap_pop_s = min(trap_pop_s + dE1_s, nTrap_s)
        trap_pop_f = min(trap_pop_f + dE1_f, nTrap_f)
        obsCounts[i] = f_i * exptime - dE1_s - dE1_f
        if dt < 5 * exptime:  # whether next exposure is in next batch of exposures
            # same orbits
            if mode == 'scanning':
                # scanning mode, no incoming flux between exposures
                dE2_s = - trap_pop_s * (1 - np.exp(-(dt - exptime)/tau_trap_s))
                dE2_f = - trap_pop_f * (1 - np.exp(-(dt - exptime)/tau_trap_f))
            elif mode == 'staring':
                # for staring mode, there is flux between exposures
                dE2_s = (eta_trap_s * f_i / c1_s - trap_pop_s) * \
                    (1 - np.exp(-c1_s * (dt - exptime)))
                dE2_f = (eta_trap_f * f_i / c1_f - trap_pop_f) * \
                    (1 - np.exp(-c1_f * (dt - exptime)))
            else:
                # others, same as scanning
                dE2_s = - trap_pop_s * (1 - np.exp(-(dt - exptime)/tau_trap_s))
                dE2_f = - trap_pop_f * (1 - np.exp(-(dt - exptime)/tau_trap_f))
            trap_pop_s = min(trap_pop_s + dE2_s, nTrap_s)
            trap_pop_f = min(trap_pop_f + dE2_f, nTrap_f)
        elif dt < 1200:
            # considering in orbit download scenario
            trap_pop_s = min(
                trap_pop_s * np.exp(-(dt-exptime)/tau_trap_s), nTrap_s)
            trap_pop_f = min(
                trap_pop_f * np.exp(-(dt-exptime)/tau_trap_f), nTrap_f)
        else:
            # switch orbit
            dt0_i = next(dt0)
            trap_pop_s = min(trap_pop_s * np.exp(-(dt-exptime-dt0_i)/tau_trap_s) +
                             next(dTrap_s), nTrap_s)
            trap_pop_f = min(trap_pop_f * np.exp(-(dt-exptime-dt0_i)/tau_trap_f) +
                             next(dTrap_f), nTrap_f)
            f_i = cRates[i + 1]
            c1_s = eta_trap_s * f_i / nTrap_s + 1 / tau_trap_s  # a key factor
            c1_f = eta_trap_f * f_i / nTrap_f + 1 / tau_trap_f
            dE3_s = (eta_trap_s * f_i / c1_s - trap_pop_s) * \
                (1 - np.exp(-c1_s * dt0_i))
            dE3_f = (eta_trap_f * f_i / c1_f - trap_pop_f) * \
                (1 - np.exp(-c1_f * dt0_i))
            dE3_s = min(trap_pop_s + dE3_s, nTrap_s) - trap_pop_s
            dE3_f = min(trap_pop_f + dE3_f, nTrap_f) - trap_pop_f
            trap_pop_s = min(trap_pop_s + dE3_s, nTrap_s)
            trap_pop_f = min(trap_pop_f + dE3_f, nTrap_f)
        trap_pop_s = max(trap_pop_s, 0)
        trap_pop_f = max(trap_pop_f, 0)

    return obsCounts


if __name__ == '__main__':
    pass
