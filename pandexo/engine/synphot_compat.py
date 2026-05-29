"""Compatibility helpers for PandExo's synphot/stsynphot usage.

This module keeps low-level ``synphot`` and ``stsynphot`` details out of the
input-building code.
"""

import numpy as np
import astropy.units as u
from synphot import SourceSpectrum, SpectralElement, units
from synphot.models import Empirical1D
import stsynphot
from stsynphot import catalog


def _wave_unit(name):
    unit_map = {
        'um': u.um,
        'micron': u.um,
        'nm': u.nm,
        'cm': u.cm,
        'Angs': u.Angstrom,
        'angstrom': u.Angstrom,
        'angstroms': u.Angstrom,
        'Hz': u.Hz,
        'hz': u.Hz,
    }
    try:
        return unit_map[name]
    except KeyError:
        raise Exception('Units are not correct. Pick um, nm, cm, hz, or Angs')


def _flux_quantity(flux, name):
    flux = np.asarray(flux, dtype=float)
    if name == 'Jy':
        return flux * u.Jy
    if name == 'jy':
        return flux * u.Jy
    if name == 'FLAM':
        return flux * units.FLAM
    if name == 'erg/cm2/s/Hz':
        return flux * 1e23 * u.Jy
    raise Exception('Units are not correct. Pick FLAM or Jy or erg/cm2/s/Hz')


def make_array_spectrum(wave, flux, wave_unit, flux_unit):
    """Create a synphot ``SourceSpectrum`` from wavelength and flux arrays.

    Parameters are the same user-facing units accepted by PandExo. Frequency
    inputs are converted onto a wavelength grid before constructing the
    spectrum so the internal spectral axis is sorted by ascending wavelength.
    """

    wave_quantity = np.asarray(wave, dtype=float) * _wave_unit(wave_unit)
    flux_quantity = _flux_quantity(flux, flux_unit)

    wave_angstrom = wave_quantity.to(u.Angstrom, equivalencies=u.spectral())
    order = np.argsort(wave_angstrom.value)
    return SourceSpectrum(
        Empirical1D,
        points=wave_angstrom[order],
        lookup_table=flux_quantity[order],
    )


def load_phoenix_spectrum(teff, metallicity, logg):
    """Load/interpolate a PHOENIX model spectrum with ``stsynphot``."""

    return catalog.grid_to_spec("phoenix", teff, metallicity, logg)


def load_bandpass_from_file(path):
    """Load a ``synphot`` bandpass from a FITS throughput file."""

    return SpectralElement.from_file(path)


def renormalize_to_vegamag(spectrum, mag, bandpass):
    """Normalize a spectrum to a VEGAMAG value through a bandpass."""

    return spectrum.normalize(
        float(mag) * units.VEGAMAG,
        band=bandpass,
        vegaspec=stsynphot.Vega,
    )


def sample_spectrum_micron_mjy(spectrum, wavelengths=None):
    """Return wavelength in microns and flux in mJy for PandExo/Pandeia."""

    spectral_axis = spectrum.waveset if wavelengths is None else wavelengths
    wave_micron = spectral_axis.to(u.micron, equivalencies=u.spectral())
    flux_mjy = spectrum(spectral_axis).to(
        u.mJy,
        equivalencies=u.spectral_density(wave_micron),
    )

    order = np.argsort(wave_micron.value)
    return wave_micron.value[order], flux_mjy.value[order]


def calculate_bin_edges(centers):
    """Calculate midpoint bin edges from a 1D array of bin centers."""

    centers = np.asarray(centers, dtype=float)
    if centers.ndim != 1 or centers.size < 2:
        raise ValueError("centers must be a 1D array with at least two values")

    edges = np.empty(centers.size + 1, dtype=float)
    edges[1:-1] = 0.5 * (centers[:-1] + centers[1:])
    edges[0] = centers[0] - (edges[1] - centers[0])
    edges[-1] = centers[-1] + (centers[-1] - edges[-2])
    return edges
