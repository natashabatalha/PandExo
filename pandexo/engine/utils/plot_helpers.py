import numpy as np


def add_wavelength_gap_breaks(x, y, gap_factor=5.0):
    """Insert visual line breaks at large gaps in a wavelength grid.

    Bokeh line and step glyphs draw a continuous segment between adjacent
    finite points. For NIRSpec spectra, the valid NRS1 and NRS2 samples can be
    adjacent in the plotted arrays even though a detector gap separates them on
    the wavelength axis. This helper inserts an extra point with a finite
    midpoint wavelength and a NaN ordinate after any isolated, unusually large
    positive wavelength step, causing Bokeh to break the rendered line at that
    location.

    The input arrays are otherwise left unchanged. This function is intended
    for plotting only; the returned arrays should not be reused for numerical
    calculations or binning.

    Parameters
    ----------
    x : array-like
        Wavelength coordinates, expected to be ordered monotonically for normal
        spectral plots. Non-finite values are ignored when estimating local
        wavelength spacings.
    y : array-like
        Values to plot at each wavelength. Must have the same shape as ``x``.
    gap_factor : float, optional
        Dimensionless multiplier on the local wavelength spacing, not a value
        in microns. A positive wavelength step larger than ``gap_factor`` times
        the median of nearby positive finite wavelength steps is treated as a
        plotting gap. For example, if the neighboring spacings are about 0.01
        microns, the default value of 5.0 breaks the plotted line at jumps
        larger than about 0.05 microns. Because the reference spacing is local,
        smooth dispersion changes across an instrument mode should not be
        treated as gaps.

    Returns
    -------
    tuple of numpy.ndarray
        The possibly expanded ``(x, y)`` arrays. If no gap is detected, or the
        inputs are too short or shape-incompatible, the arrays are returned
        after conversion to floating point NumPy arrays.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    if x.size < 3 or y.shape != x.shape:
        return x, y

    steps = np.diff(x)
    if not np.any(np.isfinite(steps) & (steps > 0.0)):
        return x, y

    def local_typical_step(step_index, radius=3):
        start = max(0, step_index - radius)
        stop = min(steps.size, step_index + radius + 1)
        local_steps = np.delete(steps[start:stop], step_index - start)
        local_steps = local_steps[np.isfinite(local_steps) & (local_steps > 0.0)]
        if local_steps.size == 0:
            return np.nan
        return np.median(local_steps)

    xout = []
    yout = []
    for i, (xi, yi) in enumerate(zip(x, y)):
        xout.append(xi)
        yout.append(yi)
        if i == x.size - 1:
            continue
        step = steps[i]
        typical_step = local_typical_step(i)
        if (
            np.isfinite(step)
            and step > 0.0
            and np.isfinite(typical_step)
            and typical_step > 0.0
            and step > gap_factor * typical_step
        ):
            xout.append(0.5 * (xi + x[i + 1]))
            yout.append(np.nan)

    return np.asarray(xout), np.asarray(yout)
