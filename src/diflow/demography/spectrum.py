"""Prepare observed frequency spectra for dadi likelihood fitting."""

from __future__ import annotations

import numpy as np

from .dadi_backend import _require_dadi


def prepare_observed_spectrum(observed_spectrum, *, polarized: bool = False):
    """Create a dadi Spectrum with correct polarization treatment.

    Ordinary VCF REF/ALT coding does not identify ancestral and derived alleles.
    DIFLOW therefore defaults to an unpolarized analysis and folds the spectrum.
    Set polarized=True only when allele orientation has been established from
    appropriate ancestral-state information upstream.
    """
    dadi = _require_dadi()
    data_array = np.asarray(observed_spectrum, dtype=float)
    if data_array.ndim != 2:
        raise ValueError("observed_spectrum must be a two-dimensional array.")
    if np.any(~np.isfinite(data_array)) or np.any(data_array < 0):
        raise ValueError("observed_spectrum must contain finite non-negative values.")
    if data_array.sum() <= 0:
        raise ValueError("observed_spectrum must contain positive total mass.")

    spectrum = dadi.Spectrum(data_array, mask_corners=True)
    if not polarized:
        spectrum = spectrum.fold()
    return spectrum
