"""
convolution_method_functions.py

This module provides functions to calculate the mean age of air from suitable
observations and a reference time series.

Main function: Conc2Age_Convolution

"""

# Standard library imports
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate


def Conc2Age_Convolution(
    t_ref,
    c_ref,
    t_obs,
    c_obs,
    rom,
    t_int=30,
    cutoff=1.0,
    positiv=False,
    comment=False,
    res="c",
):
    """Calculate mean age from tracer observations.

    based on Function con2Age_Convolution written in IDL by Harald Boenisch
    requires function calculate_lag for age values below 1 year.

    known issues:
    - For mean age below 1 year the lag time is used instead of mean age

    last modified:
    - Feb 19, 2021, Thomas Wagenhäuser
    - Apr 21, 2020, Andreas Engel

    Parameters:
    --------
    t_ref: list or np.array, floats
        reference time series time
    c_ref: list or np.array, floats
        reference time series mixing ratios
    t_obs: float
        time of observation
    c_obs: float array or float
        observed mixing ratios
    rom: float
        ratio of moments; has no default values
    t_int: int or float
        Integration period; default is 30 years
    cutoff: int or float, optional
        value below which lag time will be used. The default is 1 year.
    positiv: bool, optional
        set this keyword to replace all negative values by 0 (not tested)
    comment: bool, optional
        print some info to the console. The default is False.
    res: which rersolution is used for the forward calculation?
        c (Tracer Time series) or G (age spectrum, then use
        calculate_age_spectrum_1d and the default resolution in there).
        The default is c.

    Returns:
    --------
    mean age, array of the same dimensions as c_obs or a float

    """
    c_ref = np.asarray([c_ref]).flatten()
    c_obs = np.asarray([c_obs]).flatten()

    # a_obs (age observed) contains the sought age of air values for every
    # observed mixing ratio (c_obs). It is created as an array with the same
    # shape as "c_obs" that contains only "nan" which will be replaced later on.
    a_obs = c_obs * np.nan

    # vd (short for "valid") contains the arguments of valid values in
    vd = np.isfinite(c_obs)
    # nvd is the number of valid values in c_obs
    nvd = np.sum(vd)

    if nvd > 0:
        if comment:
            print(
                "Mean age calculation with Folding Method; lag time is used "
                "for mean age below the number of years specified as cutoff (default 1 year)"
            )

        # a_tmp defines the range and resolution with which age of air can be calculated
        # Age of air can be calculated from 1 to 10 years with an interpolation step of 1/10 years.
        a_tmp = np.linspace(cutoff, 10, int((10 - cutoff) * 10 + 1))

        c_tmp = a_tmp * np.nan

        wt = (t_ref < (t_obs)) & (t_ref > (t_obs - t_int))
        t_tmp = t_obs - t_ref[wt]

        t = np.flip(t_tmp, 0)

        if res == "c":
            # default: use resolution from observations
            for n in np.arange(len(a_tmp)):
                c_ref_rev = np.flip(c_ref[wt], 0)
                age = a_tmp[n]
                width = (rom * age) ** 0.5
                G = _inverse_gaussian(age, width, t)
                Int_G = np.trapz(G, t)
                G = G / Int_G
                # Modification to account for shifted mean AoA from normalization
                AoA_G_fit = np.trapz(G * t, t)
                time_shift = age - AoA_G_fit
                t_shift = t + time_shift
                c_ref_rev_int_shift = np.interp(t_shift, t, c_ref_rev)

                c_tmp[n] = np.trapz(c_ref_rev_int_shift * G, t)

        else:
            # use resolution from G
            for n in np.arange(len(a_tmp)):
                age = a_tmp[n]
                G = Calculate_AgeSpectrum_1D(age, rom)

                fc_int = interpolate.interp1d(
                    t_tmp, c_ref[wt], fill_value="extrapolate"
                )

                c_int = fc_int(G[:, 0])
                Int_G = np.trapz(G[:, 1], G[:, 0])
                G[:, 1] = G[:, 1] / Int_G
                AoA_G_fit = np.trapz(G[:, 1] * G[:, 0], G[:, 0])
                time_shift = age - AoA_G_fit
                t_G_shift = G[:, 0] + time_shift
                c_int_shift = np.interp(t_G_shift, G[:, 0], c_int)
                # c_tmp[n] = np.trapz(G[:,1] * c_int, G[:,0])
                c_tmp[n] = np.trapz(G[:, 1] * c_int_shift, G[:, 0])

        c2a = interpolate.interp1d(np.flip(c_tmp, 0), np.flip(a_tmp, 0))

        if nvd == 1:
            # print(c_obs, max(c_tmp), min(c_tmp))
            if (c_obs > min(c_tmp)) and (c_obs < max(c_tmp)):
                a_obs = c2a(c_obs)[0]
                if comment:
                    print("single value, folding age")
                    print(a_obs)
            else:
                a_obs = np.nan
                if comment:
                    print("single value outside of range")
                    print("a_obs:", a_obs)
                    print("c_obs:", c_obs)
            if positiv:
                a_obs = a_obs.clip(min=0)
        else:  # if array of values
            vd2 = (c_obs >= min(c_tmp)) & (c_obs <= max(c_tmp))
            a_obs = np.zeros_like(c_obs) + np.nan
            a_obs[vd2] = c2a(c_obs[vd2])
            if comment:
                print(
                    "array of values; folding age for age > 1 year; "
                    "np.nan below and for values out of range"
                )
                print(a_obs)
            if positiv:
                a_obs = a_obs.clip(min=0)

    return a_obs


# %%
def _inverse_gaussian(m_age, width, t):
    G = (m_age ** 3 / (4 * np.pi * width ** 2 * t ** 3)) ** 0.5 * np.exp(
        -m_age * (t - m_age) ** 2 / (4 * width ** 2 * t)
    )
    return G

# %%


def Calculate_AgeSpectrum_1D(m_age, rom, nyy=30, nmnth=12, per_month=30, plot=False):
    """Calculate transit time probability distribution.

    based on Calculate_AgeSpectrum_1D by Harald Boenisch written in IDL
    Age spectrum for 1-D flow advection model with diffusion :

    G(t) = m_age³/(4*pi*w_age²*t³) * exp[ -m_age*(t-m_age)²/(4*w_age²*t) ]

    Parameters:
    --------
    m_age : float
        mean age is the first moment of the age spectrum distribution
    w_age : float
        width age is the second moment of the m_age centered age
        spectrum distribution
    rom : float
        ratio of moments
    nyy : int, optional
        Number of years over which to calculate the age spectrum. The default is 30.
    nmth : int, optional
        Number of months to use (should alway be 12). The default is 12.
    per_month: int, optional
        Number of data Points per Month. The default is 10.
    plot : bool, optional
        set to true if you want a plot of the spectrum. The default is False.

    Returns:
    --------
    Agespectrum G, 2 dimensional np.array. G[:, 0] contains the transit times,
    G[:, 1] contains the probabilities

    """
    w = (rom * m_age) ** 0.5

    nt = nyy * nmnth * per_month
    dt = 1 / (nt / nyy)
    t = dt * (np.arange(nt) + 1)

    G = np.zeros((len(t) + 1, 2))

    G[1:, 0] = t
    G[1:, 1] = _inverse_gaussian(m_age, w, t)

    G[:, 1] = G[:, 1] / np.trapz(G[:, 1], x=G[:, 0])
    if plot:
        fig, ax = plt.subplots()
        ax.plot(G[:, 0], G[:, 1])
        ax.set(
            xlabel="transit time (years)",
            ylabel="Probability (1/year)",
            title="Age Spectrum, Mean Age: "
            + str(m_age)
            + ", Ratio of moments: "
            + str(rom),
        )
        ax.grid()

    return G
