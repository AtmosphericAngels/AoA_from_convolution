"""
Run convolution method to calculate mean age of air from observations.

This module is a wrapper for Conc2Age_Convolution to automatically include
trace gas mixing ratio reference time series.

"""

import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
from convolution_method_functions import Conc2Age_Convolution


# Tropical time series for CO2, SF6 and CH4 from Eric Ray
CO2_ref_file = "tropical_time_series/CO2_tropical_surface_time_series_Ray.txt"
CH4_ref_file = "tropical_time_series/CH4_tropical_surface_time_series_Ray.txt"
SF6_ref_file = "tropical_time_series/SF6_tropical_surface_time_series_Ray.txt"

# Alternative tropical time series for SF6 from Andreas Engel and Harald Boensch
SF6_ref_file = "tropical_time_series/SF6_tropical_surface_time_series_Engel.txt"


# load tropical time series for CO2
CO2_ref = pd.read_csv(
    CO2_ref_file, delim_whitespace=True, header=5, names=["time", "vmr"]
)
CO2_ref_t = np.array(CO2_ref["time"])
CO2_ref_vmr = np.array(CO2_ref["vmr"])

# load tropical time series for CH4
CH4_ref = pd.read_csv(
    CH4_ref_file, delim_whitespace=True, header=4, names=["time", "vmr"]
)
CH4_ref_t = np.array(CH4_ref["time"])
CH4_ref_vmr = np.array(CH4_ref["vmr"])

# load tropical time series for SF6
SF6_ref = pd.read_csv(
    SF6_ref_file, delim_whitespace=True, header=7, names=["time", "vmr"]
)
SF6_ref_t = np.array(SF6_ref["time"])
SF6_ref_vmr = np.array(SF6_ref["vmr"])


# %%


def SF6_to_AoA(t_obs, SF6_obs, rom, res="G"):
    """Calculate mean age from SF6 observations.

    AGE OF AIR VALUES BELOW ONE YEAR CAN NOT BE CALCULATED AND WILL RESULT IN MISSING VALUES (nan)!

    Age of Air values above 10 years can not be calculated and will result in missing values as well.

    Conc2Age_Convolution only works with one-dimensional numpy arrays!

    Parameters:
    ---------
    t_obs: float or int
        time of observation in years (full year plus fraction of year e.g. 2005.6573).
    SF6_obs: float or array
        SF6 mixing ratios at time of observation t_obs in pptV.
    rom: float or int
        ratio of first and second moment of inverse gaussian distribution (age spectrum).
        A rom of 1.2 years is a good starting point in order to get going.
    res: str
        which rersolution is used for the forward calculation?
        c (Tracer Time series) or G (age spectrum, then use
        calculate_age_spectrum_1d and the default resolution in there).
        The default is G.
    """
    SF6_AoA = Conc2Age_Convolution(
        t_ref=SF6_ref_t, c_ref=SF6_ref_vmr, t_obs=t_obs, c_obs=SF6_obs, rom=rom, res=res
    )
    return SF6_AoA


def CO2_to_AoA(t_obs, CO2_obs, rom, CH4_obs=None, res="G", verbose=False):
    """Calculate mean age from CO2 observations.

    MAKE SURE TO USE THE CORRECT UNITS FOR CO2_obs (ppmV) AND CH4_obs (ppbV)!

    AGE OF AIR VALUES BELOW ONE YEAR CAN NOT BE CALCULATED AND WILL RESULT IN MISSING VALUES (nan)!
    Age of air values below 2.5 years should be treated very carefully due to the
    seasonality of CO2.

    Age of Air values above 10 years can not be calculated and will result in missing values as well.

    AGE OF AIR CORRECTION FOR SEASONALITY of CO2 NOT YET IMPLEMENTED!

    CO_AoA_raw contains the Age of Air values calculated without (!) methane correction
    for the observed CO2 mixing ratios

    Conc2Age_Convolution only works with one-dimensional numpy arrays!

    Parameters:
    ---------
    t_obs: float or int
        time of observation in years (full year plus fraction of year e.g. 2005.6573).
    CO2_obs: float or array
        CO2 mixing ratios at time of observation t_obs in ppmV.
    rom: float or int
        ratio of first and second moment of inverse gaussian distribution (age spectrum).
        A rom of 1.2 years is a good starting point in order to get going.
    CH4_obs: float or array
        CH4 mixing ratios at time of observation t_obs in pptV. The default is None.
    res: str
        which rersolution is used for the forward calculation?
        c (Tracer Time series) or G (age spectrum, then use
        calculate_age_spectrum_1d and the default resolution in there).
        The default is G.
    verbose: bool
        Print additional information on function call. The default is False
    """
    CO2_AoA_raw = Conc2Age_Convolution(
        t_ref=CO2_ref_t, c_ref=CO2_ref_vmr, t_obs=t_obs, c_obs=CO2_obs, rom=rom, res=res
    )

    if CH4_obs is None:
        print("WARNING! No observed CH4 values were given!")
        print(
            "CH4 correction is not possible and the calculated Age of Air values will be too low!"
        )
        return CO2_AoA_raw

    else:
        # ========================================================================================
        # SIMPLIFIED METHANE CORRECTION
        # A more sophisticated methane correction ought to be implemented in the (near) future
        # ========================================================================================

        # Difference in time of observation and uncorrected Age of Air yields time in tropospheric record to take CH4 mixing ratios from
        t_CH4 = t_obs - CO2_AoA_raw

        # CH4 in ppmV transported from the troposphere to the point and time of observation in the stratosphere!
        CH4_tropo = (
            interp1d(CH4_ref_t, CH4_ref_vmr, bounds_error=None, fill_value=np.nan)(
                t_CH4
            )
            / 1000
        )

        # Difference (in ppmV) between total CH4 transported from troposphere and observed CH4 still present at point and time of observation
        CO2_from_CH4 = CH4_tropo - CH4_obs / 1000
        CO2_from_CH4 = np.asarray([CO2_from_CH4]).flatten()

        # Check if CO2_from_CH4 contains any negative i.e. nonsensical values
        CO2_from_CH4_fail_ind = np.where(CO2_from_CH4 < 0)[0]

        if len(CO2_from_CH4_fail_ind) > 0:
            print(
                f"WARNING! OBSERVED CH4 LARGER THAN TRANSPORTED CH4 FOR CH4_obs "
                f"{CH4_obs[CO2_from_CH4_fail_ind]} in {t_obs}"
            )
            print(
                "CO2_from_CH4 will be replaced with missing values (nan) "
                "for above values of t_obs"
            )
            print(
                "Please check if the correct units are used for "
                "CO2_obs (ppmV) and CH4_obs (ppbV)"
            )
            CO2_from_CH4[CO2_from_CH4_fail_ind] = np.nan

        if verbose:
            print(CO2_from_CH4, " ppmV CO2 created from CH4 for given t_obs")

        # Difference in observed CO2 and estimated CO2 created from CH4 yields corrected CO2 values
        CO2_corr = CO2_obs - CO2_from_CH4

        # Age of air from corrected CO2 mixing ratios
        CO2_AoA_corr = Conc2Age_Convolution(
            t_ref=CO2_ref_t,
            c_ref=CO2_ref_vmr,
            t_obs=t_obs,
            c_obs=CO2_corr,
            rom=rom,
            res="G",
        )
        return CO2_AoA_corr
