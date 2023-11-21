import numpy as np
import pandas as pd
from pathlib import Path
import time

from toolpac.age.Conc2Age_Convolution import age_from_conv_old, age_from_conv

from convolution_method_functions import Conc2Age_Convolution
from run_convolution_method import *


fdir = Path('tropical_time_series')
SF6_ref_file = fdir / "SF6_tropical_surface_time_series_Engel.txt"

SF6_ref = pd.read_csv(SF6_ref_file, delim_whitespace=True, header=7, names=["time", "vmr"])

SF6_ref

SF6_ref_t = np.array(SF6_ref["time"])
SF6_ref_vmr = np.array(SF6_ref["vmr"])


t_obs = 2020
SF6_obs = 9.5
rom = 1.2

SF6_to_AoA(t_obs, SF6_obs, rom)


t_obs = 2020
SF6_obs = np.ones(10) * 9.5
rom = 1.2

SF6_to_AoA(t_obs, SF6_obs, rom)

t_obs = np.ones(10) * 2020
SF6_obs = np.ones(10) * 9.5
rom = 1.2

SF6_to_AoA(t_obs, SF6_obs, rom)

t_obs = np.arange(0, 1, 0.1) + 2020
SF6_obs = np.ones(10) * 9.5
rom = 1.2

SF6_to_AoA(t_obs, SF6_obs, rom)

t_int = 2
wt = np.where(
    np.logical_and((SF6_ref_t < (t_obs)), (SF6_ref_t > (t_obs - t_int)))
)[
    0
]  # why is here a "0"?

wt.shape

wt = np.where(
    np.logical_and((SF6_ref_t < (t_obs)), (SF6_ref_t > (t_obs - t_int))))

_wt = np.logical_and((SF6_ref_t < (t_obs)), (SF6_ref_t > (t_obs - t_int)))
_wt
len(wt)

_wt_bool = (SF6_ref_t < (t_obs)) & (SF6_ref_t > (t_obs - t_int))

SF6_ref_t[wt]

SF6_ref_t[_wt]

SF6_ref_t[_wt_bool]

# %%
SF6_obs = np.ones(10) * 9.5
SF6_obs[4] = np.nan
vd = np.isfinite(SF6_obs)
vd
# nvd is the number of valid values in c_obs
nvd = np.sum(vd)
nvd

# %%
t_obs
SF6_obs
rom

age_from_conv_old(SF6_ref_t, SF6_ref_vmr, t_obs, SF6_obs, rom)

t_obs = 2018
SF6_obs = 8.5
rom

SF6_to_AoA(t_obs, SF6_obs, rom, res="G")
SF6_to_AoA(t_obs, SF6_obs, rom, res="c")

pd.Series(SF6_ref_t)
SF6_ref_vmr.shape

np.asarray([SF6_ref_vmr]).flatten().shape

t_int = 2
wt = (SF6_ref_t < (t_obs)) & (SF6_ref_t > (t_obs - t_int))
wt
SF6_ref_vmr[wt]
np.asarray([SF6_ref_vmr]).flatten()[wt]

age_from_conv(SF6_ref_t, SF6_ref_vmr, t_obs, SF6_obs, rom)
age_from_conv_old(SF6_ref_t, SF6_ref_vmr, t_obs, SF6_obs, rom)
Conc2Age_Convolution(SF6_ref_t, SF6_ref_vmr, t_obs, SF6_obs, rom, res="G")
Conc2Age_Convolution(SF6_ref_t, SF6_ref_vmr, t_obs, SF6_obs, rom, res="c")

wt[wt]

# %%

t_obs = 2018
SF6_obs = np.arange(0, 1, 0.1) + 8.3
rom = 1.2

# %%
%%timeit
SF6_to_AoA(t_obs, SF6_obs, rom, res="G")
#%%
%%timeit
SF6_to_AoA(t_obs, SF6_obs, rom, res="c")
#%%
%%timeit
age_from_conv(SF6_ref_t, SF6_ref_vmr, t_obs, SF6_obs, rom)
#%%
%%timeit
age_from_conv(SF6_ref_t, SF6_ref_vmr, t_obs, SF6_obs, rom, a_obs_init_method="quadratic_lag")
#%%
%%timeit
age_from_conv_old(SF6_ref_t, SF6_ref_vmr, t_obs, SF6_obs, rom)
#%%
%%timeit
Conc2Age_Convolution(SF6_ref_t, SF6_ref_vmr, t_obs, SF6_obs, rom, res="G")
#%%
%%timeit
Conc2Age_Convolution(SF6_ref_t, SF6_ref_vmr, t_obs, SF6_obs, rom, res="c")
# %%
