import numpy as np
import pandas as pd
from pathlib import Path

from run_convolution_method import *

fdir = Path('tropical_time_series')
SF6_ref_file = fdir / "SF6_tropical_surface_time_series_Ray.txt"

SF6_ref = pd.read_csv(SF6_ref_file, delim_whitespace=True, header=5, names=["time", "vmr"])

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
