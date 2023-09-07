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


wt = np.where(
    np.logical_and((t_ref < (t_obs[0])), (t_ref > (t_obs[0] - t_int)))
)[
    0
]  # why is here a "0"?
