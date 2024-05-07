# AoA_from_convolution
 tools to calculate the mean age of air from mixing ratios relative to a reference time series

## How to use the code
0. Make sure to have the mandatory packages installed in your python environment:
```
pandas
numpy
scipy
matplotlib
```
the code has been tested with the following package versions, but should work for a lot of other versions, too: 
```
pandas=2.1.1
numpy=1.26.2
scipy=1.11.3
matplotlib=3.8.0
```

2. Download the package
3. E.g. create an empty python file (e.g. `my_project.py`) within your downloaded package
4. within `my_project.py` you can then use the `SF6_to_AoA()` and the `CO2_to_AoA()` functions to derive mean age from observations:
```
# import everything from the run_convolution_method.py file into your current python session:

from run_convolution_method import *


# Example to calculate mean age from SF6 mixing ratio values

# Specify your observation date:
t_obs_sf6 = 2022.3  # float value in fractional year, single value only

# have your SF6 mixing ratios ready, e.g.:
sf6mr = np.arange(9.5, 11, 0.1)  # can be single value or array

# decide on what value for the ratio of moments (rom) to use. I.e. the ratio of moments of the age spectrum 
# see e.g. Garny et al. (2024), Reviews of Geophysics
ratio_of_moments = 1.2  # single value only

# calculate mean age from SF6 mixing ratios:
mean_age = SF6_to_AoA(t_obs=t_obs_sf6, SF6_obs=sf6mr, rom=ratio_of_moments)

print(mean_age)

# obviously, if you want to calculate the mean age from variable observation dates, you need to 
# write some kind of for loop

# %%
# Example to calculate mean age from CO2 mixing ratio values (and CH4 values)

# Specify your observation date:
t_obs_co2 = 2022.3  # float value in fractional year, single value only

# have your CO2 and CH4 mixing ratios ready, e.g.:
co2mr = np.arange(400, 420, 0.5)  # can be single value or array
ch4mr = np.ones(co2mr.shape) * 1100.3  # must be the same size as co2mr

# decide on what value for the ratio of moments (rom) to use. I.e. the ratio of moments of the age spectrum 
# see e.g. Garny et al. (2024), Reviews of Geophysics
ratio_of_moments = 1.2  # single value only

# calculate mean age from SF6 mixing ratios:
mean_age = CO2_to_AoA(t_obs=t_obs_co2, CO2_obs=co2mr, rom=ratio_of_moments, CH4_obs=ch4mr)

print(mean_age)

# WARNING: due to the propagation of the seasonal cycle in CO2 mixing ratios into the stratosphere,
# only mean age values above a threshold of ca. 2.5 years are expected to be physically meaningful.

# obviously, if you want to calculate the mean age from variable observation dates, you need to 
# write some kind of for loop

```
