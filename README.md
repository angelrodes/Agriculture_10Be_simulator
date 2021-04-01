# Agriculture_10Be_simulator
This MATLAB code generates a video to simulate the Be-10 accumulation in arable soils

## Model evolution

This model simulates the accumulation of Be-10 under a surface that is being eroded at a constant rate of 50 mm/ka:

* from t=0 to 50 ka, the density of the profile is constant
* from t=50 to 100 ka, the density of the profile changes to simulate soil in the upper 1.5m and saprolite until a depth of 3m. The profile is still in equilibrium: the soil production rate (SPR) is the same as the erosion rate at the top of the profile (both 50 mm/ka)
* from t=100 to 150 ka, the top 20 cm of the profile is affected by bioturbation, generating a mixing layer enriched in Be-10
* from t=150 to 200 ka, sustainable agriculture is simulated by a mixing layer of 50 cm
* from t=200 to 201 ka (last 1000 years), unsustainable agriculture is simulated by eroding the surface at a rate faster than the soil production rate (the profile is not in equilibrium!)

## Unsustainable agriculture

The effect of unsustainable agriculture is like a bulldozer that removes the top of the soil and mixes the first 0.5 m of the remaining soil.

### Effect on Be-10-derived soil production rates

At the end of the simulation, I have calculated the apparent soil production rate (a-SPR) from a sample taken in the saprolite (red star).

The a-SPR is faster than the actual SPR (110 mm/ka vs. 50 mm/ka), but much slower than the erosion rate generated by unsustainable agriculture (~3000 mm/ka).

Also, a red line illustrates the Be-10 profile we would reconstruct from that sample usign [CoSOILcal](https://github.com/angelrodes/cosoilcal), which matches well the actual Be-10 profile. This means that **even taking several samples from this profile, we probably wouldn't notice the break in the equilibrium generated by unsustainable agriculture.**

## Output

The script generates [this video](https://youtu.be/V5eV1DUaHhw):

[![Agriculture_10Be_simulator](https://img.youtube.com/vi/V5eV1DUaHhw/0.jpg)](https://www.youtube.com/watch?v=V5eV1DUaHhw)

## More info

There is a post illusrated by this model in [my blog](https://angelrodes.wordpress.com/2021/03/31/agriculture-be-10-simulator/).

## Licence

Feel free to use this code. You can cite [the CoSOILcal paper](https://doi.org/10.1016/j.mex.2019.11.026) as a reference for academic publications.
