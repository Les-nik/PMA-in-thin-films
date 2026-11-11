# Simulations

## MuMax3

The three MuMax3 scripts behind Figs. 1 and 3 are archived on Zenodo:
[10.5281/zenodo.14871425](https://doi.org/10.5281/zenodo.14871425)
(`CoFeB_2nm_pointSource_1500MHz.mx3`, `bireflection.mx3`, `trireflection.mx3`).

`mumax3/damping_contours/` holds the scripts added during the review to check
how damping changes the isofrequency contours (Supporting Information):
a 20 nm film (Q = 0.6) at 226 and 800 mT and a 2 nm film (Q = 1.2) at 372.4
and 800 mT, each with alpha = 0.001 and a large damping (0.08 or 0.1). They
excite the film with a sinc pulse in space and time and save m(t) for a
dispersion calculation; the postprocessing is the same as in
[spin-wave-dispersion](https://github.com/Les-nik/spin-wave-dispersion).

## COMSOL

The thick-film results come from a custom COMSOL model (developed at Adam Mickiewicz
University) that solves the linearized Landau-Lifshitz equation for the thickness profile
m(z) of a spin wave with in-plane wavevector k. The `.mph` files are not included (about
1.4 GB); the exported results are in `data/comsol`.
