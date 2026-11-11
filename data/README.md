# Data

## `comsol/` (in the repository)

Exports from the COMSOL Multiphysics model of a 20 nm CoFeB film
(Ms = 1344 kA/m, A = 13.6 pJ/m, Q = 0.6) and of a 2 nm film. Every file has a
COMSOL header (lines starting with `%`) and two columns. Read them with
`pma_spinwaves.comsol`.

| Folder | Columns | Content |
| --- | --- | --- |
| `thick_20nm/dispersion/dispersion_B{250,500,800}mT_Q0.6_3bands.txt` | k (rad/µm), f (GHz) | first three bands; the three rows for each k are bands 1, 2, 3 |
| `thick_20nm/dispersion/dispersion_BEC.txt` | k (rad/µm), f (GHz) | lowest band near its minimum, finer k step |
| `thick_20nm/profiles/{mx,mz}[_imag]_B{B}mT_k{k}.txt` | z (nm), value | real and imaginary parts of the dynamic magnetization across the thickness (z = 0 is the film centre) |
| `thick_20nm/torques/{tx,tz}_{exch,pma,zeeman,demag,eff}_B{B}_k{k}.txt` | z (nm), value | torque contributions across the thickness |
| `angular/thick_20nm_B250mT_angles_0_to_pi2.txt` | k (rad/µm), f (GHz) | lowest band for 17 directions of k, 0 to pi/2 in steps of pi/32 |
| `angular/ultrathin_2nm_B250mT_angles_0_to_pi2.txt` | k (rad/µm), f (GHz) | the same for the 2 nm film |

## `release/` (downloaded with `scripts/download_data.py`)

Single-frequency complex maps `Mf0` (shape ny x nx, cell sizes `cx`, `cy`
in metres) from MuMax3 runs, and one precomputed isofrequency map:

| File | Used in |
| --- | --- |
| `IFC_1500MHz.npz` | Fig. 1(e), 2 nm film, 1.5 GHz |
| `bireflection_1500MHz_25deg.npz` | Fig. 1(f), beam at 25 deg on a Q = 0 / Q = 1.2 interface |
| `thick_20nm_beam3_ABL_65.5deg_130rad_9GHz_v2.npz` | Fig. 3(b, c), 20 nm film, 9 GHz |

The data are released under CC BY 4.0, like the MuMax3 scripts on Zenodo.
