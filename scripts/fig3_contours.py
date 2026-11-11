"""Fig. 3: cowboy-hat dispersion of the 20 nm film and trireflection.

(a) f(kx, ky) of the lowest band at 250 mT from the COMSOL angular sweep.
(b) 9 GHz isofrequency contour (COMSOL, dashed) on top of the Fourier
    transform of a MuMax3 beam simulation at the same frequency.
(c) The simulated beam |m(x, y)|: one incident beam, several reflected ones.

(b) and (c) need the data release (python scripts/download_data.py).

    python scripts/fig3_contours.py
"""

from pathlib import Path

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np

from pma_spinwaves import comsol
from pma_spinwaves.fourier import isofrequency_map

ROOT = Path(__file__).resolve().parents[1]
BEAM = ROOT / "data" / "release" / "thick_20nm_beam3_ABL_65.5deg_130rad_9GHz_v2.npz"
F0 = 9.0  # GHz


def main():
    # The COMSOL sweep measures the direction of k from the field. In the
    # MuMax3 frame the field is along y, so the COMSOL "kx" is plotted as ky.
    k_par, k_perp, f = comsol.load_angular_sweep(
        comsol.DATA / "angular" / "thick_20nm_B250mT_angles_0_to_pi2.txt")
    kx, ky = k_perp, k_par

    fig = plt.figure(figsize=(12, 4))
    ax_a = fig.add_subplot(1, 3, 1, projection="3d")
    ax_a.plot_trisurf(kx, ky, f, cmap="viridis", linewidth=0, alpha=0.7)
    ax_a.tricontour(kx, ky, f, zdir="z", offset=0, levels=[7.9, 8.4, 10, 15, 20, 30],
                    cmap="gist_stern", linewidths=0.6)
    ax_a.tricontour(kx, ky, f, zdir="z", offset=0, levels=[F0], colors="m", linewidths=1.5)
    ax_a.set_zlim(0, 25)
    ax_a.set_xlabel(r"$k_x$ (rad/$\mu$m)", labelpad=-4)
    ax_a.set_ylabel(r"$k_y$ (rad/$\mu$m)", labelpad=-4)
    ax_a.set_zlabel(r"$f$ (GHz)", labelpad=-4)
    ax_a.view_init(35, 75)
    ax_a.set_title("(a) COMSOL, 20 nm, 250 mT", fontsize=9)

    ax_b = fig.add_subplot(1, 3, 2)
    ax_c = fig.add_subplot(1, 3, 3)
    if BEAM.exists():
        d = np.load(BEAM)
        _, _, amp, extent = isofrequency_map(d["Mf0"], float(d["cx"]), float(d["cy"]))
        ax_b.imshow(amp, origin="lower", extent=np.array(extent) * 1e-6, cmap="inferno",
                    norm=mcolors.Normalize(0, amp.max() / 20))
        ax_c.imshow(np.abs(d["Mf0"]), origin="lower", cmap="inferno", aspect="auto")
    else:
        for ax in (ax_b, ax_c):
            ax.text(0.5, 0.5, "run scripts/download_data.py", ha="center", va="center",
                    transform=ax.transAxes, fontsize=8)
    ax_b.tricontour(kx, ky, f, levels=[F0], colors="cyan", linewidths=0.8, linestyles="--")
    ax_b.set_xlim(-150, 150)
    ax_b.set_ylim(-150, 150)
    ax_b.set_aspect("equal")
    ax_b.set_xlabel(r"$k_x$ (rad/$\mu$m)")
    ax_b.set_ylabel(r"$k_y$ (rad/$\mu$m)")
    ax_b.set_title(f"(b) {F0:.0f} GHz: MuMax3 FFT and COMSOL contour", fontsize=9)
    ax_c.set_xticks([])
    ax_c.set_yticks([])
    ax_c.set_title("(c) MuMax3 beam, 20 nm film", fontsize=9)

    fig.tight_layout()
    out = ROOT / "figures" / "fig3_contours.png"
    fig.savefig(out, dpi=200)
    print(f"saved {out}")


if __name__ == "__main__":
    main()
