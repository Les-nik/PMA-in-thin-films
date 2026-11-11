"""Fig. 1: spin waves in a 2 nm CoFeB film with perpendicular anisotropy.

Panels (a)-(d) come from the analytical model in ``pma_spinwaves.ultrathin``.
Panels (e) and (f) show MuMax3 results and need the data release
(``python scripts/download_data.py``); without it they are left empty.

    python scripts/fig1_ultrathin.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from pma_spinwaves.ultrathin import Film, dispersion_minimum, frequency_surface, omega, omega_components

ROOT = Path(__file__).resolve().parents[1]
RELEASE = ROOT / "data" / "release"
FILM = Film()  # 2 nm CoFeB: Ms = 1344 kA/m, A = 13.6 pJ/m
GHZ = 1e-9


def components_panel(ax, B0, Q, k):
    w_oop, w_ip = omega_components(k, B0, Q, np.pi / 2, FILM)
    ax.plot(k * 1e-6, w_oop * GHZ, "r--", label=r"$\omega_\perp$")
    ax.plot(k * 1e-6, w_ip * GHZ, "g-.", label=r"$\omega_\parallel$")
    ax.plot(k * 1e-6, omega(k, B0, Q, np.pi / 2, FILM) * GHZ, "k-", lw=2, label=r"$\omega$")
    ax.axhline(0, color="k", lw=0.5)
    ax.set_title(f"Q = {Q}, $\\mu_0H_0$ = {B0 * 1e3:.0f} mT", fontsize=9)
    ax.set_xlabel(r"$k$ (rad/$\mu$m)")
    ax.set_ylabel(r"$\omega$ (rad GHz)")
    ax.legend(frameon=False, fontsize=8, ncol=3)


def main():
    k = np.linspace(0, 120e6, 1000)
    fig = plt.figure(figsize=(11, 6.5))

    ax_a = fig.add_subplot(2, 3, 1)
    components_panel(ax_a, 0.4, 0.0, k)
    ax_a.set_ylim(0, 550)

    ax_b = fig.add_subplot(2, 3, 4)
    components_panel(ax_b, 0.4, 1.2, k)
    ax_b.set_ylim(-20, 200)

    # (c) backward-volume and Damon-Eshbach branches near the softening field
    ax_c = fig.add_subplot(2, 3, 2)
    for B0, ls in [(0.4, "-"), (0.3713, "--")]:
        for phi, color, name in [(0.0, "r", "BV"), (np.pi / 2, "k", "DE")]:
            ax_c.plot(k * 1e-6, omega(k, B0, 1.2, phi, FILM) * GHZ, color=color, ls=ls,
                      label=f"{name}, {B0 * 1e3:.1f} mT")
        kmin, fmin = dispersion_minimum(B0, 1.2, np.pi / 2, FILM)
        ax_c.plot(kmin * 1e-6, 2 * np.pi * fmin * GHZ, "o", color="tab:blue", ms=4)
    ax_c.set_xlim(0, 100)
    ax_c.set_ylim(0, 40)
    ax_c.set_xlabel(r"$k$ (rad/$\mu$m)")
    ax_c.set_ylabel(r"$\omega$ (rad GHz)")
    ax_c.set_title("Q = 1.2: sombrero-like minimum", fontsize=9)
    ax_c.legend(frameon=False, fontsize=7)

    # (d) frequency surface over (kx, ky), field along y
    ax_d = fig.add_subplot(2, 3, 5, projection="3d")
    kx, ky, f = frequency_surface(80e6, 0.4, 1.2, n=161, film=FILM)
    ax_d.plot_surface(kx * 1e-6, ky * 1e-6, f * GHZ, cmap="viridis", alpha=0.6,
                      rstride=4, cstride=4, linewidth=0)
    ax_d.contour(kx * 1e-6, ky * 1e-6, f * GHZ, zdir="z", offset=np.nanmin(f) * GHZ,
                 levels=8, cmap="viridis", linewidths=0.7)
    ax_d.set_xlabel(r"$k_x$ (rad/$\mu$m)", labelpad=-2)
    ax_d.set_ylabel(r"$k_y$ (rad/$\mu$m)", labelpad=-2)
    ax_d.set_zlabel(r"$f$ (GHz)", labelpad=-2)
    ax_d.view_init(elev=30, azim=45)
    ax_d.set_title(r"$f(k_x, k_y)$, $\mathbf{H}_0 \parallel y$", fontsize=9)

    # (e), (f): MuMax3 results from the data release
    ax_e = fig.add_subplot(2, 3, 3)
    ax_f = fig.add_subplot(2, 3, 6)
    ifc = RELEASE / "IFC_1500MHz.npz"
    beam = RELEASE / "bireflection_1500MHz_25deg.npz"
    if ifc.exists():
        d = np.load(ifc)
        ext = np.array([-d["kxMax"], d["kxMax"], -d["kyMax"], d["kyMax"]]) * 1e-6
        # IFC_1500MHz.npz was saved already in display orientation, so it is
        # shown with origin="upper" exactly as in the original notebook.
        ax_e.imshow(d["IFC"], extent=ext, cmap="inferno", origin="upper")
        ax_e.set_xlim(-80, 80)
        ax_e.set_ylim(-80, 80)
        ax_e.set_xlabel(r"$k_x$ (rad/$\mu$m)")
        ax_e.set_ylabel(r"$k_y$ (rad/$\mu$m)")
        ax_e.set_title("MuMax3: isofrequency contour, 1.5 GHz", fontsize=9)
    if beam.exists():
        d = np.load(beam)
        ax_f.imshow(np.abs(d["Mf0"]), origin="lower", cmap="inferno", aspect="auto")
        ax_f.axhline(300, color="w", lw=0.5)
        ax_f.text(50, 280, "Q = 0", color="w", va="top", fontsize=8)
        ax_f.text(50, 320, "Q = 1.2", color="w", va="bottom", fontsize=8)
        ax_f.set_title("MuMax3: beam reflected at a Q interface", fontsize=9)
    for ax in (ax_e, ax_f):
        if not ax.images:
            ax.text(0.5, 0.5, "run scripts/download_data.py", ha="center", va="center",
                    transform=ax.transAxes, fontsize=8)
        if ax is ax_f:
            ax.set_xticks([])
            ax.set_yticks([])

    for ax, label in zip([ax_a, ax_c, ax_e, ax_b, ax_f], "acebf"):
        color = "w" if ax.images else "k"
        ax.text(0.02, 0.92, f"({label})", transform=ax.transAxes, color=color)
    ax_d.text2D(0.02, 0.92, "(d)", transform=ax_d.transAxes)

    fig.tight_layout()
    out = ROOT / "figures" / "fig1_ultrathin.png"
    fig.savefig(out, dpi=200)
    print(f"saved {out}")


if __name__ == "__main__":
    main()
