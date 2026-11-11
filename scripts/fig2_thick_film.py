"""Fig. 2: anti-Larmor precession in a 20 nm CoFeB film (Q = 0.6), from COMSOL data.

(a) First two bands at 250, 500 and 800 mT. The stars mark the points whose
    torques are shown in (b)-(d).
(b)-(d) x component of the torque and its contributions across the thickness.
(e) Fraction of the thickness that precesses against the Larmor sense, for
    every k at which a profile was exported.
(f) Precession orbits at five depths for selected (B, k).

    python scripts/fig2_thick_film.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from pma_spinwaves import comsol, precession

ROOT = Path(__file__).resolve().parents[1]
DISP = comsol.DATA / "thick_20nm" / "dispersion"

FIELDS = {250: "darkgreen", 500: "red", 800: "orange"}
# Points shown in panels (b)-(d): (field in mT, k of the torque export, star position k, f)
# These are the choices of the published figure. The 500 mT star sits at
# k = 93 rad/um while the torque export used for it is the k = 86 one; both
# files are in data/comsol/thick_20nm/torques if you want to compare.
TORQUE_POINTS = [(800, 80, 80, 30.55), (500, 86, 93, 19.95), (250, 106, 106, 6.1)]
TERM_STYLE = {
    "exch": ("blue", "--", r"$\tau^\mathrm{ex}$"),
    "pma": ("red", "--", r"$\tau^\mathrm{PMA}$"),
    "zeeman": ("orange", "--", r"$\tau^{0}$"),
    "demag": ("green", "--", r"$\tau^\mathrm{d}$"),
    "eff": ("black", "-", r"$\tau^\mathrm{eff}$"),
}
PROFILE_K = {250: [0, 46, 50, 53, 106, 150, 153, 159, 166, 183, 186], 500: [86, 93, 106], 800: [80, 106]}


def main():
    fig = plt.figure(figsize=(11, 7))
    gs = fig.add_gridspec(3, 4, width_ratios=[1.3, 1, 1, 1.4])

    # (a) dispersion
    ax_a = fig.add_subplot(gs[:, 0])
    for B, color in FIELDS.items():
        k, f = comsol.load_bands(DISP / f"dispersion_B{B}mT_Q0.6_3bands.txt", n_bands=3)
        ax_a.plot(k, f[0], color=color, lw=2, label=f"{B} mT")
        ax_a.plot(k, f[1], color=color, lw=1.2, ls="--", alpha=0.7)
    for _, _, ks, fs in TORQUE_POINTS:
        ax_a.plot(ks, fs, "*", color="crimson", mec="k", ms=11, zorder=5)
    ax_a.set_xlim(0, 200)
    ax_a.set_ylim(0, 40)
    ax_a.set_xlabel(r"$k$ (rad/$\mu$m)")
    ax_a.set_ylabel("frequency (GHz)")
    ax_a.legend(frameon=False, title=r"$\mu_0 H_0$")

    # (b)-(d) torques, normalised to the PMA torque as in the paper
    for row, (B, k_t, _, _) in enumerate(TORQUE_POINTS):
        ax = fig.add_subplot(gs[row, 1])
        z, terms = comsol.load_torques(B, k_t)
        scale = -1.3 * terms["pma"].max()
        for term, (c, ls, lab) in TERM_STYLE.items():
            ax.plot(terms[term] / scale, z + 10, color=c, ls=ls, lw=2 if term == "eff" else 1.4, label=lab)
        ax.axvline(0, color="maroon", ls="-.", lw=1)
        ax.set_xlim(-1, 1)
        ax.set_ylim(0, 20)
        ax.set_ylabel("z (nm)")
        ax.set_title(f"{B} mT, k = {k_t} rad/µm", fontsize=9)
        if row == 0:
            ax.legend(frameon=False, fontsize=7, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.45))
        if row == 2:
            ax.set_xlabel(r"$\tau_x$ (arb. units)")

    # (e) anti-Larmor fraction of the thickness
    ax_e = fig.add_subplot(gs[0, 2:])
    for B, ks in PROFILE_K.items():
        frac = []
        for k in ks:
            _, mx = comsol.load_profile("mx", B, k)
            _, mz = comsol.load_profile("mz", B, k)
            s = precession.handedness(mx, mz)
            # k = 0 is always Larmor; its sign fixes the convention (here -1).
            frac.append(np.mean(s > 0))
        ax_e.plot(ks, frac, "o-", color=FIELDS[B], label=f"{B} mT")
    ax_e.set_xlabel(r"$k$ (rad/$\mu$m)")
    ax_e.set_ylabel("anti-Larmor\nfraction of d")
    ax_e.set_ylim(-0.03, 0.5)
    ax_e.legend(frameon=False, fontsize=8)

    # (f) orbits at five depths
    cases = [(250, 0), (250, 46), (250, 106), (250, 166), (250, 186), (800, 80)]
    sub = gs[1:, 2:].subgridspec(1, len(cases), wspace=0.05)
    for i, (B, k) in enumerate(cases):
        ax = fig.add_subplot(sub[0, i])
        z, mx = comsol.load_profile("mx", B, k)
        _, mz = comsol.load_profile("mz", B, k)
        norm = 1.1 * np.max(np.abs(np.concatenate([mx, mz])))
        # Rotate all depths by the phase of m_x at the strongest point so the
        # ellipses are drawn with a common time origin.
        ref = np.exp(-1j * np.angle(mx[np.argmax(np.abs(mx))]))
        for j, idx in enumerate(np.linspace(0, len(z) - 1, 5).astype(int)):
            x, y = precession.orbit(mx[idx] * ref / norm, mz[idx] * ref / norm)
            s = precession.handedness(mx[idx], mz[idx])
            color = "tab:red" if s > 0 else "tab:blue"
            ax.plot(x, y + 2 * j, color=color, lw=1.2)
            ax.axhline(2 * j, color="grey", lw=0.3)
        ax.set_xlim(-1.1, 1.1)
        ax.set_ylim(-1.2, 9.2)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks(range(0, 10, 2))
        ax.set_yticklabels(["0", "5", "10", "15", "20"] if i == 0 else [])
        if i == 0:
            ax.set_ylabel("z (nm)")
        ax.set_title(f"{B} mT\nk={k}", fontsize=8)
    fig.text(0.62, 0.02, "blue: Larmor, red: anti-Larmor", fontsize=8)

    for ax, label in zip([ax_a, ax_e], "ae"):
        ax.text(0.02, 0.96, f"({label})", transform=ax.transAxes, va="top")

    fig.tight_layout()
    out = ROOT / "figures" / "fig2_thick_film.png"
    fig.savefig(out, dpi=200)
    print(f"saved {out}")


if __name__ == "__main__":
    main()
