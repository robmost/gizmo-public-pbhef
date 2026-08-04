"""Compare the donor and receiver methods on the density-jump test.

Two things are measured. The first is the energy the gas received against the analytic total, which is
also what run_matrix.sh prints. The second is the instantaneous power each method hands to the gas,
which is what separates them: the donor shares out a total it knows exactly, while the receiver has to
rebuild the dark matter density from the gas positions and loses whatever that quadrature loses.
Writes one figure beside the runs. Needs numpy, h5py and matplotlib.

    python3 scripts/test_problems/outputs/density_jump/analyse.py
"""

import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RUNS = "scripts/test_problems/outputs/density_jump/runs"
FIGURE = RUNS + "/density_jump.png"
NAMES = ["mfm_pbhef1", "mfm_pbhef2_mass", "mfm_pbhef2_iso",
         "mfv_pbhef1", "mfv_pbhef2_mass", "mfv_pbhef2_iso"]
LABELS = {"pbhef1": "receiver", "pbhef2_mass": "donor, mass-weighted", "pbhef2_iso": "donor, solid angle"}
LAST_SNAPSHOT = 4
BOXSIZE = 40.0


def energy_log(name):
    """injected over analytic energy, at the last sync point the donors deposited on"""
    d = np.loadtxt("%s/%s/output/pbhef_energy.txt" % (RUNS, name))
    d = d[d[:, 0] > 0]
    if d.shape[1] > 8:
        d = d[d[:, 8] > 0.5]  # rows where the donors deposited, so nothing lags by a step
    return d[-1, 5] / d[-1, 6]


def snapshot(name, index):
    with h5py.File("%s/%s/output/snapshot_%03d.hdf5" % (RUNS, name, index), "r") as f:
        g = f["PartType0"]
        return (f["Header"].attrs["Time"], np.asarray(g["Coordinates"]), np.asarray(g["Masses"]),
                np.asarray(g["InternalEnergy"]), np.asarray(g["PBHEF_Dtu"]))


def report():
    print("%-18s %18s %12s %12s" % ("config", "injected/expected", "power t=0", "power t=end"))
    for name in NAMES:
        totals = []
        for index in (0, LAST_SNAPSHOT):
            _, _, mass, _, dtu = snapshot(name, index)
            totals.append((dtu * mass).sum())
        expected = 2013.7  # the prefactor times the total dark matter mass, both constant here
        print("%-18s %18.9f %12.6f %12.6f"
              % (name, energy_log(name), totals[0] / expected, totals[1] / expected))


KEYS = ["pbhef2_mass", "pbhef2_iso", "pbhef1"]


def figure():
    """internal energy the way the paper plots it, then the two profiles that separate the methods"""
    fig, axes = plt.subplots(2, 3, figsize=(13, 8.4), constrained_layout=True)
    images = []
    for column, key in enumerate(KEYS):
        time, pos, mass, u, dtu = snapshot("mfm_" + key, LAST_SNAPSHOT)
        ax = axes[0][column]
        images.append(ax.scatter(pos[:, 0], pos[:, 1], c=np.log10(u), s=1.2, vmin=-2, vmax=3,
                                 cmap="inferno", rasterized=True))
        ax.axvline(0.5 * BOXSIZE, color="c", lw=0.6, ls="--")
        ax.set_xlim(0, BOXSIZE)
        ax.set_ylim(0, BOXSIZE)
        ax.set_aspect("equal")
        ax.set_title("%s, t=%.2g" % (LABELS[key], time), fontsize=10)
        ax.set_xlabel("x [kpc]")
        if column == 0:
            ax.set_ylabel("y [kpc]")
    fig.colorbar(images[-1], ax=axes[0], shrink=0.8, label="log10 internal energy [(km/s)^2]")

    edges = np.linspace(0, 0.5 * BOXSIZE, 60)
    centres = 0.5 * (edges[1:] + edges[:-1])
    for key, style in zip(KEYS, ["-", "--", "-"]):  # donor curves overlap, so vary the line
        time, pos, mass, u, dtu = snapshot("mfm_" + key, LAST_SNAPSHOT)
        radius = np.hypot(pos[:, 0] - 0.5 * BOXSIZE, pos[:, 1] - 0.5 * BOXSIZE)
        power, _ = np.histogram(radius, bins=edges, weights=dtu * mass)
        axes[1][0].plot(centres, power / np.diff(edges), style, label=LABELS[key])
        axes[1][1].plot(centres, np.cumsum(power) / (dtu * mass).sum(), style, label=LABELS[key])

        totals = [snapshot("mfm_" + key, i) for i in range(LAST_SNAPSHOT + 1)]
        axes[1][2].plot([t[0] for t in totals], [(t[4] * t[2]).sum() / 2013.7 for t in totals],
                        style, marker="o", label=LABELS[key])

    axes[1][0].set_yscale("log")
    axes[1][0].set_xlabel("radius [kpc]")
    axes[1][0].set_ylabel("power per unit radius")
    axes[1][1].set_xlabel("radius [kpc]")
    axes[1][1].set_ylabel("cumulative fraction of power")
    axes[1][2].set_xlabel("time [978.5 Myr]")
    axes[1][2].set_ylabel("power coupled / analytic")
    for ax in axes[1]:
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)

    fig.savefig(FIGURE, dpi=130)
    print("\nwrote %s" % FIGURE)


if __name__ == "__main__":
    report()
    figure()
