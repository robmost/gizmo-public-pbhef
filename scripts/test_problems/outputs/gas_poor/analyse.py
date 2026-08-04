"""Compare the donor and receiver methods as the gas is moved away from the dark matter.

The measured quantity is the fraction of the analytic evaporation power that reaches the gas. For the
receiver method that fraction is set by how much of the dark matter has gas sitting over it, so it
falls with the cavity radius. For the donor method it is one at every radius the kernel can span, and
falls only once the ceiling on that kernel is below the radius. Writes one figure beside the runs.

    python3 scripts/test_problems/outputs/gas_poor/analyse.py
"""

import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RUNS = "scripts/test_problems/outputs/gas_poor/runs"
FIGURE = RUNS + "/gas_poor.png"
RADII = [0, 2, 4, 6, 9]
LAST_SNAPSHOT = 4
BOXSIZE = 40.0
EXPECTED = 2013.7        # prefactor times the total dark matter mass, both constant here


def energy_log(name):
    """(whole-run injected/expected, lowest per-step coupled fraction, donors starved there)"""
    data = np.loadtxt("%s/%s/output/pbhef_energy.txt" % (RUNS, name))
    data = data[data[:, 0] > 0]
    ratio = data[-1, 5] / data[-1, 6]
    if data.shape[1] <= 8:
        return ratio, None, None
    live = data[data[:, 9] > 0]
    worst = live[np.argmin(live[:, 8])]
    return ratio, worst[8], int(live[:, 9].max() - worst[9])


def coupled(name, index):
    """instantaneous power the gas is being given, as a fraction of the analytic total"""
    with h5py.File("%s/%s/output/snapshot_%03d.hdf5" % (RUNS, name, index), "r") as f:
        g = f["PartType0"]
        return float((np.asarray(g["PBHEF_Dtu"]) * np.asarray(g["Masses"])).sum()) / EXPECTED


def report():
    print("%-22s %11s %11s %11s %s" % ("config", "power t=0", "power t=end", "whole-run", "worst-coupled"))
    for name in ["r%d_pbhef1" % r for r in RADII] + ["r%d_pbhef2" % r for r in RADII] \
            + ["r9_pbhef2_mass", "r9_pbhef1_shortreach", "r9_pbhef2_shortreach"]:
        ratio, worst, starved = energy_log(name)
        tail = "-" if worst is None else "%.6f (%d donors starved)" % (worst, starved)
        print("%-22s %11.6f %11.6f %11.6f %s"
              % (name, coupled(name, 0), coupled(name, LAST_SNAPSHOT), ratio, tail))


def figure():
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.6), constrained_layout=True)

    for method, label in [("pbhef1", "receiver"), ("pbhef2", "donor")]:
        axes[0].plot(RADII, [coupled("r%d_%s" % (r, method), 0) for r in RADII],
                     marker="o", label=label)
    axes[0].set_xlabel("cavity radius [kpc]")
    axes[0].set_ylabel("power coupled / analytic, at t=0")
    axes[0].set_ylim(-0.05, 1.1)

    for name, label in [("r9_pbhef2", "donor, 21.8 kpc kernel ceiling"),
                        ("r9_pbhef2_shortreach", "donor, 5.6 kpc ceiling"),
                        ("r9_pbhef1", "receiver")]:
        data = np.loadtxt("%s/%s/output/pbhef_energy.txt" % (RUNS, name))
        axes[1].plot(data[1:, 0], data[1:, 5] / data[1:, 6], label=label)
    axes[1].set_xlabel("time [978.5 Myr]")
    axes[1].set_ylabel("injected / analytic, cumulative")
    axes[1].set_title("9 kpc cavity", fontsize=10)

    for name, label in [("r6_pbhef2", "donor"), ("r6_pbhef1", "receiver")]:
        with h5py.File("%s/%s/output/snapshot_%03d.hdf5" % (RUNS, name, LAST_SNAPSHOT), "r") as f:
            pos = np.asarray(f["PartType0/Coordinates"])
            u = np.asarray(f["PartType0/InternalEnergy"])
        radius = np.hypot(pos[:, 0] - 0.5 * BOXSIZE, pos[:, 1] - 0.5 * BOXSIZE)
        edges = np.linspace(0, 0.5 * BOXSIZE, 50)
        total, _ = np.histogram(radius, bins=edges, weights=u)
        count, _ = np.histogram(radius, bins=edges)
        axes[2].plot(0.5 * (edges[1:] + edges[:-1]), total / np.maximum(count, 1), label=label)
    axes[2].set_yscale("log")
    axes[2].set_xlabel("radius [kpc]")
    axes[2].set_ylabel("mean internal energy [(km/s)^2]")
    axes[2].set_title("6 kpc cavity, t=0.25", fontsize=10)

    for ax in axes:
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
    fig.savefig(FIGURE, dpi=130)
    print("\nwrote %s" % FIGURE)


if __name__ == "__main__":
    report()
    figure()
