"""Initial conditions for the density-jump test of List et al. 2019, section 5.1.

Two dimensional, non-cosmological. The gas sits on a regular grid with a factor 1e4 jump in cell
mass across x = L/2, cold enough everywhere that nothing moves on its own. The dark matter sits on
a finer regular grid straddling that jump, with masses following a two-dimensional Gaussian, so the
mass is concentrated in a few particles at the discontinuity while every particle still has gas in
its kernel. Positions are uniform and masses are not, which is what makes the test hard for the
receiver method: it has to reconstruct a sharply peaked dark matter density at the gas positions.

Units are Gadget's: kpc, 1e10 Msun, km/s, so one code time unit is 978.5 Myr, the paper's end time.

    python3 scripts/test_problems/outputs/density_jump/make_ics.py
"""

import numpy as np
import h5py

BOXSIZE = 40.0           # kpc; the blast has to stay inside it, see the README
N_GAS = 256              # gas particles per side, as in the paper
N_DM = 513               # dark matter particles per side, as in the paper
SIGMA_DM = 3.0           # kpc, the width of the dark matter distribution
MASS_GAS_LOW = 1.5e-3    # 1.5e7 Msun, the left half
MASS_GAS_HIGH = 1.5e1    # 1.5e11 Msun, the right half
MASS_DM_TOTAL = 1.0      # 1e10 Msun over all dark matter particles
U_INIT = 0.01            # (km/s)^2, a sound speed of 0.1 km/s: 0.1 kpc of drift over the whole run

OUTPUT = "scripts/test_problems/ics/density_jump_ics.hdf5"


def grid(n):
    """cell-centred positions along one side, so nothing lands on the discontinuity or the boundary"""
    return (np.arange(n) + 0.5) * BOXSIZE / n


def make_ics():
    x_g, y_g = [a.ravel() for a in np.meshgrid(grid(N_GAS), grid(N_GAS), indexing="ij")]
    n_gas = x_g.size
    m_g = np.where(x_g < 0.5 * BOXSIZE, MASS_GAS_LOW, MASS_GAS_HIGH)
    u_g = np.full(n_gas, U_INIT)

    x_d, y_d = [a.ravel() for a in np.meshgrid(grid(N_DM), grid(N_DM), indexing="ij")]
    n_dm = x_d.size
    r2 = (x_d - 0.5 * BOXSIZE) ** 2 + (y_d - 0.5 * BOXSIZE) ** 2
    m_d = np.exp(-0.5 * r2 / SIGMA_DM ** 2)
    m_d *= MASS_DM_TOTAL / m_d.sum()

    with h5py.File(OUTPUT, "w") as f:
        npart = np.array([n_gas, n_dm, 0, 0, 0, 0])
        h = f.create_group("Header")
        h.attrs["NumPart_ThisFile"] = npart.astype(np.int32)
        h.attrs["NumPart_Total"] = npart.astype(np.uint32)
        h.attrs["NumPart_Total_HighWord"] = np.zeros(6, dtype=np.uint32)
        h.attrs["MassTable"] = np.zeros(6)
        h.attrs["BoxSize"] = BOXSIZE
        h.attrs["Time"] = 0.0
        h.attrs["Redshift"] = 0.0
        h.attrs["NumFilesPerSnapshot"] = np.int32(1)
        h.attrs["Flag_DoublePrecision"] = np.int32(0)

        p = f.create_group("PartType0")
        p.create_dataset("Coordinates", data=np.column_stack([x_g, y_g, np.zeros(n_gas)]).astype(np.float32))
        p.create_dataset("Velocities", data=np.zeros((n_gas, 3), dtype=np.float32))
        p.create_dataset("ParticleIDs", data=np.arange(1, n_gas + 1, dtype=np.uint32))
        p.create_dataset("Masses", data=m_g.astype(np.float32))
        p.create_dataset("InternalEnergy", data=u_g.astype(np.float32))

        p = f.create_group("PartType1")
        p.create_dataset("Coordinates", data=np.column_stack([x_d, y_d, np.zeros(n_dm)]).astype(np.float32))
        p.create_dataset("Velocities", data=np.zeros((n_dm, 3), dtype=np.float32))
        p.create_dataset("ParticleIDs", data=np.arange(n_gas + 1, n_gas + n_dm + 1, dtype=np.uint32))
        p.create_dataset("Masses", data=m_d.astype(np.float32))

    print("wrote %s: %d gas, %d dark matter" % (OUTPUT, n_gas, n_dm))
    print("gas mass %.4g, dark matter mass %.4g, smallest dark matter mass %.3g [1e10 Msun]"
          % (m_g.sum(), m_d.sum(), m_d.min()))


if __name__ == "__main__":
    make_ics()
