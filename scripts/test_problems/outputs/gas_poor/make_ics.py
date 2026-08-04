"""Initial conditions for the gas-poor test: a peaked dark matter distribution in a gas cavity.

Two dimensional, non-cosmological. Gas fills a regular grid everywhere except a disc of radius
R_CAVITY at the centre, and the dark matter sits on a finer regular grid with masses following a
two-dimensional Gaussian centred on the same point. Increasing R_CAVITY moves the gas away from the
dark matter without changing either distribution, so it separates the two coupling methods: a donor
widens its kernel until it finds gas, while a receiver can only heat gas that already sits where the
dark matter density is high.

One file per radius in RADII. Units are kpc, 1e10 Msun, km/s, so one code time unit is 978.5 Myr.

    python3 scripts/test_problems/outputs/gas_poor/make_ics.py
"""

import numpy as np
import h5py

BOXSIZE = 40.0           # kpc
N_GAS = 256              # gas particles per side before the cavity is cut
N_DM = 513               # dark matter particles per side
SIGMA_DM = 3.0           # kpc, the width of the dark matter distribution
MASS_GAS = 1.5e-3        # 1.5e7 Msun a cell, uniform
MASS_DM_TOTAL = 1.0      # 1e10 Msun over all dark matter particles
U_INIT = 0.01            # (km/s)^2, a sound speed of 0.1 km/s
RADII = [0.0, 2.0, 4.0, 6.0, 9.0]

OUTPUT = "scripts/test_problems/ics/gas_poor_r%g_ics.hdf5"


def grid(n):
    """cell-centred positions along one side, so nothing lands on the centre or the boundary"""
    return (np.arange(n) + 0.5) * BOXSIZE / n


def lattice(n):
    x, y = [a.ravel() for a in np.meshgrid(grid(n), grid(n), indexing="ij")]
    return x, y, np.hypot(x - 0.5 * BOXSIZE, y - 0.5 * BOXSIZE)


def write(radius, x_g, y_g, x_d, y_d, m_d):
    n_gas, n_dm = x_g.size, x_d.size
    with h5py.File(OUTPUT % radius, "w") as f:
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
        p.create_dataset("Masses", data=np.full(n_gas, MASS_GAS, dtype=np.float32))
        p.create_dataset("InternalEnergy", data=np.full(n_gas, U_INIT, dtype=np.float32))

        p = f.create_group("PartType1")
        p.create_dataset("Coordinates", data=np.column_stack([x_d, y_d, np.zeros(n_dm)]).astype(np.float32))
        p.create_dataset("Velocities", data=np.zeros((n_dm, 3), dtype=np.float32))
        p.create_dataset("ParticleIDs", data=np.arange(n_gas + 1, n_gas + n_dm + 1, dtype=np.uint32))
        p.create_dataset("Masses", data=m_d.astype(np.float32))
    return n_gas


def make_ics():
    x_d, y_d, r_d = lattice(N_DM)
    m_d = np.exp(-0.5 * r_d ** 2 / SIGMA_DM ** 2)
    m_d *= MASS_DM_TOTAL / m_d.sum()

    x_g, y_g, r_g = lattice(N_GAS)
    for radius in RADII:
        keep = r_g >= radius
        n_gas = write(radius, x_g[keep], y_g[keep], x_d, y_d, m_d)
        inside = (m_d[r_d < radius].sum() / m_d.sum()) if radius > 0 else 0.0
        print("radius %4.1f kpc: %6d gas cells, %.4f of the dark matter mass has no gas over it"
              % (radius, n_gas, inside))


if __name__ == "__main__":
    make_ics()
