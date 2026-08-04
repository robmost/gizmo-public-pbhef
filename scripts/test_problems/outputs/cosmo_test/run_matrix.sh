#!/bin/bash
# Both hydro solvers crossed with no PBHEF, the receiver method, and the donor method under each of
# its two weightings. Run from the repository root:
#
#     scripts/test_problems/outputs/cosmo_test/run_matrix.sh [ntasks] [solvers]
#
# solvers defaults to both. Pass 'mfm' or 'mfv' to run only one.
#
# The donor budget fraction is printed for each donor run and should read one. Both solvers reach z=49
# in 129 steps.
#
# Everything this generates lives under runs/<name>/, one directory per configuration, holding the
# config, the parameter file, the log and the snapshot output. Nothing is written beside this script.
#
# GIZMO_config.h is deleted before each build because config-makefile only regenerates it when the
# config file is newer than it, so a plain 'make CONFIG=...' after a different build silently keeps
# the previous configuration and reports the binary as up to date.
set -u
NTASK=${1:-2}
SOLVERS=${2:-"mfm mfv"}
DIR=scripts/test_problems/outputs/cosmo_test
EXE=GIZMO_cosmo_test

BASE="BOX_PERIODIC
PMGRID=32
ADAPTIVE_GRAVSOFT_FORGAS
EOS_GAMMA=(5.0/3.0)
USE_FFTW3
DOUBLEPRECISION_FFTW"

names=(mfm_base mfm_pbhef1 mfm_pbhef2_mass mfm_pbhef2_iso
       mfv_base mfv_pbhef1 mfv_pbhef2_mass mfv_pbhef2_iso)
flags=("HYDRO_MESHLESS_FINITE_MASS"
       "HYDRO_MESHLESS_FINITE_MASS PBHEF=1"
       "HYDRO_MESHLESS_FINITE_MASS PBHEF=2 PBHEF_MASS_WEIGHTS"
       "HYDRO_MESHLESS_FINITE_MASS PBHEF=2"
       "HYDRO_MESHLESS_FINITE_VOLUME"
       "HYDRO_MESHLESS_FINITE_VOLUME PBHEF=1"
       "HYDRO_MESHLESS_FINITE_VOLUME PBHEF=2 PBHEF_MASS_WEIGHTS"
       "HYDRO_MESHLESS_FINITE_VOLUME PBHEF=2")

printf "%-18s %-5s %-7s %-9s %-9s %s\n" config rc steps finished nan "budget fraction"
for i in "${!names[@]}"; do
    name=${names[$i]}
    case " $SOLVERS " in *" ${name%%_*} "*) ;; *) continue ;; esac
    run=$DIR/runs/$name
    rm -rf "$run"; mkdir -p "$run"

    { echo "$BASE"; for f in ${flags[$i]}; do echo "$f"; done; } > "$run/Config.sh"
    sed "s|^OutputDir .*|OutputDir  $run/output|" $DIR/cosmo_test.params > "$run/run.params"

    rm -f GIZMO_config.h compile_time_info.c
    if ! make CONFIG="$run/Config.sh" EXEC=$EXE -j8 > "$run/build.log" 2>&1; then
        printf "%-18s %s\n" "$name" "BUILD FAILED, see $run/build.log"; continue
    fi

    timeout 3600 mpirun -np "$NTASK" ./$EXE "$run/run.params" > "$run/run.log" 2>&1
    rc=$?
    steps=$(grep -c '^Step' "$run/output/cpu.txt" 2>/dev/null || echo 0)
    fin=$(grep -ac 'Simulation ends' "$run/run.log")
    nan=$(grep -ac 'nan' "$run/run.log")
    frac=$(grep -a 'budget' "$run/run.log" | tail -1 | sed 's/.*fraction=//;s/)//')
    printf "%-18s %-5s %-7s %-9s %-9s %s\n" "$name" "$rc" "$steps" "$fin" "$nan" "${frac:--}"
done
rm -f $EXE GIZMO_config.h compile_time_info.c
