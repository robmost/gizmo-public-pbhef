#!/bin/bash
# Both hydro solvers crossed with no PBHEF, the receiver method, and the donor method under each of
# its two weightings. Run from the repository root:
#
#     scripts/test_problems/outputs/density_jump/run_matrix.sh [ntasks] [solvers]
#
# solvers defaults to both. Pass 'mfm' or 'mfv' to run only one. The paper uses MFM for this test.
#
# What to compare is the ratio of the energy the gas actually received to the analytic total the
# black holes radiated, which is the last column. The donor method should read one; the receiver
# method is the number the test exists to measure. See the README.
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
DIR=scripts/test_problems/outputs/density_jump
EXE=GIZMO_density_jump

BASE="BOX_PERIODIC
BOX_SPATIAL_DIMENSION=2
SELFGRAVITY_OFF
EOS_GAMMA=(5.0/3.0)"

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

# injected/expected at the last sync point the donors deposited on; a row where they did not lags by
# up to one step, and a receiver-mode file has no such column so every row of it qualifies
ratio_of() {
    awk 'NR>1 && $1>0 && (NF<9 || $9>0.5) {r=$6/$7} END {if(r!="") {printf "%.9f", r} else {printf "-"}}' "$1"
}

printf "%-18s %-5s %-7s %-9s %-5s %s\n" config rc steps finished nan injected/expected
for i in "${!names[@]}"; do
    name=${names[$i]}
    case " $SOLVERS " in *" ${name%%_*} "*) ;; *) continue ;; esac
    run=$DIR/runs/$name
    rm -rf "$run"; mkdir -p "$run"

    { echo "$BASE"; for f in ${flags[$i]}; do echo "$f"; done; } > "$run/Config.sh"
    sed "s|^OutputDir .*|OutputDir  $run/output|" $DIR/density_jump.params > "$run/run.params"

    rm -f GIZMO_config.h compile_time_info.c
    if ! make CONFIG="$run/Config.sh" EXEC=$EXE -j8 > "$run/build.log" 2>&1; then
        printf "%-18s %s\n" "$name" "BUILD FAILED, see $run/build.log"; continue
    fi

    timeout 7200 mpirun -np "$NTASK" ./$EXE "$run/run.params" > "$run/run.log" 2>&1
    rc=$?
    rm -rf "$run/output/restartfiles"   # 170 MB a run, and nothing here resumes
    steps=$(grep -c '^Step' "$run/output/cpu.txt" 2>/dev/null || echo 0)
    fin=$(grep -ac 'Simulation ends' "$run/run.log")
    nan=$(grep -ac 'nan' "$run/run.log")
    if [ -f "$run/output/pbhef_energy.txt" ]; then r=$(ratio_of "$run/output/pbhef_energy.txt"); else r="-"; fi
    printf "%-18s %-5s %-7s %-9s %-5s %s\n" "$name" "$rc" "$steps" "$fin" "$nan" "$r"
done
rm -f $EXE GIZMO_config.h compile_time_info.c
