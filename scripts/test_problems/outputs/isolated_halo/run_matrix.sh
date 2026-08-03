#!/bin/bash
# The isolated disc galaxy with no PBHEF, the receiver method and the donor method, each without and
# with radiative cooling and star formation. Run from the repository root:
#
#     scripts/test_problems/outputs/isolated_halo/run_matrix.sh [ntasks] [name ...]
#
# Naming configurations runs only those, so a subset can be redone without discarding the rest.
#
# The runs with COOLING read the z=0 ultraviolet background from a TREECOOL file in the directory the
# code is run from. GIZMO's reader cannot skip comments, and says so when it fails to open the file, so
# a TREECOOL that still carries its header is parsed as zero entries and the run silently proceeds with
# no background at all. Check for "non-zero UVB entries" in the log: it should be 204, not 0.
#
# This is the only PBHEF problem with self-gravity and a live halo, so it is the one that tests two
# things the idealised problems cannot: donors that move appreciably between their own timesteps, and
# a dark matter distribution far more extended than the gas.
#
# Expect a coupled fraction well below one. The halo reaches 130 kpc while the gas disc is inside
# about 25 kpc, and a donor may only widen its kernel to 100 * 2.8 * SofteningHalo, which is 14 kpc
# here, so most of the halo has no gas to give its energy to. What matters is whether the two methods
# agree on the fraction they do couple.
#
# TimeMax is left at the value in the parameter file. A quarter of a code time unit is about 16
# minutes a run on four tasks; the whole unit is over an hour.
set -u
NTASK=${1:-4}
shift $(( $# > 0 ? 1 : 0 ))
WANTED="$*"
DIR=scripts/test_problems/outputs/isolated_halo
EXE=GIZMO_isolated_halo

BASE="HYDRO_MESHLESS_FINITE_MASS
ADAPTIVE_GRAVSOFT_FORGAS
EOS_GAMMA=(5.0/3.0)
MULTIPLEDOMAINS=16"

COOL="COOLING
GALSF
GALSF_EFFECTIVE_EQS
GALSF_SUBGRID_WINDS
GALSF_SUBGRID_WIND_SCALING=2"

names=(base pbhef1 pbhef2 cool_base cool_pbhef1 cool_pbhef2)
flags=("" "PBHEF=1" "PBHEF=2" "$COOL" "$COOL PBHEF=1" "$COOL PBHEF=2")

printf "%-14s %-5s %-7s %-9s %-5s %-12s %s\n" config rc steps finished nan injected/exp coupled-rate
for i in "${!names[@]}"; do
    name=${names[$i]}
    if [ -n "$WANTED" ]; then case " $WANTED " in *" $name "*) ;; *) continue ;; esac; fi
    run=$DIR/runs/$name
    rm -rf "$run"; mkdir -p "$run"

    { echo "$BASE"; for f in ${flags[$i]}; do echo "$f"; done; } > "$run/Config.sh"
    sed "s|^OutputDir .*|OutputDir  $run/output|" $DIR/isolated_halo.params > "$run/run.params"

    rm -f GIZMO_config.h compile_time_info.c
    if ! make CONFIG="$run/Config.sh" EXEC=$EXE -j8 > "$run/build.log" 2>&1; then
        printf "%-14s %s\n" "$name" "BUILD FAILED, see $run/build.log"; continue
    fi

    timeout 5400 mpirun -np "$NTASK" ./$EXE "$run/run.params" > "$run/run.log" 2>&1
    rc=$?
    rm -rf "$run/output/restartfiles"
    steps=$(grep -c '^Step' "$run/output/cpu.txt" 2>/dev/null || echo 0)
    fin=$(grep -ac 'Simulation ends' "$run/run.log")
    nan=$(grep -ac 'nan' "$run/run.log")
    log=$run/output/pbhef_energy.txt
    ratio=$(awk 'NR>1 && $1>0 && (NF<9 || $9>0) {r=$6/$7} END{if(r!="") printf "%.6f", r}' "$log" 2>/dev/null)
    frac=$(awk 'NR>1 && NF>8 && $10>0 {f=$9; n=1} END{if(n) printf "%.6f", f}' "$log" 2>/dev/null)
    printf "%-14s %-5s %-7s %-9s %-5s %-12s %s\n" "$name" "$rc" "$steps" "$fin" "$nan" "${ratio:--}" "${frac:--}"
done
rm -f $EXE GIZMO_config.h compile_time_info.c
