#!/bin/bash
# Timestep convergence for both coupling methods. Run from the repository root:
#
#     scripts/test_problems/outputs/timestep_convergence/run_matrix.sh [ntasks]
#
# The density-jump problem is run at a sequence of MaxSizeTimestep values, halving each time. That
# parameter sets the step of the dark matter, which has no acceleration in this problem, so halving it
# doubles how often the donors deposit while leaving the gas free to sub-cycle below them under the
# Courant condition. Two things should converge: the energy the gas receives against the analytic
# total, and the final state of the gas.
#
# Everything this generates lives under runs/<name>/. It reuses the density-jump IC and parameters,
# which have to exist first:
#
#     python3 scripts/test_problems/outputs/density_jump/make_ics.py
set -u
NTASK=${1:-4}
DIR=scripts/test_problems/outputs/timestep_convergence
SRC=scripts/test_problems/outputs/density_jump
EXE=GIZMO_timestep_convergence

BASE="BOX_PERIODIC
BOX_SPATIAL_DIMENSION=2
SELFGRAVITY_OFF
EOS_GAMMA=(5.0/3.0)
HYDRO_MESHLESS_FINITE_MASS"

steps=(0.02 0.01 0.005 0.0025 0.00125)
methods=(pbhef1 pbhef2)

printf "%-24s %-5s %-7s %-5s %s\n" config rc steps nan injected/expected
for method in "${methods[@]}"; do
    for step in "${steps[@]}"; do
        name="dt${step}_${method}"
        run=$DIR/runs/$name
        rm -rf "$run"; mkdir -p "$run"

        case $method in pbhef1) flag=PBHEF=1;; pbhef2) flag=PBHEF=2;; esac
        { echo "$BASE"; echo "$flag"; } > "$run/Config.sh"
        sed -e "s|^OutputDir .*|OutputDir  $run/output|" \
            -e "s|^MaxSizeTimestep .*|MaxSizeTimestep  $step|" \
            -e "s|^TimeBetSnapshot .*|TimeBetSnapshot  1|" \
            $SRC/density_jump.params > "$run/run.params"

        rm -f GIZMO_config.h compile_time_info.c
        if ! make CONFIG="$run/Config.sh" EXEC=$EXE -j8 > "$run/build.log" 2>&1; then
            printf "%-24s %s\n" "$name" "BUILD FAILED, see $run/build.log"; continue
        fi

        timeout 7200 mpirun -np "$NTASK" ./$EXE "$run/run.params" > "$run/run.log" 2>&1
        rc=$?
        rm -rf "$run/output/restartfiles"
        nsteps=$(grep -c '^Step' "$run/output/cpu.txt" 2>/dev/null || echo 0)
        nan=$(grep -ac 'nan' "$run/run.log")
        log=$run/output/pbhef_energy.txt
        ratio=$(awk 'NR>1 && $1>0 && (NF<9 || $9>0.5) {r=$6/$7} END{if(r!="") printf "%.9f", r}' "$log" 2>/dev/null)
        printf "%-24s %-5s %-7s %-5s %s\n" "$name" "$rc" "$nsteps" "$nan" "${ratio:--}"
    done
done
# the executable goes, but GIZMO_config.h stays: allvars.h includes it, so a language server cannot
# parse any file in the tree without one. Each build above deletes it first, which is what stops a
# stale configuration being reused; removing it here as well only breaks editor navigation.
rm -f $EXE
