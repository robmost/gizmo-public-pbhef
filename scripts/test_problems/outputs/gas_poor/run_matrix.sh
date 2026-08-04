#!/bin/bash
# Each cavity radius under the receiver method and the donor method, with controls. Run from the
# repository root:
#
#     scripts/test_problems/outputs/gas_poor/run_matrix.sh [ntasks]
#
# The columns are the fraction of the analytic energy the gas received over the run, and for donor
# runs the lowest per-step fraction of the intended rate that found a gas neighbour to land on. The
# lowest, not the last: a starved donor heats the gas it can reach, a little of that gas expands into
# the cavity, and once anything at all is in reach the fraction reads one again, so the final value
# says nothing about what happened. The last two entries
# lower SofteningHalo to put the ceiling on the donor kernel below the cavity radius, which is the
# only way in a box this size to make a donor fail to reach gas at all.
set -u
NTASK=${1:-4}
DIR=scripts/test_problems/outputs/gas_poor
EXE=GIZMO_gas_poor

BASE="BOX_PERIODIC
BOX_SPATIAL_DIMENSION=2
SELFGRAVITY_OFF
EOS_GAMMA=(5.0/3.0)
HYDRO_MESHLESS_FINITE_MASS"

#      name                 radius  softening  extra config flags
runs=("r0_base                 0     0.078     "
      "r0_pbhef1               0     0.078     PBHEF=1"
      "r0_pbhef2               0     0.078     PBHEF=2"
      "r2_pbhef1               2     0.078     PBHEF=1"
      "r2_pbhef2               2     0.078     PBHEF=2"
      "r4_pbhef1               4     0.078     PBHEF=1"
      "r4_pbhef2               4     0.078     PBHEF=2"
      "r6_pbhef1               6     0.078     PBHEF=1"
      "r6_pbhef2               6     0.078     PBHEF=2"
      "r9_pbhef1               9     0.078     PBHEF=1"
      "r9_pbhef2               9     0.078     PBHEF=2"
      "r9_pbhef2_mass          9     0.078     PBHEF=2 PBHEF_MASS_WEIGHTS"
      "r9_pbhef1_shortreach    9     0.020     PBHEF=1"
      "r9_pbhef2_shortreach    9     0.020     PBHEF=2")

# lowest coupled fraction over the run, and the donors that were starved when it happened. A receiver
# build writes eight columns and has neither, so it reports nothing here rather than a misleading zero
worst_coupled() {
    awk 'NF>8 && $10>0 && ($9<lo || !n) {lo=$9; miss=$10; n=1}
         END {if(n) {printf "%.6f", lo} else {printf "-"}}' "$1"
}

printf "%-22s %-5s %-7s %-5s %-10s %s\n" config rc steps nan whole-run worst-coupled
for entry in "${runs[@]}"; do
    set -- $entry
    name=$1; radius=$2; soft=$3; shift 3
    run=$DIR/runs/$name
    rm -rf "$run"; mkdir -p "$run"

    { echo "$BASE"; for f in "$@"; do echo "$f"; done; } > "$run/Config.sh"
    sed -e "s|^OutputDir .*|OutputDir  $run/output|" \
        -e "s|^InitCondFile .*|InitCondFile  scripts/test_problems/ics/gas_poor_r${radius}_ics|" \
        -e "s|^SofteningHalo .*|SofteningHalo  $soft|" \
        -e "s|^SofteningHaloMaxPhys .*|SofteningHaloMaxPhys  $soft|" \
        $DIR/gas_poor.params > "$run/run.params"

    rm -f GIZMO_config.h compile_time_info.c
    if ! make CONFIG="$run/Config.sh" EXEC=$EXE -j8 > "$run/build.log" 2>&1; then
        printf "%-22s %s\n" "$name" "BUILD FAILED, see $run/build.log"; continue
    fi

    timeout 7200 mpirun -np "$NTASK" ./$EXE "$run/run.params" > "$run/run.log" 2>&1
    rc=$?
    rm -rf "$run/output/restartfiles"
    steps=$(grep -c '^Step' "$run/output/cpu.txt" 2>/dev/null || echo 0)
    nan=$(grep -ac 'nan' "$run/run.log")
    log=$run/output/pbhef_energy.txt
    ratio=$(awk 'NR>1 && $7>0 {r=$6/$7} END{if(r!="") printf "%.6f", r}' "$log" 2>/dev/null)
    if [ -f "$log" ]; then frac=$(worst_coupled "$log"); else frac="-"; fi
    printf "%-22s %-5s %-7s %-5s %-10s %s\n" "$name" "$rc" "$steps" "$nan" "${ratio:--}" "${frac:--}"
done
rm -f $EXE GIZMO_config.h compile_time_info.c
