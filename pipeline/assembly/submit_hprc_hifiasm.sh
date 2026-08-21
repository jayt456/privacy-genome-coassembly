#!/usr/bin/env bash
# Submit HPRC 12x subsample + hifiasm array job.
# --array=0-46%5 runs at most 5 samples concurrently.
# Outputs go to /nfs/roberts/scratch/pi_hc878/jst72/hprc_results/hifiasm_12x/
#
# Usage:
#   ./submit_hprc_hifiasm.sh             # all 47 samples
#   ./submit_hprc_hifiasm.sh HG002       # single sample by name
set -euo pipefail

BASE="${COASM_BASE:-/nfs/roberts/project/pi_hc878/jst72}"
SCRATCH_BASE="${COASM_SCRATCH:-/nfs/roberts/scratch/pi_hc878/jst72}"
WORKER="$BASE/scripts/hprc_subsample_hifiasm.sh"

samples=(
    HG002    HG005    HG00438  HG00621  HG00673  HG00733  HG00735  HG00741
    HG01071  HG01106  HG01109  HG01123  HG01175  HG01243  HG01258  HG01358
    HG01361  HG01891  HG01928  HG01952  HG01978  HG02055  HG02080  HG02109
    HG02145  HG02148  HG02257  HG02486  HG02559  HG02572  HG02622  HG02630
    HG02717  HG02723  HG02818  HG02886  HG03098  HG03453  HG03486  HG03492
    HG03516  HG03540  HG03579  NA18906  NA19240  NA20129  NA21309
)

if [[ $# -eq 1 ]]; then
    # Single-sample mode: find index and submit one task
    target="$1"
    idx=-1
    for i in "${!samples[@]}"; do
        [[ "${samples[$i]}" == "$target" ]] && { idx=$i; break; }
    done
    if [[ $idx -lt 0 ]]; then
        echo "ERROR: sample '$target' not in list"
        exit 1
    fi
    # Check if already done
    gfa="$SCRATCH_BASE/hprc_results/hifiasm_12x/$target/${target}_12x.bp.hap1.p_ctg.gfa"
    if [[ -f "$gfa" ]]; then
        echo "$target: GFA exists — skipping"
        exit 0
    fi
    echo "Submitting single task: $target (array index $idx)"
    sbatch --array="$idx" "$WORKER"
else
    # Full array: skip already-done samples
    todo_indices=()
    skipped=()
    for i in "${!samples[@]}"; do
        s="${samples[$i]}"
        gfa="$SCRATCH_BASE/hprc_results/hifiasm_12x/$s/${s}_12x.bp.hap1.p_ctg.gfa"
        if [[ -f "$gfa" ]]; then
            skipped+=("$s")
        else
            todo_indices+=("$i")
        fi
    done

    if [[ ${#skipped[@]} -gt 0 ]]; then
        echo "Skipping ${#skipped[@]} already-assembled sample(s): ${skipped[*]}"
    fi

    if [[ ${#todo_indices[@]} -eq 0 ]]; then
        echo "All 47 samples already assembled."
        exit 0
    fi

    # Build comma-separated index list for --array
    array_spec=$(IFS=,; echo "${todo_indices[*]}")
    echo "Submitting ${#todo_indices[@]} sample(s) as array tasks: $array_spec"
    echo "(max 5 concurrent to manage disk usage)"
    sbatch --array="${array_spec}%5" "$WORKER"
fi
