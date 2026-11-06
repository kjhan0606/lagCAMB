#!/bin/bash
#SBATCH --job-name=bimond-history
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/bimond_scalar_history_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/bimond_scalar_history_%j.err

set -Eeuo pipefail
umask 077
: "${SLURM_JOB_ID:?Run on Grammar through Slurm}"
: "${BIMOND_EXPECTED_REVISION:?Specify the reviewed base revision}"
repo=/home/kjhan/BackUP/rar/code/lagCAMB
run=/gpfs/kjhan/RAR_QUMOND/runs/bimond_scalar_history_${SLURM_JOB_ID}
test "$(git -C "$repo" rev-parse HEAD)" = "$BIMOND_EXPECTED_REVISION"
test ! -e "$run"
mkdir -p "$run/source/camb" "$run/source/fortran/tests"
files=(
  camb/bimond_symbolic.py
  camb/bimond_scalar.py
  camb/bimond_history.py
  fortran/tests/verify_bimond_scalar_operator.py
  fortran/tests/verify_bimond_scalar_history.py
  fortran/tests/run_bimond_scalar_history.sbatch
)
for path in "${files[@]}"; do
  cp "$repo/$path" "$run/source/$path"
done
(cd "$run/source" && sha256sum "${files[@]}" > "$run/source.sha256")
git -C "$repo" rev-parse HEAD > "$run/base_revision.txt"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
# Fundamental-matrix diagnostic only; zero RAMSES snapshots or hydro outputs.
env -u LD_LIBRARY_PATH /home/kjhan/miniconda3/bin/python3 -u \
  "$run/source/fortran/tests/verify_bimond_scalar_history.py" \
  --output "$run/scalar_history.json" --dps 80
(cd "$run/source" && sha256sum -c "$run/source.sha256")
(cd "$run" && sha256sum scalar_history.json source.sha256 base_revision.txt > artifacts.sha256)
printf '%s\n' 'Finite-interval diagnostic complete; no cosmological transfer-function admission.' > "$run/COMPLETE"
