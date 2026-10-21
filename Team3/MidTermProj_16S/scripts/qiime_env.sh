# Shared helper: run QIIME2 2024.10 (ARC Apptainer image, same version as the paper)
module load apptainer/1.4.0
P=$HOME/fireblight_canker_reanalysis
export TMPDIR=/scratch/$USER/q2tmp; mkdir -p $TMPDIR $P/qiime/.q2home
SIF=/apps/common/software/QIIME2/2024.10/qiime2-amplicon-2024.10.sif
q2() { apptainer exec -B $P/qiime/.q2home:/home/qiime2 -B /scratch -B $HOME --env TMPDIR=$TMPDIR $SIF "$@"; }
qiime() { q2 qiime "$@"; }
