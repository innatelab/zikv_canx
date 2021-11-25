#! /bin/bash
#SBATCH -o /fs/pool/pool-innate-analysis/scratch/stukalov/logs/sdenolly_canxapms/%x_%j.log
#SBATCH -J sdenolly_canxapms_msglm
#SBATCH --get-user-env

#SBATCH --nodes=5
#SBATCH --ntasks-per-node=5
#SBATCH --mincpus=8
#SBATCH --cpus-per-task=8
#SBATCH --mem=24GB
#SBATCH --exclusive=user
##SBATCH --ntasks-per-socket=1
#SBATCH --mail-type=end
#SBATCH --mail-user=alexey.stukalov@tum.de
#SBATCH --export=NONE
#SBATCH --time=72:00:00
#SBATCH --wait-all-nodes=1

module load charliecloud

BIOINFO_ROOT=/fs/pool/pool-innate-bioinfo2
IMAGES_PATH=${BIOINFO_ROOT}/ccimages
CHUDIS_PATH=$(realpath $HOME/projects/chunk_dispatcher)
ANALYSIS_ROOT=$(realpath /fs/pool/pool-innate-analysis/$USER/analysis)
SCRATCH_ROOT=$(realpath /fs/pool/pool-innate-analysis/$USER/scratch)
CHARLIECLOUD_PATH=$(dirname $(which ch-run))

PROJECT_ID=sdenolly_canxapms_msglm
FIT_VERSION=20211125
CHUDIS_JOBID=${PROJECT_ID}_${FIT_VERSION}
CHUDIS_USER=ge68wan2

srun -c$SLURM_CPUS_PER_TASK --ntasks-per-node=$SLURM_NTASKS_PER_NODE --mem=$SLURM_MEM_PER_NODE --wait=0 \
     --no-kill --distribution=block --exclusive=user \
     -o "$SCRATCH_ROOT/logs/$PROJECT_ID/%x_%j_%t.log" \
${CHUDIS_PATH}/slurmstep_process_chunks.sh $CHUDIS_JOBID $CHUDIS_USER \
"${CHARLIECLOUD_PATH}/ch-run $IMAGES_PATH/archpc.msglm_v050 \
  -t --no-home --unset-env='*PATH' \
  --set-env=$HOME/projects/adhoc/$PROJECT_ID/archpc.env \
  -b $HOME/projects/adhoc:/projects/adhoc \
  -b $ANALYSIS_ROOT:/analysis \
  -b $SCRATCH_ROOT:/scratch \
  -- Rscript /projects/adhoc/$PROJECT_ID/msglm_fit_chunk.R \
  $PROJECT_ID $SLURM_JOB_NAME $FIT_VERSION"

echo "All queues job #$SLURM_JOB_ID finished"

