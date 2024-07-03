#!/usr/bin/env bash

### Inherit all current environment variables
#SBATCH --export=ALL

### Job name
#SBATCH --job-name=AP_ADalign

### Partition name
#SBATCH -p gpu

### Specify the number of tasks (mpi rank) and cpus-per-task (threads) for your job.
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --gres=gpu:0
#SBATCH --hint=multithread
#SBATCH -C broadwell

### The anticipated run-time for your job, where walltime=D-HH:MM:SS
#SBATCH --time=10-00:00:00
#SBATCH --qos=low

### Setup output and error files
#SBATCH --out=alphapulldown_features_%A_%a.out
#SBATCH --error=alphapulldown_features_%A_%a.err

### Switch to the working directory;
cd $SLURM_SUBMIT_DIR

### Setup the environmnet
__conda_setup="$('/net/prog/anaconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/net/prog/anaconda3/etc/profile.d/conda.sh" ]; then
        . "/net/prog/anaconda3/etc/profile.d/conda.sh"
    else
        export PATH="/net/prog/anaconda3/bin:$PATH"
    fi
fi

unset __conda_setup
conda activate alphafold_v2.3.1

echo AlphaPulldown 0.30.1
echo $SLURM_JOB_ID
echo $SLURM_JOB_NODELIST

## Make sure to use full paths here

export FASTA_DIR="/home/bg171/Project/ligand_fastafiles"
export RUN_DIR="/home/bg171/Project/p300/"

FASTA_FILE_PATHS=$(find $FASTA_DIR -name "*.fa" | paste -sd "," -)

#### Do not change below here ####

srun create_individual_features.py \
--fasta_paths=${FASTA_FILE_PATHS} \
--output_dir=${RUN_DIR}/ligand_Alignments \
--skip_existing=True \
--save_msa_files=True \
--max_template_date=2023-02-13 \
--data_dir=/net/common/alphafold_data/v2.3.1/ \
--seq_index=$SLURM_ARRAY_TASK_ID
