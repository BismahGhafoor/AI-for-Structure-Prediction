#!/usr/bin/env bash

###############################################################################################################################

#
# Development AlphaFold2 slurm script
# PLEASE DO NOT CHANGE ANYTHING EXCEPT THE LAST LINE
# If you change anything except the last line, I WILL NOT HELP YOU!
#
# To use:
# 1. Create a new working directory, and copy this file there (rename, if you like)
# 2. Create a FASTA file with one or multiple sequences (ex: test.fasta)
# 3. Set the required parameters on the last line of this file
#    (input and output file names, last-template date, multimer?, etc,...)
# 4. Submit your job like:
#    sbatch submit_alphafold.sh
# 5. Check the run with `squeue`, or delete with `scancel ###`, where ### is your slurm job number
# 6. Monitor the output with `tail -f alphafold.err` (normal output is in the .err file, the .out file should be empty.)
#
###############################################################################################################################

### Inherit all current environment variables
#SBATCH --export=ALL

### Job name
#SBATCH --job-name=AP_predict

### Partition name
#SBATCH -p gpu

### Specify the number of tasks (mpi rank) and cpus-per-task (threads) for your job.
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --hint=nomultithread
#SBATCH -C broadwell
### The anticipated run-time for your job, where walltime=D-HH:MM:SS
#SBATCH --time=5-00:00:00

### Setup output and error files
#SBATCH --out=alphapulldown_predict_%A_%a.out
#SBATCH --error=alphapulldown_predict_%A_%a.err

### Switch to the working directory;
cd $SLURM_SUBMIT_DIR

### Setup the environment
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

export XLA_PYTHON_CLIENT_MEM_FRACTION='2.0'
export TF_FORCE_UNIFIED_MEMORY='1'

echo AlphaPulldown 0.30.1
echo $SLURM_JOB_ID
echo $SLURM_JOB_NODELIST

###############################################################################################################################

export RUN_DIR="/home/bg171/Project/p300"
export PROTEIN_LIST_PATHS="/home/bg171/Project/txtFile/custom_filtered.txt"

# Mode can be "pulldown", "all_vs_all", "homo-oligomer", "custom"
export MODE="custom"

# Defaults, can be changed if you wish
export CYCLES=3
export PREDICTIONS_PER_MODEL=1

# Provide a default value for SLURM_ARRAY_TASK_ID if not set
JOB_INDEX=${SLURM_ARRAY_TASK_ID:-0}

#### Do not change anything below this line ####
srun run_multimer_jobs.py \
    --protein_lists="${PROTEIN_LIST_PATHS}" \
    --monomer_objects_dir="/home/dp418/af_pulldown/p300/alignments, /home/bg171/Project/p300/Alignments" \
    --output_path="${RUN_DIR}/Predictions" \
    --mode="${MODE}" \
    --num_cycle="${CYCLES}" \
    --num_predictions_per_model="${PREDICTIONS_PER_MODEL}" \
    --data_dir="/net/common/alphafold_data/v2.3.1" \
    --job_index="${JOB_INDEX}"
