#!/bin/bash
#SBATCH --job-name=f1          # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --time=168:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-user=filip.sjovall@solid.lth.se
#SBATCH --mail-type=END
#SBATCH -A lu2024-2-43
# For information on how the job went error etc, include %j to get uniqe identifier
#SBATCH -o process_%j.out
#SBATCH -e process_%j.err


# write this script to stdout-file - useful for scripting errors
cat $0
# change to the execution directory
cd $SNIC_TMP

#cp -p FictShapeOpt/ $SNIC_TMP

module purge
module load Julia

ulimit -s unlimited

cd $SNIC_TMP

date


# Run your Julia program
julia -O3 "$SLURM_SUBMIT_DIR/src/Linear/OptimizationLabFine.jl"

# Copy the output files from the compute node to the login node
cp -r "$SCRIPT_DIR/results" "$SLURM_SUBMIT_DIR"

date
