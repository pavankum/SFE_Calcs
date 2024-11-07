#! /bin/bash
#SBATCH --job-name=soleq2
#SBATCH --partition=free
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --account=DMOBLEY_LAB
#SBATCH --time=24:00:00
#SBATCH --mem=1gb
#SBATCH --array=0-19

#########################################################################
#########    Never change anything in this section              #########
#########    unless you know what you are doing of course !!!!  #########
#########################################################################

# SLURM environment variables:
# SLURM_JOB_USER             The user who started the job
# SLURM_ARRAY_TASK_ID        Job array ID (index) number.
# SLURM_ARRAY_JOB_ID         Job array’s master job ID number.
# SLURM_JOB_CPUS_PER_NODE    Count of processors available to the job on this node.
# SLURM_JOB_ID                The ID of the job allocation.
# SLURM_JOB_NAME             Name of the job.
# SLURM_JOB_NODELIST         List of nodes allocated to the job.
# SLURM_JOB_NUM_NODES        Total number of nodes in the job’s resource allocation.
# SLURM_JOB_PARTITION        Name of the partition in which the job is running.
# SLURM_NODEID               ID of the nodes allocated.
# SLURMD_NODENAME            Names of all the allocated nodes.
# SLURM_NTASKS               Same as -n, --ntasks
# SLURM_NTASKS_PER_CORE      Number of tasks requested per core. [If specified]
# SLURM_NTASKS_PER_NODE      Number of tasks requested per node. [If specified]
# SLURM_NTASKS_PER_SOCKET    Number of tasks requested per socket. [If specified]
# SLURM_PROCID               The MPI rank (or relative process ID) of the current process.
# SLURM_SUBMIT_DIR           The directory from which sbatch was invoked.
# SLURM_SUBMIT_HOST          The hostname of the computer from which sbatch was invoked.
# SLURM_TASKS_PER_NODE       Number of tasks to be initiated on each node. In the same order as SLURM_JOB_NODELIST.


# ---------------- ARRAY JOB SETTINGS ------------------ #
LAMBDA=$SLURM_ARRAY_TASK_ID
# -------------------------------------------------------#


#--------------
# Define some simpler local evironment variables:
nprocs=$SLURM_NTASKS
nnodes=$SLURM_JOB_NUM_NODES

#Function to call to run the actual code
slurm_startjob(){
#----------------- Actual calculation command goes here: ---------------------------

. ~/.bashrc
dir=${SLURM_SUBMIT_DIR}
cd $SLURM_SUBMIT_DIR
module purge
module load gromacs/2021.2/gcc.8.4.0

echo "Submitting ${ID}"
echo "Job directory: ${SLURM_SUBMIT_DIR}"


MDP=$dir/complete_files
mkdir lambda_$LAMBDA
cd lambda_$LAMBDA




##############################
# ENERGY MINIMIZATION STEP  #
##############################
echo "Starting minimization for lambda = $LAMBDA..."

mkdir EM
cd EM

gmx grompp -f $MDP/em.$LAMBDA.mdp -c $dir/0001.gro -p $dir/0001.top -o em.tpr -maxwarn 3 

gmx mdrun -deffnm em -s em.tpr -nt 1 

sleep 10

#####################
# NVT EQUILIBRATION #
#####################
echo "Starting constant volume equilibration..."

cd ../
mkdir NVT
cd NVT

gmx grompp -f $MDP/nvt.$LAMBDA.mdp -c ../EM/em.gro -p $dir/0001.top -o nvt.tpr -maxwarn 1

gmx mdrun -deffnm nvt -s nvt.tpr -nt 1

#####################
# NPT EQUILIBRATION #
#####################
echo "Starting constant pressure equilibration..."

cd ../
mkdir NPT
cd NPT

gmx grompp -f $MDP/npt.$LAMBDA.mdp -c ../NVT/nvt.gro -p $dir/0001.top -o npt.tpr -maxwarn 1

gmx mdrun -deffnm npt -s npt.tpr -nt 1


cd ../
mkdir Production_MD
cd Production_MD

gmx grompp -f $MDP/prod.$LAMBDA.mdp -c $dir/lambda_$LAMBDA/NPT/npt.gro -p $dir/0001.top -o prod.tpr -maxwarn 1



echo Job Done
}

# Function to echo informational output
slurm_info_out(){
# Informational output
echo "=================================== SLURM JOB ==================================="
date
echo
echo "The job will be started on the following node(s):"
echo $SLURM_JOB_NODELIST
echo
echo "Slurm User:         $SLURM_JOB_USER"
echo "Run Directory:      $(pwd)"
echo "Job ID:             ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "Job Name:           $SLURM_JOB_NAME"
echo "Partition:          $SLURM_JOB_PARTITION"
echo "Number of nodes:    $SLURM_JOB_NUM_NODES"
echo "Number of tasks:    $SLURM_NTASKS"
echo "Submitted From:     $SLURM_SUBMIT_HOST:$SLURM_SUBMIT_DIR"
echo "=================================== SLURM JOB ==================================="
echo
echo "--- SLURM job-script output ---"
}

slurm_info_out

slurm_startjob

