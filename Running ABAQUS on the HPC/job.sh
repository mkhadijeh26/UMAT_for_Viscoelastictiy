#!/bin/bash
#SBATCH --job-name="abaqus_visco"
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --cpus-per-task=1
#SBATCH --time=10:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --account=research-ceg-es
#SBATCH --output=abaqus_job.out

# Load required modules
module purge  # Clear all modules first
module load intel/oneapi-all
module load 2024r1
module load openmpi
module load opencoarrays
module load abaqus/2022


# First compile the UMAT code - modified for subroutine compilation
echo "Compiling UMAT code..."
gfortran -c -fcoarray=lib -fPIC -I"${ABAQUS_INC_DIR}" ViscoelasticityCode.f -o ViscoelasticityCode.o -L${OPENCOARRAYS_ROOT}/lib64 -lcaf_mpi

# Check if compilation was successful
if [ $? -eq 0 ]; then
    echo "Compilation successful, running Abaqus..."
    
    # Run ABAQUS with UMAT
    abq2022 cpus=`nproc` mp_mode=threads job=PlainStrainViscoElastic2 input=PlainStrainViscoElastic2.inp user=ViscoelasticityCode.f interactive

else
    echo "Compilation failed! Check the compiler output above."
    exit 1
fi