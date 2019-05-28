#!/bin/bash

CPUs=16
#BSUB -J "Hercules"
#BSUB -n 16
#BSUB -P infotech
#BSUB -W 10:00
#BSUB -o output.%J
#BSUB -e errors.%J

# Path to Hercules MPI4PY programs: P1(Phase 1) = Hercules-mpi.py, P2(Phase 2) = Hercules-mpi-analysis.py
_PATH_P1_="/home/infotech/jalnasir/mpi/Hercules-mpi.py"
_PATH_P2_="/home/infotech/jalnasir/mpi/Hercules-mpi-analysis.py"
# Path to GTF annotation and SAM reads files
_PATH_GTF_="/home/infotech/jalnasir/Data/Drosophila/flybase2006-exons.gtf"
_PATH_SAM_="/home/infotech/jalnasir/Data/Drosophila/sam_reads.sam"
# Working folder, where 4-mer data and intermediate results are stored (generates many files)
_PATH_WRK_="/home/infotech/jalnasir/Data/Drosophila/_working"
# Output folder - where Hercules-report.html, All-fourmers.txt is to be written
_PATH_OUT_="/home/infotech/jalnasir/mpi/"

echo "The following hosts have been allocated for this job:"
cat $LSB_DJOB_HOSTFILE

echo "Running Hercules-mpi"

# Ensure working folder exists
mkdir -p $_PATH_WRK_

# Perform PhaseI
mpirun --hostfile=$LSB_DJOB_HOSTFILE -n $CPUs python $_PATH_P1_ -g $_PATH_GTF_ -s $_PATH_SAM_ -w $_PATH_WRK_

#echo "Cleaning motif CSV files..."

# Ensure all _ are converted to , in csv output from PhaseI
find $_PATH_WRK_ -type f -name '*.csv' -exec sed -i 's/_/,/g' {} \;

# Perform PhaseII
mpirun --hostfile=$LSB_DJOB_HOSTFILE -n $CPUs python $_PATH_P2_ -g $_PATH_GTF_ -s $_PATH_SAM_ -w $_PATH_WRK_ -o $_PATH_OUT_

# Produce GC box plot of Correlation vs GC content
python ./Hercules-analysis-plot.py
                                                                                                                             
