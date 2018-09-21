#//------------------------------------------------------------------------------
# Hercules-mpi: MPI implementation of Hercules transcriptomics analysis
# Dr. Jamie Alnasir
#
# Copyright (c) 2018, Dr. Jamie Alnasir, all rights reserved
#
# (CSSB) Center for Systems and Synthetic Biology
# Royal Holloway University of London
#
#//------------------------------------------------------------------------------
#// Configuration file
#//------------------------------------------------------------------------------


#//------------------------------------------------------------------------------
# GTF Annotation file path (on local, shared filesystem)

_REMOVE_INTERMEDIATES_ = True;

_DISK_IO_WAIT_ = 0.300; # 300 ms

# REGULAR JOB FILES

# GTF annotation
_LFS_GTF_FILE_FILTERED_  = "/home/user/Data/PhD/flybase2006-exons.gtf";

# SAM reads
_LFS_SAM_READS_ = "/home/user/Data/PhD/sam_reads.sam";

# Working directory, use trailing /
_LFS_WORKING_ = "/home/user/Data/PhD/_working/";


# Path to save report of ALL fourmer correlations
_LFS_REPORT_ALL_FOURMERS = _LFS_WORKING_ + "ALL-FOURMERS.txt";

# Path to save HTML report
_LFS_HTML_REPORT = "/home/user/mpi/Hercules-report.html";

#//------------------------------------------------------------------------------
