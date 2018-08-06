#//==============================================================================
#//    _    _                     _           
#//   | |  | |                   | |          
#//   | |__| | ___ _ __ ___ _   _| | ___  ___ 
#//   |  __  |/ _ \ '__/ __| | | | |/ _ \/ __|
#//   | |  | |  __/ | | (__| |_| | |  __/\__ \
#//   |_|  |_|\___|_|  \___|\__,_|_|\___||___/ MPI
#//
#//==============================================================================

# Hercules-mpi: MPI implementation of Hercules transcriptomics analysis
# Jamie Alnasir
#
# Copyright (c) 2018, Dr. Jamie Alnasir, all rights reserved
#
# (CSSB) Center for Systems and Synthetic Biology
# Royal Holloway University of London
#
# Based on my PhD under the supervision of Dr. Hugh Shanahan
#
# Alnasir, J. & Shanahan, H. P. (2017). Transcriptomics: Quantifying non-uniform read distribution using MapReduce.. International Journal of Foundations of Computer Science (Forthcoming). [ Read (preprint) ]
#
# Alnasir, J., & Shanahan, H. (2017). A novel method to detect bias in Short Read NGS RNA-seq data. Journal of Integrative Bioinformatics, 14(3). [ Read ]
#
# Alnasir, J., & Shanahan, H. P. (2015). Transcriptomics on Spark Workshop - Introducing Hercules - an Apache Spark MapReduce algorithm for quantifying non-uniform gene expression. CloudTech'16, Marrakech, Morocco. [ Read]



#//------------------------------------------------------------------------------
#// MPI program (requires mpi4py and a configured associated mpi distribution)
#//------------------------------------------------------------------------------
from mpi4py import MPI;
from optparse import OptionParser
from time import sleep;
import sys;
import os;
import os.path;
import re;
from itertools import groupby;


# MPI initialization
comm = MPI.COMM_WORLD;			# get MPI communicator object
size = comm.Get_size();			# total number of processes
rank = comm.Get_rank();			# rank of this process
name = MPI.Get_processor_name();	# usually returns hostname
status = MPI.Status();			# get MPI status object


# Fixed configuration parameters (not specified using Option Parser)
CONF_DISK_IO_WAIT_ = 0.300; # 300 ms

# Configuration parameters (now set using Option Parser)
CONF_REMOVE_INTERMEDIATES_  = False;
CONF_LFS_WORKING_           = "";
CONF_LFS_GTF_FILE_FILTERED_ = "";
CONF_LFS_SAM_READS_         = "";



# Nucleotide bases
bases = "ATGC";


def enum(*sequential, **named):
	"""Handy way to fake an enumerated type in Python
	http://stackoverflow.com/questions/36932/how-can-i-represent-an-enum-in-python
	"""
	enums = dict(zip(sequential, range(len(sequential))), **named)
	return type('Enum', (), enums)

# Define MPI message tags
tags = enum('PING', 'WORK', 'INITIALISE', 'UNLOADGTF', 'RETIRE', 'SYNC', 'SERIALISED');

# Hercules work tags
tasks = enum('GTFSAM', 'GTFSAMCOMBINE', 'MOTIFMAP', 'MOTIFMAPCOMBINE', 'MOTIFREDUCE', 'MOTIFREDUCECOMBINE', 'GCMAP', 'GCMAPCOMBINE');


lstReceived = [];


#//------------------------------------------------------------------------------
# MPI implementation functions
#//------------------------------------------------------------------------------

def pprint(s):
# Master/Worker process print
	if (rank == 0):
		p = "Master on {}".format(name);
	else:
		p = "Process {} on {}".format(rank, name);

	print "{}: {}".format(p, s);


def ping(target):
# Ping a worker process

	pprint("pinging process {}".format(target));
	comm.send(0, dest=target, tag=tags.PING);

	# recieve the PING
	result = comm.recv(source=MPI.ANY_SOURCE, tag=tags.PING, status=status);
	source = status.Get_source();
	tag = status.Get_tag();
	return (tag == tags.PING); # assert return ping


def syncMaster():
# Send message from Worker to Master to syncronise.
# Master waits for ALL Workers messages before distributing the
# next load of work

	pprint("declaring readiness [sync] to Master");
	comm.send(0, dest=0, tag=tags.SYNC);


def syncWait():
# Wait for the syncronise (SYNC) signal to be sent from ALL worker nodes

	for i in range(1, size):
		result = comm.recv(source=MPI.ANY_SOURCE, tag=tags.SYNC, status=status);
		source = status.Get_source();
		tag = status.Get_tag();		
		pprint("received [sync] from {}".format(source));


def sendSerialised():
# Send message from Worker to Master to confirm serialised task has completed
# Master waits for EACH serialised task on Workers to complete before sending
# the next serialised task.

	pprint("succefully serialised task, sending [serialised] to Master");
	comm.send(0, dest=0, tag=tags.SERIALISED);


def serialiseWork(task, params):
# Serialise work across ALL worker nodes, waiting for each task to finish
# before allocating the next task. Necessary to prevent nodes performing I/O
# at the same time.
# Worker node must send SERIALISED when task is completed.

	for i in range(1, size):
		# Send the task
		pprint("serialising work(task={}, params={}) to worker process {}:".format(task, params, i));
		comm.send([task, params], dest=i, tag=tags.WORK);

		# MUST wait to receive SERIALISEDONE -- i.e. individual worker node to to finish it's task,
		result = comm.recv(source=MPI.ANY_SOURCE, tag=tags.SERIALISED, status=status);
		source = status.Get_source();
		tag = status.Get_tag();		
		pprint("received [serialise] from {}".format(source));
		
		# wait a little on IO before sending the next serialised task
		sleep(CONF_DISK_IO_WAIT_);


def initialise(target):
# Initialise a worker process

	print("requesting process {} to initialise".format(target));
	comm.send(0, dest=target, tag=tags.INITIALISE);

	# recieve the PING
	#result = comm.recv(source=MPI.ANY_SOURCE, tag=tags.PING, status=status);
	#source = status.Get_source();
	#tag = status.Get_tag();
	#return (tag == tags.PING); # assert return ping


def unloadGTFcache(target):
# Unload a worker process's GTF cache

	print("requesting process {} to dispose of it's GTF cache".format(target));
	comm.send(0, dest=target, tag=tags.UNLOADGTF);

def retire(target):
# Retire a worker process

	pprint("requesting process {} to retire".format(target));
	comm.send(0, dest=target, tag=tags.RETIRE);

	# recieve the PING
	#result = comm.recv(source=MPI.ANY_SOURCE, tag=tags.PING, status=status);
	#source = status.Get_source();
	#tag = status.Get_tag();
	#return (tag == tags.PING); # assert return ping


def scatter(data):
	result = comm.scatter(data, root=0);
	print 'rank',rank,'has data:',result;
	return result;

def doScatterWork(data):
	result = data + 1;
	return result;

def chunkRanges(start, stop, step):
	return [(n, min(n+step, stop)) for n in xrange(start, stop, step)];


def sendWork(target, task, params):
# send work to a target worker process

	pprint("sending work(task={}, params={}) to worker process {}:".format(task, params, target));
	comm.send([task, params], dest=target, tag=tags.WORK);


def allocWorkGTFSAM():
# send out GTFSAM work to nodes
	
	pprint("analysing SAM reads file (this may take some time): ".format(CONF_LFS_SAM_READS_));
	gtfLen = fileLineLen(CONF_LFS_SAM_READS_);
	chunkLines = gtfLen / (size - 1);
	pprint("SAM analysis done.".format(gtfLen));
	pprint("Number of reads =".format(gtfLen));
	pprint("allocating {} reads per worker process".format(chunkLines));
	params = chunkRanges(0, gtfLen, chunkLines);

	for i in range(1, size):
		sendWork(i, tasks.GTFSAM, params[i - 1]);		


def combineWorkGTFSAM():
# aggregate work from GTFSAM step

	# remove any previous result	
	pprint("attempting to remove any previous result: GTF-SAM.final.dat");
	gtfsamFile = CONF_LFS_WORKING_ + "GTF-SAM.final.dat";
	os.popen("rm -f " + gtfsamFile);

	# old method, non serialised
	#for i in range(1, size):
	#	sendWork(i, tasks.GTFSAMCOMBINE, []);

	serialiseWork(tasks.GTFSAMCOMBINE, []);
	


def allocWorkMotifMap(aMotif):
# send out MOTIF map work to nodes

	gtfsamFile      = CONF_LFS_WORKING_ + "GTF-SAM.final.dat";
	gtfsamFileMotif = CONF_LFS_WORKING_ + "GTF-SAM." + aMotif + ".dat";
	os.popen("grep " + aMotif + " " + gtfsamFile + " >" + gtfsamFileMotif);

	pprint("analysing reads with MOTIF={} (this may take some time): ".format(aMotif, gtfsamFileMotif));
	motifLen = fileLineLen(gtfsamFileMotif);
	chunkLines = motifLen / (size - 1);
	pprint("MOTIF reads analysis done.".format(motifLen));
	pprint("Number of reads =".format(motifLen));
	pprint("allocating {} reads per worker process".format(chunkLines));
	params = chunkRanges(0, motifLen, chunkLines);

	for i in range(1, size):
		sendWork(i, tasks.MOTIFMAP, [aMotif, params[i - 1]]);


def combineMotifMap(aMotif):
# aggregate work from MOTIF map step

	# remove any previous result	
	motifmapFile = CONF_LFS_WORKING_ + "MOTIF." + aMotif + ".final.dat";
	os.popen("rm -f " + motifmapFile);

	# old method, non serialised
	#for i in range(1, size):
	#	sendWork(i, tasks.MOTIFMAPCOMBINE, aMotif);

	serialiseWork(tasks.MOTIFMAPCOMBINE, aMotif);
	
					
	# remove GTF-SAM motif file from LFS - but only after serialised work completed
	gtfsamFileMotif = CONF_LFS_WORKING_ + "GTF-SAM." + aMotif + ".dat";
	if CONF_REMOVE_INTERMEDIATES_:
		os.popen("rm -f " + gtfsamFileMotif);




	
def allocWorkMotifReduce(aMotif):
# send out MOTIF reduce work to nodes

	motifmapFile = CONF_LFS_WORKING_ + "MOTIF." + aMotif + ".final.dat";
	#os.popen("rm -f " + motifmapFile);

	pprint("analysing motif-mapped reads with MOTIF={} (this may take some time): ".format(aMotif, motifmapFile));
	motifLen = fileLineLen(motifmapFile);
	chunkLines = motifLen / (size - 1);
	pprint("MOTIF reads analysis done.".format(motifLen));
	pprint("Number of reads =".format(motifLen));
	pprint("allocating {} reads per worker process".format(chunkLines));
	params = chunkRanges(0, motifLen, chunkLines);

	for i in range(1, size):
		sendWork(i, tasks.MOTIFREDUCE, [aMotif, params[i - 1]]);

		
def combineMotifReduce(aMotif):
# aggregate work from MOTIF reduce step

	# remove any previous result	
	motifreduceFile = CONF_LFS_WORKING_ + "MOTIF.REDUCE." + aMotif + ".final.dat";
	motifreduceCSVFile = CONF_LFS_WORKING_ + "herc-final-" + aMotif + ".csv";
 	
	os.popen("rm -f " + motifreduceFile);

	serialiseWork(tasks.MOTIFREDUCECOMBINE, aMotif);
	
	# Perform a final reduce on the final file because the splits may have been
	# split accross a key amd may not have been completed in order on the workers.
	motifreduceFileSorted = motifreduceFile + ".sorted";
	# Sort the entire intermediate result
	os.popen("sort " + motifreduceFile + " >" + motifreduceFileSorted);
	
	lstMR = loadText(motifreduceFileSorted);
	lstMRResult = reduceTupleList(lstMR, False); # We've already pre-sorted
	saveText(motifreduceCSVFile, lstMRResult);

	# Replace _ with , for CSV format
        os.popen("sed -i 's/_/,/g' " + motifreduceCSVFile);

	
	if CONF_REMOVE_INTERMEDIATES_:
		motifmapFile = CONF_LFS_WORKING_ + "MOTIF." + aMotif + ".final.dat";
		motifreduceFile = CONF_LFS_WORKING_ + "MOTIF.REDUCE." + aMotif + ".final.dat*"; # and .sorted

		pprint("removing intermediate data from MOTIF map and reduce steps");
		os.popen("rm -f " + motifmapFile);
		os.popen("rm -f " + motifreduceFile);



def allocWorkSAMGC(rankChunkIndex):
# send out GC computation work to nodes
# rankChunkIndex is the GTF-SAM.?.dat file where the index corresponds to the rank that generated the file
	
	wrkFileSAM = CONF_LFS_WORKING_ + "GTF-SAM." + str(rankChunkIndex) + ".dat";
	pprint("Loading partitioned SAM reads for GC Computation (this may take some time): ".format(wrkFileSAM));

	for i in range(1, size):
		sendWork(i, tasks.GCMAP, rankChunkIndex);

def combineGCMap():
# aggregate work from GTFSAM step

	# remove any previous result	
	pprint("attempting to remove any previous result: GC.final.dat");
	GCDataFile = CONF_LFS_WORKING_ + "GC.final.dat";
	os.popen("rm -f " + GCDataFile);

	# old method, non serialised
	#for i in range(1, size):
	#	sendWork(i, tasks.GTFSAMCOMBINE, []);

	serialiseWork(tasks.GCMAPCOMBINE, []);

#//------------------------------------------------------------------------------
# Hercules RNA-Seq analysis functions
#//------------------------------------------------------------------------------

def getCol(lst, col):
	return [row[col] for row in lst];
	
def getChrList(chr):
	for chromes in lstGTF_chr:
		if (chromes[0] == chr):
			return lstGTF_feat[chromes[1]:chromes[2]];
			
def getMiddle(lst):
	if lst is None:
		return -1;
	return len(lst) // 2;
	

def fileLineLen(aFile):
# Quickly, cheaply count the lines in aFile
	i=0;
	with open(aFile) as f:
		for i, l in enumerate(f):
			pass
	return i + 1

def readlinesFileSection(aFile, start, end):
# Read only lines of aFile specified by 0-based start and end (inclusive)
	lstR = [];
	fp = open(aFile);
	for i, line in enumerate(fp):
		if (i < start):
			continue;
		elif (i >= start) and (i <= end):
			lstR.append(line);
		elif (i > end):
			break;
	fp.close()
	return lstR;

def loadText(aInFile):
# Load array of strings from text file
	with open(aInFile) as f:
		data = f.readlines();
		return data;

def saveText(aOutFile, data):
# Load array of strings from text file
	with open(aOutFile, "w") as f:
		for i in data:
			f.write(i + "\n");

def CigarToTupleList(aCigarStr):
# Parse a CIGAR string into a list of tuples [(length,operation),(length,operation),(length,operation)]
	return re.findall(r'(\d+)(\w)', aCigarStr);

def getFirstCigarAlignedMatchLen(aCigarStr):
# Return the Length of the first Alignment match 
# Kludge -- there can be multiple alignment matches separated by insertions, deletions or skips. We only use
# the first portion of the read that is an alignment match, regardless of whether it is a sequence match or not.
	lstCIG = CigarToTupleList(aCigarStr);
	if (lstCIG is None): return -1;
	for item in lstCIG:
		if (item[1] == 'M'):
			return long(item[0]);
	return -1;

def isCigarSingleAlignMatch(aCigarStr):
# Returns True if the CIGAR string is a single, contiguous aligned match
	lstCIG = CigarToTupleList(aCigarStr);
	if (lstCIG is None): return False;
	if (len(lstCIG) <> 1): return False;
	return (lstCIG[0][1] == 'M');

def _isPosWithinRange(feat_start, feat_end, pos):
	return (long(pos) >= long(feat_start)) and (long(pos) <= long(feat_end));

def _isWithinRange(feat_start, feat_end, ReadPos_start, ReadPos_end):
	# Old functionality
	# Return feature key if read start and end fall totally within the feature range
	#return long(ReadPos_start) >= long(feat_start) and long(ReadPos_end) <= long(feat_end);
	# New functionality (03/05/2016)
	# Return feature key if read falls within the feature range OR straddles feature range
	bLeftStraddle  = ( long(ReadPos_start) < long(feat_start) and _isPosWithinRange(feat_start, feat_end, ReadPos_end) );
	bRightStraddle = ( _isPosWithinRange(feat_start, feat_end, ReadPos_start) and long(ReadPos_end) > long(feat_end) );
	bStraddling	= (bLeftStraddle or bRightStraddle);
	#if bStraddling:
	#	print "straddling: ", ReadPos_start, ReadPos_end;
	bWithin		= long(ReadPos_start) >= long(feat_start) and long(ReadPos_end) <= long(feat_end);
	return bStraddling or bWithin;
	

def _gtf_lookup(lst, ReadPos_start, ReadPos_end):
# Recursive GTF lookup using binary search (sequential halving of the search list)
	m = getMiddle(lst);
	if (m == -1):
		return -1;
	if _isWithinRange(lst[m][1], lst[m][2], ReadPos_start, ReadPos_end):
		return lst[m][3];	 
	if (m == 0):
		return -1;
	if (int(ReadPos_start) < lst[m][1]):
		return _gtf_lookup(lst[:m], ReadPos_start, ReadPos_end);
	if (int(ReadPos_start) > lst[m][1]):
		return _gtf_lookup(lst[m:], ReadPos_start, ReadPos_end);		
	

def gtf_lookup(chr, ReadPos_start, ReadPos_end):
# Function to look up the feature range within which a read lies using the input
# read and the GTF annotation using a binary search technique

	key = _gtf_lookup(getChrList(chr), ReadPos_start, ReadPos_end);
	if (key is None):
		return -1;
	return key;

def BaseCounts(aSeq):
# Count frequencies of bases in given sequence string
	lstBaseCounts = [];
	for base in bases:
		lstBaseCounts.append([base, len(filter(lambda x: (x == base), aSeq))]);
	return lstBaseCounts;


def GC_Content(aSeq):
# Return GC content of given sequence string
	lstB = BaseCounts(aSeq);    
	dictB = dict(lstB);
	bases_all = int( dictB['A'] + dictB['T'] + dictB['G'] + dictB['C'] );
	bases_GC = int( dictB['G'] + dictB['C'] );
	gc = float(bases_GC) / float(bases_all) * 100;
	return gc;


def motifs(aSeqStr, aMotif):
# Return an array of motif string 0-based index positions within aSeqStr
	x = 0;
	lstResult = [];
	while (x <> -1):	
		if (x == 0):
			x = aSeqStr.find(aMotif, x);
			if (x <> -1):						
				lstResult.append(x);
				x = x + 1;							
				continue;
			else:
				return None;				
		else:	
			x = aSeqStr.find(aMotif, x + len(aMotif));						
		if (x <> -1):
			lstResult.append(x);			
	return lstResult;


def gtf_sam_mapper(sam_read_line):
	# GTF/SAM READ LOOKUP Map step
	# Sorts sequence alignment (sam read fragments) (v)
	# into GTF gene/coding regions (k)
	sam_data = sam_read_line.split('\t');
	chr  = sam_data[2];
	rpos = sam_data[3];
	rlen = getFirstCigarAlignedMatchLen(sam_data[5]);
	endpos = long(rpos) + long(rlen);
	#print 'looking up read @ pos ' + rpos + '-' + endpos;
	mapkey = gtf_lookup(chr, rpos, endpos);
#	if (mapkey is None):
#		_debug('None,' + sam_read_line);

	#print 'key-gen: ' + mapkey;
	return str(str(mapkey) + ", " + sam_read_line).strip();


def motif_mapper(aReadLine, _SEARCH_MOTIF_, bReadRelative=False):
	# MOTIF Map step
	# Itemise Motif occurrence positions (offset from
	# start of read-range)
	# Set bReadRelative true for motif positions relative to read-start
	# otherwise positions are relative to exon-start
	lstMotifCorrect = None;
	#pprint('READLINE=' + aReadLine);
	key, value = aReadLine.strip('\'').split(',');
	key = key.strip('(');
	value = value.rstrip(')\n');
	if (long(key) <> -1):
		key_start = long(key[:10]);
		key_end = long(key[10:]);
	else:
		key_start = 0;
		key_end = 0;
	sam_data = value.split('\t');
	rpos = long(sam_data[3]);
	rlen = getFirstCigarAlignedMatchLen(sam_data[5]);
	endpos = rpos + rlen;		
	offset = rpos - key_start;
	keylen =  key_end - key_start;
	#print "SEQ=" + sam_data[9];
	lstMotif = motifs(sam_data[9], _SEARCH_MOTIF_);

	if (bReadRelative):
	# No offset correction as position is relative to read-start
		lstMotifCorrect = lstMotif;
	else:
		# Correct for read-offset so position is relative to exon-start
		if lstMotif is not None:					
			lstMotifOffset = map(lambda x: x+offset, lstMotif);
			lstMotifCorrect = filter(lambda x: (x > 0 and x <= keylen), lstMotifOffset);

	if lstMotifCorrect is None:
		lstMotifCorrect = [];

	#print "motifs=" + str(lstMotifCorrect);
	return str(key) + ", " + str(lstMotifCorrect);


def vector_mapper(aVectorReadLine):
	# Motif Position map step
	# Input 1 value, return 1 or more value, must be called via flatMap
	#print "LINE: " + aVectorReadLine;
	key, value = aVectorReadLine.split(',', 1);
	key = key.strip();
	value = value.strip().lstrip('[').rstrip(']');
	lstVector = list(value.split(','));
	lstVector = [v.strip().replace('L','') for v in lstVector];
	lstVector = filter(lambda x: x.isdigit(), lstVector);

	lstResult = [];
	for motif_pos in lstVector:
		lstResult.append(key + "_" + motif_pos.strip().zfill(6).replace('L', '') + ",1");

	return lstResult;
	
	
def reduceTupleList(lstT, bPreSort):
	lstR = [];
	lstTmap = map(lambda x: [x.split(",")[0], int(x.split(",")[1])], lstT);
	if bPreSort:
		lstTmapSorted = sorted(lstTmap);
	else:
		lstTmapSorted = lstTmap;
	for key, group in groupby(lstTmapSorted, lambda x: x[0]):
		lstR.append(key + "," + str(sum([x[1] for x in group])));
	return lstR; #[i.replace("_",",") for i in lstR];

def gc_mapper(aReadLine):
	# GC_MAP Map step
	# Return the GC content of the given read
	key, value = aReadLine.strip('\'').split(',');
	key = key.strip('(');
	value = value.rstrip(')\n');
	sam_data = value.split('\t');
	gc = GC_Content(sam_data[9]);
	return ",".join([str(key), str(gc)]);


def computeGCsums(lstGCchunk):
# Return a list of sums and counts for aggregation and averaging later
    lstResult = [];
    lstGC = [];
    for i in lstGCchunk:
        exon, GC = i.strip().split(",");
        GC = float(GC);
        lstGC.append([exon, GC]);
    lstGC = sorted(lstGC, key=lambda x: x[0], reverse=False);
    for key, group in groupby(lstGC, lambda x: x[0]):
        sum = 0.0;
        count = 0;
        for i in group:
            #print i[1];
            sum = sum + float(i[1]);
            count += 1;
        lstResult.append(",".join([str(key), str(count), str(sum)]));
    return lstResult;


def computeFinalGCcontent(strFinalDatFile):
# Aggregate all GC sums and counts into final mean GC for each exon
	lstResult = [];
	lstTmp    = loadText(strFinalDatFile);
	lstGCsums = [];
	for i in lstTmp:
		Exon, Count, Sum = i.strip().split(",");
		lstGCsums.append([long(Exon), long(Count), float(Sum)]);
		
	lstGCsorted = sorted(lstGCsums, key=lambda x: x[0], reverse=False);
	for key, group in groupby(lstGCsorted, lambda x: x[0]):
		GCsum = 0.0;
		CountSum = 0;
		for i in group:
			CountSum = CountSum + i[1];
			GCsum = GCsum + float(i[2]);
		lstResult.append(",".join( [str(key).rjust(20, '0'), str( float( GCsum / CountSum ) ) ] ));
	lstResult[0] = lstResult[0].replace("000000000000000000-1","-1"); # Kludge needed due to rjust
	return lstResult;


#//------------------------------------------------------------------------------
# ALL Process initialisation
#//------------------------------------------------------------------------------


# Parse input parameters
parser = OptionParser()

parser.add_option("-g", "--gtf", action="store", type="string", dest="GTF", default='ILLUMINA', help="path to GTF annotation file");
parser.add_option("-s", "--sam",     action="store", type="string", dest="SAM",   help="path to SAM reads file (single end reads)");
parser.add_option("-w", "--wrk", action="store", type="string", dest="WRK", help="path to working folder where computation is performed");
parser.add_option("-o", "--output", action="store", type="string", dest="OUT", help="path to folder where Hercules-report.html and All-fourmers.txt are written");
parser.add_option("-r", action="store_true", dest="REMOVE_INTERMEDIATES", default=False);




(options, args) = parser.parse_args()

if len(sys.argv) == 1:
    if (rank == 0):
        parser.print_help()
    exit(0);
else:    
    CONF_LFS_GTF_FILE_FILTERED_ = options.GTF;
    CONF_LFS_SAM_READS_ = options.SAM;
    CONF_LFS_WORKING_ = os.path.join(options.WRK, '');
    CONF_REMOVE_INTERMEDIATES_ = options.REMOVE_INTERMEDIATES;

    print "Using the following parameters:";
    print "GTF path: {}".format(CONF_LFS_GTF_FILE_FILTERED_);
    print "SAM path: {}".format(CONF_LFS_SAM_READS_);
    print "Working folder: {}".format(CONF_LFS_WORKING_);

#//------------------------------------------------------------------------------
# MASTER process
#//------------------------------------------------------------------------------

if rank == 0:
	#data = [(x+1)**x for x in range(size)];
	#print "scattering {}".format(data);

	#result = scatter(data);
	#result = doScatterWork(result);
	#print "RANK0, ",result;
	#lstReceived.append(result);

	for i in range(1, size):
		ping(i);

	# initialise nodes
	for i in range(1, size):
		initialise(i);
		
	gtfsamFile = CONF_LFS_WORKING_ + "GTF-SAM.final.dat";
	
	# If GTF-SAM file exists, don't perform the GTF-SAM map step
	# NB: to allow kick-off from after GTF-SAM map step
	if not os.path.exists(gtfsamFile):

		allocWorkGTFSAM();

		# unload GTF caches to free up memory on nodes
		for i in range(1, size):
			unloadGTFcache(i);

		# wait for ALL nodes to complete before combining GTF-SAM map steps
		syncWait();
		combineWorkGTFSAM(); # SERIALISED
		# combineWorkGTFSAM() is serialised so will have already waited for all workers to complete

	else:
		pprint(">> Starting off from previously performed GTF-MAP step");
		
	#syncWait();

	for a in bases:
		for b in bases:
			for c in bases:
				for d in bases:
					fourmer=a+b+c+d;

					allocWorkMotifMap(fourmer);

					# wait for ALL nodes to complete before combining MOTIF map steps
					syncWait();
					combineMotifMap(fourmer); # SERIALISED

					# combineWorkMotif() is serialised so will have already waited for all workers to complete
					
					allocWorkMotifReduce(fourmer);
					# wait for ALL nodes to complete before combining GTF-SAM map steps
					syncWait();
					
					combineMotifReduce(fourmer); # SERIALISED

					# combineWorkReduce() is serialised so will have already waited for all workers to complete


	# Perform GC content computation	
	for i in range(1, size):
		allocWorkSAMGC(i);

	# Aggregate GC content computations
	GCDatFile = CONF_LFS_WORKING_ + "GC.final.dat";
	finalGCcsv = CONF_LFS_WORKING_ + "GC-content.csv";

	combineGCMap(); # SERIALISED

	lstGCFinal = computeFinalGCcontent(GCDatFile);
	saveText(finalGCcsv, lstGCFinal);
	

	# retire nodes
	for i in range(1, size):
		retire(i);



#//------------------------------------------------------------------------------
# WORKER process
#//------------------------------------------------------------------------------

if rank <> 0:
	#data = None;
	#result = scatter(data);
	#result = doScatterWork(result);
	#print 'result',result;
	## return result to master
	#comm.send(result, dest=0, tag=tag);

	# vars for worker process(es)
	lstGTF_feat = [];
	lstGTF_chr = [];
	bWorkerRetired = False;


	while not bWorkerRetired:

		data = comm.recv(source=0, tag=MPI.ANY_TAG, status=status);
		source = status.Get_source();
		tag = status.Get_tag();

		# handle the types of messages we can recieve and perform their tasks
		if (tag == tags.PING):
			pprint("received ping from Master process, pinging back!");
			comm.send(0, dest=0, tag=tags.PING);

		if (tag == tags.INITIALISE):
			pprint("received the command to initialise.");


			# SORT ON "CHR" (str) AND "FEATURE START" (num) IS CRITICAL FOR GTF LOOK-UP'S
			# BINARY SEARCH TREE ALGORITHM! Linux sort command supports in place sorting,
			# if output (-o) is specified is the same as the input file it creates a temporary file
			# and then copies  back to the given filename.
			os.popen("sort -k1,1 -k4,4n " + CONF_LFS_GTF_FILE_FILTERED_ + " -o " + CONF_LFS_GTF_FILE_FILTERED_);

			# preload the GTF annotation
			fr = open(CONF_LFS_GTF_FILE_FILTERED_, 'r');	
			gtf_lines = fr.readlines();			

			# Load data for Binary Search method
			cur_chr = "";
			for x in range(0, len(gtf_lines)):
				chr, source, feature, feat_start, feat_end, score, strand, frame, attribute = gtf_lines[x].strip().split('\t');
				gtf_key = feat_start.zfill(10) + feat_end.zfill(10);
				lstGTF_feat.append([chr, int(feat_start), int(feat_end), gtf_key]);

				# Optimisation
				# Build list of chromosome ranges
				if (chr <> cur_chr):
						cur_chr = chr;
						lstGTF_chr.append([chr, x, -1]);
						if (x > 0):
							lstGTF_chr[len(lstGTF_chr) - 2][2] = x - 1;
			lstGTF_chr[len(lstGTF_chr) - 1][2] = len(gtf_lines) - 1;

			# Optimisation
			# Pre-sort features within chromosomes
			lstGTF_featTmp = [];
			for chromes in lstGTF_chr:
				lstTmp = lstGTF_feat[chromes[1]:chromes[2]+1];
				lstGTF_featTmp.extend( sorted(lstTmp, key=lambda x: (x[0], x[1])) );
			lstGTF_feat = [];
			lstGTF_feat.extend(lstGTF_featTmp);
			del lstGTF_featTmp;

			pprint("all GTF-cached up and ready to go!");

		if (tag == tags.UNLOADGTF):
			pprint("received the command to unload GTF cache.");
			del lstGTF_feat;
			del lstGTF_chr;
			

		if (tag == tags.WORK):
			pprint("received the command to do work!");
			#print "work instructions = ", data;

			wrkFile = CONF_LFS_WORKING_ + "GTF-SAM." + str(rank) + ".dat";
			wrkFileGC = CONF_LFS_WORKING_ + "GC." + str(rank) + ".dat";
			
			if (data[0] == tasks.GTFSAM):
				lstWrkResult = [];
				workParams = data[1];
				pprint("performing GTF-SAM partitioning on read chunk #: {}".format(workParams));
				workChunk = readlinesFileSection(CONF_LFS_SAM_READS_, workParams[0], workParams[1] - 1);

				for x in range(0, len(workChunk)):
					# Partition the read
					lstWrkResult.append(gtf_sam_mapper(workChunk[x]));
					if (x <> 0) and (x % len(workChunk) / 10 == 0):
						pprint("{}% completed".format(x * 10));

				pprint("100% completed. job done.");				
				saveText(wrkFile, lstWrkResult);
				del lstWrkResult;
				
				syncMaster();

			if (data[0] == tasks.GTFSAMCOMBINE):
				pprint("aggregating GTF-SAM data");
				pprint("input chunk: {}".format(wrkFile));

				if os.path.exists(wrkFile):
					
					gtfsamFile = CONF_LFS_WORKING_ + "GTF-SAM.final.dat";

					# filter out the un-partitionable reads
					os.popen("grep -v '\-1' " + wrkFile + " >>" + gtfsamFile);
					pprint("written GTF-SAM step file: {}".format(gtfsamFile));

					# remove wkrFile from LFS
					if CONF_REMOVE_INTERMEDIATES_:
						os.popen("rm -f " + wrkFile);
				else:
					pprint("error, input file doesn't exist {}".format(wrkFile));

				sendSerialised();

			if (data[0] == tasks.MOTIFMAP):
				lstWrkResult = [];
				lstTemp = []; # work with MOTIF map data directly without writing to disk
							  # allows performaning of MOTIF and VECTOR map steps together
							  # sequentially on each worker
				workParams = data[1];
				aMotif = workParams[0];
				wrkFile = CONF_LFS_WORKING_ + "MOTIF." + aMotif + "." + str(rank) + ".dat";
				pprint("performing MOTIF map on chunk of reads: {}".format(workParams));

				gtfsamFileMotif = CONF_LFS_WORKING_ + "GTF-SAM." + aMotif + ".dat";
				workChunk = readlinesFileSection(gtfsamFileMotif, workParams[1][0], workParams[1][1] - 1);

				# perform MOTIF MAP step
				for x in range(0, len(workChunk)):				

					mapStr = "";
					try:
						mapStr = motif_mapper(workChunk[x], aMotif);
						lstTemp.append(mapStr);
					except:
						pprint("encountered a partitioned read that could not be processed");
						# log error
						continue;
					
					if (x <> 0) and (x % len(workChunk) / 10 == 0):
						pprint("{}% completed".format(x * 10));


				pprint("100% completed. job done.");

				pprint("performing VECTOR map for chunk of data for motif {} ".format(aMotif));

				# perform VECTOR MAP step
				for x in range(0, len(lstTemp)):

					#mapStr = "";
					#try:
					lstVectorMapResult = vector_mapper(lstTemp[x]);

					for i in lstVectorMapResult:
						lstWrkResult.append(i);

					#except:
					#	pprint("encountered a MOTIF position tuple that could not be vector-mapped");
						# log error
					#	continue;
					
					if (x <> 0) and (x % len(lstTemp) / 10 == 0):
						pprint("{}% completed".format(x * 10));

				pprint("100% completed. job done.");
				
				saveText(wrkFile, lstWrkResult);
				del lstWrkResult;
				
				syncMaster();

			if (data[0] == tasks.MOTIFMAPCOMBINE):
				lstWrkResult = [];
				aMotif = data[1];
				wrkFile = CONF_LFS_WORKING_ + "MOTIF." + aMotif + "." + str(rank) + ".dat";
				pprint("aggregating MOTIF MAP data");
				pprint("input chunk: {}".format(wrkFile));

				motifmapFile = CONF_LFS_WORKING_ + "MOTIF." + aMotif + ".final.dat";
				if os.path.exists(wrkFile):
					
					os.popen("cat " + wrkFile + " >>" + motifmapFile);
					pprint("written MOTIF map step file: {}".format(motifmapFile));

					if CONF_REMOVE_INTERMEDIATES_:
						os.popen("rm -f " + wrkFile);

				else:
					pprint("error, input file doesn't exist {}".format(wrkFile));

				sendSerialised();

			if (data[0] == tasks.MOTIFREDUCE):
				lstWrkResult = [];
				workParams = data[1];
				aMotif = workParams[0];
				motifmapFile = CONF_LFS_WORKING_ + "MOTIF." + aMotif + ".final.dat";
				wrkFile = CONF_LFS_WORKING_ + "MOTIF.REDUCE." + aMotif + "." + str(rank) + ".dat";
				pprint("performing MOTIF-REDUCE on chunk of mapped reads: {}".format(workParams));
				workChunk = readlinesFileSection(motifmapFile, workParams[1][0], workParams[1][1] - 1);

				lstWrkResult = reduceTupleList(workChunk, True); # The reduce operation needs to sort it's chunk
				
				pprint("100% completed. job done.");				
				saveText(wrkFile, lstWrkResult);
				del lstWrkResult;
				
				syncMaster();

			if (data[0] == tasks.MOTIFREDUCECOMBINE):
				lstWrkResult = [];
				aMotif = data[1];
				wrkFile = CONF_LFS_WORKING_ + "MOTIF.REDUCE." + aMotif + "." + str(rank) + ".dat";
				pprint("aggregating MOTIF REDUCE data");
				pprint("input chunk: {}".format(wrkFile));

				motifreduceFile = CONF_LFS_WORKING_ + "MOTIF.REDUCE." + aMotif + ".final.dat";
				if os.path.exists(wrkFile):
					
					os.popen("cat " + wrkFile + " >>" + motifreduceFile);
					pprint("written MOTIF map step file: {}".format(motifreduceFile));

					if CONF_REMOVE_INTERMEDIATES_:
						os.popen("rm -f " + wrkFile);

				else:
					pprint("error, input file doesn't exist {}".format(wrkFile));

				sendSerialised();

			if (data[0] == tasks.GCMAP):
				lstWrkTmp = [];
				lstWrkResult = [];
				workParams = data[1];
				pprint("Computing mean exon GC content for chunk of reads: {}".format(workParams));
				workChunk = loadText(wrkFile);

				for c in workChunk:
					lstWrkTmp.append(gc_mapper(c));
				
				# Compute intermediate GC results
				lstWrkResult = computeGCsums(lstWrkTmp);

				pprint("100% completed. job done.");				
				saveText(wrkFileGC, lstWrkResult);
				del lstWrkTmp, lstWrkResult;
				
				syncMaster();
			
			if (data[0] == tasks.GCMAPCOMBINE):
				lstWrkResult = [];
				pprint("aggregating GC MAP data");
				pprint("input chunk: {}".format(wrkFileGC));

				GCDataFile = CONF_LFS_WORKING_ + "GC.final.dat";
				if os.path.exists(wrkFileGC):
					
					os.popen("cat " + wrkFileGC + " >>" + GCDataFile);
					pprint("written GC map step file: {}".format(GCDataFile));

					#if CONF_REMOVE_INTERMEDIATES_:
					#	os.popen("rm -f " + wrkFileGC);

				else:
					pprint("error, input file doesn't exist {}".format(wrkFile));

				sendSerialised();

		if (tag == tags.RETIRE):
			pprint("recieved the command to retire!")
			bWorkerRetired = True;
			exit(0); # happily retire


		sleep(0.5);
