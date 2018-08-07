#//==============================================================================
#//    _    _                     _           
#//   | |  | |				   | |		  
#//   | |__| | ___ _ __ ___ _   _| | ___  ___ 
#//   |  __  |/ _ \ '__/ __| | | | |/ _ \/ __|
#//   | |  | |  __/ | | (__| |_| | |  __/\__ \
#//   |_|  |_|\___|_|  \___|\__,_|_|\___||___/ MPI
#//
#//==============================================================================
# Hercules-mpi: MPI implementation of Hercules transcriptomics analysis
# Dr. Jamie Alnasir
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
 
 
# ** BEFORE RUNNING: MUST EDIT HerculesConf with configration **
 
 
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
import time;
import random;
import pickle;
import datetime;
  
# Stats libraries
import numpy as np
from scipy import stats;
from scipy.special import stdtr;
 
#//------------------------------------------------------------------------------
 
#// Configuration parameters and input/output file paths MUST be specified in
#// Hercules-mpi-conf.py
 
#//------------------------------------------------------------------------------
 
 
# MPI initialization
comm = MPI.COMM_WORLD;		  # get MPI communicator object
size = comm.Get_size();		 # total number of processes
rank = comm.Get_rank();		 # rank of this process
name = MPI.Get_processor_name();	# usually returns hostname
status = MPI.Status();		  # get MPI status object
 
 
# Fixed configuration parameters (not specified using Option Parser)
CONF_DISK_IO_WAIT_ = 0.300; # 300 ms
 
# Configuration parameters (now set using Option Parser)
CONF_LFS_WORKING_		   = "";
CONF_LFS_GTF_FILE_FILTERED_ = "";
CONF_LFS_SAM_READS_		 = "";
CONF_LFS_OUT_			   = "";
REPORT_HTML				 = "Hercules-report.html";
REPORT_TXT				  = "All-fourmers.txt";
 
 
 
# Nucleotide bases
bases = "ATGC";

# For exon reads median GC content lookup
dictGC = {};

# Enable for Filtering reads/counts based on exon % GC content
_GC_FILTERING_ = False;
_GC_FILTER_MIN_ = 60;
_GC_FILTER_MAX_ = 70;

_MIN_CORRELATION_COUNTS_ = 10; # NORMALLY 10
_CORREL_STAT_SPEARMAN_ = False;
FIRST_QUARTILE = 1;
 
def enum(*sequential, **named):
	"""Handy way to fake an enumerated type in Python
	http://stackoverflow.com/questions/36932/how-can-i-represent-an-enum-in-python
	"""
	enums = dict(zip(sequential, range(len(sequential))), **named)
	return type('Enum', (), enums)
 
# Define MPI message tags
tags = enum('PING', 'WORK', 'INITIALISE', 'UNLOADGTF', 'RETIRE', 'SYNC', 'SERIALISED');
 
# Hercules work tags
tasks = enum('Q1', 'MOTIFCORREL', 'PEARSONTABLE1', 'PEARSONTABLE2', 'BOXPLOT');
 
 
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
 
 
def allocWork1Q():
# send out computation to the nodes: calc 1st Inter-quartile range from the supplied motif file
	 
	pprint("Computing 1st Quartiles from motif files ".format(CONF_LFS_WORKING_));
 
	motifCount = 256;
	chunkMotifs = motifCount / (size - 1);
	pprint("Number of motifs =".format(motifCount));
	pprint("allocating {} motifs per worker process".format(chunkMotifs));
	params = chunkRanges(0, motifCount, chunkMotifs);
 
 
	for i in range(1, size):
		#print "SENDING WORK" + str(i);
		#print params;
		sendWork(i, tasks.Q1, params[i - 1]);	   
 
 
 
 
def allocWorkMotifCorrel():
# send out computation to the nodes: calc 1st Inter-quartile range from the supplied motif file
 
	pprint("Computing Correlations from motif files ".format(CONF_LFS_WORKING_));
 
	motifCount = 256;
	chunkMotifs = motifCount / (size - 1);
	pprint("Number of motifs =".format(motifCount));
	pprint("allocating {} motifs per worker process".format(chunkMotifs));
	params = chunkRanges(0, motifCount, chunkMotifs);
 
	for i in range(1, size):
		#print "SENDING WORK" + str(i);
		#print params;
		sendWork(i, tasks.MOTIFCORREL, params[i - 1]);	  
 
 
 
#//------------------------------------------------------------------------------
# Hercules RNA-Seq analysis functions
#//------------------------------------------------------------------------------
 
def prettyFloat(afloat):
	return "{0:.4f}".format(afloat);
 
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
 
 
def saveObject(objFile, obj):
# Save object using Pickle
	with open(objFile, 'wb') as f:
		pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)
 
def loadObject(objFile):
# Load object using Pickle
	with open(objFile, 'rb') as f:
		return pickle.load(f)
 
 
def loadText(aInFile):
# Load array of strings from text file
	with open(aInFile) as f:
		data = f.readlines();
		return data;
 
def saveText(aOutFile, data):
# Save array of strings to text file
	with open(aOutFile, "w") as f:
		for i in data:
			f.write(i + "\n");
 
 
def Q1(lstV):
# Use numpy to calculate Q1 (q25)   
	q75, q25 = np.percentile(lstV, [75 ,25], interpolation = 'higher');
	#return 6; # INVESTIGATION 21/03/2017 -- TEST;
	return q25;
 
def calcFirstQuartile(sorted_motifs_file):
# Calculate the First Quartile (value between lowest value and median) for ALL of the counts EXCEPT
# those which were not associated with an exon (have Exon key of -1 from the Hercules Map-Reduce Job).
	f = open(sorted_motifs_file, 'r');
	motifs_raw = f.read().splitlines();
	lstCounts = [];
	for motif_line in motifs_raw:
		exon, pos, count = motif_line.split(',');
		try:
			# Only interested in
			if (exon <> "-1"):				
				lstCounts.append(int(count));
		except:
			(1 == 1);
			#print "error: " + motif_line;
	#c_avg = float( sum(lstCounts) ) / float( len(lstCounts) );
	# Calc 1st quartile
	c_fq = Q1(lstCounts);
	return c_fq;
 
 
def _safeCalcCorrelStat(lstExonCounts, lstCounts1, lstCounts2):
# Calculate CorrelStat Co-efficient
	if (len(lstExonCounts) >= _MIN_CORRELATION_COUNTS_):
		try:
			if _CORREL_STAT_SPEARMAN_:
				Rvalue = stats.spearmanr(lstCounts1, lstCounts2)[0];
			else:
				Rvalue = np.corrcoef(lstCounts1, lstCounts2)[0, 1];		
		except Exception, e:
			Rvalue = 0;
			print "NUMPY ERR:", len(lstCounts1);
			print str(e);
	else:
		Rvalue = 0;		
	return Rvalue;
 
 
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
 
 
def removeNoise(exon_motif_reads):
	if (exon_motif_reads == []):
		return [];
  
	global _CAP_COUNT_;
	global FIRST_QUARTILE;
	lstWorking = [];
	counts = [];
	#print "FIRST QUARTILE = ", FIRST_QUARTILE;
	# itertools group iterator, pull off array
	for item in exon_motif_reads:
		counts.append(item[2]);
		lstWorking.append([item[0], item[1], item[2]]);
  
	# USE FIRST_QUARTILE GLOBAL CALCULATED WITH NUMPY
	#return filter(lambda x: (x[2] >= FIRST_QUARTILE and x[2] <= _CAP_COUNT_), lstWorking); # REMOVED CAP
	return filter(lambda x: (x[2] >= FIRST_QUARTILE), lstWorking);
 
 
def _motifsBySpace(sorted_motifs_file, space_bp, tolerance_bp):
# Return motif occurrences on the same exon within a given spacing
	global dictMotifs;  
	f = open(sorted_motifs_file, 'r');
	motifs_raw = f.read().splitlines(); 
	lstMotifs = [];
	lstResults = [];
	for motif_line in motifs_raw:
		exon, pos, count = motif_line.split(',');
		try:
			lstMotifs.append([exon, int(pos), int(count)]);
		except:
			(1 == 1);
			#print "error: " + motif_line;
	# Filter reads/counts based on exon GC content
	if _GC_FILTERING_:
		#print "GC-FILTERING OF counts for the motif by EXON GC:", _GC_FILTER_MIN_, _GC_FILTER_MAX_;
		lstMotifs = filterByGC(lstMotifs, _GC_FILTER_MIN_, _GC_FILTER_MAX_);
	for key, group in groupby(lstMotifs, lambda x: x[0]):
		if (key == "-1"):
			continue;
		filtered_group = removeNoise(group);
		filtered_group2 = pairsWithinSpace(filtered_group, space_bp, tolerance_bp);
		for item in filtered_group2:
			#print item[0] + ", " + str(item[1]) + ", " + str(item[2]) + ", " + item[3] + ", " + str(item[4]) + ", " + str(item[5]);
			lstResults.append(item);
	return lstResults;
 
 
def pairsWithinSpace(exon_motif_reads, space_bp, tolerance_bp):
	# Filter exons with motifs occurring within a given spacing (+/- tolerance)
	# Value copy exon lists for two separate filtering operations
	lstExons1 = list(exon_motif_reads);
	lstExons2 = list(lstExons1); 
	lstResult = [];
	for i in range (0, len(lstExons1)):
		for j in range (i+1, len(lstExons2)):
			#print lstExons1[i], lstExons2[j];
			if (abs(lstExons2[j][1] - lstExons1[i][1]) >= space_bp-tolerance_bp) and (abs(lstExons2[j][1] - lstExons1[i][1]) <= space_bp+tolerance_bp):
					lstResult.append(lstExons1[i] + lstExons2[j]);	 
	return lstResult;
 
 
def _calcCorrelStat(motif, motifFile, bPrint):
# Compute CorrelStat correlation co-efficients for given fourmer
# return an array of motif and correlation co-efficients: [motif, Rvalue10, Rvalue50, Rvalue100, Rvalue200];
# _MIN_CORRELATION_COUNTS_ is used to assert a minimum number of reads to compute the correlation, typically about 10.	
	dictRval = {10:0, 50:0, 100:0, 200:0}
	dictCorrelsCount = {10:0, 50:0, 100:0, 200:0}
	for spacing in [10,50,100,200]:
		if spacing in [10,50]:
			tol = 2;
		else:
			tol = 4;
		dictSp = {spacing: []}; # for holding multiple exons CorrelStat Correlations for each spacing.
								# NB: these will be averaged for each spacing.
		# n-spaced n-mers by exon
		lstMotifs = _motifsBySpace(motifFile, spacing, tol);
		lstCounts1 = getCol(lstMotifs, 2);
		lstCounts2 = getCol(lstMotifs, 5);		
		dictCorrelsCount[spacing] = len(lstCounts1);					
		if len(lstCounts1) >= _MIN_CORRELATION_COUNTS_: # we need at least this many pairs to compute correlation
			dictRval[spacing] = _safeCalcCorrelStat(lstCounts1, lstCounts1, lstCounts2);			
	if (bPrint):
		print motif + ", " + "{0:.3f}".format(FIRST_QUARTILE) + ", " + str(dictRval[10]) + " (" + str(dictCorrelsCount[10]) + "), " + str(dictRval[50]) + " (" + str(dictCorrelsCount[50]) + "), " + str(dictRval[100])  + " (" + str(dictCorrelsCount[100]) + "), " + str(dictRval[200]) + " (" + str(dictCorrelsCount[200]) + ")";
	return [motif, [dictRval[10], dictRval[50], dictRval[100], dictRval[200]], [dictCorrelsCount[10], dictCorrelsCount[50], dictCorrelsCount[100], dictCorrelsCount[200]] ];
 
 
def getMotifDictAsFlatArray():
	lstResult = []
	for k in dictCorrels.keys():
		item = dictCorrels[k][0] + dictCorrels[k][1];
		item.insert(0, k);
		lstResult.append(item);
	return lstResult;


def filterByGC(exon_motif_reads, minGC, maxGC):
# Filter exon counts by their GC content
		lstWorking = [];
		for item in exon_motif_reads:
				lstWorking.append([item[0], item[1], item[2]]);

		#print len(lstWorking), minGC, maxGC;

		# return only exon counts for exons with median GC content within minGC and maxGC
		lstR =  filter(lambda x: (GC_Lookup(x[0]) >= minGC and GC_Lookup(x[0]) < maxGC), lstWorking);
		#print len(lstR);
		return lstR;


def LoadGC(gc_file):
		global dictGC;
		fr = open(gc_file, 'r');
		lstGC = fr.read().splitlines();
		for i in lstGC:
				exon, gc_str = i.split(',');
				dictGC[exon] = float(gc_str);

def GC_Lookup(exonIDStr):
		return dictGC[exonIDStr];



def rawCorrelStatData(motif, motifFile, bPrint):
# Return array of simple tuples for motif for Mean ExonGC, MotifGC, Spacing, R-value
# There will be 4 tuples in the result, one for each spacing.

		global _GC_FILTERING_;
		global _GC_FILTER_MIN_, _GC_FILTER_MAX_; # to be writable to set filter range

		# We'll be cycling the filter range for the computations
		_GC_FILTERING_=True;

		# GC content bins
		lstREAD_GC_bins  = [[30,40], [40,50], [50,60], [60,70]];

		lstRawCorrelData = [];


		dictRval = {10:0, 50:0, 100:0, 200:0}
		for spacing in [10,50,100,200]:

				if spacing in [10,50]:
						tol = 2;
				else:
						tol = 4;

				# filter by mean GC content of the exon reads
				for gc_range in lstREAD_GC_bins:

						_GC_FILTER_MIN_ = gc_range[0];
						_GC_FILTER_MAX_ = gc_range[1];

						bin_index_exon_gc = gc_range[0];
						#print "Filtering for Median Exon GC content range: ", gc_range[0], '-', gc_range[1];
						MeanExonGC = (float(gc_range[1] - gc_range[0]) / 2) + gc_range[0];

						# n-spaced n-mers by exon
						#print "Mean GC=", MeanExonGC, " Motif GC=", GC_Content(motif), " spacing=", spacing;
						lstMotifs = _motifsBySpace(motifFile, spacing, tol); # NB filtered by Read/Exon GC

						lstCounts1 = getCol(lstMotifs, 2);
						lstCounts2 = getCol(lstMotifs, 5);


						if len(lstCounts1) >= _MIN_CORRELATION_COUNTS_: # we need at least this many pairs to compute correlation

								# Test 29/03/2017 -- use real average exon GC instead of middle of bin
								#meanExonGC = meanGC(getCol(lstMotifs, 0))
								#print "meanGC = ", meanExonGC;


								Rvalue = _safeCalcCorrelStat(lstCounts1, lstCounts1, lstCounts2);
								lstRawCorrelData.append([MeanExonGC, GC_Content(motif), spacing, Rvalue]);
								#lstRawCorrelData.append([meanExonGC, GC_Content(motif), spacing, Rvalue]);

	_GC_FILTERING_ = False;

		return lstRawCorrelData;


def saveRawCorrelationsData(aFile, lstAllCorrels):
		lines = []
		for i in lstAllCorrels:
				lineStr = ",".join(map(str, i));
				lines.append(lineStr);
		saveText(aFile, lines);


def printCorrelStat():
 
	OutlierCount = 10;
	COL_R10 = 1;
	COL_MOTIF = 0;
	COL_R10   = 1;
	COL_R50   = 2;
	COL_R100  = 3;
	COL_R200  = 4;
	CCOUNT_OFFSET = 4;
 
	lstCorrelTbl = getMotifDictAsFlatArray();
	tmp = sorted(lstCorrelTbl, key=lambda x: x[COL_R10], reverse=False)[:OutlierCount];
	lstR10min = map(lambda a,b,c:[a[COL_MOTIF],b[COL_R10],c[COL_R10 + CCOUNT_OFFSET]], tmp, tmp, tmp);
 
	lstCorrelTbl = getMotifDictAsFlatArray();
	tmp = sorted(lstCorrelTbl, key=lambda x: x[COL_R50], reverse=False)[:OutlierCount];
	lstR50min = map(lambda a,b,c:[a[COL_MOTIF],b[COL_R50],c[COL_R50 + CCOUNT_OFFSET]], tmp, tmp, tmp);
 
	lstCorrelTbl = getMotifDictAsFlatArray();
	tmp = sorted(lstCorrelTbl, key=lambda x: x[COL_R100], reverse=False)[:OutlierCount];
	lstR100min = map(lambda a,b,c:[a[COL_MOTIF],b[COL_R100],c[COL_R100 + CCOUNT_OFFSET]], tmp, tmp, tmp);
 
	lstCorrelTbl = getMotifDictAsFlatArray();
	tmp = sorted(lstCorrelTbl, key=lambda x: x[COL_R200], reverse=False)[:];
	lstR200min = map(lambda a,b,c:[a[COL_MOTIF],b[COL_R200],c[COL_R200 + CCOUNT_OFFSET]], tmp, tmp, tmp);
 
 
 
	lstCorrelTbl = getMotifDictAsFlatArray();
	tmp = sorted(lstCorrelTbl, key=lambda x: x[COL_R10], reverse=True)[:OutlierCount];
	lstR10max = map(lambda a,b,c:[a[COL_MOTIF],b[COL_R10],c[COL_R10 + CCOUNT_OFFSET]], tmp, tmp, tmp);
 
	lstCorrelTbl = getMotifDictAsFlatArray();
	tmp = sorted(lstCorrelTbl, key=lambda x: x[COL_R50], reverse=True)[:OutlierCount];
	lstR50max = map(lambda a,b,c:[a[COL_MOTIF],b[COL_R50],c[COL_R50 + CCOUNT_OFFSET]], tmp, tmp, tmp);
 
	lstCorrelTbl = getMotifDictAsFlatArray();
	tmp = sorted(lstCorrelTbl, key=lambda x: x[COL_R100], reverse=True)[:OutlierCount];
	lstR100max = map(lambda a,b,c:[a[COL_MOTIF],b[COL_R100],c[COL_R100 + CCOUNT_OFFSET]], tmp, tmp, tmp);
 
	lstCorrelTbl = getMotifDictAsFlatArray();
	tmp = sorted(lstCorrelTbl, key=lambda x: x[COL_R200], reverse=True)[:OutlierCount];
	lstR200max = map(lambda a,b,c:[a[COL_MOTIF],b[COL_R200],c[COL_R200 + CCOUNT_OFFSET]], tmp, tmp, tmp);
 
 
	lstHTMLReport.append("<br><br>");
	lstHTMLReport.append("<h2>Fourmer Outliers</h2>");
 
	lstHTMLReport.append("<table border=\"1\" cellpadding=\"3\" cellspacing=\"0\">");
	lstHTMLReport.append("<tr><td colspan=\"4\"><b>Lowest Outliers</b></td></tr>");
	lstHTMLReport.append("<tr><td><b>r<sup>2</sup> (10 bp)</b></td><td><b>r<sup>2</sup> (50 bp)</b></td><td><b>r<sup>2</sup> (100 bp)</b></td><td><b>r<sup>2</sup> (200 bp)</b></td></tr>");
 
	print "Lowest Outliers";
	print "r^2(10 bp), r^2(50 bp), r^2(100 bp), r^2(200 bp)";
 
	for x in range(0,OutlierCount):
		COL_VAL = 1;
		COL_CC  = 2;
		cc10str = " ({})".format(lstR10min[x][COL_CC]);
		cc50str = " ({})".format(lstR50min[x][COL_CC]);
		cc100str = " ({})".format(lstR100min[x][COL_CC]);
		cc200str = " ({})".format(lstR200min[x][COL_CC]);
		line10Str = lstR10min[x][COL_MOTIF]  + "=" + prettyFloat(lstR10min[x][COL_VAL]) + cc10str
		line50Str = lstR50min[x][COL_MOTIF]  + "=" + prettyFloat(lstR50min[x][COL_VAL]) + cc50str;
		line100Str = lstR100min[x][COL_MOTIF]  + "=" + prettyFloat(lstR100min[x][COL_VAL]) + cc100str;
		line200Str = lstR200min[x][COL_MOTIF]  + "=" + prettyFloat(lstR200min[x][COL_VAL]) + cc200str;
		print line10Str + "  " + line50Str + "  " + line100Str + "  " + line200Str;
 
		lstHTMLReport.append("<tr><td>" + line10Str + "</td><td>" + line50Str + "</td><td>" + line100Str + "</td><td>" + line200Str + "</td></tr>");
 
 
	print "";
	print "Highest Outliers";
	print "r^2(10 bp), r^2(50 bp), r^2(100 bp), r^2(200 bp)";
 
 
	lstHTMLReport.append("<tr><td colspan=\"4\"><b>Highest Outliers</b></td></tr>");
	lstHTMLReport.append("<tr><td><b>r<sup>2</sup> (10 bp)</b></td><td><b>r<sup>2</sup> (50 bp)</b></td><td><b>r<sup>2</sup> (100 bp)</b></td><td><b>r<sup>2</sup> (200 bp)</b></td></tr>");
 
	for x in range(0,OutlierCount):
		COL_VAL = 1;
		COL_CC  = 2;
		cc10str = " ({})".format(lstR10max[x][COL_CC]);
		cc50str = " ({})".format(lstR50max[x][COL_CC]);
		cc100str = " ({})".format(lstR100max[x][COL_CC]);
		cc200str = " ({})".format(lstR200max[x][COL_CC]);
		line10Str = lstR10max[x][COL_MOTIF]  + "=" + prettyFloat(lstR10max[x][COL_VAL]) + cc10str
		line50Str = lstR50max[x][COL_MOTIF]  + "=" + prettyFloat(lstR50max[x][COL_VAL]) + cc50str;
		line100Str = lstR100max[x][COL_MOTIF]  + "=" + prettyFloat(lstR100max[x][COL_VAL]) + cc100str;
		line200Str = lstR200max[x][COL_MOTIF]  + "=" + prettyFloat(lstR200max[x][COL_VAL]) + cc200str;
		print line10Str + "  " + line50Str + "  " + line100Str + "  " + line200Str;
 
		lstHTMLReport.append("<tr><td>" + line10Str + "</td><td>" + line50Str + "</td><td>" + line100Str + "</td><td>" + line200Str + "</td></tr>");
 
	lstHTMLReport.append("</table>");
 
 
def getHTMLStart():
	DateTimeStr = datetime.date.today().strftime("%B %d, %Y");
	return """<body>
<h1>Hercules - Fourmer Motif Analysis Report</h1>
<br>
{dt}<br>
GTF Annotation: {gtf}<br>
SAM reads: {sam}<br>
<br>
Refer to the following papers:
<ul>
  <li>
	<p style=\"line-height: 150%\" align=\"left\">Alnasir, J.  &amp; Shanahan, H. P. 
	(2017). <b> Transcriptomics: Quantifying non-uniform read distribution using MapReduce.</b>.
	<i>International Journal of Foundations of Computer Science (Forthcoming)</i>.</p>
  </li>
 
  <li>
	<p style=\"line-height: 150%\" align=\"left\">Alnasir, J., &amp; Shanahan, H. 
	(2017). <b>A novel method to detect bias in Short Read NGS RNA-seq data</b>.
	<i>Journal of Integrative Bioinformatics, 14(3).</i></p>
  </li>
</ul>
 
<h2>All fourmer motifs</h2>""".format(dt=DateTimeStr, gtf=CONF_LFS_GTF_FILE_FILTERED_, sam=CONF_LFS_SAM_READS_);
 
def getHTMLEnd():
	return """</body>""";
 
 
def doReporting():
 
	lstHTMLReport.append(getHTMLStart());
	lstHTMLReport.append("<table width=\"880\" border=\"1\" cellpadding=\"3\" cellspacing=\"0\">");
	lstHTMLReport.append("<tr><td><b>Motif</b></td><td><b>Q1</b></td><td><b>r<sup>2</sup> (10 bp)</b></td><td><b>r<sup>2</sup> (50 bp)</b></td><td><b>r<sup>2</sup> (100 bp)</b></td><td><b>r<sup>2</sup> (200 bp)</b></td></tr>");
 
	lstReportAllFourmers.append("Motif, Q1, r^2(10 bp), r^2(50 bp), r^2(100 bp), r^2(200 bp)");
 
	for fourmer in lstMotifs:	   
		fourmerQ1File = CONF_LFS_WORKING_ + "herc-final-" + fourmer + ".q1";
		fourmerCorrelFile = CONF_LFS_WORKING_ + "herc-final-" + fourmer + ".correl";
 
		# Load 1st Quartile (Q1) from file
		FIRST_QUARTILE = loadObject(fourmerQ1File);
 
		# Load fourmer motif computation results
		correlFourmer = loadObject(fourmerCorrelFile);
 
		# Store in dictionary
		dictCorrels[fourmer] = [correlFourmer[1], correlFourmer[2]];
 
		correls	  = correlFourmer[1];
		correlCounts = correlFourmer[2];
		print fourmer + ",", correls[0], "({}), ".format(correlCounts[0]),  correls[1], "({}), ".format(correlCounts[1]), correls[2], "({}),".format(correlCounts[2]), correls[3], "({}) ".format(correlCounts[3]);
 
 
		correlStr = fourmer + ", " + str(FIRST_QUARTILE) + ", "+ prettyFloat(correls[0]) + " (" + str(correlCounts[0]) + "), " + prettyFloat(correls[1]) + " (" + str(correlCounts[1]) + "), " + prettyFloat(correls[2])  + " (" + str(correlCounts[2]) + "), " + prettyFloat(correls[3]) + " (" + str(correlCounts[3]) + ")";
 
 
		HTMLcorrelStr = "<tr><td>" + correlStr.replace(",", "</td><td>") + "</td></tr>";
 
		lstHTMLReport.append(HTMLcorrelStr);
		lstReportAllFourmers.append(correlStr);
 
 
	lstHTMLReport.append("</table>");
 
	# Save the correlations report
	saveText(CONF_LFS_OUT_ + REPORT_TXT, lstReportAllFourmers);
 
	# Save the correlations dictionary
	saveObject(CONF_LFS_WORKING_ + "ALL-CORRELS.dat", dictCorrels);
 
	print "";
	printCorrelStat();
 
	lstHTMLReport.append(getHTMLEnd());
	saveText(CONF_LFS_OUT_  + REPORT_HTML, lstHTMLReport);
 
 
 
#//------------------------------------------------------------------------------
# ALL Process initialisation
#//------------------------------------------------------------------------------
 
lstMotifs = [];
dictCorrels = {}; # To store ALL of the computed correlations, in the form 'AAAA':[r,r,r,r], [c,c,c,c] -  where r=correl, c=count
lstHTMLReport = []; # To store HTML report
 
lstReportAllFourmers = [];
 
 
# Make a list of the fourmers to perform computations on. This will be processed in blocks
# by the workers.
fourmer="";
for a in bases:
	for b in bases:
		for c in bases:
			for d in bases:
				fourmer=a+b+c+d;
				lstMotifs.append(fourmer);
 
# Parse input parameters
parser = OptionParser()
 
parser.add_option("-g", "--gtf", action="store", type="string", dest="GTF", default='ILLUMINA', help="path to GTF annotation file");
parser.add_option("-s", "--sam",	 action="store", type="string", dest="SAM",   help="path to SAM reads file (single end reads)");
parser.add_option("-w", "--wrk", action="store", type="string", dest="WRK", help="path to working folder where computation is performed");
parser.add_option("-o", "--output", action="store", type="string", dest="OUT", help="path to folder where Hercules-report.html and All-fourmers.txt are written");
 
 
(options, args) = parser.parse_args()
 
if len(sys.argv) == 1:
	if (rank == 0):
		parser.print_help()
	exit(0);
else:	
	CONF_LFS_GTF_FILE_FILTERED_ = options.GTF;
	CONF_LFS_SAM_READS_ = options.SAM;
	CONF_LFS_WORKING_ = os.path.join(options.WRK, '');
	CONF_LFS_OUT_ = os.path.join(options.OUT, '');
 

# Load Exon GC dictionary
LoadGC(CONF_LFS_WORKING_ + "GC-content.csv");

 
#//------------------------------------------------------------------------------
# MASTER process
#//------------------------------------------------------------------------------
 
if rank == 0:
 
	print "Using the following parameters:";
	print "GTF path: {}".format(CONF_LFS_GTF_FILE_FILTERED_);
	print "SAM path: {}".format(CONF_LFS_SAM_READS_);
	print "Working folder: {}".format(CONF_LFS_WORKING_);
	print "Report Output folder: {}".format(CONF_LFS_OUT_);
 
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
 
 
	allocWork1Q();
 
	# wait for ALL nodes to complete
	syncWait();
 
	allocWorkMotifCorrel();
	 
	# wait for ALL nodes to complete
	syncWait();

	lstAllGCcorrels = [];
	for i in range(0, 256):
		fourmer = lstMotifs[i];
		fourmerCorrelFileGC = CONF_LFS_WORKING_ + fourmer + ".gc.correl";

		lstFourmerGCcorrel = loadObject(fourmerCorrelFileGC);
		lstAllGCcorrels.extend(lstFourmerGCcorrel);

				# Remove GC correl file as we're done with it
				#os.popen("rm " + fourmerCorrelFileGC);

	saveObject(CONF_LFS_WORKING_ + "all.gc.correl", lstAllGCcorrels);
	saveRawCorrelationsData(CONF_LFS_OUT_ + "All-correls.csv", lstAllGCcorrels);
 
	# retire nodes
	for i in range(1, size):
		retire(i);
 
 
	doReporting();
 
	# Clear files we're done with
	os.popen("rm -f " + CONF_LFS_WORKING_ + "*.q1");

 
 
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
			pprint("up and ready to go!");
 
		if (tag == tags.UNLOADGTF):
			pprint("received the command to unload GTF cache.");
			del lstGTF_feat;
			del lstGTF_chr;
			 
 
		if (tag == tags.WORK):
			pprint("received the command to do work!");
			#print "work instructions = ", data;
 
			wrkDir = CONF_LFS_WORKING_;
 
			if (data[0] == tasks.Q1):
			# Compute 1st Quartile for range of motifs (fourmers)
				workParams = data[1];
				fourmerChunk = workParams;
 
				startFourmer = fourmerChunk[0];
				endFourmer = fourmerChunk[1] + 1;
 
				pprint("Compute quartiles (Q1) for batch, worker: {}".format(rank));
				pprint("Compute quartiles (Q1), chunk: {}".format(workParams));
 
 
				 
				for i in range(startFourmer, endFourmer):
					fourmer = lstMotifs[i];
					fourmerFile = CONF_LFS_WORKING_ + "herc-final-" + fourmer + ".csv";
					fourmerQ1File = CONF_LFS_WORKING_ + "herc-final-" + fourmer + ".q1";
 
 
					pprint("Computing Q1 for fourmer: {}".format(fourmer));
					motifQ1 = calcFirstQuartile(fourmerFile);
					pprint("Q1 for {} = {}".format(fourmer, motifQ1)); 
 
					# Save result (Q1 to file)				  
					saveObject(fourmerQ1File, motifQ1);
					 
			 
				syncMaster();
 
 
			if (data[0] == tasks.MOTIFCORREL):
			# Compute Motif Correlations (for spacings of 10, 50, 100, 200 bp) for range of motifs (fourmers)
				workParams = data[1];
				fourmerChunk = workParams;
 
				startFourmer = fourmerChunk[0];
				endFourmer = fourmerChunk[1] + 1;
 
				pprint("Compute Correlations for batch, worker: {}".format(rank));
				pprint("Compute Correlations, chunk: {}".format(workParams));
 
				 
				for i in range(startFourmer, endFourmer):
					fourmer = lstMotifs[i];
					fourmerFile = CONF_LFS_WORKING_ + "herc-final-" + fourmer + ".csv";
					fourmerQ1File = CONF_LFS_WORKING_ + "herc-final-" + fourmer + ".q1";
					fourmerCorrelFile = CONF_LFS_WORKING_ + "herc-final-" + fourmer + ".correl";
					fourmerCorrelFileGC = CONF_LFS_WORKING_ + fourmer + ".gc.correl";
 
					# Load 1st Quartile (Q1) from file
					FIRST_QUARTILE = loadObject(fourmerQ1File);
 
					pprint("Computing Correlations for fourmer: {}".format(fourmer));
 
					motifCorrel = _calcCorrelStat(fourmer, fourmerFile, False);
					print motifCorrel;
 
					# Compute GC cycled correlations
					lstMotifCorrels = rawCorrelStatData(fourmer, fourmerFile, False);
 
					# Save results
					# (Q1 to file)				  
					saveObject(fourmerCorrelFile, motifCorrel);
					# Raw GC cycled correls
					saveObject(fourmerCorrelFileGC, lstMotifCorrels);
			
			 
				syncMaster();
 
				 
				 
		if (tag == tags.RETIRE):
			pprint("recieved the command to retire!")
			bWorkerRetired = True;
			exit(0); # happily retire
 
 
		sleep(0.5);
