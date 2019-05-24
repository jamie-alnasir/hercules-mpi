#//==============================================================================
#//    _    _
#//   | |  | |                    _
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
# Correlation vs GC Box plot

import numpy as np
from scipy import stats;
from scipy.special import stdtr;
 
# Plotting libraries
import matplotlib;
matplotlib.use('Agg');
import matplotlib.pyplot as plt;
import matplotlib.gridspec as gridspec;



def getCol(lst, col):
	return [row[col] for row in lst];


def loadText(aInFile):
# Load array of strings from text file
	with open(aInFile) as f:
		data = f.readlines();
		return data;


def getAllCorrels(aCorrelsFile):
# Load data from raw correlations file
    lines = loadText(aCorrelsFile);
    lstData = [];
    for l in lines:
        _exonGC, _motifGC, _spacing, _r = l.split(',');
        _exonGC = float(_exonGC);
        _motifGC = float(_motifGC);
        _spacing = float(_spacing);
        _r = float(_r);
        lstData.append([_exonGC, _motifGC, _spacing, _r]);
    return lstData;


def PlotGCAllSpacings(title, all_correls_file, plot_filename_png):
# Plot 2x Box and Whiskers (Motif GC and Exon GC) of spread of CorrelStat correlation co-efficients for ALL spacings

	# Load correlations data
	lstCorrels = getAllCorrels(all_correls_file);

	fig = plt.figure(figsize=(8, 3.4))

	ax1 = fig.add_subplot(121)
	data000 = getCol(filter(lambda x: x[1] == 0, lstCorrels),3);
	data025 = getCol(filter(lambda x: x[1] == 25, lstCorrels),3);
	data050 = getCol(filter(lambda x: x[1] == 50, lstCorrels),3);
	data075 = getCol(filter(lambda x: x[1] == 75, lstCorrels),3);
	data100 = getCol(filter(lambda x: x[1] == 100, lstCorrels),3);
	data = [data000, data025, data050, data075, data100];
	ax1.boxplot(data, 0, '',  1);
	ax1.set_ylim([-0.3,1.02]);	
	ax1.set_ylabel('correlation');   
	ax1.set_xlabel('Motif GC content');
	ax1.set_xticklabels(["0%","25%","50%","75%","100%"]);

	ax2 = fig.add_subplot(122)
	data30_40 = getCol(filter(lambda x: x[0] == 35.0, lstCorrels),3);
	data40_50 = getCol(filter(lambda x: x[0] == 45.0, lstCorrels),3);
	data50_60 = getCol(filter(lambda x: x[0] == 55.0, lstCorrels),3);
	data60_70 = getCol(filter(lambda x: x[0] == 65.0, lstCorrels),3);
	data = [data30_40, data40_50, data50_60, data60_70];
	ax2.boxplot(data, 0, '', 1);
	ax2.set_ylim([-0.3,1.02]);	
	ax2.set_ylabel('correlation');   
	ax2.set_xlabel('Mean exon GC content');
	ax2.set_xticklabels(["35%","45%","55%","65%"]);

	fig.subplots_adjust(bottom=0.15, wspace = 0.3);

	#ax = fig.add_subplot(112);	

	plt.title(title, x=-0.25, y=1.03);
	plt.savefig(plot_filename_png);


PlotGCAllSpacings("Correlation vs GC content (All spacings)", "All-correls.csv", "All-correls.png");
