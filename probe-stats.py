# extractFeatures.py
# extracts feature data from the output of the Agilent Feature Extraction software version 12.0

import genutils
import glob
import numpy as np
import os

import sys
from optparse import  OptionParser

###############################################################################
USAGE = """
python probe_stats.py 		
							--directory < Directory to write output files >

directory == Directory to write output files 

"""

parser = OptionParser(USAGE)
parser.add_option('--directory',dest='directory', help = 'Directory to write output files')
(options, args) = parser.parse_args()

if options.directory is None:
    parser.error('output probe file name not given')

###############################################################################
def write_line(outFile):
	nl = []
	nl.append(probe)
	nl.append(controlClass)
	nl.append(probe_Type)
	nl.append(chrom)
	nl.append(probeStart)
	nl.append(probeEnd)
	nl.append(mean_log2)
	nl.append(mean_red)
	nl.append(mean_green)
	nl.append(mean_female_log2)
	nl.append(mean_female_red)
	nl.append(mean_female_green)
	nl.append(mean_male_log2)
	nl.append(mean_male_red)
	nl.append(mean_male_green)
	nl = [str(j) for j in nl]
	nl = '\t'.join(nl) + '\n'
	return nl

###############################################################################

fHD_array = []
probe_array = []

outfile = options.directory + 'TotalStats/ProbeAverages.txt'
outFile = open(outfile, 'w')
#write header
outFile.write('ProbeID\tControlType\tProbeType\tChromosome\tProbeStart\tProbeEnd\tMeanLog2\tMeanRed\tMeanGreen\tFemaleMeanLog2\tFemaleMeanRed\tFemaleMeanGreen\tMaleMeanLog2\tMaleMeanRed\tMaleMeanGreen\n')

for filename in glob.glob(options.directory+"*_ProbeData.txt"):
	fHD_array.append(open(filename, 'r'))

for line_firstfile in fHD_array[0]:
	#Do something
	arr_0 = line_firstfile.split()
	if arr_0[0] != 'SampleID':
		log2 = [float(arr_0[7])]
		red_I = [float(arr_0[8])]
		green_I = [float(arr_0[9])]
		
		female_log2 = []
		female_red = []
		female_green = []

		male_log2 = []
		male_red = []
		male_green=[]
	
		#track probe
		probe = arr_0[1]
		controlClass = arr_0[2]
		probe_Type = arr_0[3]
		chrom = arr_0[4]
		probeStart = arr_0[5]
		probeEnd = arr_0[6]
	
	for i in range(1, len(fHD_array)):
		line = fHD_array[i].readline().rstrip()
		arr = line.split()
		if arr[0] == 'SampleID':
			continue
		log2.append(float(arr[7]))
		red_I.append(float(arr[8]))
		green_I.append(float(arr[9]))
		if arr[13] == 'F':
			female_log2.append(float(arr[7]))
			female_red.append(float(arr[8]))
			female_green.append(float(arr[9]))
		if arr[13] == 'M':
			male_log2.append(float(arr[7]))
			male_red.append(float(arr[8]))
			male_green.append(float(arr[9]))
	
	if arr_0[0] == 'SampleID':
		continue
	#summary statistics for this probe
	mean_log2 = np.mean(log2)
	mean_red = np.mean(red_I)
	mean_green = np.mean(green_I)

	mean_female_log2 = np.mean(female_log2)
	mean_female_red = np.mean(female_red)
	mean_female_green = np.mean(female_green)
		
	mean_male_log2 = np.mean(male_log2)
	mean_male_red = np.mean(male_red)
	mean_male_green = np.mean(male_green)

	 
	#outFile.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (probe, probe_Type, chrom, mean_log2, mean_red, mean_green))
	nl = write_line(outFile)
	outFile.write(nl)


