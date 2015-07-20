# make-tracks.py
# Takes parsed feature extraction data and makes UCSC microarray tracks for each array

import genutils
import re
import math

import sys
from optparse import  OptionParser

###############################################################################
USAGE = """
python make_tracks.py 		--input < README FILE > 
							--directory < Directory to write output files >

input == feature extraction output file
directory == Directory to write output files 

"""

parser = OptionParser(USAGE)
parser.add_option('--input',dest='input', help = 'input Agilent Feature Extraction output')
parser.add_option('--directory',dest='directory', help = 'Directory to write output files')
(options, args) = parser.parse_args()

if options.input is None:
    parser.error('input probe file name not given')
if options.directory is None:
    parser.error('output probe file name not given')

###############################################################################
###Input format:
#1		2		3			4					5					6			7			8			9		10			11					12					13						14							15
#Sample	Feature	ProbeType	ProbeID				Log2Ratio			Red			Green		Chromosome	Start	End			Log2_Red			Log2_Green			MeanCorrected_Log2Red	MeanCorrected_Log2Green		MeanCorrected_Log2Ratio
#VN37_5027	23	NonControl	chr14_794_2_DROP_R	0.00174936985345	516.0055	515.3802	chr14	58179690	58179749	9.01124263290472	9.00949330223319	0.904353586254016		0.903593034830137			1.00084169686414		


###Required output format:
#chrom	chromStart	chromEnd	name		score	strand	thickStart	thickEnd	reserved	blockCount	blockSizes	chromStarts	expCount	expIds	expScores
#chr1	159639972	159640031	2440848		500		-		159639972	159640031	0			1			59,			0,			33			0,1,2,	0.593000,1.196000,-0.190000,
#chr1	159640161	159640190	2440849		500		-		159640161	159640190	0			1			29,			0,			33			0,1,2,	-0.906000,-1.247000,0.111000,

###############################################################################
def write_multiWig_header(data):
	trackFile.write('track multiWig_CGHData\n')
	trackFile.write('type bigWig\n')
	trackFile.write('container multiWig\n')
	trackFile.write('aggregate transparent Overlay\n')
	trackFile.write('showSubtrackColoronUi on\n')
	trackFile.write('maxHeightPixels 500:100:8\n')
	trackFile.write('\n')
	
###############################################################################
def bed_to_bigWig_conversion(data):
	chromFile = '~/kidd-lab-scratch/www/track-hub/canFam3/canFam3.1-browser-chrom-sizes.fai'
	cmd = 'bedGraphToBigWig %s %s %s_MeanCorrLog2.bw' % (outFile, chromFile, data['sampleID'])
	print cmd
	genutils.runCMD(cmd)

##ultimately create data['bigWigFile']


###############################################################################
def write_track(data):
	trackFile.write('\ttrack CGHArray%s\n' % (data['sampleCount']))
	trackFile.write('\tbigDataUrl %s\n' % (data['bigWigFile'])
	trackFile.write('\tshortLabel Overlay_CGH%s\n' % data['sampleCount'])
	trackFile.write('\tlongLabel Z-scores of mean corrected Log2 ratios relative to controls for sample %s\n' % (data['sampleID'])
	trackFile.write('\tgraphTypeDefault bar\n')
	trackFile.write('\tparent multiWig_CGHData\n')
	trackFile.write('\ttype bigWig\n')
	trackFile.write('\tyLineOnOff on\n')
	trackFile.write('\tyLineMark 0.0\n')
	trackFile.write('\tyLineMark 1.0\n')
	trackFile.write('\tyLineMark -1.0\n')
	trackFile.write('\tgridDefault on\n')



###############################################################################

infile = options.input
inFile = open(options.input,'r')
print '\nReading in sample IDs from README file', options.input
print '\nReading in parsed extraction data from the following directory: ', options.directory

#Making multiwig aggregate track
trackfile = options.directory + 'TrackHub.txt'
trackFile = open(trackfile, 'w')
print '\nWriting trackhub information to', trackFile
write_multiWig_header(data)

lineNum = 0
sampleCount = 0
	
for line in inFile: #Reading through README file
	line = line.rstrip()
	line = line.split()
	lineNum+=1 
	lineCount = 0
	
	if lineNum > 1:	
		data = {}
		data['sampleID'] = line[0]		
		sampleCount += 1
		
		#Making Zscore track for UCSC
		outfile = options.directory + sampleID + '_MeanCorrectedZscores_ForUCSCTracks.bedGraph'
		outFile = open(outfile, 'w')
		
		#Read in mean corrected values from R script
		datafile = options.directory + sampleID + '_MeanCorrectedValues.txt' 
		dataFile = open(datafile, 'r')
		print 'Reading mean corrected values from', datafile
		
		for line in dataFile: #Reading through data file
			lineCount += 1
			data['sampleCount'] = sampleCount
			
			data['chrom'] = line[7]
			if data['chrom'] is 'NA':
				return
			else:
				data['start'] = line[8]
				data['end'] = line[9]
				data['name'] = line[3]
				#data['zscore'] = line[15]
				data['meanCorrLog2Ratio'] = line[14]
				outFile.write('%s\t%s\t%s\t%s' % (data['chrom'], data['start'], data['end'], data['zscore']))
		
			

		#outFile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (data['chrom'],data['start'],data['end'],data['name'],data['score'],data['strand'],data['thickStart'],data['thickEnd'],data['reserved'],data['blockCount'],data['blockSizes'],data['chromStarts'],data['expCount'],data['expIDs'],data['expScores'])