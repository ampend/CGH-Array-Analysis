# extractFeatures.py
# extracts feature data from the output of the Agilent Feature Extraction software version 12.0

import genutils
import re
#import rpy2
import math

import sys
from optparse import  OptionParser

###############################################################################
USAGE = """
python extractFeatures.py 	--input < feature extraction output file> 
							--probes <probe BED file> 
							--sample <sample ID> 
							--directory < Directory to write output files >

input == feature extraction output file
probes == probe bed file for probes that can be aligned to reference
sample == sample ID
directory == Directory to write output files 

"""

parser = OptionParser(USAGE)
parser.add_option('--input',dest='input', help = 'input Agilent Feature Extraction output')
parser.add_option('--probes',dest='probes', help = 'probe bed files')
parser.add_option('--sample',dest='sample', help = 'sample ID')
parser.add_option('--directory',dest='directory', help = 'Directory to write output files')
(options, args) = parser.parse_args()

if options.input is None:
    parser.error('input probe file name not given')
if options.probes is None:
    parser.error('input probe file name not given')
if options.sample is None:
    parser.error('output probe file name not given')
if options.directory is None:
    parser.error('output probe file name not given')


###############################################################################
def defineControls(data): # Reads Feature Extraction output and probe IDs to determine whether or not a probe is a custom or Agilent (built-in) control
	data['isControl'] = False # Default: Includes both Agilent and Custom CNV Controls
	
	#Agilent Controls
	data['isAgilentControl'] = False
	data['isCorner'] = False
	if data['ControlType'] == 1:
		if re.match(r"\S+(Corner)", data['ProbeName']) is None:
			data['isControl'] = True
			data['isAgilentControl'] = True
		if re.match(r"\S+(Corner)", data['ProbeName']) is not None:
			data['isCorner'] = True
	#Custom Controls 
	data['isCNVControl'] = False
	if re.match(r"CNVControl(\S+)", data['ProbeName']) is not None: #Probe names starting with CNVControls are custom controls (N=975 on array, replicated x2)
		match = re.match(r"(\S+)\_[RF]", data['ProbeName'])
		data['isCNVControl'] = True
		data['isControl'] = True

###############################################################################
def writeControlData(data):
	if data['isControl'] is True:
		#writes green and red signals, and log2ratio for controls only (agilent + custom controls)
		controlFile.write('%s\t%s\t%s\t%s\t%s\n' % (data['sampleID'],data['ProbeName'], data['FLOATrProcessedSignal'], data['FLOATgProcessedSignal'], data['log2LogRatio']))

###############################################################################
def writeNonControlData(data):
	if data['isControl'] is False:
		if data['isCorner'] is False:
			#writes green and red signals, and log2ratio for controls only (agilent + custom controls)
			nonControlFile.write('%s\t%s\t%s\t%s\t%s\n' % (data['sampleID'],data['ProbeName'], data['FLOATrProcessedSignal'], data['FLOATgProcessedSignal'], data['log2LogRatio']))
	 
###############################################################################
def convertLog(data): # This converts log ratios from log base 10 to log base 2
	tmpLog = float(data['LogRatio'])
	x = float(math.pow(10,tmpLog))
	data['log2LogRatio'] = math.log(x,2)

###############################################################################
def getCoordinates(data): # This extracts coordinates of probes from a bed file that contains ONLY probes mappable to canFam3, if probe is missing, then default coordinates = NA remain for the probe
	probeBedFile = open(options.probes,'r')

	for line in probeBedFile:
		line = line.rstrip()
		line = line.split()	
    
		ID = line[3]

		if ID == data['name']:
			#print 'MATCH\t%s\t%s' % (data['name'],data['ProbeName'])
			data['hasCoordinates'] = True
			data['probeChrom'] = line[0]
			data['probeStart'] = line[1]
			data['probeEnd'] = line[2]
			#shortened, non-specific probe ID
			data['probeID'] = line[3]
			#make_track_file(data)

###############################################################################
def detProbeType(data):
	if data['isControl'] is False:
		if data['isCorner'] is False:
			data['probeType'] = 'NonControl'
			data['NonControlCount'] += 1
	if data['isAgilentControl'] is True:
		data['probeType'] = 'AgilentControl'
		data['AgilentControlCount'] += 1
	if data['isCNVControl'] is True:
		data['probeType'] = 'CNVControl'
		data['CNVControlCount'] += 1
	if data['isCorner'] is True:
		data['probeType'] = 'Corner'
		data['CornerCount'] += 1
			
###############################################################################
def make_probe_stats_print_line(data):
	nl = []
	nl.append(data['sampleID'])
	nl.append(data['featureNum'])
	nl.append(data['probeType']) #Control or NonControl
	nl.append(data['ProbeName'])
	nl.append(data['log2LogRatio'])
	nl.append(data['FLOATrProcessedSignal'])
	nl.append(data['FLOATgProcessedSignal'])
	nl.append(data['probeChrom'])
	nl.append(data['probeStart'])
	nl.append(data['probeEnd'])
	
	nl = [str(j) for j in nl]
	nl = '\t'.join(nl) + '\n'
	return nl
		
###############################################################################

numProbes = 0
numPassStats = 0

inFile = open(options.input,'r')
print '\nReading in probe intensity data for Sample: ', options.sample

# Probe bed file can be found in the following directory: /home/ampend/kidd-lab/jmkidd-projects/people-projects/ampend-projects/CGH-probe-selection/results/FinalProbeSet/BEDFiles_Probes_ForUCSCBrowser/IntersectionOfProbesets
print 'Reading in probe coordinates from probe BED file', options.probes

#current directory : /home/ampend/kidd-lab/ampend-projects/CGH_Array_Design_Dog/CGH_Array_Analysis/results
print 'Writing outfiles to the following directory: ', options.directory

outfile = options.directory + options.sample + 'ProbeIntensityData' #Generating output file based on the sample ID provided 
outFile = open(outfile,'w')
print '\n\nWriting probe intensity data for controls AND noncontrols only to', outfile

#File for controls (agilent + custom) fluorescence data only
controlfile = outfile + '.controls'
controlFile = open(controlfile,'w')
print 'Writing probe intensity data for controls only to', controlfile

#File for controls (agilent + custom) fluorescence data only
noncontrolfile = outfile + '.noncontrols'
nonControlFile = open(noncontrolfile,'w')
print 'Writing probe intensity data for noncontrols only to', noncontrolfile
print '\n'


##WRITING HEADERS FOR OUTFILES
outFile.write('SampleID\tFeatureNumber\tProbeID\tProbeType\tLog2Ratio(Red/Green)\tProcessedRedSignal\tProcessedGreenSignal\tChromosome\tStart\tEnd\n') #header line
controlFile.write('SampleID\tProbeID\tRedProcessedSignal\tGreenProcessedSignal\tLog2Ratio(log2(RProcessed/GProcessed))\n')
nonControlFile.write('SampleID\tProbeID\tRedProcessedSignal\tGreenProcessedSignal\tLog2Ratio(log2(RProcessed/GProcessed))\n')

lineNum = 0
data = {}
data['NonControlCount'] = 0
data['AgilentControlCount'] = 0
data['CNVControlCount'] = 0
data['CornerCount'] = 0

for line in inFile:
	line = line.rstrip()
	line = line.split()
	lineNum +=1 
	start = line[0]

	
	#Files 
	data['inFile'] = inFile
	data['outFile'] = outFile
	data['controlFile'] = controlFile
	data['nonControlFile'] = nonControlFile
	data['sampleID'] = options.sample
	

	
#Reading in the columns from Feature Extraction output
	if start == 'DATA':
		#First seven lines are not raw data from Feature Extraction but descriptions of
		# 	slide and run info, or headers. Data begins after line 8 consistently.	
		if lineNum > 7:
			numProbes += 1
		
			#Reading probe data from Feature Extraction output file
			data['featureNum'] = int(line[1])
			data['row'] = int(line[2])
			data['col'] = int(line[3])
			data['SubTypeMask'] = int(line[4])
			data['ProbeName'] = line[6]			
			
			#### Testing for whether or not probe is an Agilent or one of our controls based
			#### 	on mrsFAST/quicKmer consensus
			
			data['ControlType'] = int(line[5]) #Binary
			defineControls(data)
			
			
			#Parsing Probe Names - some have F and R at the end of the probe IDs from 
			# 	feature extraction - This means the probe orientation was incorporated
			#Some probes only have one F or one R, some have two, but the BED coordinates
			#	would be the same for each, only orientation would matter later
			
			#Determining orientation based on Agilent Probe ID, default = Forward strand (F)
			data['probeOrientation'] = 'F'
			if re.match(r"(\S+)\_[R]", data['ProbeName']) is not None:
				 data['probeOrientation'] = 'R'
			#Parsing out probe names without 'F' and 'R' in the names to pull out 
			#	BED coordinates
			if re.match(r"(\S+)\_[RF]", data['ProbeName']) is not None:
				match = re.match(r"(\S+)\_[RF]", data['ProbeName'])
				data['name'] = match.group(1)
			else:
				data['name'] = data['ProbeName']
			
			#####################
			# Feature Extraction Data:
			#####################
			data['SystematicName'] = line[7]
			data['PositionX'] = line[8]
			data['PositionY'] = line[9]
	
			data['LogRatio'] = line[10] # Log base 10
			##Converting logratio with log base 10 to log base 2
			convertLog(data)
			data['LogRatioError'] = line[11] # of Log base 10 data
			data['PValueLogRatio'] = line[12] # of Log base 10 data
			
			######
			# Fluoresence Data
			#Green = Reference == Zoey
			#Red = Test 
			######
			#LogRatio = log10(data['FLOATgProcessedSignal']/data['FLOATgProcessedSignal'])
			data['gProcessedSignal'] = line[13]
			data['FLOATgProcessedSignal'] = float(data['gProcessedSignal'])
			data['rProcessedSignal'] = line[14]
			data['FLOATrProcessedSignal'] = float(data['rProcessedSignal'])
			
			writeControlData(data)
			writeNonControlData(data)
			
			loggreen = math.log(data['FLOATgProcessedSignal'],2)
			logred = math.log(data['FLOATrProcessedSignal'],2)
			
			ratio = data['FLOATrProcessedSignal']/data['FLOATgProcessedSignal']
			lograt = math.log(ratio,2)
			#print '%s\t%s' % (lograt,data['log2LogRatio'])
			data['gProcessedSigError'] = line[15]
			data['rProcessedSigError'] = line[16]
			data['gMedianSignal'] = float(line[17])
			data['rMedianSignal'] = float(line[18])
			data['gBGPixSDev'] = line[19]
			data['rBGPixSDev'] = line[20]
			data['gBGMedianSignal'] = line[21]
			data['rBGMedianSignal'] = line[22]
			data['gIsSaturated'] = float(line[23]) #Binary data - 0/1
			data['rIsSaturated'] = float(line[24]) #Binary data - 0/1
			data['gIsFeatNonUnifOL'] = float(line[25]) #Binary data - 0/1
			data['rIsFeatNonUnifOL'] = float(line[26]) #Binary data - 0/1
			data['gIsBGNonUnifOL'] = float(line[27]) #Binary data - 0/1
			data['rIsBGNonUnifOL'] = float(line[28]) #Binary data - 0/1
			data['gIsFeatPopnOL'] = float(line[29]) #Binary data - 0/1
			data['rIsFeatPopnOL'] = float(line[30]) #Binary data - 0/1
			data['IsManualFlag'] = float(line[31]) #Binary data - 0/1
			data['gIsPosAndSignif'] = line[36] #Binary data - 0/1
			data['rIsPosAndSignif'] = line[37] #Binary data - 0/1
			data['gIsWellAboveBG'] = line[38] #Binary data - 0/1
			data['rIsWellAboveBG'] = line[39] #Binary data - 0/1
			data['gBGSubSignal'] = line[34]
			data['rBGSubSignal'] = line[35]
			data['SpotExtentX'] = line[40]
			data['gBGMeanSignal'] = line[41]
			data['rBGMeanSignal'] = line[42]
			
			#####################
			# Default coordinates:
			# If mappable to UCSC, then the coordinates will be changed in getCoordinates
			# Else, coordinates are NA
			# This includes data from novel contigs, ChrY, and Reference Mobile Element Deletions
			#####################
			data['hasCoordinates'] = False
			data['probeChrom'] = 'NA'
			data['probeStart'] = 'NA'
			data['probeEnd'] = 'NA'
			#shortened, non-specific probe ID
			data['probeID'] = 'NA'

			#obtains probe coordinates based on data from Probe BED file
			getCoordinates(data)
			
			####################
			# Make file for R program - CGHNormaliter
			####################
			data['shortChrom'] = 'NA'
			if re.match(r"chr(\S+)", data['probeChrom']) is not None:
				match = re.match(r"chr(\S+)", data['probeChrom'])
				data['shortChrom'] = match.group(1)
			
			
			detProbeType(data) #Determines probe type : Control, NonControl, or Agilent
			
			if data['probeType'] is not 'Corner':
				nl = make_probe_stats_print_line(data)
				outFile.write(nl)
			#if numProbes > 1000:
				#break
				
				
sum = data['NonControlCount'] + data['CornerCount'] + data['AgilentControlCount'] + data['CNVControlCount']

print '%s probes were analyzed' % (numProbes)
print '%s probes are Agilent controls, %s probes are our CNV controls' % (data['AgilentControlCount'],data['CNVControlCount'])
print '%s probes are not controls' % (data['NonControlCount'])
print '%s probes are assigned to corners of arrays' % (data['CornerCount'])

	


#echo 'cd  $PBS_O_WORKDIR; python extractFeatures.py --input ../inputData/257548210013_2015-07-09_14-37_Area1_532_635_CGH_1200_Jun14_1_1.txt --probes ../inputData/ProbeBEDFiles/Mappable/TotalProbes_FINAL_Coordinates.bed --tracks ../results/out.tracks --sample test' |  qsub -l mem=40gb,walltime=10:00:00,nodes=3 -N Parse -j oe -V


