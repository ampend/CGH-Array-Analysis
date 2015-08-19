# process-array.py
# 2015-07-21
import numpy as np
import genutils
import os
import re
from optparse import  OptionParser


###############################################################################
USAGE = """
python process-array.py 	--input < README FILE > 
							--directory < Directory to write output files >

input == README file
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
def process_array_line(line):
	myData = {}
	myData['sampleID'] = data['sample']
	myData['probeName'] = line[3]
	myData['probeClass'] = line[2]
	myData['chrom'] = line[7]
	myData['startPos'] = line[8]
	myData['endPos'] = line[9]
	myData['red'] = float(line[5])
	myData['green'] = float(line[6])
	
	#Defaults
	myData['isAuto'] = True
	myData['probeType'] = 'Other'
	myData['controlClass'] = 'Noncontrol'
	#Calc log2
	myData['log2'] = np.log2(myData['red']/myData['green'])
	myData['log2Red'] = np.log2(myData['red'])
	myData['log2Green'] = np.log2(myData['green'])
		
	if myData['chrom'] == 'chrY':
		myData['isAuto'] = False
	if myData['chrom'] == 'chrX':
		myData['isAuto'] = False		
	if myData['chrom'] == 'NA':
		myData['isAuto'] = False
	if re.match(r"wolf(\S+)", myData['chrom']) is not None:
		myData['isAuto'] = True
		myData['chrom'] == 'Novel'
	if re.match(r"zoey(\S+)", myData['chrom']) is not None:
		myData['isAuto'] = True
		myData['chrom'] == 'Novel'

	if myData['probeClass'] == 'AgilentControl':
		myData['controlClass'] = 'AgilentControl'
		myData['probeType'] = 'AgilentControl'
	if 'CNVControl' in myData['probeName']:
		if myData['chrom'] == 'chrX':
			myData['controlClass'] = 'SexChromControl'
		if myData['isAuto'] is True:
			myData['controlClass'] = 'AutosomeControl'
	if 'chrY' in myData['chrom']:
		myData['probeType'] = 'chrY'
	if 'DROP' in myData['probeName']:
		myData['probeType'] = 'DROP'
	if 'CNVControl' in myData['probeName']:
		myData['probeType'] = 'CNVControl'
	if myData['chrom'] == 'chrNovel':
		myData['probeType'] = 'Novel'
	if 'SINE' in myData['probeName']:
		myData['probeType'] = 'SINE'
	if 'LINE' in myData['probeName']:
		myData['probeType'] = 'LINE'
	if 'ERV' in myData['probeName']:
		myData['probeType'] = 'ERV'
	if 'umt' in myData['probeName']:
		myData['probeType'] = 'NUMT'
		
	return myData

###############################################################################    
def make_track_file(line):
	#canFam3 track files
	trackfile = options.directory + 'UCSCTracks/' + sample + '.tracks'
	trackFile = open(trackfile, 'w')
	
	trackFile.write('track aCGH_%s_multiwig\n' % sample)
	trackFile.write('type bigWig\n')
	trackFile.write('container multiWig\n')
	trackFile.write('shortLabel aCGH_%s_multiwig\n' % sample)
	trackFile.write('longLabel Mean corrected log2 ratios relative to controls only DROP probes\n')
	trackFile.write('aggregate transparentOverlay\n')
	trackFile.write('showSubtrackColoronUi on\n')
	trackFile.write('maxHeightPixels 500:100:8\n')
	trackFile.write('browser full all\n')
	trackFile.write('viewLimits -3:3\n\n')

	trackFile.write('\ttrack aCGH_%s_grey\n' % sample)
	trackFile.write('\tparent aCGH_%s_multiwig\n' % sample)
	trackFile.write('\ttype bigWig\n')
	trackFile.write('\tbigDataUrl %s_grey.bw\n' % sample)
	trackFile.write('\tshortLabel aCGH_%s_grey\n' % sample)
	trackFile.write('\tlongLabel aCGH_%s_grey\n' % sample)
	trackFile.write('\tgraphTypeDefault bar\n')
	trackFile.write('\tLineOnOff on\n')
	trackFile.write('\tyLineMark 0.0\n')
	trackFile.write('\tcolor 163,163,163\n\n')

	trackFile.write('\ttrack aCGH_%s_black\n' % sample)
	trackFile.write('\tparent aCGH_%s_multiwig\n' % sample)
	trackFile.write('\ttype bigWig\n')
	trackFile.write('\tbigDataUrl %s_black.bw\n' % sample)
	trackFile.write('\tshortLabel aCGH_%s_black\n' % sample)
	trackFile.write('\tlongLabel aCGH_%s_black\n' % sample)
	trackFile.write('\tgraphTypeDefault bar\n')
	trackFile.write('\tLineOnOff on\n')
	trackFile.write('\tyLineMark 0.0\n')
	trackFile.write('\tcolor 0,0,0\n\n')	

	trackFile.write('\ttrack aCGH_%s_orange\n' % sample)
	trackFile.write('\tparent aCGH_%s_multiwig\n' % sample)
	trackFile.write('\ttype bigWig\n')
	trackFile.write('\tbigDataUrl %s_orange.bw\n' % sample)
	trackFile.write('\tshortLabel aCGH_%s_orange\n' % sample)
	trackFile.write('\tlongLabel aCGH_%s_orange\n' % sample)
	trackFile.write('\tgraphTypeDefault bar\n')
	trackFile.write('\tLineOnOff on\n')
	trackFile.write('\tyLineMark 0.0\n')
	trackFile.write('\tcolor 255,128,0\n\n')	
			
	trackFile.write('\ttrack aCGH_%s_blue\n' % sample)
	trackFile.write('\tparent aCGH_%s_multiwig\n' % sample)
	trackFile.write('\ttype bigWig\n')
	trackFile.write('\tbigDataUrl %s_blue.bw\n' % sample)
	trackFile.write('\tshortLabel aCGH_%s_blue\n' % sample)
	trackFile.write('\tlongLabel aCGH_%s_blue\n' % sample)
	trackFile.write('\tgraphTypeDefault bar\n')
	trackFile.write('\tLineOnOff on\n')
	trackFile.write('\tyLineMark 0.0\n')
	trackFile.write('\tcolor 0,0,255\n\n')	
	
	#Novel sequence track files
	noveltrackFile = options.directory + 'NovelTracks/' + sample + '.tracks'
	novelTrackFile = open(noveltrackFile, 'w')
	
	novelTrackFile.write('track aCGH_%s_multiwig\n' % sample)
	novelTrackFile.write('type bigWig\n')
	novelTrackFile.write('container multiWig\n')
	novelTrackFile.write('shortLabel aCGH_%s_multiwig\n' % sample)
	novelTrackFile.write('longLabel Mean corrected log2 ratios relative to controls only DROP probes\n')
	novelTrackFile.write('aggregate transparentOverlay\n')
	novelTrackFile.write('showSubtrackColoronUi on\n')
	novelTrackFile.write('maxHeightPixels 500:100:8\n')
	novelTrackFile.write('viewLimits -3:3\n\n')

	novelTrackFile.write('\ttrack aCGH_%s_grey\n' % sample)
	novelTrackFile.write('\tparent aCGH_%s_multiwig\n' % sample)
	novelTrackFile.write('\ttype bigWig\n')
	novelTrackFile.write('\tbigDataUrl %s_grey.bw\n' % sample)
	novelTrackFile.write('\tshortLabel aCGH_%s_grey\n' % sample)
	novelTrackFile.write('\tlongLabel aCGH_%s_grey\n' % sample)
	novelTrackFile.write('\tgraphTypeDefault bar\n')
	novelTrackFile.write('\tLineOnOff on\n')
	novelTrackFile.write('\tyLineMark 0.0\n')
	novelTrackFile.write('\tcolor 163,163,163\n\n')

	novelTrackFile.write('\ttrack aCGH_%s_black\n' % sample)
	novelTrackFile.write('\tparent aCGH_%s_multiwig\n' % sample)
	novelTrackFile.write('\ttype bigWig\n')
	novelTrackFile.write('\tbigDataUrl %s_black.bw\n' % sample)
	novelTrackFile.write('\tshortLabel aCGH_%s_black\n' % sample)
	novelTrackFile.write('\tlongLabel aCGH_%s_black\n' % sample)
	novelTrackFile.write('\tgraphTypeDefault bar\n')
	novelTrackFile.write('\tLineOnOff on\n')
	novelTrackFile.write('\tyLineMark 0.0\n')
	novelTrackFile.write('\tcolor 0,0,0\n\n')	

	novelTrackFile.write('\ttrack aCGH_%s_orange\n' % sample)
	novelTrackFile.write('\tparent aCGH_%s_multiwig\n' % sample)
	novelTrackFile.write('\ttype bigWig\n')
	novelTrackFile.write('\tbigDataUrl %s_orange.bw\n' % sample)
	novelTrackFile.write('\tshortLabel aCGH_%s_orange\n' % sample)
	novelTrackFile.write('\tlongLabel aCGH_%s_orange\n' % sample)
	novelTrackFile.write('\tgraphTypeDefault bar\n')
	novelTrackFile.write('\tLineOnOff on\n')
	novelTrackFile.write('\tyLineMark 0.0\n')
	novelTrackFile.write('\tcolor 255,128,0\n\n')	
			
	novelTrackFile.write('\ttrack aCGH_%s_blue\n' % sample)
	novelTrackFile.write('\tparent aCGH_%s_multiwig\n' % sample)
	novelTrackFile.write('\ttype bigWig\n')
	novelTrackFile.write('\tbigDataUrl %s_blue.bw\n' % sample)
	novelTrackFile.write('\tshortLabel aCGH_%s_blue\n' % sample)
	novelTrackFile.write('\tlongLabel aCGH_%s_blue\n' % sample)
	novelTrackFile.write('\tgraphTypeDefault bar\n')
	novelTrackFile.write('\tLineOnOff on\n')
	novelTrackFile.write('\tyLineMark 0.0\n')
	novelTrackFile.write('\tcolor 0,0,192\n\n')
###############################################################################    
def make_stats_line(data):
	nl = []
	#From README
	nl.append(data['sample']) 
	nl.append(data['DateScanned'])
	nl.append(data['ScanNumber'])
	nl.append(data['SlideNumber'])
	nl.append(data['SampleNumber'])
	nl.append(data['DLRS'])
	nl.append(data['Sex'])
	nl.append(data['GeographicOrigin'])
	nl.append(data['CanineType'])
	nl.append(data['RawDataFile'])
	#Calculated in this script
	nl.append(data['autoMean']) # autoMean
	nl.append(data['autoMed']) # autoMed
	nl.append(data['autoMin']) # autoMin
	nl.append(data['autoMax']) # autoMax
	nl.append(data['chrXMean']) # chrXMean
	nl.append(data['chrXMed']) # chrXMed
	nl.append(data['chrXMin']) # chrXMin
	nl.append(data['chrXMax']) # chrXMax
	nl.append(data['correctionFactor']) # correction factor
	nl.append(data['corr_autoMean']) # Corrected autoMean
	nl.append(data['corr_autoMed']) # Corrected autoMed
	nl.append(data['corr_autoMin']) # Corrected autoMin
	nl.append(data['corr_autoMax']) # Corrected autoMax
	nl.append(data['corr_chrXMean']) # Corrected chrXMean
	nl.append(data['corr_chrXMed']) # Corrected chrXMed
	nl.append(data['corr_chrXMin']) # Corrected chrXMin
	nl.append(data['corr_chrXMax']) # Corrected chrXMax
	nl.append(data['autoSD']) # Autosome SD
	nl.append(data['firstSDRange']) # = lower-upper (x <1.5 SD)
	nl.append(data['lowerSDRange']) # = lower1-upper1 (1.5 SD < x < 2 SD)
	nl.append(data['upperSDRange']) # = lower2-upper2 (x >2 SD)
	
	nl = [str(j) for j in nl]
	nl = '\t'.join(nl) + '\n'
	return nl	
###############################################################################
def write_R_cmds(data):
	rFile.write('require(sm)\n')
	rFile.write('require(vioplot)\n')

	rFile.write('data <- read.table("%s_ProbeData.txt", header = TRUE)\n' % (data['sample']))
	rFile.write('colnames(data) <- c("Sample", "ProbeID", "ProbeClass", "ProbeType", "ProbeChrom", "ProbeStart", "ProbeEnd", "Corr_Log2", "Corr_Red", "Corr_Green", "Uncorr_Log2", "Uncorr_Red", "Uncorr_Green", "Sex")\n')
	rFile.write('sample <- as.factor(unlist(data[1]))\n')
	rFile.write('probeID <- as.factor(unlist(data[2]))\n')
	rFile.write('probeClass <- as.factor(unlist(data[3]))\n')
	rFile.write('probeType <- as.factor(unlist(data[4]))\n')
	rFile.write('probeChrom <- as.factor(unlist(data[5]))\n')
	rFile.write('probeStart <- sapply(data[6],as.numeric)\n')
	rFile.write('probeEnd <- sapply(data[7],as.numeric)\n')
	rFile.write('corrLog2 <- sapply(data[8],as.numeric)\n')
	rFile.write('corrRed <- sapply(data[9],as.numeric)\n')
	rFile.write('corrGreen <- sapply(data[10],as.numeric)\n')
	rFile.write('uncorrLog2 <- sapply(data[11],as.numeric)\n')
	rFile.write('uncorrRed <- sapply(data[12],as.numeric)\n')
	rFile.write('uncorrGreen <- sapply(data[13],as.numeric)\n')
	rFile.write('sex <- as.factor(unlist(data[14]))\n')

	rFile.write('#Defining Dataframe Variables\n')
	rFile.write('Dataframe <- data.frame(probeClass, corrLog2)\n')
	rFile.write('sex <- Dataframe$Corr_Log2[Dataframe$probeClass==\'SexChromControl\']\n')
	rFile.write('auto <- Dataframe$Corr_Log2[Dataframe$probeClass==\'AutosomeControl\']\n')
	rFile.write('non <- Dataframe$Corr_Log2[Dataframe$probeClass==\'Noncontrol\']\n')

	rFile.write('pdf(file="%s_aCGHPlots.pdf", width = 7,height = 10)\n' % (data['sample']))
	rFile.write('par(mfrow=c(1,1), mai=c(1,0.5,0.5,0.25), pin=c(5,5))\n')
	rFile.write('#Making Log2 Boxplots\n')
	rFile.write('plot(1,1,xlim=c(0,4),ylim=range(c(sex,auto,non)),type="n", xlab="Control Types",ylab="Corrected Log2 Ratios",axes=FALSE, cex.main=0.75, cex.lab=0.75, cex.axis=0.75)\n')
	rFile.write('axis(side=1,at=1:3,labels=c("Sex Chrom.","Autosome","Noncontrol"), cex.main=0.75, cex.lab=0.75, cex.axis=0.75)\n')
	rFile.write('axis(side=2)\n')
	rFile.write('boxplot(sex, at = 1, add=TRUE)\n')
	rFile.write('boxplot(auto, at = 2, add=TRUE)\n')
	rFile.write('boxplot(non, at = 3, add = TRUE, main = "Corrected Log2 Ratios for %s (%s)")\n' % (data['sample'], data['Sex']))
	rFile.write('#Making Log2 Violin Plots\n')
	rFile.write('plot(1,1,xlim=c(0,4),ylim=range(c(sex,auto,non)),type="n", xlab="Control Types",ylab="Corrected Log2 Ratios",axes=FALSE)\n')
	rFile.write('axis(side=1,at=1:3,labels=c("Sex Chrom.","Autosome","Noncontrol"), cex.main=0.75, cex.lab=0.75, cex.axis=0.75)\n')
	rFile.write('axis(side=2)\n')
	rFile.write('vioplot(sex, at = 1, add=TRUE, col="lightskyblue2")\n')
	rFile.write('vioplot(auto, at = 2, add=TRUE, col="lightskyblue3")\n')
	rFile.write('vioplot(non, at = 3, add=TRUE, col="lightskyblue4")\n')
	rFile.write('title("Corrected Log2 Ratios for %s (%s)")\n' % (data['sample'], data['Sex']))

	rFile.write('colorDataframe <- data.frame(probeClass, corrRed, corrGreen)\n')
	rFile.write('sexRed <- colorDataframe$Corr_Red[Dataframe$probeClass==\'SexChromControl\']\n')
	rFile.write('autoRed <- colorDataframe$Corr_Red[Dataframe$probeClass==\'AutosomeControl\']\n')
	rFile.write('nonRed <- colorDataframe$Corr_Red[Dataframe$probeClass==\'Noncontrol\']\n')
	rFile.write('sexGreen <- colorDataframe$Corr_Green[Dataframe$probeClass==\'SexChromControl\']\n')
	rFile.write('autoGreen <- colorDataframe$Corr_Green[Dataframe$probeClass==\'AutosomeControl\']\n')
	rFile.write('nonGreen <- colorDataframe$Corr_Green[Dataframe$probeClass==\'Noncontrol\']\n')

	rFile.write('###Corrected Red/Green Intensities\n')
	rFile.write('# Boxplots\n')
	rFile.write('plot(2,2,xlim=c(0,7),ylim=range(c(sexRed,autoRed,nonRed,sexGreen,autoGreen,nonGreen)),type="n", xlab="Control Types",ylab="Corrected Log2 Probe Intensities",axes=FALSE)\n')
	rFile.write('axis(side=1,at=1:6,labels=c("SexChr\\n(Sample)", "Auto\\n(Sample)", "Non\\n(Sample)", "SexChr\\n(Ref.)", "Auto\\n(Ref.)", "Non\\n(Ref.)"), cex.main=0.75, cex.lab=0.55, cex.axis=0.55)\n')
	rFile.write('axis(side=2)\n')
	rFile.write('boxplot(sexRed, at = 1, col=c("darkorange"), add=TRUE)\n')
	rFile.write('boxplot(autoRed, at = 2, col=c("darkorange1"), add=TRUE)\n')
	rFile.write('boxplot(nonRed, at = 3, col=c("darkorange3"), add = TRUE)\n')
	rFile.write('boxplot(sexGreen, at = 4,col=c("dodgerblue"), add=TRUE)\n')
	rFile.write('boxplot(autoGreen, at = 5, col=c("dodgerblue3"), add=TRUE)\n')
	rFile.write('boxplot(nonGreen, at = 6, col=c("dodgerblue4"), add = TRUE, main = "Corrected Log2 Probe Intensities for %s (%s)")\n' % (data['sample'], data['Sex']))

	rFile.write('# Violin Plots\n')
	rFile.write('plot(2,2,xlim=c(0,7),range(c(sexRed,autoRed,nonRed,sexGreen,autoGreen,nonGreen)),type="n", xlab="Control Types",ylab="Corrected Log2 Probe Intensities",axes=FALSE, main = "Corrected Log2 Probe Intensities for %s (%s)")\n' % (data['sample'], data['Sex']))
	rFile.write('axis(side=1,at=1:6,labels=c("SexChr\\n(Sample)", "Auto\\n(Sample)", "Non\\n(Sample)", "SexChr\\n(Ref.)", "Auto\\n(Ref.)", "Non\\n(Ref.)"), cex.main=0.75, cex.lab=0.55, cex.axis=0.55)\n')
	rFile.write('axis(side=2)\n')
	rFile.write('vioplot(sexRed, at = 1, add=TRUE, col="darkorange")\n')
	rFile.write('vioplot(autoRed, at = 2, add=TRUE, col="darkorange1")\n')
	rFile.write('vioplot(nonRed, at = 3, add=TRUE, col="darkorange3")\n')
	rFile.write('vioplot(sexGreen, at = 4, add=TRUE, col="dodgerblue")\n')
	rFile.write('vioplot(autoGreen, at = 5, add=TRUE, col="dodgerblue3")\n')
	rFile.write('vioplot(nonGreen, at = 6, add=TRUE, col="dodgerblue4")\n')

	rFile.write('title("Corrected Probe Intensities for %s (%s)")\n' % (data['sample'], data['Sex']))
	rFile.write('dev.off()\n')	
###############################################################################        

readme_file = options.input 
readmeFile = open(readme_file, 'r')

statsfile = options.directory + 'TotalStats/ArrayStatsPerSample.txt'
statsFile = open(statsfile, 'w')
#write header
statsFile.write('SampleID\tDateScanned\tScanNumber\tSlideNumber\tSampleNumber\tSex\tOrigin\tCanineType\tDLRS\tRawDataFile\t')
statsFile.write('AutosomeMean\tAutosomeMedian\tAutosomeMin\tAutosomeMax\t')
statsFile.write('ChrXMean\tChrXMedian\tChrXMin\tChrXMax\t')
statsFile.write('CorrectionFactor\t')
statsFile.write('CorrectedAutosomeMean\tCorrectedAutosomeMedian\tCorrectedAutosomeMin\tCorrectedAutosomeMax\t')
statsFile.write('CorrectedChrXMean\tCorrectedChrXMedian\tCorrectedChrXMin\tCorrectedChrXMax\t')
statsFile.write('AutosomalSD\tLessThan1.5SDRange\tBetween1.5-2SDRange\tGreaterThan2SDRange\n')

lineNumber = 0 

trackDirectory = options.directory + 'UCSCTracks/'
novelDirectory = options.directory + 'NovelTracks/'

arrayCount = 0 

for line in readmeFile:
	line = line.rstrip()
	line = line.split()
	lineNumber += 1
	
	data = {}
	probeTable = []
	
	if lineNumber > 1:
		sample = line[0]
		arrayCount += 1
		data['sample'] = line[0]
		data['DateScanned'] = line[1]
		data['ScanNumber'] = line[2]
		data['SlideNumber'] = line[3]
		data['SampleNumber'] = line[4]
		data['Sex'] = line[5]
		data['GeographicOrigin'] = line[6]
		data['CanineType'] = line[7]
		data['DLRS'] = line[8]
		data['RawDataFile'] = line[9]

		probefile = options.directory + data['sample'] + '_ProbeData.txt'
		probeFile = open(probefile, 'w')
		probeFile.write('SampleID\tProbeName\tControlClass\tProbeType\tChrom\tStartPos\tEndPos\tCorrectedLog2Ratio\tCorrectedRedIntensity\tCorrectedGreenIntensity\tUncorrectedLog2Ratio\tUncorrectedRedIntensity\tUncorrectedGreenIntensity\tSex\tRawRed\tRawGreen\n')

		totalfile = options.directory + 'TotalProbePerformance.txt'
		totalFile = open(totalfile, 'w')

		make_track_file(line)
		
		infile = options.directory + sample + 'ProbeIntensityData'
		inFile = open(infile, 'r')
		
		print '\n\n_____________________________________\n'  
		print '\n\nReading parsed probe intensity data from', infile
		
		probeDataList = []
		for line in inFile:
			line = line.rstrip()
			line = line.split()
			if line[0] == 'SampleID':
				continue
			myProbe = process_array_line(line)
			probeDataList.append(myProbe)
		inFile.close()
		
		print '\nProcessing probe data for sample\n',  sample
		print 'Read in for %i probes' % len(probeDataList)
		
		#log2
		autoControl = []
		chrXControl = []
		#probe intensities log2 red/green
		autoControlRed = []
		autoControlGreen = []
		chrXControlRed = []
		chrXControlGreen = []
		
		for probe in probeDataList:
			if probe['probeClass'] == 'CNVControl':
				if probe['isAuto'] is True:
					autoControl.append(probe['log2'])
					autoControlRed.append(probe['log2Red'])
					autoControlGreen.append(probe['log2Green'])
				else:
					chrXControl.append(probe['log2'])
					chrXControlRed.append(probe['log2Red'])
					chrXControlGreen.append(probe['log2Green'])
		
		### AUTOSOMES
		data['autoProbeCount'] = len(autoControl)
		print 'Autosomes %i probes' % (len(autoControl))            
		#Log2
		autoMean = np.mean(autoControl)
		autoMed = np.median(autoControl)
		autoMin = min(autoControl)
		autoMax = max(autoControl)
		print '#AutoLog2\nMean %f median %f' % (autoMean,autoMed)
		print 'Min %f Max %f' % (autoMin,autoMax)
		data['autoMean'] = autoMean 
		data['autoMed'] = autoMed 
		data['autoMin'] = autoMin 
		data['autoMax'] = autoMax 
		#Red
		autoMeanRed = np.mean(autoControlRed)
		autoMedRed = np.median(autoControlRed)
		autoMinRed = min(autoControlRed)
		autoMaxRed = max(autoControlRed)
		print '#Auto_Red\nMean Red = %f Median Red = %f' % (autoMeanRed,autoMedRed)
		print 'Min Red = %f Max Red = %f' % (autoMinRed,autoMaxRed)
		data['autoMeanRed'] = autoMeanRed
		data['autoMedRed'] = autoMedRed 
		data['autoMinRed'] = autoMinRed 
		data['autoMaxRed'] = autoMaxRed 
		#Green
		autoMeanGreen = np.mean(autoControlGreen)
		autoMedGreen = np.median(autoControlGreen)
		autoMinGreen = min(autoControlGreen)
		autoMaxGreen = max(autoControlGreen)
		print '#Auto_Green\nMean Green = %f Median Green = %f' % (autoMeanGreen,autoMedGreen)
		print 'Min = %f Max = %f' % (autoMin,autoMax)
		data['autoMeanGreen'] = autoMeanGreen 
		data['autoMedGreen'] = autoMedGreen 
		data['autoMinGreen'] = autoMinGreen 
		data['autoMaxGreen'] = autoMaxGreen 
		

		### CHROMOSOME X
		data['chrXProbeCount'] = len(chrXControl)
		print 'chrX %i probes' % (len(chrXControl))            
		#Log2
		chrXMean = np.mean(chrXControl)
		chrXMed = np.median(chrXControl)
		chrXMin = min(chrXControl)
		chrXMax = max(chrXControl)		
		print '#ChrXLog2\nMean %f median %f' % (chrXMean,chrXMed)
		print 'Min %f Max %f' % (chrXMin,chrXMax)
		data['chrXMean'] = chrXMean
		data['chrXMed'] = chrXMed
		data['chrXMin'] = chrXMin
		data['chrXMax'] = chrXMax
		#Red
		chrXMeanRed = np.mean(chrXControlRed)
		chrXMedRed = np.median(chrXControlRed)
		chrXMinRed = min(chrXControlRed)
		chrXMaxRed = max(chrXControlRed)		
		print '#ChrXRed\nMean Red = %f Median Red = %f' % (chrXMeanRed,chrXMedRed)
		print 'Min Red = %f Max Red = %f' % (chrXMinRed,chrXMaxRed)
		data['chrXMeanRed'] = chrXMeanRed
		data['chrXMedRed'] = chrXMedRed
		data['chrXMinRed'] = chrXMinRed
		data['chrXMaxRed'] = chrXMaxRed
		#Green
		chrXMeanGreen = np.mean(chrXControlGreen)
		chrXMedGreen = np.median(chrXControlGreen)
		chrXMinGreen = min(chrXControlGreen)
		chrXMaxGreen = max(chrXControlGreen)		
		print '#AutoGreen\nMean Green = %f Median Green = %f' % (chrXMeanGreen,chrXMedGreen)
		print 'Min Green = %f Max Green = %f' % (chrXMinGreen,chrXMaxGreen)
		data['chrXMeanGreen'] = chrXMeanGreen
		data['chrXMedGreen'] = chrXMedGreen
		data['chrXMinGreen'] = chrXMinGreen
		data['chrXMaxGreen'] = chrXMaxGreen


		# CALCULATING CORRECTION FACTOR
		corfactor = 0.0 - autoMean
		corfactorRed = 10.0 - autoMeanRed
		corfactorGreen = 10.0 - autoMeanGreen
		
		print '\n#CorrectionFactors\nLog2 correction factor is %f' % corfactor
		print 'Red correction factor is %f' % corfactorRed
		print 'Green correction factor is %f' % corfactorGreen

		print 'Applying to all data'
		data['correctionFactor'] = corfactor
		data['correctionFactorRed'] = corfactorRed
		data['correctionFactorGreen'] = corfactorGreen

		for probe in probeDataList:
			probe['uncorr_logratio'] = probe['log2']
			probe['uncorr_Red'] = probe['log2Red']
			probe['uncorr_Green'] = probe['log2Green']
			probe['log2'] += corfactor
			probe['log2Red'] += corfactorRed
			probe['log2Green'] += corfactorGreen
			
		print 'Correction applied!\n\n'
		
		autoControl = []
		chrXControl = []
		for probe in probeDataList:
			if probe['probeClass'] == 'CNVControl':
				if probe['isAuto'] is True:
					autoControl.append(probe['log2'])
					autoControlRed.append(probe['log2Red'])
					autoControlGreen.append(probe['log2Green'])
				else:
					chrXControl.append(probe['log2'])
					chrXControlRed.append(probe['log2Red'])
					chrXControlGreen.append(probe['log2Green'])
												
		#Includes corrected values
		print '\n#Corrected Autosome Data'
		print 'Autosomes %i probes' % (len(autoControl))            
		#Log2
		autoMean = np.mean(autoControl)
		autoMed = np.median(autoControl)
		autoMin = min(autoControl)
		autoMax = max(autoControl)
		print 'Mean Corrected Log2 %f Median Corrected Log2 %f' % (autoMean,autoMed)
		print 'Min Corrected Log2 %f Max Corrected Log2 %f' % (autoMin,autoMax)
		data['corr_autoMean'] = autoMean 
		data['corr_autoMed'] = autoMed 
		data['corr_autoMin'] = autoMin 
		data['corr_autoMax'] = autoMax 
		#Red
		autoMeanRed = np.mean(autoControlRed)
		autoMedRed = np.median(autoControlRed)
		autoMinRed = min(autoControlRed)
		autoMaxRed = max(autoControlRed)
		print 'Mean Corrected Red %f Median Corrected Red %f' % (autoMeanRed,autoMedRed)
		print 'Min Corrected Red %f Max Corrected Red %f' % (autoMinRed,autoMaxRed)
		data['corr_autoMeanRed'] = autoMeanRed 
		data['corr_autoMedRed'] = autoMedRed
		data['corr_autoMinRed'] = autoMinRed 
		data['corr_autoMaxRed'] = autoMaxRed 
		#Green
		autoMeanGreen = np.mean(autoControlGreen)
		autoMedGreen = np.median(autoControlGreen)
		autoMinGreen = min(autoControlGreen)
		autoMaxGreen = max(autoControlGreen)
		print 'Mean Corrected Green %f Median Corrected Green %f' % (autoMeanGreen,autoMedGreen)
		print 'Min Corrected Green %f Max Corrected Green %f' % (autoMinGreen,autoMaxGreen)
		data['corr_autoMeanGreen'] = autoMeanGreen 
		data['corr_autoMedGreen'] = autoMedGreen 
		data['corr_autoMinGreen'] = autoMinGreen 
		data['corr_autoMaxGreen'] = autoMaxGreen 

		print '\n#Corrected chrX Data'		
		print 'chrX %i probes' % (len(chrXControl))            
		#Log 2
		chrXMean = np.mean(chrXControl)
		chrXMed = np.median(chrXControl)
		chrXMin = min(chrXControl)
		chrXMax = max(chrXControl)
		print 'Mean Corrected Log2 %f Median Corrected Log2 %f' % (chrXMean,chrXMed)
		print 'Min Corrected Log2 %f Max Corrected Log2 %f' % (chrXMin,chrXMax)
		data['corr_chrXMean'] = chrXMean
		data['corr_chrXMed'] = chrXMed
		data['corr_chrXMin'] = chrXMin
		data['corr_chrXMax'] = chrXMax
		# Red
		chrXMeanRed = np.mean(chrXControlRed)
		chrXMedRed = np.median(chrXControlRed)
		chrXMinRed = min(chrXControlRed)
		chrXMaxRed = max(chrXControlRed)
		print 'Mean Corrected Red %f Median Corrected Red %f' % (chrXMeanRed,chrXMedRed)
		print 'Min Corrected Red %f Max Corrected Red %f' % (chrXMinRed,chrXMaxRed)
		data['corr_chrXMeanRed'] = chrXMeanRed
		data['corr_chrXMedRed'] = chrXMedRed
		data['corr_chrXMinRed'] = chrXMinRed
		data['corr_chrXMaxRed'] = chrXMaxRed
		#Green
		chrXMeanGreen = np.mean(chrXControlGreen)
		chrXMedGreen = np.median(chrXControlGreen)
		chrXMinGreen = min(chrXControlGreen)
		chrXMaxGreen = max(chrXControlGreen)
		print 'Mean Corrected Green %f Median Corrected Green %f' % (chrXMeanGreen,chrXMedGreen)
		print 'Min Corrected Green %f Max Corrected Green %f' % (chrXMinGreen,chrXMaxGreen)
		data['corr_chrXMeanGreen'] = chrXMeanGreen
		data['corr_chrXMedGreen'] = chrXMedGreen
		data['corr_chrXMinGreen'] = chrXMinGreen
		data['corr_chrXMaxGreen'] = chrXMaxGreen		
			
		#Log2
		autoSd = np.std(autoControl)
		print 'Log2 auto SD %f' % autoSd
		data['autoSD'] = autoSd
		#Red
		autoSdRed = np.std(autoControlRed)
		print 'Red auto SD %f' % autoSdRed
		data['autoSDRed'] = autoSdRed
		#Green
		autoSdGreen = np.std(autoControlGreen)
		print 'Green auto SD %f\n' % autoSdGreen
		data['autoSDGreen'] = autoSdGreen		
	
		# 1.5 cutoff
		# make simple bed graph
		# make bed graph file
		outFile = open('unsorted.bedGraph','w')
		
		for probe in probeDataList:			
			###Writing probe file with corrected log2, red, and green data. For R commands to process
			p = [probe['sampleID'], probe['probeName'], probe['controlClass'], probe['probeType'], probe['chrom'],probe['startPos'],probe['endPos'],probe['log2'], probe['log2Red'], probe['log2Green'], probe['uncorr_logratio'], probe['uncorr_Red'], probe['uncorr_Green'], data['Sex'], probe['red'], probe['green']]
			p = [str(i) for i in p]			
			hd = '\t'.join(p) + '\n'
			probeFile.write(hd)
			totalFile.write(hd)
			
			#For making bedGraph
			if probe['chrom'] in ['NA','chrY']:
				continue
			if 'DROP' not in probe['probeName']:
				continue
			d = [probe['chrom'],probe['startPos'],int(probe['endPos']),probe['log2']]
			d = [str(i) for i in d]			
			nl = '\t'.join(d) + '\n'		
			outFile.write(nl)					
		outFile.close()
		probeFile.close()
		
		#sort cmd
		cmd = 'sort -k1,1 -k2,2n unsorted.bedGraph > sorted.bedGraph'
		print cmd
		genutils.runCMD(cmd)
		
		chromSizesFile = '/home/jmkidd/kidd-lab-scratch/www/track-hub/jmkidd-hubs/canFam3/canFam3.chromsizes'
		
		cmd = 'bedGraphToBigWig sorted.bedGraph %s %s%s.bw' % (chromSizesFile, trackDirectory, sample)
		print cmd
		genutils.runCMD(cmd)
		
		##############
		# now to setup for a multi track hub for DROPs		 
		# 	0-1.5 std grey
		# 	1.5 to 2 black
		# 	>2 blue/orange
		##############
		
		# GREY BIGWIG		
		lower = autoMean - 1.5 * autoSd
		upper = autoMean + 1.5 * autoSd
		outFile = open('tmp.grey.bedGraph','w')
		inFile = open('sorted.bedGraph','r')
		for line in inFile:
			ol = line
			line = line.rstrip()
			line = line.split()
			l2 = float(line[3])
			if l2 >= lower and l2 <= upper:
				outFile.write(ol)    
		inFile.close()
		outFile.close()
		cmd = 'bedGraphToBigWig tmp.grey.bedGraph %s %s%s_grey.bw' % (chromSizesFile, trackDirectory, sample)
		print cmd
		genutils.runCMD(cmd)
		
		
		lower1 = autoMean - 1.5 * autoSd
		lower2 = autoMean - 2.0 * autoSd
		
		upper1 = autoMean + 1.5 * autoSd
		upper2 = autoMean + 2.0 * autoSd
		
		##SD Cut-offs
		data['firstSDRange'] = '<' + str(upper)
		data['lowerSDRange'] = '+/-' + str(upper1)
		data['upperSDRange'] = '+/-' + str(upper2)

		# BLACK BIGWIG		
		outFile = open('tmp.black.bedGraph','w')
		inFile = open('sorted.bedGraph','r')
		for line in inFile:
			ol = line
			line = line.rstrip()
			line = line.split()
			l2 = float(line[3])
			if l2 >= lower2 and l2 < lower1:
				outFile.write(ol)    
			if l2 > upper1 and l2 <= upper2:
				outFile.write(ol)    
		inFile.close()
		outFile.close()
		cmd = 'bedGraphToBigWig tmp.black.bedGraph %s %s%s_black.bw' % (chromSizesFile, trackDirectory, sample)
		print cmd
		genutils.runCMD(cmd)		
		
		# GREEN BIGWIG
		outFile = open('tmp.blue.bedGraph','w')
		inFile = open('sorted.bedGraph','r')
		for line in inFile:
			ol = line
			line = line.rstrip()
			line = line.split()
			l2 = float(line[3])
			if l2  < lower2:
				outFile.write(ol)    
		inFile.close()
		outFile.close()
		cmd = 'bedGraphToBigWig tmp.blue.bedGraph %s %s%s_blue.bw' % (chromSizesFile, trackDirectory, sample)
		print cmd
		genutils.runCMD(cmd)

		# RED BIGWIG		
		outFile = open('tmp.orange.bedGraph','w')
		inFile = open('sorted.bedGraph','r')
		for line in inFile:
			ol = line
			line = line.rstrip()
			line = line.split()
			l2 = float(line[3])
			if l2  > upper2:
				outFile.write(ol)    
		inFile.close()
		outFile.close()
		cmd = 'bedGraphToBigWig tmp.orange.bedGraph %s %s%s_orange.bw' % (chromSizesFile, trackDirectory, sample)
		print cmd
		genutils.runCMD(cmd)
		
		nl = make_stats_line(data)
		statsFile.write(nl)
				
		os.unlink('tmp.orange.bedGraph')
		os.unlink('tmp.blue.bedGraph')
		os.unlink('tmp.black.bedGraph')
		os.unlink('tmp.grey.bedGraph')
		os.unlink('unsorted.bedGraph')		
		
		
		###NOVEL CONTIGS
		# make bed graph file
		novelFile = open('novel.unsorted.bedGraph','w')
		
		for probe in probeDataList:
			if probe['chrom'] in ['NA']:
				continue
			if 'DROP' in probe['probeName']:
				continue
			if re.match(r"chr(\d+)", probe['chrom']) is not None:
				continue
			if probe['chrom'] == 'chrX':
				continue				
			else:				
				d = [probe['chrom'],probe['startPos'],int(probe['endPos']),probe['log2']]
			
				d = [str(i) for i in d]
			
				nl = '\t'.join(d) + '\n'
				novelFile.write(nl)
		novelFile.close()
		
		#sort cmd
		cmd = 'sort -k1,1 -k2,2n novel.unsorted.bedGraph > novel.sorted.bedGraph'
		print cmd
		genutils.runCMD(cmd)
		
		#novel_chromSizesFile = '/home/jmkidd/kidd-lab-scratch/feichens-projects/kmer/canFam31/unique_kmers/canFam3.1-withnovel/canFam3.1.withUn.withNovel.fa.fai'
		novel_chromSizesFile = '/home/ampend/kidd-lab/ampend-projects/CGH_Array_Design_Dog/CGH_Array_Analysis/inputData/NovelGenome/novelGenome.chrom.sizes'
		
		cmd = 'bedGraphToBigWig novel.sorted.bedGraph %s %s%s.bw' % (novel_chromSizesFile, novelDirectory, sample)
		print cmd
		genutils.runCMD(cmd)
		
		##############
		# now to setup for a multi track hub		
		# 	0-1.5 std grey
		# 	1.5 to 2 black
		# 	>2 blue/orange
		##############
		
		# GREY BIGWIG		
		lower = autoMean - 1.5 * autoSd
		upper = autoMean + 1.5 * autoSd
		outFile = open('tmp.grey.bedGraph','w')
		inFile = open('novel.sorted.bedGraph','r')
		for line in inFile:
			ol = line
			line = line.rstrip()
			line = line.split()
			l2 = float(line[3])
			if l2 >= lower and l2 <= upper:
				outFile.write(ol)    
		inFile.close()
		outFile.close()
		cmd = 'bedGraphToBigWig tmp.grey.bedGraph %s %s%s_grey.bw' % (novel_chromSizesFile, novelDirectory, sample)
		print cmd
		genutils.runCMD(cmd)
		
		
		lower1 = autoMean - 1.5 * autoSd
		lower2 = autoMean - 2.0 * autoSd
		
		upper1 = autoMean + 1.5 * autoSd
		upper2 = autoMean + 2.0 * autoSd
		
		##SD Cut-offs
		data['firstSDRange'] = '<' + str(upper)
		data['lowerSDRange'] = '+/-' + str(upper1)
		data['upperSDRange'] = '+/-' + str(upper2)
		
		###MAKING BIGWIG FILES
		# BLACK BIGWIG		
		outFile = open('tmp.black.bedGraph','w')
		inFile = open('novel.sorted.bedGraph','r')
		for line in inFile:
			ol = line
			line = line.rstrip()
			line = line.split()
			l2 = float(line[3])
			if l2 >= lower2 and l2 < lower1:
				outFile.write(ol)    
			if l2 > upper1 and l2 <= upper2:
				outFile.write(ol)    
		inFile.close()
		outFile.close()
		cmd = 'bedGraphToBigWig tmp.black.bedGraph %s %s%s_black.bw' % (novel_chromSizesFile, novelDirectory, sample)
		print cmd
		genutils.runCMD(cmd)		
		
		# GREEN BIGWIG
		outFile = open('tmp.blue.bedGraph','w')
		inFile = open('novel.sorted.bedGraph','r')
		for line in inFile:
			ol = line
			line = line.rstrip()
			line = line.split()
			l2 = float(line[3])
			if l2  < lower2:
				outFile.write(ol)    
		inFile.close()
		outFile.close()
		cmd = 'bedGraphToBigWig tmp.blue.bedGraph %s %s%s_blue.bw' % (novel_chromSizesFile, novelDirectory, sample)
		print cmd
		genutils.runCMD(cmd)

		# ORANGE BIGWIG		
		outFile = open('tmp.orange.bedGraph','w')
		inFile = open('novel.sorted.bedGraph','r')
		for line in inFile:
			ol = line
			line = line.rstrip()
			line = line.split()
			l2 = float(line[3])
			if l2  > upper2:
				outFile.write(ol)    
		inFile.close()
		outFile.close()
		cmd = 'bedGraphToBigWig tmp.orange.bedGraph %s %s%s_orange.bw' % (novel_chromSizesFile, novelDirectory, sample)
		print cmd
		genutils.runCMD(cmd)
		
		#nl = make_stats_line(data)
		#statsFile.write(nl)		
		
		os.unlink('tmp.orange.bedGraph')
		os.unlink('tmp.blue.bedGraph')
		os.unlink('tmp.black.bedGraph')
		os.unlink('tmp.grey.bedGraph')
		os.unlink('novel.unsorted.bedGraph')
		#os.unlink('sorted.bedGraph')
		
		rfile = options.directory + data['sample'] + '_Rcmds.r'
		print 'writing R commands to', rfile
		rFile = open(rfile, 'w')
		write_R_cmds(data)
		
		#if arrayCount == 3:
			#break	
		


print '\n\nTotal arrays read...', arrayCount

cmd = 'cat %s*tracks > %sTotalTracks.txt' % (trackDirectory, trackDirectory)
print cmd
genutils.runCMD(cmd)

print 'Writing merged track files for subsequent upload to trackDb file to file: %sTotalTracks.txt' % (trackDirectory) 

cmd = 'cat %s*tracks > %sTotalTracks.txt' % (novelDirectory, novelDirectory)
print cmd
genutils.runCMD(cmd)

print 'Writing merged probe files: %sTotal_ProbeData.txt' % (options.directory) 

cmd = 'cat %s*ProbeData.txt > %sTotal_ProbeIntensities.txt' % (options.directory, options.directory)
print cmd
genutils.runCMD(cmd)

print 'DONE!!\n'