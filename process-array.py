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
    myData['probeName'] = line[3]
    myData['probeClass'] = line[2]
    myData['chrom'] = line[7]
    myData['startPos'] = line[8]
    myData['endPos'] = line[9]
    myData['red'] = float(line[5])
    myData['green'] = float(line[6])
    myData['isAuto'] = True
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
    myData['log2'] = np.log2( myData['red'] / myData['green'] )
    
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

	trackFile.write('\ttrack aCGH_%s_red\n' % sample)
	trackFile.write('\tparent aCGH_%s_multiwig\n' % sample)
	trackFile.write('\ttype bigWig\n')
	trackFile.write('\tbigDataUrl %s_red.bw\n' % sample)
	trackFile.write('\tshortLabel aCGH_%s_red\n' % sample)
	trackFile.write('\tlongLabel aCGH_%s_red\n' % sample)
	trackFile.write('\tgraphTypeDefault bar\n')
	trackFile.write('\tLineOnOff on\n')
	trackFile.write('\tyLineMark 0.0\n')
	trackFile.write('\tcolor 255,0,0\n\n')	
			
	trackFile.write('\ttrack aCGH_%s_green\n' % sample)
	trackFile.write('\tparent aCGH_%s_multiwig\n' % sample)
	trackFile.write('\ttype bigWig\n')
	trackFile.write('\tbigDataUrl %s_green.bw\n' % sample)
	trackFile.write('\tshortLabel aCGH_%s_green\n' % sample)
	trackFile.write('\tlongLabel aCGH_%s_green\n' % sample)
	trackFile.write('\tgraphTypeDefault bar\n')
	trackFile.write('\tLineOnOff on\n')
	trackFile.write('\tyLineMark 0.0\n')
	trackFile.write('\tcolor 0,255,0\n\n')	
	
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

	novelTrackFile.write('\ttrack aCGH_%s_red\n' % sample)
	novelTrackFile.write('\tparent aCGH_%s_multiwig\n' % sample)
	novelTrackFile.write('\ttype bigWig\n')
	novelTrackFile.write('\tbigDataUrl %s_red.bw\n' % sample)
	novelTrackFile.write('\tshortLabel aCGH_%s_red\n' % sample)
	novelTrackFile.write('\tlongLabel aCGH_%s_red\n' % sample)
	novelTrackFile.write('\tgraphTypeDefault bar\n')
	novelTrackFile.write('\tLineOnOff on\n')
	novelTrackFile.write('\tyLineMark 0.0\n')
	novelTrackFile.write('\tcolor 255,0,0\n\n')	
			
	novelTrackFile.write('\ttrack aCGH_%s_green\n' % sample)
	novelTrackFile.write('\tparent aCGH_%s_multiwig\n' % sample)
	novelTrackFile.write('\ttype bigWig\n')
	novelTrackFile.write('\tbigDataUrl %s_green.bw\n' % sample)
	novelTrackFile.write('\tshortLabel aCGH_%s_green\n' % sample)
	novelTrackFile.write('\tlongLabel aCGH_%s_green\n' % sample)
	novelTrackFile.write('\tgraphTypeDefault bar\n')
	novelTrackFile.write('\tLineOnOff on\n')
	novelTrackFile.write('\tyLineMark 0.0\n')
	novelTrackFile.write('\tcolor 0,255,0\n\n')
###############################################################################    
def make_stats_line(data):
	nl = []
	#From README
	nl.append(data['sample']) 
	nl.append(data['DateScanned'])
	nl.append(data['ScanNumber'])
	nl.append(data['SlideNumber'])
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

readme_file = options.input 
readmeFile = open(readme_file, 'r')

statsfile = options.input + '_STATS'
statsFile = open(statsfile, 'w')
#write header
statsFile.write('SampleID\tDateScanned\tScanNumber\tSlideNumber\tSex\tOrigin\tCanineType\tDLRS\tRawDataFile\t')
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
	
	if lineNumber > 1:
		sample = line[0]
		arrayCount += 1
		data['sample'] = line[0]
		data['DateScanned'] = line[1]
		data['ScanNumber'] = line[2]
		data['SlideNumber'] = line[3]
		data['DLRS'] = line[4]
		data['Sex'] = line[5]
		data['GeographicOrigin'] = line[6]
		data['CanineType'] = line[7]
		data['RawDataFile'] = line[8]

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
		
		autoControl = []
		chrXControl = []
		for probe in probeDataList:
			if probe['probeClass'] == 'CNVControl':
				if probe['isAuto'] is True:
					autoControl.append(probe['log2'])
				else:
					chrXControl.append(probe['log2'])
		
		# AUTOSOMES
		data['autoProbeCount'] = len(autoControl)
		print 'Autosomes %i probes' % (len(autoControl))            
		autoMean = np.mean(autoControl)
		autoMed = np.median(autoControl)
		autoMin = min(autoControl)
		autoMax = max(autoControl)
		print 'Mean %f median %f' % (autoMean,autoMed)
		print 'Min %f Max %f' % (autoMin,autoMax)
		data['autoMean'] = autoMean 
		data['autoMed'] = autoMed 
		data['autoMin'] = autoMin 
		data['autoMax'] = autoMax 

		# CHROMOSOME X
		data['chrXProbeCount'] = len(chrXControl)
		print 'chrX %i probes' % (len(chrXControl))            
		chrXMean = np.mean(chrXControl)
		chrXMed = np.median(chrXControl)
		chrXMin = min(chrXControl)
		chrXMax = max(chrXControl)		
		print 'Mean %f median %f' % (chrXMean,chrXMed)
		print 'Min %f Max %f' % (chrXMin,chrXMax)
		data['chrXMean'] = chrXMean
		data['chrXMed'] = chrXMean
		data['chrXMin'] = chrXMean
		data['chrXMax'] = chrXMean

		# CALCULATING CORRECTION FACTOR
		corfactor = 0.0 - autoMean
		print 'correction factor is %f' % corfactor
		print 'Applying to all data'
		data['correctionFactor'] = corfactor
		
		for probe in probeDataList:
			probe['log2'] += corfactor
		print 'Correction applied!\n\n'
		
		autoControl = []
		chrXControl = []
		for probe in probeDataList:
			if probe['probeClass'] == 'CNVControl':
				if probe['isAuto'] is True:
					autoControl.append(probe['log2'])
				else:
					chrXControl.append(probe['log2'])
		
		print 'Autosomes %i probes' % (len(autoControl))            
		autoMean = np.mean(autoControl)
		autoMed = np.median(autoControl)
		autoMin = min(autoControl)
		autoMax = max(autoControl)
		print 'Mean %f median %f' % (autoMean,autoMed)
		print 'Min %f Max %f' % (autoMin,autoMax)
		data['corr_autoMean'] = autoMean 
		data['corr_autoMed'] = autoMed 
		data['corr_autoMin'] = autoMin 
		data['corr_autoMax'] = autoMax 
		
		print 'chrX %i probes' % (len(chrXControl))            
		chrXMean = np.mean(chrXControl)
		chrXMed = np.median(chrXControl)
		chrXMin = min(chrXControl)
		chrXMax = max(chrXControl)
		print 'Mean %f median %f' % (chrXMean,chrXMed)
		print 'Min %f Max %f' % (chrXMin,chrXMax)
		data['corr_chrXMean'] = chrXMean
		data['corr_chrXMed'] = chrXMean
		data['corr_chrXMin'] = chrXMean
		data['corr_chrXMax'] = chrXMean
		
		print '\n___Stats___\n'
		
		autoSd = np.std(autoControl)
		print 'auto SD %f' % autoSd
		data['autoSD'] = autoSd
		
		# 1.5 cutoff
		
		
		# make simple bed graph
		
		# make bed graph file
		outFile = open('unsorted.bedGraph','w')
		
		for probe in probeDataList:
			if probe['chrom'] in ['NA','chrY']:
				continue
			if 'DROP' not in probe['probeName']:
				continue
				
			d = [probe['chrom'],probe['startPos'],int(probe['endPos']),probe['log2']]
			
			d = [str(i) for i in d]
			
			nl = '\t'.join(d) + '\n'
			outFile.write(nl)
		outFile.close()
		
		
		#sort cmd
		cmd = 'sort -k1,1 -k2,2n unsorted.bedGraph > sorted.bedGraph'
		print cmd
		genutils.runCMD(cmd)
		
		chromSizesFile = '/home/jmkidd/kidd-lab-scratch/www/track-hub/jmkidd-hubs/canFam3/canFam3.chromsizes'
		
		cmd = 'bedGraphToBigWig sorted.bedGraph %s %s%s.bw' % (chromSizesFile, trackDirectory, sample)
		print cmd
		genutils.runCMD(cmd)
		
		# now to setup for a multi track hub
		
		
		
		#0-1.5 std grey
		#1.5 to 2 black
		#>2 green/red
		
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
		outFile = open('tmp.green.bedGraph','w')
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
		cmd = 'bedGraphToBigWig tmp.green.bedGraph %s %s%s_green.bw' % (chromSizesFile, trackDirectory, sample)
		print cmd
		genutils.runCMD(cmd)

		# RED BIGWIG		
		outFile = open('tmp.red.bedGraph','w')
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
		cmd = 'bedGraphToBigWig tmp.red.bedGraph %s %s%s_red.bw' % (chromSizesFile, trackDirectory, sample)
		print cmd
		genutils.runCMD(cmd)
		
		nl = make_stats_line(data)
		statsFile.write(nl)		
		
		os.unlink('tmp.red.bedGraph')
		os.unlink('tmp.green.bedGraph')
		os.unlink('tmp.black.bedGraph')
		os.unlink('tmp.grey.bedGraph')
		os.unlink('unsorted.bedGraph')
		#os.unlink('sorted.bedGraph')


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
		
		# now to setup for a multi track hub
		
		
		
		#0-1.5 std grey
		#1.5 to 2 black
		#>2 green/red
		
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
		outFile = open('tmp.green.bedGraph','w')
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
		cmd = 'bedGraphToBigWig tmp.green.bedGraph %s %s%s_green.bw' % (novel_chromSizesFile, novelDirectory, sample)
		print cmd
		genutils.runCMD(cmd)

		# RED BIGWIG		
		outFile = open('tmp.red.bedGraph','w')
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
		cmd = 'bedGraphToBigWig tmp.red.bedGraph %s %s%s_red.bw' % (novel_chromSizesFile, novelDirectory, sample)
		print cmd
		genutils.runCMD(cmd)
		
		nl = make_stats_line(data)
		statsFile.write(nl)		
		
		os.unlink('tmp.red.bedGraph')
		os.unlink('tmp.green.bedGraph')
		os.unlink('tmp.black.bedGraph')
		os.unlink('tmp.grey.bedGraph')
		os.unlink('novel.unsorted.bedGraph')
		#os.unlink('sorted.bedGraph')
		
		#if arrayCount == 1:
			#break	


print '\n\nTotal arrays read...', arrayCount

cmd = 'cat %s*tracks > %sTotalTracks.txt' % (trackDirectory, trackDirectory)
print cmd
genutils.runCMD(cmd)

print 'Writing merged track files for subsequent upload to trackDb file to file: %sTotalTracks.txt' % (trackDirectory) 

cmd = 'cat %s*tracks > %sTotalTracks.txt' % (novelDirectory, novelDirectory)
print cmd
genutils.runCMD(cmd)

print 'DONE!!\n'




        
