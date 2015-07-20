# write-R-cmds.py
# extracts feature data from the output of the Agilent Feature Extraction software version 12.0

import genutils
import re

import sys
from optparse import  OptionParser

###############################################################################
USAGE = """
python write-R-cmds.py 		--input < Sample list of array files > 
							--directory < Directory where result files are located for processed array data >

input == Sample list of array files
directory == Directory where result files are located for processed array data

"""

parser = OptionParser(USAGE)
parser.add_option('--input',dest='input', help = 'Sample list of array files')
parser.add_option('--directory',dest='directory', help = 'Directory where result files are located for processed array data')
(options, args) = parser.parse_args()

if options.input is None:
    parser.error('input probe file name not given')
if options.directory is None:
    parser.error('output probe file name not given')

###############################################################################


inFile = open(options.input,'r')
print '\nReading in sample IDs and processed array files from: ', options.input

lineNum = 0
samples = 0

for line in inFile:
	line = line.rstrip()
	line = line.split()
	lineNum+=1 
	
	if lineNum > 1: #Skips header line
		sampleID = line[0]
		samples += 1
		outfile = options.directory + sampleID + 'Rcmds.r'
		outFile = open(outfile,'w')
		print 'Writing R commands to', outfile

		outFile.write('#############################################################################################\n')
		outFile.write('# TOTAL PROBE DATA\n')
		outFile.write('%s <- read.table("%s%sProbeIntensityData",header=TRUE)\n' % (sampleID, options.directory, sampleID))
		outFile.write('colnames(%s) <- c("Sample", "Feature", "ProbeType", "ProbeID", "Log2Ratio", "Red", "Green", "Chromosome", "Start", "End")\n'  % (sampleID))
		outFile.write('%stype<-as.factor(unlist(%s[3]))\n' % (sampleID, sampleID))
		outFile.write('%stype<-as.factor(unlist(%s[3]))\n' % (sampleID, sampleID))
		outFile.write('%slog2 <-sapply(%s[5],as.numeric)\n' % (sampleID, sampleID))
		outFile.write('%sDataframe <- data.frame(%stype, %slog2)\n' % (sampleID, sampleID, sampleID))
		outFile.write('pdf(file = "%s_ControlBoxPlot.pdf")\n' % (sampleID))
		outFile.write('boxplot(%slog2 ~ %stype, %sDataframe, main = "Boxplot of Log2 Ratios (Red:Green) for %s", ylab = "log2 Values")\n' % (sampleID, sampleID, sampleID, sampleID))
		outFile.write('dev.off()\n')
		outFile.write('%sred <-sapply(%s[6],as.numeric)\n' % (sampleID, sampleID))
		outFile.write('%sgreen <-sapply(%s[7],as.numeric)\n' % (sampleID, sampleID))
		outFile.write('log2%sred <- log2(%sred)\n' % (sampleID, sampleID))
		outFile.write('log2%sgreen <- log2(%sgreen)\n' % (sampleID, sampleID))
		outFile.write('log2ratio_total <- sapply(%s[5], as.numeric)\n' % (sampleID))
		
		outFile.write('\n# CONTROLS\n')
		outFile.write('controls <- read.table("%sProbeIntensityData.controls",header=TRUE)\n' % (sampleID))
		outFile.write('cred<-sapply(controls[3],as.numeric)\n')
		outFile.write('cgreen<-sapply(controls[4],as.numeric)\n')
		outFile.write('log2cred<-sapply((log2(controls[3])),as.numeric)\n')
		outFile.write('log2cgreen<-sapply((log2(controls[4])),as.numeric)\n')
		outFile.write('med_red%scontrols <- median(log2cred)\n' % (sampleID))
		outFile.write('mean_red%scontrols <- mean(log2cred)\n' % (sampleID))
		outFile.write('sd_red%scontrols <- sd(log2cred,na.rm=TRUE)\n' % (sampleID))
		outFile.write('med_green%scontrols <- median(log2cgreen)\n' % (sampleID))
		outFile.write('mean_green%scontrols <- mean(log2cgreen)\n' % (sampleID))
		outFile.write('sd_green%scontrols <- sd(log2cgreen,na.rm=TRUE)\n' % (sampleID))
		outFile.write('mean_corr_log2cred <- log2cred/mean_red%scontrols\n' % (sampleID))
		outFile.write('mean_corr_log2cgreen <- log2cred/mean_green%scontrols\n' % (sampleID))
		outFile.write('####Log2 ratios of controls only\n')
		outFile.write('log2controls <- sapply(controls[5],as.numeric)\n')
		outFile.write('mean_log2controls <- mean(log2controls)\n')
		outFile.write('sd_log2controls <- sd(log2controls, na.rm=TRUE)\n')
		outFile.write('mean_log2ratio_controls <- mean(sapply(controls[5], as.numeric))\n')
		outFile.write('sd_log2ratio_controls <- sd(sapply(controls[5], as.numeric))\n')
		outFile.write('mean_corr_log2ratio <- log2ratio_total/mean_log2ratio_controls\n')

		outFile.write('\n# TOTAL PROBE DATA\n')
		outFile.write('#Correcting Total Data based on mean of controls\n')
		outFile.write('mean_corr_log2%sred <- log2%sred/mean_red%scontrols\n' % (sampleID, sampleID, sampleID))
		outFile.write('mean_corr_log2%sgreen <- log2%sgreen/mean_green%scontrols\n' % (sampleID,sampleID,sampleID))
		outFile.write('#Calc Z Scores\n')
		outFile.write('zscore <- (mean_corr_log2ratio - mean_log2ratio_controls) / sd_log2ratio_controls\n')
		
		outFile.write('\n# NONCONTROLS\n')
		outFile.write('noncontrols <- read.table("%sProbeIntensityData.noncontrols",header=TRUE)\n' % (sampleID))
		outFile.write('nred<-sapply(noncontrols[3],as.numeric)\n')
		outFile.write('ngreen<-sapply(noncontrols[4],as.numeric)\n')
		outFile.write('log2nred<-sapply((log2(noncontrols[3])),as.numeric)\n')
		outFile.write('log2ngreen<-sapply((log2(noncontrols[4])),as.numeric)\n')
		outFile.write('med_corr_log2nred<-log2nred/med_red%scontrols\n' % (sampleID))
		outFile.write('mean_corr_log2nred<-log2nred/mean_red%scontrols\n' % (sampleID))
		outFile.write('med_corr_log2ngreen<-log2ngreen/med_green%scontrols\n' % (sampleID))
		outFile.write('mean_corr_log2ngreen<-log2ngreen/mean_green%scontrols\n'  % (sampleID))
		
		outFile.write('\n#############################################################################################\n')


		outFile.write('#CONTROLS ONLY\n')

		outFile.write('#Making histograms for Control Probe Intensities\n')
		outFile.write('pdf(file="%s_Control_ProbeIntensities_Histograms.pdf", width = 10,height = 7)\n' % (sampleID))
		outFile.write('par(mfrow=c(2,3), mai=c(1,0.5,0.5,0.25), pin=c(2.25,2.25))\n')
		outFile.write('hist(cred, main = "Control Processed Probe Intensities (Red Only) for %s", xlab="Probe Intensity", col="red2", cex.main=0.75, cex.lab=0.75, cex.axis=0.75)\n' % (sampleID))
		outFile.write('hist(log2cred, main = "Control Log2 Probe Intensities (Red Only) for %s", xlab="Log2 Probe Intensity", xlim=c(0,20), ylim=c(0,1200), col="red2", cex.main=0.75, cex.lab=0.75, cex.axis=0.75)\n' % (sampleID))
		outFile.write('hist(mean_corr_log2cred, main = "Mean Corrected Control Log2 Probe Intensities (Red Only) for %s", xlab="Log2 Probe Intensity", xlim=c(0,2), ylim=c(0,1200), col="red2", cex.main=0.75, cex.lab=0.75, cex.axis=0.75)\n' % (sampleID))
		outFile.write('hist(cgreen, main = "Control Processed Probe Intensities (Green Only) for %s",xlab="Probe Intensity", col="limegreen", cex.main=0.75, cex.lab=0.75, cex.axis=0.75)\n' % (sampleID))
		outFile.write('hist(log2cgreen, main = "Control Log2 Probe Intensities (Green Only) for %s", xlab="Corrected Log2 Probe Intensity", xlim=c(0,20), ylim=c(0,1200), col="limegreen", cex.main=0.75, cex.lab=0.75, cex.axis=0.75)\n' % (sampleID))
		outFile.write('hist(mean_corr_log2cgreen, main = "Mean Corrected Control Log2 Probe Intensities (Green Only) for %s", xlab="Corrected Log2 Probe Intensity", xlim=c(0,2), ylim=c(0,1200), col="limegreen", cex.main=0.75, cex.lab=0.75, cex.axis=0.75)\n' % (sampleID))
		outFile.write('dev.off()\n')

		outFile.write('#Making boxplots for Log2 Control Probe Intensities\n')
		outFile.write('pdf(file="%s_ControlProbeIntensities_BoxPlots.pdf", width = 10,height = 7)\n' % (sampleID))
		outFile.write('par(mfrow=c(1,3), mai=c(1,0.5,0.5,0.25), pin=c(2.25,5))\n')
		outFile.write('dataframe1 <- data.frame(cred,cgreen)\n')
		outFile.write('colnames(dataframe1) <- c("Red", "Green")\n')
		outFile.write('dataframe2 <- data.frame(log2cred, log2cgreen)\n')
		outFile.write('colnames(dataframe2) <- c("Log2Red", "Log2Green")\n')
		outFile.write('boxplot(dataframe1, ylab="Probe Intensity", xlab="Signals", main="Control Probe Intensities for %s", col=(c("red2","limegreen")), cex.main=0.75, cex.lab=0.75, cex.axis=0.75)\n' % (sampleID))
		outFile.write('boxplot(dataframe2, ylab="Log2 Probe Intensity", xlab="Signals", main="Control Log2 Probe Intensities for %s", col=(c("red2","limegreen")), ylim=c(0,20), cex.main=0.75, cex.lab=0.75, cex.axis=0.75)\n' % (sampleID))
		outFile.write('dataframe3 <- data.frame(log2cred, log2cgreen, mean_corr_log2cred, mean_corr_log2cgreen)\n')
		outFile.write('colnames(dataframe3) <- c("Log2 Red", "Log2 Green", "Mean Corr. \\nLog2 Red", "Mean Corr. \\nLog2 Green")\n')
		outFile.write('boxplot(dataframe3, ylab="Log2 Probe Intensity", xlab="Signals", main="Control Log2 Probe Intensities for %s", col=(c("red2","limegreen")), ylim=c(0,20), cex.main=0.75, cex.lab=0.75, cex.axis=0.75)\n' % (sampleID))
		outFile.write('dev.off()\n')

		outFile.write('\n###############################################################################################\n')
		outFile.write('#NONCONTROLS ONLY\n')
		outFile.write('#Histograms for Raw Noncontrol Probe Intensities\n')
		outFile.write('pdf(file="%s_Noncontrol_ProbeIntensities_Histograms.pdf", width = 10,height = 7)\n' % (sampleID))
		outFile.write('par(mfrow=c(2,3), mai=c(1,0.5,0.5,0.25), pin=c(2.25,2.25))\n')
		outFile.write('hist(nred, main = "Noncontrol Processed \\nProbe Intensities (Red Only) for %s", xlab="Probe Intensity", col="red2", cex.main=0.75, cex.lab=0.75, cex.axis=0.75,ylim=c(0,150000), xlim=c(0,175000), breaks=50)\n' % (sampleID))
		outFile.write('hist(log2nred, main = "Noncontrol Log2 \\nProbe Intensities (Red Only) for %s", xlab="Log2 Probe Intensity", xlim=c(0,20), ylim=c(0,30000), col="red2", cex.main=0.75, cex.lab=0.75, cex.axis=0.75, breaks=50)\n' % (sampleID))
		outFile.write('hist(mean_corr_log2nred, main = "Mean Corrected Noncontrol \\nLog2 Probe Intensities (Red Only) for %s", xlab="Log2 Probe Intensity", xlim=c(0,2), ylim=c(0,30000), col="red2", cex.main=0.75, cex.lab=0.75, cex.axis=0.75, breaks=50)\n' % (sampleID))
		outFile.write('hist(ngreen, main = "Noncontrol Processed \\nProbe Intensities (Green Only) for %s",xlab="Probe Intensity", col="limegreen", cex.main=0.75, cex.lab=0.75, cex.axis=0.75,ylim=c(0,150000), xlim=c(0,175000), breaks=50)\n' % (sampleID))
		outFile.write('hist(log2ngreen, main = "Noncontrol Log2 \\nProbe Intensities (Green Only) for %s", xlab="Corrected Log2 Probe Intensity", xlim=c(0,20), ylim=c(0,30000), col="limegreen", cex.main=0.75, cex.lab=0.75, cex.axis=0.75, breaks=50)\n' % (sampleID))
		outFile.write('hist(mean_corr_log2ngreen, main = "Mean Corrected Noncontrol \n\nLog2 Probe Intensities (Green Only) for %s", xlab="Corrected Log2 Probe Intensity", xlim=c(0,2), ylim=c(0,30000), col="limegreen", cex.main=0.75, cex.lab=0.75, cex.axis=0.75, breaks=50)\n' % (sampleID))
		outFile.write('dev.off()\n')

		outFile.write('\n#Making boxplots for Log2 Noncorrected Noncontrol Probe Intensities\n')
		outFile.write('pdf(file="%s_NoncontrolProbeIntensities_BoxPlots.pdf", width = 10,height = 7)\n' % (sampleID))
		outFile.write('par(mfrow=c(1,3), mai=c(1,0.5,0.5,0.25), pin=c(2.25,5))\n')
		outFile.write('dataframe1 <- data.frame(nred,ngreen)\n')
		outFile.write('colnames(dataframe1) <- c("Red", "Green")\n')
		outFile.write('dataframe2 <- data.frame(log2nred, log2ngreen)\n')
		outFile.write('colnames(dataframe2) <- c("Log2Red", "Log2Green")\n')
		outFile.write('boxplot(dataframe1, ylab="Probe Intensity", xlab="Signals", main="Noncontrol Probe Intensities for %s", col=(c("red2","limegreen")), cex.main=0.75, cex.lab=0.75, cex.axis=0.75)\n' % (sampleID))
		outFile.write('boxplot(dataframe2, ylab="Log2 Probe Intensity", xlab="Signals", main="Noncontrol Log2 Probe Intensities for %s", col=(c("red2","limegreen")), ylim=c(0,20), cex.main=0.75, cex.lab=0.75, cex.axis=0.75)\n' % (sampleID))
		outFile.write('dataframe3 <- data.frame(log2nred, log2ngreen, mean_corr_log2nred, mean_corr_log2ngreen)\n')
		outFile.write('colnames(dataframe3) <- c("Log2 Red", "Log2 Green", "Mean Corr. \nLog2 Red", "Mean Corr. \nLog2 Green")\n')
		outFile.write('boxplot(dataframe3, ylab="Log2 Probe Intensity", xlab="Signals", main="Noncontrol Log2 Probe Intensities for %s", col=(c("red2","limegreen")), ylim=c(0,20), cex.main=0.75, cex.lab=0.75, cex.axis=0.75)\n' % (sampleID))
		outFile.write('dev.off()\n')

		outFile.write('\n# Mean Corrected Red/Green Histograms and Boxplots for Noncontrol probes only\n')
		outFile.write('pdf(file = "%s_MeanCorrected_NonControls.pdf")\n' % (sampleID))
		outFile.write('par(mfrow=c(2,2), mai=c(1,0.5,0.5,0.25), pin=c(2.25,2.25))\n')
		outFile.write('hist <- hist(mean_corr_log2nred, main = "Mean Corrected Log2 \\nNoncontrol Probe Intensities (Red Only) for %s", xlab="Log2 (Red Intensity)", col="red2", xlim=c(0,2), ylim=c(0,30000), breaks=50, cex.main=0.75, cex.lab=0.75, cex.axis=0.75)\n' % (sampleID))
		outFile.write('hist <- hist(mean_corr_log2ngreen, main = "Mean Corrected Log2 \\nNoncontrol Probe Intensities (Green Only) for %s", xlab="Log2 (Green Intensity)", col="limegreen", xlim=c(0,2), ylim=c(0,30000), breaks=50, cex.main=0.75, cex.lab=0.75, cex.axis=0.75)\n' % (sampleID))
		outFile.write('hist <- hist(mean_corr_log2ngreen, main="Merged Mean Corrected Log2 \\nNoncontrol Probe Intensities\n(Red and Green) for %s", xlab="Log2 Probe Intensity", col=rgb(0,3/4,1/4,1/4), xlim=c(0,2), ylim=c(0,30000), breaks=50, cex.main=0.75, cex.lab=0.75, cex.axis=0.75)\n' % (sampleID))
		outFile.write('hist <- hist(mean_corr_log2nred, xlim=c(0,2), ylim=c(0,30000), breaks=50, col=rgb(1,0,0,5/8),add=T)\n')
		outFile.write('dataframe3 <- data.frame(log2nred, log2ngreen, mean_corr_log2nred, mean_corr_log2ngreen)\n')
		outFile.write('colnames(dataframe3) <- c("Log2Red", "Log2Green", "Corr\\nLog2Red", "Corr\\nLog2Green")\n')
		outFile.write('boxplot(dataframe3, ylab="Log2 Probe Intensity", xlab="Signals", main="Noncontrol Log2 \\nProbe Intensities for %s", col=(c("red2","limegreen")), ylim=c(0,20), cex.main=0.65, cex.lab=0.65, cex.axis=0.65)\n' % (sampleID))

		outFile.write('\n############################################################\n')

		outFile.write('#### MERGED DATA PLOTS\n') 
		outFile.write('# NONCONTROLS VS CONTROLS\n')
		outFile.write('pdf(file="%s_MERGED_ProbeIntensities_Histograms.pdf", paper="letter",width = 7,height = 10)\n' % (sampleID))
		outFile.write('par(mfrow=c(2,2), mai=c(1,0.5,0.5,0.25), pin=c(2.75,2.75))\n')
		outFile.write('#Controls\n')
		outFile.write('hist(log2cgreen, main = "Merged CONTROL \\nLog2 Probe Intensities for %s", xlab="Log2 Probe Intensity", xlim=c(0,20), ylim=c(0,800), col=rgb(0,3/4,1/4,1/4), cex.main=0.75, cex.lab=0.75, cex.axis=0.75, breaks=50)\n' % (sampleID))
		outFile.write('hist(log2cred, xlim=c(0,20), ylim=c(0,800), col=rgb(1,0,0,5/8), cex.main=0.75, cex.lab=0.75, cex.axis=0.75, breaks=50, add=T)\n')
		outFile.write('hist(mean_corr_log2cgreen, main = "Merged Mean Corrected CONTROL \\nLog2 Probe Intensities for %s", xlab="Corrected Log2 Probe Intensity", xlim=c(0,2), ylim=c(0,1200), col=rgb(0,3/4,1/4,1/4), cex.main=0.75, cex.lab=0.75, cex.axis=0.75, breaks=50)\n' % (sampleID))
		outFile.write('hist(mean_corr_log2cred,  xlim=c(0,2), ylim=c(0,1200), col=rgb(1,0,0,5/8),  cex.main=0.75, cex.lab=0.75, cex.axis=0.75, breaks=50, add=T)\n')
		outFile.write('#Noncontrols\n')
		outFile.write('hist(log2ngreen, main = "Merged NONCONTROL \\nLog2 Probe Intensities for %s", xlab="Log2 Probe Intensities", breaks=50, ylim=c(0,30000), xlim=c(0,20), col=rgb(0,3/4,1/4,1/4),cex.main=0.75, cex.lab=0.75, cex.axis=0.75)\n' % (sampleID))
		outFile.write('hist(log2nred, breaks=50, ylim=c(0,30000), xlim=c(0,20), col=rgb(1,0,0,5/8),cex.main=0.75, cex.lab=0.75, cex.axis=0.75,add=T)\n')
		outFile.write('hist(mean_corr_log2ngreen, main = "Merged Mean Corrected NONCONTROL \\nLog2 Probe Intensities for %s", xlab="Corrected Log2 Probe Intensity", xlim=c(0,2), ylim=c(0,30000), col=rgb(0,3/4,1/4,1/4), breaks=50, cex.main=0.75, cex.lab=0.75, cex.axis=0.75)\n' % (sampleID))
		outFile.write('hist(mean_corr_log2nred, xlim=c(0,2), ylim=c(0,30000), col=rgb(1,0,0,5/8), breaks=50, cex.main=0.75, cex.lab=0.75, cex.axis=0.75, add=T)\n')
		outFile.write('dev.off()\n')

		outFile.write('#TOTAL VERSUS CONTROLS\n')
		outFile.write('par(mfrow=c(2,2))\n')
		outFile.write('pdf(file="%s_TotalVControls_ProbeIntensities_Histograms.pdf", paper="letter",width = 7,height = 10)\n' % (sampleID))
		outFile.write('par(mfrow=c(2,2), mai=c(1,0.5,0.5,0.25), pin=c(2.75,2.75))\n')
		outFile.write('hist(log2%sred, col="indianred1", xlim=c(0,20), ylim=c(0,25000), breaks=25, main="Log2 Red Intensities of Total (Red)\n and Control (Gray) Red Probes for %s", xlab = "Log2 Probe Intensities",cex.main=0.75, cex.lab=0.75, cex.axis=0.75)\n' % (sampleID,sampleID))
		outFile.write('hist(log2cred, col="grey86", xlim=c(0,20), ylim=c(0,25000), cex.main=0.75, cex.lab=0.75, cex.axis=0.75, breaks=25, add=T)\n')
		outFile.write('hist(mean_corr_log2%sred, breaks=25, xlim=c(0,2), ylim=c(0,25000), main = "Mean Corrected Log2 Intensities \\nof Total (Red) and Control (Gray) Probes for %s", xlab="Mean Corrected Log2 Red Intensity", cex.main=0.75, cex.lab=0.75, cex.axis=0.75, col="indianred1")\n' % (sampleID,sampleID))
		outFile.write('hist(mean_corr_log2cred, breaks=25, xlim=c(0,2), ylim=c(0,25000), main = "Mean Corrected Log2 Intensities \\nof Total (Red) and Control (Gray) Probes for %s", xlab="Mean Corrected Log2 Red Intensity", cex.main=0.75, cex.lab=0.75, cex.axis=0.75, col="grey86", add=T)\n' % (sampleID))
		outFile.write('hist(log2%sgreen, col="mediumseagreen", xlim=c(0,20), ylim=c(0,25000), breaks=25, main="Log2 Red Intensities \\nof Total (Red) and Control (Gray) Red Probes for %s", xlab = "Log2 Probe Intensities",cex.main=0.75, cex.lab=0.75, cex.axis=0.75)\n' % (sampleID,sampleID))
		outFile.write('hist(log2cgreen, col="grey86", xlim=c(0,20), ylim=c(0,25000), cex.main=0.75, cex.lab=0.75, cex.axis=0.75, breaks=25, add=T)\n')
		outFile.write('hist(mean_corr_log2%sgreen, breaks=25, xlim=c(0,2), ylim=c(0,25000), main = "Mean Corrected Log2 Intensities \\nof Total (Green) and Control (Gray) Probes for %s", xlab="Mean Corrected Log2 Red Intensity", cex.main=0.75, cex.lab=0.75, cex.axis=0.75, col="mediumseagreen")\n' % (sampleID,sampleID))
		outFile.write('hist(mean_corr_log2cgreen, breaks=25, xlim=c(0,2), ylim=c(0,25000), cex.main=0.75, cex.lab=0.75, cex.axis=0.75, col="grey86", add=T)\n')
		outFile.write('dev.off()\n')


		outFile.write('%s <- read.table("%sProbeIntensityData",header=TRUE)\n' % (sampleID, sampleID))
		outFile.write('colnames(%s) <- c("Sample", "Feature", "ProbeType", "ProbeID", "Log2Ratio", "Red", "Green", "Chromosome", "Start", "End")\n' % (sampleID))
		outFile.write('sampleID<-%s[1]\n' % (sampleID))
		outFile.write('featureNum <- %s[2]\n' % (sampleID))
		outFile.write('probeType<-as.factor(unlist(%s[3]))\n' % (sampleID))
		outFile.write('probeID<-as.factor(unlist(%s[4]))\n' % (sampleID))
		outFile.write('log2 <- sapply(%s[5],as.numeric)\n' % (sampleID))
		outFile.write('rawRed <- sapply(%s[6],as.numeric)\n' % (sampleID))
		outFile.write('rawGreen <- sapply(%s[7],as.numeric)\n' % (sampleID))
		outFile.write('chrom <- as.factor(unlist(%s[8]))\n' % (sampleID))
		outFile.write('start<-sapply(%s[9],as.numeric)\n' % (sampleID))
		outFile.write('end<-sapply(%s[10],as.numeric)\n' % (sampleID))
		outFile.write('log2red <- sapply((log2(rawRed)),as.numeric)\n')
		outFile.write('log2green <- sapply((log2(rawGreen)),as.numeric)\n')
		outFile.write('mean_corr_log2red <- log2red/mean_red%scontrols\n' % (sampleID))
		outFile.write('mean_corr_log2green <- log2green/mean_green%scontrols\n' % (sampleID))
		outFile.write('mean_corr_log2ratio <- mean_corr_log2red/mean_corr_log2green\n')
		outFile.write('export_dataframe <- data.frame(sampleID, featureNum, probeType, probeID, log2, rawRed, rawGreen, chrom, start, end, log2red, log2green, mean_corr_log2red, mean_corr_log2green, mean_corr_log2ratio, zscore)\n')
		outFile.write('colnames(export_dataframe) <- c("Sample", "Feature", "ProbeType", "ProbeID", "Log2Ratio", "Red", "Green", "Chromosome", "Start", "End", "Log2_Red", "Log2_Green", "MeanCorrected_Log2Red", "MeanCorrected_Log2Green", "MeanCorrected_Log2Ratio", "Zscore_MeanCorr_Log2Ratio")\n') 		
		writeout = '/home/ampend/kidd-lab/ampend-projects/CGH_Array_Design_Dog/CGH_Array_Analysis/results/' + sampleID + '_MeanCorrectedValues.txt'
		outFile.write('write.table(export_dataframe, "%s", sep="\t", quote=FALSE, row.names=F)\n' % (writeout)) 
	
outfile = options.directory + 'TOTAL_' + 'Rcmds.r'
print '\n\nWriting commands for across-array boxplots to:', outfile
print 'Includes boxplots for samples: '
outFile = open(outfile,'w')
#max = 4 
lineNum = 0

outFile.write('pdf(file="TOTAL_Uncorrected_Boxplots_AgilentCustomCNVControlsNonControls.pdf", paper="letter", width = 7,height = 10)\n') 
outFile.write('par(mfrow=c(2,2), mai=c(.75,0.5,0.5,0.15), pin=c(3,3))\n')
infile = open(options.input, 'r')

for line in infile:
	line = line.rstrip()
	line = line.split()
	lineNum+=1 
	
	if lineNum > 1: #Skips header line
		sampleID = line[0]
		print sampleID
		outFile.write('%s <- read.table("%s%sProbeIntensityData",header=TRUE)\n' % (sampleID, options.directory, sampleID))
		outFile.write('colnames(%s) <- c("Sample", "Feature", "ProbeType", "ProbeID", "Log2Ratio", "Red", "Green", "Chromosome", "Start", "End")\n'  % (sampleID))
		outFile.write('%stype<-as.factor(unlist(%s[3]))\n' % (sampleID, sampleID))
		outFile.write('%stype<-as.factor(unlist(%s[3]))\n' % (sampleID, sampleID))
		outFile.write('%slog2 <-sapply(%s[5],as.numeric)\n' % (sampleID, sampleID))
		outFile.write('%sDataframe <- data.frame(%stype, %slog2)\n' % (sampleID, sampleID, sampleID))
		outFile.write('boxplot(%slog2 ~ %stype, %sDataframe, main = "Boxplot of Uncorrected Log2 \\nRatios (Red:Green) for %s", ylab = "log2 Values", cex.main=0.75, cex.lab=0.75, cex.axis=0.75)\n' % (sampleID, sampleID, sampleID, sampleID))
outFile.write('dev.off()\n')
		
		
		
		
		
		


		
		#cmd = 'Rscript %s' % (outfile)
		#print cmd
		#genutils.runCMD(cmd)


#cmd = 'for file in \$(ls ../results/*.r); do Rscript \$file; done;' 
#print cmd
#genutils.runCMD(cmd)

