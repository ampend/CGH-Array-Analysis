# make_commands.py
# 2015-07-22
# writes commands file (.cmds) that run-jobs.pbs will use to process the arrays

import genutils
from optparse import  OptionParser

###############################################################################
USAGE = """
python make_commands.py 	--input < README FILE > 
							--directory < Current working directory with scripts >
							--probes < Probe BED file for extractFeatures.py to use > 

input == README file
directory == Directory to write output files 
probes == Probe BED file for extractFeatures.py to use

"""

parser = OptionParser(USAGE)
parser.add_option('--input',dest='input', help = 'input Agilent Feature Extraction output')
parser.add_option('--directory',dest='directory', help = 'Directory to write output files')
parser.add_option('--probes',dest='probes', help = 'Probe BED file for extractFeatures.py to use')
(options, args) = parser.parse_args()

if options.input is None:
    parser.error('input probe file name not given')
if options.directory is None:
    parser.error('output probe file name not given')
if options.probes is None:
    parser.error('output probe file name not given')

###############################################################################

#Example command line
#python extractFeatures.py --input ../inputData/RawScanData/257548210013_2015-07-10_09-54_Area1_532_635_CGH_1200_Jun14_1_4.txt --sample VN37_5027 --directory /home/ampend/kidd-lab/ampend-projects/CGH_Array_Design_Dog/CGH_Array_Analysis/results/ --probes ../inputData/ProbeBEDFiles/Mappable/TotalProbes_FINAL_Coordinates.bed

cmdfile = options.directory + 'extractFeatures.cmds'
cmdFile = open(cmdfile, 'w')

readmeFile = open(options.input, 'r')

lineNumber = 0

for line in readmeFile:
	line = line.rstrip()
	line = line.split()
	lineNumber += 1
	
	if lineNumber > 1:
		sample = line[0]
		input = line[8]
		
		cmdFile.write('python %sextractFeatures.py --input %s --sample %s --directory ../results/ --probes %s\n' % (options.directory, input, sample, options.probes))

cmd = 'qsub %srun-jobs.pbs' % options.directory
print cmd
genutils.runCMD(cmd)

