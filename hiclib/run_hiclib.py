import os
import logging
from hiclib import mapping
from mirnylib import h5dict, genome
from hiclib import fragmentHiC
from optparse import OptionParser
import matplotlib.pyplot as plt
import numpy as np
from mirnylib import plotting
from hiclib import binnedData


# manage option and arguments processing
def main():
	global options
	global args
	usage = '''usage: %prog [options] read1.fastq read2.fastq

takes two hic fastq files and runs the hiclib pipeline on it
	'''
	parser = OptionParser(usage)
	parser.add_option("-q", "--quiet", action="store_false", dest="verbose", default=True,
					help="don't print status messages to stdout")
	parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
					help="print status messages to stdout")
	parser.add_option("-e", "--restrictionEnzyme", type="string", dest="enzyme", default="", 
					help="Name of the restriction enzyme, e.g. BglII")
	parser.add_option("-n", "--experimentName", type="string", dest="experiment", default="", 
					help="Name of the experiment")
	parser.add_option("-b", "--bowtie", type="string", dest="bowtie", default="", 
					help="location of bowtie [default: %default]")
	parser.add_option("-g", "--genome", type="string", dest="genome", default="", 
					help="genome in fasta format [default: %default]")
	parser.add_option("-i", "--index", type="string", dest="index", default="", 
					help="location of genome index")
	parser.add_option("-o", "--outputDir", type="string", dest="outputDir", default="", 
					help="output directory [default: %default]")
	parser.add_option("-c", "--cpu", type="int", dest="cpu", default=1, 
					help="number of cpus to use [default: %default]")
	parser.add_option("-t", "--tmpDir", type="string", dest="tmpDir", default="/tmp", 
					help="directory for temp files [default: %default]")
	parser.add_option("-r", "--sra-reader", type="string", dest="sra", default="fastq-dump", 
					help="location of sra reader fastq-dump in case input is SRA [default: %default]")
	
	(options, args) = parser.parse_args()
	if (len(args) != 2):
		parser.print_help()
		parser.error("incorrect number of arguments")

	if (options.genome == ""):
		print >> sys.stderr, "[ERROR] Please specify the location of the reference genome in fasta format"
		sys.exit(1)

	if (options.index == ""):
		print >> sys.stderr, "[ERROR] Please specify the location of the bowtie2 index for the reference genome"
		sys.exit(1)
		
	if (options.enzyme == ""):
		print >> sys.stderr, "[ERROR] Please specify the restriction enzyme (supported enzymes: http://www.biopython.org/DIST/docs/api/Bio.Restriction-module.html)"
		sys.exit(1)

	if (options.experiment == ""):
		print >> sys.stderr, "[ERROR] Please provide a name for the experiment, e.g. [Cellline]_[Enzymename]_[Replica]"
		sys.exit(1)
	
	if (options.output != ""): 
		options.output += "/"
		

	if (options.verbose):
		print >> sys.stderr, "read1.fastq  : %s" % (args[0])
		print >> sys.stderr, "read2.fastq  : %s" % (args[1])
			
	process()


def mapFastq(fastq):
	global options
	global args

	fileName, fileExtension = os.path.splitext(fastq)
	bamOutput = options.output+'/'+fileName.split('/')[-1]+'.bam'
	
	if (fileExtension == 'sra'):
		if (options.verbose):
			print >> sys.stdout, "Map short read archive %s utilizing %s" % (fastq, options.sra)
			
		mapping.iterative_mapping(
		    bowtie_path=options.bowtie,
		    bowtie_index_path=options.index,
		    fastq_path=fastq,
		    out_sam_path=bamOutput,
		    min_seq_len=25,
		    len_step=5,
		    seq_start=0,
		    seq_end=75,
		    nthreads=options.cpu,
		    temp_dir=options.tmpDir, 
		    bowtie_flags='--very-sensitive',
		    bash_reader=options.sra+' -Z')
	
	else:
		if (options.verbose):
			print >> sys.stdout, "Map fastq %s" % (fastq)
		
		mapping.iterative_mapping(
		    bowtie_path=options.bowtie,
		    bowtie_index_path=options.index,
		    fastq_path=fastq,
		    out_sam_path=bamOutput,
		    min_seq_len=25,
		    len_step=5,
		    seq_start=0,
		    seq_end=75,
		    nthreads=options.cpu,
		    temp_dir=options.tmpDir, 
		    bowtie_flags='--very-sensitive')
		    
	return bamOutput
		   

def collectMappedReads(bam_read1, bam_read2, mapped_reads, genome_db):
	global options
	global args
	
	mapping.parse_sam(
	    sam_basename1=bam_read1,
	    sam_basename2=bam_read2,
	    out_dict=mapped_reads,
	    genome_db=genome_db, 
	    enzyme_name=options.enzyme)

def filterFragments(genome_db):
	'''
	Filter the data at the level of individual restriction fragments

	The following reads are remove from the dataset:

	- the reads that start within the 5 bp range from the restriction site
	- the identical read pairs, with both ends starting at exactly the same positions
	- the reads coming from extremely large and extremely small restriction fragments (length > 10^5 bp or length < 100 bp)
	- the reads coming from the top 0.5% most frequently detected restriction fragments

	The rationale behind each of the filters is discussed in the hiclib publication. The API documentation contains the description of the filters.
	'''
	
	fragments = fragmentHiC.HiCdataset(
	    filename=options.output+'/fragment_dataset.hdf5',
	    genome=genome_db,
	    maximumMoleculeLength=500,
	    mode='w')
	
	# Load the parsed reads into the HiCdataset. The dangling-end filter is applied
	# at this stage, with maximumMoleculeLength specified at the initiation of the 
	# object.
	fragments.parseInputData(
	    dictLike=options.output+'/mapped_reads.hdf5')
	
	fragments.filterRsiteStart(offset=5)
	fragments.filterDuplicates()
	
	fragments.filterLarge()
	fragments.filterExtreme(cutH=0.005, cutL=0)
	
	fragments.saveHeatmap(options.output+'/heatmap-res-1M.hdf5', resolution=1000000)
	
	return fragments

def iterativeFiltering(genome_db, fragments):
	'''
	Filter the data at the binned level and perform the iterative correction.
	'''
	
	# Read resolution from the dataset.
	raw_heatmap = h5dict.h5dict(options.output+'/heatmap-res-1M.hdf5', mode='r') 
	resolution = int(raw_heatmap['resolution'])
	
	# Create a binnedData object, load the data.
	BD = binnedData.binnedData(resolution, genome_db)
	BD.simpleLoad(options.output+'/heatmap-res-1M.hdf5', options.experiment)

	# Remove the contacts between loci located within the same bin.
	BD.removeDiagonal()
	
	# Remove bins with less than half of a bin sequenced.
	BD.removeBySequencedCount(0.5)
	
	# Remove 1% of regions with low coverage.
	BD.removePoorRegions(cutoff=1)
	
	# Truncate top 0.05% of interchromosomal counts (possibly, PCR blowouts).
	BD.truncTrans(high=0.0005)
	
	# Perform iterative correction.
	BD.iterativeCorrectWithoutSS()

	# Save the iteratively corrected heatmap.
	BD.export(options.experiment, options.output+'/IC-heatmap-res-1M.hdf5')
	
	# Plot the heatmap directly.
	plotting.plot_matrix(np.log(BD.dataDict[options.experiment]))

def process():
	global options
	global args
	
	if (options.verbose):
		print >> sys.stdout, "*** START processing"

	
	logging.basicConfig(level=logging.DEBUG)
	
	if (options.verbose):
		print >> sys.stdout, "**  Create directories"

	if not os.path.exists(options.tmpDir):
	    os.mkdir(options.tmpDir)

	if not os.path.exists(options.outputDir):
	    os.mkdir(options.outputDir)
	
	if (options.verbose):
		print >> sys.stdout, "**  Create data objects"

	mapped_reads = h5dict.h5dict(options.output+'/mapped_reads.hdf5')
	genome_db    = genome.Genome(options.genome)

	if (options.verbose):
		print >> sys.stdout, "**  Map first input file"

	bam_read1 = mapFastq(args[0])

	if (options.verbose):
		print >> sys.stdout, "**  Map second input file"

	bam_read2 = mapFastq(args[1])

	if (options.verbose):
		print >> sys.stdout, "**  Collect mapped reads"
		
	collectMappedReads(bam_read1, bam_read2, mapped_reads, genome_db)
	
	if (options.verbose):
		print >> sys.stdout, "**  Filter fragments"
	
	fragments = filterFragments(genome_db)
	
	if (options.verbose):
		print >> sys.stdout, "**  Iterative filtering of fragments"

	iterativeFiltering(genome_db, fragments)
	
	if (options.verbose):
		print >> sys.stdout, "*** FINISHED processing"
	
	
	
######################################
# main
######################################
if __name__ == "__main__":
	main()
