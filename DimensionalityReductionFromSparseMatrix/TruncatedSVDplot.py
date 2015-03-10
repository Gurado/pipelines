#!/bin/python
######################################
# calculate using sklearn TruncatedSVD
#
# Author: Fabian Buske (13/01/2015)
######################################

from scipy.sparse import lil_matrix
from scipy.sparse import csr_matrix
import numpy
import regex, os, sys, errno, re
import argparse
from sklearn.decomposition import TruncatedSVD
import fileinput
import datetime
import gzip

######################################
# Timestamp
######################################
def timeStamp():
    return datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S').format()


def createFragmentResolution():
	''' 
	creates one interval tree for quick lookups
	returns 
	    fragmentsLookupTable[fragmentId] = [tuple(chrom, fragmentMidPoint)]
	'''

	if (args.verbose):
		print >> sys.stdout, "- %s START   : populate lookup table with given resolution for chromosomes matching pattern %s" % (timeStamp(), args.chromPattern)

	fragmentsCount = 0
	fragmentsLookupTable = {}
	 
	for line in fileinput.input([args.chromSizes]):
		chrom=line.split("\t")[0]
		# check if chromosome needs to be filtered out or not
		if (args.chromPattern != "" and not re.match(args.chromPattern, chrom)):
			# skip this one
			if (args.veryverbose):
				print "skipping pattern %s" % (line)
			continue
		chromlen=int(line.split("\t")[1])

		for i in range(0, chromlen, args.resolution):
			start=i
			end=min(i+ args.resolution, chromlen)
			fragmentsLookupTable[tuple([chrom, int(0.5*(start+end))])] = fragmentsCount
			fragmentsCount += 1
	
	if (args.verbose):
		print >> sys.stdout, "- %s FINISH  : counted %d fragments" % (timeStamp(), fragmentsCount)

	return [ fragmentsLookupTable, fragmentsCount ]

def readFileIntoSparseMatrix(fragmentsLookupTable, fragmentsCount):

	A = lil_matrix((len(args.contactsCountMatrices), fragmentsCount * fragmentsCount), dtype='i')
	c = 0 
	S = []
	for contactCountsFile in args.contactsCountMatrices:

		S += [os.path.basename(contactCountsFile).split(".")[0]]

		if (args.verbose):
			print >> sys.stdout, "- %s START   : reading file %s" % (timeStamp(), contactCountsFile)

		if contactCountsFile.endswith('.gz'):
			infile=gzip.open(contactCountsFile,'r')
		else:
			infile=open(contactCountsFile,'r')


		for line in infile:
			l = line.split()
			ch1,mid1,ch2,mid2,metric,qvalue = [l[0],int(l[1]),l[2],int(l[3]),float(4),float(5)]
			
			try:
				if (float(qvalue) > args.threshold):
					continue
			except:
				# probably the header line
				continue

			# skip irrelevant entries
			if (metric == 0):
				continue

			# skip trand counts if requested
			if (args.cis and ch1 != ch2):
				continue
			
			fragment1 = tuple([ch1,int(mid1)])
			fragment2 = tuple([ch2,int(mid2)])


			if (not fragmentsLookupTable.has_key(fragment1) or not fragmentsLookupTable.has_key(fragment2)):
				if (args.veryverbose):
					print "[NOTE] fragment %s or %s not detected in map, is the correct resolution specified?" % (ch1,ch2)
				continue

			# keep this symmetric matrix as sparse as possible by just filling in the top triangle
			if (fragment1 < fragment2):
				A[c, fragmentsLookupTable[fragment1] * fragmentsCount + fragmentsLookupTable[fragment2]] = metric
			else:
				A[c, fragmentsLookupTable[fragment2] * fragmentsCount +  fragmentsLookupTable[fragment1]] = metric

		c += 1 

		if (args.verbose):
			print >> sys.stdout, "- %s FINISH  : reading file" % (timeStamp())

	return (S,A.tocsr())

def explainVariance(M):

	if (args.verbose):
		print >> sys.stdout, "- %s START   : calculating SVD" % (timeStamp())

	svd = TruncatedSVD(n_components=3, random_state=42)
	svd.fit(M)
	N = svd.transform(M)
	if (args.verbose):
		print >> sys.stdout, "- %s FINISH  : calculating SVD" % (timeStamp())

	return (N,svd.explained_variance_ratio_,svd.explained_variance_ratio_.sum())


def plotVariance(S, N, explainedRatio, totalExplained):

	x = ["%.2f" % number for number in numpy.transpose(N)[0]]
	y = ["%.2f" % number for number in numpy.transpose(N)[1]]
	z = ["%.2f" % number for number in numpy.transpose(N)[2]]
	l = [args.labels.split(",")]

	jo = open("+args.outdir+'/'+args.prefix+".json',"w")

	jo.write("""{
		explained : {
			x : %.2f,
			y : %.2f,
			z : %.2f
		},
		experiments : { """ % (explainedRatio[0]*100, explainedRatio[1]*100, explainedRatio[2]*100))
	
	arr = []
	for i in range(length(l)):
		arr += [('"prefix" : "%s", x: %.2f, y: %.2f,z: %.2f ' % (l[i],x[i],y[i],z[i]))] 
	jo.write(",\n".join(arr))
	jo.write("""
		}
}	
	""")

def plotVariance(S, N, explainedRatio, totalExplained):

	f = open(args.outdir+'/'+args.imagename+".R","w")
	x = ["%.2f" % number for number in numpy.transpose(N)[0]]
	y = ["%.2f" % number for number in numpy.transpose(N)[1]]
	z = ["%.2f" % number for number in numpy.transpose(N)[2]]

	f.write("library(ggplot2)\n")
	if (args.labels != ""):
		f.write("data <- data.frame(experiment=c(\"%s\"), " % ('","'.join(args.labels.split(","))))
	else:
		f.write("data <- data.frame(experiment=c(\"%s\"), " % ('","'.join(S)))
	f.write(" x=c(%s), " % (','.join(x)))
	f.write(" y=c(%s)," % (','.join(y)))
	f.write(" z=c(%s)" % (','.join(z)))
	s=[]
	if (args.groups != ""):
		f.write(", groups=c('%s')" % ("','".join(args.groups.split(','))))
		s+=["color=groups"]
	else:
		s+=["color=experiment"]
	if (args.cutter != ""):
		f.write(", cutter=c('%s')" % ("','".join(args.cutter.split(','))))
		s+=["shape=cutter"]

	f.write(")\n")
	f.write("g1<- ggplot(data, aes(x,y))")
	f.write(" + geom_point(aes(%s), alpha = 0.5, size=5)" %(','.join(s)))
	f.write(" + xlab('Dim 1 (%.1f%% variance explained)')" % (explainedRatio[0]*100))
	f.write(" + ylab('Dim 2 (%.1f%% variance explained)')" % (explainedRatio[1]*100))
	f.write(" + geom_text(aes(label=experiment, color=experiment), size=3, vjust=3)")
	f.write(" + xlim(c(min(data$x)-0.1*(max(data$x)-min(data$x)),max(data$x)+0.1*(max(data$x)-min(data$x))))")
	f.write(" + ylim(c(min(data$y)-0.1*(max(data$y)-min(data$y)),max(data$y)+0.1*(max(data$y)-min(data$y))))")
	f.write("\n")

	f.write("g2<- ggplot(data, aes(y,z))")
	f.write(" + geom_point(aes(%s), alpha = 0.5, size=5)" %(','.join(s)))
	f.write(" + xlab('Dim 2 (%.1f%% variance explained)')" % (explainedRatio[1]*100))
	f.write(" + ylab('Dim 3 (%.1f%% variance explained)')" % (explainedRatio[2]*100))
	f.write(" + geom_text(aes(label=experiment, color=experiment), size=3, vjust=3)")
	f.write(" + xlim(c(min(data$y)-0.1*(max(data$y)-min(data$y)),max(data$y)+0.1*(max(data$y)-min(data$y))))")
	f.write(" + ylim(c(min(data$z)-0.1*(max(data$z)-min(data$z)),max(data$z)+0.1*(max(data$z)-min(data$z))))")
	f.write("\n")

	f.write("pdf('"+args.outdir+'/'+args.prefix+".pdf', width=10, height=8)\n")
	f.write("g1\n")
	f.write("g2\n")
	f.write("dev.off()\n")

	f.close()

def main():
	''' main method 
	    parsing the parameters and options
	'''
	global args
	parser = argparse.ArgumentParser(description='Reads a list of significant interaction sparse matrix files and performs TruncatedSVD')
	parser.add_argument('chromSizes', type=str, help='chomosome sizes')
	parser.add_argument('contactsCountMatrices', metavar="contactsCountMatrices", type=str, nargs='+', help='sparse interaction matrices')
	parser.add_argument("-o", '--outdir', dest='outdir', type=str, default="./",
						help='output location')
	parser.add_argument("-O", '--prefix', dest='prefix', type=str, default="truncatedSVD",
						help='prefix of the output files')
	parser.add_argument("-r", "--resolution", type=int, dest="resolution", default=1000000, 
						help="size of a fragment in bp if no genomeFragmentFile is given")
	parser.add_argument("-C", "--chrompattern", type=str, dest="chromPattern", default="", 
						help="pattern of chromosomes to filter for [default all]")
	parser.add_argument("-t", "--threshold", type=float, dest="threshold", default=0.01, 
						help="q-value threshold used to filter data [default 0.05]")
	parser.add_argument("--cis", action="store_true", help="consider cis interactions only [default all]")
	
	parser.add_argument("-g", "--groups", type=str, dest="groups", default="", 
						help="group list of count matrices via comma-separated list, e.g. 1,1,3,4")
	parser.add_argument("-l", "--labels", type=str, dest="labels", default="", 
						help="text labels for experiments used in plots, supplied via comma-separated list, e.g. PrEC_I,PrEC_II,LNCaP_I,LNCaP_II")
	parser.add_argument("-c", "--cutter", type=str, dest="cutter", default="", 
						help="restriction enzymes used in matrices, supplied via comma-separated list, e.g. HindIII,NcoI,HindIII,NcoI")
	parser.add_argument("--plot", action="store_true")
	parser.add_argument("--verbose", action="store_true")
	parser.add_argument("--veryverbose", action="store_true")

	parser.add_argument("--quiet", action="store_true")

	args = parser.parse_args()

	# try:
	# 	os.makedirs(args.outdir)
	# except OSError as exc: # Python >2.5
	# 	if exc.errno == errno.EEXIST and os.path.isdir(args.outdir):
	# 		pass
	# 	else: raise

	if (args.resolution < 1):
		parser.error("[ERROR] resolution must be a positive integer, was :"+str(args.resolution))
		sys.exit(1)
	elif (args.chromSizes == "" or not os.path.isfile(args.chromSizes)):
		parser.error("[ERROR] chromSizes not given or not existing, was :"+str(args.chromSizes))
		sys.exit(1)

	[ fragmentsLookupTable, fragmentsCount ] = createFragmentResolution()


	(S, A) = readFileIntoSparseMatrix(fragmentsLookupTable, fragmentsCount)

	(N, explainedRatio, totalExplained) = explainVariance(A)

	outputJSON(S, N, explainedRatio, totalExplained)

	if (args.plot):
		plotVariance(S, N, explainedRatio, totalExplained)

if __name__ == "__main__":
		main()
