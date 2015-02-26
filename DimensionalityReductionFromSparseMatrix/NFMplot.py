#!/bin/python
######################################
# calculate using sklearn NMF
#
# Author: Fabian Buske (13/01/2015)
######################################

from scipy.sparse import lil_matrix
from scipy.sparse import coo_matrix
from scipy.sparse import csr_matrix
import numpy
import regex, os, sys, errno, re
import argparse
from sklearn.decomposition import NMF
import fileinput
import datetime
import gzip, itertools
import matplotlib 
from time import time
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

pp = PdfPages('multipage.pdf')


######################################
# Timestamp
######################################
def timeStamp():
    return datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S').format()


def createFragmentResolution():
	''' 
	creates one interval tree for quick lookups
	returns 
	    positionLookupTable[tuple(chrom, fragmentMidPoint)] = [fragmentCount]
		fragmentsLookupTable[fragmentId] = [tuple(chrom, fragmentMidPoint)]
	'''

	if (args.verbose):
		print >> sys.stdout, "- %s START   : populate lookup table with given resolution for chromosomes matching pattern %s" % (timeStamp(), args.chromPattern)

	fragmentsCount = 0
	fragmentsLookupTable = {}
	positionLookupTable = {}
	 
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
			positionLookupTable[tuple([chrom, int(0.5*(start+end))])] = fragmentsCount
			fragmentsLookupTable[fragmentsCount] = tuple([chrom, int(0.5*(start+end))])
			fragmentsCount += 1
	
	if (args.verbose):
		print >> sys.stdout, "- %s FINISH  : counted %d fragments" % (timeStamp(), fragmentsCount)

	return [ positionLookupTable, fragmentsLookupTable, fragmentsCount ]

def readFileIntoSparseMatrix(positionLookupTable, fragmentsCount):

	A = lil_matrix((len(args.contactsCountMatrices), fragmentsCount * fragmentsCount), dtype='f')
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
			ch1,mid1,ch2,mid2,contactCounts,pvalue,qvalue,bias1,bias2,biasCorrectedContactCount=line.split()

			try:
				if (float(qvalue) > args.threshold):
					continue
			except:
				# probably the header line
				continue

			if (args.metric =="pvalue"):
				metric=1.-(float(pvalue))

			elif (args.metric == "qvalue"):
				metric=1.-(float(qvalue))
			else:
				metric=numpy.arcsinh(float(biasCorrectedContactCount))
			
			# skip irrelevant entries
			if (metric == 0.):
				continue

			# skip trand counts if requested
			if (args.cis and ch1 != ch2):
				continue
			
			if (ch1 == ch2 and mid1 == mid2 and not args.keepdiagonal):
				continue

			fragment1 = tuple([ch1,int(mid1)])
			fragment2 = tuple([ch2,int(mid2)])


			if (not positionLookupTable.has_key(fragment1) or not positionLookupTable.has_key(fragment2)):
				if (args.veryverbose):
					print "[NOTE] fragment %s or %s not detected in map, is the correct resolution specified?" % (ch1,ch2)
				continue

			# keep this symmetric matrix as sparse as possible by just filling in the top triangle
			if (fragment1 < fragment2):
				A[c, positionLookupTable[fragment1] * fragmentsCount + positionLookupTable[fragment2]] = metric
			else:
				A[c, positionLookupTable[fragment2] * fragmentsCount +  positionLookupTable[fragment1]] =  metric

		c += 1 

		if (args.verbose):
			print >> sys.stdout, "- %s FINISH  : reading file" % (timeStamp())
			
	return (S,A.tocsr())

	# # TODO Normalize columns to 1
	# if (args.verbose):
	# 	print >> sys.stdout, "- %s START   : normalize to frequencies" % (timeStamp())

	# # print A
	# colsum = A.sum(axis=1)

	# A = A.tocoo()
	
	# W = lil_matrix((len(args.contactsCountMatrices), fragmentsCount * fragmentsCount), dtype='f')
	# for i,j,v in itertools.izip(A.row, A.col, A.data):
	# 	W[i,j] = float(v)/colsum[i]

	# if (args.verbose):
	# 	print >> sys.stdout, "- %s FINISH  : normalize to frequencies" % (timeStamp())

	# return (S,W.tocsr())

def performNMF(M, fragmentsLookupTable, fragmentsCount):

	if (args.verbose):
		print >> sys.stdout, "- %s START   : calculating NMF" % (timeStamp())

	t0 = time()	
	model = NMF(n_components=args.components, init='nndsvd', beta=10000.0, max_iter=1000, tol=5e-3, sparseness='components')
	model.fit(M)
	train_time = (time() - t0)
	components_ = model.components_

	# print >> sys.stdout, components_
	N = model.transform(M)

	if (args.verbose):
		print >> sys.stdout, "- %s FINISH  : calculating NMF" % (timeStamp())

	if (args.verbose):
		print >> sys.stdout, "- %s START   : mapping components" % (timeStamp())

	# convert components into locations
	for i in xrange(args.components):
		output = gzip.open(args.outdir+"/NMF_component_"+str(i)+".txt.gz", 'wb')
		if (args.verbose):
			print >> sys.stdout, "-            : processing component %d" % (i)

	 	try:

			for j in xrange(model.components_[i].shape[0]):
					# print >> sys.stdout, model.components_[i]
					# print >> sys.stdout, "Max value %f" (numpy.max(model.components_[i]))
				# if (model.components_[i][j] != 0):
					fragment1 = j / fragmentsCount
					fragment2 = j % fragmentsCount
					
					(chr1, midpoint1) = fragmentsLookupTable[fragment1]
					(chr2, midpoint2) = fragmentsLookupTable[fragment2]
					output.write("%s\t%i\t%s\t%i\t%f\n" % (chr1, midpoint1, chr2, midpoint2, model.components_[i][j]))

		finally:
				output.close()

	if (args.verbose):
		print >> sys.stdout, "- %s FINISH  : mapping components" % (timeStamp())

	return (N,model)

def plotVariance(S, N):

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

	f.write("pdf('"+args.outdir+'/'+args.imagename+".pdf', width=10, height=8)\n")
	f.write("g1\n")
	f.write("g2\n")
	f.write("dev.off()\n")

	f.close()

def plot_gallery(title, images, image_shape, n_col=2, n_row=1):
	if (args.verbose):
		print >> sys.stdout, "- %s START   : plotting figure %s" % (timeStamp(), title)

	plt.figure(figsize=(20. * n_col, 25.2 * n_row))
	plt.suptitle(title, size=16)
	for i, comp in enumerate(images):
		plt.subplot(n_row, n_col, i + 1)
		vmax = max(comp.max(), -comp.min())
		plt.imshow(comp.reshape(image_shape), cmap=plt.cm.gray, interpolation='nearest', vmin=-vmax, vmax=vmax)
		plt.xticks(())
		plt.yticks(())
	plt.subplots_adjust(0.01, 0.05, 0.99, 0.93, 0.04, 0.)

	if (args.verbose):
		print >> sys.stdout, "- %s FINISH  : plotting figure " % (timeStamp())


def main():
	''' main method 
	    parsing the parameters and options
	'''
	global args
	parser = argparse.ArgumentParser(description='Reads a list of significant interaction sparse matrix files and performs NMF')
	parser.add_argument('chromSizes', type=str, help='chomosome sizes')
	parser.add_argument('contactsCountMatrices', metavar="contactsCountMatrices", type=str, nargs='+', help='sparse interaction matrices')
	parser.add_argument("-o", '--outdir', dest='outdir', type=str, default="./",
						help='output location')
	parser.add_argument("-r", "--resolution", type=int, dest="resolution", default=1000000, 
						help="size of a fragment in bp if no genomeFragmentFile is given")
	parser.add_argument("-n", "--components", type=int, dest="components", default=2, 
						help="number of components used in NMF, default=2")
	parser.add_argument("-C", "--chrompattern", type=str, dest="chromPattern", default="", 
						help="pattern of chromosomes to filter for [default all]")
	parser.add_argument("-t", "--threshold", type=float, dest="threshold", default=1.0, 
						help="q-value threshold used to filter data [default 1.0]")
	parser.add_argument("--cis", action="store_true", help="consider cis interactions only [default all]")
	parser.add_argument("--keepdiagonal", action="store_true", help="also consider entries on diagonal [default remove diagnoal]")
	parser.add_argument("-p", "--plotGraphic", type=str, dest="imagename", default="", 
						help="create a 2D plot with the samples")

	parser.add_argument("-g", "--groups", type=str, dest="groups", default="", 
						help="group list of count matrices via comma-separated list, e.g. 1,1,3,4")
	parser.add_argument("-l", "--labels", type=str, dest="labels", default="", 
						help="text labels for experiments used in plots, supplied via comma-separated list, e.g. PrEC_I,PrEC_II,LNCaP_I,LNCaP_II")
	parser.add_argument("-c", "--cutter", type=str, dest="cutter", default="", 
						help="restriction enzymes used in matrices, supplied via comma-separated list, e.g. HindIII,NcoI,HindIII,NcoI")
	parser.add_argument("-m", "--metric", type=str, dest="metric", default="biasCorrectedContactCount", 
						help="on of the following: pvalue, qvalue, biasCorrectedContactCount")

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

	[ positionLookupTable, fragmentsLookupTable, fragmentsCount ] = createFragmentResolution()


	(S, A) = readFileIntoSparseMatrix(positionLookupTable, fragmentsCount)

	image_shape=(fragmentsCount,fragmentsCount)
	
	plot_gallery("Default", A.todense()[:len(args.contactsCountMatrices)], image_shape,  n_col=len(args.contactsCountMatrices), n_row=1)
	plt.savefig(pp, format='pdf')

	(N, model) = performNMF(A, fragmentsLookupTable, fragmentsCount)

	plot_gallery("First chromosomes", model.components_[:args.components], image_shape,  n_col=args.components, n_row=1)
	
	plt.savefig(pp, format='pdf')
	pp.close()

	# if (args.imagename != ""):
	# 	plotVariance(S, N, explainedRatio, totalExplained)

if __name__ == "__main__":
		main()
