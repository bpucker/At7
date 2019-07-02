### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###


__usage__ = """
					python compare_ins_to_del.py
					--vcf <FULL_PATH_TO_VCF_FILE>
					--out <FULL_PATH_TO_OUTPUT_FOLDER>
					"""

import matplotlib.pyplot as plt
import numpy as np
import sys, os

# --- end of imports --- #

def main( arguments ):
	"""! @brief run everything """

	vcf = arguments[ arguments.index( '--vcf' )+1 ]
	prefix = arguments[ arguments.index( '--out' )+1 ]
	
	if prefix[-1] != "/":
		prefix += "/"
	if not os.path.exists( prefix ):
		os.makedirs( prefix )
	
	sample = []
	insertions = 0
	deletions = 0

	with open( vcf, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				sample.append( len( parts[3] ) - len( parts[4] ) )
				if len( parts[3] ) - len( parts[4] ) > 0:
					deletions += 1
				elif len( parts[3] ) - len( parts[4] ) < 0:
					insertions += 1
			line = f.readline()


	fig, ax = plt.subplots()
	ax.hist( sample, bins=1000, color="lime", alpha=0.5, label="heterozygous variants" )
	ax.set_xlim( -20, 20 )
	ax.set_xlabel( "length difference between reference and alterantive allele [bp]" )
	ax.set_ylabel( "number of variants" )
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	fig.savefig( prefix + "allele_length_differences.png", dpi=300 )


	fig, ax = plt.subplots( figsize=(2,4) )
	ax.bar( [ 1, 2 ], [ insertions, deletions ], color= ["magenta", "lime"] )
	ax.set_ylabel( "number of variants" )
	ax.xaxis.set_ticks( [ 1, 2 ] )
	ax.set_xticklabels( [ "insertions", "deletions" ] )
	ax.set_xlim( 0.5, 2.5 )
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	plt.subplots_adjust( left=0.3, top=0.99, right=0.99, bottom=0.1 )
	fig.savefig( prefix + "inversions_deletion_numbers.png", dpi=300 )

	plt.close( "all" )


if '--vcf' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
