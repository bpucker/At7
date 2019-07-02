### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python generate_cov_hist.py
					--in <COVERAGE_INPUT_FILE>
					--out <FIG_OUTPUT_FILE>
					"""

import matplotlib.pyplot as plt
import sys

# --- end of imports --- #


def load_cov( cov_file, cutoff ):
	"""! @brief load all information from coverage file """
	
	cov = []
	with open( cov_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			cov.append( min( [ float( parts[-1] ), cutoff ] ) )
			line = f.readline()
	return cov


def main( arguments ):
	
	cov_file = arguments[ arguments.index( '--in' ) + 1 ]
	fig_file = arguments[ arguments.index( '--out' ) + 1 ]
	
	if '--cut' in arguments:
		cutoff = float( arguments[ arguments.index( '--cut' ) + 1 ] )
	else:
		cutoff = 10000.0
	
	cov = load_cov( cov_file, cutoff )
	
	fig, ax = plt.subplots()
	ax.hist( cov, bins=1000, color="lime" )
	fig.savefig( fig_file, dpi=300 )
	
	plt.close( "all" )


if '--in' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
