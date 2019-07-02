### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

### WARNING: optimized for Arabidopsis thaliana Col-0 ###

__usage__ = """
					python cov_plot.py
					--in <FULL_PATH_TO_COVERAGE_FILE>
					--out <FULL_PATH_TO_OUTPUT_FILE>
					
					--res <RESOLUTION, WINDOW_SIZE_FOR_COVERAGE_CALCULATION>
					--sat <SATURATION, CUTOFF_FOR_MAX_COVERAGE_VALUE>
					--cov <AVERAGE_COVERAGE>
					"""

import sys, os
import matplotlib.pyplot as plt
import numpy as np

# --- end of imports --- #


def load_cov( cov_file ):
	"""! @brief load all information from coverage file """
	
	cov = {}
	with open( cov_file, "r" ) as f:
		line = f.readline()
		header = line.split('\t')[0]
		tmp = []
		while line:
			parts = line.strip().split('\t')
			if parts[0] != header:
				cov.update( { header: tmp } )
				header = parts[0]
				tmp = []
			tmp.append( float( parts[-1] ) )
			line = f.readline()
		cov.update( { header: tmp } )
	return cov


def generate_plot( cov, out_file, resolution, saturation, coverage ):
	"""! @brief generate figure """
	
	fig, ax = plt.subplots( figsize=( 20, 5 ) )
	
	ymax = 5	#len( cov.keys() )+1
	max_value = 0
	collected_values = {}
	
	# --- generate list for plotting --- #
	for idx, key in enumerate( sorted( cov.keys() ) ):
		y = ymax-idx-1
		x = []
		blocks = [ cov[ key ] [ i : i + resolution ] for i in xrange( 0, len( cov[ key ] ), resolution ) ]
		for block in blocks:
			x.append( np.mean( block ) )
		max_value = max( [ max_value, max( x ) ] )
		collected_values.update( { key: x } )
	
	# --- plot values --- #
	max_value = float( min( [ saturation, max_value ] ) )
	for idx, key in enumerate( sorted( cov.keys() )[:5] ):
		y = ymax - ( idx*1.3 )
		x = []
		for each in collected_values[ key ]:
			x.append( y + min( [ 1, ( each / max_value ) ] ) )
		ax.plot( x, linestyle=":", color="green" )
		
		ax.plot( [ 0, len( x ) ], [ y+( coverage / saturation ), y+( coverage / saturation ) ], color="red" , linewidth=1)
		
		ax.plot( [ 0, 0 ], [ y, y+1 ], color="black", linewidth=1, markersize=1 )
		ax.text( 0, y+1, str( int( max_value ) ), ha="right", fontsize=5 )
		ax.text( 0, y+0.5, str( int( max_value / 2 ) ), ha="right", fontsize=5 )
		ax.text( 0, y, "0", ha="right", fontsize=5 )
		
	ax.set_xlabel( "position on chromosome [ Mbp ]" )
	ax.set_ylabel( "coverage" )
	
	ax.set_xlim( 0, 30500 )
	
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.spines['left'].set_visible(False)
	ax.get_yaxis().set_ticks([])
	ax.yaxis.labelpad = 5
	
	ax.set_xticklabels( [ "0", "5", "10", "15", "20", "25", "30" ] )
	
	plt.subplots_adjust( left=0.015, right=0.999, top=0.99, bottom=0.1 )
	
	fig.savefig( out_file, dpi=300 )


def generate_hist( cov_values, outputfile, saturation ):
	"""! @brief generate coverage histogram """
	
	values = []
	for each in cov_values:
		if each > saturation:
			values.append( saturation )
		else:
			values.append( each )
	
	
	fig, ax = plt.subplots()
	
	ax.hist( values, bins=300 )
	ax.set_xlim( 0, 200 )
	
	fig.savefig( outputfile, dpi=300 )


def main( arguments ):
	"""! @brief runs everything """
	
	cov_file = arguments[ arguments.index( '--in' ) + 1 ]
	out_file = arguments[ arguments.index( '--out' ) + 1 ]
	
	if '--res' in arguments:
		resolution = int( arguments[ arguments.index( '--res' ) + 1 ] )
	else:
		resolution = 1000
	
	if '--sat' in arguments:
		saturation = int( arguments[ arguments.index( '--sat' ) + 1 ] )
	else:
		saturation = 300
	
	if '--cov' in arguments:
		coverage = float( arguments[ arguments.index( '--cov' ) + 1 ] )
	else:
		coverage = 100.0
	
	cov = load_cov( cov_file )
	
	# --- generate coverage histograms per chromosome --- #
	for key in cov.keys():
		outputfile = out_file + key + ".png"
		generate_hist( cov[ key ], outputfile, saturation )
	
	
	# --- generate per chromosome position coveage plot --- #
	generate_plot( cov, out_file, resolution, saturation, coverage )
	
	


if '--in' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
