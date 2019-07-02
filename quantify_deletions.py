### only documentation ###

import matplotlib.pyplot as plt
import numpy as np
import re

# --- end of imports --- #

def load_cov( coverage_file, cutoff ):
	"""! @brief load all information from coverage file """
	
	cov = {}
	with open( coverage_file, "r" ) as f:
		line = f.readline()
		header = line.split('\t')[0]
		tmp = []
		while line:
			parts = line.strip().split('\t')
			if parts[0] != header:
				cov.update( { header: tmp } )
				header = parts[0]
				tmp = []
			tmp.append( min( [ float( parts[-1] ), cutoff ] ) )
			line = f.readline()
		cov.update( { header: tmp } )
	return cov


def get_coverage_and_windows( coverage, window, step ):
	"""! @brief calculate coverage per window and window center """
	
	y_values = []
	for key in sorted( coverage.keys() )[:5]:
		cov = coverage[ key ]
		start = 0
		end = 0 + step
		while end < len( cov ):
			y_values.append( sum( cov[ start:end ] ) / float( window ) )
			start += step
			end += step
	return y_values


def load_BUSCO_positions( busco_position_file ):
	"""! @brief load single copy and complete BUSCO positions """
	
	busco_pos = []
	with open( busco_position_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if parts[1] == "Complete":
				busco_pos.append( { 'chr': parts[2], 'start': int( parts[3] ),  'end': int( parts[4] ) } )
			line = f.readline()
	return busco_pos


def calculate_cov_distr( busco_pos, coverage ):
	"""! @brief calculate average coverage and define lower and upper cutoff based on distribution """
	
	values = []
	for busco in busco_pos:
		values.append( np.mean( coverage[ busco['chr'] ][ busco['start']: busco['end'] ] ) )	#calculate average per BUSCO
	avg_coverage = np.mean( values )
	low_cut = avg_coverage - ( 2 * np.std( values ) )
	up_cut = avg_coverage + ( 2 * np.std( values ) )
	return avg_coverage, low_cut, up_cut


def generate_figure( avg_cov, low_cut, up_cut , windows, outputfile ):
	"""! @brief generate histogram and print result value """
	
	fig, ax = plt.subplots()
	
	low_cut = 50
	
	ax.hist( windows, bins=1000, color="lime" )
	ax.plot( [ low_cut, low_cut ], [ 0, 500 ], color="magenta" )
	ax.plot( [ avg_cov, avg_cov ], [ 0, 500 ], color="black" )
	
	ax.set_xlabel( "sequencing coverage depth" )
	ax.set_ylabel( "number of genomic blocks" )
	
	ax.set_xlim( 0, 350 )
	
	fig.savefig( outputfile, dpi=300 )
	plt.close( "all" )
	
	counter = 0
	for each in windows:
		if each < low_cut:
			counter += 1
	print "number of regions deleted in one haplotype: " + str( counter )


illumina_cov_file = "Nd1_Illumina_reads_vs_Col0.sorted.cov"
ont_cov_file = "Nd1_ONT_reads_vs_Col0.sorted.cov"
busco_position_file = "full_table_BUSCO_on_Col0.tsv"

output_dir = "./deletion_quantification/"

out_file = output_dir + "deletion_quantification.png"


window = 1000
step = 1000
cutoff = 175.0	#maximal coverage

illumina_cov = load_cov( illumina_cov_file, cutoff )
ont_cov = load_cov( ont_cov_file, cutoff )

busco_pos = load_BUSCO_positions( busco_position_file )

at7_cov = {}
for key in illumina_cov.keys():
	vals = illumina_cov[ key ]
	new_vals = []
	for idx, each in enumerate( vals ):
		new_vals.append( each + ont_cov[ key ][ idx ] )
	at7_cov.update( { key: new_vals } )

avg_cov, low_cut, up_cut = calculate_cov_distr( busco_pos, at7_cov )


at7_windows = get_coverage_and_windows( at7_cov, window, step )


generate_figure( avg_cov, low_cut, up_cut , at7_windows, out_file )
