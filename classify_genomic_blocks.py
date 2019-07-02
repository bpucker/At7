### only documentation ###

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


at7_cov_file = "At7_illumina_vs_TAIR10.cov"
fig_file = "RESULTS.png"
output_file = "GENOMIC_BLOCK_CLASSIFICATION.txt"

block_size = 5000
cutoff = 200

at7_cov = load_cov( at7_cov_file )

values = []
with open( output_file, "w" ) as out:
	for key in sorted( at7_cov.keys() ):
		a_cov = at7_cov[ key ]
		start = 0
		end = 0 + block_size
		while end < len( a_cov ):
			val = min( [ np.mean( a_cov[ start:end ] ), cutoff ] ) 
			values.append( val )
			if val < 20:
				out.write( "\t".join( map( str, [ key, start, end, "missing" ] ) ) + '\n' )
			elif val <= 40:
				out.write( "\t".join( map( str, [ key, start, end, "A" ] ) ) + '\n' )
			elif val <= 65:
				out.write( "\t".join( map( str, [ key, start, end, "B" ] ) ) + '\n' )
			elif val <= 100:
				out.write( "\t".join( map( str, [ key, start, end, "C" ] ) ) + '\n' )
			elif val <= 130:
				out.write( "\t".join( map( str, [ key, start, end, "D" ] ) ) + '\n' )
			elif val <= 160:
				out.write( "\t".join( map( str, [ key, start, end, "E" ] ) ) + '\n' )
			else:
				out.write( "\t".join( map( str, [ key, start, end, "collapsed repeat" ] ) ) + '\n' )
			start += block_size
			end += block_size


fig, ax = plt.subplots()
ax.hist( values, bins=1000, color="lime" )
ax.set_xlim( 0, cutoff )
ax.set_xlabel( "sequencing coverage depth" )
ax.set_ylabel( "number of genomic blocks" )
fig.savefig( fig_file, dpi=300 )
plt.close( "all" )
