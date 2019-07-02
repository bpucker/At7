### only documentation ###

import matplotlib.pyplot as plt
import numpy as np
import re

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


illumina_cov_file = "At7_illumina_vs_Col0.cov"
ont_cov_file = "At7_ONT_vs_Col0.cov"
col_cov_file = "GK_vs_Col0.cov"
nd1_file = "Nd1_reads_vs_Col0.sorted.cov"

output_dir = "./check_organell_cov/"

hist_file = output_dir + "chloro_mito_cov_hist.png"


illumina_cov = load_cov( illumina_cov_file )
ont_cov = load_cov( ont_cov_file )
col_cov = load_cov( col_cov_file )
nd1 = load_cov( nd1_file )

at7_cov = {}
for key in illumina_cov.keys():
	vals = illumina_cov[ key ]
	new_vals = []
	for idx, each in enumerate( vals ):
		new_vals.append( each + ont_cov[ key ][ idx ] )
	at7_cov.update( { key: new_vals } )


# --- generate histogram --- #
fig, ax = plt.subplots( figsize=( 8, 4 ) )

values = [ 	at7_cov['chloroplast'],  at7_cov['mitochondria'], col_cov['chloroplast'],  col_cov['mitochondria'],
					nd1['chloroplast'],  nd1['mitochondria']
				]

labels = [ "At7 plastome", "At7 chondrome", "Col-0 plastome", "Col-0 chondrome", "Nd-1 plastome", "Nd-1 chondrome" ]

ax.boxplot( values,
					notch=True,
					labels=labels,
					autorange = True,
					showmeans = True
					)

plt.subplots_adjust( left=0.1, right=0.999, top=0.99, bottom=0.1)
ax.set_ylabel( "sequencing coverage depth" )
fig.savefig( hist_file, dpi=300 )
plt.close( "all" )


for idx, each in enumerate( values ):
	print labels[ idx ] + ": " + str( np.median( each ) )
