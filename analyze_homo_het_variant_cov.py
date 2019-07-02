### only documentation ###

import matplotlib.pyplot as plt

# --- end of imports --- #

def load_values( vcf ):
	"""! @brief load coverage of homozygous and heterozygous variants """
	
	homo = []
	hetero = []
	with open( vcf, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if not "," in parts[4]:
				try:
					reads = parts[-1].split(':')[1].split(',')
					if len( reads ) == 2:
						reads1 = float( reads[0] )
						reads2 = float( reads[1] )
						if min( [ reads1, reads2 ] ) / (reads1 + reads2 ) < 0.05:	#homozygous
							homo.append( min( [ 1000, sum( map( int, parts[-1].split(':')[1].split(',') ) ) ] ) )
						else:	#heterozygous
							hetero.append( min( [ 1000, sum( map( int, parts[-1].split(':')[1].split(',') ) ) ] ) )
				except IndexError:
					print line
				# if parts[-1][:3] == "0/1":
					# hetero.append( min( [ 1000, sum( map( int, parts[-1].split(':')[1].split(',') ) ) ] ) )
				# elif parts[-1][:3] == "1/1":
					# homo.append( min( [ 1000, sum( map( int, parts[-1].split(':')[1].split(',') ) ) ] ) )
			line = f.readline()
	return homo, hetero



vcf_file = "At7_variants_vs_Col.vcf"
fig_file = "HOMO_HETERO_VARIANT_COV.png"

homo, hetero = load_values( vcf_file )

fig, ax = plt.subplots()
ax.hist( hetero, bins=1000, color="lime", label="variants with diviating frequencies" )
ax.hist( homo, bins=1000, color="magenta", label="homozygous variants" )
ax.set_xlim( 0, 200 )
ax.set_xlabel( "sequencing coverage depth" )
ax.set_ylabel( "number of variant positions" )
ax.legend( prop={ 'size':8 } )
fig.savefig( fig_file, dpi=300 )

plt.close( "all" )
