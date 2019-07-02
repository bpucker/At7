### only documentation ###

import matplotlib.pyplot as plt

# --- end of imports --- #

vcf_file = "At7_variants_vs_Col.vcf"
fig_file = "variant_frequencies.png"
cov_fig_file = "variant_coverages.png"

values = []
coverage = []


with open( vcf_file, "r" ) as f:
	line = f.readline()
	while line:
		if line[0] != '#':
			parts = line.strip().split('\t')
			reads = parts[-1].split(':')[1].split(',')
			if len( reads ) == 2:
				x = float( reads[0] )
				y = float( reads[1] )
				if x / y < 0.05:	#homozygous
					coverage.append( x+y )
				elif x / y > 0.95:	#homozygous
					coverage.append( x+y )
				else:	#heterozygous
					coverage.append( x+y )
					if x+y > 20:
						values.append( x / ( x+y ) )
		line = f.readline()


fig, ax = plt.subplots()
ax.hist( coverage, bins=10000, color="lime" )
ax.set_xlim( 0, 300 )
ax.set_xlabel( "sequencing coverage" )
ax.set_ylabel( "number of variants" )
fig.savefig( cov_fig_file, dpi=300 )


fig, ax = plt.subplots()
ax.hist( values, bins=100, color="lime" )
ax.set_xlabel( "allele frequency" )
ax.set_ylabel( "number of variants" )
fig.savefig( fig_file, dpi=300 )


