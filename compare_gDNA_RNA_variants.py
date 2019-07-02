### only for documentation ###

import re
import matplotlib.pyplot as plt
from scipy import stats

# --- end of imports --- #

def load_genomic_variants( gDNA_variant_file ):
	"""! @brief load genomic variants """
	
	gDNA_variants = {}
	with open( gDNA_variant_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			gDNA_variants.update( { parts[0] + '_%_' + ( parts[1].zfill( 8 ) ): parts[4] } )
			line = f.readline()
	return gDNA_variants


def load_RNA_variants( RNA_variant_file ):
	"""! @brief load RNA variants with coverage information """
	
	RNA_variants = {}
	with open( RNA_variant_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				if parts[-1][:3] == "0/1":
					numbers = map( int, re.findall( "DP4=\d+,\d+,\d+,\d+", parts[-3] )[0][4:].split(',') )
					x = sum( numbers[:2] )
					y = sum( numbers[3:] )
					RNA_variants.update( { parts[0] + '_%_' + ( parts[1].zfill( 8 ) ): { 'ref': x, 'alt': y } } )
			line = f.readline()
	return RNA_variants


def load_variants_of_interest( gDNA_variants, genomic_vcf ):
	"""! @brief load genomics variants of interest """
	
	ref_cov = []
	alt_cov = []
	
	with open( genomic_vcf, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				try:
					gDNA_variants[ parts[0] + '_%_' + ( parts[1].zfill( 8 ) )  ]
					try:
						ref, alt = map( int, parts[-1].split(':')[1].split(',') )
						ref_cov.append( ref )
						alt_cov.append( min( [ alt, 400 ] ) )
					except ValueError:
						print line
				except KeyError:
					pass
			line = f.readline()
	return ref_cov, alt_cov


gDNA_variant_file = "At7_vs_Col0_snpeff.results"
RNA_variant_file = "At7_RNA_seq_reads_vs_Col0.vcf"

genomic_vcf = "At7_DNA_seq_reads_vs_Col0.vcf"

fig_file = "ASE_at_high_impact_variants1.png"
fig_file2 = "ASE_at_high_impact_variants2.png"
fig_file3 = "ASE_at_high_impact_variants3.png"


gDNA_variants = load_genomic_variants( gDNA_variant_file )
RNA_variants = load_RNA_variants( RNA_variant_file )

print "number of affected genes: " + str( len( list( set( gDNA_variants.values() ) ) ) )

matches = 0
fails = 0

ref_values = []
alt_values = []
for key in gDNA_variants.keys():
	try:
		ref_values.append( RNA_variants[ key ]['ref'] )
		alt_values.append( RNA_variants[ key ]['alt'] )
		matches += 1
	except KeyError:
		fails += 1

print "number of high impact variants present RNA-Seq data: " + str( matches )
print "number of high impact variants NOT present RNA-Seq data: " + str( fails )


fig, ax = plt.subplots()
ax.boxplot( [ ref_values, alt_values ], labels=["reference", "alternative"] )
ax.set_xlabel( "alleles at high impact variant loci" )
ax.set_ylabel( "TPMs" )
fig.savefig( fig_file, dpi=300 )


fig, ax = plt.subplots()
ax.scatter( ref_values, alt_values, color="lime", s=3 )
ax.set_xlabel( "TPM of reference allele at high impact variant loci" )
ax.set_ylabel( "TPM of alternative allele at high impact variant loci" )
fig.savefig( fig_file2, dpi=300 )


plt.close("all")

u, p = stats.mannwhitneyu( ref_values, alt_values )
print "Mann Whitney U test: U=" + str( u ) + "\tp=" + str( p )

############################################################

ref_cov, alt_cov = load_variants_of_interest( gDNA_variants, genomic_vcf )


fig, ax = plt.subplots()
ax.boxplot( [ ref_cov, alt_cov ], labels=["reference", "alternative"] )
ax.set_xlabel( "alleles at high impact variant loci" )
ax.set_ylabel( "number of supporting DNA reads" )
fig.savefig( fig_file3, dpi=300 )


u, p = stats.mannwhitneyu( ref_cov, alt_cov )
print "Mann Whitney U test: U=" + str( u ) + "\tp=" + str( p )
