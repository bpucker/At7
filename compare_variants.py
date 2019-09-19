### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###


__usage__ = """
					python compare_variants.py
					--vcf1 <FULL_PATH_TO_VCF_FILE1>
					--vcf2 <FULL_PATH_TO_VCF_FILE2>
					--out <FULL_PATH_TO_OUTPUT_FILE>
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import numpy as np
from scipy import stats
import sys


def load_variants( vcf ):
	"""! @brief load all variants from given VCF """
	
	variants = []
	indel_counter = 0
	with open( vcf, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				variants.append( { 'ref': parts[3], 'alt': parts[4] } )
				if len( parts[3] ) != len( parts[4] ):
					indel_counter += 1
			line = f.readline()
	print "number of InDels: " + str( indel_counter )
	print "number of SNPs and MNPs: " + str( len( variants ) - indel_counter )
	return variants


def calculate_matrix( variants ):
	"""! @brief calculate matrix """
	
	SNP_matrix = {}
	for variant in variants:
		if len( variant['ref'] ) == 1:
			if len( variant['alt'] ) == 1:
				try:
					SNP_matrix[ variant['ref']+variant['alt'] ] += 1
				except KeyError:
					SNP_matrix.update( { variant['ref']+variant['alt']: 1 } )
	return SNP_matrix


def normalize_matrix( matrix ):
	"""! @brief normalize matrix """
	
	norm_mat = {}
	total = float( sum( matrix.values() ) )
	for key in matrix.keys():
		norm_mat.update( { key: 100.0 * matrix[ key ] / total } )
	return norm_mat


def check_ins_vs_del( variants ):
	"""! @brief quantify insertions and deletions among provided variants """
	
	insertion = 0
	deletion = 0
	for variant in variants:
		if len( variant['ref'] ) < len( variant['alt'] ):
			insertion += 1
		elif len( variant['ref'] ) > len( variant['alt'] ):
			deletion += 1
	print "number of insertions: " + str( insertion )
	print "number of deletions: " + str( deletion )


def main( arguments ):
	"""! @brief compare variant data sets """
	
	vcf1 = arguments[ arguments.index('--vcf1')+1 ]
	vcf2 = arguments[ arguments.index('--vcf2')+1 ]
	out_file = arguments[ arguments.index('--out')+1 ]
	
	print vcf1
	variants1 = load_variants( vcf1 )

	print vcf2
	variants2 = load_variants( vcf2 )
	
	print vcf1
	check_ins_vs_del( variants1 )
	
	print vcf2
	check_ins_vs_del( variants2 )
	
	matrix1 = calculate_matrix( variants1 )
	matrix2 = calculate_matrix( variants2 )

	data1 = []
	data2 = []
	for key in matrix1.keys():
		data1.append( matrix1[ key ] )
		data2.append( matrix2[ key ] )

	print stats.chisquare( data1, data2 )


	norm_matrix1 = normalize_matrix( matrix1 )
	norm_matrix2 = normalize_matrix( matrix2 )

	diff = []
	diff_matrix = {}
	for key in matrix1.keys():
		x = norm_matrix1[ key ] - norm_matrix2[ key ]
		diff.append( x )
		diff_matrix.update( { key: x } )

	with open( out_file, "w" ) as out:
		out.write( "\t".join( [ "", "A", "C", "G", "T" ] ) + '\n' )
		for nt1 in "ACGT":
			new_line = [ nt1 ]
			for nt2 in "ACGT":
				try:
					new_line.append( diff_matrix[ nt1+nt2 ] )
				except KeyError:
					new_line.append( 0 )
			out.write( "\t".join( map( str, new_line ) ) + '\n' )


if '--vcf1' in sys.argv and '--vcf2' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
