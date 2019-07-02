### only documentation ###

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os

# --- end of imports --- #


def load_variants_from_vcf( vcf_file ):
	"""! @brief loads the variant informaiton from a SnpEff output VCF file """
	
	snps_per_chr = [ [], [], [], [], [] ]
	indels_per_chr = [ [], [], [], [], [] ]
	
	with open( vcf_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				try:
					parts = line.strip().split('\t')
					reads = parts[-1].split(':')[1].split(',')
					if len( reads ) == 2:
						reads1 = float( reads[0] )
						reads2 = float( reads[1] )
						if min( [ reads1, reads2 ] ) / (reads1 + reads2 ) < 0.05:	#homozygous
							snps_per_chr[ int( parts[0][-1] ) - 1 ].append( int( parts[1] ) )
						else:	#heterozygous
							indels_per_chr[ int( parts[0][-1] ) - 1 ].append( int( parts[1] ) )
					# if parts[-1][:3] == "1/1":	#homozygous variants
						# snps_per_chr[ int( parts[0][-1] ) - 1 ].append( int( parts[1] ) )
					# elif parts[-1][:3] == "0/1":	#heterozygous variants
						# indels_per_chr[ int( parts[0][-1] ) - 1 ].append( int( parts[1] ) )					
				except:
					pass	#print line
			line = f.readline()
	
	return snps_per_chr, indels_per_chr


def generate_binned_values( lower_lim, upper_lim, chr_length, snps_per_chr, indels_per_chr, resolution ):
	"""! @brief group variants into bins """
	
	snp_data = []
	indel_data = []
	while True:
		if upper_lim >= chr_length:
			break
		else:
			snp_tmp = []
			indel_tmp = []
			for SNP in snps_per_chr:
				if SNP <= upper_lim and SNP > lower_lim:
					snp_tmp.append( 'X' )
			for indel in indels_per_chr:
				if indel <= upper_lim and indel > lower_lim:
					indel_tmp.append( 'X' )
			snp_data.append( len( snp_tmp ) )
			indel_data.append( len( indel_tmp ) )
		upper_lim += resolution
		lower_lim += resolution
	return max( snp_data + indel_data ), snp_data, indel_data


def construct_plot( snps_per_chr, indels_per_chr, chr_lengths, result_file, result_table, resolution ):
	"""! @brief construct variant over Col-0 genome distribution plot """
	
	fig, ax = plt.subplots()
		
	scale = 1
	snp_data = []
	indel_data = []
	for idx, chr_length in enumerate( chr_lengths ):
		max_val, snp_temp, indel_temp = generate_binned_values( 0, resolution+0, chr_length, snps_per_chr[ idx ], indels_per_chr[ idx ], resolution )
		snp_data.append( snp_temp )
		indel_data.append( indel_temp )
		scale = max( [ scale, max_val ] )
	
	scale = float( scale )
	y_max = len( chr_lengths )
	
	with open( result_table, "w" ) as out:
		for idx, chr_length in enumerate( chr_lengths ):
			chr_name = 'Chr' + str( idx + 1 )
			y = y_max-( idx*1.2 )
			x = resolution / 1000000.0
			
			ax.text( ( chr_length/ 1000000.0 ), y+0.3, chr_name, ha="right" )
			
			# --- plotting SNP and InDel distribution --- #
			for i, snps in enumerate( snp_data[ idx ] ):
				indels = indel_data[ idx ][ i ]
				
				ax.plot( [ x*i+0.5*x, x*i+0.5*x ], [ y, y+ ( snps / scale ) ], "-", color="magenta" )
				ax.plot( [ x*i+0.5*x, x*i+0.5*x ], [ y, y+ ( indels / scale ) ], "-", color="lime" )
			
			ax.plot( [ 0, 0 ], [ y, y+1 ], color="black" )
			ax.text( 0, y+1, str( int( scale ) ), ha="right", fontsize=5 )
			ax.text( 0, y+0.5, str( int( scale / 2 ) ), ha="right", fontsize=5 )
			ax.text( 0, y, "0", ha="right", fontsize=5 )
			
			# --- writing data into output table --- #
			out.write( 'Chr' + str( idx+1 ) + "homozygous variants:\t" + '\t'.join( map( str, snp_data ) ) + '\n' )
			out.write( 'Chr' + str( idx+1 ) + "heterozygous variants:\t" + '\t'.join( map( str, indel_data ) ) + '\n' )
	
	ax.set_xlabel( "genomic position [ Mbp ]" )
	ax.set_ylabel( "number of variants per interval" )	#"counts", "SNPs (black), InDel(red)"
	
	ax.set_xlim( 0, max( chr_lengths ) / 1000000.0 )	
	
	ax.legend( handles=[ mpatches.Patch(color='magenta', label='homozygous variants' ), mpatches.Patch(color='lime', label='variants with diviating frequency') ], prop={'size':8} )
	
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.spines['left'].set_visible(False)
	
	ax.get_yaxis().set_ticks([])
	ax.yaxis.labelpad = 15
	
	
	plt.subplots_adjust( left=0.07, right=0.98, top=0.99, bottom=0.1 )
	fig.savefig( result_file, dpi=600 )
	
	plt.close('all')


if __name__ == '__main__':
	vcf_file = "At7_vs_Col0_variants.vcf"
	
	output_dir = "./homo_hetero_variants/"
	
	if not output_dir[-1] == "/":
		output_dir += "/"
	
	if not os.path.exists( output_dir ):
		os.makedirs( output_dir )
	
	result_file = output_dir + "genome_wide_small_variants.png"
	result_table = output_dir + "genome_wide_small_variants.txt"
		
	snps_per_chr, indels_per_chr = load_variants_from_vcf( vcf_file )
	
	print "number of homozygous variants: " + str( len( [ x for each in snps_per_chr for x in each ] ) )
	print "number of variants with deviating coverage: " + str( len( [ x for each in indels_per_chr for x in each ]) )
	
	chr_lengths = [ 30427671, 19698289, 23459830, 18585056, 26975502 ]
	
	construct_plot( snps_per_chr, indels_per_chr, chr_lengths, result_file, result_table, resolution = 100000 )
	
	print "all done!"
