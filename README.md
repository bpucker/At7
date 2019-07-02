# At7 genome sequencing and analysis

Python scripts related to the At7 (Arabidopsis thaliana suspension cell culture) sequencing project. Some scripts are committed to customized analysis and are just included for documentation purposes. Only the usage of scripts re-useable beyond the At7 data set analyzed here will be described below.


python SnpEff_result_parser.py \
--in <FULL_PATH_TO_INPUT_VCF> \
--out <FULL_PATH_TO_OUTPUT_TEXT_FILE> \
--anno <ANNOTATION_FILE>


python get_ZCRs_single_cov.py \
--cov <FULL_PATH_TO_REFERENCE_COVERAGE_FILE> \
--out <FULL_PATH_TO_OUTPUT_FILE>


python generate_cov_hist.py \
--in <COVERAGE_INPUT_FILE> \
--out <FIG_OUTPUT_FILE>


python FASTQ_stats.py \
--in_file <FULL_PATH_TO_FASTQ_FILE> |	--in_dir <FULL_PATH_TO_DIRECTORY>


python cov_plot.py \
--in <FULL_PATH_TO_COVERAGE_FILE> \
--out <FULL_PATH_TO_OUTPUT_FILE> \
optional: \
--res <RESOLUTION, WINDOW_SIZE_FOR_COVERAGE_CALCULATION> \
--sat <SATURATION, CUTOFF_FOR_MAX_COVERAGE_VALUE> \
--cov <AVERAGE_COVERAGE>


python compare_ins_to_del.py \
--vcf <FULL_PATH_TO_VCF_FILE> \
--out <FULL_PATH_TO_OUTPUT_FOLDER>


python analyze_indel_len_CDS.py \
--vcf <FULL_PATH_TO_VCF_FILE> \
--gff <FULL_PATH_TO_GFF3_FILE> \
--out <FULL_PATH_TO_OUTPUT_FOLDER>

