#!/bin/bash


# USAGE: mapping2theMHCRef.sh sample_name project_dir sample_A_1.fastq sample_A_2.fastq  > run_pipeline.sh
# USAGE: pipeline_tophat2_mapping.sh  sample_name sample_A_1.fastq,sample_B_1.fastq sample_A_2.fastq,sample_B_2.fastq  > run_pipeline.sh
#        qsub run_pipeline.sh



echo "echo Start analysis at: \`date\`"


###################################################
#
#  SETUP PATHS
#
###################################################


# Installed tools and references

SEQTK=""
SAMTOOLS=""
BEDTOOLS=""
BAMTOOLS=""
TOPHAT2=""

GTF=""
BOWTIE2_GENOME_INDEX=""




###################################################
#
# TOPHAT2 PARAMETERS
#
###################################################

nthreads=4			# default: 1 
mate_inner_dist=150		# default: 50
mate_std_dev=100		# default: 20
read_mismatches=8		# default: 2
read_edit_dist=10		# default: 2
b2_N=1				# default: 0
library_type="fr-unstranded"	# default is unstranded (fr-unstranded)

read_gap_length=5		# default: 2
segment_mismatches=3		# default: 2


###################################################
#
# PROJECT DIRECTORY AND NAME OF DATA
#
###################################################

PROJECT_DIR="./"	# you can define your project directory
SAMPLE_NAME=$1



#OUTPUT_DIR=$PROJECT_DIR"/mapping/2.mapping/"
OUTPUT_DIR=$PROJECT_DIR

if [[ ! -e $OUTPUT_DIR ]]; then
                mkdir -p $OUTPUT_DIR
fi

SAMPLE_NAME=$1




###################################################
#
# INPUT DATA: fastq files
#
###################################################

INPUT_READ1=$2
INPUT_READ2=$3


###################################################
#
# TOPHAT2 MAPPING
#
###################################################

TOPHAT2_MAPPING(){

	echo "echo ========================================="
	echo "echo TOPHAT2 MAPPING: \`date\`"

	DIR_SAMPLE_NAME=$OUTPUT_DIR$SAMPLE_NAME
	if [[ ! -e $DIR_SAMPLE_NAME ]]; then
		echo "mkdir -p $DIR_SAMPLE_NAME.pgf"
	fi

	#tmp="cd $OUTPUT_DIR"
	#echo $tmp


	tmp="$TOPHAT2 -p $nthreads \
			--b2-very-sensitive \
			--mate-inner-dist $mate_inner_dist \
			--mate-std-dev $mate_std_dev \
			--read-mismatches $read_mismatches \
			--read-edit-dist $read_edit_dist \
			--b2-N 1 \
			--library-type $library_type \
			--read-gap-length $read_gap_length \
			--segment-mismatches $segment_mismatches \
			-G $GTF \
			-o $SAMPLE_NAME.pgf $BOWTIE2_GENOME_INDEX $2 $3 \
			2> stderr.$SAMPLE_NAME.tophat2_pgf.txt"
	echo $tmp


	echo "Manipulations on sam file:  indexing: \`date\`"

	tmp="$SAMTOOLS index $OUTPUT_DIR$SAMPLE_NAME.pgf/accepted_hits.bam \
		2> stderr.$SAMPLE_NAME.samtools_index.txt"
	echo $tmp

	tmp="$SAMTOOLS view -hb $OUTPUT_DIR$SAMPLE_NAME.pgf/accepted_hits.bam \
                        chr6:28500000-33390000 > $OUTPUT_DIR$SAMPLE_NAME.pgf/pgf.accepted_hits.bam"
	echo $tmp

	tmp="$SAMTOOLS index $OUTPUT_DIR$SAMPLE_NAME.pgf/pgf.accepted_hits.bam \
                2> stderr.$SAMPLE_NAME.samtools_index_pgf.txt"
                echo $tmp


	echo "echo INITIAL MAPPING DONE."	
}




###################################################
#
# GET READS FOR SECOND MAPPING TO THE ADDITIONAL MHC HAPLOTYPES
#
###################################################

GET_READS_FOR_SECOND_MAPPING(){

	# hg19:chr6:28502730-33359312

	echo "echo ========================================="
	echo "echo GET READS FOR SECOND MAPPING TO THE ADDITIONAL MHC HAPLOTYPES: \`date\`"



	tmp="$SAMTOOLS view -hb $OUTPUT_DIR$SAMPLE_NAME.pgf/accepted_hits.bam \
                        chr6:28500000-33400000 > $OUTPUT_DIR$SAMPLE_NAME.pgf/pgf.accepted_hits.bam"
        echo $tmp

	tmp="$SAMTOOLS index $OUTPUT_DIR$SAMPLE_NAME.pgf/pgf.accepted_hits.bam \
                2> stderr.$SAMPLE_NAME.samtools_index_pgf.txt"
                echo $tmp

	
	tmp="$SAMTOOLS view $OUTPUT_DIR$SAMPLE_NAME.pgf/accepted_hits.bam \
			chr6:28500000-33400000 | cut -f 1 | sort | uniq \
			> $OUTPUT_DIR$SAMPLE_NAME.pgf/read_ids_for_second_mapping_.txt \
			2> stderr.$SAMPLE_NAME.read_ids_for_second_mapping.txt"	
	echo $tmp

	
	tmp="$SAMTOOLS view \
			$OUTPUT_DIR$SAMPLE_NAME.pgf/unmapped.bam | cut -f 1 | sort | uniq \
			>> $OUTPUT_DIR$SAMPLE_NAME.pgf/read_ids_for_second_mapping_.txt \
			2> stderr.$SAMPLE_NAME.read_ids_for_second_mapping2.txt"
	echo $tmp

	tmp="sort $OUTPUT_DIR$SAMPLE_NAME.pgf/read_ids_for_second_mapping_.txt | uniq \
			> $OUTPUT_DIR$SAMPLE_NAME.pgf/read_ids_for_second_mapping.txt"
	echo $tmp


	READ1="$OUTPUT_DIR$SAMPLE_NAME.pgf/$SAMPLE_NAME.secondMapping_1.fastq"
	READ2="$OUTPUT_DIR$SAMPLE_NAME.pgf/$SAMPLE_NAME.secondMapping_2.fastq"


	
	tmp="$SEQTK subseq $2 $OUTPUT_DIR$SAMPLE_NAME.pgf/read_ids_for_second_mapping.txt \
			> $READ1"
	echo $tmp

	tmp="$SEQTK subseq $3 $OUTPUT_DIR$SAMPLE_NAME.pgf/read_ids_for_second_mapping.txt \
			> $READ2"
	echo $tmp

	tmp="gzip $READ1"
	echo $tmp
	tmp="gzip $READ2"
	echo $tmp

	tmp="rm $OUTPUT_DIR$SAMPLE_NAME.pgf/read_ids_for_second_mapping_.txt"
	echo $tmp
}




###################################################
#
# MAPPING TO OTHER MHC HAPLOTYPES
#
###################################################

TOPHAT2_MAPPING_TO_THE_MHC(){

	echo "echo ========================================="
	echo "echo TOPHAT2 MAPPING TO THE OTHER MHC HAPLOTYPES: \`date\`"

	for hap in cox apd dbb mann mcf qbl ssto;
	do
		dir=$OUTPUT_DIR$SAMPLE_NAME.$hap
		if [[ ! -e $dir ]]; then
			echo "mkdir -p $dir"
		fi
		ref_hap=`echo $hap | tr '[:lower:]' '[:upper:]' `
		ref_hap=$BOWTIE2_MHC$ref_hap

		READ1="$OUTPUT_DIR$SAMPLE_NAME.pgf/$SAMPLE_NAME.secondMapping_1.fastq"
		READ2="$OUTPUT_DIR$SAMPLE_NAME.pgf/$SAMPLE_NAME.secondMapping_2.fastq"

		#tmp="cd $OUTPUT_DIR"
	        #echo $tmp


		tmp="$TOPHAT2 -p $nthreads \
				--b2-very-sensitive \
				--mate-inner-dist $mate_inner_dist \
				--mate-std-dev $mate_std_dev \
                        	--read-mismatches $read_mismatches \
                        	--read-edit-dist $read_edit_dist \
                        	--b2-N 1 \
                        	--library-type $library_type \
	                        --read-gap-length $read_gap_length \
        	                --segment-mismatches $segment_mismatches \
                        	-G $GTF \
                        	-o $SAMPLE_NAME.$hap $ref_hap $READ1.gz $READ2.gz \
                        	2> stderr.$SAMPLE_NAME.tophat2_$hap.txt"
		echo $tmp
	
		tmp="$SAMTOOLS index $dir/accepted_hits.bam \
                2> stderr.$SAMPLE_NAME.samtools_index_$hap.txt"
		echo $tmp

	done

}




###################################################
#
# ENVIRONMENT SETUP
#
###################################################

ENVIRONMENT_SETUP(){
	      echo "source ~/.bashrc"
}




###################################################
#
# WRITE A SCRIPT
#
###################################################

WRITE_A_SCRIPT(){

        echo "#\!/bin/bash"
        echo "#$ -cwd -pe shmem 4 -V -N log."$1
        echo "#$ -P jknight.prjc -q long.qc"

        #ENVIRONMENT_SETUP
	      TOPHAT2_MAPPING $1 $2 $3
	      GET_READS_FOR_SECOND_MAPPING $1 $2 $3
        TOPHAT2_MAPPING_TO_THE_MHC $1 $2 $3


        echo "echo Finished analysis at: \`date\`"

}

WRITE_A_SCRIPT $1 $2 $3 $4







