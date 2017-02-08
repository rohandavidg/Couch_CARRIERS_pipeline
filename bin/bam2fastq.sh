#!/bin/sh
# generate fastqs from bams

if [ $# != 2 ];
then
	echo "usage: <input dir> <output_dir>";
else					
#	set -x
	input_dir=$1
	output_dir=$2
	picard=/projects/bsi/bictools/apps/alignment/picard/1.55/
	bams=`ls $input_dir/*.bam`
	for bam in $bams;
	do
	    echo
	    sample_name=`basename $bam | cut -d"." -f1,2`
	    qsub -N bamtofastq$sample_name -q 7-days -l h_vmem=25G -b y -M gnanaolivu.rohandavid@mayo.edu  -l h_stack=10M -pe threaded 4 -e /data2/bsi/secondary/Couch_Fergus_coucf/exome/config-170207/logs -o /data2/bsi/secondary/Couch_Fergus_coucf/exome/config-170207/raw_fastq  -wd /data2/bsi/secondary/Couch_Fergus_coucf/exome/config-170207/logs -M gnanaolivu.rohandavid@mayo.edu -V /usr/local/biotools/java/jdk1.6.0_05/bin/java -Xmx12g -Xms512m -jar $picard/SamToFastq.jar INPUT=$bam FASTQ=$output_dir/$sample_name.R1.fastq SECOND_END_FASTQ=$output_dir/$sample_name.R2.fastq VALIDATION_STRINGENCY=SILENT TMP_DIR=$output_dir/	
	done
fi
echo
#	echo gzip $output_dir/$sample.R1.fastq
#	echo gzip $output_dir/$sample.R2.fastq


