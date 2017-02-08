#! /usr/bin/env bash
### function
## updated


usage()
{
cat << EOF
######
##	script to tranfer the data to the delivery space and clean the secondary space
##
##	Script Options:
##	-r	<run_info>	-	/path/to/run info file
##	-o	<output dir>	-	/full path/ to output secondary folder
##	-m 	<multisample>	-	flag to specify if multi sample analysis.
##	-h - Display this usage/help text (No arg)
##
#################################################################################################################
EOF
}

echo "Options specified: $@"

while getopts ":o:r:mh" OPTION; do
  case $OPTION in
	m) multisample="YES" ;;
	h) usage
	exit ;;
	o) output=$OPTARG ;;
	r) run_info=$OPTARG ;;
   \?) echo "Invalid option: -$OPTARG. See output file for usage." >&2
       usage
       exit ;;
    :) echo "Option -$OPTARG requires an argument. See output file for usage." >&2
       usage
       exit ;;
  esac
done

#### start of the main wrapper
if [ -z "$run_info" ] ||[ -z "$output" ];
then
	echo "Must provide at least required options. See output file for usage." >&2
	usage
	exit 1;
fi
START=$(date +%s)
#set -x
source $run_info
source $TOOL_INFO
source $MEMORY_INFO
source $SAMPLE_INFO
source $WORKFLOW_PATH/shared_functions.sh

check_variable "$run_info:DELIVERY_FOLDER" $DELIVERY_FOLDER
check_variable "$TOOL_INFO:TAR" $TAR
check_variable "$TOOL_INFO:FIND" $FIND
check_variable "$run_info:PROJECTNAME" $PROJECTNAME 
check_variable "$run_info:DATATYPE" $DATATYPE
check_variable "$TOOL_INFO:GVCF_MODULE" $GVCF_MODULE
check_variable "$run_info:VCF_MINER_GROUP" $VCF_MINER_GROUP 
check_variable "$TOOL_INFO:VCF_MINER_PATH" $VCF_MINER_PATH
if [[ ! "$DELIVERY_FOLDER" || $DELIVERY_FOLDER == "NA" ]]
then
	echo -e "DELIVERY_FOLDER is needed to transfer the data and do clean up, it can not be blanked or NA"
	exit 1;
fi	

### check for read, write and execute permission for delivery folder
if [ ! -r ${DELIVERY_FOLDER} ]
then
	echo "Read access permission denied on ${DELIVERY_FOLDER} Directory"
	exit 1;
fi
if [ ! -w ${DELIVERY_FOLDER} ]
then
	echo "Write access permission denied on ${DELIVERY_FOLDER} Directory"
	exit 1;
fi
if [ ! -x ${DELIVERY_FOLDER} ]
then
	echo "Execute access permission denied on ${DELIVERY_FOLDER} Directory"
	exit 1;
fi

# Upload VCF file to VCF Miner
# Transfer this file first so it can prompt the user for a password right away
vcffile=`ls $output/variants/$PI.*.variants.vcf.gz`
if [ ! -f $vcffile ]
then
	echo "WARNING: Cannot find $vcffile to upload to VCF Miner. You will have to do this manually." 
else
	echo "uploading VCF file to VCF miner, user needs to enter the password once prompted to do so"
	### first create the group
	
	if [ `$VCF_MINER_PATH -f listGroups -u $USER | grep $VCF_MINER_GROUP | wc -l` -gt 0 ]
	then	
		$VCF_MINER_PATH -f createGroup -u $USER -g $VCF_MINER_GROUP -m $USERNAMES_VCF_MINER
	fi
	### add VCF 
	$VCF_MINER_PATH -f addVCF -u $USER -v $vcffile 
	### add uploaded VCF file to the group created above
	$VCF_MINER_PATH -f addVCFToGroup -u $USER -g $VCF_MINER_GROUP -v $vcffile
	echo "uploading complete"
fi

echo "now cleaning and delivery the data to delivery space specified in the run info file"
echo -e "delivery :\n$DELIVERY_FOLDER "
### 
mkdir -p $DELIVERY_FOLDER/bam
### check for BAM file existence
#for sample in `echo $SAMPLENAMES| tr ":" " "`
#do
#	if [ ! -s $output/bam/$sample.bam ]; then echo "$output/bam/$sample.bam : file doesn't exist"; exit 1; fi
#	if [ ! -s $output/bam/$sample.bam.bai ]; then echo "$output/bam/$sample.bam.bai : file doesn't exist"; exit 1; fi
#done

#echo -e "transferring BAM files ...."
#### transfer the BAM files
#for sample in `echo $SAMPLENAMES | tr ":" " "`
#do
#	mv $output/bam/$sample.bam $DELIVERY_FOLDER/bam/
#	mv $output/bam/$sample.bam.bai $DELIVERY_FOLDER/bam/
#done

#rm -Rf $output/bam
### delete files
for i in `$FIND $output/qc/fastqc -type f -name '*.txt'`
do
	rm -Rf $i
done
echo "Still working on the clean-up..."

### delete idx files
#for i in `$FIND $output/variants/gVCF/ -type f -name '*.idx'`
#do
 #   rm -Rf $i
#done

### delete the QC directory
if [ $RUN_QC == "YES" ]
then
	for i in figures inputs istats sexck summary   
	do
		rm -Rf $output/qc/ngsqc/$i
	done
fi
### cp the gVCF files to ORACLE space
if [ $GVCF_MODULE == "YES" ]
then
	check_variable "$TOOL_INFO:ORACLE_TRC_DELIVERY" $ORACLE_TRC_DELIVERY
	for sample in `echo $SAMPLENAMES | tr ":" " "`
	do
		mv $output/variants/gVCF/$sample/* $ORACLE_TRC_DELIVERY/
		rm -Rf $output/variants/gVCF/$sample
	done
fi
### make sure to delete the per chromosome folder for gvCF file creation
if [ $DO_ANALYSIS_PER_CHROMOSOME == "YES" ]
then	
	for chr in `echo $CHRINDEX | tr ":" " "`
	do
		rm -Rf $output/variants/gVCF/chr$chr
        rm -Rf $output/variants/chr$chr
	done
fi
	
echo "Still working on the clean-up..."
	
for folder in docs variants qc
do
    cp -Rf $output/$folder $DELIVERY_FOLDER/
	flag=`[ "$(ls -A $DELIVERY_FOLDER/$folder)" ] && echo "Not Empty" || echo "Empty"`
	if [[ "$flag" == "Empty" ]]
	then
		echo -e "User doesn't have access to create the folder and copy files in $DELIVERY_FOLDER/"
	else	
		if [[ "$folder" != "docs" ]]
		then
			rm -Rf $output/$folder
		fi	
	fi
done

check_variable "$run_info:VERSION" $VERSION

mv $output/index.html $DELIVERY_FOLDER/
mv $output/README $DELIVERY_FOLDER/

##update the dashboard
$WORKFLOW_PATH/dashboard.sh -T $TOOL_INFO -M $MEMORY_INFO -A $PROJECTNAME -s Delivered -d $DATATYPE

### tar ball the logs folder
cd $output/
### delete all the core files if there are 
rm -Rf .logs/core.*
echo "Still working on the clean-up and making TAR balls..."
$TAR -czf .logs.tar.gz .logs
$TAR -czf docs.tar.gz docs
rm -Rf $output/.logs
rm -Rf $output/docs
rm -Rf $output/.tmp

for i in `$FIND $DELIVERY_FOLDER -type d -name 'temp*'` 
do
	rm -Rf $i
done	
for i in `$FIND $DELIVERY_FOLDER -type d -name "$USER"`
do
	rm -Rf $i
done	

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "clean up and transfer for $output work flow execution took $DIFF seconds"