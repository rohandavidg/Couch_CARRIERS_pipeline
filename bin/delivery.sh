#! /usr/bin/env bash


usage()
{
cat << EOF
######
##	script to create final deliverbale and put everything in place adn create html page
##
##	Script Options:
##	-r	<run_info>	-	/path/to/run info file
##	-o	<output dir>	-	/full path/ to output folder
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

source $run_info
source $TOOL_INFO
source $WORKFLOW_PATH/shared_functions.sh
check_variable "$run_info:TOOL_INFO" $TOOL_INFO
check_variable "$run_info:PROJECTNAME" $PROJECTNAME
check_variable "$run_info:VERSION" $VERSION
check_variable "$run_info:MEMORY_INFO" $MEMORY_INFO
check_variable "$run_info:SAMPLE_INFO" $SAMPLE_INFO
check_variable "$TOOL_INFO:WORKFLOW_PATH" $WORKFLOW_PATH
check_variable "$TOOL_INFO:PERL" $PERL
check_variable "$TOOL_INFO:R" $R
check_variable "$run_info:DATATYPE" $DATATYPE
check_variable "$TOOL_INFO:MAIL" $MAIL
check_variable "$TOOL_INFO:FINGER" $FINGER
check_variable "$TOOL_INFO:FIND" $FIND
check_variable "$TOOL_INFO:SV_MODULE" $SV_MODULE
check_variable "$TOOL_INFO:DEBUG_MODE" $DEBUG_MODE
check_variable "$TOOL_INFO:MULTISAMPLE_REALIGN_RECAL" $MULTISAMPLE_REALIGN_RECAL
check_variable "$tool_info:PYTHON" $PYTHON
check_variable "$tool_info:PYTHONLIB" $PYTHONLIB
check_variable "$tool_info:DELIVERY_GENERATOR_PATH" $DELIVERY_GENERATOR_PATH
#set -x
### setting the environment 
export PYTHONPATH=$PYTHONLIB:$PYTHONPATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PYTHONPATH
export PATH=$PYTHON:$PYTHONPATH:$PATH

if [ $DASHBOARD == "YES" ]
then
	#### dashboad update for Results
	$WORKFLOW_PATH/dashboard.sh -T $TOOL_INFO -M $MEMORY_INFO -A $PROJECTNAME -s Results -d $DATATYPE
fi
### copy SV results
sv_results="NO"
if [[ $SV_MODULE == "YES"  && $DATATYPE == "exome" ]]
then
	if [ $SOMATIC_CALLING == "NO" ]
	then
		if [[ `echo $SAMPLENAMES | tr ":" "\n" | wc -l` -ge 3 ]]
		then
			sv_results="YES"
		fi
	else
		sv_results="YES"
	fi	
fi
if [ $sv_results == "YES" ]
then
	## copy the cnv results
	if [ "$(ls -A $output/.tmp/cnv/cnv-plot/)" ]
	then
		mv $output/.tmp/cnv/cnv-plot/* $output/qc/cnv/
	fi
	if [ "$(ls -A $output/.tmp/cnv/cnv-txt/)" ]
	then	
		mv $output/.tmp/cnv/cnv-txt/* $output/variants/cnv/
	fi
fi	

## create igv session for the workflow execution
if [ "$multisample" ]
then
	s_param="-r $run_info -o $output/docs/igv_session/ -m"
else
	s_param="-r $run_info -o $output/docs/igv_session/"
fi
echo $s_param | xargs $WORKFLOW_PATH/igv_session.sh
if [ $? -ne 0 ]
then
	$WORKFLOW_PATH/email.sh -f delivery_folder -m delivery.sh -s igv_session.sh -p "$s_param" -l $LINENO
	exit 100;
fi	
#### generate the coverage plot
samples=`echo $SAMPLENAMES | tr ":" " "` 	
### calculate the value of region covered
if [ $DATATYPE == "exome" ]
then
	file=`echo $CAPTUREKIT`
	region=`cat $file | awk '{sum+=$3-$2;print sum}'| tail -1`
else
	file=`echo $REF_GENOME | awk '{print $0".fai"}'`
	region=`cat $REF_GENOME | grep -v '>' | grep -E '^[ACTG]+$' | perl -pe 's/[[:space:]]//g' | wc -c`
	#region=`cat $file| tail -n1 | awk -F'\t' '{print $3}'` 
fi
$R/Rscript /data2/bsi/staff_analysis/m149947/couch_code/Couch_CARRIERS_pipeline/bin/coverage_plot.r $output/qc/coverage/ $output/docs/coverage.png $region $samples
if [ $? -ne 0 ]
then
	$WORKFLOW_PATH/email.sh -f $output/qc/coverage/coverage.png -m delivery.sh -s coverage_plot.r -p "$output/qc/coverage/ $output/qc/coverage/coverage.png $region $samples" -l $LINENO
	exit 100;
fi	
### generate the delivery document to deliver to PI
### generate the new documentation to deliver to the PI
export PYTHONPATH=$PYTHONPATH:$PYTHONLIB:$DELIVERY_GENERATOR_PATH
echo -e "OUTPUT_FOLDER=$output" | cat $run_info -  | sed -e 's/\"//g' > $output/run_info.txt
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PYTHONPATH $PYTHON/python $DELIVERY_GENERATOR_PATH/custom/ggps/build_ggps_document.py --output $output/ggps.json --run_info $output/run_info.txt --deliverydir $output
if [ $? -ne 0 ]
then
	$WORKFLOW_PATH/email.sh -f $output/ggps.json -m delivery.sh -s build_ggps_document.py -p "--output $output/ggps.json --run_info $run_info --deliverydir $output" -l $LINENO
	exit 100;
fi

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PYTHONPATH $PYTHON/python $DELIVERY_GENERATOR_PATH/generator.py --input $output/ggps.json
 
if [ $? -ne 0 ]
then
	$WORKFLOW_PATH/email.sh -f $output/ggps.json -m delivery.sh -s generator.py -p "--input $output/ggps.json " -l $LINENO
	exit 100;
else
	if [ -s $output/run_info.txt ]; then rm $output/run_info.txt; fi  
	if [ -s $output/ggps.json ]; then rm $output/ggps.json; fi
fi

$WORKFLOW_PATH/dashboard.sh -T $TOOL_INFO -M $MEMORY_INFO -A $PROJECTNAME -s Results -d $DATATYPE -c 

### delete the intermediate file
for sample in `echo $SAMPLENAMES | cut -d '=' -f2 | tr ":" " "`
do
	if [ -s $output/bam/$sample.mappedreads.txt ]; then rm $output/bam/$sample.mappedreads.txt; fi
done

### generate the README for the files and folders in deliverables
$WORKFLOW_PATH/generate_readme.sh $output $run_info
	
### delete the temp folders and files 
for sample in `echo $SAMPLENAMES | cut -d '=' -f2 | tr ":" " "`
do
	if [ -s $output/bam/$sample.extra.bam  ]; then rm $output/bam/$sample.extra.bam; fi
	if [ -s $output/bam/$sample.extra.bam.bai ]; then rm $output/bam/$sample.extra.bam.bai; fi
done 

if [ $DEBUG_MODE == "NO" ]
then
	### delete all the temp folders
	$FIND $output -type d -name 'temp*' | xargs rm -Rf
fi
### main document and qc is generated
rm -Rf $output/qc/numbers $output/qc/coverage
 
### cleanup the ngsqc diectory
mv $output/qc/ngsqc/summary/DataDictionary.xlsx $output/qc/ngsqc/
mv $output/qc/ngsqc/summary/ngsqc.pdf $output/qc/ngsqc/ 
mv $output/qc/ngsqc/summary/SampleMetrics.tsv $output/qc/ngsqc/ 
 
### check for any warnings and errors in the log files
##TODO


#### send completion email to the IS
email=`$FINGER $USER | awk -F ';' '{print $2}' | head -n1`
echo -e "work flow execution is completed for $run_num \nPlease check the output folder : $output\n\nThank you\nWorkflow team" | $MAIL -s "GENOMEGPS v$VERSION work flow execution Status" "$email" 
sleep 30s
echo -e "Analysis completed at:\n" >> $output/docs/log.txt 
echo `date` >> $output/docs/log.txt

$WORKFLOW_PATH/dashboard.sh -T $TOOL_INFO -M $MEMORY_INFO -A $PROJECTNAME -s Complete -d $DATATYPE -c 
END=$(date +%s)
DIFF=$(( $END - $START ))
echo "delivery document generation for workflow execution took $DIFF seconds"	