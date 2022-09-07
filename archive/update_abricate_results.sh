#!/usr/bin/env bash

help () {
    printf "\n\n"
    printf "Working directory containing only FASTAs to analyze\n"
    printf "\-c percent coverage cutoff for abricate (default 0)\n\n"
    printf "\-f percent identity cutoff for abricate (default 75)\n\n"
    printf "\-p file path for the sample directories- used to copy files after the run is completed."
    exit 1
}

hflag=
cflag=
fflag=
pflag=

while getopts ':hc:f:p:' OPTION; do
    case $OPTION in

    h) hflag=1
    ;;
    c) cflag=$OPTARG
    ;;
    f) fflag=$OPTARG
    ;;
    p) pflag=$OPTARG
    ;;

#        ?) echo "Invalid option: -$OPTARG" > &2
#        ;;
    esac
done

if [ "$hflag" ] 
then 
    help 
fi

#setup coverage and identity variables for abricate
if [[ $cflag ]]; then
    c=$cflag
else
    c=1
fi

if [[ $fflag ]]; then
    f=$fflag
else
    f=75
fi

root_dir=`pwd`

fasta_count=`ls *.fasta |wc -l`
echo "####$fasta_count files to process.####"

# Run abricate
#cd ${root_dir}/abricate
#for i in *; do name=`echo $i`; mv $i $name; done

NR_CPUS=10
# Default resfinder
for i in *fasta; do 
    (abricate --mincov $c --minid $f $i > ${i%.fasta}-resfinder.tab) & let count+=1
    [[ $((count%NR_CPUS)) -eq 0 ]] && wait
done
wait
sleep 10
abricate --summary *-resfinder.tab > resfinder_summary.txt

# #NCBI
for i in *fasta; do 
    (abricate --db ncbi --mincov $c --minid $f $i > ${i%.fasta}-ncbi.tab) & let count+=1
    [[ $((count%NR_CPUS)) -eq 0 ]] && wait
done
wait
sleep 10 
abricate --summary *-ncbi.tab > ncbi_summary.txt

sleep 10

abricate_to_excel.py

for i in *.fasta; do
    n=`echo $i | sed 's/.fasta//'`
    run_info_report.py -c $c -f $f -i ${n}
    echo "Copying ${n} to ${pflag}${n}/"
    cp abricate*.xlsx ${pflag}${n}/
    cp ${n}-run_info*.pdf ${pflag}${n}/
done



