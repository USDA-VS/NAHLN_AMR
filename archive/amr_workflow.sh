#!/usr/bin/env bash

#To be used on sever2
#Script is NOT portable. Local environment specific

: <<'COMMENT'
CHECK FOR DEPENDCIES:
~$ svn --version
~$ pip --version
~$ virtualenv --version
~$ docker run hello-world

~$ source cwl/bin/activate
(cwl) ~$ pip install -U wheel setuptools
(cwl) ~$ pip install -U cwltool[deps] PyYAML cwlref-runner

python package: svgwrite

If NAHLM samples... separate after running amr_workflow.sh on all samples.  Separate as EC (ecoli), MH, SIG (staph intermedius group) then run nahln_amr_updates.sh on each group.

COMMENT

help () {
    printf "\n\n"
    printf "Working directory containing only FASTQs to analyze\n"
    printf "Illumina pair reads only\n"
    printf "\nUsage:\n"
    printf "[wd]$ ls\n"
    printf "\tsampleA_R1.fastq.gz\n"
    printf "\tsampleA_R2.fastq.gz\n"
    printf "\tsampleB_R1.fastq.gz\n"
    printf "\tsampleB_R2.fastq.gz\n"
    printf "[wd]$ amr_workflow.sh\n\n"
    printf "\-i run ISEScan\n\n"
    printf "\-c percent coverage cutoff for abricate (default 0)\n\n"
    printf "\-f percent identity cutoff for abricate (default 75)\n\n"
    printf "\-m percent coverage cutoff for FDA pipeline (default 0.5)\n\n"
    printf "\-d percent identity cutoff for FDA pipeline (default 85)\n\n"
    exit 1
}

#Check arg1 for database type
#if [[ $1 == -h ]]; then
#    help
#elif [[ $1 == -help ]]; then
#    help
#elif [[ $1 == --help ]]; then
#    help
#fi

iflag=
hflag=
iflag=
cflag=
fflag=
mflag=
dflag=
nflag=
while getopts ':hic:f:m:d:n:' OPTION; do
    case $OPTION in
    i) iflag=1
    ;;
    h) hflag=1
    ;;
    i) iflag=1
    ;;
    c) cflag=$OPTARG
    ;;
    f) fflag=$OPTARG
    ;;
    m) mflag=$OPTARG
    ;;
    d) dflag=$OPTARG
    ;;
    n) nflag=$OPTARG
    ;;
#        ?) echo "Invalid option: -$OPTARG" > &2
#        ;;
    esac
done

if [ "$hflag" ] 
then 
    help 
fi

#activate amrfinder
source /usr/local/bin/cwl/bin/activate

root_dir=`pwd`

fasta_count=`ls *.fasta |wc -l`
fastq_count=`ls *.fastq |wc -l`

if [[ $fasta_count -gt 0 ]] && [ $fastq_count -gt 0 ]; then
    echo "FASTAs and FASTQs found. Please run different file types separately."
    exit 1
fi

if [[ $fasta_count -gt 0 ]]; then
    printf "\n\nFASTAs found!!! WARNING: SeqSero Results based on FASTAs may not be as accurate as results based on FASTQs. -> \n\n"

    for i in *.fasta; do
        cd $root_dir
        n=`echo $i | sed 's/.fasta//'`
        echo "n is : $n"
        echo $n >> sample_list
        mkdir -p $n
        mv $i $n/
        echo "working on: $i"; cd ./$n
        mlst *.fasta >> ${n}_mlst.txt
        #/usr/local/bin/amr_finder/amrfinder -n ${i}.fasta >>  ${i}_amrfinder.tab
        SeqSero2_package.py -t 4 -i *.fa*
    done

else

#Taken from packagefastq.sh
printf "\n\nPackaging FASTQs -> \n\n"
for i in *.fastq*; do 
    n=`echo $i | sed 's/_.*//' | sed 's/\..*//'`
    #n=`echo $i | sed 's/_R[1,2].fastq//'`
    echo "n is : $n"
    echo $n >> sample_list
    mkdir -p $n 
    mv $i $n/
done

sort < sample_list | uniq > sample_list_temp; mv sample_list_temp sample_list

#Taken from loop_at_terminal.txt
printf "\n\nRunning SPAdes -> \n\n"
NR_CPUS=10
for i in *; do 
    (cd $root_dir
    echo "working on: $i"; cd ./$i
    simple_spades.sh) & let count+=1
    [[ $((count%NR_CPUS)) -eq 0 ]] && wait
done
wait
cd $root_dir

#Package SPAdes
printf "\n\nPackaging SPAdes -> \n\n"
for i in *; do
    mkdir ${i}/assembly
    mv ${i}/*txt ${i}/*fasta ${i}/*log ${i}/assembly
    cd ${i}/assembly/
    mlst ${i}.fasta >> ${i}_mlst.txt
    #/usr/local/bin/amr_finder/amrfinder -n ${i}.fasta >> ${i}_amrfinder.tab
    contig_length.sh
    cd $root_dir
done

mkdir ${root_dir}/temp

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

#Run Kraken - host db
printf "\n\nRunning Kraken - host db -> \n\n"
for i in `cat sample_list`; do 
    cd $root_dir
    echo "Kraken - host db running on $f"; cd ./${i}
    cp original_reads/* ./
    simple_kraken.sh host na 
    pwd
    mkdir kraken
    mv *txt *html kraken
#    printf "Kaiju running...\n"
#    simple_kaiju.sh na
#    mkdir kaiju
#    mv *out* *html kaiju
    echo "SeroSeq2 - running on $i"
    SeqSero2_package.py -t 2 -i ${i}*_R1*.fastq.gz ${i}*_R2*.fastq.gz
    cd ./SeqSero_result*
    mv Seqsero_result.txt ${i}_Seqsero_result.txt
    cd ../
    echo "Finding FASTQ quality..."
    fastq_quality_amr.py -r1 ${i}*_R1*.fastq.gz -r2 ${i}*_R2*.fastq.gz
    wait
#    mkdir original_reads
#    mv *.fastq.gz ./original_reads/
done

pwd
echo $root_dir
cd $root_dir
pwd

#Collect stats

mkdir ./stats
find . -name "*stat_quality*.txt" -exec cp {} stats/ \;
find . -name "*assemblathon_reformat_stats.txt" -exec cp {} stats/ \;

cd ./stats
for i in `cat ../sample_list`; do
    merge_quality_files.py -a ${i}*assemblathon_reformat_stats.txt -q ${i}*stat_quality*.txt
    rm ${i}*assemblathon_reformat_stats.txt
    rm ${i}*stat_quality.txt
done

cat ./*_merged.txt | head -1 > ./batch_allmerged.txt
for i in ./*_merged.txt; do
    cat ${i} | tail -1 >> ./batch_allmerged.txt
done

stats_all_to_excel.py

cd ${root_dir}

for i in `cat sample_list`; do
    cp ${root_dir}/stats/stat-${i}*xlsx ${i}
done
#Finish FASTQ procesing
fi

cd ${root_dir}

mkdir abricate
#mkdir fastas_concatenated card-rgi abricate
#find . -name "*.fasta" -exec cp {} ./fastas_concatenated \;
#find . -name "*.fasta" -exec cp {} ./card-rgi \;
find . -name "*.fasta" -exec cp {} ./abricate \;
scaffold_number=`ls fastas_concatenated/ | wc -l`
#cd ./fastas_concatenated
#for i in *fasta; do 
#    name=`echo $i | sed s/[_.].*//`
#    echo ">${name}" > ${name}.fa; grep -v "^>" $i | tr -d "\n" >> ${name}.fa
#    printf "\n\n" >> ${name}.fa
#done
#rm *fasta
#if [[ $mflag ]]; then 
#    m=$mflag 
#else 
#    m=0.5 
#fi

#if [[ $dflag ]]; then 
#    d=$dflag 
#else 
#    d=85 
#fi

#export m d
#Run FDA pipeline on concatenated scaffolds
#echo "FDA pipeline running..."
#resistance_genes_detection_new.sh
#mkdir fastas
#mv *.fa fastas
#fda2_to_excel.py
#Copy the class.txt file to an Excel worksheet.  Name as fda_abricate_YYYY-MM-DD.xlsx
#cd $root_dir

#CARD-RGI
#NR_CPUS=10
#cd card-rgi
#echo "CARD-RGI running..."
#for i in *.fasta; do 
#    (echo "$i is running"
#    name=`echo $i | sed 's/[_.].*//'`
#    mkdir $name
#    rgi main --input_sequence $i --output_file $name --input_type contig
#    mv ${i}* $name) & let count+=1
#    [[ $((count%NR_CPUS)) -eq 0 ]] && wait
#done
#wait

#cd ${root_dir}/card-rgi
#mkdir card_summary_files
#find . -name "*.txt" -exec cp {} card_summary_files \;
#cp "${root_dir}/fastas_concatenated/class.txt" card_summary_files
#cd card_summary_files
#mv ${root_dir}/fastas_concatenated/*.xlsx .
#fda_to_excel.py class.txt
#rm ${root_dir}/card-rgi/card_summary_files/class.txt
#for i in *txt; do card_to_excel.py $i; done
#rm ${root_dir}/card-rgi/card_summary_files/*txt
#cd ${root_dir}
#for i in `cat sample_list`; do
#    cp ${root_dir}/card-rgi/card_summary_files/card-${i}*.xlsx ${root_dir}/card-rgi/card_summary_files/fda-*.xlsx ${i}
#done

#cd ${root_dir}
#mv ./card-rgi ./fastas_concatenated ${root_dir}/temp

# Run abricate
cd ${root_dir}/abricate
for i in *; do name=`echo $i`; mv $i $name; done

if [ "$iflag" ] 
then
# Make a copy for ISEScan
cp -r ${root_dir}/abricate ${root_dir}/isescan
fi

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

 #plasmidfinder
for i in *fasta; do 
    (abricate --db plasmidfinder --mincov $c --minid $f $i > ${i%.fasta}-plasmidfinder.tab) & let count+=1
    [[ $((count%NR_CPUS)) -eq 0 ]] && wait
done
wait
sleep 10
abricate --summary *-plasmidfinder.tab > plasmidfinder_summary.txt


#activate AMRFinder
#source /usr/local/bin/cwl/bin/activate

#amrfinder
for i in *.fasta; do
#   /usr/local/bin/amr_finder/amrfinder -n $i >> ${i%.fasta}_amrfind.tab
   /usr/local/bin/amrfinder/amrfinder.pl -n $i >> ${i%.fasta}_amrfind.tab
done
wait

for i in *.fasta; do
   cat ${i%.fasta}_amrfind.tab | sed 's/Contig id/SEQUENCE/g' > ${i%.fasta}-amrfinder.tab
done


abricate_to_excel.py
amrfindermerged=`echo abricate-amrfinder*.xlsx | sed 's/abricate-//'`
mv abricate-amrfinder*.xlsx ${amrfindermerged}

cd ${root_dir}
for i in `cat sample_list`; do
    cp ${root_dir}/abricate/*xlsx ${i}
    cp ${root_dir}/abricate/${i}* ${i}
done

if [ "$iflag" ]
then
#ISEScan
cd ${root_dir}/isescan
root_isescan=`pwd`

#Taken from packagefastq.sh
printf "\n\nPackaging FASTAs -> \n\n"
for i in *; do 
    n=`echo $i | sed 's/_.*//' | sed 's/\..*//'`
    echo "n is : $n"
    mkdir -p $n 
    mv $i $n/
done

printf "\n\nRunning ISEScan -> \n\n"
NR_CPUS=10
for i in *; do 
    (cd $root_isescan
    echo "working on: $i"; cd ./$i
    isescan.py ${i}.*fa* proteome hmm) & let count+=1
    [[ $((count%NR_CPUS)) -eq 0 ]] && wait
done
wait

#ISEScan results to Excel file
NR_CPUS=10
for i in *; do
    (cd $root_isescan/${i}/prediction;
    amr_update_isescan_files.sh;
    isescan_to_excel.py) & let count+=1
    [[ $((count%NR_CPUS)) -eq 0 ]] && wait
done
wait

#Copy ISEScan analysis to sample data folders
cd ${root_dir}
for i in `cat sample_list`; do
    cp ${root_dir}/isescan/${i}/prediction/*xlsx ${root_dir}/${i}
    cp -r ${root_dir}/isescan/${i} ${root_dir}/${i}/isescan
done
fi

cd $root_dir

if (( $fasta_count == 0 )) 
then
    for i in `cat sample_list`; do
        cd $root_dir
        echo "Making report for ${i}"
        cd ./$i
        mkdir original_reads
        mv *.fastq.gz ./original_reads/
        find . -name "*assemblathon_reformat_stats.txt" -exec cp {} ./ \;
        find . -name "*Seqsero_result.txt" -exec cp {} ./ \;
        find . -name "*_mlst.txt" -exec cp {} ./ \;
        find . -name "*_lk.txt" -exec cp {} ./ \;
        find . -name "*_amrfinder.txt" -exec cp {} ./ \;
        sed '1d' *-ncbi.tab | cut -f 2-6,8-10,12,13 | sed 's/_/\\_/g' | tr "\t" "&" >> ${i}-ncbi-report.tab
        sed '1d' *-resfinder.tab | cut -f 2-6,8-10,12,13 | sed 's/_/\\_/g' | tr "\t" "&"  >> ${i}-resfinder-report.tab
        sed '1d' *-amrfinder.tab | cut -f 2-4,6-14,16 | sed 's/_/\\_/g' | tr "\t" "&" >> ${i}-amrfinder-report.tab
        if [[ $nflag ]]; then
           n = $nflag
           amr_report.py -c $c -f $f -n $n -a ${i}*assemblathon_reformat_stats.txt -q ${i}*stat_quality*.txt -s ${i}*_Seqsero_result.txt -m ${i}_mlst.txt -l ${i}_lk.txt 
        else
           amr_report.py -c $c -f $f -a ${i}*assemblathon_reformat_stats.txt -q ${i}*stat_quality*.txt -s ${i}*_Seqsero_result.txt -m ${i}_mlst.txt -l ${i}_lk.txt
        fi
        mkdir report_files
        mv ${i}.aux report_files
        mv ${i}.log report_files
        mv ${i}.out report_files
        mv ${i}.tex report_files
        mv ${i}*.tab report_files
        mv ${i}_lk.txt report_files
        mv *_mlst.txt report_files
        mv *_Seqsero_result.txt report_files
        mv *_stat_quality30.txt report_files
        mv *-assemblathon_reformat_stats.txt report_files
        rm *.png
        rm *.svg
    done
fi

cd $root_dir 

#Clean-up
mv sample_list ./temp
mv ${root_dir}/abricate temp
mv ${root_dir}/isescan temp
#rm -r ${root_dir}/original_reads
rm *.txt
# tstuber 2018-04-13
