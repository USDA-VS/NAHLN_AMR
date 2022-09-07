#!/usr/bin/env bash

currentdir=`pwd`
DATE=`date +%Y-%m-%d`

ecoli=`ls EC* | wc -l`
mannh=`ls MH* | wc -l`
staph=`ls SIG* | wc -l`

#Run ecoli samples through extra steps
if [[ $ecoli -gt 0 ]]; then
    for f in EC*; do 
        cd $currentdir 
        echo $f 
        cd ./$f
        abricate --db ecoh *.fasta >> ${f}_ecoh.txt
        abricate --db ecoli_vf *.fasta >> ${f}_ecoli_vf.txt
        abricate --db vfdb *.fasta >> ${f}_vfdb.txt
        rm -r *Seq*
    done

    cd $currentdir
    mkdir copy_to_stats

    find . -name *-ncbi.tab -exec cp {} ./copy_to_stats/ \;
    find . -name *-resfinder.tab -exec cp {} ./copy_to_stats/ \;
    find . -name *-amrfinder.tab -exec cp {} ./copy_to_stats/ \;
    find . -name *-plasmidfinder.tab -exec cp {} ./copy_to_stats/ \;
    find . -name *_ecoh.txt -exec cp {} ./copy_to_stats/ \;
    find . -name *_ecoli_vf.txt -exec cp {} ./copy_to_stats/ \;
    find . -name *_vfdb.txt -exec cp {} ./copy_to_stats/ \;
    find . -name abricate*.xlsx -exec cp {} ./copy_to_stats/ \;
    find . -name amrfinder*.xlsx -exec cp {} ./copy_to_stats/ \;

    cd ./copy_to_stats/
    ncbifile=`echo $DATE"_EC_ncbi_concat.txt"`
    for i in *-ncbi.tab; do
        sed '1d' ${i} >> $ncbifile
    done
    rm *-ncbi.tab

    resfinderfile=`echo $DATE"_EC_resfinder_concat.txt"`  
    for i in *-resfinder.tab; do
        sed '1d' ${i} >> $resfinderfile
    done
    rm *-resfinder.tab

    plasmidfinderfile=`echo $DATE"_EC_plasmidfinder_concat.txt"`
    for i in *-plasmidfinder.tab; do
        sed '1d' ${i} >> $plasmidfinderfile
    done  
    rm *-plasmidfinder.tab        

    amrfinderfile=`echo $DATE"_EC_amrfinder_concat.txt"`
    for i in *-amrfinder.tab; do
        sed '1d' ${i} >> $amrfinderfile
    done
    rm *-amrfinder.tab

    ecohfile=`echo $DATE"_EC_ecoh_concat.txt"`
    for i in *_ecoh.txt; do
        sed '1d' ${i} >> $ecohfile
    done
    rm *_ecoh.txt

    ecolivffile=`echo $DATE"_EC_ecolivf_concat.txt"`
    for i in *_ecoli_vf.txt; do
        sed '1d' ${i} >> $ecolivffile
    done
    rm *_ecoli_vf.txt

    vfdbfile=`echo $DATE"_EC_vfdb_concat.txt"`
    for i in *_vfdb.txt; do
        sed '1d' ${i} >> $vfdbfile
    done
    rm *_vfdb.txt
fi

#Run Mannheimia through extra steps
if [[ $mannh -gt 0 ]]; then
    for f in MH*; do 
        cd $currentdir 
        echo $f 
        cd ./$f
        abricate --db vfdb *.fasta >> ${f}_vfdb.txt
        rm -r *Seq*
    done

    cd $currentdir
    mkdir copy_to_stats

    find . -name *-ncbi.tab -exec cp {} ./copy_to_stats/ \;
    find . -name *-resfinder.tab -exec cp {} ./copy_to_stats/ \;
    find . -name *-amrfinder.tab -exec cp {} ./copy_to_stats/ \;
    find . -name *-plasmidfinder.tab -exec cp {} ./copy_to_stats/ \;
    find . -name *_vfdb.txt -exec cp {} ./copy_to_stats/ \;
    find . -name abricate*.xlsx -exec cp {} ./copy_to_stats/ \;
    find . -name amrfinder*.xlsx -exec cp {} ./copy_to_stats/ \;

    cd ./copy_to_stats/
    ncbifile=`echo $DATE"_MH_ncbi_concat.txt"`
    for i in *-ncbi.tab; do
        sed '1d' ${i} >> $ncbifile
    done
    rm *-ncbi.tab

    resfinderfile=`echo $DATE"_MH_resfinder_concat.txt"`  
    for i in *-resfinder.tab; do
        sed '1d' ${i} >> $resfinderfile
    done
    rm *-resfinder.tab

    plasmidfinderfile=`echo $DATE"_MH_plasmidfinder_concat.txt"`
    for i in *-plasmidfinder.tab; do
        sed '1d' ${i} >> $plasmidfinderfile
    done   
    rm *-plasmidfinder.tab        

    amrfinderfile=`echo $DATE"_MH_amrfinder_concat.txt"`
    for i in *-amrfinder.tab; do
        sed '1d' ${i} >> $amrfinderfile
    done
    rm *-amrfinder.tab

    vfdbfile=`echo $DATE"_MH_vfdb_concat.txt"`
    for i in *_vfdb.txt; do
        sed '1d' ${i} >> $vfdbfile
    done
    rm *_vfdb.txt
fi

#Run SIG through extra steps
if [[ $staph -gt 0 ]]; then
    for f in SIG*; do 
        cd $currentdir 
        echo $f 
        cd ./$f
        abricate --db vfdb *.fasta >> ${f}_vfdb.txt
        rm -r *Seq*
    done

    cd $currentdir
    mkdir copy_to_stats

    find . -name *-ncbi.tab -exec cp {} ./copy_to_stats/ \;
    find . -name *-resfinder.tab -exec cp {} ./copy_to_stats/ \;
    find . -name *-amrfinder.tab -exec cp {} ./copy_to_stats/ \;
    find . -name *-plasmidfinder.tab -exec cp {} ./copy_to_stats/ \;
    find . -name *_vfdb.txt -exec cp {} ./copy_to_stats/ \;
    find . -name abricate*.xlsx -exec cp {} ./copy_to_stats/ \;
    find . -name amrfinder*.xlsx -exec cp {} ./copy_to_stats/ \;

    cd ./copy_to_stats/
    ncbifile=`echo $DATE"_SIG_ncbi_concat.txt"`
    for i in *-ncbi.tab; do
        sed '1d' ${i} >> $ncbifile
    done
    rm *-ncbi.tab

    resfinderfile=`echo $DATE"_SIG_resfinder_concat.txt"`  
    for i in *-resfinder.tab; do
        sed '1d' ${i} >> $resfinderfile
    done
    rm *-resfinder.tab

    plasmidfinderfile=`echo $DATE"_SIG_plasmidfinder_concat.txt"`
    for i in *-plasmidfinder.tab; do
        sed '1d' ${i} >> $plasmidfinderfile
    done   
    rm *-plasmidfinder.tab        

    amrfinderfile=`echo $DATE"_SIG_amrfinder_concat.txt"`
    for i in *-amrfinder.tab; do
        sed '1d' ${i} >> $amrfinderfile
    done
    rm *-amrfinder.tab

    vfdbfile=`echo $DATE"_SIG_vfdb_concat.txt"`
    for i in *_vfdb.txt; do
        sed '1d' ${i} >> $vfdbfile
    done
    rm *_vfdb.txt
fi


