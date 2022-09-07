#!/usr/bin/env bash

ts=`date "+%Y-%m-%d"`
#ab_res_old=`abricate --list|sed '7!d'|cut -f 3`
ab_res_old=`abricate --list | awk 'BEGIN{OFS="\t"} $1 == "resfinder" {print $4}'`
#ab_res_seqs_old=`abricate --list|sed '7!d'|cut -f 2`
ab_res_seqs_old=`abricate --list | awk 'BEGIN{OFS="\t"} $1 == "resfinder" {print $2}'`
#ab_ncbi_old=`abricate --list|sed '5!d'|cut -f 3`
ab_ncbi_old=`abricate --list | awk 'BEGIN{OFS="\t"} $1 == "ncbi" {print $4}'`
#ab_ncbi_seqs_old=`abricate --list|sed '5!d'|cut -f 2`
ab_ncbi_seqs_old=`abricate --list | awk 'BEGIN{OFS="\t"} $1 == "ncbi" {print $2}'`
echo "Res_old = '${ab_res_old}, Res_Seqs = ${ab_res_seqs_old}"

anaconda_dbPath=/usr/local/bin/anaconda3/db

mkdir ./old/

#mkdir /scratch/analysis/db_update/old
#cd /scratch/analysis/db_update/
cp -r $anaconda_dbPath/resfinder/ ./old/
cp -r $anaconda_dbPath/ncbi/ ./old/
echo "Resfinder databases copied"

grep "^>" ./old/resfinder/sequences  >> resfinder_headers_old
grep "^>" ./old/ncbi/sequences >> ncbi_headers_old
echo "Old Headers found."

abricate-get_db --db resfinder --force
abricate-get_db --db ncbi --force

mkdir ./new

#ab_res_new=`abricate --list|sed '7!d'|cut -f 3`
ab_res_new=`abricate --list | awk 'BEGIN{OFS="\t"} $1 == "resfinder" {print $4}'`
#ab_res_seqs_new=`abricate --list|sed '7!d'|cut -f 2`
ab_res_seqs_new=`abricate --list | awk 'BEGIN{OFS="\t"} $1 == "resfinder" {print $2}'`
#ab_ncbi_new=`abricate --list|sed '5!d'|cut -f 3`
ab_ncbi_new=`abricate --list | awk 'BEGIN{OFS="\t"} $1 == "ncbi" {print $4}'`
#ab_ncbi_seqs_new=`abricate --list|sed '5!d'|cut -f 2`
ab_ncbi_seqs_new=`abricate --list | awk 'BEGIN{OFS="\t"} $1 == "ncbi" {print $2}'`
echo "Res_new = ${ab_res_new}, Res_seqs = ${ab_res_seqs_new}"
echo "NCBI_new = ${ab_ncbi_new}, NCBI_seqs = ${ab_ncbi_seqs_new}"

cp -r $anaconda_dbPath/resfinder/ ./new/
cp -r $anaconda_dbPath/ncbi/ ./new/

grep "^>" ./new/resfinder/sequences  >> resfinder_headers_new
grep "^>" ./new/ncbi/sequences >> ncbi_headers_new
echo "New Headers found."

grep -Fxvf resfinder_headers_new resfinder_headers_old >> res_seqs_removed
grep -Fxvf resfinder_headers_old resfinder_headers_new >> res_seqs_added
res_removed=`cat res_seqs_removed |wc -l`
res_added=`cat res_seqs_added | wc -l`

grep -Fxvf ncbi_headers_new ncbi_headers_old >> ncbi_seqs_removed
grep -Fxvf ncbi_headers_old ncbi_headers_new >> ncbi_seqs_added
ncbi_removed=`cat ncbi_seqs_removed |wc -l`
ncbi_added=`cat ncbi_seqs_added | wc -l`

echo "################# Version	#####Sequences" >> abricate-resfinder_change_log_${ts}.txt
echo "Resfinder Previous: ${ab_res_old}	##### ${ab_res_seqs_old}" >> abricate-resfinder_change_log_${ts}.txt
echo "Resfinder New Version: ${ab_res_new}	##### ${ab_res_seqs_new}" >> abricate-resfinder_change_log_${ts}.txt
echo "


 ####Sequences Removed: ${res_removed}####

" >> abricate-resfinder_change_log_${ts}.txt
cat res_seqs_removed >> abricate-resfinder_change_log_${ts}.txt
echo "


 ####Sequences Added: ${res_added}####

" >> abricate-resfinder_change_log_${ts}.txt
cat res_seqs_added >> abricate-resfinder_change_log_${ts}.txt


echo "################# Version #####Sequences" >> abricate-ncbi_change_log_${ts}.txt
echo "NCBI Previous: ${ab_ncbi_old}	##### ${ab_ncbi_seqs_old}" >> abricate-ncbi_change_log_${ts}.txt
echo "NCBI New: ${ab_ncbi_new}	##### ${ab_ncbi_seqs_new}" >> abricate-ncbi_change_log_${ts}.txt
echo "


 ####Sequences Removed: ${ncbi_removed}####

" >> abricate-ncbi_change_log_${ts}.txt
cat ncbi_seqs_removed >> abricate-ncbi_change_log_${ts}.txt
echo "


 ####Sequences Added: ${ncbi_added}####

" >> abricate-ncbi_change_log_${ts}.txt
cat ncbi_seqs_added >> abricate-ncbi_change_log_${ts}.txt


