#!/usr/bin/env python3
"""
This script aggregates all pre-parsed results from the AMR pipeline,
performs final quality scaling calculations, and generates the
final PDF and summary Excel reports for a single sample.
It is designed to be called by the final Nextflow process.
"""

import argparse
import json
import os
import sys
import pandas as pd
from latex_reporter import AMR_Latex_Report
from quality_scaling import Quality_Scaling

def format_file_size(size_bytes):
    """Format file size in human-readable format."""
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if size_bytes < 1024.0:
            return f"{size_bytes:.1f}{unit}"
        size_bytes /= 1024.0
    return f"{size_bytes:.1f}PB"

def parse_seqsero_data(seqsero_json_path):
    """
    Parse SeqSero2 JSON output.
    
    Expected JSON structure from the SEQSERO_PARSER process:
    {
        "serotype": "Enteritidis",
        "antigenic_profile": "9:g,m:-",
        "subspecies": "Salmonella enterica subspecies enterica (subspecies I)",
        "note": "Detected Sdf..."
    }
    
    Returns dict with keys: serotype, antigenic_profile, subspecies, comment
    """
    try:
        with open(seqsero_json_path, 'r') as f:
            data = json.load(f)
        
        # Ensure we have the expected keys with defaults
        result = {
            'serotype': data.get('serotype', 'N/A') or 'N/A',
            'antigenic_profile': data.get('antigenic_profile', 'N/A') or 'N/A',
            'subspecies': data.get('subspecies', 'N/A') or 'N/A',
            'comment': data.get('note', '') or ''
        }
        
        print(f"DEBUG: Parsed SeqSero data successfully:")
        print(f"  - Serotype: {result['serotype']}")
        print(f"  - Antigenic profile: {result['antigenic_profile']}")
        print(f"  - Subspecies: {result['subspecies']}")
        
        return result
        
    except Exception as e:
        print(f"Warning: Error parsing SeqSero JSON: {e}")
        return {
            'serotype': None,
            'antigenic_profile': None,
            'subspecies': None,
            'comment': None
        }

def create_excel_summary(args, fastqc_data, quast_data, coverage_depth, sequencing_type):
    """Creates the final summary Excel file, replicating the original format."""
    print("Creating summary Excel file...")
    sample_name = args.meta_id
    r1_qc = fastqc_data.get('r1', {}) or {}
    r2_qc = fastqc_data.get('r2')  # This can be None for single-end

    columns = [
        'Input Type',
        'Read 1', 'R1 File Size', 'R1 Total Reads', 'R1 Mean Length', 'R1 Mean Quality', 'R1 Passing Q30',
        'Read 2', 'R2 File Size', 'R2 Total Reads', 'R2 Mean Length', 'R2 Mean Quality', 'R2 Passing Q30',
        'Assembly Contig Count', '<300bp Count', '301-999bp Count', '>1kb Count', 'Total Length', 'Longest Contig', 'N50', 'Mean Coverage ((read counts * mean read length) / total assembly length)'
    ]
    df = pd.DataFrame(index=[sample_name], columns=columns)
    
    # --- Set Input Type ---
    if sequencing_type == "Illumina Paired-End (e.g., MiSeq 2x250)":
        df.at[sample_name, 'Input Type'] = 'PE'
    elif sequencing_type == "Illumina Single-End":
        df.at[sample_name, 'Input Type'] = 'SE'
    else:
        df.at[sample_name, 'Input Type'] = 'FASTA'
    
    # --- Get original FASTQ filenames and sizes from args.reads ---
    r1_filename = 'N/A'
    r1_file_size = 'N/A'
    if args.reads:
        r1_filename = os.path.basename(args.reads[0])
        try:
            r1_size_bytes = os.path.getsize(args.reads[0])
            r1_file_size = format_file_size(r1_size_bytes)
        except (OSError, IndexError):
            print(f"Warning: Could not get file size for {args.reads[0]}")
    
    r2_filename = None
    r2_file_size = None
    if len(args.reads) > 1:
        r2_filename = os.path.basename(args.reads[1])
        try:
            r2_size_bytes = os.path.getsize(args.reads[1])
            r2_file_size = format_file_size(r2_size_bytes)
        except (OSError, IndexError):
            print(f"Warning: Could not get file size for {args.reads[1]}")
    
    # --- Populate Read 1 Stats ---
    df.at[sample_name, 'Read 1'] = r1_filename  # Use original FASTQ filename
    df.at[sample_name, 'R1 File Size'] = r1_file_size  # Use calculated file size
    df.at[sample_name, 'R1 Total Reads'] = f"{r1_qc.get('total_sequences', 0):,}"
    df.at[sample_name, 'R1 Mean Length'] = f"{r1_qc.get('mean_length', 0.0):.1f}"
    df.at[sample_name, 'R1 Mean Quality'] = f"{r1_qc.get('mean_quality', 0.0):.1f}"
    sampling_size = fastqc_data.get('sampling_size', 0)
    if sampling_size > 0:
        q30_percent = (r1_qc.get('reads_gt_q30', 0) / sampling_size) * 100
        df.at[sample_name, 'R1 Passing Q30'] = f"{q30_percent:.1f}%"
    else:
        df.at[sample_name, 'R1 Passing Q30'] = "0.0%"

    # --- Populate Read 2 Stats (if it exists) ---
    if r2_qc:
        df.at[sample_name, 'Read 2'] = r2_filename if r2_filename else 'N/A'  # Use original FASTQ filename
        df.at[sample_name, 'R2 File Size'] = r2_file_size if r2_file_size else 'N/A'  # Use calculated file size
        df.at[sample_name, 'R2 Total Reads'] = f"{r2_qc.get('total_sequences', 0):,}"
        df.at[sample_name, 'R2 Mean Length'] = f"{r2_qc.get('mean_length', 0.0):.1f}"
        df.at[sample_name, 'R2 Mean Quality'] = f"{r2_qc.get('mean_quality', 0.0):.1f}"
        if sampling_size > 0:
            q30_percent_r2 = (r2_qc.get('reads_gt_q30', 0) / sampling_size) * 100
            df.at[sample_name, 'R2 Passing Q30'] = f"{q30_percent_r2:.1f}%"
        else:
            df.at[sample_name, 'R2 Passing Q30'] = "0.0%"
    else:
        for col in ['Read 2', 'R2 File Size', 'R2 Total Reads', 'R2 Mean Length', 'R2 Mean Quality', 'R2 Passing Q30']:
            df.at[sample_name, col] = 'N/A'

    # --- Populate Assembly Stats ---
    df.at[sample_name, 'Assembly Contig Count'] = f"{quast_data.get('# contigs', 0):,}"
    df.at[sample_name, '<300bp Count'] = f"{quast_data.get('contigs_lt_300', 0):,}"
    df.at[sample_name, '301-999bp Count'] = f"{quast_data.get('contigs_300_999', 0):,}"
    df.at[sample_name, '>1kb Count'] = f"{quast_data.get('contigs_ge_1kb', 0):,}"
    df.at[sample_name, 'Total Length'] = f"{quast_data.get('Total length', 0):,}"
    df.at[sample_name, 'Longest Contig'] = f"{quast_data.get('Largest contig', 0):,}"
    df.at[sample_name, 'N50'] = f"{quast_data.get('N50', 0):,}"
    df.at[sample_name, 'Mean Coverage ((read counts * mean read length) / total assembly length)'] = f"{coverage_depth:,.1f}X"

    df.index.name = 'sample'
    output_filename = f"{args.meta_id}_stats.xlsx"
    df.to_excel(output_filename)
    print(f"Summary Excel file created: {output_filename}")

def parse_versions_file(filepath):
    """
    Parses the software_versions.txt file and returns a dictionary of software versions.
    """
    versions = {}
    keywords = {
        '  mlst: ': 'mlst',
        '  spades.py: ': 'spades',
        '  SeqSero2_package.py: ': 'seqsero2',
        '  amrfinder: ': 'amrfinder',
        '  abricate: ': 'abricate'
    }
    try:
        with open(filepath, 'r') as f:
            for line in f:
                for keyword, key_name in keywords.items():
                    if line.lstrip().startswith(keyword.lstrip()):
                        version_string = line.strip().replace(keyword, '').split()[-1]
                        print(f"DEBUG: Found version: {key_name} = {version_string}")
                        versions[key_name] = version_string
    except FileNotFoundError:
        print(f"Warning: Versions file not found at {filepath}. Versions will be 'N/A'.")
    
    print(f"DEBUG: Final versions dictionary created: {versions}")
    return versions

def main():
    INCLUDE_ABRICATE_IN_PDF = False # Set to True to include detailed Abricate results in the PDF report

    try:
        parser = argparse.ArgumentParser(description="Generate final PDF and Excel reports.")
        parser.add_argument('--meta_id', required=True, help="Sample identifier")
        parser.add_argument('--logo', required=True, help="Path to logo image file")
        parser.add_argument('--reads', nargs='*', default=[], help="R1 and optional R2 FASTQ files")
        parser.add_argument('--assembly', required=True, help="Assembled FASTA file")
        parser.add_argument('--fastqc_json', required=False, default=None, help="Parsed FastQC stats JSON (optional)")
        parser.add_argument('--quast_json', required=True, help="Parsed Quast stats JSON")
        parser.add_argument('--mlst_json', required=True, help="Parsed MLST results JSON")
        parser.add_argument('--seqsero_json', required=False, default=None, help="Optional parsed SeqSero2 results JSON")
        parser.add_argument('--abricate_json', required=True, help="Parsed Abricate stats JSON")
        parser.add_argument('--amrfinder_tab', required=True, help="Raw AMRFinder output tab file")
        parser.add_argument('--amrfinder_json', required=True, help="Parsed AMRFinder stats JSON")
        parser.add_argument('--versions', required=True, help="Software versions text file")
        parser.add_argument('--ab_ncbi_file', required=True, help="Raw Abricate NCBI tab file")
        parser.add_argument('--ab_resfinder_file', required=True, help="Raw Abricate ResFinder tab file")
        parser.add_argument('--bracken_pie', required=False, default=None, help="Optional Bracken pie chart PNG")
        args = parser.parse_args()

        # debug info for versions.txt content
        print("--- DEBUG: START OF versions.txt CONTENT ---")
        try:
            with open(args.versions, 'r') as f:
                print(f.read())
        except FileNotFoundError:
            print("--- DEBUG: versions.txt FILE NOT FOUND ---")
        print("--- DEBUG: END OF versions.txt CONTENT ---")

        software_versions = parse_versions_file(args.versions)

        fastqc_data = {}
        if args.fastqc_json and os.path.exists(args.fastqc_json):
            try:
                with open(args.fastqc_json, 'r') as f:
                    fastqc_data = json.load(f)
            except (json.JSONDecodeError, IOError):
                print("Warning: FastQC JSON not available or invalid. Proceeding without FASTQ QC data.")
        else:
            print("No FastQC data provided (FASTA-only workflow). Proceeding without FASTQ QC data.")
        
        with open(args.quast_json, 'r') as f:
            quast_data = json.load(f)
        with open(args.mlst_json, 'r') as f:
            mlst_data = json.load(f)
        with open(args.abricate_json, 'r') as f:
            abricate_data = json.load(f)
        with open(args.amrfinder_json, 'r') as f:
            amrfinder_data = json.load(f)
        
        seqsero_data = {}
        if args.seqsero_json and os.path.exists(args.seqsero_json):
            seqsero_data = parse_seqsero_data(args.seqsero_json)
        else:
            print("No SeqSero data provided (non-Salmonella sample).")
            seqsero_data = {
                'serotype': None,
                'antigenic_profile': None,
                'subspecies': None,
                'comment': None
            }

        with open(args.versions, 'r') as f:
            versions_text = f.read()

        # --- 2. Perform final scaling calculations ---
        print("Calculating quality scaling variables...")
        scaler = Quality_Scaling()
        
        try:
            genome_size = int(float(mlst_data.get('size_lookup', 0)) * 1000000)
            size_method = "Based on MLST identification"
            if genome_size == 0:
                raise ValueError
        except (ValueError, TypeError, KeyError):
            genome_size = quast_data.get('Total length', 0)
            size_method = "SPAdes total contig length"
        if genome_size == 0:
            genome_size = 1

        # Initialize variables
        r1_qc = fastqc_data.get('r1')
        r2_qc = fastqc_data.get('r2')
        
        fastq_scaling_variable = 1180
        genome_coverage_depth = quast_data.get('Mean coverage', 0)
        coverage_method = "From Assembly (SPAdes)"

        # Perform FASTQ scaling if we have R1 data
        if r1_qc:
            scaling_result = scaler.fastq_scaling(
                genome_size=genome_size,
                read1_reads_gt_q30=r1_qc.get('reads_gt_q30', 0),
                read2_reads_gt_q30=r2_qc.get('reads_gt_q30', 0) if r2_qc else 0,
                sampling_size=fastqc_data.get('sampling_size', 0),
                read1_total_read_count=r1_qc.get('total_sequences', 0),
                read2_total_read_count=r2_qc.get('total_sequences', 0) if r2_qc else 0,
                read1_read_average=r1_qc.get('mean_quality', 0.0),
                read2_read_average=r2_qc.get('mean_quality', 0.0) if r2_qc else 0.0,
                read1_length_mean=r1_qc.get('mean_length', 0.0),
                read2_length_mean=r2_qc.get('mean_length', 0.0) if r2_qc else 0.0
            )
            if scaling_result:
                fastq_scaling_variable, genome_coverage_depth = scaling_result
                coverage_method = "Calculated by (read counts * mean length) / Genome Size"

        # Perform assembly scaling
        assembly_scaling_variable = scaler.assembly_scaling(
            longest_contig=quast_data.get('Largest contig', 0),
            greater_one_kb_count=quast_data.get('# contigs (>= 1000 bp)', 0),
            contig_count=quast_data.get('# contigs', 0),
            n50=quast_data.get('N50', 0),
            l50=quast_data.get('L50', 0),
            rgl=quast_data.get('rgl', 0.0),
            stat_total_contig_lengths=quast_data.get('Total length', 0),
            stat_contig_count=quast_data.get('# contigs', 0),
            genome_coverage_depth=genome_coverage_depth
        )

        # Determine sequencing type based on FastQC data availability
        if r1_qc and r2_qc:
            sequencing_type = "Illumina Paired-End (e.g., MiSeq 2x250)"
        elif r1_qc:
            sequencing_type = "Illumina Single-End"
        else:
            sequencing_type = "Assembly-Only (FASTA input)"
        
        print(f"DEBUG: sequencing_type = {sequencing_type}")
        print(f"DEBUG: r1_qc exists: {r1_qc is not None}")
        print(f"DEBUG: r2_qc exists: {r2_qc is not None}")

        # --- Get file sizes for LaTeX report ---
        r1_file_size = 'N/A'
        r2_file_size = None
        if args.reads:
            try:
                r1_size_bytes = os.path.getsize(args.reads[0])
                r1_file_size = format_file_size(r1_size_bytes)
            except (OSError, IndexError):
                pass
        if len(args.reads) > 1:
            try:
                r2_size_bytes = os.path.getsize(args.reads[1])
                r2_file_size = format_file_size(r2_size_bytes)
            except (OSError, IndexError):
                pass

        # --- 3. Instantiate and run the LaTeX Report Generator ---
        print("Building final PDF report...")
        amr_latex_report = AMR_Latex_Report(
            fastq_scaling_variable=fastq_scaling_variable,
            assembly_scaling_variable=assembly_scaling_variable,
            genome_size=genome_size,
            genome_coverage_depth=genome_coverage_depth,
            coverage_method=coverage_method,
            size_method=size_method,
            abricate_ab_version=abricate_data.get('abricate_version', 'N/A'),
            amr_version=amrfinder_data.get('amrfinder_version', 'N/A'),
            read1_fastq=args.reads[0] if args.reads else None,
            read2_fastq=args.reads[1] if len(args.reads) > 1 else None,
            logo_path=args.logo
        )

        print(f"DEBUG: Value for spades_version being passed: {software_versions.get('spades', 'N/A')}")
        print(f"DEBUG: Value for mlst_version being passed: {software_versions.get('mlst', 'N/A')}")

        print("DEBUG: About to call latex_document()...")
        amr_latex_report.latex_document(
            sample_name=args.meta_id,
            sequencing_type=sequencing_type,
            read1_fastq=args.reads[0] if args.reads else None,
            read2_fastq=args.reads[1] if len(args.reads) > 1 else None,
            read1_file_size=r1_file_size,  # Use calculated file size
            read2_file_size=r2_file_size,  # Use calculated file size
            read1_read_average=r1_qc.get('mean_quality', 0.0) if r1_qc else 0.0,
            read2_read_average=r2_qc.get('mean_quality', 0.0) if r2_qc else None,
            read1_read_length=r1_qc.get('mean_length', 0.0) if r1_qc else 0.0,
            read2_read_length=r2_qc.get('mean_length', 0.0) if r2_qc else None,
            read1_reads_gt_q30=r1_qc.get('reads_gt_q30', 0) if r1_qc else 0,
            read2_reads_gt_q30=r2_qc.get('reads_gt_q30', 0) if r2_qc else 0,
            sampling_size=fastqc_data.get('sampling_size', 0),
            stat_contig_count=quast_data.get('# contigs', 0),
            stat_total_contig_lengths=quast_data.get('Total length', 0),
            stat_longest_contig=quast_data.get('Largest contig', 0),
            stat_greater_one_kb_count=quast_data.get('# contigs (>= 1000 bp)', 0),
            stat_n50=quast_data.get('N50', 0),
            stat_l50=quast_data.get('L50', 0),
            rgl=quast_data.get('rgl', 0.0),
            spades_version=software_versions.get('spades', 'N/A'),
            mlst_file=args.mlst_json,
            mlst_scheme=mlst_data.get('scheme', 'N/A'),
            mlst_st=mlst_data.get('type', 'N/A'),
            mlst_detail=mlst_data.get('detail', []),
            mlst_species_lookup=mlst_data.get('species_lookup', 'N/A'),
            mlst_version=software_versions.get('mlst', 'N/A'),
            seqsero2_serotype=seqsero_data.get('serotype', None),
            seqsero2_antigenic=seqsero_data.get('antigenic_profile', None),
            seqsero2_subspecies=seqsero_data.get('subspecies', None),
            seqserocomment=seqsero_data.get('comment', None),
            seqsero_file=args.seqsero_json,
            amrfinder_file=args.amrfinder_tab,
            abricate_mincov=abricate_data.get('mincov', 0),      # This is "depth" in legacy naming
            abricate_minid=abricate_data.get('minid', 75),       # This is "coverage" in legacy naming
            ab_ncbi_file=args.ab_ncbi_file,
            ab_resfinder_file=args.ab_resfinder_file,
            abricate_ncbi_version_date=abricate_data.get('db_info', {}).get('ncbi', {}).get('date', 'N/A'),
            abricate_res_version_date=abricate_data.get('db_info', {}).get('resfinder', {}).get('date', 'N/A'),
            abricate_ncbi_seq_number=abricate_data.get('db_info', {}).get('ncbi', {}).get('sequences', 0),
            abricate_res_seq_number=abricate_data.get('db_info', {}).get('resfinder', {}).get('sequences', 0),
            software_versions=versions_text,
            bracken_pie_file=args.bracken_pie,
            logo_path=args.logo,
            abricate_report=INCLUDE_ABRICATE_IN_PDF
        )
        print("DEBUG: latex_document() completed successfully")

        # --- 4. Create the final Excel summary file ---
        create_excel_summary(args, fastqc_data, quast_data, genome_coverage_depth, sequencing_type)

        print("Final report generation complete.")
        
    except Exception as e:
        print(f"FATAL ERROR: {type(e).__name__}: {str(e)}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()