log.info """\
 A M R   P I P E L I N E
 ===================================
 kraken_db    : ${params.kraken_db}
 outdir       : ${params.outdir}
 """

process FASTQC {
    tag "$meta.id"
    label 'qc_process'

    publishDir "${params.outdir}/${meta.id}/fastqc", mode: 'copy', pattern: "*.{zip,html,txt}"
    publishDir "${params.outdir}/${meta.id}", mode: 'copy', pattern: "fastq/*.{fastq,fastq.gz,fq,fq.gz}"

    input:
    tuple val(meta), path(read)

    output:
    tuple val(meta), path("*.zip"),                         emit: zip
    tuple val(meta), path("*.html"),                        emit: html
    tuple val(meta), path("*_fastqc_data.txt"),             emit: data_txt
    path "fastq/*.{fastq,fastq.gz,fq,fq.gz}",               emit: reads_out

    script:
    def prefix = read.simpleName
    """
    # Create subdirectory and copy input reads there for publishing
    mkdir -p fastq
    cp -L ${read} fastq/
    
    fastqc \\
        -o . \\
        -t ${task.cpus} \\
        ${read}

    unzip -q ${prefix}_fastqc.zip
    mv ${prefix}_fastqc/fastqc_data.txt ${prefix}_fastqc_data.txt
    """
}

process KRAKEN_PE {
    tag "$meta"
    label 'kraken_process'

    publishDir "${params.outdir}/${meta}/kraken2", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta}_report.txt"),         emit: report
    tuple val(meta), path("${meta}_output.txt"),         emit: output

    script:
    """
    kraken2 \\
        --db ${params.kraken_db} \\
        --threads ${task.cpus} \\
        --paired ${reads[0]} ${reads[1]} \\
        --report "${meta}_report.txt" \\
        --output "${meta}_output.txt"
    """
}

process KRAKEN_SE {
    tag "$meta"
    label 'kraken_process'

    publishDir "${params.outdir}/${meta}/kraken2", mode: 'copy'

    input:
    tuple val(meta), path(read)

    output:
    tuple val(meta), path("${meta}_report.txt"),         emit: report
    tuple val(meta), path("${meta}_output.txt"),         emit: output

    script:
    """
    kraken2 --db ${params.kraken_db} --threads ${task.cpus} ${read} --report "${meta}_report.txt" --output "${meta}_output.txt"
    """
}

process KRAKEN_FA {
    tag "$meta"
    label 'kraken_process'

    publishDir "${params.outdir}/${meta}/kraken2", mode: 'copy'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta}_report.txt"),         emit: report
    tuple val(meta), path("${meta}_output.txt"),         emit: output

    script:
    """
    kraken2 --db ${params.kraken_db} --threads ${task.cpus} ${fasta} --report "${meta}_report.txt" --output "${meta}_output.txt"
    """
}

process BRACKEN {
    tag "$meta"
    label 'bracken_process'

    publishDir "${params.outdir}/${meta}/kraken2", mode: 'copy'

    input:
    tuple val(meta), path(kraken_report)

    output:
    tuple val(meta), path("${meta}.bracken.xlsx"),       emit: bracken_excel

    script:
    """
    # Use a 'try/catch' block in bash to handle bracken failures gracefully
    if ! bracken -d "${params.kraken_db}" -i ${kraken_report} -o "${meta}.bracken.txt" -r 250; then
        echo "Bracken failed or produced no output. This is normal for low-yield samples. Creating a placeholder Excel file."
        
        echo "name,taxonomy_id,taxonomy_lvl,kraken_assigned_reads,added_reads,new_est_reads,fraction_total_reads" > placeholder.csv
        echo "No Bracken results - insufficient data,N/A,N/A,0,0,0,0.0" >> placeholder.csv
        
        python -c "import pandas as pd; pd.read_csv('placeholder.csv').to_excel('${meta}.bracken.xlsx', index=False)"
    else
        echo "Bracken completed successfully. Converting output to Excel."
        python -c "import pandas as pd; pd.read_csv('${meta}.bracken.txt', sep='\\t').to_excel('${meta}.bracken.xlsx', index=False)"
    fi
    """
}

process KRONA {
    tag "$meta"
    label 'krona_process'

    publishDir "${params.outdir}/${meta}/kraken2", mode: 'copy'

    input:
    tuple val(meta), path(kraken_report)

    output:
    tuple val(meta), path("${meta}.krona.html"),         emit: krona_html

    script:
    """
    kreport2krona.py -r ${kraken_report} -o "${meta}.krona.txt"

    ktImportText "${meta}.krona.txt" -o "${meta}.krona.html"
    """
}

process BRACKEN_PIE {
    tag "$meta"
    label 'bracken_pie_process'

    publishDir "${params.outdir}/${meta}/kraken2", mode: 'copy'

    input:
    tuple val(meta), path(bracken_excel)

    output:
    tuple val(meta), path("${meta}_bracken_pie.png"),    emit: pie_chart

    script:
    """
    bracken_pie.py ${meta} ${bracken_excel}
    """
}

process SPADES_PE {
    tag "$meta"
    label 'assembly_process'

    publishDir "${params.outdir}/${meta}/spades", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fasta"), emit: assembly
    tuple val(meta), path("spades_output/spades.log"), optional: true, emit: log
    tuple val(meta), path("spades_output/warnings.log"), optional: true, emit: warnings

    script:
    """
    spades.py \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        -o spades_output \\
        --threads ${task.cpus}
    
    # Prefer scaffolds.fasta, fall back to contigs.fasta if scaffolds don't exist
    if [ -f "spades_output/scaffolds.fasta" ]; then
        cp spades_output/scaffolds.fasta ${meta}.fasta
    elif [ -f "spades_output/contigs.fasta" ]; then
        cp spades_output/contigs.fasta ${meta}.fasta
    else
        echo "ERROR: Neither scaffolds.fasta nor contigs.fasta found!" >&2
        exit 1
    fi
    """
}

process SPADES_SE {
    tag "$meta"
    label 'assembly_process'

    publishDir "${params.outdir}/${meta}/spades", mode: 'copy'

    input:
    tuple val(meta), path(read)

    output:
    tuple val(meta), path("*.fasta"), emit: assembly
    tuple val(meta), path("spades_output/spades.log"), optional: true, emit: log
    tuple val(meta), path("spades_output/warnings.log"), optional: true, emit: warnings

    script:
    """
    spades.py \\
        -s ${read} \\
        -o spades_output \\
        --threads ${task.cpus}
    
    # Prefer scaffolds.fasta, fall back to contigs.fasta if scaffolds don't exist
    if [ -f "spades_output/scaffolds.fasta" ]; then
        cp spades_output/scaffolds.fasta ${meta}.fasta
    elif [ -f "spades_output/contigs.fasta" ]; then
        cp spades_output/contigs.fasta ${meta}.fasta
    else
        echo "ERROR: Neither scaffolds.fasta nor contigs.fasta found!" >&2
        exit 1
    fi
    """
}

process QUAST {
    tag "$meta"
    label 'qc_process'

    publishDir "${params.outdir}/${meta}/spades", mode: 'copy'

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("quast_output/report.tsv"),       emit: report_tsv
    path "quast_output",                                    emit: report_dir

    script:
    """
    quast.py \\
        -o quast_output \\
        -t ${task.cpus} \\
        --min-contig 0 \\
        ${assembly}
    """
}

process SEQSERO_PE {
    tag "$meta"
    label 'seqsero_process'
    
    publishDir "${params.outdir}/${meta}/seqsero", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("seqsero/SeqSero_result.*"),      emit: result_file, optional: true

    script:
    """
    SeqSero2_package.py -d seqsero -t 2 -p ${task.cpus} -m a -i ${reads[0]} ${reads[1]} || true

    if [ ! -f "seqsero/SeqSero_result.txt" ] && [ ! -f "seqsero/SeqSero_result.tsv" ]; then
        echo "Allele mode failed or produced no results. Retrying with k-mer mode..."
        rm -rf seqsero
        SeqSero2_package.py -d seqsero -t 2 -p ${task.cpus} -m k -i ${reads[0]} ${reads[1]} || true

    else
        echo "Allele mode completed."
    fi

    if [ ! -f "seqsero/SeqSero_result.txt" ] && [ ! -f "seqsero/SeqSero_result.tsv" ]; then
        echo "Both SeqSero2 modes failed to produce output. Creating an empty placeholder result file."
        mkdir -p seqsero
        touch seqsero/SeqSero_result.txt
    fi
    """
}

process SEQSERO_SE {
    tag "$meta" 
    label 'seqsero_process'
    
    publishDir "${params.outdir}/${meta}/seqsero", mode: 'copy'

    input:
    tuple val(meta), path(read)

    output:
    tuple val(meta), path("seqsero/SeqSero_result.*"),          emit: result_file, optional: true

    script:
    """
    SeqSero2_package.py -d seqsero -t 3 -p ${task.cpus} -m a -i ${read} || true

    if [ ! -f "seqsero/SeqSero_result.txt" ] && [ ! -f "seqsero/SeqSero_result.tsv" ]; then
        echo "Allele mode failed or produced no results. Retrying with k-mer mode..."
        rm -rf seqsero
        SeqSero2_package.py -d seqsero -t 3 -p ${task.cpus} -m k -i ${read} || true

    else
        echo "Allele mode completed."
    fi

    if [ ! -f "seqsero/SeqSero_result.txt" ] && [ ! -f "seqsero/SeqSero_result.tsv" ]; then
        echo "Both SeqSero2 modes failed to produce output. Creating an empty placeholder result file."
        mkdir -p seqsero
        touch seqsero/SeqSero_result.txt
    fi
    """
}

process SEQSERO_FA {
    tag "$meta"
    label 'seqsero_process'
    
    publishDir "${params.outdir}/${meta}/seqsero", mode: 'copy'

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("seqsero/SeqSero_result.*"),          emit: result_file, optional: true

    script:
    """
    # SeqSero2 assembly mode (-t 4)
    SeqSero2_package.py -d seqsero -t 4 -p ${task.cpus} -m a -i ${assembly} || true

    if [ ! -f "seqsero/SeqSero_result.txt" ] && [ ! -f "seqsero/SeqSero_result.tsv" ]; then
        echo "Allele mode failed or produced no results. Retrying with k-mer mode..."
        rm -rf seqsero
        SeqSero2_package.py -d seqsero -t 4 -p ${task.cpus} -m k -i ${assembly} || true
    else
        echo "Allele mode completed."
    fi

    if [ ! -f "seqsero/SeqSero_result.txt" ] && [ ! -f "seqsero/SeqSero_result.tsv" ]; then
        echo "Both SeqSero2 modes failed to produce output. Creating empty placeholder result file."
        mkdir -p seqsero
        touch seqsero/SeqSero_result.txt
    fi
    """

    /* SeqSero input FYI:
        -t <data_type>
        1: Interleaved paired-end
        2: Separated paired-end  ← Your SEQSERO_PE uses this
        3: Single-end            ← Your SEQSERO_SE uses this
        4: Assembly (FASTA)      ← Use this for SEQSERO_FASTA
        5: Nanopore reads
        6: PacBio reads
    */
}

process MLST {
    tag "$meta"
    label 'mlst_process'
    
    publishDir "${params.outdir}/${meta}/mlst", mode: 'copy'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta}_mlst.txt"),           emit: mlst_txt

    script:
    """
    mlst --threads ${task.cpus} ${fasta} > "${meta}_mlst.txt"
    """
}

process ABRICATE {
    tag "$meta"
    label 'abricate_process'
    
    publishDir "${params.outdir}/${meta}/abricate", mode: 'copy'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta}-resfinder.tab")   ,   emit: resfinder_tab
    tuple val(meta), path("${meta}-ncbi.tab")        ,   emit: ncbi_tab
    tuple val(meta), path("${meta}-plasmidfinder.tab"),  emit: plasmidfinder_tab

    script:
    """
    WORKDIR=\$(pwd)
    cp ${fasta} \$WORKDIR/resfinder_input.fasta
    cp ${fasta} \$WORKDIR/ncbi_input.fasta
    cp ${fasta} \$WORKDIR/plasmidfinder_input.fasta
    chmod a+r \$WORKDIR/resfinder_input.fasta \$WORKDIR/ncbi_input.fasta \$WORKDIR/plasmidfinder_input.fasta

    abricate --db resfinder     --mincov ${params.abricate_depth} --minid ${params.abricate_coverage} --threads ${task.cpus} \$WORKDIR/resfinder_input.fasta    > "\$WORKDIR/${meta}-resfinder.tab"
    abricate --db ncbi          --mincov ${params.abricate_depth} --minid ${params.abricate_coverage} --threads ${task.cpus} \$WORKDIR/ncbi_input.fasta         > "\$WORKDIR/${meta}-ncbi.tab"
    abricate --db plasmidfinder --mincov ${params.abricate_depth} --minid ${params.abricate_coverage} --threads ${task.cpus} \$WORKDIR/plasmidfinder_input.fasta > "\$WORKDIR/${meta}-plasmidfinder.tab"
    """
}

process AMRFINDER {
    tag "$meta"
    label 'amrfinder_process'
    
    publishDir "${params.outdir}/${meta}/amrfinder", mode: 'copy'

    input:
    tuple val(meta), path(fasta), path(bracken_excel)

    output:
    tuple val(meta), path("${meta}-amrfinder.tab"),      emit: amrfinder_tab

    script:
    """
    #!/bin/bash
    set -e

    DB_FLAG=""
    if [[ -n "\$AMRFINDERPLUS_DB" ]]; then
        DB_FLAG="--database \$AMRFINDERPLUS_DB"
    elif [[ -d "/opt/conda/envs/analysis_env/share/amrfinderplus/data/latest" ]]; then
        DB_FLAG="--database /opt/conda/envs/analysis_env/share/amrfinderplus/data/latest"
    fi

    TOP_ORGANISM=\$(python -c "import pandas as pd; df = pd.read_excel('${bracken_excel}'); print(df['name'].iloc[0].replace(' ', '_'))" || true)
    ORGANISM_FLAG=""
    GENUS=\$(echo "\$TOP_ORGANISM" | cut -d'_' -f1)
    
    # Logic to determine organism flag...
    case "\$TOP_ORGANISM" in
        Acinetobacter_baumannii|...)
            ORGANISM_FLAG="-O \$TOP_ORGANISM" ;;
    esac
    if [ -z "\$ORGANISM_FLAG" ]; then
        case "\$GENUS" in
            Campylobacter|Escherichia|Salmonella) ORGANISM_FLAG="-O \$GENUS" ;;
            Shigella) ORGANISM_FLAG="-O Escherichia" ;;
        esac
    fi
    if [ -z "\$ORGANISM_FLAG" ]; then
        ORGANISM_FLAG="--plus"
    fi

    amrfinder \\
        --nucleotide ${fasta} \\
        --threads ${task.cpus} \\
        \$DB_FLAG \\
        \$ORGANISM_FLAG \\
        --output "${meta}-amrfinder.tab"
    """
}

process AMRFINDER_FA {
    tag "$meta"
    label 'amrfinder_process'
    
    publishDir "${params.outdir}/${meta}/amrfinder", mode: 'copy'

    input:
    tuple val(meta), path(fasta), path(mlst_json)

    output:
    tuple val(meta), path("${meta}-amrfinder.tab"), emit: amrfinder_tab

    script:
    """
    #!/bin/bash
    set -e

    DB_FLAG=""
    if [[ -n "\$AMRFINDERPLUS_DB" ]]; then
        DB_FLAG="--database \$AMRFINDERPLUS_DB"
    elif [[ -d "/opt/conda/envs/analysis_env/share/amrfinderplus/data/latest" ]]; then
        DB_FLAG="--database /opt/conda/envs/analysis_env/share/amrfinderplus/data/latest"
    fi

    # Extract organism from MLST JSON
    ORGANISM_FLAG=""
    if [[ -f "${mlst_json}" && "${mlst_json}" != "NO_FILE" ]]; then
        ORGANISM=\$(python -c "
import json
try:
    with open('${mlst_json}') as f:
        data = json.load(f)
        species = data.get('mlst_species_lookup', '')
        if species:
            print(species.replace(' ', '_'))
except:
    print('', end='')
" || echo "")
        
        if [[ -n "\$ORGANISM" ]]; then
            case "\$ORGANISM" in
                Escherichia_coli|Salmonella_enterica|Campylobacter*)
                    GENUS=\$(echo "\$ORGANISM" | cut -d'_' -f1)
                    ORGANISM_FLAG="-O \$GENUS"
                    echo "Organism from MLST: \$GENUS"
                    ;;
            esac
        fi
    fi
    
    if [[ -z "\$ORGANISM_FLAG" ]]; then
        ORGANISM_FLAG="--plus"
        echo "Using --plus mode"
    fi

    amrfinder \\
        --nucleotide ${fasta} \\
        --threads ${task.cpus} \\
        \$DB_FLAG \\
        \$ORGANISM_FLAG \\
        --output "${meta}-amrfinder.tab"
    """
}

process PARSE_FASTQC {
    tag "$meta"
    label 'parsing_process'

    input:
    tuple val(meta), path(r1_data), path(r2_data)

    output:
    tuple val(meta), path("${meta}_fastqc.json"), emit: fastqc_json

    script:
    """
    #!/usr/bin/env python3
    import json
    import os
    import re
    import humanize

    def parse_fastqc_data(file_path, original_fastq_path):
        stats = {
            'filename': os.path.basename(original_fastq_path),
            'file_size': '0B',
            'total_sequences': 0,
            'mean_length': 0.0,
            'mean_quality': 0.0,
            'reads_gt_q30': 0
        }

        if not file_path or not os.path.exists(file_path):
            return stats

        # --- Parse Basic Stats and Length Distribution ---
        total_length_score = 0
        total_reads_for_length = 0
        in_length_module = False
        
        with open(file_path, 'r') as f:
            for line in f:
                # Get total sequences
                if line.startswith('Total Sequences'):
                    stats['total_sequences'] = int(line.strip().split('\\t')[1])
                
                # Track when we enter the length distribution module
                if '>>Sequence Length Distribution' in line:
                    in_length_module = True
                    continue
                if in_length_module and '>>END_MODULE' in line:
                    in_length_module = False
                    continue
                
                # Parse length distribution data
                if in_length_module and not line.startswith('#'):
                    try:
                        length_range, count = line.strip().split('\\t')
                        count = int(float(count))
                        
                        # Handle length ranges like "35-40" or single values like "151"
                        if '-' in length_range:
                            lengths = [int(x) for x in length_range.split('-')]
                            avg_length = sum(lengths) / len(lengths)
                        else:
                            avg_length = float(length_range)
                        
                        total_length_score += avg_length * count
                        total_reads_for_length += count
                    except (ValueError, IndexError):
                        continue
        
        # Calculate mean length from distribution
        if total_reads_for_length > 0:
            stats['mean_length'] = total_length_score / total_reads_for_length

        # Get file size
        if os.path.exists(original_fastq_path):
            stats['file_size'] = humanize.naturalsize(os.path.getsize(original_fastq_path))

        # --- Calculate Mean Quality and Reads > Q30 ---
        total_quality_score = 0
        total_reads_for_quality = 0
        reads_gt_q30 = 0
        in_quality_module = False
        
        with open(file_path, 'r') as f:
            for line in f:
                if '>>Per sequence quality scores' in line:
                    in_quality_module = True
                    continue
                if '>>END_MODULE' in line:
                    in_quality_module = False
                    continue
                if in_quality_module and not line.startswith('#'):
                    try:
                        quality, count = line.strip().split('\\t')
                        quality, count = float(quality), int(float(count))
                        total_quality_score += quality * count
                        total_reads_for_quality += count
                        if quality >= 30:
                            reads_gt_q30 += count
                    except (ValueError, IndexError):
                        continue
        
        if total_reads_for_quality > 0:
            stats['mean_quality'] = total_quality_score / total_reads_for_quality
        
        stats['reads_gt_q30'] = reads_gt_q30
        
        return stats

    # --- Main Execution Logic ---
    r1_original_fastq_zip = '${r1_data}'.replace('_fastqc/fastqc_data.txt', '.zip')
    r1_stats = parse_fastqc_data('${r1_data}', r1_original_fastq_zip)

    r2_stats = None
    if '${r2_data}' != 'null' and os.path.exists('${r2_data}'):
        r2_original_fastq_zip = '${r2_data}'.replace('_fastqc/fastqc_data.txt', '.zip')
        r2_stats = parse_fastqc_data('${r2_data}', r2_original_fastq_zip)

    sampling_size = r1_stats['total_sequences'] if r1_stats else 0

    final_data = {
        'r1': r1_stats,
        'r2': r2_stats,
        'sampling_size': sampling_size
    }

    with open('${meta}_fastqc.json', 'w') as f:
        json.dump(final_data, f, indent=4)

    print(f"Successfully parsed FastQC data for ${meta}")
    """
}

process PARSE_FASTQC_SE {
    tag "$meta"
    label 'parsing_process'

    input:
    tuple val(meta), path(data_txt)

    output:
    tuple val(meta), path("${meta}_fastqc.json"), emit: fastqc_json

    script:
    """
    #!/usr/bin/env python3
    import json
    import os
    import re
    import humanize

    def parse_fastqc_data(file_path, original_fastq_path):
        stats = {
            'filename': os.path.basename(original_fastq_path),
            'file_size': '0B',
            'total_sequences': 0,
            'mean_length': 0.0,
            'mean_quality': 0.0,
            'reads_gt_q30': 0
        }

        if not file_path or not os.path.exists(file_path):
            return stats

        # --- Parse Basic Stats and Length Distribution ---
        total_length_score = 0
        total_reads_for_length = 0
        in_length_module = False
        
        with open(file_path, 'r') as f:
            for line in f:
                # Get total sequences
                if line.startswith('Total Sequences'):
                    stats['total_sequences'] = int(line.strip().split('\\t')[1])
                
                # Track when we enter the length distribution module
                if '>>Sequence Length Distribution' in line:
                    in_length_module = True
                    continue
                if in_length_module and '>>END_MODULE' in line:
                    in_length_module = False
                    continue
                
                # Parse length distribution data
                if in_length_module and not line.startswith('#'):
                    try:
                        length_range, count = line.strip().split('\\t')
                        count = int(float(count))
                        
                        # Handle length ranges like "35-40" or single values like "151"
                        if '-' in length_range:
                            lengths = [int(x) for x in length_range.split('-')]
                            avg_length = sum(lengths) / len(lengths)
                        else:
                            avg_length = float(length_range)
                        
                        total_length_score += avg_length * count
                        total_reads_for_length += count
                    except (ValueError, IndexError):
                        continue
        
        # Calculate mean length from distribution
        if total_reads_for_length > 0:
            stats['mean_length'] = total_length_score / total_reads_for_length

        # Get file size
        if os.path.exists(original_fastq_path):
            stats['file_size'] = humanize.naturalsize(os.path.getsize(original_fastq_path))

        # --- Calculate Mean Quality and Reads > Q30 ---
        total_quality_score = 0
        total_reads_for_quality = 0
        reads_gt_q30 = 0
        in_quality_module = False
        
        with open(file_path, 'r') as f:
            for line in f:
                if '>>Per sequence quality scores' in line:
                    in_quality_module = True
                    continue
                if '>>END_MODULE' in line:
                    in_quality_module = False
                    continue
                if in_quality_module and not line.startswith('#'):
                    try:
                        quality, count = line.strip().split('\\t')
                        quality, count = float(quality), int(float(count))
                        total_quality_score += quality * count
                        total_reads_for_quality += count
                        if quality >= 30:
                            reads_gt_q30 += count
                    except (ValueError, IndexError):
                        continue
        
        if total_reads_for_quality > 0:
            stats['mean_quality'] = total_quality_score / total_reads_for_quality
        
        stats['reads_gt_q30'] = reads_gt_q30
        
        return stats

    # --- Main Execution Logic for Single-End ---
    original_fastq_zip = '${data_txt}'.replace('_fastqc_data.txt', '_fastqc.zip')
    se_stats = parse_fastqc_data('${data_txt}', original_fastq_zip)
    sampling_size = se_stats['total_sequences']

    final_data = {
        'r1': se_stats,
        'r2': None,
        'sampling_size': sampling_size
    }

    with open('${meta}_fastqc.json', 'w') as f:
        json.dump(final_data, f, indent=4)

    print(f"Successfully parsed FastQC data for ${meta}")
    """
}

process PARSE_QUAST {
    tag "$meta"
    label 'parsing_process'

    input:
    tuple val(meta), path(quast_tsv), path(assembly)

    output:
    tuple val(meta), path("${meta}_quast.json"), emit: quast_json

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import json
    from Bio import SeqIO

    # --- Parse QUAST TSV ---
    df = pd.read_csv('${quast_tsv}', sep='\\t', index_col=0)
    quast_dict = df.T.to_dict(orient='records')[0]

    # Convert to int
    for key, value in quast_dict.items():
        try:
            if isinstance(value, str) and ',' in value:
                value = value.split(',')[0]
            quast_dict[key] = int(value)
        except (ValueError, TypeError):
            pass

    # --- Count contigs by size from FASTA ---
    contigs_lt_300 = 0
    contigs_300_999 = 0
    contigs_ge_1kb = 0
    
    for record in SeqIO.parse('${assembly}', 'fasta'):
        length = len(record.seq)
        if length < 300:
            contigs_lt_300 += 1
        elif length < 1000:
            contigs_300_999 += 1
        else:
            contigs_ge_1kb += 1
    
    quast_dict['contigs_lt_300'] = contigs_lt_300
    quast_dict['contigs_300_999'] = contigs_300_999
    quast_dict['contigs_ge_1kb'] = contigs_ge_1kb

    # Custom 'rgl' calculation
    total_len = quast_dict.get('Total length', 0)
    len_gt_1k = quast_dict.get('Total length (>= 1000 bp)', 0)
    
    if total_len > 0:
        rgl = round((len_gt_1k / total_len) * 100, 2)
    else:
        rgl = 0.0
    
    quast_dict['rgl'] = rgl

    with open('${meta}_quast.json', 'w') as f:
        json.dump(quast_dict, f, indent=4)
    """
}

process PARSE_MLST {
    tag "$meta"
    label 'parsing_process'

    input:
    tuple val(meta), path(mlst_txt)
    path(lookup_json)

    output:
    tuple val(meta), path("${meta}_mlst.json"),          emit: mlst_json

    script:
    """
    python -c "
import json, os
results = {'scheme': 'N/A', 'type': 'N/A', 'detail': [], 'species_lookup': 'N/A', 'size_lookup': 'N/A'}
if os.path.exists('${mlst_txt}') and os.path.getsize('${mlst_txt}') > 0:
    with open('${mlst_txt}', 'r') as f:
        parts = f.readline().strip().split()
        if len(parts) >= 3:
            results['scheme'], results['type'], results['detail'] = parts[1], parts[2], parts[3:]
if os.path.exists('${lookup_json}') and results['scheme'] != 'N/A':
    with open('${lookup_json}', 'r') as f:
        lookup_data = json.load(f)
        scheme_data = lookup_data.get(results['scheme'], {})
        results['species_lookup'] = scheme_data.get('Scientific Name', 'N/A')
        results['size_lookup'] = scheme_data.get('Approx Genome Size (Mb)', 'N/A')
with open('${meta}_mlst.json', 'w') as f: json.dump(results, f, indent=4)
"
    """
}

process PARSE_SEQSERO {
    tag "$meta"
    label 'parsing_process'

    input:
    tuple val(meta), path(seqsero_result)  // This will be a list!

    output:
    tuple val(meta), path("${meta}_seqsero.json"),       emit: seqsero_json

    script:
    """
    python3 << 'EOF'
import json
import os
import sys

results = {
    'serotype': 'N/A',
    'antigenic_profile': 'N/A', 
    'subspecies': 'N/A',
    'note': ''
}

# When multiple files are staged, they appear as separate files in the work directory
# Find the TSV file
tsv_file = None
for fname in os.listdir('.'):
    if fname.endswith('.tsv') and 'SeqSero_result' in fname:
        tsv_file = fname
        print(f"Found TSV file: {tsv_file}")
        break

if not tsv_file:
    print("WARNING: No SeqSero TSV file found!")
    print(f"Files present: {os.listdir('.')}")
else:
    try:
        with open(tsv_file, 'r') as f:
            lines = [line.strip() for line in f.readlines() if line.strip()]
        
        print(f"Read {len(lines)} lines from {tsv_file}")
        
        if len(lines) > 1:
            header = lines[0].split('\\t')
            data = lines[1].split('\\t')
            
            print(f"Header columns: {len(header)}")
            print(f"Data columns: {len(data)}")
            
            # SeqSero2 TSV format:
            # 0: Sample, 1: Output dir, 2: Input files, 3: O antigen, 4: H1, 5: H2,
            # 6: Predicted identification, 7: Antigenic profile, 8: Serotype, 9: Note
            if len(data) >= 9:
                results['subspecies'] = data[6].strip() if data[6].strip() else 'N/A'
                results['antigenic_profile'] = data[7].strip() if data[7].strip() else 'N/A'
                results['serotype'] = data[8].strip() if data[8].strip() else 'N/A'
                results['note'] = data[9].strip() if len(data) > 9 and data[9].strip() else ''
                
                print(f"✓ Parsed: {results['serotype']} ({results['antigenic_profile']})")
            else:
                print(f"WARNING: Only {len(data)} columns, expected at least 9")
                
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()

# Write output
with open('${meta}_seqsero.json', 'w') as f:
    json.dump(results, f, indent=4)
    
print(f"\\nFinal: {json.dumps(results, indent=2)}")
EOF
    """
}

process PARSE_ABRICATE {
    tag "$meta"
    label 'parsing_process'

    input:
    tuple val(meta), path(resfinder_tab), path(ncbi_tab), path(plasmidfinder_tab)

    output:
    tuple val(meta), path("${meta}_abricate.json"),      emit: abricate_json

    script:
    """
    #!/bin/bash
    set -e
    NCBI_COUNT=\$( (grep -vc '^#' ${ncbi_tab} || true) | awk '{print (\$1 > 0) ? \$1 - 1 : 0}' )
    RESFINDER_COUNT=\$( (grep -vc '^#' ${resfinder_tab} || true) | awk '{print (\$1 > 0) ? \$1 - 1 : 0}' )
    ABRICATE_VERSION=\$(abricate --version | awk '{print \$2}')
    NCBI_DB_INFO=\$(abricate --list | grep '^ncbi'); RESFINDER_DB_INFO=\$(abricate --list | grep '^resfinder')
    NCBI_DB_SEQS=\$(echo "\$NCBI_DB_INFO" | awk '{print \$2}'); NCBI_DB_DATE=\$(echo "\$NCBI_DB_INFO" | awk '{print \$4}')
    RESFINDER_DB_SEQS=\$(echo "\$RESFINDER_DB_INFO" | awk '{print \$2}'); RESFINDER_DB_DATE=\$(echo "\$RESFINDER_DB_INFO" | awk '{print \$4}')
    python -c "
import json
data = {
    'abricate_version': '\$ABRICATE_VERSION', 
    'hit_count_ncbi': int(\$NCBI_COUNT), 
    'hit_count_resfinder': int(\$RESFINDER_COUNT),
    'mincov': ${params.abricate_depth},
    'minid': ${params.abricate_coverage},
    'db_info': {
        'ncbi': {'sequences': '\$NCBI_DB_SEQS', 'date': '\$NCBI_DB_DATE'},
        'resfinder': {'sequences': '\$RESFINDER_DB_SEQS', 'date': '\$RESFINDER_DB_DATE'}
    }
}
with open('${meta}_abricate.json', 'w') as f: json.dump(data, f, indent=4)
"
    """
}

process PARSE_AMRFINDER {
    tag "$meta"
    label 'parsing_process'

    input:
    tuple val(meta), path(amrfinder_tab)

    output:
    tuple val(meta), path(amrfinder_tab), path("${meta}_amrfinder.json"),        emit: amrfinder_parsed
    
    script:
    """
    #!/bin/bash
    set -e
    AMRFINDER_VERSION=\$(amrfinder --version | head -n 1)
    python -c "
import json
data = {'amrfinder_version': '\$AMRFINDER_VERSION', 'results_file': '${amrfinder_tab}'}
with open('${meta}_amrfinder.json', 'w') as f: json.dump(data, f, indent=4)
"
    """
}

process SOFTWARE_VERSIONS {
    
    tag "Software Versions"
    
    publishDir (
        path: "./pipeline_info", 
        mode: 'copy'
    )

    output:
    path "software_versions.txt", emit: versions_report

    script:
    """
    #!/bin/bash
    set +e
    
    exec > software_versions.txt 2>&1

    export ORIGINAL_PATH="\${PATH}"
    
    eval "\\\$(conda shell.bash hook)" 2>/dev/null || {
        if [ -f "/opt/conda/etc/profile.d/conda.sh" ]; then
            source /opt/conda/etc/profile.d/conda.sh
        elif [ -f "\${CONDA_PREFIX}/etc/profile.d/conda.sh" ]; then
            source \${CONDA_PREFIX}/etc/profile.d/conda.sh
        fi
    }
    
    echo "--- kraken_env ---"
    if conda activate kraken_env 2>/dev/null; then
        # Restore PATH utilities
        export PATH="\${CONDA_PREFIX}/bin:\${ORIGINAL_PATH}"
        
        echo "# packages in environment at \${CONDA_PREFIX}:"
        echo "#"
        echo "# Name                     Version"
        conda list 2>/dev/null | /usr/bin/grep -E "^(_openmp_mutex|abyss|bedtools|biopython|blast|bracken|bwa|fastahack|freebayes|gperftools|htslib|humanize|kraken2|krona)" | /usr/bin/awk '{print \$1, \$2}' | /usr/bin/sort
        echo ""
        
        if command -v kraken2 &> /dev/null; then
            echo -n "  kraken2: "
            kraken2 --version 2>&1 | /usr/bin/head -n 1 | /usr/bin/sed 's/^Kraken version //'
        fi
        if command -v bracken &> /dev/null; then
            echo -n "  bracken: "
            bracken -v 2>&1 | /usr/bin/head -n 1
        fi
        conda deactivate 2>/dev/null
    else
        echo "Warning: Could not activate kraken_env"
    fi
    echo ""
    echo ""
    
    echo "--- assembly_env ---"
    if conda activate assembly_env 2>/dev/null; then
        export PATH="\${CONDA_PREFIX}/bin:\${ORIGINAL_PATH}"
        
        echo "# packages in environment at \${CONDA_PREFIX}:"
        echo "#"
        echo "# Name                     Version"
        conda list 2>/dev/null | /usr/bin/grep -E "^(abyss|bedtools|biopython|bwa|fastahack|freebayes|htslib|humanize|spades|seqsero2)" | /usr/bin/awk '{print \$1, \$2}' | /usr/bin/sort
        echo ""
        
        if command -v spades.py &> /dev/null; then
            echo -n "  spades.py: "
            spades.py --version 2>&1 | /usr/bin/grep -o "SPAdes.*" | /usr/bin/head -n 1
        fi
        if command -v SeqSero2_package.py &> /dev/null; then
            echo -n "  SeqSero2_package.py: "
            SeqSero2_package.py --version 2>&1 | /usr/bin/head -n 1
        fi
        conda deactivate 2>/dev/null
    else
        echo "Warning: Could not activate assembly_env"
    fi
    echo ""
    echo ""
    
    echo "--- legacy_perl_env ---"
    if conda activate legacy_perl_env 2>/dev/null; then
        export PATH="\${CONDA_PREFIX}/bin:\${ORIGINAL_PATH}"
        
        echo "# packages in environment at \${CONDA_PREFIX}:"
        echo "#"
        echo "# Name                     Version"
        conda list 2>/dev/null | /usr/bin/grep -E "^(perl|mlst|abricate|python|pandas|openpyxl|numpy) " | /usr/bin/awk '{print \$1, \$2}' | /usr/bin/sort
        echo ""
        
        if command -v mlst &> /dev/null; then
            echo -n "  mlst: "
            mlst --version 2>&1 | /usr/bin/head -n 1
        fi
        if command -v abricate &> /dev/null; then
            echo -n "  abricate: "
            abricate --version 2>&1 | /usr/bin/head -n 1
        fi
        conda deactivate 2>/dev/null
    else
        echo "Warning: Could not activate legacy_perl_env"
    fi
    echo ""
    echo ""
    
    echo "--- analysis_env ---"
    if conda activate analysis_env 2>/dev/null; then
        export PATH="\${CONDA_PREFIX}/bin:\${ORIGINAL_PATH}"
        
        echo "# packages in environment at \${CONDA_PREFIX}:"
        echo "# Name                                     Version"
        conda list 2>/dev/null | /usr/bin/grep -E "^(bedtools|biopython|blast|bwa|cairo|cairocffi|cairosvg|amrfinder|numpy|pandas|openpyxl)" | /usr/bin/awk '{print \$1, \$2}' | /usr/bin/sort
        echo ""
        
        if command -v amrfinder &> /dev/null; then
            echo -n "  amrfinder: "
            amrfinder --version 2>&1 | /usr/bin/grep -o "[0-9]\\+\\.[0-9]\\+\\.[0-9]\\+" | /usr/bin/head -n 1
        fi
        if command -v python &> /dev/null; then
            echo -n "  python: "
            python --version 2>&1
        fi
        echo ""
        
        python -c "import pandas; print('  pandas:', pandas.__version__)" 2>/dev/null || echo "  pandas: NOT FOUND"
        python -c "import numpy; print('  numpy:', numpy.__version__)" 2>/dev/null || echo "  numpy: NOT FOUND"
        python -c "import openpyxl; print('  openpyxl:', openpyxl.__version__)" 2>/dev/null || echo "  openpyxl: NOT FOUND"
        conda deactivate 2>/dev/null
    else
        echo "Warning: Could not activate analysis_env"
    fi
    echo ""
    """
}

process GENERATE_PDF_REPORT {
    tag "$meta"
    label 'report_process'
    publishDir "${params.outdir}/${meta}", mode: 'copy', pattern: "*.pdf"
    publishDir "${params.outdir}/${meta}", mode: 'copy', pattern: "*_stats.xlsx"
    publishDir "${params.outdir}/${meta}", mode: 'copy', pattern: "*.{fastq,fastq.gz,fq,fq.gz}"
    publishDir "${params.outdir}/${meta}", mode: 'copy', pattern: "*.{fasta,fa,fna}"

    input:
    tuple val(meta), 
          path(reads, stageAs: 'reads_input/*'), 
          path(assembly, stageAs: 'assembly_input/*'),
          path(fastqc_json, stageAs: 'fastqc_input.json'), 
          path(quast_json, stageAs: 'quast_input.json'),
          path(mlst_json, stageAs: 'mlst_input.json'), 
          path(seqsero_json, stageAs: 'seqsero_input.json'),
          path(abricate_json, stageAs: 'abricate_json.json'), 
          path(amrfinder_tab, stageAs: 'amrfinder_input.tab'), 
          path(amrfinder_json, stageAs: 'amrfinder_input.json'),
          path(resfinder_tab, stageAs: 'resfinder_input.tab'), 
          path(ncbi_tab, stageAs: 'ncbi_input.tab'),
          path(bracken_pie, stageAs: 'bracken_pie.png'),
          path(software_versions, stageAs: 'versions_input.txt')
    each path(logo_file)

    output:
    tuple val(meta), path("*.pdf"),                        emit: pdf
    tuple val(meta), path("*_stats.xlsx"),                 emit: xlsx
    path "*.{fastq,fastq.gz,fq,fq.gz}", optional: true,    emit: reads_out
    path "*.{fasta,fa,fna}", optional: true,               emit: assembly_out

    script:
    // Check if we have real reads (not NO_FILE placeholder)
    def reads_files = reads instanceof List ? reads : [reads]
    def real_reads = reads_files.findAll { it.getSimpleName() != 'NO_FILE' }
    def has_real_reads = real_reads.size() > 0

    // Check if we have real fastqc json (not NO_FILE placeholder)
    def has_fastqc = fastqc_json.getSimpleName() != 'NO_FILE' && fastqc_json.size() > 0

    // Check if we have real seqsero json (not NO_FILE placeholder)
    def has_seqsero = seqsero_json.getSimpleName() != 'NO_FILE' && seqsero_json.size() > 0

    // Check if we have real bracken pie chart (not NO_FILE placeholder)
    def has_bracken_pie = bracken_pie.getSimpleName() != 'NO_FILE' && bracken_pie.size() > 0

    // Build command line arguments - only add optional inputs if they exist
    def args = []
    args << "generate_final_report.py"
    args << "--meta_id ${meta}"
    args << "--logo ${logo_file}"
    
    // Only add --reads if we have real FASTQ files
    if (has_real_reads) {
        args << "--reads reads_input/*"
    }
    
    args << "--assembly assembly_input/*"
    if (has_fastqc) args << "--fastqc_json ${fastqc_json}"
    args << "--quast_json ${quast_json}"
    args << "--mlst_json ${mlst_json}"
    if (has_seqsero) args << "--seqsero_json ${seqsero_json}"
    args << "--abricate_json ${abricate_json}"
    args << "--amrfinder_tab ${amrfinder_tab}"
    args << "--amrfinder_json ${amrfinder_json}"
    args << "--versions ${software_versions}"
    args << "--ab_ncbi_file ${ncbi_tab}"
    args << "--ab_resfinder_file ${resfinder_tab}"
    if (has_bracken_pie) args << "--bracken_pie ${bracken_pie}"

    // Join arguments with line continuation for readability
    def full_cmd = args.join(" \\\n    ")

    """
    # Copy assembly FASTA to work dir root for publishing
    cp -L assembly_input/* . 2>/dev/null || true
    
    # Copy input reads to work dir root for publishing (if they exist and are real)
    if [ -d reads_input ] && [ "\$(find reads_input -type f ! -name 'NO_FILE' | wc -l)" -gt 0 ]; then
        find reads_input -type f ! -name 'NO_FILE' -exec cp -L {} . \\;
    fi

    # Run the report generation
    ${full_cmd}
    """
}

process COMBINE_XLSX {
    tag "combined_${seq_type}"
    label 'parsing_process'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path xlsx_files
    val seq_type  // "PE", "SE", or "FA"

    output:
    path "combined_stats_${seq_type}.xlsx", emit: combined_xlsx

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import glob

    seq_type = "${seq_type}"
    
    excel_files = sorted(glob.glob("*.xlsx"))
    
    if not excel_files:
        print("ERROR: No Excel files found")
        exit(1)
    
    print(f"Combining {len(excel_files)} {seq_type} Excel files...")
    
    all_rows = []
    index_name = None  # Store the index column name
    for filename in excel_files:
        try:
            df = pd.read_excel(filename, index_col=0)
            if not df.empty:
                # Capture the index name from the first file
                if index_name is None:
                    index_name = df.index.name or 'sample'
                all_rows.append(df.iloc[0])
                print(f"  ✓ {filename}")
        except Exception as e:
            print(f"  ✗ Failed: {filename}: {e}")
    
    if not all_rows:
        print("ERROR: No valid data extracted")
        exit(1)
    
    final_df = pd.DataFrame(all_rows)
    # Restore the index name
    final_df.index.name = index_name
    
    # Generate output filename based on sequencing type
    if seq_type == "FA":
        output_file = "combined_stats_FA.xlsx"
    elif seq_type == "SE":
        output_file = "combined_stats_SE.xlsx"
    elif seq_type == "PE":
        output_file = "combined_stats_PE.xlsx"
    else:
        output_file = f"combined_stats_{seq_type}.xlsx"
    
    final_df.to_excel(output_file)
    
    print(f"\\n✓ Combined {len(final_df)} {seq_type} samples into {output_file}")
    """
}

workflow amr_pe {
    take:
        paired_fastqs

    main:
        ch_all_individual_reads = paired_fastqs
            .flatMap { meta, reads -> 
                [
                    tuple([id: meta, read_type: 'R1'], reads[0]), 
                    tuple([id: meta, read_type: 'R2'], reads[1])
                ] 
            }

        FASTQC( ch_all_individual_reads )
        
        ch_fastqc_r1 = FASTQC.out.data_txt
            .filter { meta, file -> meta.read_type == 'R1' }
            .map { meta, file -> tuple(meta.id, file) }
        
        
        ch_fastqc_r2 = FASTQC.out.data_txt
            .filter { meta, file -> meta.read_type == 'R2' }
            .map { meta, file -> tuple(meta.id, file) }
        
        ch_pe_to_parse = ch_fastqc_r1.join(ch_fastqc_r2)

        PARSE_FASTQC( ch_pe_to_parse )

        KRAKEN_PE( paired_fastqs )
        BRACKEN( KRAKEN_PE.out.report )
        BRACKEN_PIE( BRACKEN.out.bracken_excel )
        KRONA( KRAKEN_PE.out.report )

        SPADES_PE (paired_fastqs )
        spades_ch = SPADES_PE.out.assembly
        QUAST( spades_ch )
  
        SEQSERO_PE( paired_fastqs )
        MLST( spades_ch )
        ABRICATE( spades_ch )

        abricate_ch = ABRICATE.out.resfinder_tab
            .join( ABRICATE.out.ncbi_tab )
            .join( ABRICATE.out.plasmidfinder_tab )

        amrfinder_ch = spades_ch.join( BRACKEN.out.bracken_excel )
        AMRFINDER( amrfinder_ch )

        parse_quast_input = QUAST.out.report_tsv.join(SPADES_PE.out.assembly)
        PARSE_QUAST(parse_quast_input)
        PARSE_SEQSERO( SEQSERO_PE.out.result_file )
        PARSE_MLST( MLST.out.mlst_txt, params.mlst_lookup_json )
        PARSE_ABRICATE( abricate_ch )
        PARSE_AMRFINDER( AMRFINDER.out.amrfinder_tab )

        SOFTWARE_VERSIONS()

        ch_report_inputs = paired_fastqs
            .join( SPADES_PE.out.assembly )
            .join( PARSE_FASTQC.out.fastqc_json )
            .join( PARSE_QUAST.out.quast_json )
            .join( PARSE_MLST.out.mlst_json )
            .join( PARSE_ABRICATE.out.abricate_json )
            .join( PARSE_AMRFINDER.out.amrfinder_parsed )
            .join( ABRICATE.out.resfinder_tab )
            .join( ABRICATE.out.ncbi_tab )
            .join( PARSE_SEQSERO.out.seqsero_json, remainder: true)
            .join(BRACKEN_PIE.out.pie_chart)
            .map { tuple ->
                def meta = tuple[0]
                def reads = tuple[1]
                def assembly = tuple[2]
                def fastqc_json = tuple[3]
                def quast_json = tuple[4]
                def mlst_json = tuple[5]
                def abricate_json = tuple[6]
                def amrfinder_tab = tuple[7]
                def amrfinder_json = tuple[8]
                def resfinder_tab = tuple[9]
                def ncbi_tab = tuple[10]
                def seqsero_json = tuple[11]
                def braken_pie = tuple[12]

                def seqsero_file = (seqsero_json != null && seqsero_json != []) ? seqsero_json : file('NO_FILE')
                
                return [
                    meta, reads, assembly,
                    fastqc_json, quast_json, mlst_json, 
                    seqsero_file,
                    abricate_json, amrfinder_tab, amrfinder_json,
                    resfinder_tab, ncbi_tab, braken_pie
                ]
            }
            .combine( SOFTWARE_VERSIONS.out.versions_report )

        GENERATE_PDF_REPORT(ch_report_inputs, file(params.logo))

        all_xlsx = GENERATE_PDF_REPORT.out.xlsx.map { meta, xlsx -> xlsx }.collect()
        COMBINE_XLSX(all_xlsx, "PE")
}

workflow amr_se {
    take:
        single_fastqs

    main:
        ch_single_reads = single_fastqs
            .map { meta, file -> tuple([id: meta, read_type: 'SE'], file) }

        FASTQC(ch_single_reads)
        
        ch_fastqc_to_parse = FASTQC.out.data_txt
            .map { meta, file -> tuple(meta.id, file) }

        PARSE_FASTQC_SE(ch_fastqc_to_parse)

        KRAKEN_SE(single_fastqs)
        BRACKEN(KRAKEN_SE.out.report)
        BRACKEN_PIE(BRACKEN.out.bracken_excel)
        KRONA(KRAKEN_SE.out.report)

        SPADES_SE(single_fastqs)
        spades_ch = SPADES_SE.out.assembly
        QUAST(spades_ch)
  
        SEQSERO_SE(single_fastqs)
        MLST(spades_ch)
        ABRICATE(spades_ch)

        abricate_ch = ABRICATE.out.resfinder_tab
            .join(ABRICATE.out.ncbi_tab)
            .join(ABRICATE.out.plasmidfinder_tab)

        amrfinder_ch = spades_ch.join(BRACKEN.out.bracken_excel)
        AMRFINDER(amrfinder_ch)

        parse_quast_input = QUAST.out.report_tsv.join(SPADES_SE.out.assembly)
        PARSE_QUAST(parse_quast_input)
        PARSE_SEQSERO(SEQSERO_SE.out.result_file)
        PARSE_MLST(MLST.out.mlst_txt, params.mlst_lookup_json)
        PARSE_ABRICATE(abricate_ch)
        PARSE_AMRFINDER(AMRFINDER.out.amrfinder_tab)

        SOFTWARE_VERSIONS()

        ch_report_inputs = single_fastqs
            .join(QUAST.out.report_tsv)
            .join(PARSE_FASTQC_SE.out.fastqc_json)
            .map { meta, reads, quast, fastqc_json -> 
                tuple(meta, reads, file('NO_FILE'), quast, fastqc_json)
            }
            .join(PARSE_QUAST.out.quast_json)
            .join(PARSE_MLST.out.mlst_json)
            .join(PARSE_ABRICATE.out.abricate_json)
            .join(PARSE_AMRFINDER.out.amrfinder_parsed)
            .join(ABRICATE.out.resfinder_tab)
            .join(ABRICATE.out.ncbi_tab)
            .join(PARSE_SEQSERO.out.seqsero_json, remainder: true)
            .map { tuple ->
                def meta = tuple[0]
                def reads = tuple[1] 
                def assembly = tuple[2]
                def quast_report = tuple[3]
                def fastqc_json = tuple[4]
                def quast_json = tuple[5]
                def mlst_json = tuple[6]
                def abricate_json = tuple[7]
                def amrfinder_tab = tuple[8]
                def amrfinder_json = tuple[9]
                def resfinder_tab = tuple[10]
                def ncbi_tab = tuple[11]
                def seqsero_json = tuple[12]

                def seqsero_file = (seqsero_json != null && seqsero_json != []) ? seqsero_json : file('NO_FILE')
                
                return [
                    meta, 
                    reads,
                    assembly,
                    fastqc_json,
                    quast_json,      
                    mlst_json,        
                    seqsero_file,   
                    abricate_json,   
                    amrfinder_tab,
                    amrfinder_json,  
                    resfinder_tab,    
                    ncbi_tab,         
                    file('NO_FILE')
                ]
            }
            .combine(SOFTWARE_VERSIONS.out.versions_report) 

        GENERATE_PDF_REPORT(ch_report_inputs, file(params.logo))

        all_xlsx = GENERATE_PDF_REPORT.out.xlsx.map { meta, xlsx -> xlsx }.collect()
        COMBINE_XLSX(all_xlsx, "SE")
}

workflow amr_fa {
    take:
        fastas

    main:
        // FASTA files skip: FASTQC, KRAKEN, BRACKEN, SPADES
        
        QUAST(fastas)
        SEQSERO_FA(fastas)
        MLST(fastas)
        ABRICATE(fastas)

        abricate_ch = ABRICATE.out.resfinder_tab
            .join(ABRICATE.out.ncbi_tab)
            .join(ABRICATE.out.plasmidfinder_tab)

        PARSE_MLST(MLST.out.mlst_txt, params.mlst_lookup_json)

        amrfinder_fa = fastas.join(PARSE_MLST.out.mlst_json)
        AMRFINDER_FA(amrfinder_fa)
        
        parse_quast_input = QUAST.out.report_tsv.join(fastas)
        PARSE_QUAST(parse_quast_input)
        PARSE_SEQSERO(SEQSERO_FA.out.result_file)
        PARSE_ABRICATE(abricate_ch)
        PARSE_AMRFINDER(AMRFINDER_FA.out.amrfinder_tab)

        SOFTWARE_VERSIONS()

        ch_report_inputs = fastas
            .join(QUAST.out.report_tsv)
            .map { meta, fasta, quast -> 
                tuple(meta, file('NO_FILE'), fasta, quast)
            }
            .join(PARSE_QUAST.out.quast_json)
            .join(PARSE_MLST.out.mlst_json)
            .join(PARSE_ABRICATE.out.abricate_json)
            .join(PARSE_AMRFINDER.out.amrfinder_parsed)
            .join(ABRICATE.out.resfinder_tab)
            .join(ABRICATE.out.ncbi_tab)
            .join(PARSE_SEQSERO.out.seqsero_json, remainder: true)
            .map { tuple ->
                def meta = tuple[0]
                def reads = tuple[1]
                def assembly = tuple[2]
                def quast_report = tuple[3]    
                def quast_json = tuple[4]
                def mlst_json = tuple[5]
                def abricate_json = tuple[6]
                def amrfinder_tab = tuple[7]
                def amrfinder_json = tuple[8]
                def resfinder_tab = tuple[9]
                def ncbi_tab = tuple[10]
                def seqsero_json = tuple[11]

                def seqsero_file = (seqsero_json != null && seqsero_json != []) ? seqsero_json : file('NO_FILE')
                
                return [
                    meta, 
                    reads,
                    assembly,
                    file('NO_FILE'),
                    quast_json,
                    mlst_json,
                    seqsero_file,
                    abricate_json,
                    amrfinder_tab,
                    amrfinder_json,
                    resfinder_tab,
                    ncbi_tab,
                    file('NO_FILE')
                ]
            }
            .combine(SOFTWARE_VERSIONS.out.versions_report)

        GENERATE_PDF_REPORT(ch_report_inputs, file(params.logo))

        all_xlsx = GENERATE_PDF_REPORT.out.xlsx.map { meta, xlsx -> xlsx }.collect()
        COMBINE_XLSX(all_xlsx, "FA")
}

// ---MAIN WORKFLOW---
workflow {
    Channel
        .fromPath(params.input_files, checkIfExists: false)
        .set { ch_all_reads }

    ch_all_reads.branch {
        single: it.name =~ /_se\.fastq\.gz$/
        paired: it.name =~ /_R[12](_001)?\.fastq\.gz$/ &&
                !(it.name =~ /_se\.fastq\.gz$/)
        fasta: it.name =~ /\.(fasta|fa)$/
        unmatched: true
    }
    .set { sorted_reads }

    sorted_reads.unmatched.view { "UNMATCHED FILE: ${it.name}" }

    paired_fastqs = sorted_reads.paired
        .map { file ->
            def sample_id = file.name
                .replaceAll(/_S\d+(_L\d{3})?/, '')
                .replaceAll(/_R[12](_001)?\.fastq\.gz$/, '')
            tuple(sample_id, file)
        }
        .groupTuple(size: 2)

    single_fastqs = sorted_reads.single
        .map { file ->
            def sample_id = file.name
                .replaceAll(/_S\d+(_L\d{3})?/, '')
                .replaceAll(/(_R1)?_se\.fastq\.gz$/, '') // This is really a placeholder since we will implent ONT data in the future
            tuple(sample_id, file)
        }

    fastas = sorted_reads.fasta
        .map { file ->
            def sample_id = file.name
                .replaceAll(/_S\d+(_L\d{3})?/, '')
                .replaceAll(/\.(fasta|fa)$/, '')
            tuple(sample_id, file)
        }

    // Debug output
    paired_fastqs.view { "✓ PAIRED: $it" }
    single_fastqs.view { "✓ SINGLE: $it" }
    fastas.view { "✓ FASTA: $it" }

    // Call sub-workflows
    amr_pe(paired_fastqs)
    amr_se(single_fastqs)
    amr_fa(fastas)
}
