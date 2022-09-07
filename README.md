# amr_wrapper.py usage:

## Paths:

Must include path to `amr_wrapper.py` in `$PATH`.<br>
Must also include this path in `$PYTHONPATH`.

Within the directory containing `amr_wrapper.py` these other 9 files are included:

- `abricate_wrapper.py`
- `amr_wrapper.py`
- `amrfinder_wrapper.py`
- `latex_reporter.py`
- `lookup_genome_size.json`
- `mlst_wrapper.py`
- `rmlst.py`
- `seqsero2_wrapper.py`
- `sequence_score.py`

In addition to the above files the following Python files must be in `$PYTHONPATH`.

- `vsnp_fastq_quality.py`
- `spades_assembly.py`
- `spades_stats_parse.py`
- `kraken2_run.py`

For example:


## Run:

See `amr_wrapper.py` for run options.

Basic usage:

```zsh
amr_wrapper.py -r1 *_R1*fastq.gz -r2 *_R2*fastq.gz
```

# NAHLN_AMR
# NAHLN_AMR
