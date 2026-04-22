#!/usr/bin/env python
__version__ = "0.0.2"
import os
import sys
import math
import numpy as np

class Quality_Scaling:

    def fastq_scaling(self, genome_size, read1_reads_gt_q30=None, read2_reads_gt_q30=None, sampling_size=None, read1_total_read_count=None, read2_total_read_count=None, read1_read_average=None, read2_read_average=None, read1_length_mean=None, read2_length_mean=None):
        '''
        Calculates a quality score for FASTQ data and returns the score and coverage depth.
        Handles both single-end and paired-end cases.
        '''
        if not read2_total_read_count or read2_total_read_count == 0:
            # --- SINGLE-END LOGIC ---
            read_length = float(read1_length_mean) if read1_length_mean else 150
            genome_coverage_depth = (int(read1_total_read_count) * read_length) / (int(genome_size)) if genome_size > 0 else 0
            
            q30_score = ((read1_reads_gt_q30 / sampling_size) * 100) if sampling_size and sampling_size > 0 else 0
            mean_qual_score = float(read1_read_average) if read1_read_average else 0
            
            # Simple heuristic to create a normalized score
            normalized_score = ((q30_score - 70) / 30) + ((mean_qual_score - 30) / 10)
            fastq_scaling_variable = 1180 + (normalized_score * 400)
            
            # Clamp the value to the visual range of the slider
            if fastq_scaling_variable > 1980: fastq_scaling_variable = 1980
            elif fastq_scaling_variable < 1180: fastq_scaling_variable = 1180
                
            return fastq_scaling_variable, genome_coverage_depth
        
        else:
            # --- PAIRED-END LOGIC ---
            r1_length = float(read1_length_mean) if read1_length_mean else 240
            r2_length = float(read2_length_mean) if read2_length_mean else 240
            genome_coverage_depth = (int(read1_total_read_count) * r1_length + int(read2_total_read_count) * r2_length) / (int(genome_size)) if genome_size > 0 else 0

            seqcor = [[1.0000000, 0.4412253, 0.9569653, 0.4266978, 0.3393063],
                      [0.4412253, 1.0000000, 0.4446494, 0.9527558, 0.2464773],
                      [0.9569653, 0.4446494, 1.0000000, 0.4330931, 0.3315736],
                      [0.4266978, 0.9527558, 0.4330931, 1.0000000, 0.2367824],
                      [0.3393063, 0.2464773, 0.3315736, 0.2367824, 1.0000000]]

            q301 = (read1_reads_gt_q30 / sampling_size) * 100 if sampling_size and sampling_size > 0 else 0
            q302 = (read2_reads_gt_q30 / sampling_size) * 100 if sampling_size and sampling_size > 0 else 0
            
            x = [(float(read1_read_average) - 30), (float(read2_read_average) - 30), (float(q301) - 80), (float(q302) - 80), genome_coverage_depth]
            mx = np.matrix(x)
            mseqcor = np.matrix(seqcor)
            mxt = mx.transpose()
            
            # Corrected matrix multiplication
            seqv = math.sqrt(mx * mseqcor * mxt)
            fastq_scaling_variable = (seqv * 800 / 90) + 1180

            if fastq_scaling_variable > 1980:
                fastq_scaling_variable = 1980
            elif fastq_scaling_variable < 1180:
                fastq_scaling_variable = 1180
            
            return fastq_scaling_variable, genome_coverage_depth

    def assembly_scaling(self, longest_contig=None, greater_one_kb_count=None, contig_count=None, n50=None, l50=None, rgl=None, stat_total_contig_lengths=None, stat_contig_count=None, genome_coverage_depth=None):
        '''
        Calculates a quality score for the assembly.
        '''
        # Corrected indentation and ensured all inputs have default values
        longest_contig = longest_contig or 0
        greater_one_kb_count = greater_one_kb_count or 0
        contig_count = contig_count or 1 # Avoid division by zero
        n50 = n50 or 0
        l50 = l50 or 0
        rgl = rgl or 0
        stat_total_contig_lengths = stat_total_contig_lengths or 1 # Avoid division by zero
        genome_coverage_depth = genome_coverage_depth or 1 # Avoid division by zero

        assemblecor = [[1.00000000, 0.4640352, 0.39099636, -0.703376083, -0.074391092, -0.27033841, -0.26923388, -0.36627871],
                       [0.46403525, 1.0000000, 0.75927035, -0.220017666, -0.101805557, -0.56550423, 0.27504920, -0.70388508],
                       [0.39099636, 0.7592704, 1.00000000, -0.170697061, -0.054252944, -0.72295517, 0.29979708, -0.81000716],
                       [-0.70337608, -0.2200177, -0.17069706, 1.000000000, 0.003833981, 0.09220126, 0.46237217, 0.16084638],
                       [-0.07439109, -0.1018056, -0.05425294, 0.003833981, 1.000000000, 0.05500119, -0.04451447, 0.01902384],
                       [-0.27033841, -0.5655042, -0.72295517, 0.092201260, 0.055001186, 1.00000000, -0.28768966, 0.55006987],
                       [-0.26923388, 0.2750492, 0.29979708, 0.462372175, -0.044514471, -0.28768966, 1.00000000, -0.27338410],
                       [-0.36627871, -0.7038851, -0.81000716, 0.160846384, 0.019023843, 0.55006987, -0.27338410, 1.00000000]]

        massemblecor = np.matrix(assemblecor)
        
        # Corrected and made safe from division by zero errors
        perlongexpect = (float(longest_contig) / (float(genome_coverage_depth) * 1000000)) * 100 if genome_coverage_depth > 0 else 0
        perlongscaff = (float(greater_one_kb_count) / float(contig_count)) * 100 if contig_count > 0 else 0
        pern50expect = (int(n50) / float(stat_total_contig_lengths)) * 100 if stat_total_contig_lengths > 0 else 0
        perdiffexpect = abs(100 - (float(stat_total_contig_lengths) / (float(stat_total_contig_lengths) * 10000))) if stat_total_contig_lengths > 0 else 100

        quality_scale_list = [contig_count, greater_one_kb_count, l50, rgl, perdiffexpect, perlongexpect, perlongscaff, pern50expect]
        ma = np.matrix(quality_scale_list)
        mat = ma.transpose()
        
        # Use try-except for the matrix math which can fail
        try:
            assemblev = 1000 / (math.sqrt(ma * massemblecor * mat))
        except (ValueError, ZeroDivisionError):
            assemblev = 0

        assembly_scaling_variable = (assemblev * 800 / 7) + 1180

        if assembly_scaling_variable > 1980:
            assembly_scaling_variable = 1980
        elif assembly_scaling_variable < 1180:
            assembly_scaling_variable = 1180
            
        return assembly_scaling_variable

# The __main__ block is only for testing and not used by the pipeline, so it can be left as is or removed.
if __name__ == "__main__":
    print('This script is a module and should be imported, not run directly.')


# Adapted from Tod Stuber's code from 2020 by RW in 2026