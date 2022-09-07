#!/usr/bin/env python

__version__ = "0.0.1"

import os
import sys
import re
import argparse
import textwrap
import json
import math
import numpy as np

class Quality_Scaling:

    def fastq_scaling(self, genome_size, read1_reads_gt_q30=None, read2_reads_gt_q30=None, sampling_size=None, read1_total_read_count=None, read2_total_read_count=None, read1_read_average=None, read2_read_average=None):
        '''
        seqval changed to fastq_scaling_variable
        '''
        if not read2_total_read_count:
            genome_coverage_depth = int(read1_total_read_count) * 240 / (int(genome_size))
            fastq_scaling_variable = 0
        else:
            genome_coverage_depth = int(read1_total_read_count + read2_total_read_count) * 240 / (int(genome_size))

            seqcor=[[1.0000000,0.4412253,0.9569653,0.4266978,0.3393063],
            [0.4412253,1.0000000,0.4446494,0.9527558,0.2464773],
            [0.9569653,0.4446494,1.0000000,0.4330931,0.3315736],
            [0.4266978,0.9527558,0.4330931,1.0000000,0.2367824],
            [0.3393063,0.2464773,0.3315736,0.2367824,1.0000000]]

            q301 = (read1_reads_gt_q30/sampling_size) * 100
            q302 = (read2_reads_gt_q30/sampling_size) * 100
            x = [(float(read1_read_average)-30), (float(read2_read_average)-30), (float(q301)-80), (float(q302)-80), genome_coverage_depth]
            mx = np.matrix(x)
            mseqcor = np.matrix(seqcor)
            mxt = mx.transpose()
            seqv = math.sqrt(mx * mseqcor * mxt)
            fastq_scaling_variable = (seqv*800/90) + 1180

            if fastq_scaling_variable > 1980:
                fastq_scaling_variable = 1980
            elif fastq_scaling_variable < 1180:
                fastq_scaling_variable = 1180

            return fastq_scaling_variable, genome_coverage_depth

    def assembly_scaling(self, longest_contig=None, greater_one_kb_count=None, contig_count=None, n50=None, l50=None, rgl=None, stat_total_contig_lengths=None, stat_contig_count=None, genome_coverage_depth=None):
        '''
        assembleval changed to assembly_scaling_variable
        '''
        
        assemblecor=[[1.00000000,0.4640352,0.39099636,-0.703376083,-0.074391092,-0.27033841,-0.26923388,-0.36627871],
        [0.46403525,1.0000000,0.75927035,-0.220017666,-0.101805557,-0.56550423,0.27504920,-0.70388508],
        [0.39099636,0.7592704,1.00000000,-0.170697061,-0.054252944,-0.72295517,0.29979708,-0.81000716],
        [-0.70337608,-0.2200177,-0.17069706,1.000000000,0.003833981,0.09220126,0.46237217,0.16084638],
        [-0.07439109,-0.1018056,-0.05425294,0.003833981,1.000000000,0.05500119,-0.04451447,0.01902384],
        [-0.27033841,-0.5655042,-0.72295517,0.092201260,0.055001186,1.00000000,-0.28768966,0.55006987],
        [-0.26923388,0.2750492,0.29979708,0.462372175,-0.044514471,-0.28768966,1.00000000,-0.27338410],
        [-0.36627871,-0.7038851,-0.81000716,0.160846384,0.019023843,0.55006987,-0.27338410,1.00000000]]

        massemblecor = np.matrix(assemblecor)
        perlongexpect = float(longest_contig)/(int(genome_coverage_depth)*1000000)*100
        perlongscaff = float(greater_one_kb_count)/float(contig_count)*100
        pern50expect = int(n50) / (int(stat_total_contig_lengths) * 1000000) * 100
        perdiffexpect = abs(100-(float(stat_total_contig_lengths)/(int(stat_total_contig_lengths)*10000))) 

        quality_scale_list = [contig_count, greater_one_kb_count, l50, rgl, perdiffexpect, perlongexpect, perlongscaff, pern50expect] 
        ma = np.matrix(quality_scale_list)
        mat = ma.transpose()
        assemblev = 1000/(math.sqrt(ma*massemblecor*mat))
        assembly_scaling_variable = (assemblev*800/7) + 1180

        if assembly_scaling_variable > 1980 :
            assembly_scaling_variable = 1980
        elif assembly_scaling_variable < 1180:
            assembly_scaling_variable = 1180
        
        return assembly_scaling_variable
    
if __name__ == "__main__":
    '''
    '''
    print(f'### FAKE DATA BEING USED, JUST FOR TESTING, MUST USE AS PYTHON MODULE WHEN FOR REAL')
    sequence_score = Sequence_Score(mlst_scheme='senterica', total_contig_lengths=4703983, read1_total_read_count=200000, read2_total_read_count=200000, read1_reads_gt_q30=9409, read2_reads_gt_q30=8861, read1_read_average=36.179729141849165, read2_read_average=35.153855272253125)

    print(f'done')
# Created March 2020 by Tod Stuber
