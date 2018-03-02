#!/user/bin/env python3

"""Discover germline single nucleotide variants in each sample 
by using GATK & Picard. 

Three Stages are included in this script:
 
 1 Sorting the alignment SAM file using Picard. 
 2 Recalibrate base quality scores using GATK BaseRecalibrator. 
 3 Discovery variants using GATK HaplotypeCaller.
 
 Usage::
     $ python3 germline_variant_calling_GATK_zhengu_20171211.py [params]
     
  :param --picard_dir: directory to the software picard.jar 
  :param --gatk_dir: directory to the software GenomeAnalysisTK.jar 
  :param --ref_seq: directory to the FASTA file containing reference seuqences
  :param --align_dir: directory to the alignment SAM and BAM files 
  :param --sample_name: the root name of the alignment files relevant to each 
                        sample 
 :param --thres_call: the minimum phred-scaled confidence threshold at which
                      variants should be called.(10.0 by default）
 :param --out_dir: directory to the output gVCF files.
"""
__author__ = 'Maggie Ruimin Sun'
__version__ = '0.1'

import subprocess
import os
import logging
import time
import sys
import shlex
import yaml

def setup_logger(name, log_file, formatter, level=logging.DEBUG):
    handler = logging.FileHandler(log_file)
    handler.setFormatter(formatter)
    
    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)
    
    return logger

def store_logs(log_dir):
    formatter_picard_process = logging.Formatter("%(asctime)s;%(message)s")
    formatter_picard_errors = logging.Formatter("%(asctime)s;%(levelname)s;%(message)s")
    logger_picard_process = setup_logger('Picard Running Messages', 
                                      log_dir + '/picard_process.log', 
                                      formatter_picard_process)
    logger_picard_errors = setup_logger('Errors & Warnings of Picard', 
                                    log_dir + '/picard_errors.log',
                                    formatter_picard_errors)
    
    formatter_gatk_process = logging.Formatter("%(asctime)s;%(message)s")
    formatter_gatk_errors = logging.Formatter("%(asctime)s;%(levelname)s;%(message)s")
    logger_gatk_process = setup_logger('GATK Running Messages', 
                                      log_dir + '/gatk_process.log',
                                      formatter_gatk_process)
    logger_gatk_errors = setup_logger('Errors & Warnings of GATK', 
                                     log_dir + '/gatk_errors.log', 
                                     formatter_gatk_errors)
    
    
    return (logger_picard_process, logger_picard_errors, 
            logger_gatk_process, logger_gatk_errors)

def run_shell_command(command_line, logger_process, logger_errors):
    command_line_args = shlex.split(command_line)
    logger_process.info(command_line)
    command_line_process = subprocess.run(command_line_args,
                                            stdout=subprocess.PIPE, 
                                            stderr=subprocess.STDOUT)
    if command_line_process.stdout is not None:
        for stdout_line in command_line_process.stdout.splitlines():
            logger_process.info(stdout_line)
    if command_line_process.stderr is not None:
        for stderr_line in command_line_process.stderr.splitlines():
            logger_errors.error(stderr_line)
    return command_line_process.returncode

def sort_sam_picard(picard_dir, input_sam, output_bam, 
                    logger_picard_process, logger_picard_errors):
    out_dir = os.path.dirname(output_bam)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        
    if not os.path.isfile(input_sam):
        logger_picard_errors.error('%s does not exists!', input_sam)
        print('%s does not exists!', input_sam)
        return 1
    
    command_sort = 'java -jar {0} SortSam INPUT={1} OUTPUT={2} SORT_ORDER=coordinate'.format(
        picard_dir, input_sam, output_bam)
    command_index = 'java -jar {0} BuildBamIndex INPUT={1}'.format(
        picard_dir, output_bam)
    
    #logger_picard_process.info(command_sort)
    returncode_sort = run_shell_command(command_sort, logger_picard_process, 
                                        logger_picard_errors)
    logger_picard_process.info('Finished sorting SAM with return value %f.', returncode_sort)
    if not returncode_sort == 0:
        logger_picard_errors.error('Picard sorting returns non-zero value.')
        print('Picard sorting returns non-zero value! Check the input format '
              'and/or the logging files.')
        return 1

    logger_picard_process.info('Finished Picard sorting step.')
    print('Finished Picard sorting step.')
    returncode_index = run_shell_command(command_index, logger_picard_process, 
                                         logger_picard_errors)
    if not returncode_index == 0:
        logger_picard_errors.error('Picard indexing returns non-zero value.')
        print('Picard indexing returns non-zero value! Check the input format '
              'and/or the logging files.')
        return 1
    
    return 0

def recalibrate_base_quality_scores(gatk_dir, ref_fa_file, sorted_bam, 
                                    knownsites, align_dir, logger_gatk_process,
                                    logger_gatk_errors):
    
    """kownsites is a list-type argument, containing a set of known variant vcf files"""
    
    if not os.path.exists(align_dir):
        os.makedirs(align_dir)
    prefix_name = sorted_bam.split('_aligned')[0]
    sorted_bam = align_dir + sorted_bam
    if not os.path.isfile(sorted_bam):
        logger_gatk_errors.error('%s does not exists!', sorted_bam)
        print('%s does not exists!', sorted_bam)
        return 1
        
    # Step 1: analyze patterns of covariation in the sequence dataset
    command_covariation_analysis = 'java -jar {0} -T BaseRecalibrator \
    -R {1} -I {2}'.format(gatk_dir, ref_fa_file, sorted_bam)
    for knownsite in knownsites:
        command_covariation_analysis += (' -knownSites ' + knownsite)
    output_recal_table = align_dir + prefix_name + '_recal.table'
    command_covariation_analysis += (' -o ' + output_recal_table)
    returncode_covariation_analysis = run_shell_command(command_covariation_analysis, 
                                                        logger_gatk_process, 
                                                        logger_gatk_errors)
    if not returncode_covariation_analysis == 0:
        logger_gatk_errors.error('BQSR failed at covariation analysis stage.')
        print('BQSR failed at covariation analysis stage.')
        return 1
    
    # Step 2: do a second pass to analyze covariation remaining after recalibration
    command_covariation_analysis2 = 'java -jar {0} -T BaseRecalibrator -R {1} \
    -I {2}'.format(gatk_dir, ref_fa_file, sorted_bam)
    for knownsite in knownsites:
        command_covariation_analysis2 += (' -knownSites ' + knownsite)
    output_post_recal_table = align_dir + prefix_name + '_recal_post.table'
    command_covariation_analysis2 += ' -BQSR {0} -o {1}'.format(output_recal_table, 
                                                                output_post_recal_table)
    returncode_covariation_analysis2 = run_shell_command(command_covariation_analysis2, 
                                                         logger_gatk_process,
                                                         logger_gatk_errors)
    if not returncode_covariation_analysis2 == 0:
        logger_gatk_errors.error('BQSR failed at the second pass for covariation analysis.')
        print('BQSR failed at the second pass for covariation analysis.')
        return 1
    
    # Step 3: Generate before/after plots
    command_plot = 'java -jar {0} -T AnalyzeCovariates -R {1} -before {2} -after {3} \
    -plots {4}'.format(gatk_dir, ref_fa_file, output_recal_table, output_post_recal_table, 
                      align_dir + prefix_name + '_recal_plots.pdf')
    returncode_plot = run_shell_command(command_plot, logger_gatk_process, logger_gatk_errors)
    if not returncode_plot == 0:
        logger_gatk_errors.error('BQSR failed at the ploting stage.')
        print('BQSR failed at the ploting stage.')
        return 1
    
    # Step 4: Apply the recalibration to the sorted alignment data
    command_recalibrate = 'java -jar {0} -T PrintReads -R {1} -I {2} -BQSR {3} -o {4}'.format(
    gatk_dir, ref_fa_file, sorted_bam, output_recal_table, align_dir + prefix_name + '_recal.bam')
    returncode_recalibrate = run_shell_command(command_recalibrate, logger_gatk_process, 
                                               logger_gatk_errors)
    if not returncode_recalibrate == 0:
        logger_gatk_errors.error('BQSR failed at the final application stage.')
        print('BQSR failed at the final application stage.')
        return 1
    return 0

def call_variants(gatk_dir, ref_fa_file, recal_bam, threshold_call, out_dir, 
                  logger_gatk_process, logger_gatk_errors):
    
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    if not os.path.isfile(recal_bam):
        logger_gatk_errors.error('%s does not exists!', recal_bam)
        print('%s does not exists!', recal_bam)
        return 1
    samp_name = recal_bam.split('/')[-1].split('_recal.bam')[0]
    
    command_call_var = 'java -jar {0} -T HaplotypeCaller -R {1} -I {2} --genotyping_mode DISCOVERY \
    -stand_call_conf {3} --emitRefConfidence GVCF -o {4}'.format(
        gatk_dir, ref_fa_file, recal_bam, threshold_call, 
        out_dir + samp_name + '_raw_variants.g.vcf')
    returncode_call_variants = run_shell_command(command_call_var, logger_gatk_process, 
                                                 logger_gatk_errors)
    if not returncode_call_variants == 0:
        logger_gatk_errors.error('HaplotypeCaller failed.')
        print('HaplotypeCaller failed.')
        return 1

    return 0

def main():
    with open('configure.yaml') as yamlfile:
        config = yaml.load(yamlfile)
    yamlfile.close()
    picard_dir = config['sorting']['software']
    gatk_dir = config['snv_calling']['software']
    ref_seq = config['reference']['fa_file']
    align_dir = config['input_data']['input_dir'] + config['alignment']['output_dir']
    thres_call = config['snv_calling']['threshold']
    out_dir = config['snv_calling']['output_dir']
    knownsites = []
    for db in config['snv_calling']['knownsites']:
        knownsites.append(config['snv_calling']['knownsites'][db])
    log_dir = config['logging']
    
    (logger_picard_process, logger_picard_errors, logger_gatk_process, 
     logger_gatk_errors) = store_logs(log_dir)
    
    library = sys.argv[1]
    sample_name = config['input_data']['libraries'][library]['sample_name']
    input_sam = align_dir + sample_name + '_aligned.sam'
    sorted_bam = align_dir + sample_name + '_aligned_sorted.bam'
    returncode_picard = sort_sam_picard(picard_dir, input_sam, sorted_bam, 
                                        logger_picard_process, logger_picard_errors)
    
    if not returncode_picard == 0:
        return 1
    sorted_bam = sample_name + '_aligned_sorted.bam'
    returncode_bqsr = recalibrate_base_quality_scores(gatk_dir, ref_seq, sorted_bam, 
                                                      knownsites, align_dir, 
                                                      logger_gatk_process, logger_gatk_errors)
    if not returncode_bqsr == 0:
        return 2
    recal_bam = align_dir + sample_name + '_recal.bam'
    returncode_var = call_variants(gatk_dir, ref_seq, recal_bam, thres_call, out_dir, 
                                  logger_gatk_process, logger_gatk_errors)
    if not returncode_var == 0:
        return 3

if __name__ == '__main__':
    main()