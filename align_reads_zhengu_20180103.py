#!/user/bin/env python3

"""Align pair-end (PE) reads to reference genome sequences by using BWA.
 
 The following codes are inspired by the source code 'alignReads.py' 
 scripted by Martin Aryee, which can be found at GitHub
 https://github.com/aryeelab/guideseq/tree/master/guideseq.
 
 Usage::
     $ python3 align_reads_zhengu_20171211.py [params]
     
 :param --bwa_dir: directory to the software bwa
 :param --ref_index: location of the reference index files 
 :param --ref_seq: location of the FASTA file containing reference 
                   sequences                   
 :param --project_dir: directory to the corresponding project data 
 :param --read_file_name: root names of PE read files 
 :param --out_dir: directory to the alignment output files
 :param --threads: number of threads to be used when running bwa 
 :param --read_group: information of read groups
"""

__author__ = 'Maggie Ruimin Sun'
__version__ = '0.1'

import subprocess
import os
import logging
import time
import sys
import yaml

def setup_logger(name, log_file, formatter, level=logging.DEBUG):
    handler = logging.FileHandler(log_file)
    handler.setFormatter(formatter)
    
    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)
    
    return logger

def store_logs(log_dir):
    formatter_bwa_process = logging.Formatter("%(asctime)s;%(message)s")
    formatter_bwa_errors = logging.Formatter("%(asctime)s;%(levelname)s;%(message)s")
    logger_bwa_process = setup_logger('BWA Running Messages', 
                                      log_dir + '/bwa_process.log', 
                                      formatter_bwa_process)
    logger_bwa_errors = setup_logger('Errors & Warnings of BWA', 
                                    log_dir + '/bwa_errors.log',
                                    formatter_bwa_errors)
    return logger_bwa_process, logger_bwa_errors

def align_reads_bwa(bwa_dir, ref_fa_file, ref_index_name, read1, read2, 
                    out_file, n_threads, read_group, logger_bwa_process, 
                    logger_bwa_errors):
    out_dir = os.path.dirname(out_file)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    if not os.path.isfile(read1):
        logger_bwa_errors.error('%s does not exists!', read1)
        return
    
    # check the existence of each index file
    sample_alignment_dirs = {}
    index_files_extensions = ['.pac', '.amb', '.ann', '.bwt', '.sa']
    genome_indexed = True
    for extension in index_files_extensions:
        if not os.path.isfile(ref_index_name + extension):
            genome_indexed = False
            break
    if not genome_indexed:
        logger_bwa_process.info('Genome index files are not detected. '
                               'Running BWA to generate the required indices.')
        bwa_index_command = '{0} index -p {1} {2}'.format(
            bwa_dir, ref_index_name, ref_fa_file)
        logger_bwa_process.info(bwa_index_command)
        subprocess.call(bwa_index_command.split())
        logger_bwa_process.info('BWA genome index files are generated.')
    else:
        logger_bwa_process.info('BWA genome index files exist.')
    
    #Run BWA-MEM
    logger_bwa_process.info('Running paired end mapping.')
    bwa_align_command = '{0} mem -t {1} -M -R {2} {3} {4} {5} > {6}'.format(
        bwa_dir, n_threads, read_group, ref_index_name, read1, read2, out_file)
    logger_bwa_process.info(bwa_align_command)
    os.system(bwa_align_command)
    logger_bwa_process.info('Paired end mapping finished.')

def modify_sam_location(sam_file_org, sam_file_mod):
    sam_in = open(sam_file_org)
    sam_out = open(sam_file_mod, 'w')
    for row in sam_in:
        if row[0] == '@':
            sam_out.write(row)
            continue
        items = row.strip().split()
        if items[2] == '*':
            sam_out.write(row)
            continue
        chrom, start, stop = items[2].split('_')
        pos = int(v[3]) + int(start) - 1
        sam_out.write(v[0] + '\t' + v[1] + '\t' + chrom + '\t' + str(pos) 
                     + '\t' + '\t'.join(v[4:]) + '\n')
        sam_in.close()
        sam_out.close()

def main():
    with open('configure.yaml') as yamlfile:
        config = yaml.load(yamlfile)
    yamlfile.close()
    bwa_dir = config['alignment']['software']
    ref_index_name = config['alignment']['ref_index']
    ref_fa_file = config['reference']['fa_file']
    source = config['input_data']['input_dir']
    out_dir = source + config['alignment']['output_dir']
    log_dir = config['logging']
    logger_bwa_process, logger_bwa_errors = store_logs(log_dir)
    n_threads = config['alignment']['num_threads']
    
    time_start = time.time()
    library = sys.argv[1]
    #for library in config['input_data']['libraries']:
    sample_name = config['input_data']['libraries'][library]['sample_name']
    sample_id = config['input_data']['libraries'][library]['library_id']
    logger_bwa_process.info('Start alignment for data in {}.'.format(sample_name))
    filename1 = source + config['input_data']['libraries'][library]['location'] + sample_id + '-' + sample_name
    print(filename1)
    read1 = filename1 + config['input_data']['libraries'][library]['read1']
    read2 = filename1 + config['input_data']['libraries'][library]['read2']
    out_file = out_dir + sample_name + '_aligned.sam'
    read_group = '\'@RG\\tID:{0}\\tPL:Illumina\\tLB:YN\\tSM:{1}\''.format(sample_id, sample_name)
    align_reads_bwa(bwa_dir, ref_fa_file, ref_index_name, read1, read2, 
                    out_file, str(n_threads), read_group, logger_bwa_process, 
                    logger_bwa_errors)
    time_run = (time.time() - time_start) / 60
    logger_bwa_process.info('Finish alignment for library {0} after {1} min.'.format(library, str(time_run)))
    print('Finish alignment for library {0} after {1} min.'.format(library, str(time_run)))

if __name__ == '__main__':
    main()