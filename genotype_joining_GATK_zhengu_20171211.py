#!/user/bin/env python3

"""Combine all per-sample GCVFs to produce a set of joint-called SNPs and indels 
ready for filtering. The GATK tool GenotypeGCVFs is applied here.
"""
__author__ = 'Maggie Ruimin Sun'
__version__ = '0.1'

import subprocess
import os
import logging
import time
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
    formatter_gatk_process = logging.Formatter("%(asctime)s;%(message)s")
    formatter_gatk_errors = logging.Formatter("%(asctime)s;%(levelname)s;                                              %(message)s")
    logger_gatk_process = setup_logger('GATK Running Messages', 
                                      log_dir + '/gatk_process.log',
                                      formatter_gatk_process)
    logger_gatk_errors = setup_logger('Errors & Warnings of GATK', 
                                     log_dir + '/gatk_errors.log', 
                                     formatter_gatk_errors)
    
    
    return (logger_gatk_process, logger_gatk_errors)

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


def gather_gvcfs(gatk_dir, ref_seq, gvcf_list, variants_dir, out_name, thres_call, 
                 logger_gatk_process, logger_gatk_errors):
    if not os.path.exists(variants_dir):
        logger_gatk_errors.error('The directory %s does not exist!', variants_dir)
        print('ERROR:: The directory %s does not exist! Please check the correctness of '
             'input. Make sure the GVCF files are contained in the directory.')
        return 1
    
    command_joint = ('java -jar {0} -T GenotypeGVCFs -R {1} -stand_call_conf {2} --variant ' 
                     + ' --variant '.join(gvcf_per_sample for gvcf_per_sample in gvcf_list)
                     + ' -o {3}').format(gatk_dir, ref_seq, thres_call, 
                                        variants_dir+out_name+'.vcf')
    returncode_joint = run_shell_command(command_joint, logger_gatk_process, logger_gatk_errors)
    
    if not returncode_joint == 0:
        logger_gatk_errors.error('GVCFs joining returns non-zero value.')
        print('GVCFs joining returns non-zero value.')
        return 1
    
    return 0

def main():
    with open('configure.yaml') as yamlfile:
        config = yaml.load(yamlfile)
    yamlfile.close()
    gatk_dir = config['snv_calling']['software']
    ref_seq = config['reference']['fa_file']
    variants_dir = config['input_data']['input_dir'] + config['snv_calling']['output_dir']
    log_dir = config['logging']
    out_name = config['genotype_joining']['output_name']
    thres_call = config['snv_calling']['threshold']
    sample_list = []
    for library in config['input_data']['libraries']:
        sample_list.append(config['input_data']['libraries'][library]['sample_name'])
    gvcf_list = [variants_dir + x + '.g.vcf' for x in sample_list]
    
    logger_gatk_process, logger_gatk_errors = store_logs(log_dir)
    returncode = gather_gvcfs(gatk_dir, ref_seq, gvcf_list, variants_dir, out_name, thres_call,
                              logger_gatk_process, logger_gatk_errors)

if __name__ == '__main__':
    main()