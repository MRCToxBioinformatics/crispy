"""
====================
sorted CRISPR pipeline
====================


Overview
========
This pipeline performs QC, processing, quantification and statistical testing for a pooled sorted CRISPR screen


Input
-----
Reads are imported by placing files or linking to files in the :term:
`working directory`.
The following suffixes/file types are possible:

fastq.gz
   Single-end reads in fastq format.
Code
====
"""

###################################################
# load modules
###################################################

# import ruffus
from ruffus import transform, suffix, regex, merge, \
    follows, mkdir, originate, add_inputs, jobs_limit, split, \
    subdivide, formatter, collate


# import ruffus
from ruffus import *

# import useful standard python modules
import sys
import os
import inspect


# import re
# import shutil
# import sqlite3
# import glob

import numpy as np
import pandas as pd
# import cgatcore modules
#import cgatcore.experiment as E
from cgatcore import pipeline as P
import cgatcore.iotools as iotools

import ModulePipelineCrispy as Crispy

# load options from the config file
PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
    "pipeline.yml"])

conda_base_env = PARAMS['conda_base_env']

# Helper functions mapping tracks to conditions, etc
# determine the location of the input files (reads).
try:
    PARAMS["input"]
except KeyError:
    PARAMS["input"] = "."

# define input files. Here we allow single end only

READ1_SEQUENCEFILES = os.path.join(PARAMS["input"], PARAMS["basename_regex"])

SEQUENCEFILES_REGEX = regex(r"(.*\/)*(\S+).fastq.gz")


###################################################
# QC input
###################################################

@mkdir("fastqc.dir")
@transform(READ1_SEQUENCEFILES,
           SEQUENCEFILES_REGEX,
           r"fastqc.dir/\2.fastqc")
def runFastQC(infile, outfile):
    '''run FastQC on each input file.'''

    outdir=os.path.dirname(outfile)

    job_options = PARAMS['cluster_options']
    tmpdir = PARAMS['tmpdir']

    statement = '''
    fastqc
    -o %(outdir)s
    -d %(tmpdir)s
    %(infile)s
    >& %(outfile)s''' % locals()

    P.run(statement, job_condaenv=conda_base_env)


###################################################
# Build references
###################################################

@mkdir('reference.dir')
@originate('reference.dir/guides.fasta')
def makeGuidesFasta(outfile):
    with open(outfile, 'w') as outf:
        with open(PARAMS['reference_guides'], 'r') as inf:
            for line in inf:
                name, gene, sequence = line.strip().split(',')
                outf.write('>%s\n' % name)
                outf.write(sequence + '\n')


@mkdir('reference.dir')
@transform(makeGuidesFasta,
           suffix('.fasta'),
           '.1.ebwt')
def buildBowtieIndex(infile, outfile):
    '''Build bowtie index'''

    name_base = iotools.snip(outfile, '.1.ebwt')

    statement = '''
    bowtie-build %(infile)s %(name_base)s
    ''' % locals()

    P.run(statement, job_condaenv=conda_base_env)



###################################################
# Alignment
###################################################

@follows(buildBowtieIndex)
@mkdir("bowtie.dir")
@collate(READ1_SEQUENCEFILES,
         regex(os.path.join(PARAMS["input"], PARAMS["fastq_regex"])),
         r"bowtie.dir/%s.bowtie.bam" % PARAMS["fastq_pattern"],
         'reference.dir/guides.1.ebwt') # TSS: hardcoded to expected index name. Alternative, use add_inputs().
def runBowtie(infiles, outfile, index):
    '''Map reads to guides with bowtie'''

    infiles = ','.join(infiles)

    name_base = iotools.snip(index, '.1.ebwt')

    # P.run uses local variables including threads for job submission
    threads = PARAMS['bowtie_threads']

    bowtie_options = PARAMS['bowtie_options']

    tmp_file = P.get_temp_filename()


    statement = '''
    bowtie
    -x %(name_base)s
    %(bowtie_options)s
    -p %(threads)s
    -S
    %(infiles)s
    > %(tmp_file)s 2> %(outfile)s.stderr;
    samtools flagstat %(tmp_file)s > %(outfile)s.flagstat;
    samtools sort -o %(outfile)s -O BAM %(tmp_file)s;
    samtools index %(outfile)s;
    rm -f %(tmp_file)s

    ''' % locals()


    if PARAMS['cluster_queue_manager'] == "slurm":
        job_options = PARAMS['cluster_options'] + " -t 3:00:00"
    else:
        job_options = PARAMS['cluster_options']

    P.run(statement, job_condaenv=conda_base_env)


###################################################
# Alignment QC
###################################################

@transform(runBowtie,
           suffix('.bam'),
           '.bam.errors')
def countErrors(infile, outfile):
    ''' Count the number of error in the sequence reads
    relative to the expected guide sequences '''

    statement = '''
    samtools view  %(infile)s |
    cut -f14|
    sort |
    uniq -c|
    sed -e 's/^ *//g'  -e 's/[:| ]/\\t/'g |
    cut -f1,4 >
    %(outfile)s
    ''' % locals()

    P.run(statement)

@merge(countErrors,
       'bowtie.dir/all_errors.tsv')
def mergeErrorCounts(infiles, outfile):
    ''' merge the error counts across all samples'''
    sample_name = iotools.snip(os.path.basename(infiles[0]), '.bowtie.bam.errors')
    all_samples_df = pd.read_csv(infiles[0], sep='\t', header=None,
                                 names=(sample_name, 'errors')).set_index('errors')

    for infile in infiles[1:]:
        print(infile)
        print(all_samples_df)
        sample_name = iotools.snip(os.path.basename(infile), '.bowtie.bam.errors')
        samples_df = pd.read_csv(infile, sep='\t', header=None,
                                     names=(sample_name, 'errors')).set_index('errors')

        all_samples_df = all_samples_df.merge(samples_df, on='errors', how='outer')


    all_samples_df.fillna(0).to_csv(outfile, sep='\t')



###################################################
# Count guides
###################################################

@mkdir('quant.dir')
@transform(runBowtie,
           regex('bowtie.dir/(\S+).bowtie.bam'),
           r'quant.dir/\1.tsv')
def tallyGuides(infile, outfile):
    '''Count reads per gRNA'''
    # bam is sorted, so we can just use uniq -c on the contig column
    # -F 0x04 to exclude unmapped
    statement = '''
    samtools view %(infile)s -F 0x04 | cut -f3 | uniq -c | sed -e 's/^ *//;s/ /,/' >
    %(outfile)s
    ''' % locals()

    P.run(statement, job_condaenv=conda_base_env)

@merge(tallyGuides,
       'quant.dir/all_samples.tsv')
def mergeTallies(infiles, outfile):
    ''' merge the gRNA counts across all samples'''
    sample_name = iotools.snip(os.path.basename(infiles[0]), '.tsv')
    all_samples_df = pd.read_csv(infiles[0], sep=',', header=None,
                                 names=(sample_name, 'sgRNA')).set_index('sgRNA')

    for infile in infiles[1:]:
        sample_name = iotools.snip(os.path.basename(infile), '.tsv')
        samples_df = pd.read_csv(infile, sep=',', header=None,
                                 names=(sample_name, 'sgRNA')).set_index('sgRNA')
        all_samples_df = all_samples_df.merge(samples_df, on='sgRNA', how='outer')


    all_samples_df['gene'] = [x.split('_')[0] for x in all_samples_df.index]

    all_samples_df = all_samples_df.loc[
        :,['gene'] + [x for x in all_samples_df.columns.tolist() if x != 'gene']]

    all_samples_df.fillna(0).to_csv(outfile, sep='\t')


@transform(mergeTallies,
           suffix('.tsv'),
           '_norm.tsv')
def normaliseCounts(infile, outfile):

    job_options = PARAMS['cluster_options'] + " -t 0:20:00"
    job_condaenv=PARAMS['conda_base_env']

    Crispy.normaliseCounts(infile, outfile, submit=True, job_options=job_options)

###################################################
# Report
###################################################

@mkdir('report.dir')
@merge((runFastQC, runBowtie), 'report.dir/multiqc_report.html')
def runMultiQC(infiles, outfile):
    ''' Run multiqc to generate basic QC report'''

    # we'll just take the directories of the infiles to make life easy
    indirs = ' '.join(set([os.path.dirname(x) for x in infiles]))

    outdir = os.path.dirname(outfile)

    job_options = PARAMS['cluster_options']

    statement = 'multiqc %(indirs)s -o %(outdir)s' % locals()

    # run tasks locally
    P.run(statement,
          job_condaenv=conda_base_env,
          without_cluster=False)

###################################################
# Re-sampling
###################################################

@mkdir('resampled_dummy_files')
@originate(['resampled_dummy_files/%s' % round(x, 3) for x in
            np.arange(PARAMS['resample_min'], PARAMS['resample_max']+PARAMS['resample_step'], PARAMS['resample_step'])])
def create_resample_dummies(output_file):
    '''make empty files names for each level of
    resampling to use to define downstream task'''
    with open(output_file, "w"):
        pass


@mkdir('resampled_quant.dir')
@follows(mergeTallies)
@transform(create_resample_dummies,
           regex('resampled_dummy_files/(\S+)'),
           add_inputs(mergeTallies),
           r'resampled_quant.dir/\1_all_samples.tsv')
def resampleTallies(infiles, outfile):
    ''' Resample the counts table to simulate different sequencing depths '''
    dummy_infile, infile = infiles

    sample_frac = np.float(os.path.basename(dummy_infile))

    job_options = PARAMS['cluster_options'] + " -t 0:20:00"
    job_condaenv=PARAMS['conda_base_env']

    Crispy.resampleTallies(infile, sample_frac, outfile, submit=True, job_options=job_options)

###################################################
# QC
###################################################
@mkdir('qc_plots.dir')
@follows(mergeErrorCounts, mergeTallies, resampleTallies)
@originate('QC_plotting.html')
def runCrispyQC(outfile):

    this_filename = inspect.getframeinfo(inspect.currentframe()).filename
    this_dir     = os.path.dirname(os.path.abspath(this_filename))

    notebook_path = os.path.join(this_dir, 'R', 'QC_plotting.Rmd')

    if PARAMS['cluster_queue_manager'] == "slurm":
        job_options = PARAMS['cluster_options'] + " -t 1:00:00"
    else:
        job_options = PARAMS['cluster_options']

    statement = '''
    cp %(notebook_path)s . ;
    Rscript -e "rmarkdown::render('QC_plotting.Rmd')"
    ''' % locals()

    P.run(statement)


###################################################
# statistical tests
###################################################
if PARAMS['mageck_method'].lower() == 'rra':
    @mkdir('mageck.dir')
    @transform('design_*.csv',
               regex('design_(\S+).csv'),
               add_inputs(mergeTallies),
               r'mageck.dir/\1/\1.gene_summary.txt')
    def runMAGeCK(infiles, outfile):
        ''' run MAGeCK RRA to identify enriched/depleted '''

        design_inf, counts = infiles
        counts = os.path.abspath(counts)

        outfile_base = P.snip(os.path.basename(outfile), '.gene_summary.txt')

        design = pd.read_table(design_inf, sep=',')

        if not design.columns.tolist() == ['sample', 'condition']:
            raise ValueError(
                '''Unexpected design table format for file %(design_inf)s
                Table should only have sample and condition columns''')

        control_condition = design.condition[0]
        treatment_condition = [x for x in design.condition if x != control_condition][0]
        control_samples = design[design.condition==control_condition]['sample'].tolist()
        treatment_samples = design[design.condition==treatment_condition]['sample'].tolist()

        control_samples = ','.join(control_samples)
        treatment_samples = ','.join(treatment_samples)

        statement = '''
        cd mageck.dir/%(outfile_base)s;
        mageck test
        -k %(counts)s
        -c %(control_samples)s
        -t %(treatment_samples)s
        -n %(outfile_base)s
        ''' % locals()

        P.run(statement)


if PARAMS['mageck_method'].lower() == mle:
    @mkdir('mageck.dir')
    @transform(P.asList(PARAMS['mageck_designs']),
               regex('design_(\S+).txt'),
               add_inputs(mergeTallies),
               r'mageck.dir/\1/\1.gene_summary.txt')
    def runMAGeCK(infiles, outfile):
        ''' run MAGeCK MLE to identify enriched/depleted '''

        design_inf, counts = infiles
        counts = os.path.abspath(counts)

        outfile_base = P.snip(os.path.basename(outfile), '.gene_summary.txt')

        job_threads = PARAMS['mageck_mle_threads']

        statement = '''
        cd mageck.dir/%(outfile_base)s;
        mageck mle
        --norm-method none
        -k %(counts)s
        -d %(design_inf)s
        -n %(outfile_base)s
        --threads %(job_threads)s 
        ''' % locals()

        P.run(statement)

else:
    raise ValueError('mageck_method must be "rra" or "mle"')



###################################################
# targets
###################################################



# full_quant = run everything except the QC notebook
@follows(mergeErrorCounts,
         normaliseCounts)
def full_quant():
    pass


# qc = run all qc
@follows(runMultiQC,
         runCrispyQC)
def qc():
    pass

# full = run it all!
@follows(runMAGeCK,
         runCrispyQC)
def full():
    pass

###################################################
# Making pipline command-line friendly
###################################################

# Facilitate command line parsing
def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
