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
# import re
# import shutil
# import sqlite3
# import glob

import pandas as pd
# import cgatcore modules
#import cgatcore.experiment as E
from cgatcore import pipeline as P
import cgatcore.iotools as iotools


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

READ1_SEQUENCEFILES = os.path.join(PARAMS["input"], '*_R1_001.fastq.gz') # TSS: Generalise this. Config file?

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

@mkdir("bowtie.dir")
@collate(READ1_SEQUENCEFILES,
         regex('.*\/(\S+)_S\d+_L00[1-4]_R1_001.fastq.gz'), # TSS: Generalise this. Config file?
         add_inputs(buildBowtieIndex),
         r"bowtie.dir/\1.bowtie.bam")
def runBowtie(infiles, outfile):
    '''Map reads to guides with bowtie'''

    index = infiles[0][1]
    infiles = ','.join([x[0] for x in infiles])

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


    job_options = PARAMS['cluster_options'] + ' -t 1:00:00'
    P.run(statement, job_condaenv=conda_base_env)


###################################################
# Alignment QC
###################################################

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
    statement = '''
    samtools view %(infile)s | cut -f3 | uniq -c | sed -e 's/^ *//;s/ /,/' >
    %(outfile)s
    ''' % locals()

    P.run(statement, job_condaenv=conda_base_env)

@merge(tallyGuides,
       'quant.dir/all_samples.tsv')
def mergeTallies(infiles, outfile):
    ''' merge the gRNA counts across all samples'''
    sample_name = iotools.snip(os.path.basename(infiles[0]), '.tsv')
    all_samples_df = pd.read_csv(infiles[0], sep=',', header=None, names=(sample_name, 'guide')).set_index('guide')

    for infile in infiles[1:]:
        sample_name = iotools.snip(os.path.basename(infile), '.tsv')
        samples_df = pd.read_csv(infile, sep=',', header=None, names=(sample_name, 'guide')).set_index('guide')
        all_samples_df = all_samples_df.merge(samples_df, on='guide', how='outer')

    all_samples_df.fillna(0).to_csv(outfile, sep='\t')

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
# targets
###################################################



# full = run it all!
@follows(runMultiQC,
         runBowtie)
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
