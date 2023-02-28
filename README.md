# Pooled sorted CRISPR pipeline

The crispy pipeline ([pipeline_crispy.py](pipeline_crispy.py)) was build using CGAT-Core.

<a href="https://github.com/cgat-developers/cgat-core">
  <img src="https://github.com/cgat-developers/cgat-core/blob/master/docs/img/CGAT_logo.png" alt="CGAT-core" width="200">
</a>

CGAT-core is a workflow management system to build scalable data analysis pipelines. CGAT-core includes libraries and helper functions to enable researchers to quickly design and build computational workflows for the analysis of large-scale data-analysis.

[CGAT-core documentation](https://cgat-core.readthedocs.io/en/latest/ "CGAT-core read the docs")

CGAT-core is built upon the `ruffus` package. To learn more about ruffus, including the available task decorators, see [ruffus documentation](http://www.ruffus.org.uk/).

------------

### Installation

The instructions below will work for unix-like operating systems (linux or macOS). It may be possible to install and run the pipeline in windows, but I would strongly advise against trying, especially given that the optimal place to run the pipeline is on a HPC, not your personal computer.

The same installation instructions apply to run the pipeline locally and on the HPC.

See https://cgat-core.readthedocs.io/en/latest/getting_started/Installation.html
for further CGAT-core installation options if needed.

See https://github.com/mamba-org/mamba for instructions on how to install mamba.
Can alternatively only use conda to install the requirements, but mamba is quicker. Even with mamba,
you may find the installation takes 10 minutes to complete

Create a new conda environment and activate it, then install the required packages.
The R package poolr is not available via conda so must be installed via R.
```bash
mamba create --name crispr-pipeline python=3.10   
conda activate crispr-pipeline
mamba install -c conda-forge -c bioconda \
cgatcore fastqc bowtie multiqc samtools mageck pandoc \
r-base r-here r-desctools r-pheatmap r-dplyr r-tidyr r-ggplot2 r bioconductor-deseq2 r-rmarkdown r-optparse
R -e 'install.packages("poolr", repos="https://cran.ma.imperial.ac.uk/")'
```

------------

### Run the pipeline (locally)

#### Pipeline configuration
The pipeline configuration is controlled via the `pipeline.yml` config file. By default, the pipeline uses the config file in the same directory as the pipeline source code and we don't need to have a local configuration file. However, we can override paramaters using a config file in the working directory and it's good practise to have a local configuration file so we can at the very least read through it and confirm the configuration is correct.

First, create an empty directory for the pipeline output, change to this directory and run the following to generate a local pipeline configuration file.

```bash
python <PATH TO THIS REPOSITORY/pipelines/pipeline_crispr.py> config
```

The pipline needs the following input files:

1. Fastq files
2. Guide RNA sequences
3. Experimental conditions table
4. Statistical test design table(s)


#### 1. Fastq files
Create a subdirectory called `input` and add the fastq files. Better still, add softlinks to the fastqs, so the raw data can remain elsewhere on your filesystem, isolated from the pipeline files.


#### 2. Guide RNA sequences
Provide as a comma separated file with two columns, gRNA name and sequence. Examples of these files can be found in [guides](https://github.com/MRCToxBioinformatics/crispy/tree/main/guides). By default, the pipeline expects the filepath to be `guides.csv`. To specify another filepath, edit `pipeline.yml`.


#### 3. Experimental conditions table
For the QC plotting, we need to provide a map from samples names to the experimental conditions in a file called `experimental_design.tsv`.
This needs to include the columns shown in the example below.

```
sample	condition	sort	sort_by
H90_CTR_HIGH_1	Control	High	H90
H90_CTR_LOW_1	Control	Low	H90
H90_HS_HIGH_1	Heat stress	High	H90
H90_HS_LOW_1	Heat stress	Low	H90
PLASMID_LIBRARY	Control	Unsorted	Plasmid
SITA_CTR_HIGH_1	Control	High	SITA
SITA_CTR_LOW_1	Control	Low	SITA
SITA_HS_HIGH_1	Heat stress	High	SITA
SITA_HS_LOW_1	Heat stress	Low	SITA
unsorted_library_1	Control	Unsorted	Unsorted
```

Note that the samples should match the sample names as they are derived from the fastqs using the fastq_regex and fastq_pattern configuration parameters. See `pipeline.yml` for an example how the regex and pattern map from fastq filepaths to sample names.


#### 4. Statistical test design table(s)
To define which sets of samples should be compared, we need to provide a design file(s) which link the samples to the experimental conditions for the statistical test.
A basic example is shown below. This file should have two columns, sample and condition. The first condition is taken to be the reference, or control condition

```
sample,condition
H90_HS_LOW,low
H90_HS_HIGH,high
```

Note that the samples should match the sample names as they are derived from the fastqs using the fastq_regex and fastq_pattern configuration parameters. See `pipeline.yml` for an example how the regex and pattern map from fastq filepaths to sample names.

Statistical testing is performed for each design file and these *must* be named _design_xxx.csv_, where _xxx_ names the testing results.

By default, the pipeline will use the MAGeCK RRA algorithm. To use MAGeCK MLE, edit the `pipeline.yml` file (not yet supported)

#### Run the pipeline

You can now run the pipeline with:

```bash
python <PATH TO THIS REPOSITORY/pipelines/pipeline_crispr.py> --checksums=2 -p2 -v10 make full  --local
```

The options we are using are:

- `-p2`: Allow a maximum of 2 parallel tasks
- `--checksums = 2`: Re-run jobs where file timestamps are out of date, task failed, or task code has been updated. See [ruffus checksums](http://www.ruffus.org.uk/tutorials/new_tutorial/checkpointing.html?highlight=checksums)
- `-v10`: Verbose output
- `--local`: Don't submit jobs to a job scheduler, run them locally.
- `make full`: Run to the task 'full'


The final output of the pipeline are the multiqc report (in the `report.dir` directory), the gRNA counts (in the `quant.dir` directory) and the statistical test results (in the `mageck.dir` directory).

If you'd like to run the pipeline faster, change to `-p<NUMBER OF TASKS>`, replacing `<NUMBER OF TASKS>` with the number of parallel tasks you would like to run. However, you may not be able to run many parallel tasks on your desktop/laptop, in which case, use the HPC.

------------

### Run the pipeline (on the HPC)

First, log into the cambridge HPC (https://docs.hpc.cam.ac.uk/hpc/) and navigate to a suitable RDS directory. Then follow the same directions above to download the input data and install CGAT-core and the pipeline dependencies with mamba/conda.

Before we run the pipeline, we need to generate a `.cgat.yml` config file in the home directory so that CGAT-core knows how to interact with the HPC workload manager and submit jobs etc. The UoC HPC uses SLURM as the workload manager. There is a template [.cgat.file](https://github.com/MRCToxBioinformatics/Pipeline_examples/blob/main/CGATCore/.cgat.yml), which you can save to your home directory and update to provide the required details. Further details about configuring CGAT-core pipelines to run on clusters are here: https://cgat-core.readthedocs.io/en/latest/getting_started/Cluster_config.html.

After that, it's as simple as running the same commmand as previously, but without the `--local` and allowing more parallel processes.  Below, we allow a maximum of 10 tasks to be run in parallel. Note that this is different to the number of threads used in each individual task, which can be parameterised within the pipeline code and/or config file as required.

> &#x26a0;&#xfe0f; **Check you are in the appropriate working directory in your RDS before you run this command. You should have an input subdirectory with the neccessary input files and a file with gRNA sequences.**

```bash
nohup python <PATH TO THIS REPOSITORY/pipelines/pipeline_crispr.py>  --checksums=2 -p10 -v10 make full &
```

Note that we also include the command `nohup` at the start and `&` at the end of our command. `&` puts the command into the background and `nohup` ensures it will keep running even if we log out or lose connection with the HPC.

Your submitted jobs may not start immediately, but when they do, they will be much quicker and you can run many parallel tasks.  For the example data here, you are unlikely to see any significant run time benefit running on the HPC. However, if you want to analyse many samples, using the HPC becomes essential! This will also allow you to run resource-hungry jobs beyond your desktop/laptop specifications.

Above, we are running the python pipeline on the login node, which is OK for short pipelines which are not resource-hungry themselves. However, for best practise, the above command should be submitted to SLURM so that the python pipeline itself is run from a compute node, from which all jobs are subsequently submitted to SLURM. Only 1 CPU is needed to run the pipeline, though make sure to request sufficient time for the pipeline command to finish.

------------

## &#x26a0;&#xfe0f; Troubleshooting
On the HPC, you may get the following errors:

`conda Can't locate local/lib.pm in @INC (you may need to install the local::lib module)`

This can be resolved with:
`mamba install -c bioconda perl-local-lib`

------------
