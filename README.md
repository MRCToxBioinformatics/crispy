# Pooled sorted CRISPR pipeline

The pipeline was build using CGAT-Core

<a href="https://github.com/cgat-developers/cgat-core">
  <img src="https://github.com/cgat-developers/cgat-core/blob/master/docs/img/CGAT_logo.png" alt="CGAT-core" width="200">
</a>

CGAT-core is a workflow management system to build scalable data analysis pipelines. CGAT-core includes libraries and helper functions to enable researchers to quickly design and build computational workflows for the analysis of large-scale data-analysis.

[CGAT-core documentation](https://cgat-core.readthedocs.io/en/latest/ "CGAT-core read the docs")

CGAT-core is built upon the `ruffus` package. To learn more about ruffus, including the available task decorators, see [ruffus documentation](http://www.ruffus.org.uk/).


### Installation

The instructions below will work for unix-like operating systems (linux or macOS). It may be possible to install and run the pipeline in windows, but I would strongly advise against trying, given that the optimal place to run the pipeline is on a HPC, not your personal computer.

The same installation instructions apply to run the pipeline locally and on the HPC.
We're using mamba below. See https://cgat-core.readthedocs.io/en/latest/getting_started/Installation.html
for further CGAT-core installation options if needed.

See https://github.com/mamba-org/mamba for instructions on how to install mamba.
Can alternatively only use conda, but mamba is quicker for installation. Even with mamba,
you may find the final command takes a few minutes to complete

Create a new conda environment and activate it, then install the required packages
```bash
mamba create --name crispr-pipeline python=3.10   
conda activate crispr-pipeline
mamba install -c conda-forge -c bioconda \
cgatcore fastqc bowtie multiqc samtools mageck
```

### Run the pipeline (locally)
First, create an empty directory for the pipeline output, change to this directory and run the following to generate a local pipeline configuration file (`pipeline.yml`).

```bash
python <PATH TO THIS REPOSITORY/pipelines/pipeline_crispr.py> config
```

Now, create a subdirectory called `input` and add the fastq files. Better still, add softlinks to the fastqs, so the raw data can remain elsewhere on your filesystem, isolated from the pipeline files.

The only other input file the pipeline needs is a comma separated file with two columns, gRNA name and sequence. Examples of these files can be found in [guides](https://github.com/MRCToxBioinformatics/crispy/tree/main/guides). By default, the pipeline expects the filepath to be `guides.csv`. To specify another filepath, edit `pipeline.yml`.

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

If you'd like to run the pipeline faster, change to `-p<NUMBER OF TASKS>`, replacing `<NUMBER OF TASKS>` with the number of parallel tasks you would like to run. However, you may not be able to run many parallel tasks on your desktop/laptop.... which brings us to...

### Run the pipeline (on the HPC)

First, log into the cambridge HPC (https://docs.hpc.cam.ac.uk/hpc/) and navigate to a suitable RDS directory. Then follow the same directions above to download the input data and install CGAT-core and the pipeline dependencies with mamba/conda.

Before we run the pipeline, we need to generate a `.cgat.yml` config file in the home directory so that CGAT-core knows how to interact with the HPC workload manager and submit jobs etc. The UoC HPC uses SLURM as the workload manager. There is a template [.cgat.file](https://github.com/MRCToxBioinformatics/Pipeline_examples/blob/main/CGATCore/.cgat.yml), which you can save to your home directory and update to provide the required details. Further details about configuring CGAT-core pipelines to run on clusters are here: https://cgat-core.readthedocs.io/en/latest/getting_started/Cluster_config.html.

After that, it's as simple as running the same commmand as previously, but without the `--local` and allowing more parallel processes.  Below, we allow a maximum of 10 tasks to be run in parallel. Note that this is different to the number of threads used in each individual task, which can be parameterised within the pipeline code and/or config file as required.

The pipeline configuration is controlled via the `pipeline.yml` config file. By default, the pipeline uses the config file in the same directory as the pipeline source code. We can also override paramaters using a config file in the working directory.

> &#x26a0;&#xfe0f; **Check you are in the appropriate working directory in your RDS before you run this command. You should have an input subdirectory with the neccessary input files and a file with gRNA sequences.**

```bash
nohup python <PATH TO THIS REPOSITORY/pipelines/pipeline_crispr.py>  --checksums=2 -p10 -v10 make full &
```

Note that we also include the command `nohup` at the start and `&` at the end of our command. `&` puts the command into the background and `nohup` ensures it will keep running even if we log out or lose connection with the HPC.

Your submitted jobs may not start immediately, but when they do, they will be much quicker and you can run many parallel tasks.  For the example data here, you are unlikely to see any significant run time benefit running on the HPC. However, if you want to analyse many samples, using the HPC becomes essential! This will also allow you to run resource-hungry jobs beyond your desktop/laptop specifications.

Above, we are running the python pipeline on the login node, which is OK for short pipelines which are not resource-hungry themselves. However, for best practise, the above command should be submitted to SLURM so that the python pipeline itself is run from a compute node, from which all jobs are subsequently submitted to SLURM. Only 1 CPU is needed to run the pipeline, though make sure to request sufficient time for the pipeline command to finish.

## &#x26a0;&#xfe0f; Troubleshooting
On the HPC, you may get the following errors:

------------

`conda Can't locate local/lib.pm in @INC (you may need to install the local::lib module)`

This can be resolved with:
`mamba install -c bioconda perl-local-lib`

------------
