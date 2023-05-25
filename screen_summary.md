Tom Smith 
22 FEB 2023

This file documents the pooled sorted CRISPR screens to date and which samples passed QC according to Tom and Roberto Campalastri.

Currently, the input files that passed QC for all screens can be found in `input` subdirectores in the directories in `rds/rds-mrc_tox-XUr6B1Jhndg/tss38_09_CRISPR_da442_ap2260/pipeline_runs`.
These files are listed below for each screen

**POLIII - HAP1**
The following samples from the first run (GW2) could not be used:
HAP1_LIB2_R2
POL3_HS_LOW2
POL3_HS_HIGH2
The following samples from the repeat run (GW5) can be used in their place:
POL3_HS_LOW2
POL3_HS_HIGH2

_Input_files_
HAP1_LIB2_R1_001.fastq.gz
POL3_CTR_HIGH1_R1_001.fastq.gz
POL3_CTR_HIGH2_R1_001.fastq.gz
POL3_CTR_LOW1_R1_001.fastq.gz
POL3_CTR_LOW2_R1_001.fastq.gz
POL3_HS_HIGH1_R1_001.fastq.gz
POL3_HS_HIGH2_R1_001.fastq.gz
POL3_HS_LOW1_R1_001.fastq.gz
POL3_HS_LOW2_R1_001.fastq.gz


**SUMO - MOLM13**
No replicates, except for the unsorted library. Samples are from GW3:
MOLM13_LIB1_R1_001.fastq.gz
MOLM13_LIB2_R1_001.fastq.gz
sumo_CTR_HIGH1_R1_001.fastq.gz
sumo_CTR_LOW1_R1_001.fastq.gz
sumo_HS_HIGH1_R1_001.fastq.gz
sumo_HS_LOW1_R1_001.fastq.gz

**NELFE - HAP1**
NELFE_HS_LOW2 rejected (High Gini coefficient, poor correlation)

HAP1_LIB2_R1_001.fastq.gz
NELFE_CTR_HIGH1_R1_001.fastq.gz
NELFE_CTR_HIGH2_R1_001.fastq.gz
NELFE_CTR_LOW1_R1_001.fastq.gz
NELFE_CTR_LOW2_R1_001.fastq.gz
NELFE_HS_HIGH1_R1_001.fastq.gz
NELFE_HS_HIGH2_R1_001.fastq.gz
NELFE_HS_LOW1_R1_001.fastq.gz

** UBIK48 - **
Replicate 1 samples and unsorted libraries 1 & 2 are from GW3. The rest are from GW4. Clear batch
effect across the sequencing runs	

MOLM13_LIB1_R1_001.fastq.gz
MOLM13_LIB2_R1_001.fastq.gz
MOLM13_LIB3_R1_001.fastq.gz
MOLM13_LIB4_R1_001.fastq.gz
UBIK48_CTR_HIGH1_R1_001.fastq.gz
UBIK48_CTR_HIGH2_R1_001.fastq.gz
UBIK48_CTR_LOW1_R1_001.fastq.gz
UBIK48_CTR_LOW2_R1_001.fastq.gz
UBIK48_HS_HIGH1_R1_001.fastq.gz
UBIK48_HS_HIGH2_R1_001.fastq.gz
UBIK48_HS_LOW1_R1_001.fastq.gz
UBIK48_HS_LOW2_R1_001.fastq.gz

**SITA - MOLM13**
No replicates, except for MOLM13 unsorted
Two sequencing runs:

Run 1 - Low/no reads in 3/4 SITA sorted samples
Run 2 (GW5) - SITA_CTR_LOW2 has very high Gini coefficient and all 4 samples look undersampled.


**HSP90 - MOLM13**	

Sequenced twice:
First time (GW3) there were low/no reads from ¾ sorted samples, but correlations surprisingly good for all but one sample
Second time (GW5; using unsorted MOLM13 samples from GW3):
The Gini coefficient is very high for HS_LOW2 and correlations are poor
Samples are still under sequenced


**SITA - HAP1**


**HSP90 - HAP1**
