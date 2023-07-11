# kapa_nhl

## Overview
This repository contains the template scripts described in the [KAPA Bioinformatics Analysis for Longitudinal Detection of Circulating Tumor DNA](https://sequencing.roche.com/content/dam/diagnostics_microsites/sequencing/master-blueprint/en/resources/pdfs/white-papers/kapa-bioinformatics-analysis-for-longitudinal-detection-of-circulating-tumor-dna-white-paper-mc--12095.pdf) analysis white paper. 

The analysis covers the following modules:
* `main.sh`: Main analysis (single-sample, FASTQ to BAM to variant calling). 
* `longitudinal.sh`: Longitudinal analysis (multi-sample, longitudinal evaluation of mutations). Required scripts are in the `R` folder. 

The application note [Longitudinal detection of non-Hodgkin lymphoma ctDNA](https://sequencing.roche.com/content/dam/diagnostics_microsites/sequencing/master-blueprint/en/resources/pdfs/white-papers/longitudinal-detection-of-circulating-tumor-dna-white-paper-mc--11981.pdf) provides an example of the analysis using the KAPA HyperCap Design Share non-Hodgkin lymphoma panel. 

For more information on the KAPA HyperCap Design Share NHL panel, please refer to [https://sequencing.roche.com/global/en/products/group/kapa-hypercap-ds-nhl-panel.html](https://sequencing.roche.com/global/en/products/group/kapa-hypercap-ds-nhl-panel.html)

---

## Usage Tips

These scripts are meant to serve as templates only. **For more usage details, please refer to the [analysis white paper](https://sequencing.roche.com/content/dam/diagnostics_microsites/sequencing/master-blueprint/en/resources/pdfs/white-papers/kapa-bioinformatics-analysis-for-longitudinal-detection-of-circulating-tumor-dna-white-paper-mc--12095.pdf).**

The analysis scripts describes a number of steps and mini-workflows that use such third-party tools, which can be combined together into a variety of data analysis workflows for research.

Ideally, the user should develop a workflow appropriate for their experimental data using benchmark or control samples. Parameters used in the analysis should be carefully considered. 

- Note that where the text “SAMPLE” appears throughout examples shown here, the user should replace it with a unique sample name. 

- Similarly, replace “DESIGN” with the name of the target enrichment design that matches the design files supplied. 

- "NumProcessors” should be replaced with the number of CPU cores available.

- The current directory is assumed to be the location of all input files, and will also be the location of output files and report files. 

_Please note that publicly available, open-source software tools may
change and that such change is not under the control of Roche.
Therefore, Roche does not warrant and cannot be held liable for the
results obtained when using the third-party tools described herein._

---

## Prerequisites
The following tools are used:  
- BWA (0.7.17)
- FastQC (0.11.9) 
- GATK4 (4.2.0.0); Java (>1.8.0_282) 
- SAMtools (1.13)
- seqtk (1.3-r106) 
- fastp (0.20.1) 
- fgbio (1.3.0)
- VarDict java (1.8.3)
- BEDTools (2.30.0) 
- SnpSift (4.3t) 
- ctDNAtools (0.4.0)

In addition, the following R packages are required for longitudinal analysis with ctDNAtools. 
- R (≥3.6.0)
- magrittr
- dplyr (≥0.8.3)
- tidyr (≥1.0.0)
- purrr (≥0.3.2)
- Rsamtools (≥2.0.0)
- assertthat (≥0.2.1)
- GenomicRanges
- IRanges
- GenomeInfoDb
- BiocGenerics
- GenomicAlignments
- rlang (≥0.4.0)
- Biostrings
- methods
- furrr (≥0.1.0)
- ellipsis (≥0.3.0)
- VariantAnnotation (≥1.30.1)
- GenomeInfoDbData
- optparse
- BSgenome.Hsapiens (version according to genome choice)

---

## License
Distributed under Apache License 2.0. See `LICENSE.txt` for more information.

