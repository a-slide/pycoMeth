# Interval_Aggregate usage

**Bin the output of `pycoMeth CpG_Aggregate` in genomic intervals, using either an annotation file containing interval or a sliding window.**

## Example usage

* [Python API usage](https://a-slide.github.io/pycoMeth/Interval_Aggregate/API_usage/)
* [Shell CLI usage](https://a-slide.github.io/pycoMeth/Interval_Aggregate/CLI_usage/)

## Input files

### pycoMeth CpG_Aggregate output file

pycoMeth CpG_Aggregate **tsv** output file (not the bed output).

### Reference FASTA file

FASTA reference file used for read alignment and Nanopolish. This file is required and used to sort the CpG sites by coordinates 

### BED file containing intervals

Optional **sorted** and BED file containing **non-overlapping** intervals to bin CpG data into. If this file is not provided, then the program use a sliding customizable window to bim data along the entire genome.

## Output format

CpG_Aggregate can generates 2 files, a standard BED file and a tabulated file containing extra information

### Tabulated TSV file

This tabulated file contains the following fields:

* chromosome / start / end / strand: Genomic coordinates of the CpG or group of CpGs if in less than 5 bases from each other
* num_CpG_clusters: Number of CpG clusters in the interval
* median_llr: Median of log likelihood ratios for each CpG cluster aggregated
* llr_list: List of median llr values for all CpG cluster aggregated

### BED file

Standard genomic [BED6](https://genome.ucsc.edu/FAQ/FAQformat.html#format1). The score correspond to the median log likelyhood ratio.
The file is already sorted by coordinates and can be rendered with a genome browser such as IGV

The sites are color-coded as follow:

* Red: Median log likelyhood ratio >= 3 (more methylated)
* Orange: Median log likelyhood ratio >= 2 (more methylated)
* Blue: Median log likelyhood ratio <= -2 (more unmethylated)
* Light blue: Median log likelyhood ratio <= -3 (more methylated)
* Grey: Median log likelyhood ration between -2 and 2 (ambiguous methylation status)

Here is an example of multiple methylation bed files rendered with IGV

![Example Bed Files](../pictures/Interval_Aggregate.png)
