# Meth_Comp usage

**Compare methylation values for each CpG positions (`pycoMeth CpG_Aggregate`) or intervals (`pycoMeth Interval_Aggregate`) between n samples and perform a statistical test to evaluate if the positions are significantly different. For 2 samples a Mann_Withney test is performed otherwise a Kruskal Wallis test is performed for multiple samples. pValues are adjusted for multiple tests using the Benjamini & Hochberg procedure for controlling the false discovery rate.**

## Example usage

* [Python API usage](https://a-slide.github.io/pycoMeth/Meth_Comp/API_usage/)
* [Shell CLI usage](https://a-slide.github.io/pycoMeth/Meth_Comp/CLI_usage/)

## Input files

### pycoMeth CpG_Aggregate or Interval_Aggregate output file

A list of `pycoMeth CpG_Aggregate` or `pycoMeth Interval_Aggregate` **tsv** output files corresponding to different samples.

### Reference FASTA file

FASTA reference file used for read alignment and Nanopolish. This file is required and used to sort the CpG sites by coordinates.

## Output files

### Tabulated TSV file

This tabulated file contains the following fields:

* chromosome / start / end / strand: Genomic coordinates of the CpG/interval containing CpGs
* n_samples: Number of valid samples in the CpG/interval
* pvalue: pvalue of the interval obtained by Kruskal Wallis or Mann_Withney test
* adj_pvalue: FDR adjusted pvalue using the Benjamini & Hochberg procedure  
* neg_med / pos_med / ambiguous_med: Number of samples with a median below the negative llr threshold / above the positive llr threshold or with and ambiguous median between the 2 thresholds
* [optional] unique_cpg_pos: Number of CpG position for which at least samples has valid data. ONLY if Meth_Comp was ran from files generated with Interval_Aggregate.
* labels: labels of the samples tested, matching the order of values in med_llr_list and raw_llr_list
* med_llr_list: List of median llr values for each samples compared.
* raw_llr_list: List of lists of raw llr values for each samples compared
* [optional] raw_pos_list: List of lists of genomic coordinates of the center of all CpG positions/clusters found for each samples compared. The order matches values in llr_list. ONLY if Meth_Comp was ran from files generated with Interval_Aggregate.

### BED file

Standard genomic [BED9 format](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) including an RGB color field. The score correspond to the -log10(Adjusted Pvalue) capped to 1000. The file is sorted by coordinates and can be rendered with a genome browser such as IGV

The sites are color-coded as follow:

* Significant differential methylation Adjusted pValue:  Colorscale from orange (pValue=0.01) to deep purple (pValue<=0.000001)
* Non-significant: Grey

Here is an example of multiple methylation bed files rendered with IGV

![Example Bed Files](../pictures/Meth_Comp.png)
