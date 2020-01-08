# Meth_Comp usage

**Compare methylation values for each CpG positions (`pycoMeth CpG_Aggregate`) or intervals (`pycoMeth Interval_Aggregate`) between n samples and perform a statistical test to evaluate if the positions are significantly different. For 2 samples a Mann_Withney test is performed otherwise multiples samples are compared with a Kruskal Wallis test. pValues are adjusted for multiple tests using the Benjamini & Hochberg procedure for controlling the false discovery rate.**

## Example usage

* [Python API usage](https://a-slide.github.io/pycoMeth/Meth_Comp/API_usage/)
* [Shell CLI usage](https://a-slide.github.io/pycoMeth/Meth_Comp/CLI_usage/)

## Output format

Meth_Comp can generates 2 files, a standard BED file and a tabulated file containing extra information

### Tabulated TSV file

This tabulated file contains the following fields:

* chrom / start / end / strand: Genomic coordinates of the motif or group of motifs in case split_group was not selected.
* n_samples: Number of valid samples compared for position
* pvalue / statistic: pvalue /statitic for positions to have significantly different median llr (Kruskal Wallis or Mann_Withney test)
* adj_pvalue: FDR adjusted pValue
* llr_list: List of median llr values for each samples compared.

### BED file

Standard genomic BED6 (https://genome.ucsc.edu/FAQ/FAQformat.html#format1). The score correspond to the methylation frequency multiplied by 1000. The file is sorted by coordinates and can be rendered with a genome browser such as IGV

The sites are color-coded as follow:

* Red: Median log likelyhood ratio >= 2 (more methylated)
* Blue: Median log likelyhood ratio <= -2 (more unmethylated)
* Grey: Median log likelyhood ration between -2 and 2 (ambiguous methylation status)

Here is an example of multiple methylation bed files rendered with IGV

![Example Bed Files](../pictures/Meth_Comp.png)
