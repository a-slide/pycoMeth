# Meth_Comp usage

**Generate a fully responsive interactive HTML report for significant differentially methylated intervals found with `pycoMeth Meth_Comp`.**
**At the moment only data binned by intervals with `Interval_Aggregate` are supported.**

## Example usage

* [Python API usage](https://a-slide.github.io/pycoMeth/Comp_Report/API_usage/)
* [Shell CLI usage](https://a-slide.github.io/pycoMeth/Comp_Report/CLI_usage/)

## Input files

### pycoMeth CpG_Aggregate or Interval_Aggregate output file

A **TSV** output files generated with `pycoMeth Meth_comp` with data binned by intervals using `pycoMeth Interval_Aggregate`.

### Reference GFF3 file

An **Ensembl GFF3** file containing genomic annotations to extract transcripts with TSS close to the top most significant CpG Intervals.

* Main Ensembl: https://www.ensembl.org/index.html
* Ensembl Bacteria: http://bacteria.ensembl.org/index.html
* Ensembl Fungi: http://fungi.ensembl.org/index.html
* Ensembl Plants: http://plants.ensembl.org/index.html
* Ensembl Protists: http://protists.ensembl.org/index.html
* Ensembl Metazoa: http://metazoa.ensembl.org/index.html

## Output files

### Summary HTML report (pycoMeth_summary.html)

Entry point to the report containing information of methylation status at genomic level.

[Example Summary report](https://a-slide.github.io/pycoMeth/Comp_Report/medaka_html/pycoMeth_summary.html)

The report contains the following items:

* Overall summary
Simple table containing counts of intervals and CpGs significantly differentially methylated.

* Methylation category counts by sample
Interactive plotly stacked bar plot showing the number of methylated, unmethylated and ambiguous intervals for each samples.

* Methylation log-likelihood ratio by CpG interval
Interactive plotly heatmap of the median llr for all significant intervals.
Sites are ordered by genomic coordinates and samples are clustered by hierarchical clustering.

* Distribution of CpG methylation log-likelihood ratio by sample
Interactive plotly ridgeplot of median llr distribution for all sample at interval resolution.
Samples are ordered by descending overall median llr.

* Top differentially methylated intervals
Table containing details of the n top differentially methylated candidates ranked by adjusted pvalue.
The table contains links to individual intervals HTML reports.

### Top intervals HTML reports

Individual reports for the top differentially methylated interval candidates.
The right side navigation bar allows to explore all the other intervals or go back to the summary report.

[Example Interval report](https://a-slide.github.io/pycoMeth/Comp_Report/medaka_html/pycoMeth_interval_0001_chr15-13014693-13015794.html)

The report contains the following items:

* Interval details
Simple table containing details about the interval

* Methylation log-likelihood ratio by CpG position
Interactive plotly heatmap of the median llr for all CpG positions in the interval.
Sites are ordered by genomic coordinates and samples are clustered by hierarchical clustering.

* Distribution of CpG methylation log-likelihood ratio by sample
Interactive plotly ridgeplot of median llr distribution for all sample at CpG resolution.
Samples are ordered by descending overall median llr.

* Closest transcripts TSS
Table containing the list of all transcripts having their TSS within 100kb of the interval edges.
Transcripts are ranked by TSS distance to the interval
