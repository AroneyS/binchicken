# Why use Bin Chicken?

Coassembly — pooling reads from multiple samples before assembly — recovers more and better-quality genomes than assembling each sample individually.
Bin Chicken automates the selection of which samples to coassemble, making this benefit practical at scale.

## More genomes recovered

In a benchmarking study across 800 coassemblies (Aroney et al., 2025), targeted coassembly with Bin Chicken recovered **39% more species-level genome bins on average** compared to single-sample assembly of the same data.

## Efficiency at scale

Relative to SPIRE — which recovered 83,000 novel species from 500 Tbp across 99,146 samples — Bin Chicken recovered 24,000 novel species from only 17.8 Tbp across 800 coassemblies: **more than 20× the novel species yield per assembled metagenome**.

## Better genome quality

Where both approaches recovered the same species, coassembled genomes were higher quality (`completeness − 5× contamination`): **2.8% higher quality** on average (P < 0.001, n = 2,170).
This is consistent with coassembly providing greater depth over shared low-abundance organisms, producing more contiguous assemblies.

## No increase in chimerism

A common concern with coassembly is that mixing reads from multiple samples will produce chimeric bins.
In practice, Bin Chicken's use of exact matching of marker gene sequences (ensuring samples share lineages with at least species-level resolution) and exclusion of reads that map to known reference genomes substantially reduces this risk.
Chimera rates in recovered genomes (9.8%) were comparable to those found in GTDB reference databases (8.0%), showing **no discernible increase in chimerism** attributable to coassembly.

These results are shown in Extended Data Fig. 3 of the Bin Chicken paper:

![Extended Data Fig. 3 — coassembly outperforms single-sample assembly in genome recovery, quality, and chimerism](/figureS3.png)

**Extended Data Fig. 3 Benchmarking genomes recovered from coassembly against single-sample assembly with co-binning.**
(a) Comparison of genomes recovered from 30 randomly selected coassemblies and 70 single-sample assemblies from each of the samples that were coassembled. Single only vs coassembly only p = 7.86e-7, n = 30 sample-sets. (b,c,d,e,f) Comparison of genomes recovered from coassembly and single-sample from within 95% ANI species clusters. Columns represent mean ± standard deviation of (b) genome quality (CheckM2 completeness – 5 × contamination) with p = 8e-34, n = 2170 genome-pairs, (c) completeness with p = 2e-19, n = 2170 genome-pairs, (d) contamination with p = 1e-3, n = 2170 genome-pairs, (e) length-weighted contig covered fraction from one sample with p = 0.818, n = 2170 genome-pairs, and (f) proportion of genomes with GUNC detected chimerism per sample-set with p = 0.002, n = 30 sample-sets. Points represent genomes (b, c, d, e) or sample-sets (f). Data were modelled using two-sided ANOVA. NS: not significant, ** p < 0.01, *** p < 0.001.

> Aroney, S.T.N., Newell, R.J.P., Tyson, G.W. and Woodcroft B.J. *Bin Chicken: targeted metagenomic coassembly for the efficient recovery of novel genomes.* Nat Methods (2025). <https://doi.org/10.1038/s41592-025-02901-1>
