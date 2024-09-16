Unveiling RNA Structure-mediated Regulations of RNA Stability in Wheat

Introduction
This repository contains scripts and workflows for analyzing RNA structure-mediated regulations of RNA stability in wheat. This documentation is divided into five main categories, each represented by a specific set of scripts. In addition, we have included detailed steps for data preparation required for the analyses.

1.	Data preparation
The RNA-seq mapping and SNP calling workflow for the RNA-seq libraries of B028, B086, and Kronos accessions is consistent with the methods presented in Yang's study1. If there is any overlapping data, previously collected data will be utilized directly. Hisat2 v2.1.02 was used for mapping the trimmed reads to the durum wheat genome assembly (Svevo RefSeq 1.0)3 corrected by accession-specific SNVs information from previous study4.

2.  The core softwares and workflows:
The workflows of this part were developed based on two published software tools for RNA structural motifs prediction and related analyses on modeling of mRNA Decay. The workflows applied in this study for RNA structural motifs analysis and modeling of decay rates are comprehensively described in the materials and methods section of our manuscript.
2.1	Systematic identification of RNA structural motifs associated with stability was performed using the pyTEISER framework5. For more information, visit pyTEISER on GitHub website (https://github.com/goodarzilab/pyteiser).
2.2	Normalization and modeling of mRNA decay profiles were conducted using the tools and methods outlined in the RNAdecay package6. For detailed guidance on these processes, refer to the RNAdecay workflow documentation (https://bioconductor.org/packages-/release/bioc/vignettes/RNAdecay/inst/doc/RNAdecay_workflow.html). This R script modified according to the genetic characteristics of wheat have also been uploaded.
Script1: RNA_decay_rate_calculate_for_wheat.R 

3.  Parameters Calculation Python Scripts:
These Python scripts are developed for computing various parameters critical to the analysis, such as codon adaptation index (CAI), codon stabilization coefficient (CSC), frequency of codons in each transcript, GC content, AU content, sequence length, and intron numbers. These parameters are the potential factors to affect RNA stability. 
3.1 Reads count through bam file
Script2: atlas_count_rnaseq
3.2 GC content
Script3: getGCcontent.py
3.3 AU content
Script4: getATcontent
3.4 sequence length
Script5: count_length_for_fasta.py
3.5 intron numbers count through .gff3 file
Script6: get_intron_number_from_gff3.py

4. R Script for data organization and analysis:
These R scripts are primarily developed for data cleaning, organization, and classification, significance analysis, and correlation analysis, as well as for conducting GO enrichment analysis.
4.1 Data cleaning and organization
Script7: abundance_filted_0min_greater_1_all.R
Script8: merge_dada_file_by_ID_find_gene_pairs.R
4.2 Min_max_normalization
Script9: nomalized_by_Min_max_normalization.R
4.3 Significance analysis
Script10: AvsB_difference_t-test.R
4.4 Correlation analysis
Script11: Bivariate_Correlations_analysis_pearson_correlation.R

References

1.	Yang X., et al. Wheat in vivo RNA structure landscape reveals a prevalent role of RNA structure in modulating translational subgenome expression asymmetry. Genome Biol. 22, 326 (2021).
2.	Kim D., Langmead B. & Salzberg S. L. HISAT: a fast spliced aligner with low memory requirements. Nat. Methods 12, 357-360 (2015).
3.	Maccaferri M., et al. Durum wheat genome highlights past domestication signatures and future improvement targets. Nat. Genet. 51, 885-895 (2019).
4.	Zhou Y., et al. Triticum population sequencing provides insights into wheat adaptation. Nat. Genet. 52, 1412-1422 (2020).
5.	Fish L., et al. A prometastatic splicing program regulated by SNRPA1 interactions with structured RNA elements. Science 372,  (2021).
6.	Sorenson R. S., Deshotel M. J., Johnson K., Adler F. R. & Sieburth L. E. Arabidopsis mRNA decay landscape arises from specialized RNA decay substrates, decapping-mediated feedback, and redundancy. Proc. Natl. Acad. Sci. U.S.A. 115, E1485-E1494 (2018).

