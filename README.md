# Selection Signatures in CIMMYT’s International Elite Spring and Semi-Arid Wheat Yield Trials

There is a chronological order for this work, which must be followed as new R objects are saved and used for downstream pipeline analysis. Please follow this order:

1. [Hapmap](hapmap)  Transforms textual hapmap fomrat into numeric
2. [Trial Data Cleaning](trial_data_cleaning) Data preparation from field trials
3. [PCA](PCA) Principal Component Analysis
4. [FST](FST) Fixation Index Metrics over loci and populations
5. [allele_frequencies_populations](allele_frequencies_populations) Estimation of population allele frequencies. Please follow this order script:
    * [allele_frequencies_1](allele_frequencies_populations/allele_frequencies_1.R)
    * [polishing_data_2](allele_frequencies_populations/polishing_data_2.R)
    * [kernel_3](allele_frequencies_populations/kernel_3.R)
    * [plot_signatures_selection_4](allele_frequencies_populations/plot_signatures_selection_4.R)
6. [adegenet_dapc](adegenet_dapc) Population Structure Analysis
