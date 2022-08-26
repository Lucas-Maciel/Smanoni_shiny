# Schistosoma mansoni reference genes ShinyApp

http://verjolab.shinyapps.io/Reference-genes/

This site shows the expression of Schistosoma mansoni protein-coding genes (Smps) and long non-coding RNAs (SmLINC, SmLNCA, SmLNCS RNAs) that comprise the dataset which was used to determine stable reference genes (Smps) across six life-cycle stages of the parasite (Silveira et al., Assessment of reference genes at six different developmental stages of Schistosoma mansoni for quantitative RT-PCR, Scientific Reports 2021 https://doi.org/10.1038/s41598-021-96055-7 ).
As explained in the Silveira et al. 2021 paper, after reads counting with RSEM we performed a minimal pre-filtering to keep only genes that had at least 10 reads total when adding all stages. This resulted in 13,624 protein-coding genes (out of 14,548 genes in the v7.1 transcriptome (WBPS14), PMID: 27899279) and in 9,391 lncRNAs (out of 16,583 lncRNAs, PMID: 31572441) that were considered for further analyses. These genes are shown in the present site. Four different methods were used to normalize the expression across all samples from the six different life-cycle stages, namely DESeq2, Trimmed mean of M-values (TMM.CPM), Upper Quartile (UQ), and Transcripts per Million (TPM).
You can search for your gene of interest (Smp or lncRNA) by typing its gene ID name below.

Data necessary to generate the app locally is not present in this repo, just the code.

