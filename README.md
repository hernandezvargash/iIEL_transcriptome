# iIEL_transcriptome: Transcriptome analysis of mouse gamma-delta T cells of different origin


## R code corresponding to the article:

**Common and Exclusive Features of Intestinal Intraepithelial gamma-delta T Cells and Other gamma delta T Cell Subsets**, by Apostol K. Apostolov, Miriame Hamani, Hector Hernandez-Vargas, Ramdane Igalouzene, Alexandre Guyennon, Olivier Fesneau, Julien C. Marie, and Saïdi Mhoma Soudja.
*Tumor Escape Resistance Immunity Department, Cancer Research Center of Lyon, UMR INSERM 1052, CNRS 5286, Universite Claude Bernard Lyon I, Centre Leon Berard, Lyon, France*


---------------------------------

## Summary

Whole genome expression was studied in intestinal intraepithelial (IEL) TCR gamma delta cells, and 3 different subsets of TCR gamma delta cells (Type1,  Type naive-like and Type 17) that are present in lymph nodes. RNA was extracted in triplicates of each condition and hybridized into expression microarrays (genechip 2.0 ST from Affymetrix). The resulting CEL files were imported using the "oligo" package [1], and all downstream analyses were performed using this and other R/Bioconductor packages, R version 4.1.2 (2021-11-01) [https://cran.r-project.org/ ; http://www.bioconductor.org/]. The Robust Multichip Average (RMA) algorithm from the "oligo" package was used for normalization, followed by inspecion using principal component analysis (PCA) with "ggplot2" functions [2]. Differential expression was performed using "limma" functions to fit a linear model and contrasts for all pairwise comparisons [3]. Significant probes (FDR-adjusted p value < 0.05) were annotated with the "annotation mogene20sttranscriptcluster.db" package [4], and visualized with "EnhancedVolcano" for volcano plots [5], "eulerr" for genelist overlaps [6], and "NMF" for annotated supervised and unsupervised heatmaps [7]. Geneset enrichment analyses (GSEA) were performed using function from the packages "clusterprofiler" [8,9], "enrichplot" [10], "DOSE" [11], "ReactomePA" [12], "msigdbr" [13], and "fgsea" [14].

Data has been uploaded into the GEO repository, with Accession number **[GSE198703](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE198703)**.


---------------------------------

## References

1. Carvalho B. S., and Irizarry, R. A. 2010. A Framework for Oligonucleotide Microarray Preprocessing. Bioinformatics.

2. H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.

3. Ritchie, M.E., Phipson, B., Wu, D., Hu, Y., Law, C.W., Shi, W., and Smyth, G.K. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Research 43(7), e47.

4. James W. MacDonald (2021). mogene20sttranscriptcluster.db: Affymetrix mogene20 annotation data (chip mogene20sttranscriptcluster). R package version 8.8.0.

5. Kevin Blighe, Sharmila Rana and Myles Lewis (2021). EnhancedVolcano: Publication-ready volcano plots with enhanced colouring and labeling. R package version 1.10.0. https://github.com/kevinblighe/EnhancedVolcano

6. Larsson J (2020). _eulerr: Area-Proportional Euler and Venn Diagrams with Ellipses_. R package version 6.1.0, <URL: https://cran.r-project.org/package=eulerr>.

7. Renaud Gaujoux, Cathal Seoighe (2010). A flexible R package for nonnegative matrix factorization. BMC Bioinformatics 2010, 11:367. [http://www.biomedcentral.com/1471-2105/11/367]

8. T Wu, E Hu, S Xu, M Chen, P Guo, Z Dai, T Feng, L Zhou, W Tang, L Zhan, X Fu, S Liu, X Bo, and G Yu. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. The Innovation. 2021, 2(3):100141

9. Guangchuang Yu, Li-Gen Wang, Yanyan Han and Qing-Yu He. clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS: A Journal of Integrative Biology 2012, 16(5):284-287

10. Guangchuang Yu (2021). enrichplot: Visualization of Functional Enrichment Result. R package version 1.12.2. https://yulab-smu.top/biomedical-knowledge-mining-book/

11. Guangchuang Yu, Li-Gen Wang, Guang-Rong Yan, Qing-Yu He. DOSE: an R/Bioconductor package for Disease Ontology Semantic and Enrichment analysis. Bioinformatics 2015 31(4):608-609

12. Guangchuang Yu, Qing-Yu He. ReactomePA: an R/Bioconductor package for reactome pathway analysis and visualization. Molecular BioSystems 2016, 12(2):477-479

13. Igor Dolgalev (2021). msigdbr: MSigDB Gene Sets for Multiple Organisms in a Tidy Data Format. R package version 7.4.1. https://CRAN.R-project.org/package=msigdbr

14. G. Korotkevich, V. Sukhov, A. Sergushichev. Fast gene set enrichment analysis. bioRxiv (2019), doi:10.1101/060012

---------------------------------

