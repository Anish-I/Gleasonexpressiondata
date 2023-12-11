####################################################################
# Name: Anish Ivaturi
# Final project: Identification of differential expressed genes
# Comparison based on Gleason Score
####################################################################

This R script is part of a final project conducted by Anish Ivaturi for the identification of differentially expressed genes in prostate cancer based on Gleason scores. Below is an overview of the code and its key components:

### Libraries ###
The script begins by loading necessary R libraries, including dplyr, ggplot2, UCSCXenaTools, edgeR, and limma, which are essential for data processing and analysis.

### Data Retrieval and Preprocessing ###
1. The script retrieves data from the UCSC Xena browser, specifically the "GDC TCGA Prostate Cancer (PRAD)" dataset.
2. It extracts both clinical and RNA-seq data.
3. The data is preprocessed, including normalization, probe filtering, and aligning sample IDs to ensure data consistency.

### Exploratory Data Analysis ###
4. The script generates a boxplot to visualize the normalized data for the first 10 samples with different colors.

### Data Summary ###
5. The script calculates and prints the total number of samples and probes before and after data processing.

### Gleason Score Extraction ###
6. The Gleason scores are extracted from clinical data, removing non-numeric characters and converting them to numeric values.

### Differential Expression Analysis ###
7. A design matrix is created for comparing gene expression based on Gleason scores.
8. Differential expression analysis is performed using the limma package, applying a False Discovery Rate (FDR) threshold of 10%.
9. The script identifies and prints the top 30 differentially expressed probes, or fewer if the total number of probes is less than 30.

### Data Visualization ###
10. A boxplot is created to visualize the expression of a selected probe based on Gleason scores.
11. Gene names are mapped to probes, and a dataframe of results is created, including gene names, probes, log fold change, and adjusted p-values.
12. The top 5 differentially expressed probes are printed.

### Heatmap Visualization ###
13. The script generates a heatmap of the top 30 genes based on mean logCPM values, with samples color-coded by Gleason score groups (low and high).
14. The heatmap provides insights into the expression patterns of genes associated with Gleason scores.

### External Data Integration ###
15. External data from a file named "chart_F625A94F7CA01702314584249.txt" is read into the script. The specific purpose of this data integration is not explained in the code and may require additional context.

### Conclusion ###
The README concludes by summarizing the project's objectives and provides a link to the original data source on the UCSC Xena browser.

For any questions or clarifications about the code, please contact Anish Ivaturi.
