####################################################################
# Name: Anish Ivaturi
# Final project: Identification of differential expressed genes
# Comparison based on Gleason Score
####################################################################

# Load necessary libraries
library(dplyr)
library(ggplot2)
library(UCSCXenaTools)  
library(edgeR)        
library(limma)        

# Retrieve and preprocess the data
data(XenaData)
prostatecancer <- XenaData %>% filter(XenaCohorts == 'GDC TCGA Prostate Cancer (PRAD)')

# Clinical data
cli_query <- prostatecancer %>% filter(Label == "Phenotype") %>% XenaGenerate() %>% XenaQuery() %>% XenaDownload()
prostatecancer_pheno <- XenaPrepare(cli_query)

# RNA-seq data
cli_query <- prostatecancer %>% filter(Label == 'HTSeq - Counts') %>% XenaGenerate() %>% XenaQuery() %>% XenaDownload(download_probeMap = TRUE)
prostatecancer_counts <- XenaPrepare(cli_query)

# Data pre-processing
X <- data.frame(prostatecancer_counts$TCGA.PRAD.htseq_counts.tsv.gz)
rownames(X) <- X$Ensembl_ID
X <- X[,-1]  # Remove probe name column

probeMap <- prostatecancer_counts$gencode.v22.annotation.gene.probeMap
Y <- prostatecancer_pheno
colnames(X) <- gsub('\\.', '-', colnames(X))  # Format sample IDs
g <- grep('01A$', colnames(X))  # Keep '01A' tumor samples
X <- X[,g]
common_samples <- intersect(colnames(X), Y$submitter_id.samples)
mx <- match(common_samples, colnames(X))
my <- match(common_samples, Y$submitter_id.samples)
X <- X[,mx]
Y <- Y[my,]
stopifnot(all(colnames(X) == Y$submitter_id.samples))  # Check sample alignment

# Convert log2(count + 1) to count data and normalize
X <- round(2**X - 1)
dge <- DGEList(counts=X)
keep <- filterByExpr(dge, min.prop = .10)
dge <- dge[keep,,keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge, method = "TMM")
logCPM <- cpm(dge, log = TRUE, prior.count = 3)

# Generate a boxplot for the first 10 samples with colors
boxplot(logCPM[, 1:10], 
        main = "Boxplot of Normalized Data for First 10 Samples",
        ylab = "log CPM",
        col = rainbow(10))  # Use rainbow colors for the boxes


# Total number of samples with expression data before processing
total_samples_before <- ncol(dge)

# Total number of probes before processing
total_probes_before <- nrow(dge$counts)

# Total number of samples after processing
total_samples_after <- ncol(logCPM)

# Total number of probes after processing
total_probes_after <- nrow(logCPM)

# Print the counts
cat("Total number of samples with expression data before processing:", total_samples_before, "\n")
cat("Total number of probes before processing:", total_probes_before, "\n")
cat("Total number of samples with expression data after processing:", total_samples_after, "\n")
cat("Total number of probes after processing:", total_probes_after, "\n")

# Extract Gleason score column
gleason_scores <- Y$gleason_score

# Remove non-numeric characters
gleason_scores <- gsub("[^0-9]", "", gleason_scores)

# Convert to numeric if they are not already
gleason_scores <- as.numeric(gleason_scores)

gleason_scores <- na.omit(gleason_scores)
gleason_score_counts <- table(gleason_scores)

# Print the counts
print(gleason_score_counts)


library(limma)

# Create a design matrix for the comparison
design <- model.matrix(~gleason_scores)

# Perform differential expression analysis using limma
fit <- lmFit(logCPM, design)
fit <- eBayes(fit)

fdr_threshold <- 0.1
de_probes <- decideTests(fit, method = "separate", adjust.method = "fdr", p.value = fdr_threshold)
num_de_probes <- sum(de_probes != 0)

# If the number of probes is less than 30, find the top 30 probes
if (num_de_probes < 30) {
  # Get the top 30 probes based on adjusted p-values
  top_de_probes <- topTable(fit, number = 30, sort.by = "p", adjust.method = "fdr")
  
  # Print the top 30 probes
  print(top_de_probes)
} else {
  # Print the number of probes identified
  print(paste("Number of probes identified at FDR =", fdr_threshold, ":", num_de_probes))
}
class(dge)


probe_name <- rownames(logCPM)[1] 
plot_data <- data.frame(Group = gleason_scores, Expression = logCPM[probe_name, ])

# Create a boxplot 
library(ggplot2)

ggplot(plot_data, aes(x = factor(Group), y = Expression)) +
  geom_boxplot() +
  labs(
    title = paste("Probe:", probe_name, "Expression Boxplot"),
    x = "Gleason Score Group",
    y = "Expression (logCPM)"
  )

# Map gene names to probes
gene_names <- probeMap$gene[match(rownames(fit$coefficients), probeMap$id)]

# Check if gene_names has NA values and handle them
if (any(is.na(gene_names))) {
  gene_names[is.na(gene_names)] <- "Unknown"
}

# Create the Data Frame
results_df <- data.frame(
  GeneName = gene_names,
  ProbeName = rownames(fit$coefficients),
  logFC = fit$coefficients[, "gleason_scores"],  # log fold change
  AdjPValue = fit$p.value[, "gleason_scores"]    #for adjusted p-values
)

#Top 5 Probes
top_probes_df <- results_df[order(results_df$AdjPValue), ][1:5, ]
print(top_probes_df)


library(heatmaply)

# Calculate the mean logCPM
mean_logCPM <- rowMeans(logCPM)
sorted_genes <- names(sort(mean_logCPM, decreasing = TRUE))

# Define the number of top genes you want to select
top_gene_count <- 30

# Select the top genes
top_genes <- sorted_genes[1:top_gene_count]
subset_expression <- logCPM[top_genes, ]
threshold_value <- 7  #threshold for separating low and high Gleason scores
group_vector <- ifelse(gleason_scores < threshold_value, "Low Gleason", "High Gleason")
custom_color_palette <- colorRampPalette(c("yellow", "white", "blue"))(100)
heatmaply(
  subset_expression,
  colors = custom_color_palette,
  labels_col = group_vector,  # Color code samples by groups
  k_col = 2,  # Number of colors to use
  fontsize_row = 8,  # Adjust row label font size
  main = "Heatmap of Top 30 Genes Based on Mean logCPM"
)

gene_list <- rownames(logCPM)
gene_list

david_results <- read.table("C:/Users/ivatu/Downloads/Prostate/chart_F625A94F7CA01702314584249.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
print(david_results)

