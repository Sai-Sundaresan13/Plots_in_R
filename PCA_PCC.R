# PCA plot for wichers individual

counts <- read.table("counts3.txt", header = TRUE, row.names = 1)
meta <- read.csv("W.csv", header = TRUE, row.names = 1)
all(rownames(meta) == colnames(counts))

counts <- counts[, -c(1:5)]
colnames(counts) <- sub("_.*", "", colnames(counts))
all(rownames(meta) == colnames(counts))

View(meta)
meta$Time <- factor(meta$Time)   # ensure factor

dds <- DESeqDataSetFromMatrix(countData = counts,
 colData = meta,
 design = ~ Time)

vst_counts <- vst(dds, blind = TRUE)
vst_mat <- assay(vst_counts)
vst_mat <- vst_mat[apply(vst_mat, 1, var) > 0, ]
pca <- prcomp(t(vst_mat), scale. = TRUE)
pca_var <- pca$sdev^2
pca_var_percent <- pca_var / sum(pca_var) * 100

scree_df <- data.frame(
 PC = paste0("PC", 1:length(pca_var_percent)),
 Variance = pca_var_percent)

scree_df$PC <- factor(scree_df$PC, levels = paste0("PC", 1:length(pca_var_percent))) # TO ORDER THE PCA BY NUMBERS (ASCENDING ORDER)

ggplot(scree_df, aes(x = PC, y = Variance)) +
 geom_bar(stat = "identity", fill = "skyblue") +
 geom_text(aes(label = round(Variance, 1)), vjust = -0.3, size = 3) +
 theme_minimal() +
 labs(
 title = "Scree Plot of PCA",
 x = "Principal Components",
 y = "% Variance Explained"
 ) +
 theme(
 axis.text.x = element_text(angle = 45, hjust = 1)
 )

pca_df <- as.data.frame(pca$x[, 1:5])
pca_df$Sample <- rownames(pca_df)
pca_df$Condition <- meta$Time
pca_df$Sample <- rownames(pca_df)
pca_df <- cbind(pca_df, meta[rownames(pca_df), ])
pca_df

pca_df$Condition <- factor(pca_df$Time, levels = c("0 hpi", "8 hpi", "16 hpi", "24 hpi", "32 hpi", "40 hpi", "44 hpi", "48 hpi"))

pca_df$Replicate <- as.factor(pca_df$Replicate)

ggplot(pca_df, aes(PC1, PC2, color = Condition, shape = Replicate, label = Plotname)) +
     geom_point(size = 3) +
     xlab(paste0("PC1: ", round(pca_var_percent[1], 1), "% variance")) +
     ylab(paste0("PC2: ", round(pca_var_percent[2], 1), "% variance")) +
     theme_bw() +
     labs(title = "PCA of Samples") +
     theme(
         legend.title = element_text(size = 10),
         legend.text = element_text(size = 9)
     )


pcc_mat <- cor(vst_mat, method = "pearson")

head(rownames(pcc_mat))
head(rownames(meta))
head(meta$Plotname)

new_labels <- meta$Plotname[match(rownames(pcc_mat), rownames(meta))]
rownames(pcc_mat) <- new_labels
colnames(pcc_mat) <- new_labels

pheatmap(
pcc_mat,
main = "Sample-to-Sample Pearson Correlation",
clustering_distance_rows = "correlation",
clustering_distance_cols = "correlation"
)



# PCA with outline and change of legend name
ggplot(pca_df, aes(PC1, PC2, fill = Condition, shape = Replicate)) +
     geom_point(size = 3, colour = "black", stroke = 0.8) +   # black outline
     scale_fill_brewer(palette = "Set2", name = "Time Point") +                    # simple palette
     scale_shape_manual(values = c(21, 22, 23, 24)) +         # 4 easy fillable shapes
     labs(
         title = "PCA of Samples",
         x = paste0("PC1: ", round(pca_var_percent[1], 1), "% variance"),
         y = paste0("PC2: ", round(pca_var_percent[2], 1), "% variance")
     ) +
     theme_bw() +
     theme(
         legend.title = element_text(size = 10),
         legend.text  = element_text(size = 9)
     ) +
     guides(fill = guide_legend(override.aes = list(shape = 21)))

