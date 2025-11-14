# Import the library
library(tidyverse)

# ead the summary file for the counts
counts <- read.table("filename.txt", 
                      header = TRUE, 
                      sep = "\t", 
                      row.names = 1, 
                      check.names = FALSE)


View(counts)

# the below 3 lines are for data cleaning. this is personlised for each file. Omit these line if there is no need for data cleaning
counts <- counts %>% rownames_to_column(var = "Category")
colnames(counts) <- sub("_.*", "", colnames(counts))
colnames(counts)[2] <- "SRR13197315"

# create a new data frame(here a tibble) to calculate and store the percentages
counts_long <- counts %>%
     pivot_longer(-Category, names_to = "Sample", values_to = "Counts") %>%
     group_by(Sample) %>%
     mutate(Percent = Counts / sum(Counts) * 100)

View(counts_long)

# create a colour palette
n_categories <- length(unique(counts_long$Category))
my_colors <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(n_categories)

# Plot the stackd barplot 
ggplot(counts_long, aes(x = Sample, y = Percent, fill = Category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors) +
  labs(y = "Percentage of Reads", x = "Sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

