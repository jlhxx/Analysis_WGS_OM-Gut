library(phyloseq)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)

# Phyloseq object called 'ps1.com.fam.rel'

# Extract sample data
sample_data_df <- as.data.frame(as.matrix(sample_data(ps1.com.fam.rel)))

# Get unique sample names
sample_names <- unique(sample_data_df$sample2)

# Loop over each sample and create plots
for (sample in sample_names) {
  
  # Filter for current sample
  current_sample_df <- sample_data_df %>% filter(sample2 == sample)
  
  # Get sample names for the current sample
  current_sample_names <- rownames(current_sample_df)
  
  # Create vector for x-axis labels - remove prefix from sample names
  x_labels <- setNames(
    gsub("^W[0-9]+_", "", current_sample_names),
    current_sample_names
  )
  
  # Subset the phyloseq object for the current sample
  ps1_current <- prune_samples(current_sample_names, ps1.com.fam.rel)
  
  # Convert OTU table of the subsetted phyloseq object into a dataframe
  otu_table_current <- as.data.frame(otu_table(ps1_current))
  
  # Convert rownames column
  otu_table_current <- otu_table_current %>%
    rownames_to_column(var = "Taxa")
  
  # Reshape data - wide to long format
  otu_table_long <- otu_table_current %>%
    pivot_longer(
      cols = -Taxa,
      names_to = "Sample",
      values_to = "Abundance"
    )
  
  # Convert factors to determine order in the plot
  otu_table_long$Sample <- factor(otu_table_long$Sample, 
                                  levels = unique(otu_table_long$Sample))
  
  # Plot stacked barplots side-by-side
  p <- ggplot(otu_table_long, aes(x = Sample, y = Abundance, fill = Taxa)) +
    geom_bar(stat = "identity", color = "black") +
    theme_minimal() +
    labs(title = paste("Stacked Barplots of Abundance for", sample),
         x = "Sample Groups", y = "Abundance") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          axis.text.y = element_text(size = 12),
          legend.text = element_text(size = 12)) +
    scale_x_discrete(labels = x_labels)
  
  # Save plot as .png
  ggsave(filename = paste0("plot_", sample, ".png"), plot = p)
}
