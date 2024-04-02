library(dplyr)
library(tidyverse)
library("org.Hs.eg.db")
library(clusterProfiler)
library(gprofiler2)


enrichment_anaysis <- function(protein_list, plot_file = FALSE) {
  gost_object <- gost(protein_list, organism = "hsapiens", significant = TRUE, correction_method='bonferroni')
  if (is.null(gost_object)) {
    print('No gost results')
    return(NA)
  } else {
    enrich_df <- gost_object$result
    top_annot <- enrich_df %>% 
      filter(source %in% c('GO:BP','CORUM', 'KEGG', 'HPA', 'REAC', 'WP')) %>%
      group_by(source) %>% 
      slice_max(order_by = -p_value, n = 3)

    if (plot_file != F) {
      p1 <- gostplot(gost_object, interactive = FALSE)
      p1_pub <- publish_gostplot(p1, highlight_terms = top_annot$term_id, width=5, height=8)
      ggsave(plot_file, p1_pub)
    }
    return(enrich_df)
  }
}

# Prompt for user input
# e.g., For method_prefix, combinatorial, NMF, or 'hierarchical
# method_prefix = readline(prompt = 'Enter the grouping method used: ')
method_prefix = 'combinatorial'

# Define paths to read and write files
project_data_dir = './project_data/'
global_data_dir = './global_data/'

group_path = paste(project_data_dir, 'ARID1A_', '_', method_prefix, '_grouped_protein.csv', sep = '')
enrichment_path = paste(project_data_dir, 'ARID1A_', method_prefix, '_enrichment_terms.csv', sep = '')

groups <- read.csv(group_path)
group_enrichment <- data.frame()

for (i in 1:nrow(groups)){
  group <- groups[i,]
  print(group)
  protein_list <- unlist(strsplit(group$proteins, split = ','))
  print(protein_list)
  group_enrich <- enrichment_anaysis(protein_list)
  if (is.na(group_enrich)) {
    print(paste('Skipped the group', group['group']))
  }
  else {
    group_enrich['group'] <- group['group']
    group_enrichment <- rbind(group_enrichment, group_enrich)
  }
}

group_enrichment <- apply(group_enrichment, 2, as.character)
group_enrichment

write.csv(group_enrichment, row.names = FALSE, enrichment_path)

group_enrichment[duplicated(group_enrichment)]
print(enrichment_path)
