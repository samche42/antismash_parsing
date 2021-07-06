library(ggplot2)
library(gggenes)

gene_data <- read.table('genes.tab', header=TRUE, sep='\t')
domain_data <- read.table('domains.tab', header=TRUE, sep='\t')

#Draw genes with domains
domain_plot <- ggplot(gene_data, aes(xmin = start, xmax = end, y = molecule, forward=direction)) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  geom_gene_arrow(fill = "white") +
  geom_subgene_arrow(data = domain_data,
  aes(xmin = start, xmax = end, y = molecule, fill=subgene,
  xsubmin = sub_start, xsubmax = sub_end), color="black", alpha=.7) +
  theme_genes()
