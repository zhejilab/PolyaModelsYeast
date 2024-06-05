
## Import packages ------------------------------------------------------------

packages <- c("tidyverse", "RColorBrewer", "pheatmap", "gridExtra", "ggrepel", "magrittr", "matrixStats","cowplot","viridis","dendextend")
lapply(packages, library, character.only = TRUE)

## Visualize heatmaps: ARABIDOPSIS --------------------------------------------

t_rich <- c('ACTT', 'AGTT', 'ATCT', 'ATGT', 'ATTC', 'ATTG', 'ATTT', 'CATT', 'CCTT', 'CGTT', 'CTAT', 'CTCT', 'CTGT', 'CTTA', 'CTTC',
            'CTTG', 'CTTT', 'GATT', 'GCTT', 'GGTT', 'GTAT', 'GTCT', 'GTGT', 'GTTA', 'GTTC', 'GTTG', 'GTTT', 'TACT', 'TAGT', 'TATC',
            'TATG', 'TATT', 'TCAT', 'TCCT', 'TCGT', 'TCTA', 'TCTC', 'TCTG', 'TCTT', 'TGAT', 'TGCT', 'TGGT', 'TGTC', 'TGTG', 'TGTT',
            'TTAC', 'TTAG', 'TTAT', 'TTCA', 'TTCC', 'TTCG', 'TTCT', 'TTGA', 'TTGC', 'TTGG', 'TTGT', 'TTTA', 'TTTC', 'TTTG', 'TTTT')

a_rich <- c('AAAA', 'AAAC', 'AAAG', 'AAAT', 'AACA', 'AACC', 'AACG', 'AACT', 'AAGA', 'AAGC', 'AAGG', 'AAGT', 'AATA', 'AATC', 'AATG',
            'AATT', 'ACAA', 'ACAC', 'ACAG', 'ACAT', 'ACCA', 'ACGA', 'ACTA', 'AGAA', 'AGAC', 'AGAG', 'AGAT', 'AGCA', 'AGGA', 'AGTA',
            'ATAA', 'ATAC', 'ATAG', 'ATAT', 'ATCA', 'ATGA', 'ATTA', 'CAAA', 'CAAC', 'CAAG', 'CAAT', 'CACA', 'CAGA', 'CATA', 'CCAA',
            'CGAA', 'CTAA', 'GAAA', 'GAAC', 'GAAG', 'GAAT', 'GACA', 'GAGA', 'GATA', 'GCAA', 'GGAA', 'GTAA', 'TAAA', 'TAAC', 'TAAG',
            'TAAT', 'TACA', 'TAGA', 'TATA', 'TCAA', 'TGAA', 'TTAA')

tgta_rich <- c('TGTA')

dt_atha <- read_delim("./analysis/polyaclassifier_motif_importance_heatmaps/hmdf_sum.arabidopsis_thaliana.polyaclassifier_bagging3_kmers-4.txt", delim = "\t")

dt_atha['-5'] = 0
dt_atha['-4'] = 0
dt_atha['-3'] = 0
dt_atha['-2'] = 0
dt_atha['-1'] = 0
dt_atha['0'] = 0
dt_atha['1'] = 0

dt_atha_corr <- dt_atha %>% column_to_rownames('testMotif') %>% t %>% cor

paletteLength <- 50

ann <- dt_atha['testMotif'] %>% 
  mutate(motifFamily = case_when(
    testMotif %in% tgta_rich ~ 'TGTA-containing',
    testMotif %in% a_rich ~ 'A-rich',
    testMotif %in% t_rich ~ 'T-rich',
    TRUE ~ 'Other'),
    a_count = str_count(testMotif, "A"),
    c_count = str_count(testMotif, "C"),
    g_count = str_count(testMotif, "G"),
    t_count = str_count(testMotif, "T"),
  ) %>% 
  column_to_rownames('testMotif')

ann %>% group_by(motifFamily) %>% tally() %>% 
  mutate(prop = n / sum(n) * 100)

ann_colors <- list(
  'a_count' = brewer.pal(n = 7, name = "Greys"),
  'c_count' = brewer.pal(n = 7, name = "Greys"),
  'g_count' = brewer.pal(n = 7, name = "Greys"),
  't_count' = brewer.pal(n = 7, name = "Greys"),
  'motifFamily' = c('A-rich' = '#2278b5', 'T-rich' = '#2fa148', 'TGTA-containing' = '#fcb316', 'Other' = '#f7f8f8')
)

myColor <- colorRampPalette(c("seagreen", "white", "mediumpurple"))(paletteLength)
myBreaks <- c(seq(-1, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(1/paletteLength, 1, length.out=floor(paletteLength/2)))

callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = -sv)
  as.hclust(dend)
}

g1 <- dt_atha_corr %>% pheatmap(clustering_method = "ward.D2", clustering_callback = callback)

row_dend <- g1[[1]]
col_dend <- g1[[2]]

row_dend <- dendextend::rotate(row_dend, order = rev(rownames(dt_atha_corr)[seriation::get_order(col_dend)]))
col_dend <- dendextend::rotate(col_dend, order = rownames(dt_atha_corr)[seriation::get_order(row_dend)])

g2 <- dt_atha_corr %>% 
  pheatmap(color = myColor, breaks = myBreaks,
           cutree_cols = 3, cutree_rows = 3,
           show_rownames = FALSE, show_colnames = FALSE, border_color = NA,
           annotation_col = ann['motifFamily'], 
           annotation_row = ann[c('motifFamily','t_count','g_count','c_count','a_count')],
           annotation_colors = ann_colors,
           cluster_cols=as.hclust(col_dend),
           cluster_rows=as.hclust(row_dend)
           )

ggsave("./analysis/polyaclassifier_motif_importance_heatmaps/atha.4mers.hm_correlation.svg", plot = g2, device = svg())

myMin <- -300
myMax <- +300

myColor <- colorRampPalette(c("steelblue", "white", "firebrick"))(paletteLength)
myBreaks <- c(seq(myMin, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(myMax/paletteLength, myMax, length.out=floor(paletteLength/2)))

g3 <- dt_atha[match(rownames(dt_atha_corr[g2$tree_row$order,]), dt_atha$testMotif),] %>%
  column_to_rownames('testMotif') %>%
  pheatmap(cluster_rows = FALSE, cluster_cols = FALSE,
           show_rownames = TRUE, show_colnames = TRUE, border_color = NA,
           color = myColor, breaks = myBreaks,
           annotation_row = ann['motifFamily'],
           annotation_colors = ann_colors
  )

ggsave("./analysis/polyaclassifier_motif_importance_heatmaps/atha.4mers.hm_importance.reordered.svg", plot = g3, device = svg())

row_order_names <- rownames(dt_atha_corr[g2$tree_row$order,])
row_order_names


other_motifs <- ann %>% filter(motifFamily == "Other") %>% rownames %>% dput

cluster_motifs <- g2$tree_row %>% cutree(3) %>% as.data.frame %>% rownames_to_column(var = "testMotif") %>% rename("cluster" = ".")
cluster_summary <- cluster_motifs %>% group_by(cluster) %>% summarize(paste0(testMotif, collapse = ", "))
names(cluster_summary) <- c("cluster","motifs")

cluster_list <- list()
cluster_list[[1]] <- cluster_summary %>% filter(cluster == 1) %>% pull(motifs)
cluster_list[[2]] <- cluster_summary %>% filter(cluster == 2) %>% pull(motifs)
cluster_list[[3]] <- cluster_summary %>% filter(cluster == 3) %>% pull(motifs)
cluster_list


