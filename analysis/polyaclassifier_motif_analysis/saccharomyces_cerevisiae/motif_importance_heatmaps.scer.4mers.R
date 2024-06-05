
## Import packages ------------------------------------------------------------

packages <- c("tidyverse", "RColorBrewer", "pheatmap", "gridExtra", "ggrepel", "magrittr", "matrixStats","cowplot","viridis","dendextend")
lapply(packages, library, character.only = TRUE)

## Visualize heatmaps: CEREVISIAE ---------------------------------------------

a_rich_d0 <- c('AAAA')

a_rich_d1 <- c('AAAC', 'AAAG', 'AACA', 'AAGA', 'ACAA', 'AGAA', 'CAAA', 'GAAA')

a_rich_d2 <- c('AACC', 'AACG', 'AACT', 'AAGC', 'AAGG', 'AAGT', 'AATC', 'AATG', 'ACCA', 'ACGA', 'AGCA', 'AGGA', 'CAAC', 'CAAG', 'CCAA',
               'CGAA', 'CTAA', 'GAAC', 'GAAG', 'GCAA', 'GGAA', 'GTAA', 'TCAA', 'TGAA')

t_rich_d0 <- c('TTTT')

t_rich_d1 <- c('CTTT', 'GTTT', 'TCTT', 'TGTT', 'TTCT', 'TTGT', 'TTTC', 'TTTG')

t_rich_d2 <- c('ACTT', 'AGTT', 'CATT', 'CCTT', 'CGTT', 'CTTC', 'CTTG', 'GATT', 'GCTT', 'GGTT', 'GTTC', 'GTTG', 'TCCT', 'TCGT', 'TGCT',
               'TGGT', 'TTAC', 'TTAG', 'TTCA', 'TTCC', 'TTCG', 'TTGA', 'TTGC', 'TTGG')

tata_rich_d0 <- c('ATAT', 'TATA')

tata_rich_d1 <- c('AAAT', 'AATA', 'ACAT', 'AGAT', 'ATAA', 'ATAC', 'ATAG', 'ATCT', 'ATGT', 'ATTT', 'CATA', 'CTAT', 'GATA', 'GTAT', 'TAAA',
                  'TACA', 'TAGA', 'TATC', 'TATG', 'TATT', 'TCTA', 'TGTA', 'TTAT', 'TTTA')

tata_rich_d2 <- c('AATT', 'ACAC', 'ACAG', 'ACCT', 'ACGT', 'ACTA', 'AGAC', 'AGAG', 'AGCT', 'AGGT', 'AGTA', 'ATCA', 'ATCC', 'ATCG', 'ATGA',
                  'ATGC', 'ATGG', 'ATTA', 'ATTC', 'ATTG', 'CAAT', 'CACA', 'CAGA', 'CATC', 'CATG', 'CCAT', 'CCTA', 'CGAT', 'CGTA', 'CTAC',
                  'CTAG', 'CTCT', 'CTGT', 'CTTA', 'GAAT', 'GACA', 'GAGA', 'GATC', 'GATG', 'GCAT', 'GCTA', 'GGAT', 'GGTA', 'GTAC', 'GTAG',
                  'GTCT', 'GTGT', 'GTTA', 'TAAC', 'TAAG', 'TAAT', 'TACC', 'TACG', 'TACT', 'TAGC', 'TAGG', 'TAGT', 'TCAT', 'TCCA', 'TCGA',
                  'TCTC', 'TCTG', 'TGAT', 'TGCA', 'TGGA', 'TGTC', 'TGTG', 'TTAA')

dt_cer <- read_delim("./analysis/polyaclassifier_motif_importance_heatmaps/hmdf_sum.saccharomyces_cerevisiae.polyaclassifier_bagging3_kmers-4.txt", delim = "\t")

dt_cer['-5'] = 0
dt_cer['-4'] = 0
dt_cer['-3'] = 0
dt_cer['-2'] = 0
dt_cer['-1'] = 0
dt_cer['0'] = 0
dt_cer['1'] = 0

dt_cer_corr <- dt_cer %>% column_to_rownames('testMotif') %>% t %>% cor

paletteLength <- 50

ann <- dt_cer['testMotif'] %>% 
  mutate(motifFamily = case_when(
    testMotif %in% tata_rich_d0 ~ 'TA/TA-rich_d0',
    testMotif %in% tata_rich_d1 ~ 'TA/TA-rich_d1',
    testMotif %in% tata_rich_d2 ~ 'TA/TA-rich_d2',
    testMotif %in% t_rich_d0 ~ 'T-rich_d0',
    testMotif %in% t_rich_d1 ~ 'T-rich_d1',
    testMotif %in% t_rich_d2 ~ 'T-rich_d2',
    testMotif %in% a_rich_d0 ~ 'A-rich_d0',
    testMotif %in% a_rich_d1 ~ 'A-rich_d1',
    testMotif %in% a_rich_d2 ~ 'A-rich_d2',
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
  'motifFamily' = c('A-rich_d0' = '#2278b5', 'A-rich_d1' = '#6ab1e3', 'A-rich_d2' = '#9ccbec', 'T-rich_d0' = '#2fa148', 'T-rich_d1' = '#73d689', 'T-rich_d2' = '#a1e4b0', 'TA/TA-rich_d0' = '#d62a28', 'TA/TA-rich_d1' = '#e77f7e', 'TA/TA-rich_d2' = '#efaaa9', 'Other' = '#f7f8f8')
)


myColor <- colorRampPalette(c("seagreen", "white", "mediumpurple"))(paletteLength)
myBreaks <- c(seq(-1, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(1/paletteLength, 1, length.out=floor(paletteLength/2)))

callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = -sv)
  as.hclust(dend)
}

g1 <- dt_cer_corr %>% pheatmap(clustering_method = "ward.D2", clustering_callback = callback)

row_dend <- g1[[1]]
col_dend <- g1[[2]]
row_dend <- dendextend::rotate(row_dend, order = rev(rownames(dt_cer_corr)[seriation::get_order(col_dend)]))
col_dend <- dendextend::rotate(col_dend, order = rownames(dt_cer_corr)[seriation::get_order(row_dend)])

g2 <- dt_cer_corr %>% 
  pheatmap(color = myColor, breaks = myBreaks,
           cutree_cols = 3, cutree_rows = 3,
           show_rownames = TRUE, show_colnames = FALSE, border_color = NA,
           annotation_col = ann['motifFamily'], 
           annotation_row = ann[c('motifFamily','t_count','g_count','c_count','a_count')],
           annotation_colors = ann_colors,
           cluster_cols=as.hclust(col_dend),
           cluster_rows=as.hclust(row_dend),
           )

ggsave("./analysis/polyaclassifier_motif_importance_heatmaps/scer.4mers.hm_correlation.svg", plot = g2, device = svg())

myMin <- -1000
myMax <- +1000

myColor <- colorRampPalette(c("steelblue", "white", "firebrick"))(paletteLength)
myBreaks <- c(seq(myMin, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(myMax/paletteLength, myMax, length.out=floor(paletteLength/2)))

g3 <- dt_cer[match(rownames(dt_cer_corr[g2$tree_row$order,]), dt_cer$testMotif),] %>%
  column_to_rownames('testMotif') %>%
  pheatmap(cluster_rows = FALSE, cluster_cols = FALSE,
           show_rownames = FALSE, show_colnames = TRUE, border_color = NA,
           color = myColor, breaks = myBreaks,
           annotation_row = ann['motifFamily'],
           annotation_colors = ann_colors
           )

ggsave("./analysis/polyaclassifier_motif_importance_heatmaps/scer.4mers.hm_importance.reordered.svg", plot = g3, device = svg())

row_order_names <- rownames(dt_cer_corr[g2$tree_row$order,])
row_order_names


other_motifs <- ann %>% filter(motifFamily == "Other") %>% rownames %>% dput

cluster_motifs <- g2$tree_row %>% cutree(3) %>% as.data.frame %>% rownames_to_column(var = "motif") %>% rename("cluster" = ".") %>% group_by(cluster) %>% summarize(paste0(motif, collapse = ", "))
names(cluster_motifs) <- c("cluster","motifs")

cluster_list <- list()
cluster_list[[1]] <- cluster_motifs %>% filter(cluster == 1) %>% pull(motifs)
cluster_list[[2]] <- cluster_motifs %>% filter(cluster == 2) %>% pull(motifs)
cluster_list[[3]] <- cluster_motifs %>% filter(cluster == 3) %>% pull(motifs)
cluster_list


