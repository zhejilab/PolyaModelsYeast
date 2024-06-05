
## Import packages ------------------------------------------------------------

packages <- c("tidyverse", "RColorBrewer", "pheatmap", "gridExtra", "ggrepel", "magrittr", "matrixStats","cowplot","viridis","dendextend")
lapply(packages, library, character.only = TRUE)

## Visualize heatmaps: POMBE --------------------------------------------------

t_rich_d0 <- c('TTTT')

t_rich_d1 <- c('ATTT', 'CTTT', 'GTTT', 'TATT', 'TCTT', 'TGTT', 'TTAT', 'TTCT', 'TTGT', 'TTTA', 'TTTC', 'TTTG')

t_rich_d2 <- c('ACTT', 'AGTT', 'ATCT', 'ATGT', 'ATTA', 'ATTC', 'ATTG', 'CATT', 'CCTT', 'CGTT', 'CTAT', 'CTCT', 'CTGT', 'CTTA', 'CTTC',
               'CTTG', 'GATT', 'GCTT', 'GGTT', 'GTCT', 'GTGT', 'GTTA', 'GTTC', 'GTTG', 'TACT', 'TAGT', 'TATC', 'TATG', 'TCAT', 'TCCT',
               'TCGT', 'TCTA', 'TCTC', 'TCTG', 'TGAT', 'TGCT', 'TGGT', 'TGTC', 'TGTG', 'TTAC', 'TTAG', 'TTCA', 'TTCC', 'TTCG', 'TTGA',
               'TTGC', 'TTGG')

a_rich_d0 <- c('AAAA')

a_rich_d1 <- c('AAAC', 'AAAG', 'AAAT', 'AACA', 'AAGA', 'AATA', 'ACAA', 'AGAA', 'ATAA', 'CAAA', 'GAAA', 'TAAA')

a_rich_d2 <- c('AACC', 'AACG', 'AACT', 'AAGC', 'AAGG', 'AAGT', 'AATC', 'AATG', 'AATT', 'ACAC', 'ACAG', 'ACAT', 'ACCA', 'ACGA', 'ACTA',
               'AGAC', 'AGAG', 'AGAT', 'AGCA', 'AGGA', 'ATAC', 'ATAG', 'ATAT', 'ATCA', 'ATGA', 'CAAC', 'CAAG', 'CAAT', 'CACA', 'CAGA',
               'CATA', 'CCAA', 'CGAA', 'CTAA', 'GAAC', 'GAAG', 'GAAT', 'GACA', 'GAGA', 'GATA', 'GCAA', 'GGAA', 'TAAC', 'TAAG', 'TAAT',
               'TACA', 'TAGA', 'TATA', 'TCAA', 'TGAA', 'TTAA')

a_rich_d3 <- c('ACCC', 'ACCG', 'ACCT', 'ACGC', 'ACGG', 'ACGT', 'ACTC', 'ACTG', 'AGCC', 'AGCG', 'AGCT', 'AGGC', 'AGGG', 'AGGT', 'AGTC',
               'AGTG', 'ATCC', 'ATCG', 'ATGC', 'ATGG', 'CACC', 'CACG', 'CACT', 'CAGC', 'CAGG', 'CAGT', 'CATC', 'CATG', 'CCAC', 'CCAG',
               'CCAT', 'CCCA', 'CCGA', 'CCTA', 'CGAC', 'CGAG', 'CGAT', 'CGCA', 'CGGA', 'CTAC', 'CTAG', 'CTCA', 'CTGA', 'GACC', 'GACG',
               'GACT', 'GAGC', 'GAGG', 'GAGT', 'GATC', 'GATG', 'GCAC', 'GCAG', 'GCAT', 'GCCA', 'GCGA', 'GCTA', 'GGAC', 'GGAG', 'GGAT',
               'GGCA', 'GGGA', 'GTCA', 'GTGA', 'TACC', 'TACG', 'TAGC', 'TAGG', 'TCAC', 'TCAG', 'TCCA', 'TCGA', 'TGAC', 'TGAG', 'TGCA',
               'TGGA')

tag_d0 <- c('ATAG', 'CTAG', 'GTAG', 'TAGA', 'TAGC', 'TAGG', 'TAGT', 'TTAG')

gta_d0 <- c('AGTA', 'CGTA', 'GGTA', 'GTAA', 'GTAC', 'GTAG', 'GTAT', 'TGTA')

gta_tag_d0 <- c('GTAG')

dt_pom <- read_delim("./analysis/polyaclassifier_motif_importance_heatmaps/hmdf_sum.schizosaccharomyces_pombe.polyaclassifier_bagging3_kmers-4.txt", delim = "\t")

dt_pom['-5'] = 0
dt_pom['-4'] = 0
dt_pom['-3'] = 0
dt_pom['-2'] = 0
dt_pom['-1'] = 0
dt_pom['0'] = 0
dt_pom['1'] = 0

dt_pom_corr <- dt_pom %>% column_to_rownames('testMotif') %>% t %>% cor

paletteLength <- 50

ann <- dt_pom['testMotif'] %>% 
  mutate(motifFamily = case_when(
    testMotif %in% gta_tag_d0 ~ 'GTA-TAG_d0',
    testMotif %in% gta_d0 ~ 'GTA_d0',
    testMotif %in% tag_d0 ~ 'TAG_d0',
    testMotif %in% a_rich_d0 ~ 'A-rich_d0',
    testMotif %in% a_rich_d1 ~ 'A-rich_d1',
    testMotif %in% a_rich_d2 ~ 'A-rich_d2',
    testMotif %in% a_rich_d3 ~ 'A-rich_d3',
    testMotif %in% t_rich_d0 ~ 'T-rich_d0',
    testMotif %in% t_rich_d1 ~ 'T-rich_d1',
    testMotif %in% t_rich_d2 ~ 'T-rich_d2',
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
  'motifFamily' = c('A-rich_d0' = '#2278b5', 'A-rich_d1' = '#6ab1e3', 'A-rich_d2' = '#9ccbec', 'A-rich_d3' = '#cde5f5', 'T-rich_d0' = '#2fa148', 'T-rich_d1' = '#73d689', 'T-rich_d2' = '#a1e4b0', 'GTA_d0' = '#fcb316', 'TAG_d0' = '#8C0800', 'GTA-TAG_d0' = '#D55E00', 'Other' = '#f7f8f8')
)


myColor <- colorRampPalette(c("seagreen", "white", "mediumpurple"))(paletteLength)
myBreaks <- c(seq(-1, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(1/paletteLength, 1, length.out=floor(paletteLength/2)))

callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = -sv)
  as.hclust(dend)
}

g1 <- dt_pom_corr %>% pheatmap(clustering_method = "ward.D2", clustering_callback = callback)

row_dend <- g1[[1]]
col_dend <- g1[[2]]

row_dend <- dendextend::rotate(row_dend, order = rownames(dt_pom_corr)[seriation::get_order(col_dend)])
col_dend <- dendextend::rotate(col_dend, order = rownames(dt_pom_corr)[seriation::get_order(row_dend)])

g2 <- dt_pom_corr %>% 
  pheatmap(color = myColor, breaks = myBreaks,
           cutree_cols = 4, cutree_rows = 4,
           show_rownames = TRUE, show_colnames = FALSE, border_color = NA,
           annotation_col = ann['motifFamily'], 
           annotation_row = ann[c('motifFamily','t_count','g_count','c_count','a_count')],
           annotation_colors = ann_colors,
           cluster_cols=as.hclust(col_dend),
           cluster_rows=as.hclust(row_dend),
           )

ggsave("./analysis/polyaclassifier_motif_importance_heatmaps/spom.4mers.hm_correlation.svg", plot = g2, device = svg())

myMin <- -300
myMax <- +300

myColor <- colorRampPalette(c("steelblue", "white", "firebrick"))(paletteLength)
myBreaks <- c(seq(myMin, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(myMax/paletteLength, myMax, length.out=floor(paletteLength/2)))

g3 <- dt_pom[match(rownames(dt_pom_corr[g2$tree_row$order,]), dt_pom$testMotif),] %>%
  column_to_rownames('testMotif') %>%
  pheatmap(cluster_rows = FALSE, cluster_cols = FALSE,
           show_rownames = FALSE, show_colnames = TRUE, border_color = NA,
           color = myColor, breaks = myBreaks,
           annotation_row = ann['motifFamily'],
           annotation_colors = ann_colors
           )

ggsave("./analysis/polyaclassifier_motif_importance_heatmaps/spom.4mers.hm_importance.reordered.svg", plot = g3, device = svg())

row_order_names <- rownames(dt_pom_corr[g2$tree_row$order,])
row_order_names


other_motifs <- ann %>% filter(motifFamily == "Other") %>% rownames %>% dput

cluster_motifs <- g2$tree_row %>% cutree(4) %>% as.data.frame %>% rownames_to_column(var = "motif") %>% rename("cluster" = ".") %>% group_by(cluster) %>% summarize(paste0(motif, collapse = ", "))
names(cluster_motifs) <- c("cluster","motifs")

cluster_list <- list()
cluster_list[[1]] <- cluster_motifs %>% filter(cluster == 1) %>% pull(motifs)
cluster_list[[2]] <- cluster_motifs %>% filter(cluster == 2) %>% pull(motifs)
cluster_list[[3]] <- cluster_motifs %>% filter(cluster == 3) %>% pull(motifs)
cluster_list[[4]] <- cluster_motifs %>% filter(cluster == 4) %>% pull(motifs)
cluster_list


