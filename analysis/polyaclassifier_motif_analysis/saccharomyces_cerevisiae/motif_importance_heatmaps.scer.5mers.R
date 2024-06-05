
## Import packages ------------------------------------------------------------

packages <- c("tidyverse", "RColorBrewer", "pheatmap", "gridExtra", "ggrepel", "magrittr", "matrixStats","cowplot","viridis","dendextend")
lapply(packages, library, character.only = TRUE)

## Visualize heatmaps: CEREVISIAE ---------------------------------------------

a_rich_d0 <- c('AAAAA')

a_rich_d1 <- c('AAAAC', 'AAAAG', 'AAAAT', 'AAACA', 'AAAGA', 'AACAA', 'AAGAA', 'AATAA', 'ACAAA', 'AGAAA', 'CAAAA', 'GAAAA', 'TAAAA')

a_rich_d2 <- c('AAACC', 'AAACG', 'AAACT', 'AAAGC', 'AAAGG', 'AAAGT', 'AACAC', 'AACAG', 'AACAT', 'AACCA', 'AACGA', 'AACTA', 'AAGAC', 'AAGAG',
               'AAGAT', 'AAGCA', 'AAGGA', 'AAGTA', 'AATCA', 'AATGA', 'ACAAC', 'ACAAG', 'ACAAT', 'ACCAA', 'ACGAA', 'ACTAA', 'AGAAC', 'AGAAG',
               'AGAAT', 'AGCAA', 'AGGAA', 'AGTAA', 'CAAAC', 'CAAAG', 'CAAAT', 'CAACA', 'CAAGA', 'CACAA', 'CAGAA', 'CCAAA', 'CGAAA', 'CTAAA',
               'GAAAC', 'GAAAG', 'GAAAT', 'GAACA', 'GAAGA', 'GACAA', 'GAGAA', 'GCAAA', 'GGAAA', 'GTAAA', 'TAACA', 'TAAGA', 'TCAAA', 'TGAAA',
               'TTAAA')

t_rich_d0 <- c('TTTTT')

t_rich_d1 <- c('ATTTT', 'CTTTT', 'GTTTT', 'TCTTT', 'TGTTT', 'TTATT', 'TTCTT', 'TTGTT', 'TTTCT', 'TTTGT', 'TTTTA', 'TTTTC', 'TTTTG')

t_rich_d2 <- c('AATTT', 'ACTTT', 'AGTTT', 'ATTCT', 'ATTGT', 'CATTT', 'CCTTT', 'CGTTT', 'CTCTT', 'CTGTT', 'CTTCT', 'CTTGT', 'CTTTA', 'CTTTC',
               'CTTTG', 'GATTT', 'GCTTT', 'GGTTT', 'GTCTT', 'GTGTT', 'GTTCT', 'GTTGT', 'GTTTA', 'GTTTC', 'GTTTG', 'TCATT', 'TCCTT', 'TCGTT',
               'TCTTA', 'TCTTC', 'TCTTG', 'TGATT', 'TGCTT', 'TGGTT', 'TGTTA', 'TGTTC', 'TGTTG', 'TTACT', 'TTAGT', 'TTCAT', 'TTCCT', 'TTCGT',
               'TTCTA', 'TTCTC', 'TTCTG', 'TTGAT', 'TTGCT', 'TTGGT', 'TTGTA', 'TTGTC', 'TTGTG', 'TTTCA', 'TTTCC', 'TTTCG', 'TTTGA', 'TTTGC',
               'TTTGG')

tata_rich_d0 <- c('ATATA', 'TATAT')

tata_rich_d1 <- c('AAATA', 'AATAT', 'ACATA', 'AGATA', 'ATAAA', 'ATACA', 'ATAGA', 'ATATC', 'ATATG', 'ATATT', 'ATCTA', 'ATGTA', 'ATTTA', 'CATAT',
                  'CTATA', 'GATAT', 'GTATA', 'TAAAT', 'TACAT', 'TAGAT', 'TATAA', 'TATAC', 'TATAG', 'TATCT', 'TATGT', 'TATTT', 'TCTAT', 'TGTAT',
                  'TTATA', 'TTTAT')

tata_rich_d2 <- c('AAATC', 'AAATG', 'AAATT', 'AATAC', 'AATAG', 'AATCT', 'AATGT', 'AATTA', 'ACACA', 'ACAGA', 'ACATC', 'ACATG', 'ACATT', 'ACCTA',
                  'ACGTA', 'ACTAT', 'ACTTA', 'AGACA', 'AGAGA', 'AGATC', 'AGATG', 'AGATT', 'AGCTA', 'AGGTA', 'AGTAT', 'AGTTA', 'ATAAC', 'ATAAG',
                  'ATAAT', 'ATACC', 'ATACG', 'ATACT', 'ATAGC', 'ATAGG', 'ATAGT', 'ATCAA', 'ATCCA', 'ATCGA', 'ATCTC', 'ATCTG', 'ATCTT', 'ATGAA',
                  'ATGCA', 'ATGGA', 'ATGTC', 'ATGTG', 'ATGTT', 'ATTAA', 'ATTAT', 'ATTCA', 'ATTGA', 'ATTTC', 'ATTTG', 'CAATA', 'CACAT', 'CAGAT',
                  'CATAA', 'CATAC', 'CATAG', 'CATCT', 'CATGT', 'CCATA', 'CCTAT', 'CGATA', 'CGTAT', 'CTACA', 'CTAGA', 'CTATC', 'CTATG', 'CTATT',
                  'CTCTA', 'CTGTA', 'CTTAT', 'GAATA', 'GACAT', 'GAGAT', 'GATAA', 'GATAC', 'GATAG', 'GATCT', 'GATGT', 'GCATA', 'GCTAT', 'GGATA',
                  'GGTAT', 'GTACA', 'GTAGA', 'GTATC', 'GTATG', 'GTATT', 'GTCTA', 'GTGTA', 'GTTAT', 'TAAAC', 'TAAAG', 'TAACT', 'TAAGT', 'TAATA',
                  'TAATT', 'TACAA', 'TACAC', 'TACAG', 'TACCT', 'TACGT', 'TACTT', 'TAGAA', 'TAGAC', 'TAGAG', 'TAGCT', 'TAGGT', 'TAGTT', 'TATCA',
                  'TATCC', 'TATCG', 'TATGA', 'TATGC', 'TATGG', 'TATTA', 'TATTC', 'TATTG', 'TCAAT', 'TCATA', 'TCCAT', 'TCGAT', 'TCTAA', 'TCTAC',
                  'TCTAG', 'TCTCT', 'TCTGT', 'TGAAT', 'TGATA', 'TGCAT', 'TGGAT', 'TGTAA', 'TGTAC', 'TGTAG', 'TGTCT', 'TGTGT', 'TTAAT', 'TTACA',
                  'TTAGA', 'TTATC', 'TTATG', 'TTTAA', 'TTTAC', 'TTTAG')

dt_cer <- read_delim("./analysis/polyaclassifier_motif_importance_heatmaps/hmdf_sum.saccharomyces_cerevisiae.polyaclassifier_bagging3_kmers-5.txt", delim = "\t")

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
           show_rownames = FALSE, show_colnames = FALSE, border_color = NA,
           annotation_col = ann['motifFamily'], 
           annotation_row = ann[c('motifFamily','t_count','g_count','c_count','a_count')],
           annotation_colors = ann_colors,
           cluster_cols=as.hclust(col_dend),
           cluster_rows=as.hclust(row_dend),
           )

ggsave("./analysis/polyaclassifier_motif_importance_heatmaps/scer.5mers.hm_correlation.svg", plot = g2, device = svg())

myMin <- -300
myMax <- +300

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

ggsave("./analysis/polyaclassifier_motif_importance_heatmaps/scer.5mers.hm_importance.reordered.svg", plot = g3, device = svg())

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


