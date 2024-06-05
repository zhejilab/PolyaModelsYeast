
## Import packages ------------------------------------------------------------

packages <- c("tidyverse", "RColorBrewer", "pheatmap", "gridExtra", "ggrepel", "magrittr", "matrixStats","cowplot","viridis","dendextend")
lapply(packages, library, character.only = TRUE)

## Visualize heatmaps: CEREVISIAE ---------------------------------------------

a_rich_d0 <- c('AAAAAA')

a_rich_d1 <- c('AAAAAC', 'AAAAAG', 'AAAAAT', 'AAAACA', 'AAAAGA', 'AAAATA', 'AAACAA', 'AAAGAA', 'AAATAA', 'AACAAA', 'AAGAAA', 'AATAAA',
               'ACAAAA', 'AGAAAA', 'ATAAAA', 'CAAAAA', 'GAAAAA', 'TAAAAA')

a_rich_d2 <- c('AAAACC', 'AAAACG', 'AAAACT', 'AAAAGC', 'AAAAGG', 'AAAAGT', 'AAAATC', 'AAAATG', 'AAAATT', 'AAACAC', 'AAACAG', 'AAACAT',
               'AAACCA', 'AAACGA', 'AAACTA', 'AAAGAC', 'AAAGAG', 'AAAGAT', 'AAAGCA', 'AAAGGA', 'AAAGTA', 'AAATCA', 'AAATGA', 'AAATTA',
               'AACAAC', 'AACAAG', 'AACAAT', 'AACACA', 'AACAGA', 'AACCAA', 'AACGAA', 'AACTAA', 'AAGAAC', 'AAGAAG', 'AAGAAT', 'AAGACA',
               'AAGAGA', 'AAGCAA', 'AAGGAA', 'AAGTAA', 'AATAAC', 'AATAAG', 'AATAAT', 'AATCAA', 'AATGAA', 'AATTAA', 'ACAAAC', 'ACAAAG',
               'ACAAAT', 'ACAACA', 'ACAAGA', 'ACAATA', 'ACACAA', 'ACAGAA', 'ACCAAA', 'ACGAAA', 'ACTAAA', 'AGAAAC', 'AGAAAG', 'AGAAAT',
               'AGAACA', 'AGAAGA', 'AGAATA', 'AGACAA', 'AGAGAA', 'AGCAAA', 'AGGAAA', 'AGTAAA', 'ATAACA', 'ATAAGA', 'ATAATA', 'ATCAAA',
               'ATGAAA', 'ATTAAA', 'CAAAAC', 'CAAAAG', 'CAAAAT', 'CAAACA', 'CAAAGA', 'CAAATA', 'CAACAA', 'CAAGAA', 'CAATAA', 'CACAAA',
               'CAGAAA', 'CATAAA', 'CCAAAA', 'CGAAAA', 'CTAAAA', 'GAAAAC', 'GAAAAG', 'GAAAAT', 'GAAACA', 'GAAAGA', 'GAAATA', 'GAACAA',
               'GAAGAA', 'GAATAA', 'GACAAA', 'GAGAAA', 'GATAAA', 'GCAAAA', 'GGAAAA', 'GTAAAA', 'TAAAAC', 'TAAAAG', 'TAAAAT', 'TAAACA',
               'TAAAGA', 'TAACAA', 'TAAGAA', 'TAATAA', 'TACAAA', 'TAGAAA', 'TCAAAA', 'TGAAAA', 'TTAAAA')

t_rich_d0 <- c('TTTTTT')

t_rich_d1 <- c('ATTTTT', 'CTTTTT', 'GTTTTT', 'TATTTT', 'TCTTTT', 'TGTTTT', 'TTATTT', 'TTCTTT', 'TTGTTT', 'TTTATT', 'TTTCTT', 'TTTGTT',
               'TTTTAT', 'TTTTCT', 'TTTTGT', 'TTTTTA', 'TTTTTC', 'TTTTTG')

t_rich_d2 <- c('AATTTT', 'ACTTTT', 'AGTTTT', 'ATCTTT', 'ATGTTT', 'ATTATT', 'ATTCTT', 'ATTGTT', 'ATTTCT', 'ATTTGT', 'ATTTTA', 'ATTTTC',
               'ATTTTG', 'CATTTT', 'CCTTTT', 'CGTTTT', 'CTATTT', 'CTCTTT', 'CTGTTT', 'CTTATT', 'CTTCTT', 'CTTGTT', 'CTTTAT', 'CTTTCT',
               'CTTTGT', 'CTTTTA', 'CTTTTC', 'CTTTTG', 'GATTTT', 'GCTTTT', 'GGTTTT', 'GTATTT', 'GTCTTT', 'GTGTTT', 'GTTATT', 'GTTCTT',
               'GTTGTT', 'GTTTAT', 'GTTTCT', 'GTTTGT', 'GTTTTA', 'GTTTTC', 'GTTTTG', 'TAATTT', 'TACTTT', 'TAGTTT', 'TATTAT', 'TATTCT',
               'TATTGT', 'TCATTT', 'TCCTTT', 'TCGTTT', 'TCTCTT', 'TCTGTT', 'TCTTAT', 'TCTTCT', 'TCTTGT', 'TCTTTA', 'TCTTTC', 'TCTTTG',
               'TGATTT', 'TGCTTT', 'TGGTTT', 'TGTCTT', 'TGTGTT', 'TGTTAT', 'TGTTCT', 'TGTTGT', 'TGTTTA', 'TGTTTC', 'TGTTTG', 'TTAATT',
               'TTACTT', 'TTAGTT', 'TTATTA', 'TTATTC', 'TTATTG', 'TTCATT', 'TTCCTT', 'TTCGTT', 'TTCTCT', 'TTCTGT', 'TTCTTA', 'TTCTTC',
               'TTCTTG', 'TTGATT', 'TTGCTT', 'TTGGTT', 'TTGTCT', 'TTGTGT', 'TTGTTA', 'TTGTTC', 'TTGTTG', 'TTTAAT', 'TTTACT', 'TTTAGT',
               'TTTCAT', 'TTTCCT', 'TTTCGT', 'TTTCTA', 'TTTCTC', 'TTTCTG', 'TTTGAT', 'TTTGCT', 'TTTGGT', 'TTTGTA', 'TTTGTC', 'TTTGTG',
               'TTTTAA', 'TTTTAC', 'TTTTAG', 'TTTTCA', 'TTTTCC', 'TTTTCG', 'TTTTGA', 'TTTTGC', 'TTTTGG')

tata_rich_d0 <- c('ATATAT', 'TATATA')

tata_rich_d1 <- c('AAATAT', 'AATATA', 'ACATAT', 'AGATAT', 'ATAAAT', 'ATACAT', 'ATAGAT', 'ATATAA', 'ATATAC', 'ATATAG', 'ATATCT', 'ATATGT',
                  'ATATTT', 'ATCTAT', 'ATGTAT', 'ATTTAT', 'CATATA', 'CTATAT', 'GATATA', 'GTATAT', 'TAAATA', 'TACATA', 'TAGATA', 'TATAAA',
                  'TATACA', 'TATAGA', 'TATATC', 'TATATG', 'TATATT', 'TATCTA', 'TATGTA', 'TATTTA', 'TCTATA', 'TGTATA', 'TTATAT', 'TTTATA')

tata_rich_d2 <- c('AAATAC', 'AAATAG', 'AAATCT', 'AAATGT', 'AAATTT', 'AACATA', 'AACTAT', 'AAGATA', 'AAGTAT', 'AATACA', 'AATAGA', 'AATATC',
                  'AATATG', 'AATATT', 'AATCTA', 'AATGTA', 'AATTAT', 'AATTTA', 'ACACAT', 'ACAGAT', 'ACATAA', 'ACATAC', 'ACATAG', 'ACATCT',
                  'ACATGT', 'ACATTT', 'ACCTAT', 'ACGTAT', 'ACTATA', 'ACTTAT', 'AGACAT', 'AGAGAT', 'AGATAA', 'AGATAC', 'AGATAG', 'AGATCT',
                  'AGATGT', 'AGATTT', 'AGCTAT', 'AGGTAT', 'AGTATA', 'AGTTAT', 'ATAAAC', 'ATAAAG', 'ATAACT', 'ATAAGT', 'ATAATT', 'ATACAA',
                  'ATACAC', 'ATACAG', 'ATACCT', 'ATACGT', 'ATACTT', 'ATAGAA', 'ATAGAC', 'ATAGAG', 'ATAGCT', 'ATAGGT', 'ATAGTT', 'ATATCA',
                  'ATATCC', 'ATATCG', 'ATATGA', 'ATATGC', 'ATATGG', 'ATATTA', 'ATATTC', 'ATATTG', 'ATCAAT', 'ATCCAT', 'ATCGAT', 'ATCTAA',
                  'ATCTAC', 'ATCTAG', 'ATCTCT', 'ATCTGT', 'ATGAAT', 'ATGCAT', 'ATGGAT', 'ATGTAA', 'ATGTAC', 'ATGTAG', 'ATGTCT', 'ATGTGT',
                  'ATTAAT', 'ATTATA', 'ATTCAT', 'ATTGAT', 'ATTTAA', 'ATTTAC', 'ATTTAG', 'CAATAT', 'CACATA', 'CAGATA', 'CATACA', 'CATAGA',
                  'CATATC', 'CATATG', 'CATATT', 'CATCTA', 'CATGTA', 'CATTTA', 'CCATAT', 'CCTATA', 'CGATAT', 'CGTATA', 'CTAAAT', 'CTACAT',
                  'CTAGAT', 'CTATAA', 'CTATAC', 'CTATAG', 'CTATCT', 'CTATGT', 'CTCTAT', 'CTGTAT', 'CTTATA', 'GAATAT', 'GACATA', 'GAGATA',
                  'GATACA', 'GATAGA', 'GATATC', 'GATATG', 'GATATT', 'GATCTA', 'GATGTA', 'GATTTA', 'GCATAT', 'GCTATA', 'GGATAT', 'GGTATA',
                  'GTAAAT', 'GTACAT', 'GTAGAT', 'GTATAA', 'GTATAC', 'GTATAG', 'GTATCT', 'GTATGT', 'GTCTAT', 'GTGTAT', 'GTTATA', 'TAAATC',
                  'TAAATG', 'TAAATT', 'TAACTA', 'TAAGTA', 'TAATAT', 'TAATTA', 'TACACA', 'TACAGA', 'TACATC', 'TACATG', 'TACATT', 'TACCTA',
                  'TACGTA', 'TACTTA', 'TAGACA', 'TAGAGA', 'TAGATC', 'TAGATG', 'TAGATT', 'TAGCTA', 'TAGGTA', 'TAGTTA', 'TATAAC', 'TATAAG',
                  'TATAAT', 'TATACC', 'TATACG', 'TATACT', 'TATAGC', 'TATAGG', 'TATAGT', 'TATCAA', 'TATCCA', 'TATCGA', 'TATCTC', 'TATCTG',
                  'TATCTT', 'TATGAA', 'TATGCA', 'TATGGA', 'TATGTC', 'TATGTG', 'TATGTT', 'TATTAA', 'TATTCA', 'TATTGA', 'TATTTC', 'TATTTG',
                  'TCAATA', 'TCATAT', 'TCCATA', 'TCGATA', 'TCTAAA', 'TCTACA', 'TCTAGA', 'TCTATC', 'TCTATG', 'TCTATT', 'TCTCTA', 'TCTGTA',
                  'TGAATA', 'TGATAT', 'TGCATA', 'TGGATA', 'TGTAAA', 'TGTACA', 'TGTAGA', 'TGTATC', 'TGTATG', 'TGTATT', 'TGTCTA', 'TGTGTA',
                  'TTAAAT', 'TTAATA', 'TTACAT', 'TTAGAT', 'TTATAA', 'TTATAC', 'TTATAG', 'TTATCT', 'TTATGT', 'TTCATA', 'TTCTAT', 'TTGATA',
                  'TTGTAT', 'TTTAAA', 'TTTACA', 'TTTAGA', 'TTTATC', 'TTTATG')

dt_cer <- read_delim("./analysis/polyaclassifier_motif_importance_heatmaps/hmdf_sum.saccharomyces_cerevisiae.polyaclassifier_bagging3_kmers-6.txt", delim = "\t")

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

g1 <- dt_cer_corr %>% pheatmap(clustering_method = "ward.D2", ## options: average, single, complete, median, centroid, mcquitty, ward.D, ward.D2
                               clustering_callback = callback
                               )

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

ggsave("./analysis/polyaclassifier_motif_importance_heatmaps/scer.6mers.hm_correlation.svg", plot = g2, device = svg())

myMin <- -100
myMax <- +100

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

ggsave("./analysis/polyaclassifier_motif_importance_heatmaps/scer.6mers.hm_importance.reordered.svg", plot = g3, device = svg())

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

cluster_families <- g2$tree_row %>% 
  cutree(3) %>% as.data.frame %>% 
  rownames_to_column(var = "motif") %>% rename("cluster" = ".") %>% 
  mutate(motifFamily = case_when(
    motif %in% tata_rich_d0 ~ 'TA/TA-rich_d0',
    motif %in% tata_rich_d1 ~ 'TA/TA-rich_d1',
    motif %in% tata_rich_d2 ~ 'TA/TA-rich_d2',
    motif %in% t_rich_d0 ~ 'T-rich_d0',
    motif %in% t_rich_d1 ~ 'T-rich_d1',
    motif %in% t_rich_d2 ~ 'T-rich_d2',
    motif %in% a_rich_d0 ~ 'A-rich_d0',
    motif %in% a_rich_d1 ~ 'A-rich_d1',
    motif %in% a_rich_d2 ~ 'A-rich_d2',
    TRUE ~ 'Other'))
  
cluster_families %>% arrange(cluster, motifFamily)



