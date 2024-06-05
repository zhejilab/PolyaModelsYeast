
## Import packages ------------------------------------------------------------

packages <- c("tidyverse", "RColorBrewer", "pheatmap", "gridExtra", "ggrepel", "magrittr", "matrixStats","cowplot","viridis","dendextend")
lapply(packages, library, character.only = TRUE)

## Visualize heatmaps: POMBE --------------------------------------------------

t_rich_d0 <- c('TTTTT')

t_rich_d1 <- c('ATTTT', 'CTTTT', 'GTTTT', 'TATTT', 'TCTTT', 'TGTTT', 'TTATT', 'TTCTT', 'TTGTT', 'TTTAT', 'TTTCT', 'TTTGT', 'TTTTA', 'TTTTC',
               'TTTTG')

t_rich_d2 <- c('AATTT', 'ACTTT', 'AGTTT', 'ATATT', 'ATCTT', 'ATGTT', 'ATTAT', 'ATTCT', 'ATTGT', 'ATTTA', 'ATTTC', 'ATTTG', 'CATTT', 'CCTTT',
               'CGTTT', 'CTATT', 'CTCTT', 'CTGTT', 'CTTAT', 'CTTCT', 'CTTGT', 'CTTTA', 'CTTTC', 'CTTTG', 'GATTT', 'GCTTT', 'GGTTT', 'GTCTT',
               'GTGTT', 'GTTAT', 'GTTCT', 'GTTGT', 'GTTTA', 'GTTTC', 'GTTTG', 'TAATT', 'TACTT', 'TAGTT', 'TATAT', 'TATCT', 'TATGT', 'TATTA',
               'TATTC', 'TATTG', 'TCATT', 'TCCTT', 'TCGTT', 'TCTAT', 'TCTCT', 'TCTGT', 'TCTTA', 'TCTTC', 'TCTTG', 'TGATT', 'TGCTT', 'TGGTT',
               'TGTCT', 'TGTGT', 'TGTTA', 'TGTTC', 'TGTTG', 'TTAAT', 'TTACT', 'TTAGT', 'TTATA', 'TTATC', 'TTATG', 'TTCAT', 'TTCCT', 'TTCGT',
               'TTCTA', 'TTCTC', 'TTCTG', 'TTGAT', 'TTGCT', 'TTGGT', 'TTGTC', 'TTGTG', 'TTTAA', 'TTTAC', 'TTTAG', 'TTTCA', 'TTTCC', 'TTTCG',
               'TTTGA', 'TTTGC', 'TTTGG')

a_rich_d0 <- c('AAAAA')

a_rich_d1 <- c('AAAAC', 'AAAAG', 'AAAAT', 'AAACA', 'AAAGA', 'AAATA', 'AACAA', 'AAGAA', 'AATAA', 'ACAAA', 'AGAAA', 'ATAAA', 'CAAAA', 'GAAAA',
               'TAAAA')

a_rich_d2 <- c('AAACC', 'AAACG', 'AAACT', 'AAAGC', 'AAAGG', 'AAAGT', 'AAATC', 'AAATG', 'AAATT', 'AACAC', 'AACAG', 'AACAT', 'AACCA', 'AACGA',
               'AACTA', 'AAGAC', 'AAGAG', 'AAGAT', 'AAGCA', 'AAGGA', 'AATAC', 'AATAG', 'AATAT', 'AATCA', 'AATGA', 'AATTA', 'ACAAC', 'ACAAG',
               'ACAAT', 'ACACA', 'ACAGA', 'ACATA', 'ACCAA', 'ACGAA', 'ACTAA', 'AGAAC', 'AGAAG', 'AGAAT', 'AGACA', 'AGAGA', 'AGATA', 'AGCAA',
               'AGGAA', 'ATAAC', 'ATAAG', 'ATAAT', 'ATACA', 'ATAGA', 'ATATA', 'ATCAA', 'ATGAA', 'ATTAA', 'CAAAC', 'CAAAG', 'CAAAT', 'CAACA',
               'CAAGA', 'CAATA', 'CACAA', 'CAGAA', 'CATAA', 'CCAAA', 'CGAAA', 'CTAAA', 'GAAAC', 'GAAAG', 'GAAAT', 'GAACA', 'GAAGA', 'GAATA',
               'GACAA', 'GAGAA', 'GATAA', 'GCAAA', 'GGAAA', 'TAAAC', 'TAAAG', 'TAAAT', 'TAACA', 'TAAGA', 'TAATA', 'TACAA', 'TAGAA', 'TATAA',
               'TCAAA', 'TGAAA', 'TTAAA')

a_rich_d3 <- c('AACCC', 'AACCG', 'AACCT', 'AACGC', 'AACGG', 'AACGT', 'AACTC', 'AACTG', 'AACTT', 'AAGCC', 'AAGCG', 'AAGCT', 'AAGGC', 'AAGGG',
               'AAGGT', 'AAGTC', 'AAGTG', 'AAGTT', 'AATCC', 'AATCG', 'AATCT', 'AATGC', 'AATGG', 'AATGT', 'AATTC', 'AATTG', 'ACACC', 'ACACG',
               'ACACT', 'ACAGC', 'ACAGG', 'ACAGT', 'ACATC', 'ACATG', 'ACATT', 'ACCAC', 'ACCAG', 'ACCAT', 'ACCCA', 'ACCGA', 'ACCTA', 'ACGAC',
               'ACGAG', 'ACGAT', 'ACGCA', 'ACGGA', 'ACTAC', 'ACTAG', 'ACTAT', 'ACTCA', 'ACTGA', 'ACTTA', 'AGACC', 'AGACG', 'AGACT', 'AGAGC',
               'AGAGG', 'AGAGT', 'AGATC', 'AGATG', 'AGATT', 'AGCAC', 'AGCAG', 'AGCAT', 'AGCCA', 'AGCGA', 'AGCTA', 'AGGAC', 'AGGAG', 'AGGAT',
               'AGGCA', 'AGGGA', 'AGTCA', 'AGTGA', 'AGTTA', 'ATACC', 'ATACG', 'ATACT', 'ATAGC', 'ATAGG', 'ATAGT', 'ATATC', 'ATATG', 'ATCAC',
               'ATCAG', 'ATCAT', 'ATCCA', 'ATCGA', 'ATCTA', 'ATGAC', 'ATGAG', 'ATGAT', 'ATGCA', 'ATGGA', 'ATTAC', 'ATTAG', 'ATTCA', 'ATTGA',
               'CAACC', 'CAACG', 'CAACT', 'CAAGC', 'CAAGG', 'CAAGT', 'CAATC', 'CAATG', 'CAATT', 'CACAC', 'CACAG', 'CACAT', 'CACCA', 'CACGA',
               'CACTA', 'CAGAC', 'CAGAG', 'CAGAT', 'CAGCA', 'CAGGA', 'CATAC', 'CATAG', 'CATAT', 'CATCA', 'CATGA', 'CATTA', 'CCAAC', 'CCAAG',
               'CCAAT', 'CCACA', 'CCAGA', 'CCATA', 'CCCAA', 'CCGAA', 'CCTAA', 'CGAAC', 'CGAAG', 'CGAAT', 'CGACA', 'CGAGA', 'CGATA', 'CGCAA',
               'CGGAA', 'CTAAC', 'CTAAG', 'CTAAT', 'CTACA', 'CTAGA', 'CTATA', 'CTCAA', 'CTGAA', 'CTTAA', 'GAACC', 'GAACG', 'GAACT', 'GAAGC',
               'GAAGG', 'GAAGT', 'GAATC', 'GAATG', 'GAATT', 'GACAC', 'GACAG', 'GACAT', 'GACCA', 'GACGA', 'GACTA', 'GAGAC', 'GAGAG', 'GAGAT',
               'GAGCA', 'GAGGA', 'GATAC', 'GATAG', 'GATAT', 'GATCA', 'GATGA', 'GATTA', 'GCAAC', 'GCAAG', 'GCAAT', 'GCACA', 'GCAGA', 'GCATA',
               'GCCAA', 'GCGAA', 'GCTAA', 'GGAAC', 'GGAAG', 'GGAAT', 'GGACA', 'GGAGA', 'GGATA', 'GGCAA', 'GGGAA', 'GTCAA', 'GTGAA', 'GTTAA',
               'TAACC', 'TAACG', 'TAACT', 'TAAGC', 'TAAGG', 'TAAGT', 'TAATC', 'TAATG', 'TACAC', 'TACAG', 'TACAT', 'TACCA', 'TACGA', 'TACTA',
               'TAGAC', 'TAGAG', 'TAGAT', 'TAGCA', 'TAGGA', 'TATAC', 'TATAG', 'TATCA', 'TATGA', 'TCAAC', 'TCAAG', 'TCAAT', 'TCACA', 'TCAGA',
               'TCATA', 'TCCAA', 'TCGAA', 'TCTAA', 'TGAAC', 'TGAAG', 'TGAAT', 'TGACA', 'TGAGA', 'TGATA', 'TGCAA', 'TGGAA', 'TTAAC', 'TTAAG',
               'TTACA', 'TTAGA', 'TTCAA', 'TTGAA')

tag_d0 <- c('AATAG', 'ACTAG', 'AGTAG', 'ATAGA', 'ATAGC', 'ATAGG', 'ATAGT', 'ATTAG', 'CATAG', 'CCTAG', 'CGTAG', 'CTAGA', 'CTAGC', 'CTAGG',
            'CTAGT', 'CTTAG', 'GATAG', 'GCTAG', 'GGTAG', 'GTAGA', 'GTAGC', 'GTAGG', 'GTAGT', 'GTTAG', 'TAGAA', 'TAGAC', 'TAGAG', 'TAGAT',
            'TAGCA', 'TAGCC', 'TAGCG', 'TAGCT', 'TAGGA', 'TAGGC', 'TAGGG', 'TAGGT', 'TAGTA', 'TAGTC', 'TAGTG', 'TAGTT', 'TATAG', 'TCTAG',
            'TGTAG', 'TTAGA', 'TTAGC', 'TTAGG', 'TTAGT', 'TTTAG')

gta_d0 <- c('AAGTA', 'ACGTA', 'AGGTA', 'AGTAA', 'AGTAC', 'AGTAG', 'AGTAT', 'ATGTA', 'CAGTA', 'CCGTA', 'CGGTA', 'CGTAA', 'CGTAC', 'CGTAG',
            'CGTAT', 'CTGTA', 'GAGTA', 'GCGTA', 'GGGTA', 'GGTAA', 'GGTAC', 'GGTAG', 'GGTAT', 'GTAAA', 'GTAAC', 'GTAAG', 'GTAAT', 'GTACA',
            'GTACC', 'GTACG', 'GTACT', 'GTAGA', 'GTAGC', 'GTAGG', 'GTAGT', 'GTATA', 'GTATC', 'GTATG', 'GTATT', 'GTGTA', 'TAGTA', 'TCGTA',
            'TGGTA', 'TGTAA', 'TGTAC', 'TGTAG', 'TGTAT', 'TTGTA')

gta_tag_d0 <- c('AGTAG', 'CGTAG', 'GGTAG', 'GTAGA', 'GTAGC', 'GTAGG', 'GTAGT', 'TAGTA', 'TGTAG')

dt_pom <- read_delim("./analysis/polyaclassifier_motif_importance_heatmaps/hmdf_sum.schizosaccharomyces_pombe.polyaclassifier_bagging3_kmers-5.txt", delim = "\t")

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
  pheatmap(
    # clustering_method = "average",
           color = myColor, breaks = myBreaks,
           cutree_cols = 4, cutree_rows = 4,
           show_rownames = FALSE, show_colnames = FALSE, border_color = NA,
           annotation_col = ann['motifFamily'], 
           annotation_row = ann[c('motifFamily','t_count','g_count','c_count','a_count')],
           annotation_colors = ann_colors,
           # clustering_callback = callback,
           cluster_cols=as.hclust(col_dend),
           cluster_rows=as.hclust(row_dend),
           )

ggsave("./analysis/polyaclassifier_motif_importance_heatmaps/spom.5mers.hm_correlation.svg", plot = g2, device = svg())

myMin <- -200
myMax <- +200

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

ggsave("./analysis/polyaclassifier_motif_importance_heatmaps/spom.5mers.hm_importance.reordered.svg", plot = g3, device = svg())

row_order_names <- rownames(dt_pom_corr[g2$tree_row$order,])
row_order_names


other_motifs <- ann %>% filter(motifFamily == "Other") %>% rownames %>% dput

cluster_motifs <- g2$tree_row %>% cutree(3) %>% as.data.frame %>% rownames_to_column(var = "motif") %>% rename("cluster" = ".") %>% group_by(cluster) %>% summarize(paste0(motif, collapse = ", "))
names(cluster_motifs) <- c("cluster","motifs")

cluster_list <- list()
cluster_list[[1]] <- cluster_motifs %>% filter(cluster == 1) %>% pull(motifs)
cluster_list[[2]] <- cluster_motifs %>% filter(cluster == 2) %>% pull(motifs)
cluster_list[[3]] <- cluster_motifs %>% filter(cluster == 3) %>% pull(motifs)
cluster_list


