
## Import packages ------------------------------------------------------------

packages <- c("tidyverse", "RColorBrewer", "pheatmap", "gridExtra", "ggrepel", "magrittr", "matrixStats","cowplot","viridis","dendextend")
lapply(packages, library, character.only = TRUE)

## Visualize heatmaps: POMBE --------------------------------------------------

t_rich_d0 <- c('TTTTTT')

t_rich_d1 <- c('ATTTTT', 'CTTTTT', 'GTTTTT', 'TATTTT', 'TCTTTT', 'TGTTTT', 'TTATTT', 'TTCTTT', 'TTGTTT', 'TTTATT', 'TTTCTT', 'TTTGTT',
               'TTTTAT', 'TTTTCT', 'TTTTGT', 'TTTTTA', 'TTTTTC', 'TTTTTG')

t_rich_d2 <- c('AATTTT', 'ACTTTT', 'AGTTTT', 'ATATTT', 'ATCTTT', 'ATGTTT', 'ATTATT', 'ATTCTT', 'ATTGTT', 'ATTTAT', 'ATTTCT', 'ATTTGT',
               'ATTTTA', 'ATTTTC', 'ATTTTG', 'CATTTT', 'CCTTTT', 'CGTTTT', 'CTATTT', 'CTCTTT', 'CTGTTT', 'CTTATT', 'CTTCTT', 'CTTGTT',
               'CTTTAT', 'CTTTCT', 'CTTTGT', 'CTTTTA', 'CTTTTC', 'CTTTTG', 'GATTTT', 'GCTTTT', 'GGTTTT', 'GTCTTT', 'GTGTTT', 'GTTATT',
               'GTTCTT', 'GTTGTT', 'GTTTAT', 'GTTTCT', 'GTTTGT', 'GTTTTA', 'GTTTTC', 'GTTTTG', 'TAATTT', 'TACTTT', 'TAGTTT', 'TATATT',
               'TATCTT', 'TATGTT', 'TATTAT', 'TATTCT', 'TATTGT', 'TATTTA', 'TATTTC', 'TATTTG', 'TCATTT', 'TCCTTT', 'TCGTTT', 'TCTATT',
               'TCTCTT', 'TCTGTT', 'TCTTAT', 'TCTTCT', 'TCTTGT', 'TCTTTA', 'TCTTTC', 'TCTTTG', 'TGATTT', 'TGCTTT', 'TGGTTT', 'TGTCTT',
               'TGTGTT', 'TGTTAT', 'TGTTCT', 'TGTTGT', 'TGTTTA', 'TGTTTC', 'TGTTTG', 'TTAATT', 'TTACTT', 'TTAGTT', 'TTATAT', 'TTATCT',
               'TTATGT', 'TTATTA', 'TTATTC', 'TTATTG', 'TTCATT', 'TTCCTT', 'TTCGTT', 'TTCTAT', 'TTCTCT', 'TTCTGT', 'TTCTTA', 'TTCTTC',
               'TTCTTG', 'TTGATT', 'TTGCTT', 'TTGGTT', 'TTGTCT', 'TTGTGT', 'TTGTTA', 'TTGTTC', 'TTGTTG', 'TTTAAT', 'TTTACT', 'TTTAGT',
               'TTTATA', 'TTTATC', 'TTTATG', 'TTTCAT', 'TTTCCT', 'TTTCGT', 'TTTCTA', 'TTTCTC', 'TTTCTG', 'TTTGAT', 'TTTGCT', 'TTTGGT',
               'TTTGTC', 'TTTGTG', 'TTTTAA', 'TTTTAC', 'TTTTAG', 'TTTTCA', 'TTTTCC', 'TTTTCG', 'TTTTGA', 'TTTTGC', 'TTTTGG')

a_rich_d0 <- c('AAAAAA')

a_rich_d1 <- c('AAAAAC', 'AAAAAG', 'AAAAAT', 'AAAACA', 'AAAAGA', 'AAAATA', 'AAACAA', 'AAAGAA', 'AAATAA', 'AACAAA', 'AAGAAA', 'AATAAA',
               'ACAAAA', 'AGAAAA', 'ATAAAA', 'CAAAAA', 'GAAAAA', 'TAAAAA')

a_rich_d2 <- c('AAAACC', 'AAAACG', 'AAAACT', 'AAAAGC', 'AAAAGG', 'AAAAGT', 'AAAATC', 'AAAATG', 'AAAATT', 'AAACAC', 'AAACAG', 'AAACAT',
               'AAACCA', 'AAACGA', 'AAACTA', 'AAAGAC', 'AAAGAG', 'AAAGAT', 'AAAGCA', 'AAAGGA', 'AAATAC', 'AAATAG', 'AAATAT', 'AAATCA',
               'AAATGA', 'AAATTA', 'AACAAC', 'AACAAG', 'AACAAT', 'AACACA', 'AACAGA', 'AACATA', 'AACCAA', 'AACGAA', 'AACTAA', 'AAGAAC',
               'AAGAAG', 'AAGAAT', 'AAGACA', 'AAGAGA', 'AAGATA', 'AAGCAA', 'AAGGAA', 'AATAAC', 'AATAAG', 'AATAAT', 'AATACA', 'AATAGA',
               'AATATA', 'AATCAA', 'AATGAA', 'AATTAA', 'ACAAAC', 'ACAAAG', 'ACAAAT', 'ACAACA', 'ACAAGA', 'ACAATA', 'ACACAA', 'ACAGAA',
               'ACATAA', 'ACCAAA', 'ACGAAA', 'ACTAAA', 'AGAAAC', 'AGAAAG', 'AGAAAT', 'AGAACA', 'AGAAGA', 'AGAATA', 'AGACAA', 'AGAGAA',
               'AGATAA', 'AGCAAA', 'AGGAAA', 'ATAAAC', 'ATAAAG', 'ATAAAT', 'ATAACA', 'ATAAGA', 'ATAATA', 'ATACAA', 'ATAGAA', 'ATATAA',
               'ATCAAA', 'ATGAAA', 'ATTAAA', 'CAAAAC', 'CAAAAG', 'CAAAAT', 'CAAACA', 'CAAAGA', 'CAAATA', 'CAACAA', 'CAAGAA', 'CAATAA',
               'CACAAA', 'CAGAAA', 'CATAAA', 'CCAAAA', 'CGAAAA', 'CTAAAA', 'GAAAAC', 'GAAAAG', 'GAAAAT', 'GAAACA', 'GAAAGA', 'GAAATA',
               'GAACAA', 'GAAGAA', 'GAATAA', 'GACAAA', 'GAGAAA', 'GATAAA', 'GCAAAA', 'GGAAAA', 'TAAAAC', 'TAAAAG', 'TAAAAT', 'TAAACA',
               'TAAAGA', 'TAAATA', 'TAACAA', 'TAAGAA', 'TAATAA', 'TACAAA', 'TAGAAA', 'TATAAA', 'TCAAAA', 'TGAAAA', 'TTAAAA')

a_rich_d3 <- c('AAACCC', 'AAACCG', 'AAACCT', 'AAACGC', 'AAACGG', 'AAACGT', 'AAACTC', 'AAACTG', 'AAACTT', 'AAAGCC', 'AAAGCG', 'AAAGCT',
               'AAAGGC', 'AAAGGG', 'AAAGGT', 'AAAGTC', 'AAAGTG', 'AAAGTT', 'AAATCC', 'AAATCG', 'AAATCT', 'AAATGC', 'AAATGG', 'AAATGT',
               'AAATTC', 'AAATTG', 'AAATTT', 'AACACC', 'AACACG', 'AACACT', 'AACAGC', 'AACAGG', 'AACAGT', 'AACATC', 'AACATG', 'AACATT',
               'AACCAC', 'AACCAG', 'AACCAT', 'AACCCA', 'AACCGA', 'AACCTA', 'AACGAC', 'AACGAG', 'AACGAT', 'AACGCA', 'AACGGA', 'AACTAC',
               'AACTAG', 'AACTAT', 'AACTCA', 'AACTGA', 'AACTTA', 'AAGACC', 'AAGACG', 'AAGACT', 'AAGAGC', 'AAGAGG', 'AAGAGT', 'AAGATC',
               'AAGATG', 'AAGATT', 'AAGCAC', 'AAGCAG', 'AAGCAT', 'AAGCCA', 'AAGCGA', 'AAGCTA', 'AAGGAC', 'AAGGAG', 'AAGGAT', 'AAGGCA',
               'AAGGGA', 'AAGTCA', 'AAGTGA', 'AAGTTA', 'AATACC', 'AATACG', 'AATACT', 'AATAGC', 'AATAGG', 'AATAGT', 'AATATC', 'AATATG',
               'AATATT', 'AATCAC', 'AATCAG', 'AATCAT', 'AATCCA', 'AATCGA', 'AATCTA', 'AATGAC', 'AATGAG', 'AATGAT', 'AATGCA', 'AATGGA',
               'AATTAC', 'AATTAG', 'AATTAT', 'AATTCA', 'AATTGA', 'AATTTA', 'ACAACC', 'ACAACG', 'ACAACT', 'ACAAGC', 'ACAAGG', 'ACAAGT',
               'ACAATC', 'ACAATG', 'ACAATT', 'ACACAC', 'ACACAG', 'ACACAT', 'ACACCA', 'ACACGA', 'ACACTA', 'ACAGAC', 'ACAGAG', 'ACAGAT',
               'ACAGCA', 'ACAGGA', 'ACATAC', 'ACATAG', 'ACATAT', 'ACATCA', 'ACATGA', 'ACATTA', 'ACCAAC', 'ACCAAG', 'ACCAAT', 'ACCACA',
               'ACCAGA', 'ACCATA', 'ACCCAA', 'ACCGAA', 'ACCTAA', 'ACGAAC', 'ACGAAG', 'ACGAAT', 'ACGACA', 'ACGAGA', 'ACGATA', 'ACGCAA',
               'ACGGAA', 'ACTAAC', 'ACTAAG', 'ACTAAT', 'ACTACA', 'ACTAGA', 'ACTATA', 'ACTCAA', 'ACTGAA', 'ACTTAA', 'AGAACC', 'AGAACG',
               'AGAACT', 'AGAAGC', 'AGAAGG', 'AGAAGT', 'AGAATC', 'AGAATG', 'AGAATT', 'AGACAC', 'AGACAG', 'AGACAT', 'AGACCA', 'AGACGA',
               'AGACTA', 'AGAGAC', 'AGAGAG', 'AGAGAT', 'AGAGCA', 'AGAGGA', 'AGATAC', 'AGATAG', 'AGATAT', 'AGATCA', 'AGATGA', 'AGATTA',
               'AGCAAC', 'AGCAAG', 'AGCAAT', 'AGCACA', 'AGCAGA', 'AGCATA', 'AGCCAA', 'AGCGAA', 'AGCTAA', 'AGGAAC', 'AGGAAG', 'AGGAAT',
               'AGGACA', 'AGGAGA', 'AGGATA', 'AGGCAA', 'AGGGAA', 'AGTCAA', 'AGTGAA', 'AGTTAA', 'ATAACC', 'ATAACG', 'ATAACT', 'ATAAGC',
               'ATAAGG', 'ATAAGT', 'ATAATC', 'ATAATG', 'ATAATT', 'ATACAC', 'ATACAG', 'ATACAT', 'ATACCA', 'ATACGA', 'ATACTA', 'ATAGAC',
               'ATAGAG', 'ATAGAT', 'ATAGCA', 'ATAGGA', 'ATATAC', 'ATATAG', 'ATATAT', 'ATATCA', 'ATATGA', 'ATATTA', 'ATCAAC', 'ATCAAG',
               'ATCAAT', 'ATCACA', 'ATCAGA', 'ATCATA', 'ATCCAA', 'ATCGAA', 'ATCTAA', 'ATGAAC', 'ATGAAG', 'ATGAAT', 'ATGACA', 'ATGAGA',
               'ATGATA', 'ATGCAA', 'ATGGAA', 'ATTAAC', 'ATTAAG', 'ATTAAT', 'ATTACA', 'ATTAGA', 'ATTATA', 'ATTCAA', 'ATTGAA', 'ATTTAA',
               'CAAACC', 'CAAACG', 'CAAACT', 'CAAAGC', 'CAAAGG', 'CAAAGT', 'CAAATC', 'CAAATG', 'CAAATT', 'CAACAC', 'CAACAG', 'CAACAT',
               'CAACCA', 'CAACGA', 'CAACTA', 'CAAGAC', 'CAAGAG', 'CAAGAT', 'CAAGCA', 'CAAGGA', 'CAATAC', 'CAATAG', 'CAATAT', 'CAATCA',
               'CAATGA', 'CAATTA', 'CACAAC', 'CACAAG', 'CACAAT', 'CACACA', 'CACAGA', 'CACATA', 'CACCAA', 'CACGAA', 'CACTAA', 'CAGAAC',
               'CAGAAG', 'CAGAAT', 'CAGACA', 'CAGAGA', 'CAGATA', 'CAGCAA', 'CAGGAA', 'CATAAC', 'CATAAG', 'CATAAT', 'CATACA', 'CATAGA',
               'CATATA', 'CATCAA', 'CATGAA', 'CATTAA', 'CCAAAC', 'CCAAAG', 'CCAAAT', 'CCAACA', 'CCAAGA', 'CCAATA', 'CCACAA', 'CCAGAA',
               'CCATAA', 'CCCAAA', 'CCGAAA', 'CCTAAA', 'CGAAAC', 'CGAAAG', 'CGAAAT', 'CGAACA', 'CGAAGA', 'CGAATA', 'CGACAA', 'CGAGAA',
               'CGATAA', 'CGCAAA', 'CGGAAA', 'CTAAAC', 'CTAAAG', 'CTAAAT', 'CTAACA', 'CTAAGA', 'CTAATA', 'CTACAA', 'CTAGAA', 'CTATAA',
               'CTCAAA', 'CTGAAA', 'CTTAAA', 'GAAACC', 'GAAACG', 'GAAACT', 'GAAAGC', 'GAAAGG', 'GAAAGT', 'GAAATC', 'GAAATG', 'GAAATT',
               'GAACAC', 'GAACAG', 'GAACAT', 'GAACCA', 'GAACGA', 'GAACTA', 'GAAGAC', 'GAAGAG', 'GAAGAT', 'GAAGCA', 'GAAGGA', 'GAATAC',
               'GAATAG', 'GAATAT', 'GAATCA', 'GAATGA', 'GAATTA', 'GACAAC', 'GACAAG', 'GACAAT', 'GACACA', 'GACAGA', 'GACATA', 'GACCAA',
               'GACGAA', 'GACTAA', 'GAGAAC', 'GAGAAG', 'GAGAAT', 'GAGACA', 'GAGAGA', 'GAGATA', 'GAGCAA', 'GAGGAA', 'GATAAC', 'GATAAG',
               'GATAAT', 'GATACA', 'GATAGA', 'GATATA', 'GATCAA', 'GATGAA', 'GATTAA', 'GCAAAC', 'GCAAAG', 'GCAAAT', 'GCAACA', 'GCAAGA',
               'GCAATA', 'GCACAA', 'GCAGAA', 'GCATAA', 'GCCAAA', 'GCGAAA', 'GCTAAA', 'GGAAAC', 'GGAAAG', 'GGAAAT', 'GGAACA', 'GGAAGA',
               'GGAATA', 'GGACAA', 'GGAGAA', 'GGATAA', 'GGCAAA', 'GGGAAA', 'GTCAAA', 'GTGAAA', 'GTTAAA', 'TAAACC', 'TAAACG', 'TAAACT',
               'TAAAGC', 'TAAAGG', 'TAAAGT', 'TAAATC', 'TAAATG', 'TAAATT', 'TAACAC', 'TAACAG', 'TAACAT', 'TAACCA', 'TAACGA', 'TAACTA',
               'TAAGAC', 'TAAGAG', 'TAAGAT', 'TAAGCA', 'TAAGGA', 'TAATAC', 'TAATAG', 'TAATAT', 'TAATCA', 'TAATGA', 'TAATTA', 'TACAAC',
               'TACAAG', 'TACAAT', 'TACACA', 'TACAGA', 'TACATA', 'TACCAA', 'TACGAA', 'TACTAA', 'TAGAAC', 'TAGAAG', 'TAGAAT', 'TAGACA',
               'TAGAGA', 'TAGATA', 'TAGCAA', 'TAGGAA', 'TATAAC', 'TATAAG', 'TATAAT', 'TATACA', 'TATAGA', 'TATATA', 'TATCAA', 'TATGAA',
               'TATTAA', 'TCAAAC', 'TCAAAG', 'TCAAAT', 'TCAACA', 'TCAAGA', 'TCAATA', 'TCACAA', 'TCAGAA', 'TCATAA', 'TCCAAA', 'TCGAAA',
               'TCTAAA', 'TGAAAC', 'TGAAAG', 'TGAAAT', 'TGAACA', 'TGAAGA', 'TGAATA', 'TGACAA', 'TGAGAA', 'TGATAA', 'TGCAAA', 'TGGAAA',
               'TTAAAC', 'TTAAAG', 'TTAAAT', 'TTAACA', 'TTAAGA', 'TTAATA', 'TTACAA', 'TTAGAA', 'TTATAA', 'TTCAAA', 'TTGAAA', 'TTTAAA')

tag_d0 <- c('AAATAG', 'AACTAG', 'AAGTAG', 'AATAGA', 'AATAGC', 'AATAGG', 'AATAGT', 'AATTAG', 'ACATAG', 'ACCTAG', 'ACGTAG', 'ACTAGA',
            'ACTAGC', 'ACTAGG', 'ACTAGT', 'ACTTAG', 'AGATAG', 'AGCTAG', 'AGGTAG', 'AGTAGA', 'AGTAGC', 'AGTAGG', 'AGTAGT', 'AGTTAG',
            'ATAGAA', 'ATAGAC', 'ATAGAG', 'ATAGAT', 'ATAGCA', 'ATAGCC', 'ATAGCG', 'ATAGCT', 'ATAGGA', 'ATAGGC', 'ATAGGG', 'ATAGGT',
            'ATAGTA', 'ATAGTC', 'ATAGTG', 'ATAGTT', 'ATATAG', 'ATCTAG', 'ATGTAG', 'ATTAGA', 'ATTAGC', 'ATTAGG', 'ATTAGT', 'ATTTAG',
            'CAATAG', 'CACTAG', 'CAGTAG', 'CATAGA', 'CATAGC', 'CATAGG', 'CATAGT', 'CATTAG', 'CCATAG', 'CCCTAG', 'CCGTAG', 'CCTAGA',
            'CCTAGC', 'CCTAGG', 'CCTAGT', 'CCTTAG', 'CGATAG', 'CGCTAG', 'CGGTAG', 'CGTAGA', 'CGTAGC', 'CGTAGG', 'CGTAGT', 'CGTTAG',
            'CTAGAA', 'CTAGAC', 'CTAGAG', 'CTAGAT', 'CTAGCA', 'CTAGCC', 'CTAGCG', 'CTAGCT', 'CTAGGA', 'CTAGGC', 'CTAGGG', 'CTAGGT',
            'CTAGTA', 'CTAGTC', 'CTAGTG', 'CTAGTT', 'CTATAG', 'CTCTAG', 'CTGTAG', 'CTTAGA', 'CTTAGC', 'CTTAGG', 'CTTAGT', 'CTTTAG',
            'GAATAG', 'GACTAG', 'GAGTAG', 'GATAGA', 'GATAGC', 'GATAGG', 'GATAGT', 'GATTAG', 'GCATAG', 'GCCTAG', 'GCGTAG', 'GCTAGA',
            'GCTAGC', 'GCTAGG', 'GCTAGT', 'GCTTAG', 'GGATAG', 'GGCTAG', 'GGGTAG', 'GGTAGA', 'GGTAGC', 'GGTAGG', 'GGTAGT', 'GGTTAG',
            'GTAGAA', 'GTAGAC', 'GTAGAG', 'GTAGAT', 'GTAGCA', 'GTAGCC', 'GTAGCG', 'GTAGCT', 'GTAGGA', 'GTAGGC', 'GTAGGG', 'GTAGGT',
            'GTAGTA', 'GTAGTC', 'GTAGTG', 'GTAGTT', 'GTATAG', 'GTCTAG', 'GTGTAG', 'GTTAGA', 'GTTAGC', 'GTTAGG', 'GTTAGT', 'GTTTAG',
            'TAATAG', 'TACTAG', 'TAGAAA', 'TAGAAC', 'TAGAAG', 'TAGAAT', 'TAGACA', 'TAGACC', 'TAGACG', 'TAGACT', 'TAGAGA', 'TAGAGC',
            'TAGAGG', 'TAGAGT', 'TAGATA', 'TAGATC', 'TAGATG', 'TAGATT', 'TAGCAA', 'TAGCAC', 'TAGCAG', 'TAGCAT', 'TAGCCA', 'TAGCCC',
            'TAGCCG', 'TAGCCT', 'TAGCGA', 'TAGCGC', 'TAGCGG', 'TAGCGT', 'TAGCTA', 'TAGCTC', 'TAGCTG', 'TAGCTT', 'TAGGAA', 'TAGGAC',
            'TAGGAG', 'TAGGAT', 'TAGGCA', 'TAGGCC', 'TAGGCG', 'TAGGCT', 'TAGGGA', 'TAGGGC', 'TAGGGG', 'TAGGGT', 'TAGGTA', 'TAGGTC',
            'TAGGTG', 'TAGGTT', 'TAGTAA', 'TAGTAC', 'TAGTAG', 'TAGTAT', 'TAGTCA', 'TAGTCC', 'TAGTCG', 'TAGTCT', 'TAGTGA', 'TAGTGC',
            'TAGTGG', 'TAGTGT', 'TAGTTA', 'TAGTTC', 'TAGTTG', 'TAGTTT', 'TATAGA', 'TATAGC', 'TATAGG', 'TATAGT', 'TATTAG', 'TCATAG',
            'TCCTAG', 'TCGTAG', 'TCTAGA', 'TCTAGC', 'TCTAGG', 'TCTAGT', 'TCTTAG', 'TGATAG', 'TGCTAG', 'TGGTAG', 'TGTAGA', 'TGTAGC',
            'TGTAGG', 'TGTAGT', 'TGTTAG', 'TTAGAA', 'TTAGAC', 'TTAGAG', 'TTAGAT', 'TTAGCA', 'TTAGCC', 'TTAGCG', 'TTAGCT', 'TTAGGA',
            'TTAGGC', 'TTAGGG', 'TTAGGT', 'TTAGTA', 'TTAGTC', 'TTAGTG', 'TTAGTT', 'TTATAG', 'TTCTAG', 'TTGTAG', 'TTTAGA', 'TTTAGC',
            'TTTAGG', 'TTTAGT', 'TTTTAG')

gta_d0 <- c('AAAGTA', 'AACGTA', 'AAGGTA', 'AAGTAA', 'AAGTAC', 'AAGTAG', 'AAGTAT', 'AATGTA', 'ACAGTA', 'ACCGTA', 'ACGGTA', 'ACGTAA',
            'ACGTAC', 'ACGTAG', 'ACGTAT', 'ACTGTA', 'AGAGTA', 'AGCGTA', 'AGGGTA', 'AGGTAA', 'AGGTAC', 'AGGTAG', 'AGGTAT', 'AGTAAA',
            'AGTAAC', 'AGTAAG', 'AGTAAT', 'AGTACA', 'AGTACC', 'AGTACG', 'AGTACT', 'AGTAGA', 'AGTAGC', 'AGTAGG', 'AGTAGT', 'AGTATA',
            'AGTATC', 'AGTATG', 'AGTATT', 'AGTGTA', 'ATAGTA', 'ATCGTA', 'ATGGTA', 'ATGTAA', 'ATGTAC', 'ATGTAG', 'ATGTAT', 'ATTGTA',
            'CAAGTA', 'CACGTA', 'CAGGTA', 'CAGTAA', 'CAGTAC', 'CAGTAG', 'CAGTAT', 'CATGTA', 'CCAGTA', 'CCCGTA', 'CCGGTA', 'CCGTAA',
            'CCGTAC', 'CCGTAG', 'CCGTAT', 'CCTGTA', 'CGAGTA', 'CGCGTA', 'CGGGTA', 'CGGTAA', 'CGGTAC', 'CGGTAG', 'CGGTAT', 'CGTAAA',
            'CGTAAC', 'CGTAAG', 'CGTAAT', 'CGTACA', 'CGTACC', 'CGTACG', 'CGTACT', 'CGTAGA', 'CGTAGC', 'CGTAGG', 'CGTAGT', 'CGTATA',
            'CGTATC', 'CGTATG', 'CGTATT', 'CGTGTA', 'CTAGTA', 'CTCGTA', 'CTGGTA', 'CTGTAA', 'CTGTAC', 'CTGTAG', 'CTGTAT', 'CTTGTA',
            'GAAGTA', 'GACGTA', 'GAGGTA', 'GAGTAA', 'GAGTAC', 'GAGTAG', 'GAGTAT', 'GATGTA', 'GCAGTA', 'GCCGTA', 'GCGGTA', 'GCGTAA',
            'GCGTAC', 'GCGTAG', 'GCGTAT', 'GCTGTA', 'GGAGTA', 'GGCGTA', 'GGGGTA', 'GGGTAA', 'GGGTAC', 'GGGTAG', 'GGGTAT', 'GGTAAA',
            'GGTAAC', 'GGTAAG', 'GGTAAT', 'GGTACA', 'GGTACC', 'GGTACG', 'GGTACT', 'GGTAGA', 'GGTAGC', 'GGTAGG', 'GGTAGT', 'GGTATA',
            'GGTATC', 'GGTATG', 'GGTATT', 'GGTGTA', 'GTAAAA', 'GTAAAC', 'GTAAAG', 'GTAAAT', 'GTAACA', 'GTAACC', 'GTAACG', 'GTAACT',
            'GTAAGA', 'GTAAGC', 'GTAAGG', 'GTAAGT', 'GTAATA', 'GTAATC', 'GTAATG', 'GTAATT', 'GTACAA', 'GTACAC', 'GTACAG', 'GTACAT',
            'GTACCA', 'GTACCC', 'GTACCG', 'GTACCT', 'GTACGA', 'GTACGC', 'GTACGG', 'GTACGT', 'GTACTA', 'GTACTC', 'GTACTG', 'GTACTT',
            'GTAGAA', 'GTAGAC', 'GTAGAG', 'GTAGAT', 'GTAGCA', 'GTAGCC', 'GTAGCG', 'GTAGCT', 'GTAGGA', 'GTAGGC', 'GTAGGG', 'GTAGGT',
            'GTAGTA', 'GTAGTC', 'GTAGTG', 'GTAGTT', 'GTATAA', 'GTATAC', 'GTATAG', 'GTATAT', 'GTATCA', 'GTATCC', 'GTATCG', 'GTATCT',
            'GTATGA', 'GTATGC', 'GTATGG', 'GTATGT', 'GTATTA', 'GTATTC', 'GTATTG', 'GTATTT', 'GTCGTA', 'GTGGTA', 'GTGTAA', 'GTGTAC',
            'GTGTAG', 'GTGTAT', 'GTTGTA', 'TAAGTA', 'TACGTA', 'TAGGTA', 'TAGTAA', 'TAGTAC', 'TAGTAG', 'TAGTAT', 'TATGTA', 'TCAGTA',
            'TCCGTA', 'TCGGTA', 'TCGTAA', 'TCGTAC', 'TCGTAG', 'TCGTAT', 'TCTGTA', 'TGAGTA', 'TGCGTA', 'TGGGTA', 'TGGTAA', 'TGGTAC',
            'TGGTAG', 'TGGTAT', 'TGTAAA', 'TGTAAC', 'TGTAAG', 'TGTAAT', 'TGTACA', 'TGTACC', 'TGTACG', 'TGTACT', 'TGTAGA', 'TGTAGC',
            'TGTAGG', 'TGTAGT', 'TGTATA', 'TGTATC', 'TGTATG', 'TGTATT', 'TGTGTA', 'TTAGTA', 'TTCGTA', 'TTGGTA', 'TTGTAA', 'TTGTAC',
            'TTGTAG', 'TTGTAT', 'TTTGTA')

gta_tag_d0 <- c('AAGTAG', 'ACGTAG', 'AGGTAG', 'AGTAGA', 'AGTAGC', 'AGTAGG', 'AGTAGT', 'ATAGTA', 'ATGTAG', 'CAGTAG', 'CCGTAG', 'CGGTAG',
                'CGTAGA', 'CGTAGC', 'CGTAGG', 'CGTAGT', 'CTAGTA', 'CTGTAG', 'GAGTAG', 'GCGTAG', 'GGGTAG', 'GGTAGA', 'GGTAGC', 'GGTAGG',
                'GGTAGT', 'GTAGAA', 'GTAGAC', 'GTAGAG', 'GTAGAT', 'GTAGCA', 'GTAGCC', 'GTAGCG', 'GTAGCT', 'GTAGGA', 'GTAGGC', 'GTAGGG',
                'GTAGGT', 'GTAGTA', 'GTAGTC', 'GTAGTG', 'GTAGTT', 'GTATAG', 'GTGTAG', 'TAGGTA', 'TAGTAA', 'TAGTAC', 'TAGTAG', 'TAGTAT',
                'TCGTAG', 'TGGTAG', 'TGTAGA', 'TGTAGC', 'TGTAGG', 'TGTAGT', 'TTAGTA', 'TTGTAG')

dt_pom <- read_delim("./analysis/polyaclassifier_motif_importance_heatmaps/hmdf_sum.schizosaccharomyces_pombe.polyaclassifier_bagging3_kmers-6.txt", delim = "\t")

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
           show_rownames = FALSE, show_colnames = FALSE,
           annotation_col = ann['motifFamily'], 
           annotation_row = ann[c('motifFamily','t_count','g_count','c_count','a_count')],
           annotation_colors = ann_colors,
           cluster_cols=as.hclust(col_dend),
           cluster_rows=as.hclust(row_dend),
           )

ggsave("./analysis/polyaclassifier_motif_importance_heatmaps/spom.6mers.hm_correlation.svg", plot = g2, device = svg())

myMin <- -100
myMax <- +100

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

ggsave("./analysis/polyaclassifier_motif_importance_heatmaps/spom.6mers.hm_importance.reordered.svg", plot = g3, device = svg())

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


