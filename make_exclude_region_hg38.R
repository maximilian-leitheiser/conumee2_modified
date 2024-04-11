exclude_region_hg38 = GenomicRanges::GRanges(seqnames = c("chr6", "chr19"), 
                       ranges = IRanges(start = c(32432224, 42495849), 
                                        end = c(32832223, 43395848)), 
                       name = c("HLA", "CEACAM-PSG"))
  

save(list = c("exclude_region_hg38"), 
     file = "/home/leitheim/CNV_analysis/conumee2/conumee2_modified/data/exclude_region_hg38.rda")


  