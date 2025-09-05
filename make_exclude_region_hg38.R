exclude_regions_hg38 = GenomicRanges::GRanges(seqnames = c("chr6", "chr19", "chr22"), 
                       ranges = IRanges(start = c(32432224, 42495849, 23997807), 
                                        end = c(32832223, 43395848, 24400000)), 
                       name = c("HLA", "CEACAM-PSG", "GSTT1"))

save(list = c("exclude_regions_hg38"), 
     file = "/mnt/ssd/private/leitheim/leitheim_nas/conumee_current/conumee2_modified/data/exclude_regions_hg38.rda")


  