Cosmic_CancerGeneCensus_v99_GRCh37 = read.csv("/mnt/nas/private/leitheim/CNV_analysis/CGC_annotation/Cosmic_CancerGeneCensus_v99_GRCh37.tsv", sep = "\t")
Cosmic_CancerGeneCensus_v99_GRCh38 = read.csv("/mnt/nas/private/leitheim/CNV_analysis/CGC_annotation/Cosmic_CancerGeneCensus_v99_GRCh38.tsv", sep = "\t")


save(list = c("Cosmic_CancerGeneCensus_v99_GRCh37"), 
     file = "/mnt/nas/private/leitheim/conumee_current/conumee2_modified/data/Cosmic_CancerGeneCensus_v99_GRCh37.rda")
save(list = c("Cosmic_CancerGeneCensus_v99_GRCh38"), 
     file = "/mnt/nas/private/leitheim/conumee_current/conumee2_modified/data/Cosmic_CancerGeneCensus_v99_GRCh38.rda")
