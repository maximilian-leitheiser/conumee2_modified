chromInfo = read.csv("/media/Leia/shared/Scripts/Methylation/max_scripts/manifest_files/UCSC_tables/chromInfo_hg38.tsv", sep = "\t")
gap = read.csv("/media/Leia/shared/Scripts/Methylation/max_scripts/manifest_files/UCSC_tables/gap_hg38.tsv", sep = "\t")
cytoBand = read.csv("/media/Leia/shared/Scripts/Methylation/max_scripts/manifest_files/UCSC_tables/cyto_band_hg38.tsv", sep = "\t")
centromeres = read.csv("/media/Leia/shared/Scripts/Methylation/max_scripts/manifest_files/UCSC_tables/centromere_hg38.tsv", sep = "\t")


tbl_UCSC_hg38 = list(chromInfo = chromInfo, gap = gap, cytoBand = cytoBand, centromeres = centromeres)
tbl_UCSC_hg38 = lapply(tbl_UCSC_hg38, function(tbl){
  names(tbl)[[1]] = str_sub(names(tbl)[[1]], start = 3)
  tbl
})
save(list = "tbl_UCSC_hg38", file = "/home/leitheim/CNV_analysis/conumee2/conumee2_modified/data/tbl_UCSC_hg38.rda")
