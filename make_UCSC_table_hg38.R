chromInfo = read.csv("/media/Leia/shared/Scripts/Methylation/max_scripts/manifest_files/UCSC_tables/chromInfo.tsv", sep = "\t")
gap = read.csv("/media/Leia/shared/Scripts/Methylation/max_scripts/manifest_files/UCSC_tables/gap.tsv", sep = "\t")
cytoBand = read.csv("/media/Leia/shared/Scripts/Methylation/max_scripts/manifest_files/UCSC_tables/cyto_band.tsv", sep = "\t")

tbl_UCSC_hg38 = list(chromInfo = chromInfo, gap = gap, cytoBand = cytoBand)

save(list = "tbl_UCSC_hg38", file = "/home/leitheim/CNV_analysis/conumee2/conumee2_modified/data/tbl_UCSC_hg38.rda")