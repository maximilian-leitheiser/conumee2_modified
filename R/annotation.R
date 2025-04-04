#### aux functions ####
# TODO: deal more clearly with genome information, especially XY chromosome in focal detection
# TODO: make .create_focal_ranges use genome info from anno_object instead of 'Seqinfo(hg19)'
.create_focal_ranges = function(genome_info){
  # load Cancer Genome Census information
  CGC_df = Cosmic_CancerGeneCensus_v99_GRCh37
  rownames(CGC_df) = CGC_df$GENE_SYMBOL
  CGC_df = CGC_df[!(is.na(CGC_df$GENOME_START) | is.na(CGC_df$GENOME_STOP)), ]
  
  # adapt to GRCh37 sequence names to hg19 (as used in rest of package)
  CGC_df$CHROMOSOME = paste0("chr", CGC_df$CHROMOSOME)
  
  # remove genes that are not on chromosomes currently included (relevant in particular for chrX, chrY)
  CGC_df = CGC_df[CGC_df$CHROMOSOME %in% genome_info$chr, ]
  
  # create GRanges object
  non_meta_colnames_CGC = c("CHROMOSOME", "GENOME_START", "GENOME_STOP")
  CGC_ranges = GRanges(seqnames = CGC_df$CHROMOSOME,
                       ranges = IRanges(start = CGC_df$GENOME_START,
                                        end = CGC_df$GENOME_STOP,
                                        names = CGC_df$GENE_SYMBOL),
                       strand = "*",
                       seqinfo = Seqinfo(seqnames = genome_info$chr, 
                                         seqlengths = genome_info$size),
                       mcols = CGC_df[, -match(non_meta_colnames_CGC, colnames(CGC_df))])
  
  colnames(mcols(CGC_ranges)) = str_remove(colnames(mcols(CGC_ranges)), pattern = "mcols.")
  
  # return
  return(CGC_ranges)
}




.create_detail_ranges = function(genome_info, detail_regions){
  
  ## empty detail ranges
  if(is.null(detail_regions)){
    detail_ranges = GRanges(seqinfo = Seqinfo(genome_info$chr, genome_info$size))
    return(detail_ranges)
  }
  
  ## filled detail ranges
  if (class(detail_regions) == "GRanges") {
    detail_ranges <- GRanges(as.vector(seqnames(detail_regions)),
                             ranges(detail_regions), seqinfo = Seqinfo(genome_info$chr,
                                                                       genome_info$size))
    if (any(grepl("name", names(values(detail_regions))))) {
      values(detail_ranges)$name <- values(detail_regions)[[grep("name",
                                                                 names(values(detail_regions)))[1]]]
    }
    if (any(grepl("IRanges", sapply(values(detail_regions), class)))) {
      values(detail_ranges)$thick <- values(detail_regions)[[grep("IRanges",
                                                                  sapply(values(detail_regions), class))[1]]]
    }
    detail_ranges <- sort(detail_ranges)
  } else {
    detail_ranges <- sort(rtracklayer::import(detail_regions, seqinfo = Seqinfo(genome_info$chr,
                                                                                genome_info$size)))
  }
  if (!is.element("name", names(values(detail_ranges)))) {
    stop("detailed region bed file must contain name column.")
  }
  if (!all(table(values(detail_ranges)$name) == 1)) {
    stop("detailed region names must be unique.")
  }
  
  if (!is.element("thick", names(values(detail_ranges)))) {
    values(detail_ranges)$thick <- resize(ranges(detail_ranges), fix = "center",
                                          1e+06)
  }
  
  return(detail_ranges)
}



.create_anno_human = function(anno_object, array_type, chrXY, detail_regions, exclude_regions, 
                              bin_minprobes = 15, bin_minsize = 50000, bin_maxsize = 5e+06, 
                              custom_probe_set = NULL, ref_gene = c("hg19", "hg38"), cg_probes_only = FALSE){
  
  # get genome info (chromosome dataframe and gap ranges )
  tmp_list = .get_genome_info_and_gap(ref_gene = ref_gene, chrXY = chrXY)
  anno_object@genome = tmp_list$genome
  anno_object@gap = tmp_list$gap
  
  
  if(ref_gene == "hg19"){
    # TODO: this is simply the old code
    # TODO: EPICv2 is in facto NOT hg19. THis should be corrected/harmonized
    
    probes450k <- probesEPIC <- probesEPICv2 <- GRanges()
    if (is.element(array_type, c("450k", "overlap"))) {
      message("getting 450k annotations")
      data("UCSC_RefGene_Name_450k")
      probes450k <- minfi::getLocations(IlluminaHumanMethylation450kanno.ilmn12.hg19::IlluminaHumanMethylation450kanno.ilmn12.hg19)
      probes450k$genes <- UCSC_RefGene_Name_450k
      probes450k <- sort(probes450k)
    }
    if (is.element(array_type, c("EPIC", "overlap"))) {
      message("getting EPIC annotations")
      data("UCSC_RefGene_Name_EPIC")
      probesEPIC <- minfi::getLocations(IlluminaHumanMethylationEPICanno.ilm10b4.hg19::IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
      probesEPIC$genes <- UCSC_RefGene_Name_EPIC
      probesEPIC <- sort(probesEPIC)
    }
    
    if (is.element(array_type, "EPICv2")) {
      message("getting EPICv2 annotations")
      data("EPICv2_hg19_probes")
      probesEPICv2 <- sort(EPICv2_hg19_probes)
    }
    
  } else if(ref_gene == "hg38") {
    # TODO: 'genes' are still missing
    # TODO: sesameDataGet should maybe not be used in this package that way
    
    probes450k <- probesEPIC <- probesEPICv2 <- GRanges()
    
    if (is.element(array_type, c("450k", "overlap"))) {
      message("getting 450k annotations")
      probes450k <- sesameData::sesameDataGet("HM450.address")$hg38
      probes450k <- sort(probes450k, ignore.strand = TRUE)
      probes450k$genes = "NA"
    }
    if (is.element(array_type, c("EPIC", "overlap"))) {
      message("getting EPIC annotations")
      probesEPIC <- sesameData::sesameDataGet("EPIC.address")$hg38
      probesEPIC <- sort(probesEPIC, ignore.strand = TRUE)
      probesEPIC$genes = "NA"
    }
    
    if (is.element(array_type, "EPICv2")) {
      message("getting EPICv2 annotations")
      probesEPICv2 <- sesameData::sesameDataGet("EPICv2.address")$hg38
      probesEPICv2 <- sort(probesEPICv2, ignore.strand = TRUE)
      probesEPICv2$genes = "NA"
      
      # TODO: collapsing is done inefficiently
      is_ctl_probe = grepl(pattern = "ctl", x = names(probesEPICv2))
      prefix_vec = names(probesEPICv2)
      prefix_vec[!is_ctl_probe]  = sapply(str_split(names(probesEPICv2[!is_ctl_probe]), "_"), head, 1)
      
      probes_EPICv2_coll = probesEPICv2[!duplicated(prefix_vec)]
      is_ctl_probe = grepl(pattern = "ctl", x = names(probes_EPICv2_coll))
      names(probes_EPICv2_coll)[!is_ctl_probe] = sapply(str_split(names(probes_EPICv2_coll)[!is_ctl_probe], "_"), head, 1)
      
      probesEPICv2 <- sort(probesEPICv2, ignore.strand = TRUE)
    }
    
  } else {
    stop("reference genome not recognized")
  }
  
  
  
  if (array_type == "overlap") {
    probes <- sort(subsetByOverlaps(probes450k, probesEPIC))
  } else {
    probes <- sort(c(probes450k, probesEPIC, probesEPICv2))
    probes = probes[!duplicated(names(probes))]
  }
  
  if(!is.null(custom_probe_set)){
    # check for unmatched probes
    unmatched_probe_vec = setdiff(custom_probe_set, names(probes))
    if(length(unmatched_probe_vec) > 0){
      stop(length(unmatched_probe_vec), " probes in custom_probe_set are not present in the probe set of the chosen chip type.")
    }
    
    # reduce chip probe set to custom probe set
    probes = probes[custom_probe_set]
    probes = sort(unique(probes))
  }
  
  
  # CpG probes only
  if(cg_probes_only){
    allowed_prefix_set = c("cg")
  } else {
    allowed_prefix_set = c("cg", "ch")
  }
  
  ao_probes <- probes[substr(names(probes), 1, 2) %in% allowed_prefix_set & is.element(as.vector(seqnames(probes)), anno_object@genome$chr)]
  anno_object@probes <- sort(GRanges(as.vector(seqnames(ao_probes)), ranges(ao_probes),
                                seqinfo = Seqinfo(anno_object@genome$chr, anno_object@genome$size)))
  anno_object@probes$genes <- ao_probes$genes
  
  message(" - ", length(anno_object@probes), " probes used")
  
  ## exclude regions
  if (!is.null(exclude_regions)) {
    message("importing regions to exclude from analysis")
    if (class(exclude_regions) == "GRanges") {
      anno_object@exclude <- GRanges(as.vector(seqnames(exclude_regions)),
                                ranges(exclude_regions), seqinfo = Seqinfo(anno_object@genome$chr,
                                                                           anno_object@genome$size))
      values(anno_object@exclude) <- values(exclude_regions)
      anno_object@exclude <- sort(anno_object@exclude)
    } else {
      anno_object@exclude <- sort(rtracklayer::import(exclude_regions,
                                                 seqinfo = Seqinfo(anno_object@genome$chr, anno_object@genome$size)))
    }
  } else {
    anno_object@exclude <- GRanges(seqinfo = Seqinfo(anno_object@genome$chr,
                                                anno_object@genome$size))
  }
  
  ## detail region
  # TODO: does not work yet for hg38
  message("importing regions for detailed analysis")
  anno_object@detail = .create_detail_ranges(genome_info = anno_object@genome, 
                                             detail_regions = detail_regions)
  
  ## cancer genes
  # TODO: does not work yet for hg38
  message("importing cancer-related genes for focal analysis")
  anno_object@cancer_genes = .create_focal_ranges(genome_info = anno_object@genome)
  
  ## creating bins
  message("creating bins")
  anno.tile <- CNV.create_bins(hg19.anno = anno_object@genome, bin_minsize = bin_minsize,
                               hg19.gap = anno_object@gap, hg19.exclude = anno_object@exclude)
  message(" - ", length(anno.tile), " bins created")
  
  message("merging bins")
  anno_object@bins <- CNV.merge_bins(hg19.anno = anno_object@genome, hg19.tile = anno.tile,
                                bin_minprobes = bin_minprobes, hg19.probes = anno_object@probes, bin_maxsize = bin_maxsize)
  message(" - ", length(anno_object@bins), " bins remaining")
  
  
  # TODO: does not work yet for hg38 (because genes annotation is not right yet)
  message("getting the gene annotations for each bin")
  o <- findOverlaps(anno_object@probes, anno_object@bins)
  bin_genes <- sapply(lapply(sapply(split(anno_object@probes$genes[queryHits(o)], names(anno_object@bins)[subjectHits(o)]),
                                    function(x) na.omit(unlist(strsplit(x,split = ";")))), unique), paste, collapse = ";")
  
  
  
  anno_object@bins$genes <- bin_genes
  
  return(anno_object)
}



.create_anno_mouse = function(anno_object, chrXY, exclude_regions, detail_regions, array_type,
                              bin_minprobes = 15, bin_minsize = 50000, bin_maxsize = 5e+06){
  
  data("mouse_annotation")
  
  if (chrXY) {
    anno_object@genome <- data.frame(chr = paste("chr", c(1:19, "X", "Y"),
                                            sep = ""), stringsAsFactors = FALSE)
  } else {
    anno_object@genome <- data.frame(chr = paste("chr", 1:19, sep = ""),
                                stringsAsFactors = FALSE)
  }
  
  rownames(anno_object@genome) <- anno_object@genome$chr
  
  message("using mm10 genome annotations from UCSC")
  
  anno_object@genome$size <- mouse_annotation[[1]]$size[1:nrow(anno_object@genome)]
  
  tbl.gap <- mouse_annotation[[2]][is.element(mouse_annotation[[2]]$chrom, anno_object@genome$chr),]
  
  anno_object@gap <- sort(GRanges(as.vector(tbl.gap$chrom), IRanges(tbl.gap$chromStart + 1,
                                                               tbl.gap$chromEnd), seqinfo = Seqinfo(anno_object@genome$chr, anno_object@genome$size)))
  
  
  mouse_probes <- GRanges(as.vector(paste("chr",mouse_annotation[[3]]$CHR, sep = "")),
                          IRanges(start = mouse_annotation[[3]]$MAPINFO,
                                  end = mouse_annotation[[3]]$MAPINFO), seqinfo = Seqinfo(anno_object@genome$chr, anno_object@genome$size, genome = "mm10"))
  
  names(mouse_probes) <- mouse_annotation[[3]]$Name
  mouse_probes$genes <- mouse_annotation[[3]]$genes
  
  
  # CpG probes only
  mouse_probes <- mouse_probes[substr(names(mouse_probes),1, 2) == "cg" & is.element(as.vector(seqnames(mouse_probes)), anno_object@genome$chr)]
  seqlevels(mouse_probes)<- c(paste("chr", 1:19, sep = ""), "chrY", "chrX")
  mouse_probes <- sort(mouse_probes)
  
  anno_object@probes <- mouse_probes
  
  message(" - ", length(anno_object@probes), " probes used")
  
  if (!is.null(exclude_regions)) {
    message("importing regions to exclude from analysis")
    if (class(exclude_regions) == "GRanges") {
      anno_object@exclude <- GRanges(as.vector(seqnames(exclude_regions)),
                                ranges(exclude_regions), seqinfo = Seqinfo(anno_object@genome$chr,
                                                                           anno_object@genome$size))
      values(anno_object@exclude) <- values(exclude_regions)
      anno_object@exclude <- sort(anno_object@exclude)
    } else {
      anno_object@exclude <- sort(rtracklayer::import(exclude_regions,
                                                 seqinfo = Seqinfo(anno_object@genome$chr, anno_object@genome$size)))
    }
  } else {
    anno_object@exclude <- GRanges(seqinfo = Seqinfo(anno_object@genome$chr,
                                                anno_object@genome$size))
  }
  
  if (!is.null(detail_regions)) {
    message("importing regions for detailed analysis")
    if (class(detail_regions) == "GRanges") {
      anno_object@detail <- GRanges(as.vector(seqnames(detail_regions)),
                               ranges(detail_regions), seqinfo = Seqinfo(anno_object@genome$chr,
                                                                         anno_object@genome$size))
      if (any(grepl("name", names(values(detail_regions))))) {
        values(anno_object@detail)$name <- values(detail_regions)[[grep("name",
                                                                   names(values(detail_regions)))[1]]]
      }
      if (any(grepl("IRanges", sapply(values(detail_regions), class)))) {
        values(anno_object@detail)$thick <- values(detail_regions)[[grep("IRanges",
                                                                    sapply(values(detail_regions), class))[1]]]
      }
      anno_object@detail <- sort(anno_object@detail)
    } else {
      anno_object@detail <- sort(rtracklayer::import(detail_regions, seqinfo = Seqinfo(anno_object@genome$chr,
                                                                                  anno_object@genome$size)))
    }
    if (!is.element("name", names(values(anno_object@detail)))) {
      stop("detailed region bed file must contain name column.")
    }
    if (!all(table(values(anno_object@detail)$name) == 1)) {
      stop("detailed region names must be unique.")
    }
  } else {
    anno_object@detail <- GRanges(seqinfo = Seqinfo(anno_object@genome$chr, anno_object@genome$size))
  }
  if (!is.element("thick", names(values(anno_object@detail)))) {
    values(anno_object@detail)$thick <- resize(ranges(anno_object@detail), fix = "center",
                                          1e+06)
  }
  
  message("creating bins")
  anno.tile <- CNV.create_bins(hg19.anno = anno_object@genome, bin_minsize = bin_minsize,
                               hg19.gap = anno_object@gap, hg19.exclude = anno_object@exclude)
  message(" - ", length(anno.tile), " bins created")
  
  message("merging bins")
  anno_object@bins <- CNV.merge_bins(hg19.anno = anno_object@genome, hg19.tile = anno.tile,
                                bin_minprobes = bin_minprobes, hg19.probes = anno_object@probes, bin_maxsize = bin_maxsize)
  message(" - ", length(anno_object@bins), " bins remaining")
  
  message("getting the gene annotations for each bin")
  
  o <- findOverlaps(anno_object@probes, anno_object@bins)
  bin_genes <- sapply(lapply(sapply(split(anno_object@probes$genes[queryHits(o)], names(anno_object@bins)[subjectHits(o)]),
                                    function(x) na.omit(unlist(strsplit(x,split = ";")))), unique), paste, collapse = ";")
  
  
  anno_object@bins$genes <- bin_genes
  
  return(anno_object)
}



.get_genome_info_and_gap = function(ref_gene = c("hg19", "hg38"), chrXY = FALSE){

  # parse arguments
  ref_gene = match.arg(ref_gene)
  
  # set sequence names
  if (chrXY) {
    genome_info <- data.frame(chr = paste("chr", c(1:22, "X", "Y"), sep = ""), 
                              stringsAsFactors = FALSE)
  } else {
    genome_info <- data.frame(chr = paste("chr", 1:22, sep = ""),
                              stringsAsFactors = FALSE)
  }
  
  rownames(genome_info) <- genome_info$chr
  
  if(ref_gene == "hg19"){
    cur_tbl_UCSC = tbl_ucsc
    message("using hg19 genome annotations from UCSC")
  } else {
    # TODO: How to access tbl_UCSC_hg38 properly?
    cur_tbl_UCSC = conumee2.modified::tbl_UCSC_hg38
    message("using hg38 genome annotations from UCSC")
  }
  
  
  tbl.chromInfo <- cur_tbl_UCSC$chromInfo[match(genome_info$chr, cur_tbl_UCSC$chromInfo$chrom),
                                          "size"]
  genome_info$size <- tbl.chromInfo
  
  tbl.gap <- cur_tbl_UCSC$gap[is.element(cur_tbl_UCSC$gap$chrom, genome_info$chr),]
  
  
  gap_ranges <- sort(GRanges(as.vector(tbl.gap$chrom), IRanges(tbl.gap$chromStart +
                                                                 1, tbl.gap$chromEnd),
                             seqinfo = Seqinfo(genome_info$chr, genome_info$size)))
  
  
  
  if(ref_gene == "hg19"){
    # find the gap that is overlapping with the end of the last p-band, use
    # center of that gap for indicating centromers in the genome plots
    
    tbl.cytoBand <- cur_tbl_UCSC$cytoBand[is.element(cur_tbl_UCSC$cytoBand$chrom,
                                                     genome_info$chr), ]
    
    pq <- sapply(split(tbl.cytoBand$chromEnd[grepl("p", tbl.cytoBand$name)],
                       as.vector(tbl.cytoBand$chrom[grepl("p", tbl.cytoBand$name)])),
                 max)
    
    genome_info$pq <- start(resize(subsetByOverlaps(gap_ranges, 
                                                    GRanges(names(pq),
                                                            IRanges(pq, pq))), 
                                   1, 
                                   fix = "center"))
    
  } else if (ref_gene == "hg38"){
    # use centromers table directly (centromers are not included in gap table for hg38 anymore)
    # combine multiple centromere entries per chromosome and get most extreme values
    centromere_table = cur_tbl_UCSC$centromeres %>% 
      group_by(chrom) %>% 
      dplyr::summarise(start = min(chromStart), end = max(chromEnd)) %>% 
      mutate(middle = (start+end)/2)
    
    genome_info$pq = centromere_table$middle[match(genome_info$chr, centromere_table$chrom)]
    
    
  } else {
    stop("reference genome not recognized")
  }
  
  out = list(genome = genome_info, gap = gap_ranges)
  return(out)
}



#### main functions ####
#' @import minfi
#' @import IlluminaHumanMethylation450kmanifest
#' @import IlluminaHumanMethylation450kanno.ilmn12.hg19
#' @import IlluminaHumanMethylationEPICmanifest
#' @import IlluminaHumanMethylationEPICanno.ilm10b4.hg19
#' @import GenomicRanges
#' @import IRanges
#' @import GenomeInfoDb
#' @importFrom rtracklayer import
NULL



#' CNV.create_anno
#' @description Create annotations for CNV analysis.
#' @param bin_minprobes numeric. Minimum number of probes per bin. Bins are iteratively merged with neighboring bin until minimum number is reached.
#' @param bin_minsize numeric. Minimum size of a bin.
#' @param bin_maxsize numeric. Maximum size of a bin. Merged bins that are larger are filtered out.
#' @param array_type character. One of \code{450k}, \code{EPIC}, \code{EPICv2}, \code{mouse} or \code{overlap}. Defaults to \code{450k}.
#' @param exclude_regions GRanges object or path to bed file containing genomic regions to be excluded.
#' @param detail_regions GRanges object or path to bed file containing genomic regions to be examined in detail.
#' @param chrXY logical. Should chromosome X and Y be included in the analysis?
#' @return \code{CNV.anno} object.
#' @details This function collects all annotations required for CNV analysis using Illumina 450k, EPIC or Mouse arrays. The output \code{CNV.anno} object is not editable. Rerun \code{CNV.create_anno} to change parameters.
#' @examples
#' # create annotation object
#' anno <- CNV.create_anno(array_type = "450k", detail_regions = detail_regions, exclude_regions = exclude_regions)
#' anno
#' @author Volker Hovestadt \email{conumee@@hovestadt.bio}, Bjarne Daenekas
#' @export
CNV.create_anno <- function(bin_minprobes = 15, bin_minsize = 50000, bin_maxsize = 5e+06,
    array_type = "450k", ref_gene = c("hg19", "hg38"), exclude_regions = NULL, detail_regions = NULL, chrXY = FALSE, 
    custom_probe_set = NULL, cg_probes_only = FALSE) {
    object <- new("CNV.anno")
    object@date <- date()

    ref_gene = match.arg(ref_gene)
    
    a1 <- formals()
    a2 <- as.list(match.call())[-1]
    object@args <- as.list(sapply(unique(names(c(a1, a2))), function(an) if (is.element(an,
        names(a2)))
        a2[[an]] else a1[[an]], simplify = FALSE))

    if (is.null(array_type)) {
      array_type <- "450k"
    }
    if (!is.element(array_type, c("450k", "EPIC","EPICv2", "mouse","overlap"))) {
      stop("array_type must be on of 450k, EPIC, EPICv2, mouse or overlap")
    }


    if(array_type == "mouse") {
      object = .create_anno_mouse(anno_object = object, 
                                  array_type = array_type,
                                  chrXY = chrXY, 
                                  exclude_regions = exclude_regions, 
                                  detail_regions = detail_regions,
                                  bin_minprobes = bin_minprobes, 
                                  bin_minsize = bin_minsize, 
                                  bin_maxsize = bin_maxsize)
    } else {
      object = .create_anno_human(anno_object = object, 
                                  array_type = array_type,
                                  chrXY = chrXY, 
                                  exclude_regions = exclude_regions,
                                  detail_regions = detail_regions,
                                  bin_minprobes = bin_minprobes, 
                                  bin_minsize = bin_minsize, 
                                  bin_maxsize = bin_maxsize,
                                  custom_probe_set = custom_probe_set,
                                  ref_gene = ref_gene,
                                  cg_probes_only = cg_probes_only)
    }
    
    return(object)
}


#' CNV.create_bins
#' @description Split genome into bins of defined size.
#' @param hg19.anno foo
#' @param bin_minsize foo
#' @param hg19.gap foo
#' @param hg19.exclude foo
#' @return \code{GRanges} object.
CNV.create_bins <- function(hg19.anno, bin_minsize = 50000, hg19.gap, hg19.exclude) {
  hg19.tile <- sort(tileGenome(Seqinfo(hg19.anno$chr, hg19.anno$size),
                               tilewidth = bin_minsize, cut.last.tile.in.chrom = TRUE))
  # setdiff for gaps (on every second window to avoid merging)
  hg19.tile <- sort(c(setdiff(hg19.tile[seq(1, length(hg19.tile), 2)],
                              hg19.gap), setdiff(hg19.tile[seq(2, length(hg19.tile), 2)], hg19.gap)))
  # setdiff for exluded regions
  hg19.tile <- sort(c(setdiff(hg19.tile[seq(1, length(hg19.tile), 2)],
                              hg19.exclude), setdiff(hg19.tile[seq(2, length(hg19.tile), 2)],
        hg19.exclude)))
    return(hg19.tile)
}

#' CNV.merge_bins
#' @description Merge bins containing less than the defined number probes with neighboring bin containing fewer probes.
#' @param hg19.anno foo
#' @param hg19.tile foo
#' @param bin_minprobes foo
#' @param hg19.probes foo
#' @param bin_maxsize foo
#' @param verbose foo
#' @return \code{GRanges} object.
CNV.merge_bins <- function(hg19.anno, hg19.tile, bin_minprobes = 15, hg19.probes,
    bin_maxsize = 5e+06, verbose = FALSE) {
    values(hg19.tile)$probes <- countOverlaps(hg19.tile, hg19.probes)
    hg19.tile.df <- as.data.frame(hg19.tile)[, c("seqnames", "start", "end",
        "probes")]
    hg19.tile.df$seqnames <- as.vector(hg19.tile.df$seqnames)  # not factor

    hg19.tile.df.bin <- do.call(rbind, lapply(split(hg19.tile.df, hg19.tile.df$seqnames),
        function(mdf) {
            while (min(mdf$probes) < bin_minprobes) {
                mw <- which(mdf$probes == min(mdf$probes))[1]
                mwn <- NA
                mwns <- Inf
                # left
                if (is.element(mdf[mw, "start"] - 1, mdf[, "end"])) {
                  mwn <- mw - 1
                  mwns <- mdf[mw - 1, "probes"]
                  # }
                }
                # right
                if (is.element(mdf[mw, "end"] + 1, mdf[, "start"])) {
                  if (mdf[mw + 1, "probes"] < mwns) {
                    mwn <- mw + 1
                    mwns <- mdf[mw + 1, "probes"]
                  }
                }
                if (is.na(mwn)) {
                  if (verbose)
                    message(paste(mdf[mw, 1:3], collapse = "-"), " has only ",
                      mdf[mw, "probes"], " probes and cannot be merged - remove")
                  mdf <- mdf[-mw, ]
                } else {
                  # merge
                  mdf[mwn, "start"] <- min(mdf[c(mwn, mw), "start"])
                  mdf[mwn, "end"] <- max(mdf[c(mwn, mw), "end"])
                  mdf[mwn, "probes"] <- sum(mdf[c(mwn, mw), "probes"])
                  mdf <- mdf[-mw, ]
                }
            }
            return(mdf)
        }))

    hg19.tile.bin <- sort(GRanges(hg19.tile.df.bin$seqnames, IRanges(hg19.tile.df.bin$start,
        hg19.tile.df.bin$end), seqinfo = seqinfo(hg19.tile)))
    hg19.tile.bin <- hg19.tile.bin[width(hg19.tile.bin) <= bin_maxsize]

    values(hg19.tile.bin)$probes <- countOverlaps(hg19.tile.bin, hg19.probes)
    values(hg19.tile.bin)$midpoint <- as.integer(start(hg19.tile.bin) +
        (end(hg19.tile.bin) - start(hg19.tile.bin))/2)
    # values(hg19.tile.bin)$offset <-
    # hg19.anno[as.vector(seqnames(hg19.tile.bin)),
    # 'offset']+values(hg19.tile.bin)$midpoint

    names(hg19.tile.bin) <- paste(as.vector(seqnames(hg19.tile.bin)), formatC(unlist(lapply(table(seqnames(hg19.tile.bin)),
        seq)), width = nchar(max(table(seqnames(hg19.tile.bin)))), format = "d",
        flag = "0"), sep = "-")
    return(hg19.tile.bin)
}

