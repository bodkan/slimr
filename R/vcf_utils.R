# Filter positions based on given mutation criteria.
filter_muts <- function(vcf, mut_type = NULL, pop_origin = NULL, t_min = -Inf, t_max = Inf) {
  vcf_info <- VariantAnnotation::info(vcf)

  mut_pos <- vcf_info$GO >= t_min & vcf_info$GO <= t_max

  if (!is.null(mut_type)) {
    mut_pos <- mut_pos & vcf_info$MT == mut_type
  }

  if (!is.null(pop_origin)) {
    mut_pos <- mut_pos & (vcf_info$PO == pop_origin)
  }

  mut_pos
}


# Split the diploid GT matrix into haploid GT matrix.
# (I.e. split diploid "0/1"-like GTs into two haploid "0" and "1" GTs.)
split_haplotypes <- function(gt_mat) {
  hap_mat <- matrix(nrow = nrow(gt_mat), ncol = ncol(gt_mat) * 2)
  colnames(hap_mat) <- paste0("chr", 1:ncol(hap_mat))

  for (i in seq_len(ncol(gt_mat))) {
    k <- (2 * i) - 1
    hap_mat[, c(k, k + 1)] <-
      as.integer(stringr::str_split_fixed(gt_mat[, i], "\\|", 2))
  }

  hap_mat
}


#' Download the coordinates of centromeres.
#' These are used for "splitting" admixture tracts that overlap a centromere
#' and would therefore appear extensively long.
#' http://rohsdb.usc.edu/GBshape/cgi-bin/hgTables?db=hg19&hgta_group=map&hgta_track=gap&hgta_table=gap&hgta_doSchema=describe+table+schema
#'
#' @export
get_gaps <- function() {
  session <- rtracklayer::browserSession("UCSC")
  rtracklayer::genome(session) <- "hg19"

  query <- rtracklayer::ucscTableQuery(
    session,
    "gap",
    rtracklayer::GRangesForUCSCGenome("hg19", chrom = paste0("chr", 1:22))
  )
  tbl <- rtracklayer::getTable(query)

  gaps_df <- dplyr::filter(tbl,
                           chrom %in% paste0("chr", 1:22))
  gaps_gr <- GenomicRanges::makeGRangesFromDataFrame(gaps_df, keep.extra.columns = TRUE)

  IRanges::reduce(gaps_gr)
}
