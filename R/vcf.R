#' Read a simulated VCF file.
#'
#' @param path Path to a VCF file.
#'
#' @return VCF object.
#'
#' @seealso VariantAnnotation::readVCF
#'
#' @export
read_vcf <- function(path) {
  VariantAnnotation::readVcf(path)
}


#' Load mutation information (ignoring genotypes).
#'
#' @param vcf VCF object.
#' @param mut_type Mutation type as an integer.
#' @param pop_origin Integer ID of a population of origin.
#' @param t_min Lower bound on time of origin.
#' @param t_max Upper bound on time of origin.
#'
#' @return GRanges object.
#'
#' @importFrom magrittr %>%
mut_info <- function(vcf, mut_type = NULL, pop_origin = NULL, t_min = -Inf, t_max = Inf) {
  mut_pos <- filter_muts(vcf, mut_type, pop_origin, t_min, t_max)

  gr <- GenomicRanges::granges(vcf)[mut_pos]
  names(gr) <- NULL

  sample_ids <- VariantAnnotation::samples(VariantAnnotation::header(vcf))
  GenomicRanges::mcols(gr) <- as.data.frame(VariantAnnotation::info(vcf)) %>%
    dplyr::filter(mut_pos) %>%
    dplyr::mutate(freq = AC / (2 * length(sample_ids)))

  # shift VCF coordinates back to the SLiM 0-based system
  gr <- GenomicRanges::shift(gr, shift = -1)

  sort(gr)
}


#' Load mutation genotypes.
#'
#' @param vcf VCF object.
#' @param mut_type Mutation type as an integer.
#' @param pop_origin Integer ID of a population of origin.
#' @param t_min Lower bound on time of origin.
#' @param t_max Upper bound on time of origin.
#'
#' @return GRanges object.
#'
#' @export
#'
#' @importFrom magrittr %>%
mut_gt <- function(vcf, mut_type = NULL, pop_origin = NULL, t_min = -Inf,
                   t_max = Inf) {
  mut_pos <- filter_muts(vcf, mut_type, pop_origin, t_min, t_max)

  gr <- GenomicRanges::granges(vcf)[mut_pos]

  # get the GTs at given sites
  gt_mat <- VariantAnnotation::geno(vcf)$GT[mut_pos, ]
  hap_mat <- as.data.frame(split_haplotypes(gt_mat))

  # get the INFO columns
  info_df <- mut_info(vcf, mut_type, pop_origin, t_min, t_max) %>%
    GenomicRanges::mcols() %>%
    as.data.frame

  GenomicRanges::mcols(gr) <- dplyr::bind_cols(info_df, hap_mat)
  names(gr) <- NULL

  # shift VCF coordinates back to the SLiM 0-based system
  gr <- GenomicRanges::shift(gr, shift = -1)

  sort(gr)
}


#' Transpose the simulated sites in the SLiM 0-based coordinate system into
#' realistic coordinate system.
#'
#' SLiM simulates all mutations as part of a single continuous "chromosome".
#' Therefore, simulating fixed mutations that are present in a real genome
#' requires their transposition into a 0-based continuous coordinate system.
#' This function performs this transposition based on the given positions of
#' the sites in their original coordinate system.
#'
#' @param sim_sites GRanges object with SLiM VCF output data.
#' @param real_sites GRanges object with real coordinates of sim_sites.
#'
#' @return GRanges object with transposed coordinates.
#'
#' @seealso read_sites
#'
#' @export
transpose_sites <- function(sim_sites, real_sites) {
  # for each simulated 0-based SLiM site, find a corresponding real coordinate
  hits <- GenomicRanges::findOverlaps(real_sites, sim_sites)

  # get the coordinates of the realistic sites from the DataFrame portion
  # of the GRanges object
  transposed_df <- as.data.frame(GenomicRanges::mcols(real_sites))
  names(transposed_df) <- c("chrom", "start", "end")

  # convert the realistic coordinates into a GRanges object
  transposed_gr <- GenomicRanges::makeGRangesFromDataFrame(
      transposed_df, starts.in.df.are.0based = TRUE
  )[IRanges::from(hits)]

  # assign the INFO/GT DataFrame object of the simulated sites to the newly
  # generated realistic coordinates GRanges object
  GenomicRanges::mcols(transposed_gr) <- GenomicRanges::mcols(sim_sites)

  transposed_gr
}


#' Read positions of simulated sites in their original coordinate system.
#'
#' SLiM simulates all mutations as part of a single continuous "chromosome".
#' Therefore, simulating fixed mutations that are present in a real genome
#' requires their transposition into a 0-based continuous coordinate system.
#' This function reads a table in a BED-like format that makes it possible
#' to "reverse transpose" the simulated 0-based coordinates into their
#' original form using the "transpose_sites" function.
#'
#' The file has a strict (tab-separated) format and must contain the following
#' columns in the same order: original_chromosome_id original_start_pos
#' original_end_pos SLiM_start_pos SLiM_end_pos.
#'
#' @param file Path to a BED-like file. See the required format bellow.
#'
#' @return GRanges object.
#'
#' @seealso admixr::transpose_sites
#'
#' @export
#'
#' @importFrom magrittr %>%
read_sites <- function(file) {
  gr <- readr::read_tsv(file, col_types = "ciiiic") %>%
    dplyr::select(
      real_chrom = chrom,
      real_start = start,
      real_end = end,
      start = slim_start,
      end = slim_end
    ) %>%
    dplyr::mutate(chrom = 1) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
  gr
}


#' Given a chromosome with fixed ancestry markers, infer coordinates of
#' introgressed tracks.
#'
#' @param markers GRanges object with marker genotypes.
#' @param chroms Identifiers of chromosomes to analyze.
#'
#' @return GRanges object with coordinates of detected tracks.
#'
#' @export
#'
#' @importFrom magrittr %>%
ancestry_tracks <- function(markers, chroms = NULL) {
  # if no chromosomes were specified, look for ancestry tracks in all
  if (is.null(chroms)) {
    cols <- colnames(GenomicRanges::mcols(markers))
    chroms <- cols[grep("^chr", cols)]
  }

  lapply(chroms, function(chrom) {

  marker_runs <- rle(as.integer(GenomicRanges::mcols(markers)[[chrom]] == 1))

  # if there are no informative markers on this chromosome, there are no
  # introgressed tracks
  if (!any(marker_runs$values)) return(GenomicRanges::GRanges())

  block_idx <- c(0, cumsum(marker_runs$lengths))
  block_start <- block_idx[-length(block_idx)]
  block_end <- block_idx[2:length(block_idx)]

  # a single run of markers implies a 100% of ancestry
  anc_pos <- if (length(marker_runs$values) > 1) marker_runs$values == 1 else 1
  track_start <- block_start[anc_pos]
  track_end <- block_end[anc_pos]

  anc_haps <- GenomicRanges::GRanges(
    unique(GenomicRanges::seqnames(markers)),
    IRanges::IRanges(start = GenomicRanges::start(markers[track_start + 1, ]),
                     end   = GenomicRanges::start(markers[track_end,       ]))
  )

  anc_haps

  }) %>% GenomicRanges::GRangesList() %>% setNames(chroms)
}
