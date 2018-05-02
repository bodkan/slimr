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
#' @export
#'
#' @importFrom magrittr %>%
mut_info <- function(vcf, mut_type = NULL, pop_origin = NULL, t_min = -Inf,
                     t_max = Inf) {
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


#' Transpose coordinates in the SLiM 0-based coordinate system into
#' realistic coordinate system.
#'
#' SLiM simulates all mutations as part of a single continuous "chromosome".
#' Therefore, simulating fixed mutations that are present in a real genome
#' requires their transposition into a 0-based continuous coordinate system.
#' This function performs this transposition based on the given positions of
#' the sites in their original coordinate system.
#'
#' The real_coords GRanges object must contain a DataFrame object with three
#' columns specifying a realistic chromosome, start and end coordinate of
#' a region.
#'
#' @param sim_coords GRanges object with SLiM, 0-based coordinates.
#' @param real_coords GRanges object with real coordinate.
#'
#' @return GRanges object with transposed coordinates.
#'
#' @seealso read_coordinates
#'
#' @export
transpose_coordinates <- function(sim_coords, real_coords) {
  # for each simulated 0-based SLiM position, find a corresponding
  # real coordinate
  hits <- GenomicRanges::findOverlaps(real_coords, sim_coords)

  # get the realistic coordinates from a DataFrame portion of a GRanges object
  transposed_df <- as.data.frame(GenomicRanges::mcols(real_coords))
  names(transposed_df) <- c("chrom", "start", "end")

  # convert the realistic coordinates into a GRanges object
  transposed_gr <- GenomicRanges::makeGRangesFromDataFrame(
      transposed_df, starts.in.df.are.0based = TRUE
  )[IRanges::from(hits)]

  # assign the INFO/GT DataFrame object of the simulated coordinates to the
  # newly generated realistic coordinates GRanges object
  GenomicRanges::mcols(transposed_gr) <- GenomicRanges::mcols(sim_coords)

  transposed_gr
}


#' Read positions of simulated coordinates (regions or sites) in their
#' original coordinate system.
#'
#' SLiM simulates everything as part of a single continuous "chromosome".
#' Therefore, matching simulated coordinates to their real coordinates
#' requires their transposition into a 0-based continuous coordinate system.
#' This function reads a table in a BED-like format that makes it possible
#' to "reverse transpose" the simulated 0-based coordinates into their
#' original form using the "transpose_sites" function.
#'
#' The coordinate file must be in a TSV format and must contain columns
#' chrom, start, end, slim_start and slim_end.
#'
#' @param file Path to a BED-like file. See the required format bellow.
#'
#' @return GRanges object.
#'
#' @seealso admixr::transpose_coordinates
#'
#' @export
#'
#' @importFrom magrittr %>%
read_coordinates <- function(file) {
  gr <- read.table(file, header = TRUE, stringsAsFactors = FALSE) %>%
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

#' Extract coordinates of deserts from a given VCF file.
#'
#' @param markers GRanges object with marker genotypes.
#' @param cutoff Frequency cutoff for admixture markers.
#'
#' @return GRanges object with ancestry deserts.
#'
#' @export
ancestry_deserts <- function(markers, cutoff = 0) {
  all_chrom <- list()

  for (chrom in unique(GenomicRanges::seqnames(markers))) {
    chrom_markers <- markers[GenomicRanges::seqnames(markers) == chrom, ]
    desert_runs <- rle(as.integer(chrom_markers$freq > cutoff))

    # if there's not desert on this chromosome (i.e., the whole chromosome
    # comes from an introgressing population), return nothing
    if (length(desert_runs$values) < 2) {
      all_chrom[[chrom]] <- NULL
      next
    }

    block_idx <- c(0, cumsum(desert_runs$lengths))
    block_start <- block_idx[-length(block_idx)]
    block_end <- block_idx[2:length(block_idx)]

    desert_start <- block_start[desert_runs$values == 0]
    desert_end <- block_end[desert_runs$values == 0]

    all_chrom[[chrom]] <- GenomicRanges::GRanges(
        chrom,
        IRanges(start = GenomicRanges::start(chrom_markers[desert_start + 1, ]),
                end = GenomicRanges::start(chrom_markers[desert_end, ]))
    )
  }
  Reduce(c, GRangesList(all_chrom))
}


#' Fill in frequencies of alleles lost during the simulation (to prevent a SFS
#' to be "skewed" at the end of the simulation).
#'
#' @param sim_sites GRanges object of simulated sites.
#' @param all_sites GRanges object of all sites at the start of a simulation.
#'
#' @return GRanges object with frequencies of lost sites set to zero. Otherwise
#'   contains the same data as given by the sim_sites argument.
#'
#' @export
#'
#' @importFrom magrittr %>%
fill_lost <- function(sim_sites, all_sites) {
  sim_df <- dplyr::select(as.data.frame(sim_sites), -strand)
  all_df <- as.data.frame(all_sites)[c("seqnames", "start", "end")]

  joined <- dplyr::full_join(sim_df, all_df,
                             by = c("seqnames", "start", "end")) %>%
    dplyr::mutate(freq=ifelse(is.na(freq), 0, freq)) %>%
    dplyr::mutate_at(dplyr::vars(dplyr::starts_with("chr")),
                     dplyr::funs(ifelse(is.na(.), 0, .)))

  sort(GenomicRanges::makeGRangesFromDataFrame(joined,
                                               keep.extra.columns = TRUE))
}
