#' Get Mutation Context Matrix
#' @description
#' from Grange objects that lists all the mutations.
#' This is a wrapper function to use the modified mut_matrix_UI function
#' This function also allows for assymetrical context to be generated
#'
#' @param input GRange object conatining the mutations
#' @param context defines the extent for with the context of a mutation at C:G position is determined
#' @param UI if set to TRUE, the UI values will be calculated instead of counting the mutations
#' @param custom_name for when there is only one sample, you can set the name to be displayed
#'
#' @references https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-018-0539-0
#'
#' @export


get_mut_context_matrix <- function(input , context = "NNCN", UI = FALSE, custom_name = NULL ){
  bases <- c("A", "C", "G", "T")
  base_subs <- c("[C>A]", "[C>G]", "[C>T]", "[T>A]", "[T>C]", "[T>G]")

  nchar(sub("C.*$", "", context)) -> n1
  nchar(sub("^.*C", "", context)) -> n2
  regex_pattern <- paste0(".{", n1, '}\\[.>.\\].{', n2, "}")

  if(UI){
    mut_mat_ext_context <- mut_matrix_UI(input, ref_genome, extension = max(c(n1,n2)),UI=TRUE)
  } else {
    mut_mat_ext_context <- mut_matrix_UI(input, ref_genome, extension = max(c(n1,n2)),UI=FALSE)
  }

  if(ncol(mut_mat_ext_context)==1 & !is.null(custom_name)){
    colnames(mut_mat_ext_context) <- custom_name
  }

  data.frame(matrix(rep(bases,n1),nrow = 4*6,ncol = n1),base_subs,
             matrix(rep(bases,n2),nrow = 4*6,ncol = n2)) %>%
    mutate_all(as.factor) %>%
    do.call(tidyr::crossing,.) %>% do.call(paste0, .) -> combi_tb
  mut_mat_ext_context %>% as.data.frame %>%
    mutate(context =stringr::str_extract(rownames(.),stringr::regex(regex_pattern) ) ) %>%
    filter(context %in% combi_tb) %>%
    group_by(context) %>% summarise_all(sum) %>%
    filter(.,rowSums(select(.,-context))>0) %>%
    tibble::column_to_rownames("context") %>% as.matrix
}

############################################################################################
#' Make mutation count matrix
#'
#' @description Make mutation count matrix
#'this is a modification of the mut_matrix() function from MutationPatterns package
#' @param vcf_list GRangesList or GRanges object.
#' @param ref_genome BSgenome reference genome object
#' @param extension The number of bases, that's extracted upstream and
#' downstream of the base substitutions. (Default: 1).
#' @param UI if set to TRUE the mean UI value will be calculated
#' @return mutation count matrix
#' @references https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-018-0539-0
#' @export


mut_matrix_UI <- function(vcf_list, ref_genome, extension = 1,UI=FALSE) {

  # Convert list to grl if necessary
  if (inherits(vcf_list, "list")) {
    vcf_list <- GenomicRanges::GRangesList(vcf_list)
  }

  # Determine nr mutations per sample
  if (inherits(vcf_list, "CompressedGRangesList")) {
    gr_sizes <- S4Vectors::elementNROWS(vcf_list)
    gr <- BiocGenerics::unlist(vcf_list)
  } else if (inherits(vcf_list, "GRanges")) {
    gr <- vcf_list
    gr_sizes <- length(gr)
    names(gr_sizes) <- "My_sample"
  } else {
    .not_gr_or_grl(vcf_list)
  }
  # Determine type and context of all mutations
  type_context <- type_context(gr, ref_genome, extension)

  if(UI){
  # get U values
  # Count the type and context to create the mut_mat
  type_context$Uvals <- gr$U
  mut_mat <- mut_96_UI(type_context, gr_sizes,UI=TRUE)
  }else{
  # Count the type and context to create the mut_mat
  mut_mat <- mut_96_UI(type_context, gr_sizes,UI=FALSE)
  }
  return(mut_mat)
}
############################################################################################

#' Count mutation occurrences in the provided contexts
#'
#'  @details
#'  This function is called by mut_matrix_UI. It is a modification of the mut_96_occurrences() function from the MutationalPatterns package.
#'  It calculates the context for all variants and then splits these per GRanges (samples). It then calculates how often each context occurs.
#'  This function is compatible with assymetric context and UI value
#'
#'
#' @param type_context result from type_context function
#' @param gr_sizes A vector indicating the number of variants per GRanges
#' @param UI if TRUE, it will calculated the mean UI instead of counting the mutations
#' @return Mutation matrix with 96 trinucleotide mutation occurrences
#' @references https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-018-0539-0
#' @importFrom magrittr %>%
#'


mut_96_UI <- function(type_context, gr_sizes,UI=FALSE) {

  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  categories <- count <- NULL

  # Determine nr of bases
  nr_bases <- nchar(type_context$context[[1]])
  middle_base <- ceiling(nr_bases / 2)


  # Determine all possible contexts
  bases_left <- c("A", "C", "G", "T")
  bases_right <- c("A", "C", "G", "T")
  base_subs <- c("[C>A]", "[C>G]", "[C>T]", "[T>A]", "[T>C]", "[T>G]")

  # Loop over each base substitution
  full_context_poss <- vector("list", length(base_subs))
  for (i in seq_along(base_subs)) {
    sub <- base_subs[[i]]
    sub_context <- sub
    # Repeatedly add bases left and right
    for (j in seq_len(middle_base - 1)) {
      combi_tb <- tidyr::crossing(bases_left, sub_context, bases_right)
      sub_context <- paste0(combi_tb$bases_left, combi_tb$sub_context, bases_right)
    }
    full_context_poss[[i]] <- sub_context
  }
  full_context_poss <- do.call(c, full_context_poss)


  # Determine 96 context for all variants
  full_context <- stringr::str_c(
    substr(type_context$context, 1, middle_base - 1),
    "[",
    type_context$types,
    "]",
    substr(type_context$context, middle_base + 1, nr_bases)
  ) %>%
    factor(levels = full_context_poss)

  # Set names if they are not yet present
  if (is.null(names(gr_sizes))) {
    names(gr_sizes) <- seq_along(gr_sizes)
  }

  # Create vector describing the sample of each variant
  sample_vector <- rep(names(gr_sizes), gr_sizes) %>%
    factor(levels = names(gr_sizes))

  if(UI){
    # calculated the mean UI
    counts <- tibble::tibble("categories" = full_context, "sample" = sample_vector, U = type_context$U) %>%
      dplyr::filter(!is.na(categories)) %>%
      dplyr::group_by(categories, sample, .drop = FALSE) %>%
      dplyr::summarise(UI = mean(U)) %>%
      dplyr::mutate(UI = ifelse(is.na(UI), 0, UI))
      counts <- tidyr::spread(counts, key = sample, value = UI, fill = 0)
  }else{
    # Count the mutations per type and per sample
    counts <- tibble::tibble("categories" = full_context, "sample" = sample_vector) %>%
      dplyr::filter(!is.na(categories)) %>%
      dplyr::group_by(categories, sample, .drop = FALSE) %>%
      dplyr::summarise(count = dplyr::n())
      counts <- tidyr::spread(counts, key = sample, value = count, fill = 0)
  }
  # Transform the data into a mutation matrix
  unnecesary_cols <- which(colnames(counts) == "<NA>")
  mut_mat <- as.matrix(counts[, -c(1, unnecesary_cols)])
  rownames(mut_mat) <- counts$categories
  return(mut_mat)
}

