#' River Plots
#' @description
#' RivPlot is a wrapper for the modified plot_river2 function.
#' the original plot_river is taken from the MutationalPatterns Package
#' @param input mutation count matrix from the get_mut_context_matrix() function
#' @param fontsize font size of the text of the nucleotide/mutation labels
#' @param is_ui are you plotting UI values or mutation count?
#' @param is_norm is the data normalized?
#' @param title Plot's title
#' @param AltCol using alternative colors?
#' @param facet a character vector to be used for facet labels

RivPlot <- function(input, fontsize=4,is_ui=FALSE,is_norm=FALSE,title="",AltCol=T, facet=NULL){
  plot_river2(input ,UI=is_ui, is_normalized=is_norm, fontsize=fontsize,alt_colors=AltCol, FacetNames=facet) +
    theme(legend.position = "none",
          panel.spacing.y = unit(1, "line"),
          strip.background = element_blank(),
          strip.text = element_text(size=18),
          plot.title = element_text(size=22),
          axis.title = element_text(size=20),
          axis.text.y = element_text(size=17),
          axis.text.x = element_text(size=19))+
    labs(title=title ) +
    scale_y_continuous(expand = c(0.01,0.01)) +
    scale_x_discrete(expand = c(0.03,0.03)) -> g
  return(g)
}


########################
#' plot_river2
#' @description
#' plot_river2 is modified version of plot_river from MutationalPatterns
#' @param name description
#' @references https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-018-0539-0
#'
#' #' Plot a riverplot
#'
#' Function to plot a SNV mutation matrix as a riverplot.
#' This is especially useful when looking at a wide
#' mutational context
#'
#' @param mut_matrix Matrix containing mutation counts.
#' @param condensed More condensed plotting format. Default = F.
#' @param UI set to TRUE if ploting Uracilation Index Values
#' @param is_normalized set to TRUE if the input data is normalized, to change labels
#' @param fontsize font size of the text of the nucleotide/mutation labels
#' @param alt_colors if set tot TRUE, alternative colors will be used
#' @param FacetNames a character vector to be used for facet labels
#'
#' @return A ggplot object
#'
#' @export
#' @import ggplot2
#' @import ggalluvial
#'
#' @seealso
#' \code{\link{mut_matrix}},
#' \code{\link{plot_96_profile}},
#' \code{\link{plot_profile_heatmap}}
#'
#' @examples
#'
#' ## See the 'mut_matrix()' examples for how we obtained the
#' ## mutation matrix information:
#' ## Get regular matrix
#' mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' ## Create heatmap of profile
#' plot_river(mut_mat[,c(1,4)])
#'
#' ## Get extended matrix
#' mut_mat_extended <- readRDS(system.file("states/mut_mat_data_extended.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' ## Create heatmap of extended profile
#' plot_river(mut_mat_extended[,c(1,4)])
#'
#' ## Create condensed version of riverplot
#' plot_river(mut_mat_extended[,c(1,4)], condensed = TRUE)
#'
plot_river2 = function(mut_matrix, condensed = FALSE, UI=FALSE, is_normalized=FALSE,fontsize = 3, alt_colors = FALSE, FacetNames=NULL){

  context_m <- .river_create_context_matrix(mut_matrix)
  lodes_tb <- .river_create_lodes(mut_matrix, context_m)
  fig <- .river_plot_river(lodes_tb, mut_matrix, condensed,UI,is_normalized,fontsize,alt_colors,FacetNames)

  return(fig)
}

#' Split the contexts of a mutation matrix
#'
#' Function to split the contexts of a mutation matrix
#' in multiple columns.
#'
#' @param mut_matrix Matrix containing mutation counts.
#'
#' @return Matrix containing the contexts of a mutation matrix.
#' @noRd
#' @importFrom magrittr %>%
#'
.river_create_context_matrix = function(mut_matrix){

  # Split left and right context from mutations
  context_m <- mut_matrix %>%
    rownames() %>%
    stringr::str_split("\\[|\\]", simplify = TRUE)

  # Split contexts into single bases
  left_m <- stringr::str_split(context_m[,1], "", simplify = TRUE)
  right_m <- stringr::str_split(context_m[,3], "", simplify = TRUE)

  # Combine all
  context_m <- cbind(left_m, context_m[,2], right_m)

  # Set column names
  mut_positions <- c((-1*rev(seq_len(ncol(left_m)))), 0, seq_len(ncol(right_m)))
  colnames(context_m) <- paste0("pos_", mut_positions)

  return(context_m)
}


#' Create a long format lodes tibble.
#'
#' @param mut_matrix Matrix containing mutation counts.
#' @param context_m Matrix containing the contexts of a mutation matrix.
#'
#' @return A tibble in lodes format containing the mutation contexts.
#' @noRd
#'
#' @importFrom magrittr %>%
#'
.river_create_lodes = function(mut_matrix, context_m){

  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  pos <- type <- sample_name <- rowid <- nrmuts <- NULL

  mut_df <- as.data.frame(mut_matrix)
  context_df <- as.data.frame(context_m)

  # Transform data into long (lodes) format.
  lodes_tb <- cbind(mut_df, context_df) %>%
    tibble::rowid_to_column() %>%
    tidyr::pivot_longer(cols = dplyr::contains("pos_"), names_to = "pos", values_to = "type") %>%
    tidyr::pivot_longer(cols = c(-pos, -type, -rowid),
                        names_to = "sample_name",
                        values_to = "nrmuts") %>%
    dplyr::mutate(pos = stringr::str_remove(pos, "pos_"),
                  pos = factor(pos, levels = unique(pos)),
                  type = factor(type, levels = c( "A","C","G","T","C>A", "C>G","C>T","C>U", "T>A", "T>C", "T>G")),
                  sample_name = factor(sample_name, levels = unique(sample_name)))
  return(lodes_tb)
}

#' Plot a riverplot with a lodes tibble as input
#'
#' @param lodes_tb A tibble in lodes format containing the mutation contexts.
#' @param condensed More condensed plotting format.
#' @param UI If TRUE, plot UI instead of mutation counts.
#'
#' @return A ggplot object
#' @noRd
#'
#' @import ggplot2
#' @import ggalluvial
#'
.river_plot_river = function(lodes_tb, mut_matrix, condensed,UI,is_normalized,fontsize,alt_colors,FacetNames){

  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  rowid <- nrmuts <- NULL


  # Set colours
  if (alt_colors == TRUE) {
    colours_v <- c("#41BCFF","#fe3b6b","#1F64FF","#e21f00",
                   "#ffd0d4","#d3c7d5","#B4B4B4","#B4B4B4",
                   "#d2ead2","#9dc19b","#c4c5ab")

    lodes_tb$type <- factor(lodes_tb$type, levels = c("A","G","C","T","C>A", "C>G","C>T","C>U", "T>A", "T>C", "T>G"))
  } else {
    colours_v <- c("#EC1456","#71AC94","#FEBC64","#37B9ED",
                   "#4464F6","#000000","#DB0016","#DB0016",
                   "#D4D2D2","#AECF55","#EFCFCE")
  }


  TYPES <- c( "A","C","G","T","C>A", "C>G", "C>T", "C>U", "T>A", "T>C", "T>G")
  names(colours_v) <- TYPES
  used_colours <- colours_v[which(names(colours_v) %in% unique(as.character(lodes_tb$type)))]
  USED_TYPES <- TYPES[which(names(colours_v) %in% unique(as.character(lodes_tb$type)))]

  # Change plotting parameters based on whether plot should be condensed.
  if (condensed == TRUE) {
    spacing <- 0
  } else {
    spacing <- 0.5
  }

  # Create facet texts
  if(UI==TRUE){
    uimean <- round(colMeans(mut_matrix),2)
    facet_labs_y <- stringr::str_c(colnames(mut_matrix), " (mean UI = ", uimean, ")")
    names(facet_labs_y) <- colnames(mut_matrix)
  } else if(is_normalized == TRUE){
    facet_labs_y <- stringr::str_c(colnames(mut_matrix))
    names(facet_labs_y) <- colnames(mut_matrix)
  }else{
    nr_muts <- colSums(mut_matrix)
    facet_labs_y <- stringr::str_c(colnames(mut_matrix), " (n = ", nr_muts, ")")
    names(facet_labs_y) <- colnames(mut_matrix)
  }

  if(!is.null(FacetNames)){
    facet_labs_y <- FacetNames
    names(facet_labs_y) <- colnames(mut_matrix)
  }

  #Add stratum stat.
  StatStratum <- ggalluvial::StatStratum

  # Create plot
  fig <- ggplot(lodes_tb, aes(x = pos,
                              stratum = type,
                              y = nrmuts,
                              alluvium = rowid,
                              fill =  type,
                              label = type)) +
    ggalluvial::geom_stratum() +
    ggalluvial::geom_flow() +
    geom_text(stat = StatStratum, size = fontsize, colour = "white") +
    facet_grid(sample_name ~ ., scales = "free_y",
               labeller = labeller(sample_name = facet_labs_y)) +
    scale_fill_manual(values = used_colours, breaks=USED_TYPES) +
    labs(fill="",x = "Position", y = "Nr mutations") +
    theme_bw() +
    theme(panel.spacing.y = unit(spacing, "lines"))
  if(UI==TRUE){ fig <- fig + labs(y= "UI")}
  if(is_normalized==TRUE){fig <- fig + labs(y= "Normalized mutation rate")}

  return(fig)
}



