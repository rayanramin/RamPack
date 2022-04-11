#'DESeq2 result table
#'@description reformats the output of DESeq2 to a more useful format.
#' @export

DEResult_table <- function(input , alpha = 0.05){
  input %>% data.frame() %>%
    rownames_to_column(var="gene") %>%
    as_tibble() %>% filter(!is.na(padj)) %>%  mutate(pval=pvalue,fold_change=2^log2FoldChange , sig = ifelse( padj <=  alpha  , "Significant" , "N.S.")  , diff_exp = ifelse( log2FoldChange  > 1 , "UP" , ifelse( log2FoldChange  < -1 , "DOWN", "NO" )) )
}
