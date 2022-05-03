#' Detect sequence context (stranded)
#' @description  This function detects contexts and returns a vector of -1, 1 and 0 showing the presense and absense of the selected context. -1 and +1 shows the strand that cytosine is on.
#' This function is not optimized for speed.
#' @param ref_fasta  fasta file
#' @param chr chromosome(s) to be looked up
#' @param pos position(s) to be looked up
#' @param context context to be detected, accepts IUPAC nucleotide code,
#'  if more than 1 C in the context, the target C must be denoted by using an X,
#'  for example, "TCX" looks for cytosines preceded by TC. And "XC" looks for cytosines followed by a C.
#'
#' @note Corrently this function would return an error for out of bound ranges,
#'  such as positions at the edges of the sequence.
#'
#' @export
#' @examples
#' Detect_Context("./ref.fa",chr="chr2",pos=c(1234,2345),"WRCY") -> out
#'

Detect_Context <- function(ref_fasta,chr,pos,context){
  if (missing(ref_fasta)) {
    stop("ERROR: input ref_fasta is required")}
  if (missing(chr)) {
    stop("ERROR: input chr is required")}
  if (missing(pos)) {
    stop("ERROR: input pos is required")}
  if (missing(context)) {
    stop("ERROR: input context is required")}
  if(stringr::str_count(context,"C")>1 & stringr::str_count(context,"X")==0){
    stop('ERROR: wrong context, if more than 1 "C" in the context, you should specify the target "C" by an "X"')}
  if(stringr::str_count(context,"X")>1){
    stop('ERROR: wrong context, you can only use 1 "X" to denote the target "C"')}
  if(stringr::str_count(context,"X")==1){
  context <- toupper(context)
    l <- nchar(context)
    xp <- stringr::str_locate(context,"X")[1]
    ba = l-xp; bb = xp-1;
    if(bb < ba){ context <- paste0(paste0(rep("N",ba-bb),collapse = ""),context)}
    if(bb > ba){ context <- paste0(context, paste0(rep("N",bb-ba),collapse = ""))}
    context <- stringr::str_replace(context,"X","C")
  }
  if(stringr::str_count(context,"X")==0){
    l <- nchar(context)
    xp <- stringr::str_locate(context,"C")[1]
    ba = l-xp; bb = xp-1;
    if(bb < ba){ context <- paste0(paste0(rep("N",ba-bb),collapse = ""),context)}
    if(bb > ba){ context <- paste0(context, paste0(rep("N",bb-ba),collapse = ""))}
  }
  flank <-  max(ba,bb)
  context <- bioseq::dna(context)
  revcom_context <-  bioseq::seq_complement(bioseq::seq_reverse(context))
  GenomicRanges::GRanges(paste0(chr,":",pos-flank,"-",pos+flank)) -> g
  Rsamtools::scanFa(ref_fasta, param=g)  %>%
  as.character %>% unname %>% bioseq::dna() %>%
    {case_when(bioseq::seq_detect_pattern(.,pattern = context) ~ 1,
               bioseq::seq_detect_pattern(.,pattern =  revcom_context) ~ -1,
               TRUE ~ 0)}
}
