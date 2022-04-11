#' read readcount output with correct column names
#' @param Input char. path to RC table
#' @keywords readcount RC BRC
#' @export
#' @examples
#' read.readcount("./path/to/RC.txt") -> RC
#'
#'@description use the folloing line to genereate readcount table:
#' $ bam-readcount -w0 -d 1000000000 -f ref.fa alignment.sorted.bam | awk -F ":|\t|=" 'BEGIN {OFS = "\t"}; {print $1, $2, $3 , $4, $21 , $35, $49 , $63}' > BRC.txt

read.readcount <- function(X){
read.table(X, sep="\t", col.names=c("Chr.","POS","REF","COV","A","C","G","T"))
}