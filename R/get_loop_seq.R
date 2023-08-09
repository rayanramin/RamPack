#' Get hairpin loop sequence
#' @description gets the sequence of hairpin loops from a data.frame
#' @param X input must be a data.frame of Hairpin_survey output with columns looplen,looppos, minus2,minus1,minus0,plus1,plus2,plus3
#' output will be of character strings like "(NATC.C.C)" with the targeted C shown as .C.
#' @example get_loop_seq(Data) -> Data$loopseq
#' @export



get_loop_seq <- function(X){
  ACGT <- c("A","C","G","T")
  tmp <- paste0("NNNNNNNNNN",ACGT[match(X$minus2,1:4)],ACGT[match(X$minus1,1:4)],ACGT[match(X$minus0,1:4)],".C.",ACGT[match(X$plus1,1:4)],ACGT[match(X$plus2,1:4)],ACGT[match(X$plus3,1:4)],"NNNNNNNNNN")
  out <- rep(NA,nrow(X))
  ifelse(X$looppos==0,paste0(substr(tmp,14,16),"(",paste0(substr(tmp,17,17+X$looplen-1),")")),
         ifelse(X$looppos==X$looplen+1,paste0("(",substr(tmp,14-X$looplen,13),")",substr(tmp,14,16)),
                paste0("(",substr(tmp,15-X$looppos,16+X$looplen-X$looppos),")")))
}
duplicated_all <- function(x){
  tmp <- x[duplicated(x)]
  x %in% tmp
}


get_hairpin_seq <- function(D,stem_length=1, onlyLoop=FALSE){

  ## load the genome
  stopifnot(!require("BSgenome.Ecoli.WayneState.BH214V4", quietly = TRUE))
  library(BSgenome.Ecoli.WayneState.BH214V4)

  D %>% mutate(roc = looplen-looppos, loc = looppos-1) %>%
    mutate(st = ifelse(ref==2, (POS-loc), (POS-roc)) , en = ifelse(ref==2, POS+roc, POS+loc)) %>%
    mutate(stem_st = st - stem_length , stem_en = en + stem_length ) -> D

  D$seq <- NA
  ## to deal with circular genome
  # beginning
  D$seq[which(D$stem_st <=0)] <- lapply(which(D$stem_st <=0), function(i) paste0(BH214$BH214V4[(length(BH214$BH214V4)+D$stem_st[i]):length(BH214$BH214V4)] , BH214$BH214V4[1:D$stem_en[i]])) %>% unlist %>% toupper
  ## end
  D$seq[which(D$stem_en >length(BH214$BH214V4))] <- lapply(which(D$stem_en >length(BH214$BH214V4)), function(i) paste0(BH214$BH214V4[(D$stem_st[i]):length(BH214$BH214V4)] , BH214$BH214V4[1:(D$stem_en[i]-length(BH214$BH214V4))])) %>% unlist %>% toupper
  ## middle of the sequence
  D$seq[ which(D$stem_st >0 & D$stem_en <=length(BH214$BH214V4))] <- lapply(which(D$stem_st >0 & D$stem_en <=length(BH214$BH214V4)), function(i) paste0(BH214$BH214V4[(D$stem_st[i]):D$stem_en[i]])) %>% unlist %>% toupper

  D$seq[D$ref==3] <- lapply(which(D$ref==3), function(i) seqinr::c2s(rev(seqinr::comp(seqinr::s2c(D$seq[i])))) ) %>% unlist %>% toupper

  if (any(!(substr(D$seq,stem_length+D$looppos,stem_length+D$looppos)=="C"))) {
    stop("something went wrong")
  } else {
    SEQS <-  lapply(D$seq, seqinr::s2c)
    lapply(1:length(SEQS), function(i) {
      SEQS[[i]][stem_length+D$looppos[i]] <- ".C." ;
      SEQS[[i]][1+stem_length] <-  paste0('(',SEQS[[i]][1+stem_length] ) ;
      SEQS[[i]][stem_length+D$looplen[i]] <-  paste0(SEQS[[i]][stem_length+D$looplen[i]],')' );
      seqinr::c2s(SEQS[[i]]) }) %>%  unlist -> out
  }
  if(onlyLoop){
    stringr::str_extract_all(out,"\\(.*?(.C.)?\\)|\\.C\\.",simplify = TRUE) %>%
      apply(. ,1,paste,collapse="") -> out
  }
  return(out)
}
