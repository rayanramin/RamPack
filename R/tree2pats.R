#' Tree to Pattern
#' @description turn edges(branches) of phylogenetic trees into binary patterns of present/absent of terminal nodes
#' @param infile phylogentic tree
#' @export
#'
#' @example tree2pats(tree) -> out

tree2pats <- function(infile) {
  tree <- infile
  tmp <- as.data.frame(tree$edge)
  tree$tip.label -> labels
  tmp$labels <- NA
  nnum <- tree$Nnode
  tipnum <- length(labels)
  pats <- matrix(F,nrow(tmp),tipnum)
  X <- list(edges=tmp,pats=pats)
  for(i in 1:tipnum){
    j<- which(X$edges$V2==i)
    X$edges$labels[j] <- labels[i]
    k <- as.integer(labels[i])
    X$pats[j,k] <- T
  }
  X$edges$id <- row.names(X$edges)
  intnode <- NULL
  while(any(is.na(X$edges$labels))){
    tmp <- filter(X$edges, !is.na(X$edges$labels) & !(X$edges$V1 %in% intnode))
    # find nodes directly connected to tips
    a <- unique(tmp$V1[duplicated(tmp$V1)])
    intnode <- c(intnode , a)
    for(i in a){
      s <- tmp$V2[tmp$V1==i]
      X$pats[ which(X$edges$V2==i) ,] <- (colSums(X$pats[X$edges$V2 %in% s,])>=1)

      l <- X$edges$labels[X$edges$V2 %in% s] %>% as.character() %>% strsplit(  split = "_") %>% unlist
      X$edges$labels[ which(X$edges$V2==i)]<- paste(unique(sort(as.numeric(l))),collapse="_")
    }
  }
  X$edges <- rbind(X$edges, data.frame(V1=0,V2=tipnum+1,labels=paste(1:tipnum,collapse="_"),id=nrow(X$edges)+1))
  X$pats <-  rbind(X$pats, rep(T,1,tipnum))

  X$edges$patterns <- NA
  for(i in 1:nrow(X$edges)){
    X$edges$patterns[i] <- as.numeric(X$pats[i,]) %>% paste(collapse="")
  }
  X$edges$samp <- names(samps)[match(X$edges$labels , samps)]
  X
}
