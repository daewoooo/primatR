#' Function to test significance of ...
#'
#' @param query.gr A \code{\link{GRanges-class}} object.
#' @param feature.g A \code{\link{GRanges-class}} object.
#' @param nperm ...
#' @param randomize ...
#' @param mask.gr A \code{\link{GRanges-class}} object.
#' @param bsgenome ...
#' @return A \code{list} object ...
#' @author David Porubsky
#' @export
#' 
permDistanceToFeature <- function(query.gr=NULL, feature.gr=NULL, nperm=100, randomize='query', mask.gr=NULL, bsgenome=NULL) {

  observed.dist <- getMinDist(gr = query.gr, userTrack = feature.gr)
  
  rand.dist.l <- list()
  for (iter in 1:nperm) {
    message("Working on permutation: ", iter)
    
    if (randomize == 'query') {
      rand.gr <- randomizeRanges(gr = query.gr[,0], mask.gr = mask.gr, bsgenome = bsgenome)
      rand.dist <- getMinDist(gr = rand.gr, userTrack = feature.gr)
    } else if (randomize == 'feature') {
      rand.gr <- randomizeRanges(gr = feature.gr[,0], mask.gr = mask.gr, bsgenome = bsgenome)
      rand.dist <- getMinDist(gr = rand.gr, userTrack = feature.gr)
    } else {
      message('Please set what ranges to randomize: query or feature?')
    }
    rand.dist.l[[iter]] <- rand.dist
  }
  rand.dist.all <- unlist(rand.dist.l, use.names = FALSE)
  
  ## Stat test
  ttest.p.val <- t.test(observed.dist, rand.dist.all, paired = FALSE, alternative = 'less') #paired false since all ranges are the same
  wilcox.p.val <- wilcox.test(observed.dist, rand.dist.all, paired = FALSE, alternative = 'less')
  enrich.analysis <- list(ttest.p.val = ttest.p.val$p.value,  
                          wilcox.p.val = wilcox.p.val$p.value,
                          observed.dist = observed.dist,
                          rand.dist = rand.dist.all)
  return(enrich.analysis)
}