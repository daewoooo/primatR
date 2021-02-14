#' Prepare NJ tree based on data matrix
#'
#' @param data.matrix ...
#' @param boot.iter ...
#' @return A \code{ggplot} object.
#' @importFrom ape nj boot.phylo
#' @author David Porubsky
#' @export
#'
plotDistanceTree <- function(data.matrix, boot.iter=10000) {
  ## Construct a tree
  tree <- ape::nj(X = dist(data.matrix))
  if (boot.iter > 0) {
    boot <- ape::boot.phylo(tree, data.matrix, function(x) nj(dist(x)), B = boot.iter)
    boot <- (boot/boot.iter)*100
    tree$node.label <- boot
  }
  offset <- max(tree$edge.length) + 100
  ## Plot phylogenetic tree
  plt <- ggplot(tree) + 
    geom_tree() + 
    theme_tree2() + 
    geom_tiplab() + 
    geom_nodelab(hjust = -0.1) +
    ggplot2::xlim(0, offset)
  return(plt)
}