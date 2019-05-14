bind_single <- function(tree= tree, species = species, sample = c(1), root.age = c()){
	
###add dummy tip
tip <-list(edge=matrix(c(2,1),1,2),
           tip.label= species,
           edge.length=1.0,
           Nnode=1)

class(tip)<-"phylo"
	
tip_tree <- bind.tip(tip, "dummy", edge.length = 0.5)

trees <- tree.bind(tip_tree, tree, sample, root.age)

for(i in 1:length(trees)){

trees[[i]] <- drop.tip(trees[[i]], "dummy")

}

return(trees)	
}



###the jiggle.bind function joins two trees together at a random point between two potential time points

jiggle.bind <- function (x, y, sample, min.age, max.age) 
{
  
  x_class <- mulTree:::check.class(x, c("multiPhylo", "phylo"))
  y_class <- mulTree:::check.class(y, c("multiPhylo", "phylo"))
  if (x_class == "phylo") 
    x <- list(x)
  class(x) <- "multiPhylo"
  if (y_class == "phylo") 
    y <- list(y)
  class(y) <- "multiPhylo"
  if (missing(sample)) {
    sample <- 1
  }
  else {
    mulTree:::check.class(sample, "numeric")
  }
  if (missing(min.age)) {
    min.age <- 0
  }
  else {
    mulTree:::check.class(min.age, "numeric")
  }
  if (missing(max.age)) {
    max.age <- 0
  }
  else {
    mulTree:::check.class(max.age, "numeric")
  }
  rand_x <- mulTree:::sample.trees(x, sample, mulTree:::get.replace(x, sample, 
                                                                    TRUE))
  rand_y <- mulTree:::sample.trees(y, sample, mulTree:::get.replace(y, sample, 
                                                                    TRUE))
  sample_list <- as.list(seq(1:sample))
  age_list <- runif(sample, min.age, max.age)
  
  ########this still selects the first input I think
  lapply.bind.tree.jiggle <- function (element, x, y, rand_x, rand_y, age_list) 
  {
    return(mulTree:::add.root.edge(x[[rand_x[element]]], age_list[[element]]) + mulTree:::add.root.edge(y[[rand_y[element]]], 
                                                                                                        age_list[[element]]))
  }
  
  binded_trees <- lapply(sample_list, lapply.bind.tree.jiggle, x, 
                         y, rand_x, rand_y, age_list)
  if (all(unlist(lapply(binded_trees, Ntip)) == unlist(lapply(binded_trees, 
                                                              function(x) length(unique(x$tip.label)))))) {
    if (length(binded_trees) == 1) {
      binded_trees <- binded_trees[[1]]
      class(binded_trees) <- "phylo"
    }
    else {
      class(binded_trees) <- "multiPhylo"
    }
    return(list(binded_trees, age_list))
  }
  else {
    warning("Some trees have duplicated tip labels.\nThe output can not be converted into phylo or multiPhylo objects.")
    return(list(binded_trees, age_list))
  }
}
