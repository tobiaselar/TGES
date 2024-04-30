library(RBGL)
library(causalDisco)

randtierDAG <- function(incpar, accpar, tino, lB = 0, uB = 1){
  numno <- sum(tino)
  adjmat <- matrix(0,nrow = numno,ncol = numno)
  #tier order
  tior <- c()
  for (o in 1:length(tino)){
    tior <- c(tior,rep(o,tino[o]))
  }
  
  emptyDAG <- T
  while (emptyDAG) {
      #In tier DAGs
    for (i in 1:length(tino)){
      if (tino[i] > 1){
        tempDAG <- r.gauss.pardag(p = tino[i], prob = incpar, lbe = lB, ube = uB, neg.coef = F)
        adjmat[which(tior == i),which(tior == i)] <- t(tempDAG$weight.mat())
        #createAdjMatrixFromList(tempDAG$.in.edges) 
      }
    }
    # Make across tier edges
    for (n in 1:(length(tino)-1)){
      plustier <- which(tior > n)
      tier <- which(tior == n)
      adjmat[plustier,tier] <- rbinom(n = length(plustier)*length(tier), size = 1,prob = accpar)*runif(length(plustier)*length(tier), min = lB, max = uB)
    }
    emptyDAG <- (sum(adjmat) == 0)
    
  }
  
  return(adjmat)
  
}

rmvTDAG <- function(n, DAG){
  graphnelDAG <- as(t(DAG), "graphNEL")
  graphnelDAG <- subGraph(tsort(graphnelDAG),graphnelDAG)
  data <- rmvDAG(n,graphnelDAG, errDist = "normal")
  data <- data[,as.character(sort(as.integer(colnames(data))))]
  return(data)
}

ges_to_simple_tges_adj <- function(ges_fit_f,ord_f){
  cpadjmat <- t(as(ges_fit_f$essgraph,"matrix")) #t(as(dag2cpdag(ges_fit_f$repr), "matrix")) 
  num_of_nodes <- length(ord_f)
  node.numbers <- 1:num_of_nodes
  Forbidden.edges.matrix <- matrix(0,ncol = num_of_nodes, nrow = num_of_nodes)
  for (nod in node.numbers){
    Forbidden.edges.matrix[nod,] <- ord_f[nod]<ord_f
  }
  adjmat_return <- addBgKnowledge(cpadjmat*!Forbidden.edges.matrix, checkInput = F)
  return(adjmat_return)  
}

# Overwrite causalDisco confusion "directional" function
dir_confusion <- function(est_amat, true_amat) {
  est_edges <- causalDisco::edges(est_amat)
  true_edges <- causalDisco::edges(true_amat)
  
  true_adj <- c(true_edges$undir, true_edges$dir)
  true_adj <- c(true_adj, lapply(true_adj, rev))
  
  true_dir <- true_edges$dir
  true_revdir <- lapply(true_dir, rev)
  true_undir <- true_edges$undir
  true_undir <- c(true_undir, lapply(true_undir, rev)) #Added to fix bug
  
  est_dir <- est_edges$dir
  est_undir <- est_edges$undir
  
  dir_fp <- 0
  dir_fn <- 0
  dir_tp <- 0
  dir_tn <- 0
  
  #count metrics for undirected edges
  if (length(est_undir) > 0) {
    for (i in 1:length(est_undir)) {
      thisedge <- est_undir[i]
      if (thisedge %in% true_adj) {
        if (thisedge %in% true_undir) { #is correctly undirected
          dir_tn <- dir_tn + 1
        } else if (thisedge %in% c(true_dir, true_revdir)) { #is undirected, should be directed
          dir_fn <- dir_fn + 1
        }
      }
    }
  }
  
  #count metrics for directed edges
  if (length(est_dir) > 0) {
    for (i in 1:length(est_dir)) {
      thisedge <- est_dir[i]
      if (thisedge %in% true_adj) {
        if (thisedge %in% true_undir) { #is directed, should be undirected
          dir_fp <- dir_fp + 1 
        } else if (thisedge %in% true_dir) { #is directed in correct direction
          dir_tp <- dir_tp + 1
        } 
        if (thisedge %in% true_revdir) { #is directed in incorrect direction
          dir_fp <- dir_fp + 1
          dir_fn <- dir_fn + 1
        }
      }
    }
  }
  list(tp = dir_tp, tn = dir_tn,
       fp = dir_fp, fn = dir_fn)
}

confusion <- function(est_amat, true_amat, type = "adj") {
  #UseMethod("confusion")
  
  est_class <- class(est_amat)
  true_class <- class(true_amat)
  
  if (any(est_class %in% c("tpdag", "cpdag"))) { 
    est_amat <- amat(est_amat)
  }
  
  if (any(true_class %in% c("tpdag", "cpdag"))) {
    true_amat <- amat(true_amat)
  }
  
  if (type == "adj") {
    adj_confusion(est_amat, true_amat)
  } else if (type == "dir") {
    dir_confusion(est_amat, true_amat)
  } else {
    stop("Type must be either adj or dir.")
  }
}


intier_confusion <- function(est_amat,true_amat,type = "adj",tierorder){
  conf <- confusion(1,1) #Confusion matrix with all zeroes
  for (clu in unique(tierorder)){
    amat_index <- which(tierorder == clu)
    if (length(amat_index) > 1){
      conf <- mapply("+",confusion(est_amat = est_amat[amat_index,amat_index],
                                   true_amat = true_amat[amat_index,amat_index],
                                   type = type) ,
                     conf , SIMPLIFY = FALSE)
    }
    
  }
  if (all(conf == 0)){
    conf$tp <- NA
  }
  return(conf)
}

dag2tmpdag <- function(adjmatrix,odr){
  cpdag <- dag2cpdag(as(t(adjmatrix), "graphNEL"))
  cpadjmat <- t(as(cpdag,"matrix"))
  number_nodes <- length(odr)
  node_numbers <- 1:number_nodes
  Forbidden.edges.matrix <- matrix(0,ncol = number_nodes, nrow = number_nodes)  #list of nodes all with integer(0) entry
  for (nod in node_numbers){
    Forbidden.edges.matrix[nod,] <- odr[nod]<odr
  }
  tmadjmat <- addBgKnowledge(cpadjmat*!Forbidden.edges.matrix, checkInput = F)
  return(tmadjmat)
}
