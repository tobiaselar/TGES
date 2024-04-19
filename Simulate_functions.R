library(RBGL)

randtierDAG <- function(incpar, accpar, tino, lB = 0, uB = 1){
  numno <- sum(tino)
  adjmat <- matrix(0,nrow = numno,ncol = numno)
  #tier order
  tior <- c()
  for (o in 1:length(tino)){
    tior <- c(tior,rep(o,tino[o]))
  }
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
  return(adjmat)
  
}

rmvTDAG <- function(n, DAG){
  graphnelDAG <- as(t(DAG), "graphNEL")
  graphnelDAG <- subGraph(tsort(graphnelDAG),graphnelDAG)
  data <- rmvDAG(n,graphnelDAG, errDist = "normal")
  data <- data[,as.character(sort(as.integer(colnames(data))))]
  return(data)
}

ges_to_simple_tges_adj <- function(ges_fit_f,order){
  num_of_nodes <- length(order)
  node.numbers <- 1:num_of_nodes
  Forbidden.edges <-  ges_fit_f$essgraph$.in.edges
  node.names <- ges_fit_f$essgraph$.nodes
  for (nu in node.numbers){
    Forbidden.edges[[nu]] <- node.numbers[order[nu]<order]
  }
  adjmat_return <- addBgKnowledge(createAdjMatrixFromList(ges_fit_f$essgraph$.in.edges)*!createAdjMatrixFromList(Forbidden.edges), checkInput = F)
  return(adjmat_return)  
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
