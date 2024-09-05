# Load packages
library(pcalg)
# Define new classes

## Define new EssGraph class (TEssgraph) with new greedy step that also return new edges added

TEssGraph <- setRefClass("TEssGraph",
                         contains = "EssGraph",
                         methods = list(
                           #' Performs one greedy step
                           greedy.step = function(direction = c("forward", "backward", "turning"), verbose = FALSE, ...) {
                             stopifnot(!is.null(score <- getScore()))
                             
                             ## Cast direction
                             direction <- match.arg(direction)
                             alg.name <- switch(direction,
                                                forward = "GIES-F",
                                                backward = "GIES-B",
                                                turning = "GIES-T")
                             
                             new.graph <- .Call(causalInference,
                                                .in.edges,
                                                score$pp.dat,
                                                alg.name,
                                                score$c.fcn,
                                                causal.inf.options(caching = FALSE,
                                                                   maxSteps = 1, 
                                                                   verbose = verbose,
                                                                   #adaptive = "none", #added by TOBIAS
                                                                   ...))
                             if (identical(new.graph, "interrupt"))
                               return(FALSE)
                             
                             if (new.graph$steps > 0) {
                               last.edges <- .in.edges
                               .in.edges <<- new.graph$in.edges
                               names(.in.edges) <<- .nodes
                               
                               new.in.edges <- .in.edges[sapply(names(.in.edges), function(x) !identical(.in.edges[[x]], last.edges[[x]]))]
                               #returns if any new edges have been added. (Also returns if an edge have been removed and 
                               #there still is an edge going in to this node.)
                             }
                             else
                               new.in.edges <- list()
                             
                             
                             return(c((new.graph$steps == 1),new.in.edges))
                           }), inheritPackage = TRUE
)


## Define score for tges (TBIC)

setRefClass("GaussL0penIntScoreORDER",
            contains = "GaussL0penIntScore",
            
            fields = list( 
              .order = "vector"),
            
            methods = list(
              initialize = function(data = matrix(1, 1, 1),
                                    nodes = colnames(data),
                                    lambda = 0.5*log(nrow(data)),
                                    intercept = TRUE,
                                    format = c("raw", "scatter"),
                                    use.cpp = FALSE, 
                                    order = rep(1,ncol(data)),
                                    ...) {
                .order <<- order
                callSuper(data = data,
                          targets = list(integer(0)),
                          target.index = rep(as.integer(1), nrow(data)),
                          nodes = nodes,
                          lambda = lambda,
                          intercept = intercept,
                          format = format,
                          use.cpp = use.cpp,
                          ...)},#Same as GaussL0penObsScore
              
              
              
              
              #' Calculates the local score of a vertex and its parents
              local.score = function(vertex, parents,...) {
                ## Check validity of arguments
                validate.vertex(vertex)
                validate.parents(parents)
                order <- .order
                if (order[vertex] >= max(c(order[parents],-Inf))){
                  #Checks if the tier of parents are before or same as node
                  
                  if (c.fcn == "none") {
                    ## Calculate score in R
                    if (.format == "raw") {
                      ## calculate score from raw data matrix
                      ## Response vector for linear regression
                      Y <- pp.dat$data[pp.dat$non.int[[vertex]], vertex]
                      sigma2 <- sum(Y^2)
                      
                      if (length(parents) + pp.dat$intercept != 0) {
                        ## Get data matrix on which linear regression is based
                        Z <- pp.dat$data[pp.dat$non.int[[vertex]], parents, drop = FALSE]
                        if (pp.dat$intercept)
                          Z <- cbind(1, Z)
                        
                        ## Calculate the scaled error covariance using QR decomposition
                        Q <- qr.Q(qr(Z))
                        sigma2 <- sigma2 - sum((Y %*% Q)^2)
                      }
                    }
                    else if (.format == "scatter") {
                      ## Calculate the score based on pre-calculated scatter matrices
                      ## If an intercept is allowed, add a fake parent node
                      parents <- sort(parents)
                      if (pp.dat$intercept)
                        parents <- c(pp.dat$vertex.count + 1, parents)
                      
                      pd.scMat <- pp.dat$scatter[[pp.dat$scatter.index[vertex]]]
                      sigma2 <- pd.scMat[vertex, vertex]
                      if (length(parents) != 0) {
                        b <- pd.scMat[vertex, parents]
                        sigma2 <- sigma2 - as.numeric(b %*% solve(pd.scMat[parents, parents], b))
                      }
                    }
                    
                    ## Return local score
                    lscore <- -0.5*pp.dat$data.count[vertex]*(1 + log(sigma2/pp.dat$data.count[vertex])) -
                      pp.dat$lambda*(1 + length(parents))
                    #print(c(lscore,"v",vertex,"p",parents))
                    return(lscore)
                  } else {
                    ## Calculate score with the C++ library (NOT ABLE TO DO THIS YET)
                    #.Call(localScore, c.fcn, pp.dat, vertex, parents, c.fcn.options(...))
                    stop("Not able to compute using C++. Set use.cpp = FALSE")
                  } 
                }
                else { skip <- -Inf
                #print(c(skip,"v",vertex,"p",parents))
                return(skip)}#set score to minus infinity if vertex earlier than parents
              }),
            inheritPackage = TRUE
)


# Define new functions

## transform list of type .in.edges to adjancency matrix

createAdjMatrixFromList <- function(inputList) {
  # Determine the number of rows and columns (n x n matrix)
  n <- length(inputList)
  
  # Initialize a matrix of 0s
  resultMatrix <- matrix(0, nrow = n, ncol = n)
  
  # Iterate over the list to update the matrix
  for (i in seq_along(inputList)) {
    indices <- inputList[[i]]
    # Ensure indices are within the bounds of the matrix
    validIndices <- indices[indices <= n]
    if (length(validIndices) > 0) {
      resultMatrix[i, validIndices] <- 1
    }
  }
  rownames(resultMatrix) <- names(inputList)
  colnames(resultMatrix) <- names(inputList)
  return(resultMatrix)
}

## Adjancency matrix to .in.edges list for an essgraph

createListFromAdjMatrix <- function(adjMatrix) {
  # Get the number of rows (or columns, since it's n x n) in the adjacency matrix
  n <- nrow(adjMatrix)
  
  # Initialize an empty list
  resultList <- vector("list",n)
  names(resultList) <- rownames(adjMatrix)
  
  # Iterate over the rows of the matrix
  for (i in 1:n) {
    # Find the indices (column numbers) where there is an in edge (value == 1)
    connectedIndices <- as.integer(which(adjMatrix[i, ] == 1))
    
    
    if (length(connectedIndices) > 0) {
      resultList[[i]] <- connectedIndices
    } else {
      # If no connections, assign NULL or an empty vector
      resultList[[i]] <- integer(0)
    }
  }
  
  return(resultList)
}

tges <- function(score, verbose = FALSE){
  node.numbers <- 1:score$pp.dat$vertex.count
  essgraph <- new("TEssGraph", nodes = as.character(node.numbers), score = score)
  Forbidden.edges <- essgraph$.in.edges #list of nodes all with integer(0) entry
  node.names <- score$.nodes
  num.bidir <- 0
  num.directed <- 0
  order <- score$.order
  for (n in node.numbers){
    Forbidden.edges[[n]] <- node.numbers[order[n]<order]
  }
  
  cont <- TRUE
  while(cont) {
    cont <- FALSE
    
    #Forward phase
    runwhile <- TRUE
    while(runwhile){
      tempstep <- essgraph$greedy.step("forward", verbose = verbose)
      runwhile <- as.logical(tempstep[1])
      if (runwhile){cont <- TRUE}
      else break
      
      for (i in names(tempstep[-1])){ #Run through the nodes that have been changed
        in.node.edges <- tempstep[-1][[i]] #save the in.node edges of node i
        forbidden.node.edges <- Forbidden.edges[[as.numeric(i)]]
        removed.edges <- in.node.edges[in.node.edges %in% forbidden.node.edges] #List of edges to be removed
        if (length(removed.edges) > 0){
          bgx <- rep(as.numeric(i),length(removed.edges))
          bgy <- removed.edges
          amatbg <- addBgKnowledge(gInput = createAdjMatrixFromList(essgraph$.in.edges), x = bgx, y = bgy,verbose = verbose)
          no.forbidden.edges <- createListFromAdjMatrix(amatbg)
          essgraph$.in.edges <- no.forbidden.edges
        }
      }
    }
    
    #Backward phase
    runwhile <- TRUE
    while(runwhile){
      tempstep <- essgraph$greedy.step("backward", verbose = verbose)
      runwhile <- as.logical(tempstep[1]) 
      if (runwhile){cont <- TRUE}
      else break
      
      for (i in names(tempstep[-1])){ #Run through the nodes that have been changed
        in.node.edges <- tempstep[-1][[i]] #save the in.node edges of node i
        forbidden.node.edges <- Forbidden.edges[[as.numeric(i)]]
        removed.edges <- in.node.edges[in.node.edges %in% forbidden.node.edges] #List of edges to be removed
        if (length(removed.edges) > 0){
          bgx <- rep(as.numeric(i),length(removed.edges))
          bgy <- removed.edges
          amatbg <- addBgKnowledge(gInput = createAdjMatrixFromList(essgraph$.in.edges), x = bgx, y = bgy, verbose = verbose)
          no.forbidden.edges <- createListFromAdjMatrix(amatbg)
          essgraph$.in.edges <- no.forbidden.edges
        }
      }
    }
    
    #Turning phase
    runwhile <- TRUE
    while(runwhile){
      tempstep <- essgraph$greedy.step("turning", verbose = verbose)
      runwhile <- as.logical(tempstep[1]) 
      if (runwhile){cont <- TRUE}
      else break
      
      for (i in names(tempstep[-1])){ #Run through the nodes that have been changed
        in.node.edges <- tempstep[-1][[i]] #save the in.node edges of node i
        forbidden.node.edges <- Forbidden.edges[[as.numeric(i)]]
        removed.edges <- in.node.edges[in.node.edges %in% forbidden.node.edges] #List of edges to be removed
        if (length(removed.edges) > 0){
          bgx <- rep(as.numeric(i),length(removed.edges))
          bgy <- removed.edges
          amatbg <- addBgKnowledge(gInput = createAdjMatrixFromList(essgraph$.in.edges), x = bgx, y = bgy, verbose = verbose)
          no.forbidden.edges <- createListFromAdjMatrix(amatbg)
          essgraph$.in.edges <- no.forbidden.edges
        }
        
      }
    }
  }
  essgraph$.nodes <- node.names # Save names of nodes
  names(essgraph$.in.edges) <- node.names # Save names of nodes
  if (verbose){cat("Number of edges directed", num.bidir,"\nNumber of directed edges removed", num.directed, "\n")}
  return(essgraph)
}