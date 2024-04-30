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
                                                causal.inf.options(caching = FALSE, #Why not true (TOBIAS NOTE)
                                                                   maxSteps = 1, #Could be more? (Tobias NOTE). If no maxSteps algorithm continue untill no further improvement (no need for while loop)
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
                               
                               #new.in.edges <- .in.edges[!(.in.edges %in% last.edges)]
                               #returns if any new edges have been added. (Also returns if a edge have been removed and theres still an edge going in to this note. Can be fixed by not doing it for "backward")
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
            
            #validity = function(object) {
            #Define something to check validity of "ORDER" TODO
            #},
            
            methods = list(
              #' Constructor
              initialize = function(data = matrix(1, 1, 1),
                                    nodes = colnames(data),
                                    lambda = 0.5*log(nrow(data)),
                                    intercept = TRUE,
                                    format = c("raw", "scatter"),
                                    use.cpp = TRUE,
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
                  #Checks if parents are before or same
                  
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
                    #return(0)
                  } else {
                    ## Calculate score with the C++ library
                    .Call(localScore, c.fcn, pp.dat, vertex, parents, c.fcn.options(...))
                  } # IF c.fcn
                }
                else { skip <- -Inf
                #print(c(skip,"v",vertex,"p",parents))
                return(skip)}#set score to minus infinity if vertex earlier than parents
              }),
            inheritPackage = TRUE
)