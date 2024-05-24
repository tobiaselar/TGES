source("Simulate_functions.R") #Import functions needed for simulationg tiered DAGs

# Simulate graph

## Choose parameters
n <- 10000
lB <- 0.5
uB <- 1
incpar <- 0.33
accpar <- 0.5
num_no <- 5
num_ti <- 3

set.seed(123)
tiered_ordering <- sort(c(sample(c(1:num_ti),(num_no-num_ti),replace = T),1:num_ti))
tino <- c()
for (o in 1:num_ti){
  tino <- c(tino,sum(tiered_ordering == o))
}
# Simulate DAG from chosen tiers, nodes and density parameters
weightadjmat <- randtierDAG(incpar,accpar,tino, lB = lB, uB = uB)

# Introduce a v-structure
#weightadjmat <- weightadjmat*matrix(c(rep(1,4),0,rep(1,20)), nrow = 5) 

# Generate data from DAG
data <- rmvTDAG(n, weightadjmat)

# Save adjmat for DAG, CPDAG an tiered MPDAG
adjmat <- ceiling(weightadjmat/uB)
DAGadjmat <- adjmat
CPDAG <- dag2cpdag(as(t(adjmat), "graphNEL"))
CPDAGadjmat <- t(as(CPDAG,"matrix"))
adjmat <- dag2tmpdag(adjmatrix = adjmat,odr = tiered_ordering)
#plot(CPDAG)

# Fit TGES
source("TGES_function.R") #Import function and class needed
t_score <- new("GaussL0penIntScoreORDER", order = tiered_ordering
               , data = data 
               , use.cpp = F)
tges_fit <- tges(score = t_score, order = tiered_ordering)
TGES_est_adjmat <- t(as(tges_fit, "matrix"))

plot(tges_fit)







# Fit GES
ges_score <- new("GaussL0penIntScore"
                 , data = data )
ges_fit <- ges(score = ges_score)
GES_est_adjmat <- createAdjMatrixFromList(ges_fit$essgraph$.in.edges)

plot(ges_fit$essgraph)

# Fit Simple TGES

SimpleTGES_est_adjmat <- ges_to_simple_tges_adj(ges_fit_f = ges_fit,ord_f = tiered_ordering)

#Fit TPC
tpc_fit <- tpc::tpc(suffStat = list(C = cor(data), n = n),
                    indepTest = gaussCItest, alpha = 0.01, tiers = tiered_ordering, p = num_no)
TPC_est_adjmat <- as(tpc_fit,"matrix")
plot(tpc_fit@graph)

# TGES fully connected cross tier
#tges_fully_fit <- tges_fully(t_score,tiered_ordering,verbose = F)
#plot(tges_fully_fit)
