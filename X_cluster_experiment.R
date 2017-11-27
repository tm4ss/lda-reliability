# CLUSTERING EXPERIMENT
# ---------------------
# Which clustering method performs best with which 
# cooccurrence measure?

require(lda)
require(topicmodels)
require(htmltools)
require(xtable)
require(igraph)
require(cluster)
source("functions_LDA_evaluation.R")
options(stringsAsFactors = F)

# ================================================
sampleName <- "online_climate_uk"
# ================================================
sampleDirectory <- paste0("../", sampleName)


# ------------------------------------------------
blacklistClimate <- readLines("backlist_climate.txt", encoding = "UTF-8")
blacklistBadTokens <- readLines("backlist_badtokens.txt", encoding = "UTF-8")
tokensToIgnore <- union(blacklistClimate, blacklistBadTokens)

# -----------------------------------------------
# EXPERIMENTS ON CLUSTERINGS
K <- 30
iterations <- 200
repetitions <- 2 # repeat experiment three times
alphaPriorsToTest <- rep(0.002, repetitions) 
etaPriorsToTest <- c(1 / K)

coocMethods <- c("DICE", "LL", "POISSON")
clusterMethods <- c("INFOMAP", "KMST", "PAM")

reliabilityMatrix <- matrix(0, nrow = length(coocMethods), ncol = length(clusterMethods), dimnames = list(coocMethods, clusterMethods))
coherenceMatrix <- matrix(0, nrow = length(coocMethods), ncol = length(clusterMethods), dimnames = list(coocMethods, clusterMethods))
completeEvalData <- vector("list", length(coocMethods) * length(clusterMethods))

i <- 0
for (coocMethod in coocMethods) {
  for (clusterMethod in clusterMethods) {
    i <- i + 1
    modelEvaluationData <- runLDAModelEvalution(sampleDirectory, K, iterations, alphaPriorsToTest, etaPriorsToTest, blacklist = tokensToIgnore, 
                                                initTopicAssignments = T,
                                                clusterMethod = clusterMethod,
                                                cooccurrenceMeasure = coocMethod,
                                                nCoocs = 1500)
    completeEvalData[[i]] <- list(evalData = modelEvaluationData, params = paste0(coocMethod, " ", clusterMethod))
    reliabilityMatrix[coocMethod, clusterMethod] <- getModelAlignments(modelEvaluationData)
    coherenceMatrix[coocMethod, clusterMethod] <- mean(modelEvaluationData$modelCoherence)
  }
}
save(completeEvalData, file=paste0("RESULT_cluster_experiment_", sampleName, "_K", K, "_i", iterations, ".RData"))
print(reliabilityMatrix)
print(coherenceMatrix)



# GET SAVED RESULTS
# -----------------
load("RESULT_cluster_experiment_online_climate_de_K30_i500.RData")
i <- 0
reliabilityMatrix <- matrix(0, nrow = length(coocMethods), ncol = length(clusterMethods), dimnames = list(coocMethods, clusterMethods))
coherenceMatrix <- matrix(0, nrow = length(coocMethods), ncol = length(clusterMethods), dimnames = list(coocMethods, clusterMethods))
for (coocMethod in coocMethods) {
  for (clusterMethod in clusterMethods) {
    i <- i + 1
    reliabilityMatrix[coocMethod, clusterMethod] <- getModelAlignments(completeEvalData[[i]]$evalData, verbose = T)
    coherenceMatrix[coocMethod, clusterMethod] <- mean(completeEvalData[[i]]$evalData$modelCoherence)
  }
}
print(reliabilityMatrix)
print(coherenceMatrix)




# PARALLELIZATION
# --------------------------------
K <- 30
iterations <- 1000
repetitions <- 10 # repeat experiment three times
alphaPriorsToTest <- rep(0.002, repetitions) 
etaPriorsToTest <- c(1 / K)

coocMethods <- c("DICE", "LL", "POISSON")
clusterMethods <- c("INFOMAP", "KMST", "PAM")

reliabilityMatrix <- matrix(0, nrow = length(coocMethods), ncol = length(clusterMethods), dimnames = list(coocMethods, clusterMethods))
coherenceMatrix <- matrix(0, nrow = length(coocMethods), ncol = length(clusterMethods), dimnames = list(coocMethods, clusterMethods))
completeEvalData <- vector("list", length(coocMethods) * length(clusterMethods))

require(parallel)
cl <- makeCluster(9, methods = T, outfile = "")
clusterExport(cl, c("runLDAModelEvalution", "topicCoherence", "ldaformat2Matrix", "getModelAlignments"))


i <- 0
for (coocMethod in coocMethods) {
  for (clusterMethod in clusterMethods) {
    i <- i + 1
    completeEvalData[[i]] <- list(evalData = NULL, params = list(coocMethod = coocMethod, clusterMethod = clusterMethod))
  }
}

runLDAModelEvalutionPara <- function(i, paramList, sampleDirectory, K, iterations, alphaPriorsToTest, etaPriorsToTest, blacklist, 
                                     initTopicAssignments, nCoocs) {
  runLDAModelEvalution(sampleDirectory, K, iterations, alphaPriorsToTest, etaPriorsToTest, blacklist = blacklist, 
                       initTopicAssignments = initTopicAssignments,
                       clusterMethod = paramList[[i]]$params$clusterMethod,
                       cooccurrenceMeasure = paramList[[i]]$params$coocMethod,
                       nCoocs = nCoocs)
}

# MAP: Parallel computation of parameter modes
parallelEvalData <- parSapply(cl, 1:length(completeEvalData), runLDAModelEvalutionPara, paramList = completeEvalData,
                              sampleDirectory = sampleDirectory, K = K, iterations = iterations, 
                              alphaPriorsToTest = alphaPriorsToTest, etaPriorsToTest = etaPriorsToTest, 
                              blacklist = tokensToIgnore, initTopicAssignments = T, nCoocs = 3000)

for (i in 1:ncol(parallelEvalData)) {
  completeEvalData[[i]]$evalData = as.data.frame(parallelEvalData[, i])
}

stopCluster(cl)


# REDUCE: Combine results from single map processes
i <- 0
reliabilityMatrix <- matrix(0, nrow = length(coocMethods), ncol = length(clusterMethods), dimnames = list(coocMethods, clusterMethods))
coherenceMatrix <- matrix(0, nrow = length(coocMethods), ncol = length(clusterMethods), dimnames = list(coocMethods, clusterMethods))
for (coocMethod in coocMethods) {
  for (clusterMethod in clusterMethods) {
    i <- i + 1
    reliabilityMatrix[coocMethod, clusterMethod] <- mean(getModelAlignments(completeEvalData[[i]]$evalData, verbose = T))
    coherenceMatrix[coocMethod, clusterMethod] <- mean(completeEvalData[[i]]$evalData$modelCoherence)
  }
}
save(completeEvalData, file=paste0("RESULT_cluster_experiment_", sampleName, "_K", K, "_i", iterations, "_rep", repetitions, ".RData"))
print(reliabilityMatrix)
print(coherenceMatrix)




# Get statistical significance of difference 
# ------------------------------------------

#load("RESULT_cluster_experiment_online_climate_uk_K30_i1000_rep10.RData")
#load("RESULT_cluster_experiment_online_climate_de_K30_i1000_rep10.RData")
#load("RESULT_cluster_experiment_online_food_uk_K30_i1000_rep10.RData")
load("RESULT_cluster_experiment_online_food_de_K30_i1000_rep10.RData")


# RELIABILITY
worstResultId <- 9
reliabilitiesOfWorstModel <- getModelAlignments(completeEvalData[[worstResultId]]$evalData, verbose = F)
i <- 0
pMatrix <- matrix(0, nrow = length(coocMethods), ncol = length(clusterMethods), dimnames = list(coocMethods, clusterMethods))
improvementMatrix  <- matrix(0, nrow = length(coocMethods), ncol = length(clusterMethods), dimnames = list(coocMethods, clusterMethods))
for (coocMethod in coocMethods) {
  for (clusterMethod in clusterMethods) {
    i <- i + 1
    allReliabilities <- getModelAlignments(completeEvalData[[i]]$evalData, verbose = F)
    pMatrix[coocMethod, clusterMethod] <- t.test(reliabilitiesOfWorstModel, allReliabilities, alternative = "greater")$p.value
    improvementMatrix[coocMethod, clusterMethod] <- mean(allReliabilities) / mean(reliabilitiesOfWorstModel) - 1
  }
}
print(improvementMatrix)
print(pMatrix)


# COHERENCE
#worstResultId <- which.min(as.vector(t(coherenceMatrix)))
worstResultId <- 7
i <- 0
pMatrix <- matrix(0, nrow = length(coocMethods), ncol = length(clusterMethods), dimnames = list(coocMethods, clusterMethods))
improvementMatrix  <- matrix(0, nrow = length(coocMethods), ncol = length(clusterMethods), dimnames = list(coocMethods, clusterMethods))
for (coocMethod in coocMethods) {
  for (clusterMethod in clusterMethods) {
    i <- i + 1
    pMatrix[coocMethod, clusterMethod] <- t.test(completeEvalData[[worstResultId]]$evalData$modelCoherence, completeEvalData[[i]]$evalData$modelCoherence, alternative = "less")$p.value
    improvementMatrix[coocMethod, clusterMethod] <- 1 - mean(completeEvalData[[i]]$evalData$modelCoherence) / mean(completeEvalData[[worstResultId]]$evalData$modelCoherence)
  }
}
print(improvementMatrix)
print(pMatrix)