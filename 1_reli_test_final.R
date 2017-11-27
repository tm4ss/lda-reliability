require(lda)
require(topicmodels)
require(htmltools)
require(xtable)
require(igraph)
require(cluster)
source("util/functions_LDA_evaluation.R")
options(stringsAsFactors = F)

require(doMC)
registerDoMC(6)
require(foreach)

# Configurations
tokensToIgnore <- NULL
iterations <- 1000
repetitions <- 10 # repeat experiment three times
alphaPriorsToTest <- rep(0.01, repetitions) 
etaPriorsToTest <- 0.2
initTopicAssignments <- T
coocMethod <- "LL"
clusterMethod <- "PAM"

source("experiments.R")

data_sets <- list(
  list(dir = "data/all_food_de", id = "DE", K = 50),
  list(dir = "data/online_food_us", id = "US", K = 50),
  list(dir = "data/all_food_uk", id = "UK", K = 30)
)

# Create experimental data
for (data_set in data_sets) {
  sampleDirectory <- data_set$dir
  K <- data_set$K
  
  progress_monitor_file <- paste0(data_set$id, "progress.txt")
  write(Sys.time(), progress_monitor_file)
  completeEvalData <- foreach(i = 1:length(experiments)) %dopar% {
    write(paste0("experiment ", i), progress_monitor_file, append = T)
    modelEvaluationData <- runLDAModelEvalution(sampleDirectory, K, iterations = experiments[[i]]$iter, 
                                                alphaPriorsToTest, etaPriorsToTest, blacklist = tokensToIgnore, 
                                                initTopicAssignments = experiments[[i]]$initTopicAssignments,
                                                clusterMethod = clusterMethod,
                                                cooccurrenceMeasure = coocMethod,
                                                nCoocs = 2500)
    save(modelEvaluationData, file = paste0("tmp_eval_data_", data_set$id, "_exp", i, ".RData"))
    return(list(evalData = modelEvaluationData, params = experiments[[i]]))
  }
  write(Sys.time(), progress_monitor_file, append = T)
  save(completeEvalData, file = paste0("final_eval_data_", data_set$id, "_nIterations.RData"))
}


# visualize by running 2_eval_plots.R


