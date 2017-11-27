require(ggplot2)
source("util/functions_LDA_evaluation.R")

get_final_plot_data <- function(completeEvalData, dataset_id) {
  
  completeReliabilty <- NULL
  completeCoherence <- NULL
  for (i in 1:length(completeEvalData)) {
    modelEvaluationData <- completeEvalData[[i]]$evalData
    
    current_res_r <- data.frame(
      reliability = getModelAlignments(modelEvaluationData, topWordsToMatch = 20, similarityThreshold = 0.3),
      strategy = completeEvalData[[i]]$params$initTopicAssignments, 
      iteration = completeEvalData[[i]]$params$iter)
    completeReliabilty <- rbind(completeReliabilty, current_res_r)
    
    current_res_c <- data.frame(
      coherence = modelEvaluationData$modelCoherence,
      strategy = completeEvalData[[i]]$params$initTopicAssignments, 
      iteration = completeEvalData[[i]]$params$iter)
    completeCoherence <- rbind(completeCoherence, current_res_c)
  }
  
  print(ggplot(completeReliabilty, aes(iteration, reliability, group = strategy, color = strategy)) +
    stat_summary(fun.data = mean_se, geom = "errorbar") +
    stat_summary(fun.y = mean, geom = "line", size = 1))
  
  
  print(ggplot(completeCoherence, aes(iteration, coherence, group = strategy, color = strategy)) +
    stat_summary(fun.data = mean_se, geom = "errorbar") +
    stat_summary(fun.y = mean, geom = "line", size = 1))
  
  completeReliabilty <- cbind(completeReliabilty, dataset = dataset_id)
  completeCoherence <- cbind(completeCoherence,  dataset = dataset_id) 
  
  return(list(completeReliabilty = completeReliabilty, completeCoherence = completeCoherence))
}



# --------------
load("final_eval_data_DE_nIterations.RData")
de_res <- get_final_plot_data(completeEvalData, "DE, K = 50")

load("final_eval_data_US_nIterations.RData")
us_res <- get_final_plot_data(completeEvalData, "US, K = 50")

load("final_eval_data_UK_nIterations.RData")
uk_res <- get_final_plot_data(completeEvalData, "UK, K = 30")


# --------------------------------------
final_reli <- rbind(us_res$completeReliabilty, de_res$completeReliabilty, uk_res$completeReliabilty)
final_cohe <- rbind(us_res$completeCoherence, de_res$completeCoherence, uk_res$completeCoherence)

final_reli$strategy <- factor(final_reli$strategy)
levels(final_reli$strategy) <- c("RANDOM", "SEED", "CLUSTER")

final_cohe$strategy <- factor(final_cohe$strategy)
levels(final_cohe$strategy) <- c("RANDOM", "SEED", "CLUSTER")

# save(final_reli, final_cohe, file = "final_plot_data.RData")
# load("final_plot_data.RData")

ggplot(final_reli, aes(iteration, reliability, group = strategy, color = strategy)) +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_summary(fun.y = mean, geom = "line") +
  facet_wrap(~dataset) + theme(legend.position="bottom")


ggplot(final_cohe, aes(iteration, coherence, group = strategy, color = strategy)) +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_summary(fun.y = mean, geom = "line") +
  facet_wrap(~dataset) + theme(legend.position="bottom")

# Without 1st iteration
df <- final_cohe[final_cohe$iteration > 1, ]
ggplot(df, aes(iteration, coherence, group = strategy, color = strategy)) +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_summary(fun.y = mean, geom = "line") +
  facet_wrap(~dataset) + theme(legend.position="bottom")
