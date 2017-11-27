# Authors: 
# - Gregor Wiedemann, gregor.wiedemann@uni-leipzig.de
# - Andreas Niekler, aniekler@informatik.uni-leipzig.de
# ASV, Universität Leipzig
# Leipzig, November 2017

runLDAModelEvalution <- function(
  sampleDirectory, 
  K, 
  iterations, 
  alphaPriorsToTest, 
  etaPriorsToTest, 
  blacklist = NULL, 
  initTopicAssignments = T,
  clusterMethod = "PAM",
  cooccurrenceMeasure = "LL",
  nCoocs = 2500,
  seed = 9127
) {
  require(xtable)
  require(Matrix)
  require(lda)
  require(topicmodels)
  # Read input data
  corpusData <- read.documents(filename = paste0(sampleDirectory, "/corpus/corpus"))
  corpusVocab <- readLines(paste0(sampleDirectory, "/corpus/vocabulary.dat"), encoding = "UTF-8")
  
  if (!is.null(blacklist)) {
    corpusDTM <- ldaformat2dtm(corpusData, corpusVocab)
    corpusDTM <- corpusDTM[, !(colnames(corpusDTM) %in% blacklist)]
    corpusLDA <- dtm2ldaformat(corpusDTM)
    corpusData <- corpusLDA$documents
    corpusVocab <- corpusLDA$vocab
  }
  
  if (initTopicAssignments == "SEED") {
    
    # set fixed random seed for initialization of neglected terms in coccurrence computation
    set.seed(seed)
    vocabClusterAssignments <- sample.int(K, size = length(corpusVocab), replace=TRUE) - 1
    # set seed to random system time again
    st <- as.numeric(Sys.time()); set.seed((st - floor(st)) * 1e8)
    initAssignments <- lapply(corpusData, FUN = function(x) {vocabClusterAssignments[(x[1, ] + 1)]})
    initAssignments <- lapply(initAssignments, as.integer) 
    initAssignments <- list(assignments = initAssignments)
    initStrategy <- "seed"
    
  } else if (initTopicAssignments == TRUE) {
    source("util/topicMapping.R")
    initAssignments <- computeInitialAssignments(corpusData, corpusVocab, K, nCoocs = nCoocs, clusterMethod = clusterMethod, cooccurrenceMeasure = cooccurrenceMeasure)
    initAssignments <- list(assignments = initAssignments)
    initStrategy <- "fix"
    #t <- lda.collapsed.gibbs.sampler(corpusData, K, vocab = corpusVocab, num.iterations = iterations, alpha = alpha, eta = eta, trace = 1L, compute.log.likelihood = T)
    #s <- data.frame(sapply(initAssignments, length), sapply(t$assignments, length))
    #s[,1] == s[,2]
  } else if (initTopicAssignments == FALSE) {
    initAssignments <- NULL
    initStrategy <- "rnd"
  }

  # Prepare output data
  modelEvaluationData <- data.frame(modelID = integer(), K = integer(), alpha = double(), eta = double(), modelLikelihood = double(), modelCoherence = integer(), topicsHTML = character(), modelFileName = character())
  modelEvaluationHTML <- data.frame()
  # Run evaluation
  modelID <- 0
  for (alpha in alphaPriorsToTest) {
    for (eta in etaPriorsToTest) {
      modelID <- modelID + 1
      topicModel <- lda.collapsed.gibbs.sampler(corpusData, K, vocab = corpusVocab, num.iterations = iterations, alpha = alpha, eta = eta, initial = initAssignments, trace = 1L, compute.log.likelihood = T)
      topicTerms <- top.topic.words(topicModel$topics, 25, by.score=TRUE)
      topicProportions <- t(topicModel$document_sums) / colSums(topicModel$document_sums)
      modelLikelihood <- tail(as.vector(topicModel$log.likelihoods[2, ]), 1)
      topicCoherenceForAllTopics <- topicCoherence(ldaformat2Matrix(corpusData, corpusVocab), topicModel)
      modelCoherence <- mean(topicCoherenceForAllTopics)
      
      # Prob
      tProportions <- colSums(topicProportions) / nrow(topicProportions)
      oProportions <- order(tProportions, decreasing = T)
      # Rank 1
      firstDocTopics <- apply(topicProportions, 1, FUN=function(x) order(x, decreasing=TRUE)[1])
      primaryDocTopics <- factor(firstDocTopics, 1:K)
      nRanks1 <- table(primaryDocTopics)
      tRanks1 <- as.integer(nRanks1)
      oRanks1 <- order(tRanks1, decreasing = T)
      # Coherence
      tCoherence <- topicCoherenceForAllTopics
      oCoherence <- order(topicCoherenceForAllTopics, decreasing = T)
      topics <- data.frame(
        TopicID = sprintf("%02d", 1:K),
        Rank1 = tRanks1,
        Prob = tProportions,
        Coherence = tCoherence,
        Terms = paste0("<pre>", apply(topicTerms, 2, paste0, collapse=" "), "</pre>")
      )
      topics <- topics[oRanks1, ]
      topicsHTML <- print(xtable(topics, digits=5), print.results = F, type="html", sanitize.text.function=function(x){x}, include.rownames = F, html.table.attributes = paste0('id="T', modelID, '" class="sortable"'))

      modelFileName <- paste0(sampleDirectory, "/model-", modelID, 
                              "_K", K, "_a", alpha, "_e", eta, "_i", iterations, initStrategy,
                              "_d", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), ".RData")
      save(topicModel, file = modelFileName)
      
      modelEvaluationHTML <- rbind(modelEvaluationHTML, data.frame(modelID, K, alpha, eta, modelLikelihood, modelCoherence, topicsHTML, modelFileName))
      print(paste0("LDA Model ", modelID, " has been saved to ", modelFileName))
    }
  }
  evalFileName <- paste0(sampleDirectory, "/modelEvaluationData_K", K, "_i", iterations, initStrategy, ".html")
  file.copy("html_template.html", evalFileName, overwrite = T)
  print(xtable(modelEvaluationHTML, digits=4), type="html", file=evalFileName, append=T, sanitize.text.function=function(x){x}, include.rownames = F)
  
  print(paste0("Final LDA model evaluation CSV data has been has been saved to ", modelFileName))
  
  return(modelEvaluationHTML)
}




# Berechnet topic coehrence  nach(Mimno et al. 2011) 
# - params: Dokument-Term-Matrix
# - LDA Modell
# - Anzahl an höchstwahrscheinlichen Termen zur Berechnung der Dokument-Kookkurrenz

topicCoherence <- function(DTM, ldaModel, N = 25) {
  
  # Ensure matrix or Matrix-format (convert if SparseM)
  require(Matrix)
  require(slam)
  if (is.simple_triplet_matrix(DTM)) {
    DTM <- sparseMatrix(i=DTM$i, j=DTM$j, x=DTM$v, dims=c(DTM$nrow, DTM$ncol), dimnames = dimnames(DTM))
  }
  
  DTMBIN <- DTM
  DTMBIN[DTMBIN > 0] <- 1
  
  documentFrequency <- colSums(DTMBIN)
  names(documentFrequency) <- colnames(DTMBIN)
  
  K <- nrow(ldaModel$topics)
  
  topNtermsPerTopic <- top.topic.words(ldaModel$topics, N, by.score=TRUE)
  allTopicModelTerms <- unique(as.vector(topNtermsPerTopic))
  
  DTMpreprocessed <- DTMBIN[, allTopicModelTerms]
  DTMpreprocessedCooc <- t(DTMpreprocessed) %*% DTMpreprocessed
  DTMpreprocessedCooc <- t((DTMpreprocessedCooc + 1) / colSums(DTMpreprocessed))
  DTMpreprocessedCooc <- log(DTMpreprocessedCooc)
  DTMpreprocessedCooc <- as.matrix(DTMpreprocessedCooc)
  
  coherence <- rep(0, K)
  pb <- txtProgressBar(max = K)
  for (topicIdx in 1:K) {
    setTxtProgressBar(pb, topicIdx)
    topWordsOfTopic <- topNtermsPerTopic[, topicIdx]
    coherence[topicIdx] <- 0
    for (m in 2:length(topWordsOfTopic)) {
      for (l in 1:(m-1)) {
        mTerm <- topWordsOfTopic[m]
        lTerm <- topWordsOfTopic[l]
        coherence[topicIdx] <- coherence[topicIdx] + DTMpreprocessedCooc[mTerm, lTerm]
      }
    }
  }
  close(pb)
  return(coherence)
}


ldaformat2Matrix <- function (documents, vocab) {
  require(slam)
  stm <- simple_triplet_matrix(
    i = rep(seq_along(documents), sapply(documents, ncol)), 
    j = as.integer(unlist(lapply(documents,"[", 1, )) + 1L), 
    v = as.integer(unlist(lapply(documents,"[", 2, ))), 
    nrow = length(documents), 
    ncol = length(vocab), 
    dimnames = list(names(documents), vocab))
  dtm <- sparseMatrix(i=stm$i, j=stm$j, x=stm$v, dims=c(stm$nrow, stm$ncol), dimnames = dimnames(stm))
}


peekIntoModelCorpus <- function(sampleCorpusFulltext, topicModel, topicToInvestigate, topicThresholdInDocument = NULL, n = 1) {
  topicProportions <- t(topicModel$document_sums) / colSums(topicModel$document_sums)
  if (!is.null(topicThresholdInDocument)) {
    idx <- which(topicProportions[, topicToInvestigate] > topicThresholdInDocument)
    docSampledIds <- sample(idx, min(n, length(idx)), n)
  } else {
    docSampledIds <- order(topicProportions[, topicToInvestigate], decreasing = T)[1:n]
  }
  sampleTexts <- sampleCorpusFulltext[docSampledIds]
  html_print(HTML(paste0(sampleTexts, collapse = "<hr><br><br/><br/>")))
}


getModelAlignments <- function(modelEvaluationData, topWordsToMatch = 100, similarityThreshold = 0.2, verbose = F) {
  source("util/compare.R")
  numModels <- nrow(modelEvaluationData)
  if (numModels < 2) stop("Nothing to compare, got just than one model!")
  cat(c("Parameters: nWords", topWordsToMatch, "| threshold", similarityThreshold, "\n"))
  pairs <- combn(as.character(modelEvaluationData$modelFileName), 2, simplify = F)
  allReliabilities <- rep(0, length(pairs))
  i <- 0
  for (pair in pairs) {
    i <- i + 1
    cat(c("------", "\n"))
    cat(c(pair[1], "\n"))
    cat(c(pair[2], "\n"))
    tm1 <- get(load(file = pair[1]))
    tm2 <- get(load(file = pair[2]))
    alignment <- alignTopicModels(tm1, tm2, topWordsToMatch, similarityThreshold)
    printAlignedTopics(alignment, verbose = verbose)
    allReliabilities[i] <- alignment$reliability
  }
  return(allReliabilities)
}
