# Topic mapping (Lancichinetti et al. 2015)

source("util/co-occurrence.R")

computeInitialAssignments <- function(corpusData, corpusVocab, K, nCoocs = 1000, clusterMethod = "PAM", cooccurrenceMeasure = "LL") {
  
  # compute term cooccurrences for LDA initialization
  binDTM <- ldaformat2Matrix(corpusData, corpusVocab)
  binDTM[binDTM > 1] <- 1
  topNTerms <- colnames(binDTM)[order(colSums(binDTM), decreasing = T)[1:nCoocs]]
  binDTM <- binDTM[, topNTerms]
  corpusTTM <- t(binDTM) %*% binDTM
  corpusTTM[corpusTTM < 3] <- 0
  
  if (cooccurrenceMeasure == "DICE") {
    corpusTTMsig <- computeSignificanceDice(as.matrix(corpusTTM))
    corpusTTMsig[corpusTTMsig < 0.005] <- 0
  } else if (cooccurrenceMeasure == "LL") {
    corpusTTMsig <- computeSignificanceLL(as.matrix(corpusTTM), nrow(binDTM))
    corpusTTMsig[corpusTTMsig < 6.63] <- 0
  } else if (cooccurrenceMeasure == "POISSON") {
    DTM <- ldaformat2Matrix(corpusData, corpusVocab)
    topNTerms <- colnames(DTM)[order(colSums(DTM), decreasing = T)[1:nCoocs]]
    DTM <- DTM[, topNTerms]
    corpusTTM <- t(DTM) %*% DTM
    diag(corpusTTM) <- colSums(DTM)
    corpusTTM[corpusTTM < 3] <- 0
    corpusTTMsig <- computeSignificancePoisson(as.matrix(corpusTTM), rowSums(DTM), pLevel = 0.05)
    #getTopCoocs(corpusTTMsig)
  } else {
    stop(paste0("cooccurrence significance measure ", cooccurrenceMeasure, " is undefined"))
  }
  # print(getTopCoocs(corpusTTMsig))
  
  clusteredTerms <- clusterTerms(corpusTTMsig, K, method = clusterMethod)
  print(sapply(1:K, FUN=function(x) topNTerms[clusteredTerms == x]))
  
  # set fixed random seed for initialization of neglected terms in coccurrence computation
  set.seed(1000)
  vocabClusterAssignments <- sample.int(K, size = length(corpusVocab), replace=TRUE) - 1
  as.numeric(Sys.time())-> t; set.seed((t - floor(t)) * 1e8 -> seed)
  
  for (k in 1:K) {
    terms <- topNTerms[clusteredTerms == k]
    vocabClusterAssignments[corpusVocab %in% terms] <- k - 1
  }
  initialAssignments <- lapply(corpusData, FUN = function(x) {vocabClusterAssignments[(x[1, ] + 1)]})
  initialAssignments <- lapply(initialAssignments, as.integer) 
  return(initialAssignments)

}


clusterTerms <- function(ttm, K, method = "PAM"){
  print("Clustering terms (may take a while)")
  
  topNTerms <- colnames(ttm)
  
  if (method == "PAM") {
    
    # partitioning around mediods
    require(cluster)
    pamResult <- pam(ttm, K)
    t <- table(pamResult$clustering)
    t[t > 1]
    return(pamResult$clustering)
  
  } else if (method == "KMST") {
    
    # k minimum spanning tree
    require("igraph")
    g <- as.undirected(graph.adjacency(ttm * -1, mode="undirected", weighted=TRUE))
    g_mst <- mst(g, method="prim")
    g_mst <- delete_edges(g_mst, order(E(g_mst)$weight)[1:(K-1)])
    #count_components(g_mst)
    return(components(g_mst)$membership)
  
  } else if (method == "INFOMAP") {
    
    # infomap
    require("igraph")
    g <- as.undirected(graph.adjacency(ttm, mode="undirected", weighted = TRUE))
    g_imap <- infomap.community(g)
    t <- table(g_imap$membership)
    t[t > 1]
    sapply(1:K, FUN=function(x) topNTerms[g_imap$membership == x])
    return(g_imap$membership)
    
  } else {
    stop(paste0("clustering method ", method, " is undefined"))
  }
  
}
