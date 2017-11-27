computeSignificancePoisson = function(ttmatrix, wordsInDocs, pLevel = 0.05) {
  
  # Poisson weight (Lancichinetti et al 2015)
  if (FALSE) {
    a <- 10
    b <- 200
    wordsInDocs <- round(runif(1000, 10, 100))
    Lc <- sum(wordsInDocs) ^ 2
    Ld <- sum(wordsInDocs ^ 2)
    lambda <- (a * b / Lc) * Ld
    plot(dpois(0:9, lambda), type="o")
    abline(h = 0.05)
    z <- qpois(pLevel, lambda, lower.tail = F)
    sum(dpois(z:100000, lambda))
    z
  }
  
  Lc <- sum(wordsInDocs) ^ 2
  Ld <- sum(wordsInDocs ^ 2)
  
  m <- nrow(ttmatrix)
  n <- ncol(ttmatrix)
  
  sigMatrix <- matrix(0,nrow=m,ncol=n)
  rownames(sigMatrix) <- rownames(ttmatrix)
  colnames(sigMatrix) <- colnames(ttmatrix)
  print("Computing poisson weights for co-occurrences")
  pb <- txtProgressBar(max=m)
  for (i_a in 1:m) {
    setTxtProgressBar(pb, value=i_a)
    for (i_b in 1:n) {
      # sentences containing a
      a = ttmatrix[i_a, i_a]
      
      # sentences containing b
      b = ttmatrix[i_b, i_b]
      
      # sentences containing a and b
      k = ttmatrix[i_a, i_b]
      
      if (a == 0 || b == 0 || k == 0 || i_a == i_b) {
        sigMatrix[i_a, i_b] = 0;
      } else {
        lambda <- (a * b / Lc) * Ld
        sigMatrix[i_a, i_b] = k - qpois(pLevel, lambda, lower.tail = F)
      }
    }
  }
  close(pb)
  sigMatrix[sigMatrix < 0] <- 0
  # getTopCoocs(sigMatrix)
  return(sigMatrix)
}


computeSignificanceLL = function(ttmatrix, sizeOfCorpus) {
  
  c <- sizeOfCorpus
  
  m <- nrow(ttmatrix)
  n <- ncol(ttmatrix)
  
  sigMatrix <- matrix(0,nrow=m,ncol=n)
  rownames(sigMatrix) <- rownames(ttmatrix)
  colnames(sigMatrix) <- colnames(ttmatrix)
  print("Computing log-likelihood statistic for co-occurrences")
  pb <- txtProgressBar(max=m)
  for (i_a in seq(1,m)) {
    setTxtProgressBar(pb, value=i_a)
    for (i_b in seq(1,n)) {
      # sentences containing a
      a = ttmatrix[i_a, i_a]
      
      # sentences containing b
      b = ttmatrix[i_b,i_b]
      
      # sentences containing a and b
      k = ttmatrix[i_a,i_b]
      
      # Log Likelihood Significance
      if (a == 0 || b == 0 || k == 0 || i_a == i_b) {
        sigMatrix[i_a,i_b] = 0;
      } else {
        # in case a or b equals k we need to increase k a tiny bit, to not produce NaNs
        if (k==a) a = a + 0.001;
        if (k==b) b = b + 0.001;
        sigMatrix[i_a,i_b] = 2 * ((c * log(c)) - (a * log(a)) - (b * log(b)) + (k * log(k)) 
                                  + (c - a - b + k) * log(c - a - b + k) 
                                  + (a - k) * log(a - k) + (b - k) * log(b - k) 
                                  - (c - a) * log(c - a) - (c - b) * log(c - b))
      }      
    }
  }
  close(pb)
  return(sigMatrix)
}


computeSignificanceDice = function(ttmatrix, listOfTerms = NULL) {
  
  m <- nrow(ttmatrix)
  n <- ncol(ttmatrix)
  
  if (!is.null(listOfTerms)) {
    ttmatrix <- ttmatrix[listOfTerms,]
    m <- length(listOfTerms)
  }
  
  sigMatrix <- matrix(0,nrow=m,ncol=n)
  rownames(sigMatrix) <- rownames(ttmatrix)
  colnames(sigMatrix) <- colnames(ttmatrix)
  print("Computing dice statistic for co-occurrences")
  pb <- txtProgressBar(max=m)
  for (i_a in seq(1,m)) {
    setTxtProgressBar(pb, value=i_a)
    for (i_b in seq(1,n)) {
      # sentences containing a
      a = ttmatrix[i_a, i_a]
      
      # sentences containing b
      b = ttmatrix[i_b,i_b]
      
      # sentences containing a and b
      k = ttmatrix[i_a,i_b]
      
      #print(paste0(c(a,b,k), collapse=" "))
      
      # Dice Koeffizient
      if ((a + b) == 0) {
        sigMatrix[i_a,i_b] = 0;
      } else {
        if (i_a == i_b) {
          sigMatrix[i_a,i_b] = 0;
        } else {
          sigMatrix[i_a,i_b] = (2 * k) / (a + b);
        }
      }
    }
  }
  close(pb)
  return(sigMatrix)
}


getTopCoocs = function(coocSentenceSignificance, decr = TRUE, dtm2 = NULL, N = 300) {
  dictTermsDTM <- colnames(coocSentenceSignificance)
  topCoocsIdx <- arrayInd(order(coocSentenceSignificance, decreasing=decr)[1:(N*2)], dim(coocSentenceSignificance))
  
  if (is.null(dtm2)) {
    topCoocs <- t(apply(topCoocsIdx, 1, FUN=function(x) c(dictTermsDTM[x],coocSentenceSignificance[x[1],x[2]])))
  } else {
    topCoocs <- t(apply(topCoocsIdx, 1, FUN=function(x) c(dictTermsDTM[x],coocSentenceSignificance[x[1],x[2]],dtm2[x[1],x[2]],   dtm2[x[1],x[2]] * 100 / coocSentenceSignificance[x[1],x[2]])))
  }
  
  topCoocs <- topCoocs[seq(1,(N*2),2),]
  return(topCoocs)  
}