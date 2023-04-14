#---------------------------------------
#----- Define R functions --------------
#---------------------------------------


## Generate binary survey data
generateSurveyResponse <- function (n1, n2, p1, p2, summary=FALSE){
  # n1: sample size in the control group
  # n2: sample size in the treatment group
  # p1: positive response rate in the control group
  # p2: positive response rate in the treatment group
  # summary: if True return 2x2 contingency table; if False return raw data
  # return: simulated survey response data
  
  response1 <- rbinom(n1, 1, p1)
  response2 <- rbinom(n2, 1, p2)
  response <- list(response1=response1, response2=response2)
  
  if(summary){
    npr1 <- sum(response1)
    npr2 <- sum(response2)
    response <- matrix(c(n1, n2, npr1, npr2), 2, 2)
  }
  
  return(response)
}

## Get the normalized Z-statistic
normalizedz <- function(data){
  # data:2x2 contingency table
  # |-----------------------|-----------------------------------|
  # | control sample size   | control positive response count   |
  # |-----------------------|-----------------------------------|
  # | treatment sample size | treatment positive response count |
  # |-----------------------|-----------------------------------|
  # return: Z-statistic = sample mean difference / its standard deviation
  
  n1 <- data[1,1]
  n2 <- data[2,1]
  m1 <- data[1,2]
  m2 <- data[2,2]
  
  p1 <- m1/n1
  p2 <- m2/n2
  p <- (m1+m2)/(n1+n2)
  
  if(p==0 | p==1){
    z_stat <- 0
  } else {
    z_stat <- (p2-p1)/sqrt(p*(1-p)*(1/n1 + 1/n2))
  }
  
  return(z_stat)
}


## Calculate the absolute lift and relative lift
liftCalculator <- function(data){
  # data:2x2 contingency table
  # |-----------------------|-----------------------------------|
  # | control sample size   | control positive response count   |
  # |-----------------------|-----------------------------------|
  # | treatment sample size | treatment positive response count |
  # |-----------------------|-----------------------------------|
  # return: absolute lift and relative lift
  
  n1 <- data[1,1]
  n2 <- data[2,1]
  m1 <- data[1,2]
  m2 <- data[2,2]
  prr1 <- m1/n1
  prr2 <- m2/n2
  absolute_lift <- prr2 - prr1
  
  if(m1==0 & m2==0){
    relative_lift <- 0
  } else if (m1==0) {
    relative_lift <- 100000
  } else {
    relative_lift <- (prr2 - prr1)/prr1
  }
  
  return(list(absolute_lift=absolute_lift, relative_lift=relative_lift))
  
}


liftCalculator2 <- function(x1, x2){
  
  n1 <- length(x1)
  n2 <- length(x2)
  p1b <- mean(x1)
  p2b <- mean(x2)
  v1b <- var(x1)*(n1-1)/n1
  v2b <- var(x2)*(n2-1)/n2
  
  absolute_lift <- p2b - p1b
  
  if (p1b==0 & p2b==0){
    relative_lift <- 0
  } else if (p1b==0){
    relative_lift <- 10000
  } else {
    relative_lift <- (p2b - p1b)/p1b
  }
  
  if(p1b==p2b & v1b==0 & v2b==0){
    std_absolute_lift <- 0
  } else if (v1b==0 & v2b==0){
    std_absolute_lift <- 10000
  } else {
    std_absolute_lift <- (p2b-p1b)/sqrt(v1b/n1+v2b/n2)
  }
  
  return(list(absolute_lift=absolute_lift, 
              relative_lift=relative_lift, 
              std_absolute_lift=std_absolute_lift))
}


matrix2long <- function(data, rowgrid, colgrid, threshold){
  
  nrow_mat <- length(rowgrid)
  ncol_mat <- length(colgrid)
  nrow_long <- nrow_mat * ncol_mat
  sample_size <- rep(colgrid, nrow_mat)
  prr <- rep(rowgrid, each=ncol_mat)
  pvalue <- as.vector(data)
  label <- ifelse(pvalue > threshold, 'invalid', 'valid')
  pval <- data.frame(sample_size, prr, pvalue, label)
  return(pval)
  
}


## Perform two sample proportion test with different methods
proportionTest <- function(data, method, B){
  # data:2x2 contingency table
  # |-----------------------|-----------------------------------|
  # | control sample size   | control positive response count   |
  # |-----------------------|-----------------------------------|
  # | treatment sample size | treatment positive response count |
  # |-----------------------|-----------------------------------|
  # return: absolute lift and relative lift
  # method = ("clt", "bootstrap", "bootstrapmean", "permutation")
  #   -clt: standard proportion test with normal approximation
  #   -bootstrap: bootstrap test to test whether two distributions are the same
  #   -bootstrapmean: bootstrap test to test whether the means of two populations are the same
  #   -permutation: permutation test to test whether two distributions are the same
  # B: number of replications for bootstrap test and permutation test
  # return: absolute lift and relative lift together with p-values
  
  lift <- liftCalculator(data)
  absolute_lift <- lift$absolute_lift
  relative_lift <- lift$relative_lift
  normalized_zstat <- normalizedz(data)
  
  n1 <- data[1,1]
  n2 <- data[2,1]
  m1 <- data[1,2]
  m2 <- data[2,2]
  
  if(method == "clt"){
    
    p1 <- m1/n1
    p2 <- m2/n2
    v1 <- p1*(1-p1)
    v2 <- p2*(1-p2)
    
    if(m1==0 & m2==0){
      t_absolute <- 0
      t_relative <- 0
    } else if (m1==0){
      t_absolute <- (p2-p1)/sqrt(v1/n1 + v2/n2)
      t_relative <- 0
    } else {
      t_absolute <- (p2-p1)/sqrt(v1/n1 + v2/n2)
      vr <- v2/(n2*p1^2) + p2^2*v1/(n1*p1^4)
      t_relative <- relative_lift/sqrt(vr)
      
    }
    pvalue_absolute <- 2*(1-pnorm(abs(t_absolute)))
    pvalue_relative <- 2*(1-pnorm(abs(t_relative)))
    lift <- c(absolute_lift, relative_lift)
    pval <- c(pvalue_absolute, pvalue_relative)
    
  }
  
  if(method == 'bootstrap'){
    
    absolute_lift_B <- rep(NA, B)
    relative_lift_B <- rep(NA, B)
    normalizedz_B <- rep(NA, B)
    
    for (b in 1:B){
      n1b <- n1
      n2b <- n2
      m1b <- rbinom(1, n1b, (m1+m2)/(n1+n2))
      m2b <- rbinom(1, n2b, (m1+m2)/(n1+n2))
      datab <- matrix(c(n1b, n2b, m1b, m2b), 2, 2)
      liftb <- liftCalculator(datab)
      absolute_lift_B[b] <- liftb$absolute_lift
      relative_lift_B[b] <- liftb$relative_lift
      normalizedz_B[b] <- normalizedz(datab)
    }
    pvalue_absolute <- sum(abs(absolute_lift_B) > abs(absolute_lift))/B
    pvalue_relative <- sum(abs(relative_lift_B) > abs(relative_lift))/B  
    pvalue_normalizedz <- sum(abs(normalizedz_B) > abs(normalized_zstat))/B  
    
    lift <- c(absolute_lift, relative_lift)
    pval <- c(pvalue_absolute, pvalue_relative, pvalue_normalizedz)
  }
  
  
  if(method == 'bootstrapmean'){
    
    absolute_lift_B <- rep(NA, B)
    relative_lift_B <- rep(NA, B)
    normalizedz_B <- rep(NA, B)
    
    x1 <- c(rep(1, m1), rep(0, n1-m1))
    x2 <- c(rep(1, m2), rep(0, n2-m2))
    p1 <- mean(x1)
    p2 <- mean(x2)
    p12 <- mean(c(x1, x2))
    x1n <- x1 - p1 + p12
    x2n <- x2 - p2 + p12
    
    for (b in 1:B){
      
      x1b <- sample(x1n, n1, replace=T)
      x2b <- sample(x2n, n2, replace=T)
      liftb <- liftCalculator2(x1b, x2b)
      absolute_lift_B[b] <- liftb$absolute_lift
      relative_lift_B[b] <- liftb$relative_lift
      normalizedz_B[b] <- liftb$std_absolute_lift
    }
    pvalue_absolute <- sum(abs(absolute_lift_B) > abs(absolute_lift))/B
    pvalue_relative <- sum(abs(relative_lift_B) > abs(relative_lift))/B  
    pvalue_normalizedz <- sum(abs(normalizedz_B) > abs(normalized_zstat))/B  
    
    lift <- c(absolute_lift, relative_lift)
    pval <- c(pvalue_absolute, pvalue_relative, pvalue_normalizedz)
  }
  
  
  if(method == 'permutation'){
    
    absolute_lift_B <- rep(NA, B)
    relative_lift_B <- rep(NA, B)
    normalizedz_B <- rep(NA, B)
    x_pooled <- c(rep(1, m1+m2), rep(0, n1+n2-m1-m2))
    
    for (b in 1:B){
      xb <- sample(x_pooled)
      n1b <- n1
      n2b <- n2
      m1b <- sum(xb[1:n1])
      m2b <- sum(xb[(n1+1):(n1+n2)])
      datab <- matrix(c(n1b, n2b, m1b, m2b), 2, 2)
      liftb <- liftCalculator(datab)
      absolute_lift_B[b] <- liftb$absolute_lift
      relative_lift_B[b] <- liftb$relative_lift
      normalizedz_B[b] <- normalizedz(datab)
    }
    pvalue_absolute <- sum(abs(absolute_lift_B) > abs(absolute_lift))/B
    pvalue_relative <- sum(abs(relative_lift_B) > abs(relative_lift))/B  
    pvalue_normalizedz <- sum(abs(normalizedz_B) > abs(normalized_zstat))/B
    
    lift <- c(absolute_lift, relative_lift)
    pval <- c(pvalue_absolute, pvalue_relative, pvalue_normalizedz)
  }
  
  return(list(lift=lift, pval=pval))
  
}