### R script for Wald test for the randomly missing datasets
### June 12th. 
### For the cluster

library(dplyr)

Em_wald<- function(data, M=100, K, beta, theta1){
  
  #First obtain genotype frequencies: 
  geno_comp<- data[,3]
  Xstar<- as.matrix(cbind(x0=1, x1=geno_comp)) # Complete genotype vector with intercept
  geno_prop<- c(0.3834,0.4705,0.2057)
  
  ###Initial estimates from GASTON
  theta0<- c()
  theta0[1:ncol(Xstar)]<- c(26.949, beta) 
  theta0[3]<- 10.7414 ###sigmasqub
  theta0[4]<- 7.00284 ##sigmasque
  theta0[5:7] <- geno_prop
  id <- data[,1]
  u.id <- unique(id)
  n <- length(u.id)
  y<- data[,2]
  beta0<- theta0[1:ncol(Xstar)]
  sigmasqub0<- theta0[3] ##Variance for the random effect
  sigmasque0<- theta0[4] ## Variance estimates of the residual error
  initial.theta <- theta0
  theta.hist <- theta0
  
  #for (it in 1:niter) {
  #D <- sigmasqub0 * K
  D <- c(theta1[3]) * K
  #Vary<- D + diag(sigmasque0, nrow=n, ncol=n) 
  Vary<- D + diag(c(theta1[4]), nrow=n, ncol=n) 
  Vy<- solve(solve(D) + diag(rep(theta1[4], n))) 
 
  #-----||-----||-----||-----||-----||-----||-----||-----||-----
  #STEP2: Compute weights for the individuals with missing genotype data at thetahat
  
  ## For randomly missing individuals, we first randomly delete these observations
  ## maybe 5%, 10% or 20 % from the complete dataset.
  ## Then we obtain the IDS of these individuals so that we can insert the weights
  ## we obtain back into the full dataset after computing. 
  
  #Step 2(i): Define the missing and observed data:
  
  # Initialize data frames
  rand_miss <- data.frame()    # Data frame to hold the missing genotype data
  clean_data <- data.frame()
  dat_missing <- data.frame()
  Obs_data1 <- data.frame()
  non_NA <- data.frame()       # Should just be the whole dataset if there's no NA
  
  # Sample function to randomly delete data (for demonstration purposes)
  randomly_delete_data <- function(df) {
    n <- nrow(df)
    delete_count <- ceiling(0.05 * n)
    delete_indices <- sample(1:n, delete_count)
    
    deleted_data <- df[delete_indices, ]
    remaining_data <- df[-delete_indices, ]
    
    return(list(remaining_data = remaining_data, deleted_data = deleted_data))
  }
  
  # Check if there are any missing values in the Geno column
  if (any(is.na(data$Geno))) {
    # Loop through each row
    for (i in 1:nrow(data)) {
      # Check if the value in the Geno column is NA
      if (is.na(data$Geno[i])) {
        # If NA is encountered, copy the entire row to 'rand_miss'
        rand_miss <- rbind(rand_miss, data[i, ])
      } else {
        # Copy the non-NA row to 'non_NA'
        non_NA <- rbind(non_NA, data[i, ])
      }
    }
    
    # Process non_NA data to create clean_data
    clean_data <- non_NA
    
    # Apply the randomly_delete_data function to clean_data
    result1 <- randomly_delete_data(clean_data)
    Obs_data1 <- result1$remaining_data
    miss_data1 <- result1$deleted_data
    
    # Combine rand_miss and miss_data1 to form dat_missing
    dat_missing <- rbind(rand_miss, miss_data1)
  } else {
    # If no NA values are found, process the entire dataset
    result1 <- randomly_delete_data(data)
    Obs_data1 <- result1$remaining_data
    dat_missing <- result1$deleted_data
  }
  
  ##Now we compute weights for individuals with missing genotypes
  ## Using the phenotype information for the individuals with missing genotypes:
  
  missing_geno_weights<- matrix(rep(0, nrow(dat_missing) * 3), ncol = 3)
  #Ind_var<- Vary[row(Vy) == col(Vy)]
  Ind_var<- Vary[row(Vary) == col(Vary)] ##Obtain the diagonal elements of the variance matrix 
  # Subset to obtain those for just the missing individuals
  Missing_ind_var<- Ind_var[as.integer(rownames(dat_missing))] 
  y_mis<- dat_missing$BMI
  
  miss_geno_probs<- matrix(rep(0, nrow(dat_missing) * 3), ncol = 3)
  for (i in seq_along(y_mis)) {
    for (j in 0:2){
      mean_value <- theta1[1] + j * theta1[2]
      mean_value <- as.numeric(mean_value[, 1])
      #y_imis <- y_mis[id_mis == u.idMis[i]]
      # miss_geno_probs[i, j+1]<- dnorm(y_mis[i], (c(theta1[1])) + j*(c(theta1[2])), sqrt(Missing_ind_var[i])) * geno_prop[j+1]
      miss_geno_probs[i, j+1]<- dnorm(y_mis[i], mean_value, sqrt(Missing_ind_var[i])) * geno_prop[j+1]
    }
  }
  
  missing_geno_weights<- matrix(rep(0, nrow(dat_missing) * 3), ncol = 3)
  for (i in seq_along(y_mis)){
    for (j in seq_along(1:3)){
      missing_geno_weights[i,j]<- miss_geno_probs[i,j]/sum(miss_geno_probs[i,])
    }
  }
  
  #-----||-----||-----||-----||-----||-----||-----||-----||-----
  #STEP2b: Prepare the weights for the individuals with complete genotype data 
  # Recall the weight is 1 for the actual genotype of the individual with complete 
  # genotypes and 0 for the other two genotypes. We show this below:
  #-----||-----||-----||-----||-----||-----||-----||-----||----
  #First initialise 
  Weights_complete<- matrix(rep(0, nrow(Obs_data1) * 3), ncol = 3)
  for (i in 1:nrow(Obs_data1)){
    if (Obs_data1$Geno[i] =='0'){
      Weights_complete[i,1]<- 1
      Weights_complete[i,2]<- 0
      Weights_complete[i,3]<- 0
    } else if (Obs_data1$Geno[i] =='1'){
      Weights_complete[i,1]<- 0
      Weights_complete[i,2]<- 1
      Weights_complete[i,3]<- 0
    }
    else{
      Weights_complete[i,1]<- 0
      Weights_complete[i,2]<- 0
      Weights_complete[i,3]<- 1
    }
  }
  
  
  #-----||-----||-----||-----||-----||-----||-----||-----||-----
  #STEP2c: Monte Carlo Sampling:
  # A genotype with level 0,1,2 is sampled for each individual 
  # with missing genotype based vectors from the distribution.
  # missing_geno_weights.
  #-----||-----||-----||-----||-----||-----||-----||-----||----
  #-----||-----||-----||-----||-----||-----||-----||-----||----
  M=100
  sample_outcomes<- matrix(ncol=nrow(missing_geno_weights), nrow = M)
  genof= c(0,1,2)
  for (mis_ind in 1:nrow(missing_geno_weights)){
    for (l in seq(1,M)){
      sample_outcomes[l,mis_ind]<- sample(genof, 1, replace = TRUE, prob = missing_geno_weights[mis_ind,])
    }
  }
  
  #-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||--
  # STEP2d: Obtain the g_l i.e the lth genotype vector over all individuals.
  # This is a combination of gl= g_obs, g_(mis,l) where g_obs is the vector of 
  # complete genotypes and g_(mis_l) is the lth sampled genotype vector for the 
  # missing individuals. I want to extract every line of the matrix 
  # sample_outcomes and combine it with g_obs to form g_l.
  #-----||-----||-----||-----||-----||-----||-----||-----||----||-----||--
  #-----||-----||-----||-----||-----||-----||-----||-----||----||----||--
  
  g_obs<- c(Obs_data1$Geno)
  gl <-list()
  for (i in 1:nrow(sample_outcomes)){
    gl[[i]] <- c(g_obs, sample_outcomes[i,])
  }
  
  #-----//----------//----------//--------
  ## Now for each lth sampled genotype vector, combine it with the intercept column
  X_lsampled<- list()
  for (i in 1:nrow(sample_outcomes)){
    X_lsampled[[i]] <- as.matrix(cbind(Xstar[,1], gl[[i]]))
  }
  
  #-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||--
  # STEP2e: Obtain b_l and V_l at the t+1th iteration
  # At the t=0 iteration: 
  #-----||-----||-----||-----||-----||-----||-----||-----||----||-----||--
  #-----||-----||-----||-----||-----||-----||-----||-----||----||----||--
  # This code below works but the next one is better on the eyes. 
  # beta1<-theta1[1:ncol(Xstar)]
  # gll  <- list()
  # b_l1 <- list()
  # for (l in 1:M){
  #   gll[[l]]<- (y - (X_lsampled[[l]]%*%matrix(beta1)))
  #   b_l1[[l]] <-  (1/theta1[4]) * Vy %*% gll[[l]]
  # }
  # 
  beta1<- as.numeric(theta1[1:ncol(Xstar)])
  b_l1 <- list()
  for (l in seq_along(X_lsampled)) {
    X_l_current <- X_lsampled[[l]]
    
    # Ensure X_l_current is a numeric matrix
    if (!is.matrix(X_l_current)) {
      X_l_current <- as.matrix(X_l_current)
    }
    
    # Compute the term y - g_l * beta_g - X * beta_X
    term <- y -  X_l_current %*% beta1 
    
    # Compute b_l_next
    b_l_next <- (1 / theta1[,4]) * Vy %*% term
    
    # Store the result in the list
    b_l1[[l]] <- b_l_next
  }
  
  
  #-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||--
  # STEP2f:Compute betastar at the t+1th iteration.
  # At the t=0 iteration: 
  #-----||-----||-----||-----||-----||-----||-----||-----||----||-----||--
  #-----||-----||-----||-----||-----||-----||-----||-----||----||----||--
  # 
  D<- list()
  for (l in 1:M){
    D[[l]]<- t(X_lsampled[[l]]) %*% X_lsampled[[l]]
  }
  # sum_Xl<- Reduce('+', D)
  # XX<-solve(Reduce('+', D)) ###Taking the inverse of the sum of lists
  # 
  # F1<- list()
  # F2<- list()
  # for (l in 1:M){
  #   F1[[l]]<- t(y - b_l1[[l]])
  #   F2[[l]] <- F1[[l]] %*% X_lsampled[[l]]
  # }
  # XZ<- t(Reduce('+', F2))
  # betastar1<- ((XX) %*% XZ)
  
  
  # Compute the sum of X_lsampled' X_lsampled for all l
  if (!is.matrix(X_l_current)) {
    X_l_current <- as.matrix(X_l_current)
  }
  
  #X_lsampled <- as.matrix(X_lsampled) 
  # sum_Xl_Xl is what I have previously called D, it gives the same result.
  sum_Xl_Xl <- Reduce("+", lapply(X_lsampled, function(Xl) t(Xl) %*% Xl))
  
  # Compute the quantity
  # Previously XX gives the same result
  beta_star_next <- t(solve(sum_Xl_Xl))
  
  ### Maybe a better way of computing the sum in the list: 
  sum_term <- Reduce("+", lapply(1:M, function(l) {
    t(y - b_l1[[l]]) %*% X_lsampled[[l]]
  }))
  
  # Compute the quantity
  result <- t(sum_term)
  
  betastar1<- beta_star_next %*% result
  
  
  #####
  ####Dimension check##
  # Assuming the following variables are already defined and properly initialized:
  # X_lsampled: List of matrices or data frames
  # y: Numeric vector
  # b_l1: List of numeric vectors from the previous step
  # M: Number of matrices in X_lsampled
  
  #-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||--
  # STEP2g:Estimating covariance parameters at the t+1th iteration.
  # At the t=0 iteration: 
  # For sigmaque
  #-----||-----||-----||-----||-----||-----||-----||-----||----||-----||--
  #-----||-----||-----||-----||-----||-----||-----||-----||----||----||--
  
  ee <- c()
  e<- list()
  #ee<- list()
  for (l in 1:M){
    e[[l]] <- (y - (X_lsampled[[l]] %*% betastar1) - b_l1[[l]])
    #ee[[l]]<- t(e[[l]]) %*% e[[l]]
    ee<- sum(t(e[[l]]) %*% e[[l]])
    
  }
  # E<- sum(ee)
  #E<- Reduce('+', ee)
  
  #sigmasque1<-(1/n) * (sum(diag(Vy))) + ((1/M) * ee)
  #sigmasque1<- (1/n) * ( (sum(diag(Vy))) + ((1/M) * ee))
  
  ######First computation of this parameter (gives 34580).
  ###checked through again: not correct
  # #E<-0
  # ee <- c()
  # e<- list()
  # #ee<- list()
  # for (l in 1:M){
  #   e[[l]] <- (y - (X_lsampled[[l]] %*% betastar1) - b_l1[[l]])
  #   #ee[[l]]<- t(e[[l]]) %*% e[[l]]
  #   ee<- sum(t(e[[l]]) %*% e[[l]]) ##this is not really the sum of the list, only giving the last product
  # }
  # 
  # # E<- sum(ee)
  # #E<- Reduce('+', ee)
  # 
  # sigmasque1<-(sum(diag(Vy))) + ((1/M) * ee)
  
  # Compute the trace of Vy
  trace_V <- sum(diag(Vy))
  
  # Compute the sum term
  sum_term <- 0
  for (l in 1:M) {
    diff <- y - X_lsampled[[l]] %*% betastar1- b_l1[[l]]
    sum_term <- sum_term + t(diff) %*% diff
  }
  
  # Compute sigma_e2
  sigmasque1 <- (1/n)*(trace_V + (1 / M) * sum_term)
  
  #-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||--
  # STEP2f:Estimating covariance parameters at the t+1th iteration.
  # At the t=0 iteration: 
  # For sigmasqub
  ##-----||-----||-----||-----||-----||-----||-----||-----||----||-----||--
  #-----||-----||-----||-----||-----||-----||-----||-----||----||----||--
  #B<-0
  ##Different computation giving different results?? Will follow the second computation. 
  # invK<- solve(K)
  # bk<- list()
  # for (l in seq_along(M)){
  #   bk[[l]]<- t(b_l1[[l]]) %*% invK %*% b_l1[[l]]
  # }
  # B <- Reduce('+', bk)
  # kv<- sum(diag(solve(K) %*% Vy)) ##trace of k inverse and vy
  # sigmasqub1<- (1/n) * ( kv + (1/M)*(B))
  
  # Compute K^-1
  K_inv <- solve(K)
  # Compute the trace of K^-1 V
  trace_term <- sum(diag(K_inv %*% Vy))
  
  # Compute the sum term
  sum_term <- 0
  for (l in 1:M) {
    b_l <- b_l1[[l]]
    sum_term <- sum_term + t(b_l) %*% K_inv %*% b_l
  }
  
  # Compute sigma_b2
  sigmasqub1<- (1 / n) * (trace_term + (1 / M) * sum_term)
  
  #-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||--
  # STEP2g:Estimating genetic covariate parameters at the t+1th iteration.
  # At the t=0 iteration: 
  #-----||-----||-----||-----||-----||-----||-----||-----||----||-----||--
  #-----||-----||-----||-----||-----||-----||-----||-----||----||----||--
  
  All_ind_weights<- as.data.frame(rbind(Weights_complete, missing_geno_weights))
  ## Now to compute the estimates of the genotype parameters:
  Obs_freq1<- length(which(All_ind_weights$V1 =="1"))
  miss_freq1<- All_ind_weights$V1[All_ind_weights$V1 < 1]
  gamma0<-  (Obs_freq1/n) + (sum(miss_freq1)/n)
  gamma1<- (length(which(All_ind_weights$V2 =="1"))/n) + sum(All_ind_weights$V2[All_ind_weights$V2 <1])/n
  gamma2<- 1-gamma0-gamma1
  geno_prop1 <- c(gamma0,gamma1,gamma2)
  
  ######## Collate all estimates #####
  thetahat <- c(betastar1,sigmasqub1,sigmasque1, geno_prop1)
  #theta.hist <- rbind(theta.hist, theta1)
  
  ### Monte Carlo Approach to the Wald test.
  ### Second derivative for betag
  ## Gives a two by two matrix: To resuse this in the information matrix,
  ## we pick the diagonal elements which correspond to the second derivatives of 
  ## the intercept and betag
  t1<- c(theta1[,1], theta1[,2], theta1[,3], theta1[,4])
  t2 <- t1[1:2]
  
  ## Same result of computing second derivative. Just different method
  #Xll<- Reduce('+', D) ##
  
  ##An alternative to computing Xll:
  ## same results??
  # Assuming the following variables are already defined and properly initialized:
  # X_l: List of matrices, each of dimension 2539 x 2
  # M: Number of matrices in X_l (length of the list X_l)
  # sigma_e2: The value of sigma_e squared (a numeric scalar)
  
  # Initialize the sum of products matrix
  sum_Xl_Xl <- matrix(0, nrow = ncol(X_lsampled[[1]]), ncol = ncol(X_lsampled[[1]]))
  
  # Iterate over each matrix in the list X_l
  for (l in 1:M) {
    X_l_current <- X_lsampled[[l]]
    
    # Compute the product X_l' * X_l
    Xl_prime_Xl <- t(X_l_current) %*% X_l_current
    
    # Sum the products
    sum_Xl_Xl <- sum_Xl_Xl + Xl_prime_Xl
  }
  
  # Compute the final result
  betagsqu <- (1 / (M * t1[4])) * sum_Xl_Xl
  
  ### First derivative computation:
  ##Initialize list to store components of the bracket
  # FF1<- list()
  # for (l in 1:M){
  #   FF1[[l]]<- (t(y - b_l1[[l]]) %*% X_lsampled[[l]] - (t(t2) %*% t(X_lsampled[[l]]) %*% X_lsampled[[l]]))^2
  # }
  # FF2<- (1/(M * t1[4])) * (Reduce('+', FF1))
  
  # More detailed calculation but same results: 
  # Initialize the sum
  sum_term <- 0
  
  # Iterate over each matrix/vector in the lists
  for (l in 1:M) {
    X_l_current <- X_lsampled[[l]]
    b_l_current <- b_l1[[l]]
    
    # Compute the term (y - b_l) %*% X_l
    term1 <- t(y - b_l_current) %*% X_l_current
    
    # Compute the term beta_star_prime %*% t(X_l) %*% X_l
    term2 <- t(t2) %*% t(X_l_current) %*% X_l_current
    
    # Compute the difference and square it
    diff <- (term1 - term2)
    diff_squared <- diff %*% t(diff)
    
    # Sum the terms
    sum_term <- sum_term + diff_squared
  }
  
  
  #  Compute the final result
  #  Same as FF2
  result <- (1 / (M * t1[4])) * sum_term
  
  ## Monte Carlo Approximation for the standard error of betag
  ## Since the betagsqu is block diagonal 
  varBeta<- betagsqu[2,2] + result
  varBeta<- solve(infor_mat) #variance of betahat is the inverse of the information matrix
  standard_error <- sqrt(varBeta)
  estimated_parameter <- t1[2]
  
  ### Wald test statistic
  compute_wald_statistic <- function(est_parameter, varBeta) {
    wald_statistic <- (est_parameter)^2 / varBeta
    return(wald_statistic)
  }
  wald_value <- compute_wald_statistic(estimated_parameter, varBeta)
  p_value1 <- pchisq(wald_value, df =1, lower.tail = FALSE)
  thetap <- c(t1, standard_error, p_value1)
  
  # ### Wald test statistic
  # compute_wald_statistic <- function(est_parameter, standard_error) {
  #   wald_statistic <- est_parameter / standard_error
  #   return(wald_statistic)
  # }
  # wald_value <- compute_wald_statistic(estimated_parameter, standard_error)
  # degrees_of_freedom <- 1
  # p_value1 <- pchisq(wald_value, df = degrees_of_freedom, lower.tail = FALSE)
  # thetap <- c(t1, standard_error, p_value1)
  
  return(list(theta.hat=thetahat, theta.hist = thetap))
  
}


#---#---#---# Call the GWAS_Wald function #---#---#---#---
###Use array job for this:
###Prepare for GWAS
Comp_data<-read.table("/home/monif064/projects/def-kburkett/monif064/projects/Maryam/Prostrate_randomMissing/Comp_data1.txt", header=FALSE)
## Read in SNP effect sizes
beta<-read.table("/home/monif064/projects/def-kburkett/monif064/projects/Maryam/Prostrate_randomMissing/selected_betas.txt")
##Read in GRM
kin2<-read.table("/home/monif064/scratch/Prostrate_QC_files/New_Prostrate_QC/Files_to_graham/my_prostrate_data.sXX.txt")
##Concatenated version of the EM estimates of the SNPs in Comp_data
Em_estimate<- read.table("EM_theta1.txt")

##Set up array job
nSNPs<- (ncol(Comp_data)-2)
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
task_id <- as.numeric(slurm_arrayid)

total_files=seq(from=1, to=nSNPs, by=50)

starting = total_files[task_id]

#Compute starting index

ending = total_files[task_id +1] - 1

####Loop through each SNP
for (g in 1:nSNPs){
  my_data <- data.frame(Comp_data[,1], Comp_data[,2], Comp_data[, (g+2)])
  colnames(my_data) <- c("ID", "BMI", "Geno")
  theta1<- Em_estimate[g,]
  ModelEM<- Em_wald(my_data, M=100, niter=50, kin2, beta$V1[g], theta1, tol = .00001)
  modelEM.out<- ModelEM$theta.hist
  write(modelEM.out, file = paste("EM_theta_rand", g, ".txt", sep=""), sep= " ", append=FALSE, ncolumns=length(modelEM.out))
}


