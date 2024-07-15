###  A refined working version of Prostrate_Params.R
###  Script for randomly missing on GWAS on prostrate data
###  Just randomly missing 5%

setwd("/Users/maryamonifade/Documents/Prostrate_Testing")

data <- data.frame(Comp_data[,1], Comp_data[,2], Comp_data[, 5])
colnames(data) <- c("ID", "BMI", "Geno")


Em_miss<- function(data=my_data, M=100, niter=50, K=kin2, beta=beta$V1[3], tol = .00001){

  geno_comp<- data[,3]
  Xstar<- as.matrix(cbind(x0=1, x1=geno_comp)) # Complete genotype vector with intercept
  ##Use genotype frequencies of a complete SNP for initial genotype frequencies
  geno_prop<- c(0.3834,0.4705,0.2057)
  ###Initial estimates from GASTON
  theta0<- c()
  theta0[1:ncol(Xstar)]<- c(26.949, beta)
  theta0[3]<- 10.7414  ##variance estimates for the random effect
  theta0[4]<- 7.00284  ## variance estimate for the residual
  theta0[5:7] <- geno_prop
  id <- data[,1]
  u.id <- unique(id)
  n <- length(u.id)
  y<-data[,2]
  beta0<-theta0[1:ncol(Xstar)]
  sigmasqub0<-theta0[3] ##Variance for the random effect
  sigmasque0<-theta0[4] ## Variance estimates of the residual error
  initial.theta<-theta0
  theta.hist<-theta0
  
  for (it in 1:niter) {
    D <- sigmasqub0 * K
    ## Variance of y
    Vary<- D + diag(sigmasque0, nrow=n, ncol=n) 
    ## variance due to the sufficient statistics: V^{t} in our document
    Vy<- solve(solve(D) + diag(rep(sigmasque0, n))) 
    
    #Step 2(i): Randomly delete 5% of the genotypes: this is equal to about 125 observations

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
        #y_imis <- y_mis[id_mis == u.idMis[i]]
        miss_geno_probs[i, j+1]<- dnorm(y_mis[i], (theta0[1] + j*theta0[2]), sqrt(Missing_ind_var[i])) * geno_prop[j+1]
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
    print(length(g_obs))
    gl <-list()
    for (i in 1:nrow(sample_outcomes)){
      gl[[i]] <- c(g_obs, sample_outcomes[i,])
    }
    
    #-----//----------//----------//--------
    ## Now for each lth sampled genotype vector, combine it with the intercept column
    X_lsampled<- list()
    for (i in 1:nrow(sample_outcomes)){
      X_lsampled[[i]]<-   as.matrix(cbind(Xstar[,1], gl[[i]]))
    }
    
    
    #-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||--
    # STEP2e: Obtain b_l and V_l at the t+1th iteration
    # At the t=0 iteration: 
    #-----||-----||-----||-----||-----||-----||-----||-----||----||-----||--
    #-----||-----||-----||-----||-----||-----||-----||-----||----||----||--
  
    b_l1 <- list()
    
    for (l in seq_along(X_lsampled)) {
      X_l_current <- X_lsampled[[l]]
      
      # Compute the term y - g_l * beta_g - X * beta_X
      term <- y -  X_l_current %*% beta0 
      
      # Compute b_l_next
      b_l_next <- (1 / sigmasque0) * Vy %*% term
      
      # Store the result in the list
      b_l1[[l]] <- b_l_next
    }
    
    
    #-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||--
    # STEP2f:Compute betastar at the t+1th iteration.
    # At the t=0 iteration: 
    #-----||-----||-----||-----||-----||-----||-----||-----||----||-----||--
    #-----||-----||-----||-----||-----||-----||-----||-----||----||----||--
    
    # Compute the sum of X_lsampled' X_lsampled for all l
    sum_Xl_Xl <- Reduce("+", lapply(X_lsampled, function(Xl) t(Xl) %*% Xl))
    
    # Compute the quantity
    beta_star_next <- t(solve(sum_Xl_Xl))
   
    ### Maybe a better way of computing the sum in the list: 
    sum_term <- Reduce("+", lapply(1:M, function(l) {
      t(y - b_l1[[l]]) %*% X_lsampled[[l]]
    }))
    
    # Compute the quantity
    result <- t(sum_term)
    
    betastar1<- beta_star_next %*% result
    
    #-----||-----||-----||-----||-----||-----||-----||-----||-----||-----||--
    # STEP2g:Estimating covariance parameters at the t+1th iteration.
    # At the t=0 iteration: 
    # For sigmaque
    #-----||-----||-----||-----||-----||-----||-----||-----||----||-----||--
    #-----||-----||-----||-----||-----||-----||-----||-----||----||----||--
    
    # Compute the trace of V
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
    #-----||-----||-----||-----||-----||-----||-----||-----||----||-----||--
    #-----||-----||-----||-----||-----||-----||-----||-----||----||----||--
    
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
    theta1 <- c(betastar1, sigmasque1, sigmasqub1, geno_prop1)
    theta.hist <- rbind(theta.hist, theta1)
    
  }
  
  #################################################################################
  theta.hist<-as.data.frame(theta.hist)
  #thetaComb<- as.data.frame(thetaComb)
  colnames(theta.hist) <- c("beta0", "beta1", "sigsque", "sigsqub", "gamma0", "gamma1", "gamma2")
  #rownames(theta.hist) <- paste("it", 0:50, sep = "")
  return(list(theta0 = initial.theta, theta.hat = round(theta1,6), theta.hist = theta.hist))
}

###Prepare for GWAS   
    
Comp_data <- read.table("Comp_data1.txt", header = FALSE)
beta<-read.table("selected_betas.txt")
kin2<-read.table("my_prostrate_data.sXX.txt")

###Use array job for this:
nSNPs<- (ncol(Comp_data)-2)
print(nSNPs)
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
task_id <- as.numeric(slurm_arrayid)

total_files=seq(from=1, to=nSNPs, by=50)

starting = total_files[task_id]

#Compute starting index
ending = total_files[task_id +1] - 1

####Loop through each SNP
for (g in 1:nSNPs){
  my_data <- data.frame(Comp_data[,1], Comp_data[,2], Comp_data[, (g+2)])
  #my_data <- data.frame(Comp_data[,1], Comp_data[,2], Comp_data[, 5])
  colnames(my_data) <- c("ID", "BMI", "Geno")
  #ModelEM<- Em_miss(my_data, M=100, niter=50, kin2, beta$V1[1], tol = .00001)
  ModelEM<- Em_miss(my_data, M=100, niter=50, kin2, beta$V1[g], tol = .00001)
  modelEM.out<- ModelEM$thetaComb
  write(modelEM.out, file = paste("EM_theta_rand", g, ".txt", sep=""), sep= " ", append=FALSE, ncolumns=length(modelEM.out))
}

###Analyse some betas obtained from 1000 SNPs in comparison with betas from GEMMA used in also starting the algorithm

Em1<- read.table("EM_theta1.txt") ##select 186 from selected beta
Em2<- read.table("EM_theta3.txt") ## select from 1000-1500 
Em3<- read.table("EM_theta4.txt") ## select from 1500-2000

# Define the ranges
ranges <- list(1:186, 1000:1500, 1500:2000)

# Initialize an empty vector to store the selected entries
selected_entries <- c()

# Loop through each range and append the selected entries to the vector
for (range in ranges) {
  selected_entries <- c(selected_entries, beta$V1[range])
}

EMM1<- Em1[,2]
EMM2<- Em2[,2]
EMM3<- Em3[,2]
Em_total<- c(EMM1,EMM2,EMM3)

# Load necessary library
library(ggplot2)
iterations<- 1:length(Em_total)
Gema_beta<- selected_entries[3:1188]
EM_beta<- Em_total

data <- data.frame(param1, param2)

# Melt the data frame to long format for ggplot2
#library(reshape2)
#data_long <- melt(data, id = "iterations")

# Plot the data using ggplot2
ggplot(data, aes(x = param1, y = param2)) +
  geom_point(size = 3) +
  geom_line() +
  labs(title = "Scatter Plot of Two Regression Parameters",
       x = "Parameter 1",
       y = "Parameter 2") +
  theme_minimal()


