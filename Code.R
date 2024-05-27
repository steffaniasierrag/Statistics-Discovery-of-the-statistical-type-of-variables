#Authors:
#Steffania Sierra Galvis. s.ssierragalvis@studenti.unipi.it

#install.packages('MCMCpack')

#### Libraries ####
library(MCMCpack) # for distributions
library(MASS)  #for multivariate distributions
library(numDeriv) #to compute derivatives

#### Variables: ####
# Z: low-rank representation
# {b_ld} vectors of weight of contribution of the latent variables to the original variables
# {y_nld} pseudo-observations associated to x_nd
# w_d vectors of weights of each attribute of the original data
# S matrix of likelihood assignments
# U matrix of noise variables.

#### Algorithm ####
datatypes <- function(input_X, K, iterations){
  
  #### Constants ####
  N = length(input_X[,1]) #number of rows of X  
  D = length(input_X[1,]) #number of columns of X
  s_Ld = 3 #Number of data types. In this case, only continuous data types are considered
  # l = 1,2,3 correspond to the datatypes RealValuedData, PositiveRealValuedData, IntervalData respectively.
  
  #### Standard deviation parameters ####
  sd_z = 1
  sd_b = 1
  sd_y = 1
  sd_u = 0.001
  
  #### Parameters for rescaling #### 
  mu_param = 0
  epsilon_param = (max(input_X) - min(input_X))/10000
  thetaL = min(input_X)-epsilon_param
  thetaH = max(input_X) - epsilon_param
#  w_param = 2/ max(input_X) #for the continuous datatypes
  w_param <- c()  #parameter used in the transformation functions
  for (d in 1:D){
    w_param <- c(w_param, 2/max(input_X[,d]))
  }
  
  #### Continuous variables transformations ####
  RealValuedData <- function(x,d){return(w_param[d]*x+mu_param)}
  PositiveRealValuedData <- function(x,d){return(log(1+exp(w_param[d]*x)))}
  IntervalData <- function(x,d){return((thetaH - thetaL)/(1+exp(-w_param[d]*x))+thetaL)}
  
  #### inverse function of the transformations
  inverse_RVD <- function(y,d) {return((y-mu_param)/w_param[d])}
  inverse_PRVD <- function(y,d) {return(log(exp(1/w_param[d]*y)-1))}
  #inverse_PRVD <- function(y,d) {return(log((1/w_param[d])*exp(y)-1))}
  inverse_IVD <- function(y,d) {return(-1/w_param[d] * log((thetaH-y)/(y-thetaL)))}
  
  #### derivatives of the inverse of the transformations
  D_I_RVD <- function(x,d){return(1/w_param[d])}
  D_I_PRVD <- function(x,d){return(exp(x/w_param[d])/(w_param[d]*(exp(x/w_param[d])-1)))} 
  #D_I_PRVD <- function(x,d){return(exp(x)/w_param[d]*(exp(x)-1))} 
  D_I_IVD <- function(x,d){return((-thetaH+thetaL)/(w_param[d]*(x-thetaL)*(thetaH-x)))}  
  
  #### Initialization of all the variables ####
  
  ### Vector of weights W
  alpha_dir <- rep(1, times = s_Ld) #uniform dirichtlet distribution. This vector will change every time S is sampled. 
  W_matrix <-  matrix(rdirichlet(D, alpha_dir), ncol = s_Ld, byrow = FALSE) #Matrix of size D x |L^d|. Every wd is a row of W. Result: D samples (vectors) of length the same length than the vector of parameters alpha.   
  
  ### Matrix of likelihood assignments S ###
  
  # The entry snd is a variable with multinomial distribution with vector of probabilities wd.
  # rmultinom(n = number of samples, size = sum of the entries of every sample, vector of probabilities p). The result will be n vectors of size len(p) with sum(entries of each vector) = sample.
  # In our case, every snd represents a likelihood assignment corresponding to the observation xnd, which can correspond only to one data type. 
  #Therefore, we want that only one of the entries of the output vector(which has length the number of data types) be non-zero.
  # Since snd depends only of d for every n, we will generate for every d, N sample vectors, and to extract the index of the non-zero component.
  
  S <- list()
  for (d in 1:D) {
    sd <- c()
    multnom_output_d <- matrix(rmultinom(n = N, size = 1, W_matrix[d,]), nrow = length(W_matrix[d,]), byrow = FALSE)
    for (i in 1:N){
      sd <- c(sd, which(multnom_output_d[,i] == 1))
    }
    S <- c(S, list(sd))
  }
  S <- do.call(cbind, S) #Matrix of size N x D.
  
  ### Vectors b_ld
  #Generates s_Ld vectors b_l^d for every d. Every b_ld vector has K-components.
  B <- mvrnorm(n= s_Ld*D, mu = rep(0, length = K), Sigma = sd_b * diag(K))
  B <- t(B) #Vectors b_ld are column vectors. B is a matrix of size K x (s_Ld * D)
  
  ### Matrix Z
  # Z represents the low-rank representation of our original data set X. |Z| = N x K 
  Z <- mvrnorm(n= N, mu = rep(0, length = K), Sigma = sd_z * diag(K))
  
  ###Pseudo-observations y_nld
  Y <- list()  #list of vectors y_id, ordered first by d and n.
  for (d in 1:D){
    for (n in 1:N){
      y_nd <- c() 
      for (l in 1:s_Ld){
        mean_v <- Z[n,,drop = FALSE]%*%B[,s_Ld*(d-1)+l, drop = FALSE]
        y_nld <- rnorm(n = 1, mean = mean_v, sd = sd_y) # component l of the vector y_nd. 
        y_nd <- c(y_nd, y_nld) #add each component l to the vector y_nd
      }
      Y <- c(Y, list(y_nd))
    }
  }
  Y <- do.call(rbind,Y) #|Y| = (N*D) x s_Ld. Each row is a vector y_nd associated to the observation x_nd. 
  
  ### Noise variable U
  U <- matrix(rnorm(n = N*D, mean = 0, sd = sd_u), nrow = N, byrow = TRUE)
  
  ### Auxiliary functions for sampling
  # function used to update the vectors of probabilities used when sampling s_nd as multinomial variable
  prob_xnd_zn_bld_l <- function(d,n,l){
    c1 <- 1/(sqrt(2*pi*(sd_y+sd_u)))
    c2 <- -1 / 2*(sd_y+sd_u)
    if (l == 1){
      der <- abs(D_I_RVD(input_X[n,d],d))
      exp_term <- (inverse_RVD(input_X[n,d],d)-Z[n,, drop = FALSE]%*%B[,(d-1)*s_Ld + l, drop = FALSE])^2
      return(c1*der*exp(c2*exp_term))}
    if (l == 2){
      der <- abs(D_I_PRVD(input_X[n,d],d))
      exp_term <- (inverse_PRVD(input_X[n,d],d)-Z[n,, drop = FALSE]%*%B[,(d-1)*s_Ld + l, drop = FALSE])^2
      return(c1*der*exp(c2*exp_term))}
    if (l == 3){
      der <- abs(D_I_IVD(input_X[n,d],d))
      exp_term <- (inverse_IVD(input_X[n,d],d)-Z[n,, drop = FALSE]%*%B[,(d-1)*s_Ld + l, drop = FALSE])^2
      return(c1*der*exp(c2*exp_term))}
  }
  
  #### Body of the code ####
  iter <- 0
  while (iter <= iterations){
    ### Sampling Z given B and Y ###
    # 1st compute the covariance matrix sigma_z.
    prod_BBt <- list()   #List of matrices b_ld*t(b_ld)
    for (d in 1:D){
      for (l in 1:s_Ld){
        pr_ld <- B[,(d-1)*s_Ld+l,drop = FALSE] %*% t(B[,(d-1)*s_Ld+l,drop = FALSE])   #matrix of size K*K
        prod_BBt <- c(prod_BBt, list(pr_ld))}}
    sigma_z <- Reduce('+', prod_BBt) + (1/sd_z)*diag(K) #it doesn't depends on n.
    sigma_z <- solve(sigma_z)
    
    # Update Z[n,]
    for (n in 1:N){
      mu_z_n <- list()  # vector with mean values
      for (d in 1:D){
        for (i in n:N){
          for (l in 1:s_Ld){
            B_y <- B[,(d-1)*s_Ld + l,drop = FALSE]*Y[(d-1)*N + i, l]
            mu_z_n <- c(mu_z_n, list(B_y))}
        }
      }
      mu_z_n <- sigma_z %*% Reduce('+', mu_z_n) #it depends on n. Multiplication of a matrix KXK by a vector K x 1
      Z[n,] <- mvrnorm(n = 1, mu = mu_z_n, Sigma = sigma_z)
    }
    
    ### Sample y's given x,z,b,s (sampling pseudo-observations) ### 
    ## To update y_ndl there are two cases. If s_nd != l, use the its prior distribution. 
    #Otherwise, the posterior distribution is defined as normal distribution, with parameters mu_hat_y(d,l,n) and sd2y_hat
    
    sd2y_hat <- (1/sd_y + 1/sd_u)^(-1) #it doesn't depend of N, D or s_Ld
    
    # Function that computes the mean parameter of the normal distribution of the y's.
    mu_hat_y <- function(d,l,n){
      term1 <- (Z[n,, drop = FALSE]%*% B[, (d-1)*s_Ld + l, drop = FALSE])/sd_y  
      if (l == 1){term2 <- inverse_RVD(input_X[n,d], d)/sd_u}
      else if (l == 2){term2 <- inverse_PRVD(input_X[n,d], d)/sd_u}
      else if (l == 3){term2 <- inverse_IVD(input_X[n,d], d)/sd_u}
      return((term1+term2)*sd2y_hat)}
    
    sigma_b <- solve( (1/sd_y) * t(Z)%*%Z + (1/sd_b)*diag(K)) #Covariance matrix used to update b_ld. It depends only of Z.
    
    for (d in 1:D){
      for (l in 1:s_Ld){
        #Sample Y_ndl
        for (n in 1:N){ if (S[n,d] == l){
          mu_hat_value <- mu_hat_y(d,l,n)
          if (is.nan(mu_hat_value) | mu_hat_value == 0){next}
          else {Y[(d-1)*N + n, l] <- rnorm(n = 1, mean = Y[(d-1)*N + n, l]/ mu_hat_value, sd = sd2y_hat)}
          }
        }
        
        # Sample b's given Z and y's
        mu_ld_1 <- matrix(numeric(K), ncol = 1)  #First term of mu_ld (mean vector of the normal distr.)
        for (n in 1:N){
          Zn_Y_nld <- t(Z[n,, drop = FALSE])*Y[(d-1)*N + n,l]
          mu_ld_1 <- mu_ld_1 + Zn_Y_nld}
        mu_ld <- sigma_b %*% mu_ld_1
        B[,(d-1)*s_Ld + l] <- mvrnorm(n = 1, mu = mu_ld, Sigma = sigma_b)
      }
      
      # Sample S given X, Z, and B
      for (n in 1:N){
        p_snd <- c() #vector of probabilities to use when sampling snd as multinomial variable
        for (l in 1:s_Ld){
          new_prob_vector <- prob_xnd_zn_bld_l(d,n,l)
          p_snd <- c(p_snd, W_matrix[d,l]*new_prob_vector)
        }
        if (any(is.nan(p_snd)) | all(p_snd==0)){p_snd <- W_matrix[d,]} #the likelihood functions can generate NaN values because of the log.
        else{p_snd <- p_snd/sum(p_snd)}
        S[n,d] <- which(rmultinom(n = 1, size = 1, prob = p_snd)==1)
      }
      
      #sample w's given S
      for (l in 1:s_Ld){alpha_dir[l] <- alpha_dir[l] + length(which(S[,d] == l))}
      W_matrix[d,] <-  rdirichlet(1, alpha_dir)
    }
    iter <- iter + 1
  }
  return(W_matrix)
  }
