#Authors:
#Steffania Sierra Galvis. s.ssierragalvis@studenti.unipi.it
#Alla Usova. a.usova@studenti.unipi.it

#### Function that generates random datasets ####

input_X_generation <- function(){
  num_rows <- 100
  num_col <- 1
  real_data <- data.frame(
    matrix(rnorm(num_col*num_rows), nrow = num_rows)
  )
  
  positive_real_data <- data.frame(
    matrix(runif(num_col*num_rows, 0, 1000), nrow = num_rows)
  )
  
  distinct_values=sample(5:20, num_col)
  
  generate_interval_data <- function(num_rows, num_distinct) {
    intervals <- lapply(num_distinct, function(distinct) {
      seq(0, 100, length.out = distinct)
    })
    unlist(lapply(intervals, function(interval) {
      sample(interval, num_rows, replace = TRUE)
    }))
  }
  
  interval_data <- data.frame(
    matrix(generate_interval_data(num_rows, distinct_values),
           nrow = num_rows)
  )
  
  full_dataset <- cbind(real_data, positive_real_data, interval_data)
  
  # Add column names
  colnames(full_dataset) <- c(
    paste("Real_", 1:num_col, sep = ""),
    paste("Positive_Real_", 1:num_col, sep = ""),
    paste("Interval_", 1:num_col, sep = "")
  )
  return(full_dataset)
}

#### Simulations ####

#First simulation. Run the algorithm 100 times with different datasets and with parameters: 
#data = random dataset with 100 rows and 3 columns. K = 3, Iterations = 10

#Second simulation. Run the algorithm 100 times with the same datasets than simulation 1, and with K = 3, and Iterations = 30.

simulation_1 <- list()
simulation_2 <- list()
for (i in 1:100){
  X_i <- input_X_generation()
  result_sim_1 <- datatypes(X_i,K=3,iterations = 10)
  result_sim_2 <- datatypes(X_i,K=3,iterations = 30)
  simulation_1 <- c(simulation_1, list(result_sim_1))
  simulation_2 <- c(simulation_2, list(result_sim_2))
}

#### Results of the simulations ####

#First simulation
results <- list()
for (attrib in 1:3){
  results_byattrib <- c()
  for (i in 1:100){
    results_byattrib <- c(results_byattrib, which.max(simulation_1[[i]][attrib,])) 
  }
  results <- c(results, list(results_byattrib))
}


categories <- c('Real-valued', 'Positive-valued', 'Interval value')
counts <- c(length(which(results[[1]] == 1)), length(which(results[[2]] == 2)), length(which(results[[3]] == 3)))
barplot(counts, names.arg = categories, col = "lightblue", main = "1st simulation with nrow = 100, K= 3, iter = 10", xlab = 'Data type',ylab = "Percentage of prediction", ylim = c(0, 100), width = 0.5)

#Second simulation
results_s2 <- list()
for (attrib in 1:3){
  results_byattrib <- c()
  for (i in 1:100){
    results_byattrib <- c(results_byattrib, which.max(simulation_2[[i]][attrib,])) 
  }
  results_s2 <- c(results_s2, list(results_byattrib))
}

counts_s2 <- c(length(which(results_s2[[1]] == 1)), length(which(results_s2[[2]] == 2)), length(which(results_s2[[3]] == 3)))
barplot(counts_s2, names.arg = categories, col = "lightblue", main = "2nd simulation with nrow = 100, K= 3, iter = 30", xlab = 'Data type',ylab = "Percentage of prediction",ylim = c(0, 100))

