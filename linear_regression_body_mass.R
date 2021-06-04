# *****************************************************************
# Linear Regression to find min population size for any body mass *
# *****************************************************************
population_size <- function(bm, rm){
  
  #Convert the given body mass in log, into a dataframe
  bm_df <- data.frame(log(bm))
  names(bm_df) <- c('Body.mass')
  
  #Select the required data according to the rm value
  if (rm == 0.6) {
    data <- as.data.frame(cbind(log(mean_min_pop$Population[1:31]), log(mean_min_pop$Body.mass[1:31])))
    names(data) <- c('Population', 'Body.mass')}
  
  if (rm == 0.8) {
    data <- as.data.frame(cbind(log(mean_min_pop$Population[32:62]), log(mean_min_pop$Body.mass[32:62])))
    names(data) <- c('Population', 'Body.mass')}
  
  if (rm == 1) {
    data <- as.data.frame(cbind(log(mean_min_pop$Population[63:93]), log(mean_min_pop$Body.mass[63:93])))
    names(data) <- c('Population', 'Body.mass')}
  
  #Build Linear Model
  lmMod <- lm(Population ~ Body.mass, data = data)
  
  #Check the quality of the linear regression
  r2 <- summary(lmMod)$r.squared
  pvalue <- summary(lmMod)$coefficients[2,4]
  
  #Continue only if the model correctly represents the data
  if (r2 < 0.6 & pvalue > 0.05) stop("The linear regression model does not represent correctly the data")

  #Predict the Minimum Population Size for the given body mass from the Linear Model
  return(as.integer(exp(predict(lmMod, bm_df))))
}

population_size(bm, rm)