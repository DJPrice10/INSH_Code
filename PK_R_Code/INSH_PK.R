# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Script to implement INSH algorithm for the Pharmacokinetic example,       #
# using the utility evaluation from ACE (Overstall & Woods 2016)            #
#                                                                           #
# Author: David J. Price                                                    #
# Date: 31st January, 2017                                                  #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

require(acebayes)
require(truncnorm)
require(foreach)
require(doParallel)


# Function to propose discrete designs according to a truncated MVN distribution centred at a previous design
discrete_trunc_sample <- function(current_d, minmaxT, sd_T, stepsize, m, mindiff){
  new_design <- matrix(-1, nrow = m, ncol=length(current_d))
  for (i in 1:m){
    while( any(diff(new_design[i,])<mindiff) ){
      # Use the top line to propose designs across a grid of spacing "stepsize"
      # new_design[i,] = round(matrix(rtruncnorm(n = 1, a = minmaxT[1], b = minmaxT[2], mean = as.numeric(current_d), sd = sd_T), nrow = 1, ncol=ntimes)/stepsize)*stepsize
      new_design[i,] = matrix(rtruncnorm(n = 1, a = minmaxT[1], b = minmaxT[2], mean = as.numeric(current_d), sd = sd_T), nrow = 1, ncol=ntimes)
    }
  }
  return(new_design)
}


tic <- Sys.time()
set.seed(1)

# Number of simulations/samples used to evaluate utility (with utilcomp15sig)
B <- 5000

# Number of waves to run INSH
W <- 60

# Number of samples to retain in INSH at each iteration
# ntokeep <- c(rep(150,6), rep(75, 6),  rep(50, 6))
nper <- 12
ntokeep <- c(rep(150, nper), rep(75, nper),  rep(50, nper), rep(25,nper), rep(10,nper))

# Number of new designs to sample about each retained design
m <- c(rep(2, nper), rep(4,nper), rep(6,nper), rep(12,nper), rep(30,nper))

if((length(m)+length(ntokeep))/2 != W){
  stop('Check that m and ntokeep are same length as W')
}


# Number of observations, and range for those times
ntimes <- 15
minmaxT <- c(0, 24)

# Standard deviation of perturbation kernel
pert_sd <- 0.2
sample_sd <- rep(pert_sd, times = ntimes)

# Grid over which to search for obs times
gridsize <- 0.05


# Minimum time allowed between subsequent samples
min_difference <- 0.25

# Sample initial designs uniformly/randomly from grid over 0:gridsize:24
n_init <- 1200
current_wave_designs <- matrix(NA, nrow = n_init, ncol = ntimes);

for (i in 1:n_init){
  dd <- rep(1,ntimes)
  # Make sure subsequent times are atleast 15 minutes apart
  while (min(diff(dd))<0.25){
   # dd <- sort(sample(x = seq(0,24, by=gridsize),size = 15, replace = FALSE))
    dd <- sort(runif(ntimes, min = minmaxT[1], max = minmaxT[2]))
   }
  current_wave_designs[i,] <- dd
}


keep_all <- NULL


cl <- makeCluster(4)
registerDoParallel(cl)

for (w in 1:W){
  print(w)
  
  # Evaluate utility for each of the designs in the current wave using utilcomp15sig from acebayes package
  util <- vector("numeric", length = nrow(current_wave_designs))
  mydesigns <- current_wave_designs/12 -1 # Convert to [-1,1] for utility evaluation

  util <- foreach(d=1:nrow(current_wave_designs), .packages = "acebayes", .combine = c) %dopar% {
    mean(utilcomp15sig(mydesigns[d,], B))
  }
  

  # If utility returns Inf for any design, remove that design
  if(any(is.na(util) | is.infinite(util))){
    omit <- which(is.na(util) | is.infinite(util))
    current_wave_designs <- current_wave_designs[-omit,]
    util <- util[-omit]
  }
  
  # Establish cutoff for which designs to retain
  tmp_util <- sort(util, decreasing = TRUE)
  cutoff <- tmp_util[ntokeep[w]]
  
  # Retain designs
  update_designs <- current_wave_designs[util >= cutoff,]
  
  # Store ALL designs considered, the corresponding utility, and the wave
  keep_all <- rbind(keep_all, cbind(current_wave_designs,util, w))
  
  # IF the optimal design occured in a previous iteration, add back into the designs for consideration so new designs are sampled around this region
  if (keep_all[which.max(keep_all[,ntimes+1]),ntimes+2]!=w){
    update_designs <- rbind(update_designs, keep_all[which.max(keep_all[,ntimes+1]), 1:ntimes])
  }
  
  # Sample new designs
  current_wave_designs <- foreach(r=1:nrow(update_designs), .combine = "rbind", .packages = "truncnorm") %dopar% {
    discrete_trunc_sample(update_designs[r,], minmaxT, sample_sd, gridsize, m[w], min_difference)
  }

}

stopCluster(cl)
(toc <- Sys.time() - tic)

# Find optimal
opt_loc <- which.max(keep_all[,ntimes+1])
(Opt_design <- keep_all[opt_loc,1:ntimes])
(Opt_util <- keep_all[opt_loc,ntimes+1])
(Total_N_designs <- nrow(keep_all))


write.csv(x = keep_all, file = "INSH_PK_run.csv")





