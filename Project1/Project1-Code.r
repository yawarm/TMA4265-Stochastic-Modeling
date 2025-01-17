# Amandus Omholt Nygaard, Yawar Mahmood

# Problem 1

# Task 1 c)
set.seed(1)
alfa <- 0.005
beta <- 0.01
gamma <- 0.1
days <- 7300
X_initial <- 0

P <- matrix(c(1-beta, beta, 0, 0, 1-gamma, gamma, alfa, 0, 1-alfa), nrow=3, byrow=TRUE)

MC_simulate <- function(X0, P, steps){
  U <- 1:dim(P)[1]
  X <- rep(NA, steps+1)
  X[1] <- X0 + 1
  for(i in 1:steps){
    X[i+1] <- sample(U, 1, prob = P[X[i], ])
  }
  X <- X - 1
  
  return(X)
}

X = MC_simulate(X_initial, P, days)

# plot with options
plot(X,
     main = 'Simulation of Measles over 20 years (7300 time steps)',
     xlab = 'Time step #',
     ylab = 'State #',)


# set n to the number of times we want to run the realisation
set.seed(1)
n = 30
pi0 <- rep(NA, n)
pi1 <- rep(NA, n)
pi2 <- rep(NA, n)
for(i in 1:n){
  res <- MC_simulate(X_initial, P, days)
  pi <- table(tail(res, days/2))
  pi0[i] <- pi[1]/(days/2)
  pi1[i] <- pi[2]/(days/2)
  pi2[i] <- pi[3]/(days/2)
}

# upper bound
upper_pi0 <- mean(pi0) + qt(0.975, length(pi0)-1) * sd(pi0) / sqrt(length(pi0))
# lower bound
lower_pi0 <- mean(pi0) - qt(0.975, length(pi0)-1) * sd(pi0) / sqrt(length(pi0))
# upper bound
upper_pi1 <- mean(pi1) + qt(0.975, length(pi1)-1) * sd(pi1) / sqrt(length(pi1))
# lower bound
lower_pi1 <- mean(pi1) - qt(0.975, length(pi1)-1) * sd(pi1) / sqrt(length(pi1))
# upper bound
upper_pi2 <- mean(pi2) + qt(0.975, length(pi2)-1) * sd(pi2) / sqrt(length(pi2))
# lower bound
lower_pi2 <- mean(pi2) - qt(0.975, length(pi2)-1) * sd(pi2) / sqrt(length(pi2))

plot(pi0,
     main = 'limiting distribution state 0',
     xlab = 'Realisation #',
     ylab = 'limiting distribution',)
abline(h=upper_pi0, col="blue")
abline(h=lower_pi0, col="blue")

plot(pi1,
     main = 'limiting distribution state 1',
     xlab = 'Realisation #',
     ylab = 'limiting distribution',)
abline(h=upper_pi1, col="blue")
abline(h=lower_pi1, col="blue")

plot(pi2,
     main = 'limiting distribution state 2',
     xlab = 'Realisation #',
     ylab = 'limiting distribution',)
abline(h=upper_pi2, col="blue")
abline(h=lower_pi2, col="blue")

# Task 1 e)
set.seed(1)
N <- 1000
n = 300
Y0 <- c(950, 50, 0)

Y_n_simulate <- function(alfa, gamma, N, n, Y0){
  Y <- list()
  Y[[1]] <- Y0
  
  for(i in 2:n){
    beta_n <- 0.5*Y[[i-1]][2]/N
    syke_syke <- rbinom(1, Y[[i-1]][2], 1-gamma)
    sus_syke <- rbinom(1, Y[[i-1]][1], beta_n)
    syke_recover <- Y[[i-1]][2] - syke_syke 
    recover_recover <- rbinom(1, Y[[i-1]][3], 1-alfa)
    
    I_n <- syke_syke + sus_syke
    R_n <- syke_recover + recover_recover
    S_n <- N - I_n - R_n
    
    temp <- c(S_n, I_n, R_n)
    Y[[i]] <- temp
  }
  return(Y = Y)
}

Y_n <- Y_n_simulate(alfa, gamma, N, n, Y0)
Y_n

I_n <- rep(NA, n)
S_n <- rep(NA, n)
R_n <- rep(NA, n)

for(i in 1:n){
  S_n[i] <- Y_n[[i]][1]
  I_n[i] <- Y_n[[i]][2]
  R_n[i] <- Y_n[[i]][3]
}

plot(S_n,
     main = 'Infected (red), Recovered (green), Subceptible (blue) \n individuals in one realisation',
     xlab = 'time step #',
     ylab = '# of individuals',
     xlim = c(0, 300),
     ylim = c(0, 1000),
     col="blue",)
par(new = TRUE)
plot(I_n,
     main = '',
     xlab = '',
     ylab = '',
     yaxt="n",
     xaxt="n",
     xlim = c(0, 300),
     ylim = c(0, 1000),
     col="red",)
par(new = TRUE)
plot(R_n,
     main = '',
     xlab = '',
     ylab = '',
     yaxt="n",
     xaxt="n",
     xlim = c(0, 300),
     ylim = c(0, 1000),
     col="green")

# Task 1 f)
set.seed(1)
Sim_Num <- 1000
E_maxI <- rep(NA, Sim_Num) 
E_min_argmaxI <- rep(NA, Sim_Num)  
for(i in 1:Sim_Num){
  I <- rep(NA, n)
  res <- Y_n_simulate(alfa, gamma, N, n, Y0)
  for(j in 1:n){
    I[j] <- res[[j]][2]
  }
  E_maxI[i] <- max(I)
  E_min_argmaxI[i] <- which.max(I)
}
mean_E_maxI <- mean(E_maxI)
mean_E_min_argmaxI <- mean(E_min_argmaxI)

# upper bound
mean_E_maxI + qnorm(0.975) * sd(E_maxI) / sqrt(Sim_Num)
# lower bound
mean_E_maxI - qnorm(0.975) * sd(E_maxI) / sqrt(Sim_Num)
# upper bound
mean_E_min_argmaxI + qnorm(0.975) * sd(E_min_argmaxI) / sqrt(Sim_Num)
# lower bound
mean_E_min_argmaxI - qnorm(0.975) * sd(E_min_argmaxI) / sqrt(Sim_Num)

# Task 1 g)

Y_n_simulate_Immune <- function(alfa, gamma, N, immune, n, Y0){
  Y <- list()
  Y[[1]] <- Y0
  
  for(i in 2:n){
    beta_n <- 0.5*Y[[i-1]][2]/N
    syke_syke <- rbinom(1, Y[[i-1]][2], 1-gamma)
    sus_syke <- rbinom(1, Y[[i-1]][1], beta_n)
    syke_recover <- Y[[i-1]][2] - syke_syke 
    recover_recover <- rbinom(1, Y[[i-1]][3]-immune, 1-alfa) + immune
    
    I_n <- syke_syke + sus_syke
    R_n <- syke_recover + recover_recover
    S_n <- N - I_n - R_n
    
    Y[[i]] <- c(S_n, I_n, R_n)
  }
  return(Y = Y)
}
n = 300
N <- 1000
Sim_Num <- 1000

# 100 immune
set.seed(1)
immune <- 100
Y0 <- c(950-immune, 50, 0+immune)

Y_n <- Y_n_simulate_Immune(alfa, gamma, N, immune, n, Y0)

I_n <- rep(NA, n)
S_n <- rep(NA, n)
R_n <- rep(NA, n)

for(i in 1:n){
  S_n[i] <- Y_n[[i]][1]
  I_n[i] <- Y_n[[i]][2]
  R_n[i] <- Y_n[[i]][3]
}

plot(S_n,
     main = 'Infected (red), Recovered (green), Subceptible (blue) \n individuals in one realisation \n with 100 immune',
     xlab = 'time step #',
     ylab = '# of individuals',
     xlim = c(0, 300),
     ylim = c(0, 1000),
     col="blue",)
par(new = TRUE)
plot(I_n,
     main = '',
     xlab = '',
     ylab = '',
     yaxt="n",
     xaxt="n",
     xlim = c(0, 300),
     ylim = c(0, 1000),
     col="red",)
par(new = TRUE)
plot(R_n,
     main = '',
     xlab = '',
     ylab = '',
     yaxt="n",
     xaxt="n",
     xlim = c(0, 300),
     ylim = c(0, 1000),
     col="green")

E_maxI_100 <- rep(NA, Sim_Num) 
E_min_argmaxI_100 <- rep(NA, Sim_Num)  
for(i in 1:Sim_Num){
  I <- rep(NA, n)
  res <- Y_n_simulate_Immune(alfa, gamma, N, immune, n, Y0)
  for(j in 1:n){
    I[j] <- res[[j]][2]
  }
  E_maxI_100[i] <- max(I)
  E_min_argmaxI_100[i] <- which.max(I)
}

mean_E_maxI_100 <- mean(E_maxI_100)
mean_E_min_argmaxI_100 <- mean(E_min_argmaxI_100)

# upper bound
mean_E_maxI_100 + qnorm(0.975) * sd(E_maxI_100) / sqrt(Sim_Num)
# lower bound
mean_E_maxI_100 - qnorm(0.975) * sd(E_maxI_100) / sqrt(Sim_Num)
# upper bound
mean_E_min_argmaxI_100 + qnorm(0.975) * sd(E_min_argmaxI_100) / sqrt(Sim_Num)
# lower bound
mean_E_min_argmaxI_100 - qnorm(0.975) * sd(E_min_argmaxI_100) / sqrt(Sim_Num)

# 600 immune
set.seed(1)
immune <- 600
Y0 <- c(950-immune, 50, 0+immune)

Y_n <- Y_n_simulate_Immune(alfa, gamma, N, immune, n, Y0)

I_n <- rep(NA, n)
S_n <- rep(NA, n)
R_n <- rep(NA, n)

for(i in 1:n){
  S_n[i] <- Y_n[[i]][1]
  I_n[i] <- Y_n[[i]][2]
  R_n[i] <- Y_n[[i]][3]
}

plot(S_n,
     main = 'Infected (red), Recovered (green), Subceptible (blue) \n individuals in one realisation \n with 600 immune',
     xlab = 'time step #',
     ylab = '# of individuals',
     xlim = c(0, 300),
     ylim = c(0, 1000),
     col="blue",)
par(new = TRUE)
plot(I_n,
     main = '',
     xlab = '',
     ylab = '',
     yaxt="n",
     xaxt="n",
     xlim = c(0, 300),
     ylim = c(0, 1000),
     col="red",)
par(new = TRUE)
plot(R_n,
     main = '',
     xlab = '',
     ylab = '',
     yaxt="n",
     xaxt="n",
     xlim = c(0, 300),
     ylim = c(0, 1000),
     col="green")

E_maxI_600 <- rep(NA, Sim_Num) 
E_min_argmaxI_600 <- rep(NA, Sim_Num)  
for(i in 1:Sim_Num){
  I <- rep(NA, n)
  res <- Y_n_simulate_Immune(alfa, gamma, N, immune, n, Y0)
  for(j in 1:n){
    I[j] <- res[[j]][2]
  }
  E_maxI_600[i] <- max(I)
  E_min_argmaxI_600[i] <- which.max(I)
}

mean_E_maxI_600 <- mean(E_maxI_600)
mean_E_min_argmaxI_600 <- mean(E_min_argmaxI_600)

# upper bound
mean_E_maxI_600 + qnorm(0.975) * sd(E_maxI_600) / sqrt(Sim_Num)
# lower bound
mean_E_maxI_600 - qnorm(0.975) * sd(E_maxI_600) / sqrt(Sim_Num)
# upper bound
mean_E_min_argmaxI_600 + qnorm(0.975) * sd(E_min_argmaxI_600) / sqrt(Sim_Num)
# lower bound
mean_E_min_argmaxI_600 - qnorm(0.975) * sd(E_min_argmaxI_600) / sqrt(Sim_Num)

# 800 immune
set.seed(1)
immune <- 800
Y0 <- c(950-immune, 50, 0+immune)

Y_n <- Y_n_simulate_Immune(alfa, gamma, N, immune, n, Y0)

I_n <- rep(NA, n)
S_n <- rep(NA, n)
R_n <- rep(NA, n)

for(i in 1:n){
  S_n[i] <- Y_n[[i]][1]
  I_n[i] <- Y_n[[i]][2]
  R_n[i] <- Y_n[[i]][3]
}

plot(S_n,
     main = 'Infected (red), Recovered (green), Subceptible (blue) \n individuals in one realisation \n with 800 immune',
     xlab = 'time step #',
     ylab = '# of individuals',
     xlim = c(0, 300),
     ylim = c(0, 1000),
     col="blue",)
par(new = TRUE)
plot(I_n,
     main = '',
     xlab = '',
     ylab = '',
     yaxt="n",
     xaxt="n",
     xlim = c(0, 300),
     ylim = c(0, 1000),
     col="red",)
par(new = TRUE)
plot(R_n,
     main = '',
     xlab = '',
     ylab = '',
     yaxt="n",
     xaxt="n",
     xlim = c(0, 300),
     ylim = c(0, 1000),
     col="green")

E_maxI_800 <- rep(NA, Sim_Num) 
E_min_argmaxI_800 <- rep(NA, Sim_Num)  
for(i in 1:Sim_Num){
  I <- rep(NA, n)
  res <- Y_n_simulate_Immune(alfa, gamma, N, immune, n, Y0)
  for(j in 1:n){
    I[j] <- res[[j]][2]
  }
  E_maxI_800[i] <- max(I)
  E_min_argmaxI_800[i] <- which.max(I)
}

mean_E_maxI_800 <- mean(E_maxI_800)
mean_E_min_argmaxI_800 <- mean(E_min_argmaxI_800)

# upper bound
mean_E_maxI_800 + qnorm(0.975) * sd(E_maxI_800) / sqrt(Sim_Num)
# lower bound
mean_E_maxI_800 - qnorm(0.975) * sd(E_maxI_800) / sqrt(Sim_Num)
# upper bound
mean_E_min_argmaxI_800 + qnorm(0.975) * sd(E_min_argmaxI_800) / sqrt(Sim_Num)
# lower bound
mean_E_min_argmaxI_800 - qnorm(0.975) * sd(E_min_argmaxI_800) / sqrt(Sim_Num)



# Problem 2

set.seed(1)

rate = 1.5
t = 59

N = 1000

num_insurance = function(rate,t) {
  n <- rpois(1,rate*t)
  
  X = matrix(0,nrow = 2, ncol = n+2)
  
  X[2,2:(n+1)] <- 1:n
  X[2,n+2] <- n
  
  X[1,0] <- 0
  X[1,n+2] <- t
  X[1,2:(n+1)] <- sort(runif(n,min = 0,max = t))
  
  return(X)
}


colors_plot = c('#009999','#0000FF','#FF3399','#000000','#CC0CC0','#00FFFF','#00CC99','#FF9933','#999999','#00C00C')

claims = vector('numeric',length = N)


plot(0,0,xlim = c(0,t), ylim = c(0,120), ylab = "Number of claims", xlab = "Time in days",
     main = "10 realizations of insurance claims")
for (i in 1:N) {
  
  X <- num_insurance(rate,t)
  
  claims[i] <- X[2,ncol(X)]
  
  if (i <= 10) {
    segments(X[1,1:(ncol(X)-1)],X[2,1:(ncol(X)-1)],X[1,2:ncol(X)],X[2,1:(ncol(X)-1)], col = colors_plot[i])
  }
}

# Printer estimering for sannsyneligheten at antall claims > 100
print(sum(claims > 100)/N)


ratem <- 10

monetary_insurance = function(rate,t,ratem) {
  X <- num_insurance(rate,t)
  
  Z <- matrix(0,nrow = nrow(X),ncol = ncol(X))
  
  Z[1,] <- X[1,]
  
  for (i in 2:ncol(Z)) {
    Z[2,i] = sum(rexp(X[2,i]-X[2,i-1],ratem)) + Z[2,i-1]
  }
  
  return(Z)
}

monetary_claims = vector('numeric',length = N)


plot(NULL,NULL,xlim = c(0,t), ylim = c(0,13), xlab = "Time in days", ylab = "Total claim amount in mill. kr",
     main = "10 realizations of claim amount.")
for (i in 1:N) {
  
  Z <- monetary_insurance(rate,t,ratem)
  
  monetary_claims[i] <- Z[2,ncol(Z)]
  
  if (i <= 10) {
    segments(Z[1,1:(ncol(Z)-1)],Z[2,1:(ncol(Z)-1)],Z[1,2:ncol(Z)],Z[2,1:(ncol(Z)-1)], col = colors_plot[i])
  }
}

print(sum(monetary_claims > 8)/N)