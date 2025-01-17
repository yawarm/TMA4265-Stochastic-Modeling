# Task 1 b)

# i)

lam <- 5 # arrival rate [patients per hour]
mu <- 1/10 * 60 # service rate [per hour]
N <- 1200 # total duration of simulation [hour]
set.seed(24)

MM1_simulate <- function(N, lam, mu) {
  x_0 <- 0 # initial condition
  t_0 <- 0 # initial time in condition
  x <- c(x_0) # vector of states
  s <- c(t_0) # vector of transition times
  while( max(s) < N ) {
    state = tail(x, 1) # current state
    if (state == 0) {
      s <- c(s, tail(s, 1) + rexp(1, lam))
      x <- c(x, 1)
    } else {
      s <- c(s, tail(s,1) + rexp(1, lam+mu))
      if (runif(1) < lam/(lam + mu)) {
        x <- c(x, state + 1)
      } else {
        x <- c(x, state - 1)
      }
    }
  }
  res <- list(x, s)
  return(res)
}
# ii)

W_UCC <- function(x, s) {
  
  numStates <- max(x) # number of states
  
  #time spent at each state [hour]
  numT <- rep(0, numStates + 1)
  for(i in 1:(length(x) - 1)) {
    numT[x[i] + 1] <- numT[x[i]+1]+s[i+1]-s[i]
  }
  
  pi <- numT/s[length(x)-1] # long term proportion of time spent in each state
  
  L <- sum((0:numStates)*pi) # average number of patients in the system
  
  W <- L/lam # average time spent in the system [hours]
  
  return(W)
}

# iii)
n <- 30
W_vec <- c()
for (i in 1:n) {
  res <- MM1_simulate(N, lam, mu)
  x <- res[[1]]
  s <- res[[2]]
  W <- W_UCC(x, s)
  print(W)
  W_vec <- c(W_vec, W)
}

# upper bound
upper_bound <- mean(W_vec) + qt(0.975, length(W_vec)-1) * sd(W_vec) / sqrt(length(W_vec))
# lower bound
lower_bound <- mean(W_vec) - qt(0.975, length(W_vec)-1) * sd(W_vec) / sqrt(length(W_vec))


# iv)

res <- MM1_simulate(N, lam, mu)
x <- res[[1]]
s <- res[[2]]

index = match(12, floor(s)) # index of the value where 12 hours have passed

plot(NULL,
     NULL,
     xlim = c(0, 12),
     ylim = c(0, 12),
     xlab = "Time (h)",
     ylab = "Number of patients",
     cex.axis = 1.5, cex.lab = 1.5,
     main = 'Number of patients in the UCC as a function of time')
for(i in 1:index){
  lines(s[i:(i+1)], rep(x[i],2), lwd = 2)
}




# g)

# i

lam <- 5 # arrival rate [patients per hour]
mu <- 1/10 * 60 # service rate [per hour]
Num <- 1200 # total duration of simulation [hour]
p <- 0.8 # probability that the patient is urgent
set.seed(24)

U_and_N_simulate <- function(Num, lam, mu) {
  t_0 <- 0 # initial condition time
  x_0 <- 0 # initial condition state
  U <- c(x_0) # vector to hold urgent patients
  U_s <- c(t_0) # vector of transition times for urgent patients
  N_s <- c(t_0) # vector of transition times for normal patients
  N <- c(x_0) # vector to hold normal patients
  s <- c(t_0) # total time

  while( tail(s,1) < Num ) {
      if (tail(U,1) > 0) {
        if (runif(1) > lam/(lam + mu)) {
          exp_u <- rexp(1, p*lam+mu)
          U_s <- c(U_s, tail(U_s,1) + exp_u)
          s <- c(s, tail(s,1) + exp_u)
          U <- c(U, tail(U,1) - 1)
        } else {
          if (runif(1) < p) {
            exp_u <- rexp(1, p*lam)
            U_s <- c(U_s, tail(U_s, 1) + exp_u)
            s <- c(s, tail(s,1) + exp_u)
            U <- c(U, tail(U,1) + 1)
          } else {
            exp_n <- rexp(1, (1-p)*lam)
            N_s <- c(N_s, tail(N_s, 1) + exp_n)
            s <- c(s, tail(s,1) + exp_n)
            N <- c(N, tail(N,1) + 1)
          }
        }
      } else if (tail(N,1) > 0) {
        if (runif(1) > lam/(lam + mu)) {
          exp_n <- rexp(1, (1-p)*lam+mu)
          N_s <- c(N_s, tail(N_s,1) + exp_n)
          s <- c(s, tail(s,1) + exp_n)
          N <- c(N, tail(N,1) - 1)
        } else {
          if (runif(1) < p) {
            exp_u <- rexp(1, p*lam)
            U_s <- c(U_s, tail(U_s, 1) + exp_u)
            s <- c(s, tail(s,1) + exp_u)
            U <- c(U, tail(U,1) + 1)
          } else {
            exp_n <- rexp(1, (1-p)*lam)
            N_s <- c(N_s, tail(N_s, 1) + exp_n)
            s <- c(s, tail(s,1) + exp_n)
            N <- c(N, tail(N,1) + 1)
          }
        }
      } else {
        if (runif(1) < p) {
          exp_u <- rexp(1, p*lam)
          U_s <- c(U_s, tail(U_s, 1) + exp_u)
          s <- c(s, tail(s,1) + exp_u)
          U <- c(U, tail(U,1) + 1)
        } else {
          exp_n <- rexp(1, (1-p)*lam)
          N_s <- c(N_s, tail(N_s, 1) + exp_n)
          s <- c(s, tail(s,1) + exp_n)
          N <- c(N, tail(N,1) + 1)
        }
      }
  }
  res <- list(U, N, U_s, N_s, s)
  return(res)
}
  
res <- U_and_N_simulate(Num, lam, mu)

U <- res[[1]]
N <- res[[2]]
U_s <- res[[3]]
N_s <- res[[4]]
s <- res[[5]]

# ii

W_N <- function(N, N_s) {
  
  numStates <- max(N) # number of states
  
  #time spent at each state [hour]
  numT <- rep(0, numStates + 1)
  for(i in 1:(length(N) - 1)) {
    numT[N[i] + 1] <- numT[N[i]+1]+N_s[i+1]-N_s[i]
  }
  
  pi <- numT/N_s[length(N)] # long term proportion of time spent in each state
  
  L_N <- sum((0:numStates)*pi) # average number of patients in the system
  
  W_N <- L_N/((1-p)*lam) # average time spent in the system [hours]
  
  return(W_N)
}

n <- 30
W_N_vec <- c()
for (i in 1:n) {
  res <- U_and_N_simulate(Num, lam, mu)
  U <- res[[1]]
  N <- res[[2]]
  U_s <- res[[3]]
  N_s <- res[[4]]
  W <- W_N(N, N_s)
  W_N_vec <- c(W_N_vec, W)
}

# upper bound
upper_bound_W_N <- mean(W_N_vec) + qt(0.975, length(W_N_vec)-1) * sd(W_N_vec) / sqrt(length(W_N_vec))
# lower bound
lower_bound_W_N <- mean(W_N_vec) - qt(0.975, length(W_N_vec)-1) * sd(W_N_vec) / sqrt(length(W_N_vec))


W_U <- function(U, U_s) {
  
  numStates <- max(U) # number of states
  
  #time spent at each state [hour]
  numT <- rep(0, numStates + 1)
  for(i in 1:(length(U) - 1)) {
    numT[U[i] + 1] <- numT[U[i]+1]+U_s[i+1]-U_s[i]
  }
  
  pi <- numT/U_s[length(U)] # long term proportion of time spent in each state
  
  L_U <- sum((0:numStates)*pi) # average number of patients in the system
  
  W_U <- L_U/(p*lam) # average time spent in the system [hours]
  
  return(W_U)
}

n <- 30
W_U_vec <- c()
for (i in 1:n) {
  res <- U_and_N_simulate(Num, lam, mu)
  U <- res[[1]]
  N <- res[[2]]
  U_s <- res[[3]]
  N_s <- res[[4]]
  W <- W_U(U, U_s)
  W_U_vec <- c(W_U_vec, W)
}

# upper bound
upper_bound_W_U <- mean(W_U_vec) + qt(0.975, length(W_U_vec)-1) * sd(W_U_vec) / sqrt(length(W_U_vec))
# lower bound
lower_bound_W_U <- mean(W_U_vec) - qt(0.975, length(W_U_vec)-1) * sd(W_U_vec) / sqrt(length(W_U_vec))



# iii)

res <- U_and_N_simulate(Num, lam, mu)
U <- res[[1]]
N <- res[[2]]
U_s <- res[[3]]
N_s <- res[[4]]
s <- res[[5]]

plot(NULL,
     NULL,
     xlim = c(0, 12),
     ylim = c(0, 12),
     xlab = "Time (h)",
     ylab = "Number of patients",
     cex.axis = 1.5, cex.lab = 1.5,
     main = 'Number of Urgent patients in the UCC as a function of time')
for(i in 1:index){
  lines(s[i:(i+1)], rep(U[i],2), lwd = 2, col = 'blue')
}
par(new = TRUE)
plot(NULL,
     NULL,
     col = 'green', 
     xlim = c(0, 12),
     ylim = c(0, 12),
     xlab = "Time (h)",
     ylab = "Number of patients",
     cex.axis = 1.5, cex.lab = 1.5,
     main = 'Number of Urgent patients in the UCC as a function of time')
for(i in 1:index){
  lines(s[i:(i+1)], rep(N[i],2), lwd = 2, col = 'green')
}

# Task 2

set.seed(1)

# a)

weather_points <- c(0.3,0.35,0.39,0.41,0.45)
weather_values <- c(0.5,0.32,0.4,0.35,0.6)

n <- 51
theta_vec <- round(seq(from = 0.25, to = 0.5, by = 0.005),3)

weather_index = which(theta_vec %in% weather_points)

mu = rep(0.5,n)

sigma2 = 0.5**2
Cov = matrix(0,nrow = n,ncol = n)

corr = function(theta1,theta2) {return((1+15*abs(theta1-theta2))*exp(-15*abs(theta1-theta2)))}

for (i in 1:n) {
  for (j in 1:(i)) {
    Cov[i,j] = sigma2*corr(theta_vec[i],theta_vec[j])
    Cov[j,i] = Cov[i,j]
  } 
}

CovC = Cov-Cov[,weather_index, drop = FALSE]%*%solve(Cov[weather_index,weather_index, drop = FALSE], Cov[weather_index,,drop = FALSE])
muC = mu + Cov[,weather_index,drop = FALSE]%*%solve(Cov[weather_index,weather_index,drop = FALSE], weather_values-rep(0.5,5))

y_vec = rnorm(n,muC,diag(CovC))

plot(theta_vec,y_vec, main = "Weather value prediction \n 5 measurments", xlab = "theta", ylab = "y(theta)",
     ylim = c(0,1.5))
lines(theta_vec, muC+1.64*sqrt(diag(CovC)))
lines(theta_vec, muC-1.64*sqrt(diag(CovC)))

#b)

Y_vec = pnorm(rep(0.3,n),mean = muC,sd = sqrt(diag(CovC)))

plot(theta_vec,Y_vec, main = "Pr(Y(theta)<0.3)\n 5 measurments", xlab = "Theta", ylab = "Pr(Y(theta)<0.3)")


#c)

weather_points <- c(weather_points,0.33)
weather_values <- c(weather_values,0.4)

weather_index = which(theta_vec %in% weather_points)

CovC = Cov-Cov[,weather_index, drop = FALSE]%*%solve(Cov[weather_index,weather_index, drop = FALSE], Cov[weather_index,,drop = FALSE])
muC = mu + Cov[,weather_index,drop = FALSE]%*%solve(Cov[weather_index,weather_index,drop = FALSE], weather_values-rep(0.5,6))

y_vec = rnorm(n,muC,diag(CovC))

plot(theta_vec,y_vec, main = "Weather value prediction\n 6 measurements", xlab = "theta", ylab = "y(theta)",
     ylim = c(0,1))
lines(theta_vec, muC+1.64*sqrt(diag(CovC)))
lines(theta_vec, muC-1.64*sqrt(diag(CovC)))

Y_vec = pnorm(rep(0.3,n),mean = muC,sd = sqrt(diag(CovC)))

plot(theta_vec,Y_vec, main = "Pr(Y(theta)<0.3) \n 6 measurments", xlab = "Theta", ylab = "Pr(Y(theta)<0.3)")