#################################
#Epsilon-strong with rough paths
#################################

################
#1. Algorithm 1
################

#maps integer l>=1 to a pair (n,k) with l = 2^(n-1) + k; 1 -> (1,0)
unzip <- function(l){

  n <- floor(log2(l)) + 1
  k <- l - 2^(n-1)

  out <- matrix(c(n,k), ncol = 2)
  return(out)
}

#samples (n,k)-indices of all record breakers of first type up to final one
record_breakers1 <- function(){

  R <- 0
  S <- c()

  while(TRUE){
    D <- 0; U <- 1
    V <- runif(1)

    while((D < V) & (V < U)){
      R <- R+1
      n <- floor(log2(R)) + 1
      U <- (2*pnorm(4 * sqrt(n+1)) - 1) * U
      D <- (1 - R^(-7)) * U
    }

    if(V > U){
      S <- c(S,R)
    } else if (V < D){
      if(length(S) == 0){
        pairs <- matrix(0,ncol=2,nrow=0)
      } else {
        pairs <- unzip(S)
      }
      out <- pairs
      return(out)
    }
  }
}

#sample from tail of gaussian
normal_tail_sample <- function(mu_min){

  while(TRUE){
    alpha <- 1/2 * (mu_min + sqrt(mu_min^2 + 4))

    z_prop <- rexp(1,rate = alpha) + mu_min
    rho <- exp(-1/2 * (z_prop - alpha)^2)
    U <- runif(1)

    if(U < rho){
      return(z_prop)
    }
  }
}

#sample outside of tail of gaussian
normal_body_sample <- function(mu_plus){

  while(TRUE){
    prop <- rnorm(1)
    if(prop < mu_plus){
      return(prop)
    }
  }
}

#samples Brownian motion conditional on RB1 pairs up to given resolution n;
#also tracks Vm = max|Wm| and gives this as vector
ZV_sample <- function(n,pairs){

  Z <- c(0,rnorm(1))
  V <- abs(Z[2])

  #if no record breakers, immediately return
  if(n==0){
    return(list(Z=Z,V=V))
  } else {
    r <- nrow(pairs)
    current_pair <- 1
    for(m in 1:n){
      W_new <- numeric(2^(m-1))
      #CHECK SCALE
      scale <- 2^(-(m+1)/2)

      for(k in 1:2^(m-1)){
        if(current_pair <= r){
          m_curr <- pairs[current_pair,1]; k_curr <- pairs[current_pair,2]
          if((m==m_curr) & (k==k_curr)){
            W_new[k] <- normal_tail_sample(4*sqrt(m+1))
            current_pair <- current_pair+1
          } else {
            W_new[k] <- normal_body_sample(4*sqrt(m+1))
          }
        } else {
          W_new[k] <- normal_body_sample(4*sqrt(m+1))
        }
      }

      V <- c(V,max(abs(W_new)))
      Z_midpoints <- apply(matrix(c(Z[1:2^(m-1)],Z[2:(2^(m-1)+1)]), ncol=2),1,mean)
      Z_new <- Z_midpoints + scale*W_new
      Z <-  c(rbind(Z[1:2^(m-1)], Z_new), Z[2^(m-1)+1])
    }
    return(list(Z=Z,V=V))
  }
}

###################
#2. Procedure A
###################

#################################
#2a) Finding continuity constants: K_alpha, K_2alpha, Gamma_R
#################################

#########
#K_alpha
#########
C_K_alpha <- function(n,alpha){
  out <- 2^(-n/2 * (1/2 - alpha))*sqrt(n+1)
  return(out)
}

C_max <- function(N1, alpha){

  n1 <- floor(2/(1-2*alpha)*1/log(2) - 1)

  if(N1 >= n1){
    out <- C_func(N1+1,alpha)
  } else {
    C1 <- C_K_alpha(n1,alpha); C2 <- C_K_alpha(n1+1,alpha)
    if(C1>=C2){
      return(C1)
    } else {
      return(C2)
    }
  }
}

K_alpha <- function(N_1, alpha, V){

  scales <- 2^(-(1/2-alpha)*0:N_1)
  t1 <- 2^(2*alpha+1)*sum(scales*V)

  C <- C_max(N_1,alpha)
  num <- 2^(-1/2 * (N_1+1)*(1/2-alpha))
  denom <- 1 - 2^(-1/2 * (1/2 - alpha))
  t2 <- 2^(2*alpha + 3)*C*num/denom

  out <- t1+t2
  return(out)

}

###############
#Gamma_L
###############

zoom_out_once <- function(Z){
  n <- log2(ncol(Z)-1)
  indices <- seq(1,2^n+1,by=2)
  out <- Z[,indices]
  return(out)
}

#given multidim Z at resolution n, recover resolution m<=n
zoom_out <- function(Z,m){
  n <- log2(ncol(Z)-1)
  if(m > n){
    return("must have m <= n")
  } else if(m == n){
    return(Z)
  } else {
    indices <- seq(1, 2^n+1, by=2^(n-m))
    out <- Z[,indices]
    return(out)
  }
}

#find for given Z,i,j, max (L(m)-L(l))/((m-l)^beta * delta_n^2alpha)
Gamma_Ln_ij <- function(Z_i,Z_j,alpha,beta){

  n <- log2(length(Z_i)-1)

  Z11 <- Z_i[2*(1:2^(n-1))-1]
  Z12 <- Z_i[2*(1:2^(n-1))]
  Lambda_i <- Z12-Z11

  Z21 <- Z_j[2*(1:2^(n-1))]
  Z22 <- Z_j[2*(1:2^(n-1))+1]
  Lambda_j <- Z22 - Z21

  Lambda_prod <- Lambda_i*Lambda_j

  scale <- 2^(-n)
  denoms <- ((1:2^(n-1))^beta)*(scale^(2*alpha))
  out <- 1

  for(m in 0:(2^(n-1)-1)){
    sum_vals <- Lambda_prod[(m+1):2^(n-1)]
    s <- abs(cumsum(sum_vals))
    print(s)
    scaled_Ldiffs <- s/denoms[1:length(sum_vals)]
    out <- max(out,scaled_Ldiffs)
  }
  return(out)
}

Gamma_Ln_ij2 <- function(Z_i,Z_j,alpha,beta){

  n <- log2(length(Z_i)-1)

  Z11 <- Z_i[2*(1:2^(n-1))-1]
  Z12 <- Z_i[2*(1:2^(n-1))]
  Lambda_i <- Z12-Z11

  Z21 <- Z_j[2*(1:2^(n-1))]
  Z22 <- Z_j[2*(1:2^(n-1))+1]
  Lambda_j <- Z22 - Z21

  Lambda_prod <- Lambda_i*Lambda_j

  scale <- 2^(-n)
  denoms <- ((1:2^(n-1))^beta)*(scale^(2*alpha))

  tmp1 <- lapply(2^(n-1):1, function(i) cumsum(Lambda_prod[i:2^(n-1)]))
  nums <- abs(mapply(function(y, l) c(rep(NA, 2^(n-1)-l), y), tmp1, lengths(tmp1)))

  tmp2 <- matrix(NA,2^(n-1),2^(n-1))
  denoms <- sequence((2^(n-1):1))^beta * (scale^(2*alpha))
  tmp2[lower.tri(tmp2,diag=TRUE)] <- denoms

  out <- max(nums/tmp2, na.rm = TRUE)

  return(out)
}

#find Gamma_Ln for multidim Z
Gamma_Ln <- function(Z,alpha,beta){
  dprime <- nrow(Z)
  out <- 1
  for(i in 1:(dprime-1)){
    Z_i <- Z[i,]
    for(j in (i+1):dprime){
      Z_j <- Z[j,]
      candidate <- Gamma_Ln_ij(Z_i,Z_j,alpha,beta)
      out <- max(out, candidate)
    }
  }
  return(out)
}

Gamma_L <- function(Z,alpha,beta){
  n <- log2(ncol(Z)-1)
  Gamma_L <- 1
  for(i in 1:(n-1)){
    print(i)
    Z <- zoom_out_once(Z)
    candidate <- Gamma_Ln(Z,alpha,beta)
    Gamma_L <- max(Gamma_L,candidate)
  }
  return(Gamma_L)
}

Gamma_L_partial <- function(Z,m,alpha,beta){
  Gamma_L <- 1
  for(i in 1:m){
    Z <- zoom_out_once(Z)
    candidate <- Gamma_Ln(Z,alpha,beta)
    Gamma_L <- max(Gamma_L,candidate)
  }
  return(Gamma_L)
}

##############
#Gamma_R
##############
Gamma_R <- function(Gamma_L, alpha, beta){

  num <- 2^(-(2*alpha - beta))
  denom <- 1 - num
  out <- num/denom * Gamma_L
  return(out)
}

##############
#K_2alpha
##############
K_2alpha <- function(K_alpha, alpha, Gamma_R){

  l <- 2*Gamma_R/(1 - 2^(-2*alpha))
  r <- K_alpha^2 * 2^(1-alpha)/(1-2^(-alpha))
  out <- l+r
  return(out)
}

##############
#2b) find G given continuity constants
##############
find_C <- function(alpha, K_alpha, K_2alpha, K_R, M, d, dprime){

  d <- max(d, dprime)

  C1_delta <- d*M*K_alpha + 1/2
  C2_delta <- (d^3)*(M^2)*K_2alpha + 1/2
  scale <- 2/(1 - 2^(1-3*alpha))
  C3_delta <- scale*(M*C1_delta + d*M*(C1_delta^2)*K_alpha + (d^2)*M*C2_delta*K_alpha + (d^2)*(M^2)*K_alpha + 2*(d^3)*(M^2)*C1_delta*K_2alpha)

  print(c(C1_delta,C2_delta,C3_delta))

  delta <- 1
  cv1 <- 1; cv2 <- 1

  while((cv1 >= 1/2) | (cv2 >= 1/2)){

    delta <- delta/2
    cv1 <- C3_delta*(delta^(2*alpha)) + M*(delta^(1-alpha)) + (d^3)*(M^2)*K_2alpha*(delta^alpha)
    cv2 <- C3_delta*(delta^alpha)
  }

  C1 <- 2/delta * C1_delta
  C2 <- 2/delta * (C2_delta + M*C1 + d*M*C1*K_alpha)
  C3 <- scale * (M*C1 + d*M*(C1^2)*K_alpha + (d^2)*M*C2*K_alpha + 2*(d^3)*(M^2)*C1*K_2alpha)

  out <- list(C1=C1, C2=C2, C3=C3)
  return(out)
}

find_B <- function(alpha, K_alpha, K_2alpha, K_R, M){

  Bs <- find_C(alpha, K_alpha, K_2alpha, K_R, 2*M, 1, 1)
  B <- Bs$C1
  return(B)
}

find_G <- function(alpha, beta, K_alpha, K_2alpha, K_R, M, d, dprime){

  Cs <- find_C(alpha, K_alpha, K_2alpha, K_R, M, d, dprime)
  C1 <- Cs$C1; C3 <- Cs$C3
  B <- find_B(alpha, K_alpha, K_2alpha, K_R, M)

  G1 <- (1+B)*C3

  delta <- 1
  threshold <- 2^(alpha + beta) - 2
  while(B*delta > threshold){
    delta <- delta/2
  }

  scale <- 1 - (2 + B*(delta^alpha))/(2^(alpha+beta))
  C4_delta <- 2/scale * (B*(d^3)*(M^2)*(K_R) + 2*(d^3)*(M^2)*C1*K_R)
  C4 <- (1 + B*(delta^alpha))*C4_delta + 2/delta * (B*(d^3)*(M^2)*K_R + 2*(d^3)*(M^2)*C1*K_R)
  G2 <- C4 + (d^3)*(M^2)*K_R
  G <- G1 + G2
  return(G)
}

########################
#3. Procedure Aux & Procedure B
########################

##############
#3a) Recursions
##############

theta_plus<- function(n,m,l,r,k,kprime,theta){

  if(l==1){
    if((k < r) & (r <= kprime)){
      return(theta)
    } else {
      return(0)
    }
  } else {
    out <- theta_plus(n,m,l-1,2*r-1,k,kprime,theta)+theta_plus(n,m,l-1,2*r,k,kprime,theta)
    return(out)
  }
}

theta_minus<- function(n,m,l,r,k,kprime,theta){

  if(l==1){
    if((k < r) & (r <= kprime)){
      return(theta)
    } else {
      return(0)
    }
  } else {
    out <- theta_minus(n,m,l-1,2*r-1,k,kprime,theta)-theta_minus(n,m,l-1,2*r,k,kprime,theta)
    return(out)
  }
}

#note: eta_plus/minus(l, eta1) is found using theta_plus/minus(l, eta1)

rho_fn <- function(n,m,l,tp,ep){

  delta <- 2^(-(n+m-l+2))
  num <- delta*tp
  denom <- 1 - 2*delta*ep
  out <- num/denom
  return(out)
}

g_fn <- function(n,m,l,tp,ep){

  delta <- 2^(-(n+m-l+2))
  num <- delta*tp
  denom <- 1 - 2*delta*ep
  out <- num/(denom^2)
  return(out)
}

h_fn <- function(n,m,l,ep,rho){

  delta <- 2^(-(n+m-l+2))
  num <- delta
  denom <- (1 - 2*delta*ep)*(1-rho^2)
  out <- num/denom
  return(out)
}

eta_fn <- function(tm,ep,em,rho,h){

  out <- ep/4 + h/8 * (tm^2 + 4*em^2 + 4*tm*em*rho)
  return(out)
}

theta_fn <- function(tp,tm,em,rho,g,h){

  out <- tp/4 + h*(tm*em + (tm^2)*g/4 + (em^2)*rho)
  return(out)
}

C_fn <- function(n,m,l,r,k,kprime,theta_1,eta_1){

  delta <- 2^(-(n+m-l+1))
  tp <- theta_plus(n,m,l,r,k,kprime,theta_1)
  ep <- theta_plus(n,m,l,r,k,kprime,eta_1)
  rho <- rho_fn(n,m,l,tp,ep)

  denom <- (1 - 2*delta*ep)*sqrt(1-rho^2)
  out <- 1/denom
  return(out)
}

###############
#3b) Procdure Aux
###############

############
#Step 1: trivial
############

############
#Step 2: sample k, kprime
############
aux_2 <- function(n,m,alpha,alpha_prime,beta,gamma){

  klim <- 2^(n+m-2)

  C1 <- 2^(-2*n*(alpha - alpha_prime))
  C2 <- 2^(-m*(2*alpha - 1))
  v <- 8*exp(-gamma/2*C1*C2*(1:klim)^(beta - 1/2))
  multiplicities <- klim:1

  v2 <- v*multiplicities
  q <- v2/sum(v2)
  diff <- sample(1:klim, 1, prob = q)
  k <- sample(0:(klim-diff),1)
  kprime <- k+diff
  q_val <- v[diff]/sum(v2)

  out <- list(k=k,kprime=kprime,q=q_val)
  return(out)
}

#############
#Step 3: sample proposal Z from tilted distribution
#############

#sampling from P'
#sample from distribution defined in Lem. 5.1
five_one_sample <- function(n,m,Lambda_i,Lambda_j,theta0){

  delta <- 2^(-(n+m+1))

  cv12 <- -theta0*delta
  scale <- 1/(1-(theta0^2)*(delta^2))
  cov <- scale*matrix(c(1,cv12,cv12,1), nrow=2)

  m <- c(Lambda_j,-Lambda_i)
  mean <- (1/2 * theta0*sqrt(delta))*(cov %*% m)

  out <- mvrnorm(1, mu = mean, Sigma = cov)
  return(out)
}

#sample from distribution in Lem. 5.2
five_two_sample <- function(n,m,l,r,k,kprime,Lambda_i,Lambda_j,theta,eta){

  delta <- 2^(-(n+m-l+1))
  tp <- theta_plus(n,m,l,r,k,kprime,theta)
  tm <- theta_minus(n,m,l,r,k,kprime,theta)
  ep <- theta_plus(n,m,l,r,k,kprime,eta)
  em <- theta_minus(n,m,l,r,k,kprime,eta)
  rho <- rho_fn(n,m,l,tp,ep)
  g <- g_fn(n,m,l,tp,ep)

  cv11 <- 1/(1 - 2*delta*ep);
  cv12 <- g; cv21 <- g
  cv22 <- cv11
  scale <- 1/(1-rho^2)

  cov <- scale*matrix(c(cv11,cv21,cv12,cv22), ncol = 2)

  m1 <- Lambda_i*em + 1/2 * Lambda_j*tm
  m2 <- Lambda_j*em + 1/2 * Lambda_i*tm
  m <- c(m1,m2)
  mean <- sqrt(delta)*(cov %*% m)

  out <- mvrnorm(1, mu = mean, Sigma = cov)
  return(out)
}

#Given Z^n and proposed W_new, find Z^(n+1)
BM_fill <- function(Z, W){

  n <- log2(ncol(W))
  dprime <- nrow(W)
  delta <- 2^(-(n+2))

  Z_midpoints <- 1/2 * (Z[,1:2^n] + Z[,2:(2^n+1)])
  Z_new <- Z_midpoints + sqrt(delta)*W

  out <- matrix(0,nrow=dprime,ncol=2^(n+1)+1)
  out[,seq(1,2^(n+1)+1,by=2)] <- Z
  out[,seq(2,2^(n+1),by=2)] <- Z_new
  return(out)
}
#check if an RB2 occurs
#COULD IMPROVE THIS - HALT AS SOON AS A RB2 IS FOUND, GAMMA CHECKS ALL FIRST
RB2_check <- function(Z,alpha,beta){

  candidate <- Gamma_Ln(Z,alpha,beta)
  if(candidate > 1){
    return(1)
  } else {
    return(0)
  }
}

#Resolve final step
aux_32 <- function(Z,n,m,i,j,k,kprime,theta_0){

  dprime <- nrow(Z)
  W <- matrix(rnorm(dprime*(2^(n+m-1))), nrow = dprime)

  Z_i1 <- Z[i,1:2^(n+m-1)]; Z_i2 <- Z[i,2:(2^(n+m-1)+1)]
  Z_j1 <- Z[j,1:2^(n+m-1)]; Z_j2 <- Z[j,2:(2^(n+m-1)+1)]
  Lambda_i <- Z_i2 - Z_i1; Lambda_j <- Z_j2 - Z_j1

  for(r in 1:(2^(n+m-1))){
    Lambda_ir <- Lambda_i[r];Lambda_jr <- Lambda_j[r]
    W[c(i,j),r] <- five_one_sample(n,m,Lambda_ir,Lambda_jr,theta_0)
  }

  Z <- BM_fill(Z, W)
  out <- Z
  return(out)
}

#Resolve other steps
aux_31 <- function(Z,n,m,l,i,j,k,kprime,theta_1,eta_1){

  dprime <- nrow(Z)
  W <- matrix(rnorm(dprime*2^(n+m-l-1)), nrow = dprime)

  Z_i1 <- Z[i,1:2^(n+m-l-1)]; Z_i2 <- Z[i,2:(2^(n+m-l-1)+1)]
  Z_j1 <- Z[j,1:2^(n+m-l-1)]; Z_j2 <- Z[j,2:(2^(n+m-l-1)+1)]
  Lambda_i <- Z_i2 - Z_i1; Lambda_j <- Z_j2 - Z_j1

  for(r in 1:(2^(n+m-l-1))){
    Lambda_ir <- Lambda_i[r];Lambda_jr <- Lambda_j[r]
    W[c(i,j),r] <- five_two_sample(n,m,l,r,k,kprime,Lambda_ir,Lambda_jr,theta_1,eta_1)
  }

  Z <- BM_fill(Z, W)
  out <- Z
  return(out)
}

#Whole of step 3
aux_3 <- function(Z,n,m,i,j,k,kprime,alpha,beta,theta_0,theta_1,eta_1){

  if(m>1){
    for(l in (m-1):1){
      Z <- aux_31(Z,n,m,l,i,j,k,kprime,theta_1,eta_1)
      check <- RB2_check(Z,alpha,beta)
      if(check==1){
        return(0)
      }
    }
  }

  Z <- aux_32(Z,n,m,i,j,k,kprime,theta_0)
  return(Z)
}

#################
#Step 4: find Xi and N
################

#calculate L(kprime)-L(k) given Z
L_diff <- function(Z,i,j,k,kprime){

  Z11 <- Z[i,2*((k+1):kprime)-1]
  Z12 <- Z[i,2*((k+1):kprime)]
  Lambda_i <- Z12-Z11

  Z21 <- Z[j,2*((k+1):kprime)]
  Z22 <- Z[j,2*((k+1):kprime)+1]
  Lambda_j <- Z22 - Z21

  out <- sum(Lambda_i*Lambda_j)
  return(out)
}

psi_fn <- function(Z,n,m,i,j,k,kprime,theta_0,theta_1,eta_1){

  delta <- 2^(-2*(n+m+1))
  scale <- (1 - (theta_0^2)*(delta^2))^(-(kprime-k)/2)

  C <- 1
  if(m > 1){
    for(l in 2:m){
      for(r in 1:(2^(n+m-l))){
        C <- C_fn(n,m,l,r,k,kprime,theta_1,eta_1)*C
      }
    }
  }

  sum1 <- 0
  sum2 <- 0

  for(r in 1:(2^n)){

    Lambda_i <- Z[i,r+1] - Z[i,r]
    Lambda_j <- Z[j,r+1] - Z[j,r]

    tp <- theta_plus(n,m,m,r,k,kprime,theta_1)
    tm <- theta_minus(n,m,m,r,k,kprime,theta_1)
    ep <- theta_plus(n,m,m,r,k,kprime,eta_1)
    em <- theta_minus(n,m,m,r,k,kprime,eta_1)
    rho <- rho_fn(n,m,m,tp,ep)
    g <- g_fn(n,m,m,tp,ep)
    h <- h_fn(n,m,m,ep,rho)
    theta <- theta_fn(tp,tm,em,rho,g,h)
    eta <- eta_fn(tm,ep,em,rho,h)

    sum1 <- sum1 + theta*Lambda_i*Lambda_j
    sum2 <- sum2 + eta*(Lambda_i^2 + Lambda_j^2)
  }

  exp_val <- scale * C * exp(sum1) * exp(sum2)
  out <- log(exp_val)
  return(out)
}

Xi_fn <- function(Z,n,m,i,j,k,kprime,q,theta_0,theta_1,eta_1){

  dprime <- nrow(Z)
  num <- exp(1)*factorial(m-1)*dprime*(dprime-1)
  L_diff <- L_diff(Z,i,j,k,kprime)
  psi <- psi_fn(Z,n,m,i,j,k,kprime,theta_0,theta_1,eta_1)
  denom <- q*exp(theta_0*L_diff - psi)

  out <- num/denom
  return(out)
}

#count the number of RB2s which occur for given i,j
RB2_count_ij <- function(Z,i,j,alpha,beta){

  Z_i <- Z[i,]; Z_j <- Z[j,]
  n <- log2(length(Z_i)-1)

  Z11 <- Z[i,2*(1:2^(n-1))-1]
  Z12 <- Z[i,2*(1:2^(n-1))]
  Lambda_i <- Z12-Z11

  Z21 <- Z[j,2*(1:2^(n-1))]
  Z22 <- Z[j,2*(1:2^(n-1))+1]
  Lambda_j <- Z22 - Z21

  Lambda_prod <- Lambda_i*Lambda_j

  scale <- 2^(-n)
  denoms <- ((1:2^(n-1))^beta)*(scale^(2*alpha))
  num_RB2s <- 0
  for(m in 0:(2^(n-1)-1)){
    s <- abs(cumsum(Lambda_prod[(m+1):2^(n-1)]))
    scaled_Ldiffs <- s/denoms[(m+1):2^(n-1)]
    num_RB2s<- num_RB2s + sum(scaled_Ldiffs > 1)
  }
  return(num_RB2s)
}

#count RB2s for all pairs i,j
RB2_count <- function(Z,alpha,beta){
  dprime <- nrow(Z)
  num <- 0
  for(i in 1:(dprime-1)){
    for(j in (i+1):dprime){
      num <- RB2_count_ij(Z,i,j,alpha,beta)
    }
  }
  return(num)
}

#################
#Procedure Aux, full implementation
#################

aux <- function(Z,alpha,alpha_prime,beta,gamma){

  if(gamma > 1/4){
    return("value of gamma must be <= 1/4")
  }

  if((alpha <= 1/3) | (alpha >= 1/2)){
    return("value of alpha must lie in (1/3,1/2)")
  }

  if((alpha_prime <= alpha) | (alpha_prime >= 1/2)){
    return("value of alpha_prime must lie in (alpha,1/2)")
  }

  if((beta <= 1-alpha) | (beta >= 2*alpha)){
    return("value of beta must lie in (1-alpha,2*alpha")
  }

  dprime <- nrow(Z)
  n <- log2(ncol(Z)-1)

  #step 1: sample m, i, j
  m <- rpois(1,1) + 1

  indices <- sample(1:dprime, 2, replace=FALSE)
  i <- indices[1]; j <- indices[2]

  #step 2: sample k, kprime
  ks <- aux_2(n,m,alpha,alpha_prime,beta,dprime)
  k <- ks$k; kprime <- ks$kprime; q <- ks$q

  #define corresponding parameters for step 3
  theta0_denom <- sqrt(kprime - k)*(2^(-2*alpha_prime*n))*(2^(-m))
  theta_0 <- gamma/theta0_denom

  theta1_denom <- 4 * ( 1 -((theta_0^2) * (2^(-2*(m+n+1)))) )
  theta_1 <- theta_0/theta1_denom
  eta_1 <- (theta_0^2)*(2^(-(n+m)))/(2*theta1_denom)

  #step 3: sample up to resolution n+m
  Z <- aux_3(Z,n,m,i,j,k,kprime,alpha,beta,theta_0,theta_1,eta_1)

  #if an RB2 is detected too early, return F=0
  if(identical(Z,0)){
    return(list(F_bernoulli=0))
  }

  #then check RB2 occurs at proposal; if not, again F=0
  check <- abs(L_diff(Z,i,j,k,kprime))
  crit <- (kprime-k)^beta * 2^(-2*alpha*(n+m))
  if(check < crit){
    return(list(F_bernoulli=0))
  }

  #step 4: calculate Xi, N
  Xi <- Xi_fn(Z,n,m,i,j,k,kprime,q,theta_0,theta_1,eta_1)
  N <- RB2_count(Z,alpha,beta)

  #step 5: sample indicator for new record breaker
  U <- runif(1)
  if(U < Xi/N){
    out <- list(F_bernoulli=1, Z = Z)
    return(out)
  } else {
    return(list(F_bernoulli=0))
  }
}

########################
#3c) Procedure B
########################

################
#Procedure B
################

#step 3 altered for proc B
procB_32 <- function(Z,n,m,i,j,k,kprime,theta_0){

  dprime <- nrow(Z)
  W <- matrix(rnorm(dprime*(2^(n+m-1))), nrow = dprime)

  Z_i1 <- Z[i,1:2^(n+m-1)]; Z_i2 <- Z[i,2:(2^(n+m-1)+1)]
  Z_j1 <- Z[j,1:2^(n+m-1)]; Z_j2 <- Z[j,2:(2^(n+m-1)+1)]
  Lambda_i <- Z_i2 - Z_i1; Lambda_j <- Z_j2 - Z_j1

  for(r in 1:(2^(n+m-1))){
    Lambda_ir <- Lambda_i[r];Lambda_jr <- Lambda_j[r]
    W[c(i,j),r] <- five_one_sample(n,m,Lambda_ir,Lambda_jr,theta_0)
  }

  if(max(abs(W)) > 4*sqrt(m+n)){
    return(list(H_bernoulli=1))
  } else {
    Z <- BM_fill(Z, W)
    return(list(H_bernoulli=0,Z=Z))
  }
}

procB_31 <- function(Z,n,m,l,i,j,k,kprime,theta_1,eta_1){

  dprime <- nrow(Z)
  W <- matrix(rnorm(dprime*2^(n+m-l-1)), nrow = dprime,ncol=2^(n+m-l-1))

  Z_i1 <- Z[i,1:2^(n+m-l-1)]; Z_i2 <- Z[i,2:(2^(n+m-l-1)+1)]
  Z_j1 <- Z[j,1:2^(n+m-l-1)]; Z_j2 <- Z[j,2:(2^(n+m-l-1)+1)]
  Lambda_i <- Z_i2 - Z_i1; Lambda_j <- Z_j2 - Z_j1

  for(r in 1:(2^(n+m-l-1))){
    Lambda_ir <- Lambda_i[r];Lambda_jr <- Lambda_j[r]
    W[c(i,j),r] <- five_two_sample(n,m,l,r,k,kprime,Lambda_ir,Lambda_jr,theta_1,eta_1)
  }

  if(max(abs(W)) > 4*sqrt(m+n-l)){
    return(list(H_bernoulli=1))
  } else {
    Z <- BM_fill(Z, W)
    return(list(H_bernoulli=0,Z=Z))
  }
}

procB_3 <- function(Z,n,m,i,j,k,kprime,alpha,beta,theta_0,theta_1,eta_1){

  #assume no early RB2
  Cbar_nl = 0

  if(m>1){
    for(l in (m-1):1){
      info <- procB_31(Z,n,m,l,i,j,k,kprime,theta_1,eta_1)
      if(info$H_bernoulli == 1){
        return(list(H_bernoulli=1))
      }
      Z <- info$Z
      check <- RB2_check(Z,alpha,beta)
      if(check==1){
        Cbar_nl <- 1
      }
    }
  }

  info <- procB_32(Z,n,m,i,j,k,kprime,theta_0)
  if(info$H_bernoulli == 1){
    return(list(H_bernoulli=1))
  } else {
    Z <- info$Z
    return(list(H_bernoulli=0,Z=Z,Cbar_nl=Cbar_nl))
  }
}
#extra weighting factor for acceptance ratio
#CAREFUL ABOUT W/Z indexing with n - check length(W^n) & length(Z^n)
H_prob <- function(n,m){
  out <- 1
  for(l in 1:m){
    out <- out*pnorm(4*sqrt(n+l+1))^(2^(n+l-1))
  }
  return(out)
}

proc_B <- function(Z,alpha,alpha_prime,beta,gamma,epsilon_0){

  if(gamma > 1/4){
    return("value of gamma must be <= 1/4")
  }

  if((alpha <= 1/3) | (alpha >= 1/2)){
    return("value of alpha must lie in (1/3,1/2)")
  }

  if((alpha_prime <= alpha) | (alpha_prime >= 1/2)){
    return("value of alpha_prime must lie in (alpha,1/2)")
  }

  if((beta <= 1-alpha) | (beta >= 2*alpha)){
    return("value of beta must lie in (1-alpha,2*alpha")
  }

  dprime <- nrow(Z)
  n <- log2(ncol(Z)-1)

  H <- 1
  while(H==1){
    #step 1: sample m, i, j
    m <- rpois(1,1) + 1
    indices <- sample(1:dprime, 2, replace=FALSE)
    i <- indices[1]; j <- indices[2]

    #step 2: sample k, kprime
    ks <- aux_2(n,m,alpha,alpha_prime,beta,gamma)
    k <- ks$k; kprime <- ks$kprime; q <- ks$q

    #define corresponding parameters for step 3
    theta0_denom <- sqrt(kprime - k)*(2^(n))*(2^(-m))
    theta_0 <- gamma/theta0_denom
    theta1_denom <- 4 * ( 1 -((theta_0^2) * (2^(-2*(m+n+1)))) )
    theta_1 <- theta_0/theta1_denom
    eta_1 <- (theta_0^2)*(2^(-(n+m)))/(2*theta1_denom)

    #step 3: sample up to resolution n+m
    info <- procB_3(Z,n,m,i,j,k,kprime,alpha,beta,theta_0,theta_1,eta_1)
    H <- info$H_bernoulli
  }

  #if an RB2 is detected too early, return F=0
  Cbar_nl <- info$Cbar_nl
  if(Cbar_nl==1){
    return(list(F_bernoulli=0))
  }

  #if no ealy RB2, update Z and proceed
  Z <- info$Z
  n <- n+m

  #then check RB2 occurs at proposal; if not, again F=0
  check <- abs(L_diff(Z,i,j,k,kprime))
  crit <- (kprime-k)^beta * 2^(-2*alpha*(n+m))
  if(check < crit){
    return(list(F_bernoulli=0))
  }

  #step 4: calculate Xi, N, P_H
  Xi <- Xi_fn(Z,n,m,i,j,k,kprime,q,theta_0,theta_1,eta_1)
  N <- RB2_count(Z,alpha,beta)
  P_H <- H_prob(n,m)

  #step 5: sample indicator for new record breaker
  U <- runif(1)
  if(U < P_H*Xi/N){
    out <- list(F_bernoulli=1, Z = Z)
    return(out)
  } else {
    return(list(F_bernoulli=0))
  }
}

###############################
#4. Algorithm 2
###############################

############
##Step 1
############

#sample a list of recordbreakers for each dimension, to get N1
alg_211 <- function(dprime){

  pairs_list <- vector("list",dprime)
  N1 <- 0
  #first sample lists of pairs to find N1
  for(i in 1:dprime){
    pairs_i <- record_breakers1()
    pairs_list[[i]] <- pairs_i
    if(nrow(pairs_i)>0){
      N1 <- max(N1, pairs_i[nrow(pairs_i)],1)
    }
  }
  return(list(N1=N1, pairs_list = pairs_list))
}

#conditional on above, sample Zn
alg_212 <- function(dprime, N1, pairs_list){

  Z <- matrix(0,nrow=dprime,ncol=(2^N1+1))
  V_matrix <- matrix(0,nrow=dprime,ncol=N1+1)
  #sample each dimension up to resolution N1
  for(i in 1:dprime){
    pairs_i <- pairs_list[[i]]
    info_i <- ZV_sample(N1,pairs_i)
    Z[i,] <- info_i$Z
    V_matrix[i,] <- info_i$V
  }

  V <- apply(V_matrix,2,max)
  return(list(Z=Z,V=V))
}

#step 1 wrapped into single function
alg_21 <- function(dprime){

  info <-alg_211(dprime)
  N1 <- info$N1
  pairs_list <- info$pairs_list
  out <- alg_212(dprime, N1, pairs_list)
  return(out)
}

################
#Step 2
################

#given i&j, check whether conditions 5.8 & 5.9 are met
five89_check <- function(Z_i,Z_j, alpha_prime, beta, epsilon_0){

  #condition 5.8
  n <- log2(length(Z_i)-1)
  Lambda_i <- diff(Z_i); Lambda_j <- diff(Z_j)
  L_i <- max(abs(Lambda_i)); L_j <- max(abs(Lambda_j))
  crit1 <- 2^(-n*alpha_prime)
  if(max(L_i,L_j) > crit1){
    return(0)
  }

  #condition 5.9
  Lambda_prod <- Lambda_i*Lambda_j
  crit2 <- epsilon_0*((1:2^n)^beta)*(crit1^2)
  for(l in 0:(2^n-1)){
    s <- abs(cumsum(Lambda_prod[(l+1):2^n]))
    crit2_l <- crit2[1:(2^n-l)]
    if(any(s>crit2_l)){
      return(0)
    }
  }
  return(1)
}

five89_check_full <- function(Z,alpha_prime,beta,epsilon_0){
  dprime <- nrow(Z)
  for(i in 1:(dprime-1)){
    Z_i <- Z[i,]
    for(j in (i+1):dprime){
      Z_j <- Z[j,]
      check <- five89_check(Z_i,Z_j,alpha_prime,beta,epsilon_0)
      if(check==0){
        return(0)
      }
    }
  }
  return(1)
}

#resolves multidim Z one further step condtional on no more RB1s
resolve_once <- function(Z){

  dprime <- nrow(Z)
  n <- log2(ncol(Z)-1)
  #CHECK SCALE
  scale <- 2^(-(n+2)/2)

  Z_out <- matrix(0,nrow=dprime,ncol=2^(n+1)+1)

  for(i in 1:dprime){
    W_new <- numeric(2^n)
    Z_midpoints <- apply(matrix(c(Z[i,1:2^n],Z[i,2:(2^n+1)]), ncol=2),1,mean)
    for(j in 1:(2^n)){
      W_new[j] <- normal_body_sample(4*sqrt(n+1))
    }
    Z_new <- Z_midpoints + scale*W_new
    Z_out[i,] <-  c(rbind(Z[i,1:2^n], Z_new), Z[i,2^n+1])
  }
  return(Z_out)
}

#resolve Z until conditions 5.8 & 5.9 are met
alg_22 <- function(Z,alpha_prime,beta,epsilon_0){

  check <- 0
  while(check==0){
    check <- five89_check_full(Z, alpha_prime, beta, epsilon_0)
    if(check==0){
      Z <- resolve_once(Z)
    }
  }
  return(Z)
}

################
#step 3: Procedure B
################

################
#Step 4: no code necessary
################

################
#Step 5: Procedure A
################

################
#Step 6
################

Pn_sample_once <- function(Z){
  n <- log2(ncol(Z)-1)
  dprime <- nrow(Z)
  delta <- 2^(-(n+2))

  W <- matrix(rnorm(dprime*(2^n)),nrow=dprime)
  if(max(abs(W)) > 4*sqrt(n+2)){
    return(list(H_bernoulli=1))
  }

  Z_midpoints <- 1/2 * (Z[,1:2^n] + Z[,2:(2^n+1)])
  Z_new <- Z_midpoints + sqrt(delta)*W
  out <- matrix(0,nrow=dprime,ncol=2^(n+1)+1)
  out[,seq(1,2^(n+1)+1,by=2)] <- Z
  out[,seq(2,2^(n+1),by=2)] <- Z_new
  return(list(H_bernoulli=0,Z=out))
}

Pn_sample <- function(Z,m,alpha,beta){

  for(l in 1:m){
    info <- Pn_sample_once(Z)
    H <- info$H_bernoulli
    if(H==1){
      return(list(P_bernoulli=1))
    }
    Z <- info$Z
    check <- RB2_check(Z,alpha,beta)
    if(check==1){
      return(list(P_bernoulli=1))
    }
  }
  return(list(P_bernoulli=0,Z=Z))
}

alg_26_proposal <- function(Z,m,alpha,beta){

  P_bernoulli <- 1
  while(P_bernoulli==1){
    info <- Pn_sample(Z,m)
    P_bernoulli <- info$P_bernoulli
    Z <- info$Z
  }
  return(Z)
}

alg_26 <- function(Z,m,alpha,alpha_prime,beta,gamma,epsilon_0){

  F_bernoulli <- 1
  while(F_bernoulli==1){
    prop <- alg_26_proposal(Z,m,alpha,beta)
    B <- proc_B(prop,alpha,alpha_prime,beta,gamma,epsilon_0)
    F_bernoulli <- B$F_bernoulli
  }

  out <- B$Z
  return(Z)
}

################
#Algorithm 2, full implementation
################

BM_rp <- function(dprime,d,epsilon,M,alpha=1/3+1/1000,alpha_prime=1/3+2/1000,beta=2/3,gamma=1/4,epsilon_0=0.49){

  if((alpha <= 1/3) | (alpha >= 1/2)){
    return("value of alpha must lie in (1/3,1/2)")
  }
  if((alpha_prime <= alpha) | (alpha_prime >= 1/2)){
    return("value of alpha_prime must lie in (alpha,1/2)")
  }
  if((beta <= 1-alpha) | (beta >= 2*alpha)){
    return("value of beta must lie in (1-alpha,2*alpha)")
  }
  if(gamma > 1/4){
    return("value of gamma must be <= 1/4")
  }
  if((epsilon_0 <= 0) | (epsilon_0 >= 1/2)){
    return("value of epsilon_0 must lie in (0,1/2)")
  }

  #step one: sample BM up to level containing final RB1
  info <- alg_21(dprime)
  Z <- info$Z; V <- info$V; N1 <- log2(ncol(Z)-1) #record V,N1 to get constants for step 5
  print(c("N1",log2(ncol(Z)-1)))

  #step two: ensure conditions 5.8 & 5.9 are met
  Z <- alg_22(Z,alpha_prime,beta,epsilon_0)
  print(c("N1",log2(ncol(Z)-1)))
  #set indicator for existence of another RB2
  more_RB2 <- 1

  #loop steps 3-4
  while(more_RB2 == 1){

    #step three: decide if any RB2 occurs, and if so where
    B <- proc_B(Z,alpha,alpha_prime,beta,gamma,epsilon_0)
    print(c("N2",log2(ncol(Z)-1)))

    if(B$F_bernoulli == 1){
      #return to step 2
      Z <- B$Z
    } else {
      #got to step 5
      more_RB2 <- 0
    }
  }

  N2 <- log2(ncol(Z)-1)
  print(c("N2",N2))

  #step 5: find G and corresponding N3
  K_alpha <- K_alpha(N1, alpha, V)
  #Gamma_L <- Gamma_L(Z,alpha,beta)
  Gamma_L <- 5
  Gamma_R <- Gamma_R(Gamma_L, alpha, beta)
  K_2alpha <- K_2alpha(K_alpha, alpha, Gamma_R)

  G <- find_G(alpha, beta, K_alpha, K_2alpha, Gamma_R, M, d, dprime)
  print(G)
  N3 <- ceiling(-log2(epsilon/G)/(2*alpha - beta))

  return(N3)

  #step 6
  if(N3 > N2){
    m <- N3-n
    Z <- alg_26(Z,m,alpha,alpha_prime,beta,gamma,epsilon_0)
  }

  return(Z)
}

####################################
#5. Euler scheme
####################################

rpeuler_step <- function(n,k,d,dprime,i,mu,sigma,grad_sigma,X,Z_1,Z_2){

  drift <- k*2^(-n)*mu(X)[i]

  vol <- 0
  s <- sigma(X)

  for(j in 1:dprime){
    vol <- s[i,j]*(Z_2[j] - Z_1[j])
  }

  A_diag <- 1/2 * ((Z_1 - Z_2)^2 - 2^(-n))

  rp_correction <- 0

  ds <- grad_sigma(X)
  s1 <- s[i,i]

  for(l in 1:d){
    ds1 <- ds[i,i,l]
    for(m in 1:dprime){
      s2 <- s[l,m]
      ds2 <- ds[l,m,l]
      A <- A_diag[m]
      rp_correction <- rp_correction + A*(s1*ds2 + s2*ds1)
    }
  }

  out <- X[i] + drift + vol + rp_correction
  return(out)
}

rp_euler <- function(n,d,dprime,mu,sigma,grad_sigma,x_0,Z){

  out <- matrix(0, ncol = 2^n + 1, nrow = d)
  X_old <-x_0
  out[,1] <- X_old

  for(k in 1:(2^n)){
    X_new <- numeric(d)
    Z_1 <- Z[,k]; Z_2 <- Z[,k+1]
    for(i in 1:d){
      X_new[i] <- rpeuler_step(n,k,d,dprime,i,mu,sigma,grad_sigma,X_old,Z_1,Z_2)
    }
    out[,k+1] <- X_new
  }

  return(out)
}

