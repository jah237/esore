####################################
#Epsilon-strong for Brownian motion
####################################

########################
#0. Bernoulli samplers (general purpose)
########################

##given a function f to calculate an alternating series f_n in [0,1], carries
##out a retrospective Bernoulli sample

retro_bernoulli <- function(f, ...){

  U <- runif(1)
  n <- 1

  while(TRUE){

    upper_bound <- f(..., n = n)
    lower_bound <- f(..., n = n+1)

    if ((0 <= upper_bound) & (upper_bound < U)){

      return(0)
    } else if (lower_bound > U){

      return(1)
    } else {

      n <- n+2
    }
  }
}

#for non-binary samples
retro_categorical <- function(f, param_range, ...){

  U <- runif(1)

  n <- 1

  while(TRUE){

    old_lb <- 0
    old_ub <- 0

    for(p in seq_along(param_range)){

      i <- param_range[p]

      lb <- max(0, old_lb + f(..., i = i, n = n+1))
      ub <- min(old_ub + f(..., i = i, n = n), 1)

      #technical fix
      if(ub < 0){

        n <- n+2
        break

      } else if((old_ub < U) & (U < lb)){

        #current category determined
        return(p)

      } else if((lb < U) & (U < ub)){

        #undetermined - increase precision from 1st category
        n <- n+2
        break

      } else {

        #check following categories
        old_lb <- lb
        old_ub <- ub
      }
    }
  }
}

########################
#1. Primitive functions
########################

############
#1a) numerical versions, for validation
############

sigma_bar_j <- function(x, y, l, delta, psi, j){

  out <- exp(-2/l * (delta*j + psi - x) * (delta*j + psi - y))
  return(out)
}

tau_bar_j <- function(x, y, l, delta, j){

  out <- exp(-2*j/l * (delta^2 * j + delta*(x-y)))
  return(out)
}

sigma_j <- function(L, U, l, x, y, j){

  out <- sigma_bar_j(x, y, l, (U-L), L, j) + sigma_bar_j(-x, -y, l, (U-L), -U, j)
  return(out)
}

tau_j <- function(L, U, l, x, y, j){

  out <- tau_bar_j(x, y, l, (U-L), j) + tau_bar_j(-x, -y, l, (U-L), j)
  return(out)
}

sigma_sum <- function(L, U, l, x, y, n){

    out <- sum(sigma_j(L, U, l, x, y, 1:n))
    return(out)
}

tau_sum <- function(L, U, l, x, y, n){

  out <- sum(tau_j(L, U, l, x, y, 1:n))
  return(out)
}

############
#1b) coefficient versions (exact), for final use
############

#store the five coeffs for the expression c_i * exp(a_i + b_i*w) * pi(w) * Ind_{L_i < w < U_i}
sigma_bar_j2 <- function(z, l, delta, psi, j){

  coeffs <- -2/l * (delta*j + psi - z) * c(delta*j + psi, -1)
  out <- c(coeffs, 1)
  return(out)
}

tau_bar_jx <- function(z, l, delta, j){

  coeffs <- -2*j*delta/l * c((delta*j - z), 1)
  out <- c(coeffs, -1)
  return(out)
}

tau_bar_jy <- function(z, l, delta, j){

  coeffs <- -2*j*delta/l * c((delta*j + z), -1)
  out <- c(coeffs, -1)
  return(out)
}

#produce two rows of coeffs
sigma_j2 <- function(L, U, l, z, j){

  out <- rbind(sigma_bar_j2(z, l, (U-L), L, j), c(1, -1, 1)*sigma_bar_j2(-z, l, (U-L), -U, j))
  out <- cbind(out, cbind(rep(L, 2), rep(U, 2)))
  return(out)
}

tau_j2 <- function(L, U, l, z, w_is, j){

  if (w_is == "x"){

    out <- rbind(tau_bar_jx(z, l, (U-L), j), c(1, -1, 1)*tau_bar_jx(-z, l, (U-L), j))

  } else if (w_is == "y"){

    out <- rbind(tau_bar_jy(z, l, (U-L), j), c(1, -1, 1)*tau_bar_jy(-z, l, (U-L), j))
  }

  out <- cbind(out, cbind(rep(L, 2), rep(U, 2)))
  return(out)
}

#for convenience in Sn_zeta
sigma_sum2 <- function(L, U, l, z, n){

  out <- matrix(0, nrow = 2*n, ncol = 5)

  for (j in 1:n){

    out[(2*j-1):(2*j), ] <- sigma_j2(L, U, l, z, j)
  }

  return(out)
}

tau_sum2 <- function(L, U, l, z, w_is, n){

  out <- matrix(0, nrow = 2*n, ncol = 5)

  for (j in 1:n){

    out[(2*j-1):(2*j), ] <- tau_j2(L, U, l, z, w_is, j)
  }

  return(out)
}

#get numerical values from coeffs, for validation
coeffs_eval <- function(coeffs, w){

  L <- coeffs[4]
  U <- coeffs[5]

  if ((L < w) & (w < U)){

    out <- coeffs[3]*exp(coeffs[1] + coeffs[2]*w)
    return(out)
  } else {

    return(0)
  }
}

get_value <- function(array, w){

  out <- sum(apply(array, 1, coeffs_eval, w=w))
  return(out)
}

#################
#2. Alternating series
#################

########
#2a) Numerical versions, for validation
########

#zeta
Sn_zeta <- function(L, U, l, x, y, n){

  if ((L < x) & (x < U) & (L < y) & (y < U)){

    if (n == 1){

      out <- sigma_j(L, U, l, x, y, 1)
      return(out)

    } else {

      n_sigma <- floor((n+1)/2)
      n_tau <- floor(n/2)

      nn_sigma <<- n_sigma
      nn_tau <<- n_tau

      out <- sigma_sum(L, U, l, x, y, n_sigma) - tau_sum(L, U, l, x, y, n_tau)
      return(out)

    }
  } else {

    return(1)
  }
}

#beta
Sn_beta <- function(L_down, L_up, U_down, U_up, l, x, y, n){

  Ld <<- L_down; Lu <<- L_up; Ud <<- U_down; Uu <<- U_up
  ll <<- l; xx <<- x; yy <<- y; nn <<- n

  out <- (Sn_zeta(L_down, U_down, l, x, y, n) - Sn_zeta(L_down, U_up, l, x, y, n+1)
          + Sn_zeta(L_up, U_up, l, x, y, n) - Sn_zeta(L_up, U_down, l, x, y, n+1))

  return(out)
}

#rho
rho_term <- function(L, U, q, r, x, w, y, n, k){

  if (k == 1){

    ns <- c(n+1,n)
  } else {

    ns <- c(n, n+1)
  }

  out <- (1 - Sn_zeta(L, U, q, x, w, ns[1]) - Sn_zeta(L, U, r, w, y, ns[1])
          +  Sn_zeta(L, U, q, x, w, ns[2])*Sn_zeta(L, U, r, w, y, ns[2]))
  return(out)
}

Sn_rho <- function(L_down, L_up, U_down, U_up, q, r, x, w, y, n){

  S1 <- rho_term(L_down, U_up, q, r, x, w, y, n, 1)
  S2 <- rho_term(L_up, U_up, q, r, x, w, y, n, 2)
  S3 <- rho_term(L_down, U_down, q, r, x, w, y, n, 2)
  S4 <- rho_term(L_up, U_down, q, r, x, w, y, n, 1)

  out <- S1 - S2 - S3 + S4
  return(out)
}

##########
#2b) coeffiecient versions (exact), for final use
##########
Sn_zeta2 <- function(L, U, l, z, w_is, n){

  if (n==1){

    out <- sigma_j2(L, U, l, z, 1)

  } else {

    n_sigma <- floor((n+1)/2)
    n_tau <- floor(n/2)

    nrow <- 2*(n_sigma + n_tau)
    out <- matrix(0, nrow = nrow, ncol = 5)

    out[1:(2*n_sigma), ] <- sigma_sum2(L, U, l, z, n_sigma)
    out[1:(2*n_tau) + (2*n_sigma), ] <- tau_sum2(L, U, l, z, w_is, n_tau)

  }

  return(out)
}

Sn_zeta_i2 <- function(L_down, L_up, U_down, U_up, q, r, x, y, n, i){

  if (i%%2 == 1){

    l <- q
    z <- x
    w_is <- "y"
  } else {

    l <- r
    z <- y
    w_is <- "x"
  }

  if (i <= 2){
    L <- L_down; U <- U_up
  } else if (i <= 4){
    L <- L_up; U <- U_up
  } else if (i <= 6){
    L <- L_down; U <- U_down
  } else {
    L <- L_up; U <- U_down
  }

  out <- Sn_zeta2(L, U, l, z, w_is, n)
  return(out)
}


#changes sign of ci
sign_flip2 <- function(A){

  A[, 3] <- -A[, 3]

  return(A)
}

#finds coefficient matrix for a product of two expressions
coeff_product2 <- function(A, B){

  out <- matrix(0, nrow = nrow(A)*nrow(B), ncol = 5)

  out[, 1] <- apply(expand.grid(A[, 1], B[, 1]), 1, sum)
  out[, 2] <- apply(expand.grid(A[, 2], B[, 2]), 1, sum)
  out[, 3] <- apply(expand.grid(A[, 3], B[, 3]), 1, prod)
  out[, 4] <- apply(expand.grid(A[, 4], B[, 4]), 1, max)
  out[, 5] <- apply(expand.grid(A[ ,5], B[, 5]), 1, min)

  return(out)
}

rho_term2 <- function(L, U, q, r, x, y, n, k){

  if (k == 1){

    ns <- c(n+1,n)

    out <- matrix(0, nrow = 4*n^2 + 4*n + 4, ncol = 5)
  } else {

    ns <- c(n, n+1)

    out <- matrix(0, nrow = 4*n^2 + 12*n + 4, ncol = 5)
  }

  out[1:(2*ns[1]), ] <- sign_flip2(Sn_zeta2(L, U, q, x, "y", ns[1]))
  out[1:(2*ns[1]) + (2*ns[1]), ] <- sign_flip2(Sn_zeta2(L, U, r, y, "x", ns[1]))
  out[(4*ns[1] + 1):nrow(out), ] <- coeff_product2(Sn_zeta2(L, U, q, x, "y", ns[2]), Sn_zeta2(L, U, r, y, "x", ns[2]))

  return(out)
}

Sn_rho2 <- function(L_down, L_up, U_down, U_up, q, r, x, y, n){

  out <- matrix(0, nrow = 16*(n+1)^2, ncol = 5)

  S1 <- rho_term2(L_down, U_up, q, r, x, y, n, 1)
  n1 <- nrow(S1)
  S2 <- rho_term2(L_up, U_up, q, r, x, y, n, 2)
  n2 <- nrow(S2)
  S3 <- rho_term2(L_down, U_down, q, r, x, y, n, 2)
  S4 <- rho_term2(L_up, U_down, q, r, x, y, n, 1)

  out[1:n1, ] <- S1
  out[1:n2 + n1, ] <- sign_flip2(S2)
  out[1:n2 + n1 + n2, ] <- sign_flip2(S3)
  out[1:n1 + n1 + 2*n2, ] <- S4

  return(out)
}

######################
#3. Sampling procedures
######################

###########
#3a) Initial layers
###########

sample_initial_layers <- function(as, l, x, y){

  U <- runif(1)

  x_bar <- min(x, y)
  y_bar <- max(x, y)

  n <- 1

  while(TRUE){

    s1 <- seq_along(as)

    pairs <- expand.grid(s1, s1)
    param_range <- seq_len(nrow(pairs))

    as <- c(0, as)

    count <- 0

    while(TRUE){

      old_ub <- 0
      old_lb <- 0

      for(i in param_range){

        a1 <- as[pairs[i, 1]]; a2 <- as[pairs[i, 1] + 1]
        b1 <- as[pairs[i, 2]]; b2 <- as[pairs[i, 2] + 1]

        ub <- min(old_ub + Sn_beta(x_bar - a2, x_bar - a1, y_bar + b1, y_bar + b2, l, x, y, n), 1)
        lb <- old_lb + Sn_beta(x_bar - a2, x_bar - a1, y_bar + b1, y_bar + b2, l, x, y, n+1)

        if ((U < lb) & (old_ub < U)){
          return(c(x_bar - a2, x_bar - a1, y_bar + b1, y_bar + b2))

        } else if ((lb < U) & (U < ub)) {

          n <- n+2
          break

        } else {

          old_lb <- lb
          old_ub <- ub
        }
      }

      if(i == nrow(pairs)){

        as <- c(as[2:length(as)], 2^(1:5) * as[length(as)])
        break
      }
    }
  }
}

sil_test <- function(as, l, x, y){

  U <- runif(1)

  x_bar <- min(x, y)
  y_bar <- max(x, y)

  n <- 1

  while(TRUE){

    s1 <- seq_along(as)

    pairs <- expand.grid(s1, s1)
    param_range <- seq_len(nrow(pairs))

    as <- c(0, as)

    while(TRUE){

      old_ub <- 0
      old_lb <- 0

      for(i in param_range){

        a1 <- as[pairs[i, 1]]; a2 <- as[pairs[i, 1] + 1]
        b1 <- as[pairs[i, 2]]; b2 <- as[pairs[i, 2] + 1]

        ub <- min(old_ub + Sn_beta(x_bar - a2, x_bar - a1, y_bar + b1, y_bar + b2, l, x, y, n), 1)
        lb <- old_lb + Sn_beta(x_bar - a2, x_bar - a1, y_bar + b1, y_bar + b2, l, x, y, n+1)

        if ((U < lb) & (old_ub < U)){

          return(i)

        } else if ((lb < U) & (U < ub)) {

          n <- n+2
          break

        } else {

          old_lb <- lb
          old_ub <- ub
        }
      }

      if(i == nrow(pairs)){

        as <- c(as[2:length(as)], 2^(1:5) * as[length(as)])
        break
      }
    }
  }
}

###########
#3b) Midpoint
###########

###
#3bi) Numerical versions, for validation
###

pi_fun <- function(q, r, x, w, y){

  l <- q + r

  mean <- r/l * x + q/l * y
  var <- q*r/l

  out <- exp(-(1/(2*var))*(w-mean)^2)
  return(out)
}

fn <- function(L_down, L_up, U_down, U_up, q, r, x, w, y, n){

  out <- Sn_rho(L_down, L_up, U_down, U_up, q, r, x, w, y, n)*pi_fun(q, r, x, w, y)
  return(out)
}

fn_vec <- function(L_down, L_up, U_down, U_up, q, r, x, w, y, n){

  out <- sapply(w, fn, L_down = L_down, L_up = L_up, U_down = U_down, U_up = U_up,
                q = q, r = r, x = x, y = y, n = n)
  return(out)
}

Fn_unnormalised <- function(L_down, L_up, U_down, U_up, q, r, x, w, y, n){

  if(w <= L_down){

    return(0)
  } else{

    out <- integrate(fn_vec, lower = L_down, upper = min(w, U_up), L_down = L_down, L_up = L_up, U_down = U_down, U_up = U_up,
                     q = q, r = r, x = x, y = y, n = n)$value

    return(out)
  }
}

Fn <- function(L_down, L_up, U_down, U_up, q, r, x, w, y, n){

  if (w <= L_down){

    return(0)
  } else if (w >= U_up){

    return(1)
  } else{

    denom <- Fn_unnormalised(L_down, L_up, U_down, U_up, q, r, x, U_up, y, n)

    out <- Fn_unnormalised(L_down, L_up, U_down, U_up, q, r, x, w, y, n)/denom
    return(out)
  }
}

propose_midpoint <- function(L_down, L_up, U_down, U_up, q, r, x, y, n, tol){

  U <- runif(1)

  l <- L_down
  u <- U_up

  error <- 2*tol

  while (error > tol){

    candidate <- (l+u)/2
    cum_prob <- Fn(L_down, L_up, U_down, U_up, q, r, x, l, y, n)

    error <- abs(cum_prob - U)

    if (cum_prob < U){

      l <- candidate
    } else if (cum_prob > U){

      u <- candidate
    }
  }

  return(candidate)
}

###
#3bii) Coefficient versions, for final use
###

pi_params <- function(q, r, x, y){

  l <- q + r

  mu <- x * r/l + y * q/l
  sigma2 <- q*r/l

  return(list(mu = mu, sigma2 = sigma2))
}

fn2 <- function(L_down, L_up, U_down, U_up, q, r, x, y, n){

  old_coeffs <- Sn_rho2(L_down, L_up, U_down, U_up, q, r, x, y, n)
  ai <- old_coeffs[, 1]
  bi <- old_coeffs[, 2]
  ci <- old_coeffs[, 3]

  params <- pi_params(q, r, x, y)
  mu <- params$mu
  sigma2 <- params$sigma2

  new_coeffs <- matrix(0, nrow = nrow(old_coeffs), ncol = 4)

  weights <- ci * exp(ai + bi*mu + 1/2 * sigma2*bi^2)
  means <- mu + sigma2*bi

  new_coeffs[, 1] <- weights
  new_coeffs[, 2] <- means
  new_coeffs[, 3:4] <- old_coeffs[, 4:5]

  return(new_coeffs)
}

eval_fn2 <- function(L_down, L_up, U_down, U_up, q, r, x, w, y, n){

  coeffs <- fn2(L_down, L_up, U_down, U_up, q, r, x, y, n)
  sigma2 <- q*r/(q+r)

  out_fn <- function(coeff_row, sigma2, w){

    L <- coeff_row[3]
    U <- coeff_row[4]

    if ((L <= w) & (w <= U)){

      out <- coeff_row[1] * exp(-(1/(2*sigma2)) * (w - coeff_row[2])^2)
      return(out)
    } else {

      return(0)
    }
  }

  out <- sum(apply(coeffs, 1, out_fn, sigma2=sigma2, w=w))
  return(out)
}

#get numerical values from coeffs
cdf_eval <- function(coeff_row, sigma, w){

  L <- coeff_row[3]
  U <- coeff_row[4]

  if (w <= L){

    return(0)
  } else {

    up_val <- pnorm(min(w, coeff_row[4]), mean = coeff_row[2], sd = sigma)
    lw_val <- pnorm(coeff_row[3], mean = coeff_row[2], sd = sigma)

    out <- coeff_row[1]*(up_val - lw_val)
    return(out)
  }
}

Fn2_unnormalised <- function(L_down, L_up, U_down, U_up, q, r, x, w, y, n){

  if ((x < L_down) | (x > U_up) | (y < L_down) | (y > U_up)){

    return("x/y out of bounds")
  } else {

    coeffs <- fn2(L_down, L_up, U_down, U_up, q, r, x, y, n)

    sigma <- sqrt(q*r/(q+r))

    Ls <- coeffs[, 3]
    Us <- coeffs[, 4]

    out <- sqrt(2*pi) * sigma * sum(apply(coeffs, 1, cdf_eval, w=w, sigma = sigma))
  }

  return(out)
}

Fn2 <- function(L_down, L_up, U_down, U_up, q, r, x, w, y, n){

  if (w <= L_down){

    return(0)
  } else if (w >= U_up){

    return(1)
  } else{

    denom <- Fn2_unnormalised(L_down, L_up, U_down, U_up, q, r, x, U_up, y, n)
    out <- Fn2_unnormalised(L_down, L_up, U_down, U_up, q, r, x, w, y, n)/denom

    return(out)
  }
}

propose_midpoint2 <- function(L_down, L_up, U_down, U_up, q, r, x, y, tol){

  U <- runif(1)

  l <- L_down
  u <- U_up

  range <- (u - l)

  error <- 1

  while (error > tol){

    candidate <- (l+u)/2
    cum_prob <- Fn2(L_down, L_up, U_down, U_up, q, r, x, candidate, y, 1)

    if (cum_prob < U){

      l <- candidate
    } else if (cum_prob > U){

      u <- candidate
    }

    error <- 0.5 * error
  }

  return(candidate)
}

acc_ratio <- function(L_down, L_up, U_down, U_up, q, r, x, w, y, n){

  denom <- Sn_rho(L_down, L_up, U_down, U_up, q, r, x, w, y, 1)
  out <- Sn_rho(L_down, L_up, U_down, U_up, q, r, x, w, y, n)/denom

  return(out)
}

sample_midpoint <- function(L_down, L_up, U_down, U_up, q, r, x, y){

  tol <- 0.0001

  while(TRUE){

    w <- propose_midpoint2(L_down, L_up, U_down, U_up, q, r, x, y, tol)

    decision <- retro_bernoulli(acc_ratio, L_down = L_down, L_up = L_up, U_down = U_down, U_up = U_up,
                                q = q, r = r, x = x, w = w, y = y)

    if(decision == 1){

      return(w)
    }
  }
}

sample_midpoint_tst <- function(L_down, L_up, U_down, U_up, q, r, x, y){

  tol <- 0.0001

  w <- propose_midpoint2(L_down, L_up, U_down, U_up, q, r, x, y, tol)


  decision <- retro_bernoulli(acc_ratio, L_down = L_down, L_up = L_up, U_down = U_down, U_up = U_up,
                              q = q, r = r, x = x, w = w, y = y)

  print(c(w, decision))
}

############
#3c) New layers
############

#define alternating series for layer sampling
Sn_beta_i <- function(L_down, L_up, U_down, U_up, q, r, x, w, y, n, i){

  w_x <- min(x, w); w_y <- min(w, y)
  wx <- max(x, w); wy <- max(w, y)

  denom <- Sn_rho(L_down, L_up, U_down, U_up, q, r, x, w, y, n+1)

  if(i == 1){
    out <- Sn_beta(L_down, L_up, U_down, U_up, q, x, w, n)*Sn_beta(L_down, L_up, U_down, U_up, r, w, y, n)
  } else if (i ==2){
    out <- Sn_beta(L_down, L_up, U_down, U_up, q, x, w, n)*Sn_beta(L_up, w_y, U_down, U_up, r, w, y, n)
  } else if (i == 3){
    out <- Sn_beta(L_down, L_up, U_down, U_up, q, x, w, n)*Sn_beta(L_down, L_up, wy, U_down, r, w, y, n)
  } else if (i == 4){
    out <- Sn_beta(L_down, L_up, U_down, U_up, q, x, w, n)*Sn_beta(L_up, w_y, wy, U_down, r, w, y, n)
  } else if (i == 5){
    out <- Sn_beta(L_up, w_x, U_down, U_up, q, x, w, n)*Sn_beta(L_down, L_up, U_down, U_up, r, w, y, n)
  } else if (i == 6){
    out <- Sn_beta(L_up, w_x, U_down, U_up, q, x, w, n)*Sn_beta(L_down, L_up, wy, U_down, r, w, y, n)
  } else if (i == 7){
    out <- Sn_beta(L_down, L_up, wx, U_down, q, x, w, n)*Sn_beta(L_down, L_up, U_down, U_up, r, w, y, n)
  } else if (i == 8){
    out <- Sn_beta(L_down, L_up, wx, U_down, q, x, w, n)*Sn_beta(L_up, w_y, U_down, U_up, r, w, y, n)
  } else if (i == 9){
    out <- Sn_beta(L_up, w_x, wx, U_down, q, x, w, n)*Sn_beta(L_down, L_up, U_down, U_up, r, w, y, n)
  }

  out <- out/denom
  return(out)
}

update_layers <- function(L_down, L_up, U_down, U_up, q, r, x, w, y){

  cat <- retro_categorical(Sn_beta_i, param_range = seq_len(9), L_down = L_down, L_up = L_up,
                           U_down = U_down, U_up = U_up, q = q, r = r, x = x, w = w, y = y)

  w_x <- min(x, w); w_y <- min(w, y)
  wx <- max(x, w); wy <- max(w, y)

  if (cat == 1){
    out <- c(L_down, L_up, U_down, U_up, L_down, L_up, U_down, U_up)
  } else if (cat == 2){
    out <- c(L_down, L_up, U_down, U_up, L_up, w_y, U_down, U_up)
  } else if (cat == 3){
    out <- c(L_down, L_up, U_down, U_up, L_down, L_up, wy, U_down)
  } else if (cat == 4){
    out <- c(L_down, L_up, U_down, U_up, L_up, w_y, wy, U_down)
  } else if (cat == 5){
    out <- c(L_up, w_x, U_down, U_up, L_down, L_up, U_down, U_up)
  } else if (cat == 6){
    out <- c(L_up, w_x, U_down, U_up, L_down, L_up, wy, U_down)
  } else if (cat == 7){
    out <- c(L_down, L_up, wx, U_down, L_down, L_up, U_down, U_up)
  } else if (cat == 8){
    out <- c(L_down, L_up, wx, U_down, L_up, w_y, U_down, U_up)
  } else if (cat == 9){
    out <- c(L_up, w_x, wx, U_down, L_down, L_up, U_down, U_up)
  }

  return(out)
}

update_layers_test <- function(L_down, L_up, U_down, U_up, q, r, x, w, y){

  cat <- retro_categorical(Sn_beta_i, param_range = seq_len(9), L_down = L_down, L_up = L_up,
                           U_down = U_down, U_up = U_up, q = q, r = r, x = x, w = w, y = y)

  return(cat)
}

#############
#3d) Refine layers
#############

upper_ratio <- function(L_down, L_up, U_down, U_up, U_star=NA, l, x, y, n){

  if(is.na(U_star) == TRUE){

    U_star <- 1/2 * (U_down + U_up)
  }

  denom <- Sn_beta(L_down, L_up, U_down,  U_up, l, x, y, n+1)
  out <- Sn_beta(L_down, L_up, U_star, U_up, l, x, y, n)/denom

  return(out)
}

lower_ratio <- function(L_down, L_up, U_down, U_up, L_star=NA, l, x, y, n){

  if(is.na(L_star) == TRUE){

    L_star <- 1/2 * (L_down + L_up)
  }

  denom <- Sn_beta(L_down, L_up, U_down,  U_up, l, x, y, n+1)
  out <- Sn_beta(L_down, L_star, U_down,  U_up, l, x, y, n)/denom

  return(out)
}

refine_upper_layer <- function(info, U_star=NA){

  L_down <- info[5]; L_up <- info[6]
  U_down <- info[7]; U_up <- info[8]

  l <- info[4]-info[3]
  x <- info[1]; y <- info[2]

  decision <- retro_bernoulli(upper_ratio2, L_down = L_down, L_up = L_up,
                              U_down = U_down,  U_up = U_up, U_star = U_star,
                              l = l, x = x, y = y)

  if(decision == 1){

    info[7] <- U_star
  } else {

    info[8] <- U_star
  }

  return(info)
}

refine_lower_layer <- function(info, L_star=NA){

  L_down <- info[5]; L_up <- info[6]
  U_down <- info[7]; U_up <- info[8]

  l <- info[4]-info[3]
  x <- info[1]; y <- info[2]

  decision <- retro_bernoulli(lower_ratio2, L_down = L_down, L_up = L_up,
                              U_down = U_down,  U_up = U_up, L_star = L_star,
                              l = l, x = x, y = y)

  if(decision == 1){

    info[6] <- L_star
  } else {

    info[5] <- L_star
  }

  return(info)
}

refine_layers_generic <- function(info){

  l <- info[4]-info[3]
  thresh <- sqrt(l/2)

  while(info[6] - info[5] > thresh){
    info <- refine_lower_layer(info)
  }

  while(info[8] - info[7] > thresh){
    info <- refine_upper_layer(info)
  }

  return(info)
}

#refine the pair of layers resulting from a new midpoint (convenience fn)
refine_layer_pair <- function(info_pair){

  out <- t(apply(info_pair,1,refine_layers_generic))
  return(out)
}


##################################
#5. Epsilon-strong algorithms (fixed epsilon)
##################################

####################
#5a) One-dimensional
####################

esbm_1 <- function(x,t, epsilon){

  y <- rnorm(1,x,t)

  as <- c(abs(x-y)/2,abs(x-y),2*abs(x-y))
  init_layers <- sample_initial_layers(as, t, x, y)
  L_down <- init_layers[1];
  U_up <- init_layers[4]

  info_list <- c(x,y,0,t,init_layers)
  max_err <- U_up - L_down
  err_list <- as.vector(max_err)

  while(max_err > epsilon){

    length <- nrow(info_list)
    new_length <- length + sum(err_list > epsilon)
    count <- 1

    new_info_list <- matrix(0,nrow=new_length,ncol=8)
    new_err_list <- numeric(new_length)

    for(i in nrow(info_list)){

      if(err_list[i] > epsilon){

        info <- info_list[i,]
        #WRITE A SPLIT AND REFINE FUNCTION
        new_info_list[c(count,count+1),] <- refine_layer_pair(transect(info))

        err_1 <- new_info_list[count,8] - new_info_list[count,5]
        err_2 <- new_info_list[count+1,8] - new_info_list[count+1,5]
        new_err_list[c(count,count+1)] <- c(err_1,err_2)
        count <- count+2

      } else {

        new_info_list[count,] <- info
        new_err_list[count] <- err_list[i]
        count <- count+1
      }
    }

    info_list <- new_info_list
    err_list <- new_err_list
    max_err <- max(err_list)
  }

  return(info_list)
}

earl2 <- function(info){
  out <- earl(info[1],info[2],info[3],info[4],info[5],info[6],info[7],info[8])$ref
  return(out)
}

esbm_12 <- function(x, t, epsilon){

  y <- rnorm(1,x,t)
  info_list <- matrix(eadl(0,t,x,y)$layer, nrow=1)

  L_down <- info_list[1,5];
  U_up <- info_list[1,8]
  max_err <- U_up - L_down
  err_list <- as.vector(max_err)

  while(max_err > epsilon){

    length <- nrow(info_list)
    new_length <- length + sum(err_list > epsilon)
    count <- 1

    new_info_list <- matrix(0,nrow=new_length,ncol=8)
    new_err_list <- numeric(new_length)

    for(i in 1:nrow(info_list)){

      info <- info_list[i,]

      if(err_list[i] > epsilon){

        midpoint <- earp(info[1],0.5*(info[1]+info[2]),info[2],info[3],info[4],info[5],info[6],info[7],info[8])$draw
        split_info <- eabl(info[1],0.5*(info[1]+info[2]),info[2],info[3],midpoint,info[4],info[5],info[6],info[7],info[8])$lyr2
        new_info_list[c(count,count+1),] <- t(apply(split_info,1,earl2))


        err_1 <- new_info_list[count,8] - new_info_list[count,5]
        err_2 <- new_info_list[count+1,8] - new_info_list[count+1,5]
        new_err_list[c(count,count+1)] <- c(err_1,err_2)
        count <- count+2

      } else {

        new_info_list[count,] <- info
        new_err_list[count] <- err_list[i]
        count <- count+1
      }
    }

    info_list <- new_info_list
    err_list <- new_err_list
    max_err <- max(err_list)
  }

  return(info_list)
}

################################
#5b) Multiple dimensions
################################

esbm <- function(x, s, t, epsilon){

  d <- length(x)
  y <- rnorm(d,x,t-s)

  info_list <- matrix(0,nrow=d,ncol=8)
  max_err <- 0
  for(j in 1:d){
    info_list[j,] <- eadl(s,t,x[j],y[j])$layer
    max_err <- max(max_err, info_list[j,8] - info_list[j,5])
  }

  err_list <- as.vector(max_err)

  while(max_err > 2*epsilon){


    length <- nrow(info_list)/d
    new_length <- length + sum(err_list > epsilon)
    count <- 1

    new_info_list <- matrix(0,nrow=d*new_length,ncol=8)
    new_err_list <- numeric(new_length)

    for(i in 1:length){
      info <- info_list[((i-1)*d+1):(i*d),]

      if(err_list[i] > epsilon){

        max_err_1 <- max_err_2 <- 0
        for(j in 1:d){

          midpoint <- earp(0,0.5*(info[j,2]-info[j,1]),info[j,2]-info[j,1],info[j,3],info[j,4],info[j,5],info[j,6],info[j,7],info[j,8])$draw
          split_info <- eabl(0,0.5*(info[j,2]-info[j,1]),info[j,2]-info[j,1],info[j,3],midpoint,info[j,4],info[j,5],info[j,6],info[j,7],info[j,8])$lyr2
          new_info_list[c((count-1)*d+j,count*d+j),] <- t(apply(split_info,1,earl2))
          new_info_list[c((count-1)*d+j,count*d+j),1:2] <- new_info_list[c((count-1)*d+j,count*d+j),1:2] + info[j,1]

          err_1 <- new_info_list[(count-1)*d+j,8] - new_info_list[(count-1)*d+j,5]
          err_2 <- new_info_list[count*d+j,8] - new_info_list[count*d+j,5]
          max_err_1 <- max(max_err_1, err_1)
          max_err_2 <- max(max_err_2, err_2)
        }

        new_err_list[c(count,count+1)] <- c(max_err_1,max_err_2)
        count <- count+2

      } else {

        new_info_list[((count-1)*d+1):(count*d),] <- info
        new_err_list[count] <- err_list[i]
        count <- count+1
      }
    }

    info_list <- new_info_list
    err_list <- new_err_list
    max_err <- max(err_list)
  }

  ss <- info_list[d*1:(nrow(info_list)/d),1]
  ts <- info_list[d*1:(nrow(info_list)/d),2]
  xs <- matrix(info_list[,3],byrow=TRUE,ncol=2)
  ys <- matrix(info_list[,4],byrow=TRUE,ncol=2)
  centres <- matrix(0.5*(info_list[,8]+info_list[,5]),byrow=TRUE,ncol=d)
  errors <- err_list/2
  extra_info <- matrix(t(info_list[,3:8]),byrow=TRUE,ncol=6*d)
  fin_loc <- info_list[(nrow(info_list)-d+1):nrow(info_list),4]

  out <- list(ss=ss,ts=ts,xs=xs,ys=ys,centres=centres,errors=errors,
              extra_info=extra_info,fin_loc=fin_loc)
  return(out)
}

esbm_cond <- function(ss, ts, d, extra_info, epsilon){

  info_list <- matrix(0,nrow=d*length(as.vector(ss)),ncol=8)
  info_list[,1] <- rep(ss,rep(d,length(as.vector(ss))))
  info_list[,2] <- rep(ts,rep(d,length(as.vector(ss))))
  info_list[,3:8] <- matrix(extra_info,byrow=TRUE,ncol=6)

  err_list <- apply(matrix(info_list[,8]-info_list[,5],byrow=TRUE,ncol=d),1,max)
  max_err <- max(err_list)

  while(max_err > 2*epsilon){

    length <- nrow(info_list)/d
    new_length <- length + sum(err_list > epsilon)
    count <- 1

    new_info_list <- matrix(0,nrow=d*new_length,ncol=8)
    new_err_list <- numeric(new_length)

    for(i in 1:length){

      info <- info_list[((i-1)*d+1):(i*d),]
      info_check <<- info[,1:2]
      if(err_list[i] > epsilon){

        max_err_1 <- max_err_2 <- 0
        for(j in 1:d){

          midpoint <- earp(0,0.5*(info[j,2]-info[j,1]),(info[j,2]-info[j,1]),info[j,3],info[j,4],info[j,5],info[j,6],info[j,7],info[j,8])$draw
          split_info <- eabl(0,0.5*(info[j,2]-info[j,1]),(info[j,2]-info[j,1]),info[j,3],midpoint,info[j,4],info[j,5],info[j,6],info[j,7],info[j,8])$lyr2
          new_info_list[c((count-1)*d+j,count*d+j),] <- t(apply(split_info,1,earl2))
          new_info_list[c((count-1)*d+j,count*d+j),1:2] <- new_info_list[c((count-1)*d+j,count*d+j),1:2] + info[j,1]

          err_1 <- new_info_list[(count-1)*d+j,8] - new_info_list[(count-1)*d+j,5]
          err_2 <- new_info_list[count*d+j,8] - new_info_list[count*d+j,5]
          max_err_1 <- max(max_err_1, err_1)
          max_err_2 <- max(max_err_2, err_2)
        }

        new_err_list[c(count,count+1)] <- c(max_err_1,max_err_2)
        count <- count+2

      } else {

        new_info_list[((count-1)*d+1):(count*d),] <- info
        new_err_list[count] <- err_list[i]
        count <- count+1
      }
    }

    info_list <- new_info_list
    err_list <- new_err_list
    max_err <- max(err_list)
  }
  ss <- info_list[d*1:(nrow(info_list)/d),1]
  ts <- info_list[d*1:(nrow(info_list)/d),2]
  xs <- matrix(info_list[,3],byrow=TRUE,ncol=2)
  ys <- matrix(info_list[,4],byrow=TRUE,ncol=2)
  centres <- matrix(0.5*(info_list[,8]+info_list[,5]),byrow=TRUE,ncol=d)
  errors <- err_list/2
  extra_info <- matrix(t(info_list[,3:8]),byrow=TRUE,ncol=6*d)
  fin_loc <- info_list[(nrow(info_list)-d+1):nrow(info_list),4]

  out <- list(ss=ss,ts=ts,xs=xs,ys=ys,centres=centres,errors=errors,fin_loc=fin_loc,extra_info=extra_info)
  return(out)
}

esbm_print <- function(x, s, t, epsilon){

  d <- length(x)
  y <- rnorm(d,x,t-s)

  info_list <- matrix(0,nrow=d,ncol=8)
  max_err <- 0
  for(j in 1:d){
    info_list[j,] <- eadl(s,t,x[j],y[j])$layer
    max_err <- max(max_err, info_list[j,8] - info_list[j,5])
  }

  err_list <- as.vector(max_err)

  while(max_err > 2*epsilon){

    length <- nrow(info_list)/d
    new_length <- length + sum(err_list > epsilon)
    count <- 1

    new_info_list <- matrix(0,nrow=d*new_length,ncol=8)
    new_err_list <- numeric(new_length)

    for(i in 1:length){

      info <- info_list[((i-1)*d+1):(i*d),]

      if(err_list[i] > epsilon){

        max_err_1 <- max_err_2 <- 0
        for(j in 1:d){
          midpoint <- earp(info[j,1],0.5*(info[j,1]+info[j,2]),info[j,2],info[j,3],info[j,4],info[j,5],info[j,6],info[j,7],info[j,8])$draw
          split_info <- eabl(info[j,1],0.5*(info[j,1]+info[j,2]),info[j,2],info[j,3],midpoint,info[j,4],info[j,5],info[j,6],info[j,7],info[j,8])$lyr2
          new_info_list[c((count-1)*d+j,count*d+j),] <- t(apply(split_info,1,earl2))

          err_1 <- new_info_list[(count-1)*d+j,8] - new_info_list[(count-1)*d+j,5]
          err_2 <- new_info_list[count*d+j,8] - new_info_list[count*d+j,5]
          max_err_1 <- max(max_err_1, err_1)
          max_err_2 <- max(max_err_2, err_2)
        }

        new_err_list[c(count,count+1)] <- c(max_err_1,max_err_2)
        count <- count+2

      } else {

        new_info_list[((count-1)*d+1):(count*d),] <- info
        new_err_list[count] <- err_list[i]
        count <- count+1
      }
    }

    info_list <- new_info_list
    err_list <- new_err_list
    max_err <- max(err_list)
  }

  centres <- matrix(0.5*(info_list[,8]+info_list[,5]),byrow=TRUE,ncol=d)
  errors <- err_list

  minx <- min(info_list[seq(1,(nrow(info_list)-1),by=2),3])
  miny <- min(info_list[seq(2,nrow(info_list),by=2),3])
  maxx <- max(info_list[seq(1,(nrow(info_list)-1),by=2),3])
  maxy <- max(info_list[seq(2,nrow(info_list),by=2),3])
  plot(c(minx-1, maxx+1), c(miny-1,maxy+1), type = "n",asp=1)
  rect(centres[,1]-err_list,centres[,2]-err_list,
       centres[,1]+err_list,centres[,2]+err_list)
  abline(h=0,col="red"); abline(v=0,col="red")
  lines(info_list[seq(1,(nrow(info_list)-1),by=2),3],info_list[seq(2,nrow(info_list),by=2),3],col="blue")

  return(err_list)
}

esbm_test <- function(x, s, t, epsilon){

  d <- length(x)
  y <- rnorm(d,x,t-s)

  info_list <- matrix(0,nrow=d,ncol=8)
  max_err <- 0
  for(j in 1:d){
    info_list[j,] <- eadl(s,t,x[j],y[j])$layer
    max_err <- max(max_err, info_list[j,8] - info_list[j,5])
  }

  err_list <- as.vector(max_err)

  while(max_err > 2*epsilon){

    length <- nrow(info_list)/d
    new_length <- length + sum(err_list > epsilon)
    count <- 1

    new_info_list <- matrix(0,nrow=d*new_length,ncol=8)
    new_err_list <- numeric(new_length)

    for(i in 1:length){

      info <- info_list[((i-1)*d+1):(i*d),]

      if(err_list[i] > epsilon){

        max_err_1 <- max_err_2 <- 0
        for(j in 1:d){
          midpoint <- earp(info[j,1],0.5*(info[j,1]+info[j,2]),info[j,2],info[j,3],info[j,4],info[j,5],info[j,6],info[j,7],info[j,8])$draw
          split_info <- eabl(info[j,1],0.5*(info[j,1]+info[j,2]),info[j,2],info[j,3],midpoint,info[j,4],info[j,5],info[j,6],info[j,7],info[j,8])$lyr2
          new_info_list[c((count-1)*d+j,count*d+j),] <- t(apply(split_info,1,earl2))

          err_1 <- new_info_list[(count-1)*d+j,8] - new_info_list[(count-1)*d+j,5]
          err_2 <- new_info_list[count*d+j,8] - new_info_list[count*d+j,5]
          max_err_1 <- max(max_err_1, err_1)
          max_err_2 <- max(max_err_2, err_2)
        }

        new_err_list[c(count,count+1)] <- c(max_err_1,max_err_2)
        count <- count+2

      } else {

        new_info_list[((count-1)*d+1):(count*d),] <- info
        new_err_list[count] <- err_list[i]
        count <- count+1
      }
    }

    info_list <- new_info_list
    err_list <- new_err_list
    max_err <- max(err_list)
  }

  check1 <- sum((info_list[,3] > info_list[,5]) & (info_list[,3]<info_list[,8]))
  check2 <- sum((info_list[,4] > info_list[,5]) & (info_list[,4]<info_list[,8]))
  N <- nrow(info_list)

  return(N-c(check1,check2))
}


