##################################
#Euler schemes
##################################

###########################
#1. Brownian motion
###########################

###########################
#1a. One-dimensional
###########################

cd_euler <- function(x, delta, a, b){

  while((a <= x) & (x <= b)){
    x <- x + rnorm(1, mean = 0, sd = sqrt(delta))
  }

  if(x < a){
    return(0)
  } else {
    return(1)
  }
}

splitting_euler <- function(n, levels, replicates, x, delta){

  survivors <- n

  for (i in 1:length(levels)){
    print(c("i",i))
    count <- 0
    for (j in 1:survivors){
      print(c("j",j))
      count <- count + cd_euler(x, delta, 0, levels[i])
    }
    if(count == 0){
      return(0)
    }
    if(i < length(levels)){
      survivors <- replicates[i] * count
      delta <- replicates[i]*delta
    }
    x <- levels[i]
  }
  denom <- prod(replicates) * n
  return(count/denom)
}

###########################
#1b. Multi-dimensional
###########################

cdm_euler <- function(x, delta, a, b, reac_coord){

  d <- length(x)

  while((a <= reac_coord(x)) & (reac_coord(x) <= b)){
    x <- x + rnorm(d, mean = 0, sd = sqrt(delta))
  }

  if(reac_coord(x) < a){
    return(list(survived=0))
  } else {
    return(list(survived=1,x=x))
  }
}

msplitting_euler <- function(n, z_A, levels, x, reac_coord, delta, delta_scale){

  d <- length(x)
  m <- length(levels)
  survivors <- matrix(rep(x,n),byrow=TRUE,nrow=n)
  n_surv <- integer(m)

  for (i in 1:m){
    new_survivors <- matrix(0,nrow=0,ncol=d)

    for (j in 1:n){
      trial <- cdm_euler(survivors[j,], delta, z_A, levels[i], reac_coord)
      if(trial$survived==1){
        new_survivors <- rbind(new_survivors,trial$x)
      }
    }

    n_surv[i] <- nrow(new_survivors)

    if(n_surv[i] == 0){
      return(0)
    }

    if(i < m){
      indices <- sample(1:n_surv[i], n, replace=TRUE)
      survivors <- new_survivors[indices,]
      delta <- delta*delta_scale[i]
    }
  }

  p_est <- prod(n_surv) / n^m
  return(list(p_est=p_est, survival_record=n_surv))
}

msplitting_euler2 <- function(n, z_A, levels, x, reac_coord, delta, delta_scale){

  d <- length(x)
  m <- length(levels)
  survivors <- matrix(rep(x,n),byrow=TRUE,nrow=n)
  n_surv <- integer(m)

  for (i in 1:m){
    new_survivors <- matrix(0,nrow=0,ncol=d)

    for (j in 1:n){
      trial <- cdm_euler(survivors[j,], delta, z_A, levels[i], reac_coord)
      if(trial$survived==1){
        new_survivors <- rbind(new_survivors,trial$x)
      }
    }

    n_surv[i] <- nrow(new_survivors)

    if(n_surv[i] == 0){
      return(0)
    }

    if(i < m){
      indices <- sample(1:n_surv[i], n, replace=TRUE)
      survivors <- new_survivors[indices,]
      delta <- delta*delta_scale[i]
    }
  }
  out <- prod(n_surv) / cumprod(rep(n,m))
  return(out)
}


