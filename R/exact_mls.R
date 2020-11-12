##################
#1. Plotting functions
##################

plot_balls <- function(centres,errors,r,m){

  plot(c(0,0),xlim=c(-2,2^(0.5*(m+2))),ylim=c(-2,2^(0.5*(m+2))),type="n")
  rect(centres[,1]-errors,centres[,2]-errors,
       centres[,1]+errors,centres[,2]+errors)
  abline(h=0,col="red"); abline(v=0,col="red")
  for(l in r:m){
    segments(x0=2^(0.5*(l+1)),y0=2^(0.5*(l+1)),x1=100000,col="green")
    segments(x0=2^(0.5*(l+1)),y0=2^(0.5*(l+1)),y1=100000,col="green")
  }
}

plot_all <- function(xs,centres,errors,r,m){

  plot(xs[,1],xs[,2],xlim=c(-2,2^(0.5*(m+2))), ylim = c(-2,2^(0.5*(m+2))),type="l",col="blue")
  rect(centres[,1]-errors,centres[,2]-errors,
       centres[,1]+errors,centres[,2]+errors)
  abline(h=0,col="red"); abline(v=0,col="red")
  for(l in r:m){
    segments(x0=2^(0.5*(l+1)),y0=2^(0.5*(l+1)),x1=100000,col="green")
    segments(x0=2^(0.5*(l+1)),y0=2^(0.5*(l+1)),y1=100000,col="green")
  }
}

###################
#2. Crossing decisions
###################

#positions for one-sided barrier
mark_pos_basic2 <- function(inf_coord,sup_coord,centre,error,z,dir){

  inf <- inf_coord(centre,error)
  sup <- sup_coord(centre,error)
  degree <- 0

  if(sup < z){
    position <- 1
  } else if(inf > z){
    position <- 3
  } else {
    position <- 2
    if(dir=="up"){
      degree <- (sup-z)/(sup-inf)
    } else if(dir == "down") {
      degree <- (z-inf)/(sup-inf)
    }
  }
  return(list(position=position,degree=degree))
}

#positions for two sided barrier
mark_positions2 <- function(inf_coord,sup_coord,centres,errors,z_A,z_B){

  if(is.matrix(centres)==0){
    centres <- matrix(centres, nrow=1, byrow=TRUE)
  }

  n <- nrow(centres)

  positions <- numeric(n)
  degrees <- numeric(n)

  for(j in 1:n){

    inf <- inf_coord(centres[j,],errors[j])
    sup <- sup_coord(centres[j,],errors[j])

    if(sup < z_A){
      positions[j] <- 1
    } else if(inf < z_A){
      positions[j] <- 2
      degrees[j] <- (z_A - inf)/(sup-inf)
    } else if(sup < z_B){
      positions[j] <- 3
    } else if(inf < z_B){
      positions[j] <- 4
      degrees[j] <- (sup - z_B)/(sup-inf)
    } else {
      positions[j] <- 5
    }
  }
  return(list(positions=positions,degrees=degrees))
}

split_into_blocks<- function(positions, degrees, info_list){

  centres <- info_list$centres
  errors <- info_list$errors
  ss <- info_list$ss
  ts <- info_list$ts
  extra_info <- info_list$extra_info
  K <- length(positions)

  block_list <- vector(length=0,mode='list')
  degree_list <- vector(length=0,mode='list')
  centre_list <- vector(length=0,mode='list')
  error_list <- vector(length=0,mode='list')
  ss_list <- vector(length=0,mode='list')
  ts_list <- vector(length=0,mode='list')
  extra_info_list <- vector(length=0,mode='list')

  no_blocks <- 0

  if((any(positions<=2)) | (any(positions>=4))){

    i <- min(which((positions<=2) | (positions>=4)))
    type <- positions[i]

    while(TRUE){

      no_blocks <- no_blocks+1

      if(type==2){
        if(any(positions>=4)){
          i <- min(which(positions>=4))
          block_list[[no_blocks]] <- positions[1:(i-1)]
          degree_list[[no_blocks]] <- degrees[1:(i-1)]
          centre_list[[no_blocks]] <- centres[1:(i-1),]
          error_list[[no_blocks]] <- errors[1:(i-1)]
          ss_list[[no_blocks]] <- ss[1:(i-1)]
          ts_list[[no_blocks]] <- ts[1:(i-1)]
          extra_info_list[[no_blocks]] <- extra_info[1:(i-1),]

          positions <- positions[i:K]
          degrees <- degrees[i:K]
          centres <- centres[i:K,]
          errors <- errors[i:K]
          ss <- ss[i:K]
          ts <- ts[i:K]
          extra_info <- extra_info[i:K,]
          K < length(positions)

          type <- 4

        } else {

          block_list[[no_blocks]] <- positions
          degree_list[[no_blocks]] <- degrees
          centre_list[[no_blocks]] <- centres
          error_list[[no_blocks]] <- errors
          ss_list[[no_blocks]] <- ss
          ts_list[[no_blocks]] <- ts
          extra_info_list[[no_blocks]] <- extra_info

          return(list(block_list=block_list, degree_list=degree_list,
                      centre_list=centre_list, error_list=error_list,
                      ss_list=ss_list, ts_list=ts_list,
                      extra_info_list=extra_info_list))
        }
      } else if(any(positions<=2)){

        i <- min(which(positions<=2))
        block_list[[no_blocks]] <- positions[1:(i-1)]
        degree_list[[no_blocks]] <- degrees[i:(i-1)]
        centre_list[[no_blocks]] <- centres[1:(i-1),]
        error_list[[no_blocks]] <- errors[1:(i-1)]
        ss_list[[no_blocks]] <- ss[1:(i-1)]
        ts_list[[no_blocks]] <- ts[1:(i-1)]
        extra_info_list[[no_blocks]] <- extra_info[1:(i-1),]

        positions <- positions[i:K]
        degrees <- degrees[i:K]
        centres <- centres[i:K,]
        errors <- errors[i:K]
        ss <- ss[i:K]
        ts <- ts[i:K]
        extra_info <- extra_info[i:K,]
        K < length(positions)

        type <- 2

      } else {
        block_list[[no_blocks]] <- positions
        degree_list[[no_blocks]] <- degrees
        centre_list[[no_blocks]] <- centres
        error_list[[no_blocks]] <- errors
        ss_list[[no_blocks]] <- ss
        ts_list[[no_blocks]] <- ts
        extra_info_list[[no_blocks]] <- extra_info

        return(list(block_list=block_list, degree_list=degree_list,
                    centre_list=centre_list, error_list=error_list,
                    ss_list=ss_list, ts_list=ts_list,
                    extra_info_list=extra_info_list))
      }
    }
  } else {
    no_blocks <- 1
    block_list[[no_blocks]] <- positions
    degree_list[[no_blocks]] <- degrees
    centre_list[[no_blocks]] <- centres
    error_list[[no_blocks]] <- errors
    ss_list[[no_blocks]] <- ss
    ts_list[[no_blocks]] <- ts
    extra_info_list[[no_blocks]] <- extra_info

    return(list(block_list=block_list, degree_list=degree_list,
                centre_list=centre_list, error_list=error_list,
                ss_list=ss_list, ts_list=ts_list,
                extra_info_list=extra_info_list))
  }
}


crossing_below3 <- function(pos,degrees,centres,errors,ss,ts,extra_info,
                            method_cond,z_A,inf_coord,sup_coord, babble=0,r,m){

  K <- length(pos)
  if(K==0){
    return(-1)
  }

  if(is.matrix(centres)==0){
    centres <- matrix(centres,nrow=1,byrow=TRUE)
  }

  d <- ncol(centres)

  while(TRUE){

    if(is.matrix(centres)==0){
      centres <- matrix(centres,nrow=1,byrow=TRUE)
    }

    if(is.matrix(extra_info)==0){
      extra_info <- matrix(extra_info,nrow=1)
    }

    if(babble==1){

      plot_balls(centres,errors,r,m)
    }

    if(K==0){
      return(-1)

    } else if(sum(pos==1)>0){
      return(1)

    } else {
      i <- which.max(degrees)
      epsilon_new <- errors[i]/2
      new_skeleton <- method_cond(ss[i],ts[i],d,extra_info[i,],epsilon_new)
      new_centres <- as.matrix(new_skeleton$centres,ncol=d)
      K_new <- nrow(new_centres)
      new_errors <- new_skeleton$errors
      new_ss <- new_skeleton$ss
      new_ts <- new_skeleton$ts
      new_extra_info <- new_skeleton$extra_info

      new_positions <- numeric(K_new)
      new_degrees <- numeric(K_new)
      #move loop into mark position function
      for(j in 1:K_new){
        pos_update <- mark_pos_basic2(inf_coord,sup_coord,new_centres[j,],new_errors[j],z_A,dir="down")
        new_positions[j] <- pos_update$position
        new_degrees[j] <- pos_update$degree
      }
      indices <- which(new_positions < 3)

      if(K==1){
        centres <- new_centres[indices,]
        errors <- new_errors[indices]
        ss <- new_ss[indices]
        ts <- new_ts[indices]
        extra_info <- new_extra_info[indices,]
        pos <- new_positions[indices]
        degrees <- new_degrees[indices]
        K <- length(indices)

      } else {
        if(i == 1){
          centres <- rbind(new_centres[indices,], centres[2:K,])
          errors <- c(new_errors[indices],errors[2:K])
          ss <- c(new_ss[indices],ss[2:K])
          ts <- c(new_ts[indices],ts[2:K])
          extra_info <- rbind(new_extra_info[indices,],extra_info[2:K,])
          pos <- c(new_positions[indices],pos[2:K])
          degrees <- c(new_degrees[indices],degrees[2:K])
          K <- nrow(centres)
        } else if (i == K){
          centres <- rbind(centres[1:(K-1),],new_centres[indices,])
          errors <- c(errors[1:(K-1)], new_errors[indices])
          ss <- c(ss[1:(K-1)], new_ss[indices])
          ts <- c(ts[1:(K-1)], new_ts[indices])
          extra_info <- rbind(extra_info[1:(K-1),], new_extra_info[indices,])
          pos <- c(pos[1:(K-1)], new_positions[indices])
          degrees <- c(degrees[1:(K-1)], new_degrees[indices])
          K <- nrow(centres)
        } else {
          centres <- rbind(rbind(centres[1:(i-1),],new_centres[indices,]),centres[(i+1):K,])
          errors <- c(errors[1:(i-1)], new_errors[indices],errors[(i+1):K])
          ss <- c(ss[1:(i-1)], new_ss[indices], ss[(i+1):K])
          ts <- c(ts[1:(i-1)], new_ts[indices], ts[(i+1):K])
          extra_info <- rbind(rbind(extra_info[1:(i-1),], new_extra_info[indices,]),extra_info[(i+1):K,])
          pos <- c(pos[1:(i-1)], new_positions[indices], pos[(i+1):K])
          degrees <- c(degrees[1:(i-1)], new_degrees[indices], degrees[(i+1):K])
          K <- nrow(centres)
        }
      }
    }
  }
}

crossing_above4 <- function(pos,degrees,centres,errors,ss,ts,extra_info,fin_loc,
                            method_cond,z_B,inf_coord,sup_coord, babble=0,r,m){

  K <- length(pos)

  if(is.matrix(centres)==0){
    centres <- matrix(centres,nrow=1,byrow=TRUE)
  }

  d <- ncol(centres)

  while(TRUE){

    if(is.matrix(centres)==0){
      centres <- matrix(centres,nrow=1,byrow=TRUE)
    }

    if(is.matrix(extra_info)==0){
      extra_info <- matrix(extra_info,nrow=1)
    }

    if(babble==1){

      plot_balls(centres,errors,r,m)
    }

    if(K==0){
      return(list(crossed=-1))

    } else if(any(pos==3)){
      ind <- which(pos==3)[1]
      K <- length(pos)
      remaining_centres <- centres[ind:K,]
      remaining_errors <- errors[ind:K]
      remaining_ss <- ss[ind:K]
      remaining_ts <- ts[ind:K]
      remaining_extra_info <- extra_info[ind:K,]
      return(list(crossed=1,remaining_info = list(centres=remaining_centres,
                                                  errors=remaining_errors,
                                                  ss=remaining_ss,
                                                  ts=remaining_ts,
                                                  extra_info=remaining_extra_info,
                                                  fin_loc=fin_loc)))

    } else {
      i <- which.max(degrees)
      epsilon_new <- errors[i]/2
      new_skeleton <- method_cond(ss[i],ts[i],d,extra_info[i,],epsilon_new)
      new_centres <- as.matrix(new_skeleton$centres,ncol=d)
      K_new <- nrow(new_centres)
      new_errors <- new_skeleton$errors
      new_ss <- new_skeleton$ss
      new_ts <- new_skeleton$ts
      new_extra_info <- new_skeleton$extra_info

      new_positions <- numeric(K_new)
      new_degrees <- numeric(K_new)
      #move loop into mark position function
      for(j in 1:K_new){
        pos_update <- mark_pos_basic2(inf_coord,sup_coord,new_centres[j,],new_errors[j],z_B,dir="up")
        new_positions[j] <- pos_update$position
        new_degrees[j] <- pos_update$degree
      }
      indices <- which(new_positions > 1)

      if(K==1){
        centres <- new_centres[indices,]
        errors <- new_errors[indices]
        ss <- new_ss[indices]
        ts <- new_ts[indices]
        extra_info <- new_extra_info[indices,]
        pos <- new_positions[indices]
        degrees <- new_degrees[indices]
        K <- length(indices)

      } else {
        if(i == 1){
          centres <- rbind(new_centres[indices,], centres[2:K,])
          errors <- c(new_errors[indices],errors[2:K])
          ss <- c(new_ss[indices],ss[2:K])
          ts <- c(new_ts[indices],ts[2:K])
          extra_info <- rbind(new_extra_info[indices,],extra_info[2:K,])
          pos <- c(new_positions[indices],pos[2:K])
          degrees <- c(new_degrees[indices],degrees[2:K])
          K <- nrow(centres)
        } else if (i == K){
          centres <- rbind(centres[1:(K-1),],new_centres[indices,])
          errors <- c(errors[1:(K-1)], new_errors[indices])
          ss <- c(ss[1:(K-1)], new_ss[indices])
          ts <- c(ts[1:(K-1)], new_ts[indices])
          extra_info <- rbind(extra_info[1:(K-1),], new_extra_info[indices,])
          pos <- c(pos[1:(K-1)], new_positions[indices])
          degrees <- c(degrees[1:(K-1)], new_degrees[indices])
          K <- nrow(centres)
        } else {
          centres <- rbind(rbind(centres[1:(i-1),],new_centres[indices,]),centres[(i+1):K,])
          errors <- c(errors[1:(i-1)], new_errors[indices],errors[(i+1):K])
          ss <- c(ss[1:(i-1)], new_ss[indices], ss[(i+1):K])
          ts <- c(ts[1:(i-1)], new_ts[indices], ts[(i+1):K])
          extra_info <- rbind(rbind(extra_info[1:(i-1),], new_extra_info[indices,]),extra_info[(i+1):K,])
          pos <- c(pos[1:(i-1)], new_positions[indices], pos[(i+1):K])
          degrees <- c(degrees[1:(i-1)], new_degrees[indices], degrees[(i+1):K])
          K <- nrow(centres)
        }
      }
    }
  }
}

#Algorithm 4
crossing_direction_cond2 <- function(info_list, epsilon_0, method, method_cond,
                                    z_A, z_B, inf_coord, sup_coord, babble=0, r, m0){

  centres <- info_list$centres
  errors <- info_list$errors
  ss <- info_list$ss
  ts <- info_list$ts
  extra_info <- info_list$extra_info
  fin_loc <- info_list$fin_loc

  pos_info <- mark_positions2(inf_coord,sup_coord,centres,errors,z_A,z_B)
  positions <- pos_info$positions
  degrees <- pos_info$degrees

  d <- ncol(centres)
  K <- nrow(centres)
  t_fin <- ts[length(ts)]
  increment <- t_fin-ss[1]

  while(TRUE){

    if(is.matrix(centres)==0){
      centres <- matrix(centres,nrow=1,byrow=TRUE)
    }

    if(is.matrix(extra_info)==0){
      extra_info <- matrix(extra_info,nrow=1)
    }

    block_info <- split_into_blocks(positions, degrees, info_list)
    block_list <- block_info$block_list
    degree_list <- block_info$degree_list
    centre_list <- block_info$centre_list
    error_list <- block_info$error_list
    ss_list <- block_info$ss_list
    ts_list <- block_info$ts_list
    extra_info_list <- block_info$extra_info_list

    L <- length(block_list)

    for (j in 1:L){

      block_positions <- block_list[[j]]
      block_degrees <- degree_list[[j]]
      block_centres <- centre_list[[j]]
      block_errors <- error_list[[j]]
      block_ss <- ss_list[[j]]
      block_ts <- ts_list[[j]]
      block_extra_info <- extra_info_list[[j]]

      if(any(block_positions==1)){
        return(list(dir=-1))
      } else if(any(block_positions==5)){
        return(list(dir=+1,info_list=list(centres=centres, errors=errors,
                                          ss=ss, ts=ts, extra_info=extra_info, fin_loc=fin_loc)))
      } else if(any(block_positions==2)){

        trial <- crossing_below3(block_positions, block_degrees, block_centres,
                                 block_errors, block_ss, block_ts, block_extra_info,
                                 method_cond, z_A, inf_coord, sup_coord,
                                 babble, r, m0)

        if(trial==1){
          return(list(dir=-1))
        }
      } else if(any(block_positions==4)){

        trial <- crossing_above4(block_positions, block_degrees, block_centres,
                                 block_errors, block_ss, block_ts, block_extra_info,
                                 fin_loc, method_cond, z_B, inf_coord, sup_coord,
                                 babble, r, m0)

        outcome <- trial$crossed

        if(outcome==1){
          return(list(dir=1,info=trial$remaining_info))
        }
      }
    }

    x <- fin_loc
    s <- 0
    t <- max(increment,1)

    info_list <- method(x,s,t,epsilon_0)

    centres <- info_list$centres
    errors <- info_list$errors
    ss <- info_list$ss
    ts <- info_list$ts
    extra_info <- info_list$extra_info
    fin_loc <- info_list$fin_loc
    pos_info <- mark_positions2(inf_coord,sup_coord,centres,errors,z_A,z_B)
    positions <- pos_info$positions
    degrees <- pos_info$degrees

    K <- nrow(centres)
    t_fin <- ts[length(ts)]
  }
}

#######################
#3. Splitting
#######################

split_smc3 <- function(N, x_0, t, epsilon_0, method, method_cond,
                       z_A, levels, inf_coord, sup_coord, epsilon_scale, babble=0, save_seed=FALSE, id){

  d <- length(x_0)
  m <- length(levels)

  #initial paths
  continuing_paths <- vector(length=N,mode="list")
  for(j in 1:N){
    continuing_paths[[j]] <- method(x_0,0,t,epsilon_0)
  }

  N_surv <- numeric(m)

  for(i in 1:m){
    print(c("i",i))
    z_B <- levels[i]
    surviving_paths <- vector(length=N,mode="list")
    survivors <- 0

    for(j in 1:N){

      print(c("j",j))

      if(!exists(".Random.seed")){
        set.seed(NULL)
      }

      if(save_seed==TRUE){
        saved_seed <- .Random.seed
        if((j==1) & (i>1)){
          unlink(paste(paste("data/task",id,"level",i-1,"particle",N,sep="_"),".RData",sep=""))
        } else if(j>1){
          unlink(paste(paste("data/task",id,"level",i,"particle",j-1,sep="_"),".RData",sep=""))
        }

        name <- paste(paste("data/task",id,"level",i,"particle",j,sep="_"),".RData",sep="")
        save(list=ls(), file=name)
      }

      info_list <- continuing_paths[[j]]

      #crossing decision
      trial <- crossing_direction_cond2(info_list, epsilon_0, method, method_cond,
                                        z_A, z_B, inf_coord, sup_coord, babble, i, m)

      if(trial$dir== 1){
        survivors <- survivors + 1
        surviving_paths[[survivors]] <- trial$info
      }
    }

    print(c("i",i)); print(N_surv[i])
    N_surv[i] <- survivors

    #resampling
    if(i < m){
      if(N_surv[i] == 0){
        return(0)
      } else {
        indices <- sample(N_surv[i], N, replace=TRUE)
        continuing_paths <- surviving_paths[indices]
      }
      epsilon_0 <- epsilon_scale[i]*epsilon_0
    }
  }

  print("final N_surv"); print(N_surv)
  out <- prod(N_surv)/N^m
  return(out)
}


