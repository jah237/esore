################
#Extend finished simulation runs, and complete halted ones
################

###############
#1. Complete
###############

complete_splitting_smc <- function(id){

  task <- id

  filename <- paste("pick-up-points/task_", task, ".RData", sep="")
  load(filename)

  .Random.seed <<- saved_seed
  jx <- j
  ix <- i

  if(jx < N){
    for(j in jx:N){

      print(c("j",j))

      if(max_levels[j] > ix){

        dbl_crossings <- dbl_crossings+1

        survivors_location <- rbind(survivors_location, starting_locations[j,])
        survivors_max_level <- c(survivors_max_level, max_levels[j])
        survivors_time <- c(survivors_time, times[j])
      } else {

        x <- starting_locations[j,]
        s <- times[j]


        if(save_seed==TRUE){
          saved_seed <- .Random.seed
          name<- paste(paste("more-data/task",id,"level",ix,"particle",j,sep="_"),".RData",sep="")
          save(list=ls(), file=name)
        }
        #here method has been specified as "exsp" to avoid conflict with
        #version of "method" etc in loaded data file (else changes to esbm
        #don't register) - should be a neater fix...

        trial <- extended_crossing2(x, s, s+t, epsilon_0, esbm, esbm_cond,
                                    z_A, levels[i:m], inf_coord, sup_coord,
                                    babble, ix, m)

        if(trial$survived == 1){
          survivors_location <- rbind(survivors_location, trial$fin_loc)
          survivors_max_level <- c(survivors_max_level, trial$max_level + ix)
          survivors_time <- c(survivors_time, trial$t_fin)
        }
      }
    }
  }


  if(ix < m){
    for(i in (ix+1):m){

      print(c("i",i))
      z_B <- levels[i]
      survivors_location <- matrix(0, ncol = d, nrow = 0)
      survivors_max_level <- numeric(0)
      survivors_time <- numeric(0)
      N_new <- 1

      for(j in 1:N){

        print(c("j",j))

        if(max_levels[j] > i){

          dbl_crossings <- dbl_crossings+1

          survivors_location <- rbind(survivors_location, starting_locations[j,])
          survivors_max_level <- c(survivors_max_level, max_levels[j])
          survivors_time <- c(survivors_time, times[j])
        } else {

          x <- starting_locations[j,]
          s <- times[j]


          if(save_seed==TRUE){
            saved_seed <- .Random.seed
            name<- paste(paste("data/task",id,"level",i,"particle",j,sep="_"),".RData",sep="")
            save(list=ls(), file=name)
          }

          #SEE analogous line above
          trial <- extended_crossing2(x, s, s+t, epsilon_0, esbm, esbm_cond,
                                      z_A, levels[i:m], inf_coord, sup_coord,
                                      babble, i, m)

          if(trial$survived == 1){
            survivors_location <- rbind(survivors_location, trial$fin_loc)
            survivors_max_level <- c(survivors_max_level, trial$max_level + i)
            survivors_time <- c(survivors_time, trial$t_fin)
          }
        }
      }

      N_surv[i] <- length(survivors_time)
      print("n_surv");print(N_surv[i])

      if(i < m){
        if(N_surv[i] == 0){
          return(0)
        } else {
          indices <- sample(N_surv[i], N, replace=TRUE)
          starting_locations <- survivors_location[indices,]
          times <- survivors_time[indices]
          max_levels <- survivors_max_level[indices]
        }
        epsilon_0 <- epsilon_scale[i]*epsilon_0
      }
    }
  }

  print("final N_surv"); print(N_surv)
  print("dbl_crossings");   print(dbl_crossings)

  out <- prod(N_surv)/N^m

  name<- paste(paste("more-results/task",task,sep="_"),".RData",sep="")
  save(list=ls(), file=name)

  return(out)
}

#####################
#2. Extend
#####################

extend_splitting_smc <- function(id, prev_level,target_level){

  task <- id

  filename <- paste(paste("data/task",task,"level",prev_level,"particle",100,sep="_"),".RData", sep="")
  load(filename)

  .Random.seed <<- saved_seed
  ix <- prev_level
  m <- target_level
  levels <- 2^(0.5*(2:(m+1)))

  for(i in (ix+1):m){

    print(c("i",i))
    z_B <- levels[i]
    survivors_location <- matrix(0, ncol = d, nrow = 0)
    survivors_max_level <- numeric(0)
    survivors_time <- numeric(0)
    N_new <- 1

    for(j in 1:N){

      print(c("j",j))

      if(max_levels[j] > i){

        dbl_crossings <- dbl_crossings+1

        survivors_location <- rbind(survivors_location, starting_locations[j,])
        survivors_max_level <- c(survivors_max_level, max_levels[j])
        survivors_time <- c(survivors_time, times[j])
      } else {

        x <- starting_locations[j,]
        s <- times[j]


        if(save_seed==TRUE){

          if((j==1) & (i>prev_level+1)){
            unlink(paste(paste("more-data/task",task,"level",i-1,"particle",N,sep="_"),".RData",sep=""))
          } else if(j>1){
            unlink(paste(paste("more-data/task",task,"level",i,"particle",j-1,sep="_"),".RData",sep=""))
          }

          saved_seed <- .Random.seed
          name<- paste(paste("more-data/task",id,"level",i,"particle",j,sep="_"),".RData",sep="")
          save(list=ls(), file=name)
        }

        #SEE analogous line above
        trial <- extended_crossing2(x, s, s+t, epsilon_0, esbm, esbm_cond,
                                    z_A, levels[i:m], inf_coord, sup_coord,
                                    babble, i, m)

        if(trial$survived == 1){
          print("success!")
          survivors_location <- rbind(survivors_location, trial$fin_loc)
          survivors_max_level <- c(survivors_max_level, trial$max_level + i)
          survivors_time <- c(survivors_time, trial$t_fin)
        }
      }
    }

    N_surv[i] <- length(survivors_time)
    print("n_surv");print(N_surv[i])

    if(i < m){
      if(N_surv[i] == 0){
        return(0)
      } else {
        indices <- sample(N_surv[i], N, replace=TRUE)
        starting_locations <- survivors_location[indices,]
        times <- survivors_time[indices]
        max_levels <- survivors_max_level[indices]
      }
      epsilon_0 <- epsilon_scale[i]*epsilon_0
    }
  }

  print("final N_surv"); print(N_surv)
  print("dbl_crossings");   print(dbl_crossings)

  out <- prod(N_surv)/N^m

  name<- paste(paste("more-results/task",task,sep="_"),".RData",sep="")
  save(list=ls(), file=name)

  return(out)
}


