###############
#COllecting data from output files
###############

#for example 2
gather <- function(folder_name){

  filenames <- list.files(folder_name, full.names=TRUE)
  L <- length(filenames)
  record <- numeric(L)

  for(i in 1:L){

    load(filenames[[i]])
    record[i] <- v[1]
  }

  return(record)
}

gather_record <- function(folder_name, final_level){

  filenames <- list.files(folder_name, full.names=TRUE)
  L <- length(filenames)
  p_est<- numeric(L)

  for(k in 1:L){

    load(filenames[[k]])
    p_est[k] <- prod(N_surv[1:final_level]/100)  }

  return(p_est)
}

#for euler schemes
gather_euler <- function(no_trials,step_size){

  record <- numeric(no_trials)
  for(j in 1:no_trials){
    load(paste(paste("euler-results/step_size",step_size,"task",j,sep="_"),".RData",sep=""))
    record[j] <- p_estimates
  }

  return(record)
}

gather_euler_record <- function(no_trials,step_size){

  record <- matrix(0,row=no_trials,ncol=17)

  for(j in 1:no_trials){
    load(paste(paste("euler-results/step_size",step_size,"task",j,sep="_"),".RData",sep=""))
    record[j,] <- survival_record
  }

  return(record)
}

#for example 1
gather_ex2 <- function(folder_name){

  filenames <- list.files(folder_name, full.names=TRUE)
  L <- length(filenames)
  record <- numeric(L)

  for(k in 1:L){

    load(filenames[[k]])
    record[k] <- out
  }

  return(record)
}

