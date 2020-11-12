#To be written as functions

ex1_euler <- function(no_trials,no_particles,no_levels,step_size,random_seed){
  
  set.seed(random_seed)

  levels <- 2^seq(1,by=0.5,length=no_levels)
  x = c(1.5,1.5)
  reac_coord <- function(x){min(x)}
  delta_scale <- rep(sqrt(2),no_levels-1)

  p_estimates <- numeric(no_trials)
  survival_record <- matrix(0,ncol=no_levels,nrow=no_trials)

  for(i in 1:no_trials){
    print(c("i",i))
    trial <- msplitting_euler(no_particles,0,levels,x,reac_coord,step_size,delta_scale)
    p_estimates[i] <- trial$p_est
    survival_record[i,] <- trial$survival_record
  }

  mean_survivals <- apply(survival_record,2,mean)

  name <- paste(paste("euler-results/step_size",step_size,"task",random_seed,sep="_"),".RData",sep="")
  save(list=ls(),file=name)

  return(list(survival_record = survival_record, mean_survivals=mean_survivals, p_estimates=p_estimates))
}

ex1 <- function(no_trials=1,no_particles,no_levels,save_seed,random_seed=1){

  n <- no_trials
  z_A <- 0
  levels <- 2^(seq(1,length=no_levels,by=0.5))
  x <- c(1.5,1.5)
  epsilon_scale <- rep(sqrt(2),times=no_levels-1)
  reac_coord <- function(x){
    return(min(x))
  }
  inf_coord <- function(centre,error){

    return(min(centre-error))
  }
  sup_coord <- function(centre,error){

    return(min(centre+error))
  }

  set.seed(random_seed)
  v <- numeric(no_trials)
  for(i in 1:no_trials){
    print(c("trial",i))
    v[i] <- split_smc3(no_particles,x,1,0.4,esbm,esbm_cond,0,levels,inf_coord,sup_coord,epsilon_scale,babble=0,save_seed,random_seed)
  }

  name<- paste(paste("results/task",random_seed,sep="_"),".RData",sep="")
  save(v, file=name)

  return(v)
}
