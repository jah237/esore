##########################
#1D Multilevel Splitting
##########################

###############
#0. Alternating series function
#############

#eabeta in simplest form (cf `eabetaC' in murrays_code.R)
eabeta <- function(m,s,t,x,y,Ll,Lu,Ul,Uu){
  if(m>=mpfrthr){
    pbn<-m*mpfrpbn
    s<-mpfr(s,precBits=pbn);t<-mpfr(t,precBits=pbn)
    x<-mpfr(x,precBits=pbn);y<-mpfr(y,precBits=pbn)
    Ll<-mpfr(Ll,precBits=pbn);Lu<-mpfr(Lu,precBits=pbn)
    Ul<-mpfr(Ul,precBits=pbn);Uu<-mpfr(Uu,precBits=pbn)
  }

  z1<-eazeta(m+1,s,t,x,y,Ll,Uu); z2<-eazeta(m,s,t,x,y,Lu,Uu)
  z3<-eazeta(m,s,t,x,y,Ll,Ul); z4<-eazeta(m+1,s,t,x,y,Lu,Ul)

  out<-as.numeric(-z1+z2+z3-z4)

  return(out)
}

#################
#1. Brownian bridge operations
#################

#midpoint and new layers
transect <- function(info){

  s <- info[1]; t <- info[2]; q <- (s+t)/2
  x <- info[3]; y <- info[4]
  Ll <- info[5]; Lu <- info[6]
  Ul <- info[7]; Uu <- info[8]

  w <- earp(s,q,t,x,y,Ll,Lu,Ul,Uu)$draw

  Lu <- min(w,Lu); Ul <- max(w,Ul)

  out <- eabl(s,q,t,x,w,y,Ll,Lu,Ul,Uu)$lyr2
  return(out)
}

#alternating series for refining upper layer
upper_ratio_m <- function(info, U_star, n){

  s <- info[1]; t <- info[2]
  x <- info[3]; y <- info[4]
  Ll <- info[5]; Lu <- info[6]
  Ul <- info[7]; Uu <- info[8]

  #check this is right
  denom <- eabeta(n+1,s,t,x,y,Ll,Lu,Ul,Uu)
  out <- eabeta(n,s,t,x,y,Ll,Lu,U_star,Uu)/denom

  return(out)
}

#alternating series for refining lower layer
lower_ratio_m <- function(info, L_star, n){

  s <- info[1]; t <- info[2]
  x <- info[3]; y <- info[4]
  Ll <- info[5]; Lu <- info[6]
  Ul <- info[7]; Uu <- info[8]

  denom <- eabeta(n+1,s,t,x,y,Ll,Lu,Ul,Uu)
  out <- eabeta(n,s,t,x,y,Ll,L_star,Ul,Uu)/denom

  return(out)
}

#refine upper layer once at chosen target
refine_upper_layer_m <- function(info, U_star){

  s <- info[1]; t <- info[2]
  x <- info[3]; y <- info[4]
  L_down <- info[5]; L_up <- info[6]
  U_down <- info[7]; U_up <- info[8]

  decision <- retro_bernoulli(upper_ratio_m, info = info, U_star = U_star)

  if(decision == 1){

    info[7] <- U_star
  } else {

    info[8] <- U_star
  }

  return(info)
}

#refine lower layer once at chosen target
refine_lower_layer_m <- function(info, L_star){

  decision <- retro_bernoulli(lower_ratio_m, info = info, L_star = L_star)

  if(decision == 1){

    info[6] <- L_star
  } else {

    info[5] <- L_star
  }

  return(info)
}

###############
#2. Multilevel splitting
##############

#check crossing for single barrier over given time-interval
crossing_decision <- function(info, a, b){

  Ll <- info[5]; Lu <- info[6]
  Ul <- info[7]; Uu <- info[8]

  if((Ll < a) & (a < Lu)){
    info <- refine_lower_layer_m(info, L_star = a)
  }

  if((Ul < b) & (b < Uu)){
    info <- refine_upper_layer_m(info, U_star = b)
  }

  Ll <- info[5]; Lu <- info[6]
  Ul <- info[7]; Uu <- info[8]

  if((Lu <= a) & (Uu <= b)){

    return(list(decision=-1))
  } else if((a <= Ll) & (b <= Ul)){

    return(list(decision=1))
  } else if((a <= Ll) & (Uu <= b)){

    return(list(decision=0))
  } else if((Lu <= a) & ((b <= Ul))){

    return(list(decision=2,info=info))
  } else {

    return("disaster!")
  }
}

#simulate crossing for single barrier until ocurrence
resolve <- function(x, a, b, plot = FALSE){

  s <- 0
  #set length over time interval over which to extend Brownian path
  t_hor <- min((x-a)^2, (x-b)^2)
  t <- t_hor
  y <- rnorm(1, x, sd=sqrt(t))
  info <- matrix(eadl(s,t,x,y)$layer, nrow = 1)

  while(TRUE){

    if(nrow(info) == 1){

      if(plot == TRUE){
        plot(info[1:2],info[3:4], xlim = c(info[1]-1,info[2]+1), ylim = c(min(info[5],a)-1,max(info[8],b)+1),
           main ="", type = "l")

        points(info[1:2], info[3:4], col = "orange", pch = 16)
        abline(h=c(a,b), col="blue")
        rect(info[1], info[5], info[2], info[6], density=15, col = "purple", border = TRUE)
        rect(info[1], info[7], info[2], info[8], density=15, col = "purple", border = TRUE)
      }

      decision_list <- crossing_decision(info, a, b)
      decision <- decision_list$decision

      if((decision == 1) | (decision == -1)){

        return(decision)
      } else if(decision == 0){

        s <- info[2]
        x <- info[4]
        t <- s + t_hor
        y <- rnorm(1, x, sqrt(t-s))
        info <- matrix(eadl(s,t,x,y)$layer, nrow = 1)

      } else {

        info <- decision_list$info
        info <- transect(info)


      }
    } else if(nrow(info) == 2){

      times <- c(info[1,1:2],info[2,2])
      locs <- c(info[1,3:4],info[2,4])

      if(plot == TRUE){
        plot(times,locs, xlim = c(info[1,1]-1,info[2,2]+1), ylim = c(min(c(info[,5],a))-1,max(c(info[,8],b))+1),
           main ="", type = "l")

        points(times, locs, col = "orange", pch = 16)
        abline(h=c(a,b), col="blue")
        rect(info[1,1], info[1,5], info[1,2], info[1,6], density=15, col = "purple", border = TRUE)
        rect(info[1,1], info[1,7], info[1,2], info[1,8], density=15, col = "purple", border = TRUE)
        rect(info[2,1], info[2,5], info[2,2], info[2,6], density=15, col = "purple", border = TRUE)
        rect(info[2,1], info[2,7], info[2,2], info[2,8], density=15, col = "purple", border = TRUE)
      }

      for(j in 1:2){

        info_slice <- info[j, ]

        decision_list <- crossing_decision(info_slice, a, b)
        decision <- decision_list$decision

        if((decision == 1) | (decision == -1)){

          return(decision)
        } else if (decision == 0){

          if(j == 1){
            #run on to j=2
          } else if(j == 2){

            s <- info_slice[2]
            x <- info_slice[4]
            t <- s + t_hor
            y <- rnorm(1, x, sqrt(t-s))
            info <- matrix(eadl(s,t,x,y)$layer, nrow = 1)

          }
        } else {

          info <- decision_list$info
          info <- transect(info)
          break
        }
      }
    }
  }
}

#Multilevel splitting
simple_splitting <- function(n, levels, replicates, x){

  N <- length(levels)
  survivors <- n

  for (i in 1:N){

    print(c("i",i))

    b <- levels[i]
    count <- 0

    for (j in 1:survivors){

      if(j%%10 == 0){
        print(c("j",j))
      }

      #note "resolve" returns +/- 1
      count <- count + 1/2 * (resolve(x, 0, b) + 1)
    }

    if(count == 0){
      return(0)
    }

    if (i < N){
      survivors <- replicates[i]*count
      x <- b
    }
  }

  denom <- prod(replicates)*n
  return(count/denom)
}

#MLS calibrated for paths starting x_0=1
ex2 <- function(no_particles, base, no_levels){

  levels <- base^(1:no_levels)
  replicates <- rep(base,no_levels-1)
  x_0 <- 1

  out <- simple_splitting(no_particles, levels, replicates, x_0)
  return(out)
}


