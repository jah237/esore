matsort <- function(mat,n) {mat[rank(mat[,n]),]<- mat[c(1:nrow(mat)),];return(mat)}
#### 0.1 - Euler Discretisation
euler 		<- function(mesh,x,s,t,drift,lambda,nu){deltat <- (t-s)/mesh; skel <- matrix(0,2,mesh+1); skel[1,] <- seq(s,t,deltat); skel[2,1] <- x; for(i in 2:(mesh+1)){skel[2,i]<-skel[2,i-1]+drift(skel[2,i-1])*deltat+rnorm(1,0,sd=(deltat)^(1/2))+rbinom(1,1,1-exp(-lambda(skel[2,i-1])*deltat))*nu(skel[2,i-1])}; y <- skel[2,dim(skel)[2]]
list(skel=skel,y=y)} # Euler Discretisation
##### 0.1.1 - Euler Discretisation - Multiples
eulersims	<- function(sims,mesh,x,s,t,drift,lambda,nu){output<-numeric(sims);for(i in 1:sims){output[i]<-euler(mesh,x,s,t,drift,lambda,nu)$y};output}

## 1 - Primary Functions, C functions - m odd only
### 1.0 - mpfr - high  precision kick in levels
mpfrthr		<- 9 # m threshold for high precision
mpfrpbn		<- 10 # m muliptlier for precision
### 1.1 - Brownian Bridge Functions (used by earpreg (no mpfr required as earpreg superceded at high precision))
easiga 		<- function(P,z,A,L){P*(A*(A+2*L-z)+L*(L-z))}; easigb<-function(P,z,A,L){P*(z-A-L)}
eaphia		<- function(P,z,A,si){P*(A^2-si*A*z)}; eaphib<-function(P,z,A,si){P*si*A}
### 1.2 - Zeta Functions (self contained in 1 - no mpfr used)
eazeta		<- function(n,s,t,x,y,L,U){if(max(x-U,y-U,L-x,L-y)>=0){1}else{j<-1:(ceiling(n/2));P<--2/(t-s);D<-U-L;D1<-D*j+L;D2<-D*j-U;z<-y-x;if(even(n)){sum(exp(P*(D1-x)*(D1-y))+exp(P*(D2+x)*(D2+y))-exp(P*j^2*D^2-P*j*D*z)-exp(P*j^2*D^2+P*j*D*z))}else{sum(exp(P*(D1-x)*(D1-y))+exp(P*(D2+x)*(D2+y)))-sum(exp(P*j[1:length(j)-1]^2*D^2-P*j[1:length(j)-1]*D*z)+exp(P*j[1:length(j)-1]^2*D^2+P*j[1:length(j)-1]*D*z))}}}
eazetaC		<- function(m,s,t,x,y,L,U){if(max(x-U,y-U,L-x,L-y)>=0){s1<-s2<-1}else{j<-1:((m+1)/2);P<--2/(t-s);D<-U-L;D1<-D*j+L;D2<-D*j-U;z<-y-x;s2<-sum(exp(P*(D1-x)*(D1-y))+exp(P*(D2+x)*(D2+y))-exp(P*j^2*D^2-P*j*D*z)-exp(P*j^2*D^2+P*j*D*z));s1<-s2+exp(P*((m+1)/2)^2*D^2-P*((m+1)/2)*D*z)+exp(P*((m+1)/2)^2*D^2+P*((m+1)/2)*D*z)};c(s1=s2,s2=s2)}
eaze3C		<- function(m,s,t,x,y,L,U){if(max(x-U,y-U,L-x,L-y)>=0){s1<-s2<-s3<-1}else{j<-1:((m+1)/2);P<--2/(t-s);D<-U-L;D1<-D*j+L;D2<-D*j-U;z<-y-x;s2<-sum(exp(P*(D1-x)*(D1-y))+exp(P*(D2+x)*(D2+y))-exp(P*j^2*D^2-P*j*D*z)-exp(P*j^2*D^2+P*j*D*z));s1<-s2+exp(P*((m+1)/2)^2*D^2-P*((m+1)/2)*D*z)+exp(P*((m+1)/2)^2*D^2+P*((m+1)/2)*D*z);s3<-s2+exp(P*(D*((m+3)/2)+L-x)*(D*((m+3)/2)+L-y))+exp(P*(D*((m+3)/2)-U+x)*(D*((m+3)/2)-U+y))};c(s1=s1,s2=s2,s3=s3)}
eaze4C		<- function(m,s,t,x,y,L,U){if(max(x-U,y-U,L-x,L-y)>=0){s1<-1;s2<-1;s3<-1;s4<-1}else{j<-1:((m+1)/2);P<--2/(t-s);D<-U-L;D1<-D*j+L;D2<-D*j-U;z<-y-x;s2<-sum(exp(P*(D1-x)*(D1-y))+exp(P*(D2+x)*(D2+y))-exp(P*j^2*D^2-P*j*D*z)-exp(P*j^2*D^2+P*j*D*z));s1<-s2+exp(P*((m+1)/2)^2*D^2-P*((m+1)/2)*D*z)+exp(P*((m+1)/2)^2*D^2+P*((m+1)/2)*D*z);s3<-s2+exp(P*(D*((m+3)/2)+L-x)*(D*((m+3)/2)+L-y))+exp(P*(D*((m+3)/2)-U+x)*(D*((m+3)/2)-U+y));s4<-s3-exp(P*((m+3)/2)^2*D^2-P*((m+3)/2)*D*z)-exp(P*((m+3)/2)^2*D^2+P*((m+3)/2)*D*z)};c(s1=s1,s2=s2,s3=s3,s4=s4)}
### 1.3 - Gamma Functions (mpfr used for eagammaC only as eagamma used within high precision functions only)
eagamma		<- function(n,s,t,x,y,L,U){1-eazeta(n,s,t,x,y,L,U)}
eagammaC	<- function(m,s,t,x,y,L,U){if(m>=mpfrthr){pbn<-m*mpfrpbn;s<-mpfr(s,precBits=pbn);t<-mpfr(t,precBits=pbn);x<-mpfr(x,precBits=pbn);y<-mpfr(y,precBits=pbn);L<-mpfr(L,precBits=pbn);U<-mpfr(U,precBits=pbn)};z<-eazetaC(m,s,t,x,y,L,U);c(s1=as.numeric(1-z[1]),s2=as.numeric(1-z[2]))}
### 1.4 - Bessel Functions (mpfr used for eadelC only)
eapsi		<- function(j,s,t,m,xoy,u){P<--2*abs(u-m)*j/(t-s);(2*abs(u-m)*j-(xoy-m))*exp(P*(abs(u-m)*j-(xoy-m)))}
eachi		<- function(j,s,t,m,xoy,u){P<--2*abs(u-m)*j/(t-s);(2*abs(u-m)*j+(xoy-m))*exp(P*(abs(u-m)*j+(xoy-m)))}
eadel2  	<- function(n,s,t,m,xoy,u){if(max(xoy-u,m-xoy)>=0){0}else{if(even(n)){j<-1:(n/2);1-(sum(eapsi(j,s,t,m,xoy,u)-eachi(j,s,t,m,xoy,u)))/(xoy-m)}else{if(n>1){j<-1:((n-1)/2);1-(sum(eapsi(j,s,t,m,xoy,u)-eachi(j,s,t,m,xoy,u)))/(xoy-m)-eapsi(max(j)+1,s,t,m,xoy,u)/(xoy-m)}else{1-eapsi(1,s,t,m,xoy,u)/(xoy-m)}}}}
eadel		<- function(n,s,t,x,y,m,u){if(x==m){xI<-1}else{xI<-0};if(y==m){yI<-1}else{yI<-0};if(max(xI,yI)==1){delT<-2}else{delT<-1};if(m>max(x,y)){x<--x;y<--y;m<--m;u<--u};if(max(x-u,y-u,m-x,m-y)>=0){out<-0};if(delT==1){out<-eagamma(n,s,t,x,y,m,u)/(1-exp(-2*(x-m)*(y-m)/(t-s)))};if(delT==2){if(xI*yI==0){xoy<-max(x,y);out<-eadel2(n,s,t,m,xoy,u)}else{out<-0}};if(out<0){out<-0};if(out>1){out<-1};out}
eadelC 		<- function(mt,s,t,x,y,m,u){if(mt>=mpfrthr){pbn<-mt*mpfrpbn;s<-mpfr(s,precBits=pbn);t<-mpfr(t,precBits=pbn);x<-mpfr(x,precBits=pbn);y<-mpfr(y,precBits=pbn);m<-mpfr(m,precBits=pbn);u<-mpfr(u,precBits=pbn)};c(s1=eadel(mt,s,t,x,y,m,u),s2=eadel(mt+1,s,t,x,y,m,u))}
### 1.5 - Beta Functions (mpfr used)
eabetaC		<- function(m,s,t,x,y,Ll,Lu,Ul,Uu){if(m>=mpfrthr){pbn<-m*mpfrpbn;s<-mpfr(s,precBits=pbn);t<-mpfr(t,precBits=pbn);x<-mpfr(x,precBits=pbn);y<-mpfr(y,precBits=pbn);Ll<-mpfr(Ll,precBits=pbn);Lu<-mpfr(Lu,precBits=pbn);Ul<-mpfr(Ul,precBits=pbn);Uu<-mpfr(Uu,precBits=pbn)};z1<-eazetaC(m,s,t,x,y,Ll,Uu); z2<-c(eazetaC(m,s,t,x,y,Lu,Uu)[2],eazetaC(m+2,s,t,x,y,Lu,Uu)[1]);z3<-c(eazetaC(m,s,t,x,y,Ll,Ul)[2],eazetaC(m+2,s,t,x,y,Ll,Ul)[1]);z4<-eazetaC(m,s,t,x,y,Lu,Ul);c(s1=as.numeric(-z1[1]+z2[1]+z3[1]-z4[1]),s2=as.numeric(-z1[2]+z2[2]+z3[2]-z4[2]))}
eabe3C		<- function(m,s,t,x,y,Ll,Lu,Ul,Uu){if(m>=mpfrthr){pbn<-m*mpfrpbn;s<-mpfr(s,precBits=pbn);t<-mpfr(t,precBits=pbn);x<-mpfr(x,precBits=pbn);y<-mpfr(y,precBits=pbn);Ll<-mpfr(Ll,precBits=pbn);Lu<-mpfr(Lu,precBits=pbn);Ul<-mpfr(Ul,precBits=pbn);Uu<-mpfr(Uu,precBits=pbn)};z1<-eaze3C(m,s,t,x,y,Ll,Uu); z2<-c(eaze3C(m,s,t,x,y,Lu,Uu)[2:3],eaze3C(m+2,s,t,x,y,Lu,Uu)[2]);z3<-c(eaze3C(m,s,t,x,y,Ll,Ul)[2:3],eaze3C(m+2,s,t,x,y,Ll,Ul)[2]);z4<-eaze3C(m,s,t,x,y,Lu,Ul);c(s1=as.numeric(-z1[1]+z2[1]+z3[1]-z4[1]),s2=as.numeric(-z1[2]+z2[2]+z3[2]-z4[2]),s3=as.numeric(-z1[3]+z2[3]+z3[3]-z4[3]))}
### 1.6 - Rho Functions (mpfr used)
earhoC		<- function(m,s,q,t,x,w,y,Ll,Lu,Ul,Uu){if(m>=mpfrthr){pbn<-m*mpfrpbn;s<-mpfr(s,precBits=pbn);q<-mpfr(q,precBits=pbn);t<-mpfr(t,precBits=pbn);x<-mpfr(x,precBits=pbn);w<-mpfr(w,precBits=pbn);y<-mpfr(y,precBits=pbn);Ll<-mpfr(Ll,precBits=pbn);Lu<-mpfr(Lu,precBits=pbn);Ul<-mpfr(Ul,precBits=pbn);Uu<-mpfr(Uu,precBits=pbn)}; z1L<-eaze3C(m,s,q,x,w,Ll,Uu); z1R<-eaze3C(m,q,t,w,y,Ll,Uu); z2L<-eaze3C(m,s,q,x,w,Lu,Uu); z2R<-eaze3C(m,q,t,w,y,Lu,Uu); z3L<-eaze3C(m,s,q,x,w,Ll,Ul); z3R<-eaze3C(m,q,t,w,y,Ll,Ul); z4L<-eaze3C(m,s,q,x,w,Lu,Ul); z4R<-eaze3C(m,q,t,w,y,Lu,Ul);c(s1=as.numeric(-z1L[2]-z1R[2]+z1L[1]*z1R[1]+z2L[1]+z2R[1]-z2L[2]*z2R[2]+z3L[1]+z3R[1]-z3L[2]*z3R[2]-z4L[2]-z4R[2]+z4L[1]*z4R[1]),s2=as.numeric(-z1L[3]-z1R[3]+z1L[2]*z1R[2]+z2L[2]+z2R[2]-z2L[3]*z2R[3]+z3L[2]+z3R[2]-z3L[3]*z3R[3]-z4L[3]-z4R[3]+z4L[2]*z4R[2]))}
earh3C		<- function(m,s,q,t,x,w,y,Ll,Lu,Ul,Uu){if(m>=mpfrthr){pbn<-m*mpfrpbn;s<-mpfr(s,precBits=pbn);q<-mpfr(q,precBits=pbn);t<-mpfr(t,precBits=pbn);x<-mpfr(x,precBits=pbn);w<-mpfr(w,precBits=pbn);y<-mpfr(y,precBits=pbn);Ll<-mpfr(Ll,precBits=pbn);Lu<-mpfr(Lu,precBits=pbn);Ul<-mpfr(Ul,precBits=pbn);Uu<-mpfr(Uu,precBits=pbn)}; z1L<-eaze4C(m,s,q,x,w,Ll,Uu); z1R<-eaze4C(m,q,t,w,y,Ll,Uu); z2L<-eaze4C(m,s,q,x,w,Lu,Uu); z2R<-eaze4C(m,q,t,w,y,Lu,Uu); z3L<-eaze4C(m,s,q,x,w,Ll,Ul); z3R<-eaze4C(m,q,t,w,y,Ll,Ul); z4L<-eaze4C(m,s,q,x,w,Lu,Ul); z4R<-eaze4C(m,q,t,w,y,Lu,Ul); c(s1=as.numeric(-z1L[2]-z1R[2]+z1L[1]*z1R[1]+z2L[1]+z2R[1]-z2L[2]*z2R[2]+z3L[1]+z3R[1]-z3L[2]*z3R[2]-z4L[2]-z4R[2]+z4L[1]*z4R[1]),s2=as.numeric(-z1L[3]-z1R[3]+z1L[2]*z1R[2]+z2L[2]+z2R[2]-z2L[3]*z2R[3]+z3L[2]+z3R[2]-z3L[3]*z3R[3]-z4L[3]-z4R[3]+z4L[2]*z4R[2]),s3=as.numeric(-z1L[4]-z1R[4]+z1L[3]*z1R[3]+z2L[3]+z2R[3]-z2L[4]*z2R[4]+z3L[3]+z3R[3]-z3L[4]*z3R[4]-z4L[4]-z4R[4]+z4L[3]*z4R[3]))}

### 2 - Biased Brownian Motion & Brownian Motion Functions
eabbm		<- function(x,s,t,bbmsim,bbmtar,bbmdom){ind<-0;while(ind==0){draw<-bbmsim(x,s,t,1);if(runif(1,0,1)<=bbmtar(draw,x,s,t)/bbmdom(draw,x,s,t)){ind <- 1}};list(draw=draw)} # Biased Brownian Motion
eabrbr		<- function(q,s,t,x,y){list(w=rnorm(1,x+(q-s)*(y-x)/(t-s),sd=((t-q)*(q-s)/(t-s))^(1/2)))} # Brownian Bridge
eabrbrm		<- function(q,s,t,x,y){qs<-sort(q);qout<-numeric(length(qs));qout[1]<-eabrbr(qs[1],s,t,x,y)$w;for(i in 2:length(qs)){qout[i]<-eabrbr(qs[i],qs[i-1],t,qout[i-1],y)$w};list(w=qout)} # Brownian Bridge Mutliple

### 3 - EA3 Related Functions
#### 3.1 - EA3 Single Layer
easl		<- function(s,t,x,y){
  xb <- min(x,y); yb<-max(x,y); adf<-sqrt((t-s)/2); act<-1
  gind <- 0; u1<-runif(1,0,1); m1<-3; il<-eagammaC(m1,s,t,x,y,xb-adf*act,yb+adf*act)
  while(gind==0){if(u1>=il[2]){act<-act+1;m1<-3;il<-eagammaC(m1,s,t,x,y,xb-adf*act,yb+adf*act)}else{if(u1<=il[1]){gind<-1}else{m1<-m1+2;il<-eagammaC(m1,s,t,x,y,xb-adf*act,yb+adf*act)}}}
  list(xb=xb,yb=yb,adf=adf,act=act,al=(act-1)*adf,au=act*adf)}
#### 3.2 - EA Min/Max
eaminmax	<- function(s,t,x,y,Ll,Lu,Ul,Uu){
  if(rbinom(1,1,0.5)==1){minI<-1}else{minI<--1; x<--x;y<--y;Ll<--Uu;Lu<--Ul}
  u1<-runif(1,exp(-2*(Ll-x)*(Ll-y)/(t-s)),exp(-2*(Lu-x)*(Lu-y)/(t-s))); u2<-runif(1,0,1); m<-x-0.5*(sqrt((y-x)^2-2*(t-s)*log(u1))-(y-x))
  c1 <- (y-m)^2/(2*(t-s)); c2 <- (m-x)^2/(2*(t-s)); I1<-rinvgauss(1,sqrt(c1/c2),2*c1);I2<-1/rinvgauss(1,sqrt(c2/c1),2*c2);V<-if(runif(1,0,1)<=(1/(1+sqrt(c1/c2)))){I1}else{I2}; tau <- (s*V+t)/(1+V)
  if(minI==-1){m<--m};list(m=m,tau=tau,minI=minI)}
#### 3.3 - Bessel Mid Points
eabesmid	<- function(q,s,tau,t,x,m,y,minI){
  if(minI==-1){x<--x;y<--y;m<--m} # Reflection
  bI<-0;if(q==s){bI<-1;w<-x};if(q==tau){bI<-1;w<-m};if(q==t){bI<-1;w<-y} # Boundary
  t<-t-s;tau<-tau-s;q<-q-s # Rescale time
  if(bI==0){if(q<tau){Ra1<-sqrt(tau);Ra2<-(x-m)*(tau-q)/((tau)^(3/2));Ra3<-(tau-q)/(tau)}else{Ra1<-sqrt(t-tau);Ra2<-(y-m)*(q-tau)/((t-tau)^(3/2));Ra3<-(q-tau)/(t-tau)};BB3<-rnorm(3,0,sqrt(Ra3*(1-Ra3)));w<-m+Ra1*sqrt((Ra2+BB3[1])^2+(BB3[2])^2+(BB3[3])^2)}
  list(w=minI*w)}
eabesmidm	<- function(Q,s,tau,t,x,m,y,minI){ ## Note - ensure Q includes Xs=x and Xy =y (omitted)
  wVn<-length(Q);wV<-numeric(wVn);Qlhs<-sort(Q[Q<=tau],decreasing=FALSE);Qrhs<-sort(Q[Q>tau],decreasing=TRUE)
  wV[1]<-eabesmid(Qlhs[1],s,tau,t,x,m,y,minI)$w;wV[wVn]<-eabesmid(Qrhs[1],s,tau,t,x,m,y,minI)$w
  if(length(Qlhs)>1){for(i in 2:length(Qlhs)){wV[i]<-eabesmid(Qlhs[i],Qlhs[i-1],tau,t,wV[i-1],m,y,minI)$w}}
  if(length(Qrhs)>1){for(i in 2:length(Qrhs)){wV[wVn-i+1]<-eabesmid(Qrhs[i],s,tau,Qrhs[i-1],x,m,wV[wVn-i+2],minI)$w}}
  list(wV=wV)}
#### 3.4 - Bessel Exceedance Evaluation
eabesex		<- function(sV,tV,xV,yV,m,B1,B2,minI){ # Vectors of equal length
  if(minI==-1){xV<--xV;yV<--yV;m<--m;B1<--B1;B2<--B2}
  u <- runif(1,0,1); mt<-3;  em<-matrix(0,length(sV),8); em[,1]<-sV; em[,2]<-tV; em[,3]<-xV; em[,4]<-yV
  B1evI<-B2evI<-0; while(B1evI==0){for(i in 1:dim(em)[1]){em[i,5:6]<-eadelC(mt,em[i,1],em[i,2],em[i,3],em[i,4],m,B1)};if(u<=prod(em[,5])){B1evI<-B2evI<-1;con1I<-1;con2I<-1;ex1I<-0;ex2I<-0}else{if(u>prod(em[,6])){B1evI<-1;con1I<-0;ex1I<-1}else{B1evI<-0;con1I<-0;ex1I<-0;mt<-mt+2}}}
  while(B2evI==0){for(i in 1:dim(em)[1]){em[i,7:8]<-eadelC(mt,em[i,1],em[i,2],em[i,3],em[i,4],m,B2)};if(u<=prod(em[,7])){B2evI<-1;con2I<-1;ex1I<-0}else{if(u>prod(em[,8])){B2evI<-1;con2I<-0;ex2I<-1}else{B2evI<-0;con2I<-0;ex2I<-0;mt<-mt+2}}}
  if(minI==-1){em[,3]<--em[,3];em[,4]<--em[,4]}; accI<-0; if(con1I==1){accI<-1}else{if(con2I==1){if(rbinom(1,1,0.5)==1){accI<-1}}}
  list(accI=accI,u=u,con1I=con1I,con2I=con2I,em=em)}

### 4 - AUEA Related Functions
#### 4.1 - Double Layer
eadl		<- function(s,t,x,y){
  # Initialise
  xb <- min(x,y); yb<-max(x,y); adf<-((t-s)/2)^(1/2); act<-1; Ll<-xb-adf; Lu<-xb; Ul<-yb; Uu<-yb+adf; Di<-0
  # Gamma initial single layer simulation
  gind <- 0; u1<-runif(1,0,1); m1<-3; il<-eagammaC(m1,s,t,x,y,Ll,Uu); while(gind==0){if(u1>=il[2]){act<-act+1;Lu<-Ll;Ll<-Ll-adf;Ul<-Uu;Uu<-Uu+adf;m1<-3;il<-eagammaC(m1,s,t,x,y,Ll,Uu)}else{if(u1<=il[1]){gind<-1}else{m1<-m1+2;il<-eagammaC(m1,s,t,x,y,Ll,Uu)}}}
  # Additional layer simulation
  if(act==1){Di<-1}else{bbind<-deind<-0; m2<-m3<-3; u2<-runif(1,0,1)
  while(deind==0){debd<-eabe3C(m2,s,t,x,y,Ll,xb,yb,Uu)[2:3];if(debd[2]<=0){m2<-m2+2}else{deind<-1}}
  while(bbind==0){nubd<-eabetaC(m3,s,t,x,y,Ll,Lu,Ul,Uu);if(u2<=(nubd[1]/debd[1])){bbind<-1;Di<-1}else{if(u2>=(nubd[2]/debd[2])){bbind<-1;u3<-runif(1,0,1);if(u3<1/2){Uu<-Ul;Ul<-yb;Di<-2}else{Ll<-Lu;Lu<-xb;Di<-3}}else{m2<-m2+2;m3<-m3+2;debd<-eabe3C(m2,s,t,x,y,Ll,xb,yb,Uu)[2:3]}}}}
  list(layer=c(s,t,x,y,Ll,Lu,Ul,Uu))}
#### 4.2 - Simulation of Double Layer Intermediary Points (Simulation of normal mixture approximation matrix where n=1 - primary function)
##### 4.2.1 - Primary Functions
earplays 	<- function(mwt,mmu){ # Nett matrix type 1
  nmat<-matrix(0,0,2); dmat<-cbind(mmu,mwt); colnames(dmat)<-c("mean","weight"); colnames(nmat)<-c("weight","mean")
  while(dim(dmat)[1]>0){
    eqmat<-rbind(1:dim(dmat)[1],mapply(function(ub,b){if(isTRUE(all.equal(ub,b))==TRUE){1}else{0}},rep(dmat[1,"mean"],dim(dmat)[1]),dmat[,"mean"])); commat<-dmat[eqmat[1,eqmat[2,]==1],,drop=FALSE]; dmat<-dmat[setdiff(1:dim(dmat)[1],eqmat[1,eqmat[2,]==1]),,drop=FALSE]
    wnew<-sum(commat[,"weight"]);if(wnew!=0){nmat<-rbind(nmat,c(wnew,commat[1,"mean"]))}};nmat}
##### 4.2.2 - Regular

earpreg	<- function(s,q,t,x,y,Ll,Lu,Ul,Uu){
  bbm <- x + (q-s)*(y-x)/(t-s); bbsd <- ((t-q)*(q-s)/(t-s))^(0.5);P1<--2/(q-s);P2<--2/(t-q);Alu<-Uu-Ll;Auu<-Uu-Lu;All<-Ul-Ll;Aul<-Ul-Lu
  mat <- matrix(c(-1,easiga(P1,x,Alu,Ll),easigb(P1,x,Alu,Ll),-1,easiga(P1,-x,Alu,-Uu),-easigb(P1,-x,Alu,-Uu),1,eaphia(P1,x,Alu,-1),eaphib(P1,x,Alu,-1),1,eaphia(P1,-x,Alu,-1),-eaphib(P1,-x,Alu,-1),-1,easiga(P2,y,Alu,Ll),easigb(P2,y,Alu,Ll),-1,easiga(P2,-y,Alu,-Uu),-easigb(P2,-y,Alu,-Uu),1,eaphia(P2,y,Alu,1),eaphib(P2,y,Alu,1),1,eaphia(P2,-y,Alu,1),-eaphib(P2,-y,Alu,1),1,easiga(P1,x,Alu,Ll)+easiga(P2,y,Alu,Ll),easigb(P1,x,Alu,Ll)+easigb(P2,y,Alu,Ll),1,easiga(P1,x,Alu,Ll)+easiga(P2,-y,Alu,-Uu),easigb(P1,x,Alu,Ll)-easigb(P2,-y,Alu,-Uu),1,easiga(P1,-x,Alu,-Uu)+easiga(P2,y,Alu,Ll),-easigb(P1,-x,Alu,-Uu)+easigb(P2,y,Alu,Ll),1,easiga(P1,-x,Alu,-Uu)+easiga(P2,-y,Alu,-Uu),-easigb(P1,-x,Alu,-Uu)-easigb(P2,-y,Alu,-Uu),1,easiga(P1,x,Auu,Lu),easigb(P1,x,Auu,Lu),1,easiga(P1,-x,Auu,-Uu),-easigb(P1,-x,Auu,-Uu),1,easiga(P2,y,Auu,Lu),easigb(P2,y,Auu,Lu),1,easiga(P2,-y,Auu,-Uu),-easigb(P2,-y,Auu,-Uu),-1,easiga(P1,x,Auu,Lu)+easiga(P2,y,Auu,Lu),easigb(P1,x,Auu,Lu)+easigb(P2,y,Auu,Lu),-1,easiga(P1,x,Auu,Lu)+easiga(P2,-y,Auu,-Uu),easigb(P1,x,Auu,Lu)-easigb(P2,-y,Auu,-Uu),1,easiga(P1,x,Auu,Lu)+eaphia(P2,y,Auu,1),easigb(P1,x,Auu,Lu)+eaphib(P2,y,Auu,1),1,easiga(P1,x,Auu,Lu)+eaphia(P2,-y,Auu,1),easigb(P1,x,Auu,Lu)-eaphib(P2,-y,Auu,1),-1,easiga(P1,-x,Auu,-Uu)+easiga(P2,y,Auu,Lu),-easigb(P1,-x,Auu,-Uu)+easigb(P2,y,Auu,Lu),-1,easiga(P1,-x,Auu,-Uu)+easiga(P2,-y,Auu,-Uu),-easigb(P1,-x,Auu,-Uu)-easigb(P2,-y,Auu,-Uu),1,easiga(P1,-x,Auu,-Uu)+eaphia(P2,y,Auu,1),-easigb(P1,-x,Auu,-Uu)+eaphib(P2,y,Auu,1),1,easiga(P1,-x,Auu,-Uu)+eaphia(P2,-y,Auu,1),-easigb(P1,-x,Auu,-Uu)-eaphib(P2,-y,Auu,1),1,eaphia(P1,x,Auu,-1)+easiga(P2,y,Auu,Lu),eaphib(P1,x,Auu,-1)+easigb(P2,y,Auu,Lu),1,eaphia(P1,x,Auu,-1)+easiga(P2,-y,Auu,-Uu),eaphib(P1,x,Auu,-1)-easigb(P2,-y,Auu,-Uu),-1,eaphia(P1,x,Auu,-1)+eaphia(P2,y,Auu,1),eaphib(P1,x,Auu,-1)+eaphib(P2,y,Auu,1),-1,eaphia(P1,x,Auu,-1)+eaphia(P2,-y,Auu,1),eaphib(P1,x,Auu,-1)-eaphib(P2,-y,Auu,1),1,eaphia(P1,-x,Auu,-1)+easiga(P2,y,Auu,Lu),-eaphib(P1,-x,Auu,-1)+easigb(P2,y,Auu,Lu),1,eaphia(P1,-x,Auu,-1)+easiga(P2,-y,Auu,-Uu),-eaphib(P1,-x,Auu,-1)-easigb(P2,-y,Auu,-Uu),-1,eaphia(P1,-x,Auu,-1)+eaphia(P2,y,Auu,1),-eaphib(P1,-x,Auu,-1)+eaphib(P2,y,Auu,1),-1,eaphia(P1,-x,Auu,-1)+eaphia(P2,-y,Auu,1),
                  -eaphib(P1,-x,Auu,-1)-eaphib(P2,-y,Auu,1),1,easiga(P1,x,All,Ll),easigb(P1,x,All,Ll),1,easiga(P1,-x,All,-Ul),-easigb(P1,-x,All,-Ul),1,easiga(P2,y,All,Ll),easigb(P2,y,All,Ll),1,easiga(P2,-y,All,-Ul),-easigb(P2,-y,All,-Ul),-1,easiga(P1,x,All,Ll)+easiga(P2,y,All,Ll),easigb(P1,x,All,Ll)+easigb(P2,y,All,Ll),-1,easiga(P1,x,All,Ll)+easiga(P2,-y,All,-Ul),easigb(P1,x,All,Ll)-easigb(P2,-y,All,-Ul),1,easiga(P1,x,All,Ll)+eaphia(P2,y,All,1),easigb(P1,x,All,Ll)+eaphib(P2,y,All,1),1,easiga(P1,x,All,Ll)+eaphia(P2,-y,All,1),easigb(P1,x,All,Ll)-eaphib(P2,-y,All,1),-1,easiga(P1,-x,All,-Ul)+easiga(P2,y,All,Ll),-easigb(P1,-x,All,-Ul)+easigb(P2,y,All,Ll),-1,easiga(P1,-x,All,-Ul)+easiga(P2,-y,All,-Ul),-easigb(P1,-x,All,-Ul)-easigb(P2,-y,All,-Ul),1,easiga(P1,-x,All,-Ul)+eaphia(P2,y,All,1),-easigb(P1,-x,All,-Ul)+eaphib(P2,y,All,1),1,easiga(P1,-x,All,-Ul)+eaphia(P2,-y,All,1),-easigb(P1,-x,All,-Ul)-eaphib(P2,-y,All,1),1,eaphia(P1,x,All,-1)+easiga(P2,y,All,Ll),eaphib(P1,x,All,-1)+easigb(P2,y,All,Ll),1,eaphia(P1,x,All,-1)+easiga(P2,-y,All,-Ul),eaphib(P1,x,All,-1)-easigb(P2,-y,All,-Ul),-1,eaphia(P1,x,All,-1)+eaphia(P2,y,All,1),eaphib(P1,x,All,-1)+eaphib(P2,y,All,1),-1,eaphia(P1,x,All,-1)+eaphia(P2,-y,All,1),eaphib(P1,x,All,-1)-eaphib(P2,-y,All,1),1,eaphia(P1,-x,All,-1)+easiga(P2,y,All,Ll),-eaphib(P1,-x,All,-1)+easigb(P2,y,All,Ll),1,eaphia(P1,-x,All,-1)+easiga(P2,-y,All,-Ul),-eaphib(P1,-x,All,-1)-easigb(P2,-y,All,-Ul),-1,eaphia(P1,-x,All,-1)+eaphia(P2,y,All,1),-eaphib(P1,-x,All,-1)+eaphib(P2,y,All,1),-1,eaphia(P1,-x,All,-1)+eaphia(P2,-y,All,1),-eaphib(P1,-x,All,-1)-eaphib(P2,-y,All,1),-1,easiga(P1,x,Aul,Lu),easigb(P1,x,Aul,Lu),-1,easiga(P1,-x,Aul,-Ul),-easigb(P1,-x,Aul,-Ul),1,eaphia(P1,x,Aul,-1),eaphib(P1,x,Aul,-1),1,eaphia(P1,-x,Aul,-1),-eaphib(P1,-x,Aul,-1),-1,easiga(P2,y,Aul,Lu),easigb(P2,y,Aul,Lu),-1,easiga(P2,-y,Aul,-Ul),-easigb(P2,-y,Aul,-Ul),1,eaphia(P2,y,Aul,1),eaphib(P2,y,Aul,1),1,eaphia(P2,-y,Aul,1),-eaphib(P2,-y,Aul,1),1,easiga(P1,x,Aul,Lu)+easiga(P2,y,Aul,Lu),easigb(P1,x,Aul,Lu)+easigb(P2,y,Aul,Lu),1,easiga(P1,x,Aul,Lu)+easiga(P2,-y,Aul,-Ul),easigb(P1,x,Aul,Lu)-easigb(P2,-y,Aul,-Ul),1,easiga(P1,-x,Aul,-Ul)+easiga(P2,y,Aul,Lu),-easigb(P1,-x,Aul,-Ul)+easigb(P2,y,Aul,Lu),1,easiga(P1,-x,Aul,-Ul)+easiga(P2,-y,Aul,-Ul),-easigb(P1,-x,Aul,-Ul)-easigb(P2,-y,Aul,-Ul)),64,3,byrow=TRUE)
  mwt <- mat[,1]*exp(0.5*(bbsd*mat[,3])^2+bbm*mat[,3]+mat[,2]); mwt <- mwt[c(1:12,33:52,1:64,1:32)]; mmu <- bbm+(bbsd^2)*mat[,3]; mmu<-mmu[c(1:12,33:52,1:64,1:32)]
  umat1 <- earplays(mwt[1:32],mmu[1:32]); umat2<-earplays(mwt[33:96],mmu[33:96]); umat3<-earplays(mwt[97:128],mmu[97:128])
  nmat <- rbind(cbind(L=Ll,U=Lu,B=1,umat1),cbind(L=Lu,U=Ul,B=2,umat2),cbind(L=Ul,U=Uu,B=3,umat3)); nmat<-nmat[nmat[,4]!=0,];
  eb <- (nmat[,5]-bbm)/(bbsd^2); ea<-log(abs(nmat[,4]))-0.5*(bbsd*eb)^2-bbm*eb;esi<-sign(nmat[,4])
  Lp <- pnorm(nmat[,1],nmat[,5],sd=bbsd); Up<-pnorm(nmat[,2],nmat[,5],sd=bbsd); itwt<-nmat[,4]*(Up-Lp)

  nmat <- cbind(esi,ea,eb,wt=nmat[,4],itwt,L=nmat[,1],Lp,U=nmat[,2],Up,mmu=nmat[,5],B=nmat[,3]);  pmat<-nmat[nmat[,1]==1,]; pmat[,5] <- pmat[,5]/sum(pmat[,5])
  dind1<-0; while(dind1==0){
    dind2<-0; while(dind2==0){
      sp  <- pmat[sample(1:(dim(pmat)[1]),1,replace=TRUE,prob=pmat[,5]),]; dr	<- qnorm(runif(1,sp[7],sp[9]),sp[10],bbsd); pmt <- pmat[pmat[,11]==sp[11], ,drop=FALSE]; nmt <- nmat[nmat[,11]==sp[11],]
      u <- runif(1,0,sum(pmt[,1]*exp(pmt[,2]+pmt[,3]*dr))); if(u<=sum(nmt[,1]*exp(nmt[,2]+nmt[,3]*dr))){dind2<-1}}
    dind3<-0; m<-3; while(dind3==0){cauc<-earh3C(m,s,q,t,x,dr,y,Ll,Lu,Ul,Uu);if(u<=cauc[2]){dind1<-dind3<-1}else{if(u>=cauc[3]){dind3<-1}else{m<-m+2}}}}
  list(draw=dr)
}
##### 4.2.3 - Hybrid (Bessel)
eahybbesex	<- function(sV,tV,xV,yV,m,B1,B2,minI){ # Vectors of equal length
  if(minI==-1){xV<--xV;yV<--yV;m<--m;B1<--B1;B2<--B2}
  u <- runif(1,0,1); mt<-3;  em<-matrix(0,length(sV),8); em[,1]<-sV; em[,2]<-tV; em[,3]<-xV; em[,4]<-yV
  B1evI<-B2evI<-0; while(B1evI==0){for(i in 1:dim(em)[1]){em[i,5:6]<-eadelC(mt,em[i,1],em[i,2],em[i,3],em[i,4],m,B1)};if(u<=prod(em[,5])){B1evI<-B2evI<-1;con1I<-1;con2I<-1;ex1I<-0;ex2I<-0}else{if(u>prod(em[,6])){B1evI<-1;con1I<-0;ex1I<-1}else{B1evI<-0;con1I<-0;ex1I<-0;mt<-mt+2}}}
  while(B2evI==0){for(i in 1:dim(em)[1]){em[i,7:8]<-eadelC(mt,em[i,1],em[i,2],em[i,3],em[i,4],m,B2)};if(u<=prod(em[,7])){B2evI<-1;con2I<-1;ex1I<-0}else{if(u>prod(em[,8])){B2evI<-1;con2I<-0;ex2I<-1}else{B2evI<-0;con2I<-0;ex2I<-0;mt<-mt+2}}}
  if(minI==-1){em[,3]<--em[,3];em[,4]<--em[,4]}; accI<-0; if(con1I==0){if(con2I==1){accI<-1}}
  list(accI=accI,u=u,con1I=con1I,con2I=con2I,em=em)}
earphyb	  	<- function(s,q,t,x,y,Ll,Lu,Ul,Uu){
  accI <-0; while(accI==0){
    mnmx <- eaminmax(s,t,x,y,Ll,Lu,Ul,Uu); m<-mnmx$m; tau<-mnmx$tau; minI<-mnmx$minI; if(minI==1){B1<-Ul;B2<-Uu}else{B1<--Lu;B2<--Ll} # Determine minimum and barriers
    dr <- eabesmid(q,s,tau,t,x,m,y,minI)$w # Make draw at q conditional on m
    skel <- matrix(0,4,2); skel[,1] <- c(s,t,tau,q); skel[,2] <- c(x,y,m,dr); skel <- matsort(skel,1)
    accI <- eahybbesex(skel[1:3,1],skel[2:4,1],skel[1:3,2],skel[2:4,2],m,B1,B2,minI)$accI}
  list(draw=dr)}
##### 4.2.4 - Mixture (Reg with Bessel cut out)
earpmix	  	<- function(s,q,t,x,y,Ll,Lu,Ul,Uu){
  cnt1I <- 0; cnt1M <- 250; cnt2M <- 25 # earpreg breakout count
  bbm <- x + (q-s)*(y-x)/(t-s); bbsd <- ((t-q)*(q-s)/(t-s))^(0.5);P1<--2/(q-s);P2<--2/(t-q);Alu<-Uu-Ll;Auu<-Uu-Lu;All<-Ul-Ll;Aul<-Ul-Lu
  mat <- matrix(c(-1,easiga(P1,x,Alu,Ll),easigb(P1,x,Alu,Ll),-1,easiga(P1,-x,Alu,-Uu),-easigb(P1,-x,Alu,-Uu),1,eaphia(P1,x,Alu,-1),eaphib(P1,x,Alu,-1),1,eaphia(P1,-x,Alu,-1),-eaphib(P1,-x,Alu,-1),-1,easiga(P2,y,Alu,Ll),easigb(P2,y,Alu,Ll),-1,easiga(P2,-y,Alu,-Uu),-easigb(P2,-y,Alu,-Uu),1,eaphia(P2,y,Alu,1),eaphib(P2,y,Alu,1),1,eaphia(P2,-y,Alu,1),-eaphib(P2,-y,Alu,1),1,easiga(P1,x,Alu,Ll)+easiga(P2,y,Alu,Ll),easigb(P1,x,Alu,Ll)+easigb(P2,y,Alu,Ll),1,easiga(P1,x,Alu,Ll)+easiga(P2,-y,Alu,-Uu),easigb(P1,x,Alu,Ll)-easigb(P2,-y,Alu,-Uu),1,easiga(P1,-x,Alu,-Uu)+easiga(P2,y,Alu,Ll),-easigb(P1,-x,Alu,-Uu)+easigb(P2,y,Alu,Ll),1,easiga(P1,-x,Alu,-Uu)+easiga(P2,-y,Alu,-Uu),-easigb(P1,-x,Alu,-Uu)-easigb(P2,-y,Alu,-Uu),1,easiga(P1,x,Auu,Lu),easigb(P1,x,Auu,Lu),1,easiga(P1,-x,Auu,-Uu),-easigb(P1,-x,Auu,-Uu),1,easiga(P2,y,Auu,Lu),easigb(P2,y,Auu,Lu),1,easiga(P2,-y,Auu,-Uu),-easigb(P2,-y,Auu,-Uu),-1,easiga(P1,x,Auu,Lu)+easiga(P2,y,Auu,Lu),easigb(P1,x,Auu,Lu)+easigb(P2,y,Auu,Lu),-1,easiga(P1,x,Auu,Lu)+easiga(P2,-y,Auu,-Uu),easigb(P1,x,Auu,Lu)-easigb(P2,-y,Auu,-Uu),1,easiga(P1,x,Auu,Lu)+eaphia(P2,y,Auu,1),easigb(P1,x,Auu,Lu)+eaphib(P2,y,Auu,1),1,easiga(P1,x,Auu,Lu)+eaphia(P2,-y,Auu,1),easigb(P1,x,Auu,Lu)-eaphib(P2,-y,Auu,1),-1,easiga(P1,-x,Auu,-Uu)+easiga(P2,y,Auu,Lu),-easigb(P1,-x,Auu,-Uu)+easigb(P2,y,Auu,Lu),-1,easiga(P1,-x,Auu,-Uu)+easiga(P2,-y,Auu,-Uu),-easigb(P1,-x,Auu,-Uu)-easigb(P2,-y,Auu,-Uu),1,easiga(P1,-x,Auu,-Uu)+eaphia(P2,y,Auu,1),-easigb(P1,-x,Auu,-Uu)+eaphib(P2,y,Auu,1),1,easiga(P1,-x,Auu,-Uu)+eaphia(P2,-y,Auu,1),-easigb(P1,-x,Auu,-Uu)-eaphib(P2,-y,Auu,1),1,eaphia(P1,x,Auu,-1)+easiga(P2,y,Auu,Lu),eaphib(P1,x,Auu,-1)+easigb(P2,y,Auu,Lu),1,eaphia(P1,x,Auu,-1)+easiga(P2,-y,Auu,-Uu),eaphib(P1,x,Auu,-1)-easigb(P2,-y,Auu,-Uu),-1,eaphia(P1,x,Auu,-1)+eaphia(P2,y,Auu,1),eaphib(P1,x,Auu,-1)+eaphib(P2,y,Auu,1),-1,eaphia(P1,x,Auu,-1)+eaphia(P2,-y,Auu,1),eaphib(P1,x,Auu,-1)-eaphib(P2,-y,Auu,1),1,eaphia(P1,-x,Auu,-1)+easiga(P2,y,Auu,Lu),-eaphib(P1,-x,Auu,-1)+easigb(P2,y,Auu,Lu),1,eaphia(P1,-x,Auu,-1)+easiga(P2,-y,Auu,-Uu),-eaphib(P1,-x,Auu,-1)-easigb(P2,-y,Auu,-Uu),-1,eaphia(P1,-x,Auu,-1)+eaphia(P2,y,Auu,1),-eaphib(P1,-x,Auu,-1)+eaphib(P2,y,Auu,1),-1,eaphia(P1,-x,Auu,-1)+eaphia(P2,-y,Auu,1),-eaphib(P1,-x,Auu,-1)-eaphib(P2,-y,Auu,1),1,
                  easiga(P1,x,All,Ll),easigb(P1,x,All,Ll),1,easiga(P1,-x,All,-Ul),-easigb(P1,-x,All,-Ul),1,easiga(P2,y,All,Ll),easigb(P2,y,All,Ll),1,easiga(P2,-y,All,-Ul),-easigb(P2,-y,All,-Ul),-1,easiga(P1,x,All,Ll)+easiga(P2,y,All,Ll),easigb(P1,x,All,Ll)+easigb(P2,y,All,Ll),-1,easiga(P1,x,All,Ll)+easiga(P2,-y,All,-Ul),easigb(P1,x,All,Ll)-easigb(P2,-y,All,-Ul),1,easiga(P1,x,All,Ll)+eaphia(P2,y,All,1),easigb(P1,x,All,Ll)+eaphib(P2,y,All,1),1,easiga(P1,x,All,Ll)+eaphia(P2,-y,All,1),easigb(P1,x,All,Ll)-eaphib(P2,-y,All,1),-1,easiga(P1,-x,All,-Ul)+easiga(P2,y,All,Ll),-easigb(P1,-x,All,-Ul)+easigb(P2,y,All,Ll),-1,easiga(P1,-x,All,-Ul)+easiga(P2,-y,All,-Ul),-easigb(P1,-x,All,-Ul)-easigb(P2,-y,All,-Ul),1,easiga(P1,-x,All,-Ul)+eaphia(P2,y,All,1),-easigb(P1,-x,All,-Ul)+eaphib(P2,y,All,1),1,easiga(P1,-x,All,-Ul)+eaphia(P2,-y,All,1),-easigb(P1,-x,All,-Ul)-eaphib(P2,-y,All,1),1,eaphia(P1,x,All,-1)+easiga(P2,y,All,Ll),eaphib(P1,x,All,-1)+easigb(P2,y,All,Ll),1,eaphia(P1,x,All,-1)+easiga(P2,-y,All,-Ul),eaphib(P1,x,All,-1)-easigb(P2,-y,All,-Ul),-1,eaphia(P1,x,All,-1)+eaphia(P2,y,All,1),eaphib(P1,x,All,-1)+eaphib(P2,y,All,1),-1,eaphia(P1,x,All,-1)+eaphia(P2,-y,All,1),eaphib(P1,x,All,-1)-eaphib(P2,-y,All,1),1,eaphia(P1,-x,All,-1)+easiga(P2,y,All,Ll),-eaphib(P1,-x,All,-1)+easigb(P2,y,All,Ll),1,eaphia(P1,-x,All,-1)+easiga(P2,-y,All,-Ul),-eaphib(P1,-x,All,-1)-easigb(P2,-y,All,-Ul),-1,eaphia(P1,-x,All,-1)+eaphia(P2,y,All,1),-eaphib(P1,-x,All,-1)+eaphib(P2,y,All,1),-1,eaphia(P1,-x,All,-1)+eaphia(P2,-y,All,1),-eaphib(P1,-x,All,-1)-eaphib(P2,-y,All,1),-1,easiga(P1,x,Aul,Lu),easigb(P1,x,Aul,Lu),-1,easiga(P1,-x,Aul,-Ul),-easigb(P1,-x,Aul,-Ul),1,eaphia(P1,x,Aul,-1),eaphib(P1,x,Aul,-1),1,eaphia(P1,-x,Aul,-1),-eaphib(P1,-x,Aul,-1),-1,easiga(P2,y,Aul,Lu),easigb(P2,y,Aul,Lu),-1,easiga(P2,-y,Aul,-Ul),-easigb(P2,-y,Aul,-Ul),1,eaphia(P2,y,Aul,1),eaphib(P2,y,Aul,1),1,eaphia(P2,-y,Aul,1),-eaphib(P2,-y,Aul,1),1,easiga(P1,x,Aul,Lu)+easiga(P2,y,Aul,Lu),easigb(P1,x,Aul,Lu)+easigb(P2,y,Aul,Lu),1,easiga(P1,x,Aul,Lu)+easiga(P2,-y,Aul,-Ul),easigb(P1,x,Aul,Lu)-easigb(P2,-y,Aul,-Ul),1,easiga(P1,-x,Aul,-Ul)+easiga(P2,y,Aul,Lu),-easigb(P1,-x,Aul,-Ul)+easigb(P2,y,Aul,Lu),1,easiga(P1,-x,Aul,-Ul)+easiga(P2,-y,Aul,-Ul),-easigb(P1,-x,Aul,-Ul)-easigb(P2,-y,Aul,-Ul)),64,3,byrow=TRUE)
  mwt <- mat[,1]*exp(0.5*(bbsd*mat[,3])^2+bbm*mat[,3]+mat[,2]); mwt <- mwt[c(1:12,33:52,1:64,1:32)]; mmu <- bbm+(bbsd^2)*mat[,3]; mmu<-mmu[c(1:12,33:52,1:64,1:32)]
  umat1 <- earplays(mwt[1:32],mmu[1:32]); umat2<-earplays(mwt[33:96],mmu[33:96]); umat3<-earplays(mwt[97:128],mmu[97:128])
  nmat <- rbind(cbind(L=Ll,U=Lu,B=1,umat1),cbind(L=Lu,U=Ul,B=2,umat2),cbind(L=Ul,U=Uu,B=3,umat3)); nmat <- nmat[nmat[,4]!=0,]
  eb <- (nmat[,5]-bbm)/(bbsd^2); ea<-log(abs(nmat[,4]))-0.5*(bbsd*eb)^2-bbm*eb;esi<-sign(nmat[,4])
  Lp <- pnorm(nmat[,1],nmat[,5],sd=bbsd); Up<-pnorm(nmat[,2],nmat[,5],sd=bbsd); itwt<-nmat[,4]*(Up-Lp)
  nmat <- cbind(esi,ea,eb,wt=nmat[,4],itwt=(itwt)/sum(itwt),L=nmat[,1],Lp,U=nmat[,2],Up,mmu=nmat[,5],B=nmat[,3]); pmat<-nmat[nmat[,1]==1,]
  dind1<-0; while(dind1==0){
    dind2<-0;cnt2I<-0;while(dind2==0){
      cnt2I<-cnt2I+1
      sp  <- pmat[sample(1:(dim(pmat)[1]),1,replace=TRUE,prob=pmat[,5]),]; dr	<- qnorm(runif(1,sp[7],sp[9]),sp[10],bbsd); pmt <- pmat[pmat[,11]==sp[11],]; nmt <- nmat[nmat[,11]==sp[11],]
      u <- runif(1,0,sum(pmt[,1]*exp(pmt[,2]+pmt[,3]*dr))); if(u<=sum(nmt[,1]*exp(nmt[,2]+nmt[,3]*dr))){dind2<-1}
      if(dind2==0){if(cnt2I>=cnt2M){dind2<--1}}}
    if(dind2==1){dind3<-0; m<-3; while(dind3==0){cauc<-earh3C(m,s,q,t,x,dr,y,Ll,Lu,Ul,Uu);if(u<=cauc[2]){dind1<-dind3<-1}else{if(u>=cauc[3]){dind3<-1}else{m<-m+2}}}}
    if(dind1==0){cnt1I<-cnt1I+1; if(cnt1I>=cnt1M){dind1<-1;dr<-earphyb(s,q,t,x,y,Ll,Lu,Ul,Uu)$draw}}} # Index breakout counter and breakout with draw
  list(draw=dr)}
##### 4.2.5 - Brute Force Rejection Sampler (Lipschitz Approach)
brutepbn<-200
brutesig	<- function(j,s,t,x,y,L,U,side){ # Calculate sigma
  pbn<-brutepbn;s<-mpfr(s,precBits=pbn);t<-mpfr(t,precBits=pbn);x<-mpfr(x,precBits=pbn);y<-mpfr(y,precBits=pbn);L<-mpfr(L,precBits=pbn);U<-mpfr(U,precBits=pbn)
  if(side==-1){a<-(-2/(t-s))*((abs(U-L)*j)^2+2*abs(U-L)*j*min(L,U)+(min(L,U))^2-abs(U-L)*j*x-min(L,U)*x); b<-(-2/(t-s))*(-abs(U-L)*j-min(L,U)+x)}
  if(side==1){a<-(-2/(t-s))*((abs(U-L)*j)^2+2*abs(U-L)*j*min(L,U)+(min(L,U))^2-abs(U-L)*j*y-min(L,U)*y);b<-(-2/(t-s))*(-abs(U-L)*j-min(L,U)+y)}
  list(a=as.numeric(a),b=as.numeric(b))}
brutetau	<- function(j,s,t,x,y,L,U,side){ # Calculate tau
  pbn<-brutepbn;s<-mpfr(s,precBits=pbn);t<-mpfr(t,precBits=pbn);x<-mpfr(x,precBits=pbn);y<-mpfr(y,precBits=pbn);L<-mpfr(L,precBits=pbn);U<-mpfr(U,precBits=pbn)
  if(side==-1){a<-(-2*j/(t-s))*(j*abs(U-L)^2+abs(U-L)*x);b<-(-2*j/(t-s))*(-abs(U-L))}
  if(side==1){a<-(-2*j/(t-s))*(j*abs(U-L)^2-abs(U-L)*y);b<-(-2*j/(t-s))*(abs(U-L))}
  list(a=as.numeric(a),b=as.numeric(b))}
brutezet	<- function(n,s,t,x,y,L,U,side){ # Calculate zeta
  zetm<-matrix(0,0,3); colnames(zetm)<-c("sign","a","b")
  if(even(n)){for(i in 1:(n/2)){spos<-brutesig(i,s,t,x,y,L,U,side); sneg<-brutesig(i,s,t,-x,-y,-L,-U,side); tpos<-brutetau(i,s,t,x,y,L,U,side); tneg<-brutetau(i,s,t,-x,-y,-L,-U,side); zetm<-rbind(zetm,c(1,spos$a,spos$b),c(1,sneg$a,-sneg$b),c(-1,tpos$a,tpos$b),c(-1,tneg$a,-tneg$b))}}
  if(odd(n)){if(n>1){for(i in 1:((n-1)/2)){spos<-brutesig(i,s,t,x,y,L,U,side); sneg<-brutesig(i,s,t,-x,-y,-L,-U,side); tpos<-brutetau(i,s,t,x,y,L,U,side); tneg<-brutetau(i,s,t,-x,-y,-L,-U,side); zetm<-rbind(zetm,c(1,spos$a,spos$b),c(1,sneg$a,-sneg$b),c(-1,tpos$a,tpos$b),c(-1,tneg$a,-tneg$b))}};spos<-brutesig((n+1)/2,s,t,x,y,L,U,side); sneg<-brutesig((n+1)/2,s,t,-x,-y,-L,-U,side); zetm<-rbind(zetm,c(1,spos$a,spos$b),c(1,sneg$a,-sneg$b))}
  list(zetm=zetm)}
brutezetmut <- function(n,s,q,t,x,y,L,U){w <- 1 # Dummy w	# Calculate zeta multiples
zetm<-matrix(0,0,3); colnames(zetm)<-c("sign","a","b"); zetmL<-brutezet(n,s,q,x,w,L,U,-1)$zetm; zetmR<-brutezet(n,q,t,w,y,L,U,1)$zetm; for(i in 1:(dim(zetmL)[1])){for(j in 1:(dim(zetmR)[1])){zetm<-rbind(zetm,c(zetmL[i,"sign"]*zetmR[j,"sign"],zetmL[i,"a"]+zetmR[j,"a"],zetmL[i,"b"]+zetmR[j,"b"]))}}
list(zetm=zetm)}
bruterho	<- function(n,s,q,t,x,y,Ll,Lu,Ul,Uu){w <- 1 # Dummy w	# Calculate rho
gamm1<-matrix(0,0,3); zetmL1<-brutezet(n+1,s,q,x,w,Ll,Uu,-1)$zetm; zetmL1[,"sign"]<--zetmL1[,"sign"]; zetmR1<-brutezet(n+1,q,t,w,y,Ll,Uu,1)$zetm; zetmR1[,"sign"]<--zetmR1[,"sign"]; zetmLR1<-brutezetmut(n,s,q,t,x,y,Ll,Uu)$zetm; gamm1<-rbind(gamm1,zetmL1,zetmR1,zetmLR1); gamm1<-cbind(gamm1,rep(Ll,dim(gamm1)[1]),rep(Uu,dim(gamm1)[1]))
gamm2<-matrix(0,0,3); zetmL2<-brutezet(n,s,q,x,w,Lu,Uu,-1)$zetm; zetmR2<-brutezet(n,q,t,w,y,Lu,Uu,1)$zetm; zetmLR2<-brutezetmut(n+1,s,q,t,x,y,Lu,Uu)$zetm; zetmLR2[,"sign"]<--zetmLR2[,"sign"]; gamm2<-rbind(gamm2,zetmL2,zetmR2,zetmLR2); gamm2<-cbind(gamm2,rep(Lu,dim(gamm2)[1]),rep(Uu,dim(gamm2)[1]))
gamm3<-matrix(0,0,3); zetmL3<-brutezet(n,s,q,x,w,Ll,Ul,-1)$zetm; zetmR3<-brutezet(n,q,t,w,y,Ll,Ul,1)$zetm; zetmLR3<-brutezetmut(n+1,s,q,t,x,y,Ll,Ul)$zetm; zetmLR3[,"sign"]<--zetmLR3[,"sign"]; gamm3<-rbind(gamm3,zetmL3,zetmR3,zetmLR3); gamm3<-cbind(gamm3,rep(Ll,dim(gamm3)[1]),rep(Ul,dim(gamm3)[1]))
gamm4<-matrix(0,0,3); zetmL4<-brutezet(n+1,s,q,x,w,Lu,Ul,-1)$zetm; zetmL4[,"sign"]<--zetmL4[,"sign"]; zetmR4<-brutezet(n+1,q,t,w,y,Lu,Ul,1)$zetm; zetmR4[,"sign"]<--zetmR4[,"sign"]; zetmLR4<-brutezetmut(n,s,q,t,x,y,Lu,Ul)$zetm; gamm4<-rbind(gamm4,zetmL4,zetmR4,zetmLR4); gamm4<-cbind(gamm4,rep(Lu,dim(gamm4)[1]),rep(Ul,dim(gamm4)[1]))
rhom<-matrix(0,0,5); rhom<-rbind(rhom,gamm1,gamm2,gamm3,gamm4); colnames(rhom)<-c("sign","a","b","L","U")
list(rhom=rhom,gamm1=gamm1,gamm2=gamm2,gamm3=gamm3,gamm4=gamm4)}
highevalab	<- function(mat,w){ # Calculate rho at point
  helmat <<- mat; helw <<- w
  rmat<-mat[mat[,5]>=w,,drop=FALSE]; rmat<-rmat[rmat[,4]<w,,drop=FALSE]
  helrmat <<- rmat
  if(dim(rmat)[1]==0){ev<-0}else{pbn<-ceiling(100+10*(log(max(abs(rmat[,2]+rmat[,3]*w)))/log(2))); helpbn <<- pbn; s<-mpfr(rmat[,1],pbn); a<-mpfr(rmat[,2],pbn);b<-mpfr(rmat[,3],pbn);w<-mpfr(w,pbn); ev<-sum(s*exp(a+b*w))}
  c(ev=as.numeric(ev))}
highevalab2	<- function(mat,w){ # Calculate bound for rho using only positive terms (backup if highevalab negative)
  helmat <<- mat; helw <<- w
  rmat<-mat[mat[,5]>=w,,drop=FALSE]; rmat<-rmat[rmat[,4]<w,,drop=FALSE]
  rmat <- rmat[rmat[,1]==1, ,drop=FALSE]
  if(dim(rmat)[1]==0){ev<-0}else{pbn<-ceiling(100+10*(log(max(abs(rmat[,2]+rmat[,3]*w)))/log(2))); helpbn <<- pbn; s<-mpfr(rmat[,1],pbn); a<-mpfr(rmat[,2],pbn);b<-mpfr(rmat[,3],pbn);w<-mpfr(w,pbn); ev<-sum(s*exp(a+b*w))}
  c(ev=as.numeric(ev))}
brutelip	<- function(nmu,bbsd,L,U){ # Calculate normal lipschitz constant in interval
  conI<-0; if(min(U-(nmu+bbsd),(nmu+bbsd)-L)>=0){conI<-1}; if(min(U-(nmu-bbsd),(nmu-bbsd)-L)>=0){conI<-1}
  if(conI==1){lip<-abs((1/bbsd)*exp(-1/2))}else{lip<-max(abs(((nmu-L)/(bbsd^2))*exp(-(1/(2*bbsd^2))*(L-nmu)^2)),abs(((nmu-U)/(bbsd^2))*exp(-(1/(2*bbsd^2))*(U-nmu)^2)))}
  c(lip=lip)}
brutenett1 	<- function(mat,bbm){ # Nett matrix type 1
  nmat<-matrix(0,0,dim(mat)[2]); dmat<-mat; colnames(nmat)<-colnames(dmat)<-colnames(mat)
  while(dim(dmat)[1]>0){
    eqmat<-rbind(1:dim(dmat)[1],mapply(function(ub,b){if(isTRUE(all.equal(ub,b))==TRUE){1}else{0}},rep(dmat[1,"b"],dim(dmat)[1]),dmat[,"b"])); commat<-dmat[eqmat[1,eqmat[2,]==1],,drop=FALSE]; dmat<-dmat[setdiff(1:dim(dmat)[1],eqmat[1,eqmat[2,]==1]),,drop=FALSE]
    pbn<-brutepbn; sm<-mpfr(commat[,"sign"],pbn); am<-mpfr(commat[,"a"],pbn); bm<-mpfr(commat[1,"b"],pbn); Lm<-mpfr(commat[1,"L"],pbn); Um<-mpfr(commat[1,"U"],pbn); bbmm<-mpfr(bbm,pbn); bbsdm<-mpfr(commat[1,"sd"],pbn);
    wnew<-sum(sm*exp(am+bm*bbmm+0.5*(bm*bbsdm)^2)); if(wnew!=0){
      snew<-sign(wnew); anew<-log(wnew*snew)-bm*bbmm-0.5*(bm*bbsdm)^2; trnew<-pnorm(Um,mean=bbmm+bm*bbsdm^2,sd=bbsdm)-pnorm(Lm,mean=bbmm+bm*bbsdm^2,sd=bbsdm); trwnew<-wnew*trnew; lipnew<-brutelip(bbmm+bm*bbsdm^2,bbsdm,Lm,Um)
      nmat<-rbind(nmat,c(snew,as.numeric(anew),commat[1,"b"],commat[1,"L"],commat[1,"U"],as.numeric(wnew),as.numeric(trnew),as.numeric(trwnew),commat[1,"mean"],commat[1,"sd"],as.numeric(lipnew)))}}
  list(nmat=nmat)}
brutenett2 	<- function(mat){ # Nett matrix type 2
  nmat<-matrix(0,0,dim(mat)[2]); dmat<-mat; colnames(nmat)<-colnames(dmat)<-colnames(mat)
  while(dim(dmat)[1]>0){
    eqmat<-rbind(1:dim(dmat)[1],mapply(function(ub,b){if(isTRUE(all.equal(ub,b))==TRUE){1}else{0}},rep(dmat[1,"b"],dim(dmat)[1]),dmat[,"b"])); commat<-dmat[eqmat[1,eqmat[2,]==1],,drop=FALSE]; dmat<-dmat[setdiff(1:dim(dmat)[1],eqmat[1,eqmat[2,]==1]),,drop=FALSE]
    a1<-sum(commat[,"a1"]);w1<-sum(commat[,"w1"]);t1<-sum(commat[,"t1"]);wt1<-sum(commat[,"wt1"]);l1<-sum(commat[,"l1"]);a2<-sum(commat[,"a2"]);w2<-sum(commat[,"w2"]);t2<-sum(commat[,"t2"]);wt2<-sum(commat[,"wt2"]);l2<-sum(commat[,"l2"]);a3<-sum(commat[,"a3"]);w3<-sum(commat[,"w3"]);t3<-sum(commat[,"t3"]);wt3<-sum(commat[,"wt3"]);l3<-sum(commat[,"l3"])
    wt<-wt1+wt2+wt3; lip<-max(l1,l2,l3)
    nmat<-rbind(nmat,c(commat[1,"b"],commat[1,"mean"],commat[1,"sd"],a1,w1,t1,wt1,l1,a2,w2,t2,wt2,l2,a3,w3,t3,wt3,l3,wt,lip))}
  list(nmat=nmat)}
brutenorm 	<- function(mat,bbm,bbsd){ # Calculate normalising costant, Lipshitz constant and normalized Lipshitz constant
  pbn<-brutepbn; sm<-mpfr(mat[,"sign"],pbn);am<-mpfr(mat[,"a"],pbn);bm<-mpfr(mat[,"b"],pbn);Lm<-mpfr(mat[,"L"],pbn);Um<-mpfr(mat[,"U"],pbn);trm<-mpfr(mat[,"trunc"],pbn);lipm<-mpfr(mat[,"lip"],pbn);bbmm<-mpfr(bbm,pbn);bbsdm<-mpfr(bbsd,pbn)
  wei<-sm*exp(am+bm*bbmm+0.5*(bm*bbsdm)^2); mat[,"wei"]<-as.numeric(wei); trwei<-wei*trm; mat[,"trwei"]<-as.numeric(trwei); wlip<-abs(wei)*lipm; mat[,"wlip"]<-as.numeric(wlip); Z<-sum(trwei); glip<-sum(wlip);
  nwei<-wei/Z; mat[,"nwei"]<-as.numeric(nwei); ntrwei<-nwei*trm; mat[,"ntrwei"]<-as.numeric(ntrwei); nwlip<-abs(nwei)*lipm; mat[,"nwlip"]<-as.numeric(nwlip); nglip<-sum(nwlip)
  list(mat=mat,Z=as.numeric(Z),glip=as.numeric(glip),nglip=as.numeric(nglip))}
bruterhopi	<- function(n,s,q,t,x,y,Ll,Lu,Ul,Uu){ # Calculate full density matrix, normalising constant and Lipschitz constant for given parameters
  vars_brp <<- c(n,s,q,t,x,y,Ll,Lu,Ul,Uu)
  # Required fields
  mat <- bruterho(n,s,q,t,x,y,Ll,Lu,Ul,Uu); rhom <- mat$rhom; gamm1 <- mat$gamm1; gamm2 <- mat$gamm2; gamm3 <- mat$gamm3; gamm4 <- mat$gamm4 # Rho matrix
  bbm	<- x + (q-s)*(y-x)/(t-s); bbsd <- ((t-q)*(q-s)/(t-s))^(0.5) # Mean and Variance
  # Matrix
  prmat <- matrix(0,(dim(rhom)[1]),11); colnames(prmat) <- c("sign","a","b","L","U","weight","trunc","trweight","mean","sd","lip")
  prmat[,"sign"] <- rhom[,"sign"]; prmat[,"a"] <- rhom[,"a"]; prmat[,"b"] <- rhom[,"b"]; prmat[,"L"] <- rhom[,"L"]; prmat[,"U"] <- rhom[,"U"]; prmat[,"sd"] <- bbsd
  prmat[,"mean"] <- bbm+prmat[,"b"]*bbsd^2; prmat[,"weight"] <- prmat[,"sign"]*exp(prmat[,"a"]+prmat[,"b"]*bbm+0.5*(prmat[,"b"]*bbsd)^2); prmat[,"trunc"] <- pnorm(prmat[,"U"],mean=prmat[,"mean"],sd=bbsd)-pnorm(prmat[,"L"],mean=prmat[,"mean"],sd=bbsd);
  prmat[,"trweight"] <- prmat[,"weight"]*prmat[,"trunc"]; for(i in 1:dim(prmat)[1]){prmat[i,"lip"]<-brutelip(prmat[i,"mean"],bbsd,prmat[i,"L"],prmat[i,"U"])}
  m1 <- dim(gamm1)[1]; m2 <- m1 + dim(gamm2)[1]; m3 <- m2 + dim(gamm3)[1]; m4 <- m3 + dim(gamm4)[1]; prmat1 <- prmat[1:m1,]; prmat2 <- prmat[(m1+1):m2,]; prmat3 <- prmat[(m2+1):m3,]; prmat4 <- prmat[(m3+1):m4,]
  # Band split matrix
  splmat1<-splmat2<-splmat3<-matrix(0,0,11); colnames(splmat1)<-colnames(splmat2)<-colnames(splmat3)<-colnames(prmat);
  splmat1<-rbind(splmat1,prmat1,prmat3); splmat1[,"L"]<-Ll; splmat1[,"U"]<-Lu; nmat1<-brutenett1(splmat1,bbm)$nmat
  splmat2<-rbind(splmat2,prmat1,prmat2,prmat3,prmat4); splmat2[,"L"]<-Lu; splmat2[,"U"]<-Ul; nmat2<-brutenett1(splmat2,bbm)$nmat
  splmat3<-rbind(splmat3,prmat1,prmat2); splmat3[,"L"]<-Ul; splmat3[,"U"]<-Uu; nmat3<-brutenett1(splmat3,bbm)$nmat
  # Net matrix
  nmata<-rbind(nmat1,nmat2,nmat3)
  nmatb<-matrix(0,0,20); colnames(nmatb)<-c("b","mean","sd","a1","w1","t1","wt1","l1","a2","w2","t2","wt2","l2","a3","w3","t3","wt3","l3","wt","l")
  nmatb<-rbind(nmatb,cbind(nmat1[,c(3,9,10),drop=FALSE],nmat1[,c(2,6,7,8,11),drop=FALSE],rep(0,dim(nmat1)[1]),rep(0,dim(nmat1)[1]),rep(0,dim(nmat1)[1]),rep(0,dim(nmat1)[1]),rep(0,dim(nmat1)[1]),rep(0,dim(nmat1)[1]),rep(0,dim(nmat1)[1]),rep(0,dim(nmat1)[1]),rep(0,dim(nmat1)[1]),rep(0,dim(nmat1)[1]),rep(0,dim(nmat1)[1]),rep(0,dim(nmat1)[1])))
  nmatb<-rbind(nmatb,cbind(nmat2[,c(3,9,10),drop=FALSE],rep(0,dim(nmat2)[1]),rep(0,dim(nmat2)[1]),rep(0,dim(nmat2)[1]),rep(0,dim(nmat2)[1]),rep(0,dim(nmat2)[1]),nmat2[,c(2,6,7,8,11),drop=FALSE],rep(0,dim(nmat2)[1]),rep(0,dim(nmat2)[1]),rep(0,dim(nmat2)[1]),rep(0,dim(nmat2)[1]),rep(0,dim(nmat2)[1]),rep(0,dim(nmat2)[1]),rep(0,dim(nmat2)[1])))
  nmatb<-rbind(nmatb,cbind(nmat3[,c(3,9,10),drop=FALSE],rep(0,dim(nmat3)[1]),rep(0,dim(nmat3)[1]),rep(0,dim(nmat3)[1]),rep(0,dim(nmat3)[1]),rep(0,dim(nmat3)[1]),rep(0,dim(nmat3)[1]),rep(0,dim(nmat3)[1]),rep(0,dim(nmat3)[1]),rep(0,dim(nmat3)[1]),rep(0,dim(nmat3)[1]),nmat3[,c(2,6,7,8,11),drop=FALSE],rep(0,dim(nmat3)[1]),rep(0,dim(nmat3)[1])))
  nmatb<-brutenett2(nmatb)
  # Normalised matrix
  normm<-cbind(nmata[,c(1:6,8,6,8,7,9,10,11,11,11)]); colnames(normm)<-c("sign","a","b","L","U","wei","trwei","nwei","ntrwei","trunc","mean","sd","lip","wlip","nwlip")
  normmfn<-brutenorm(normm,bbm,bbsd); normm<-normmfn$mat; Z<-normmfn$Z; glip<-normmfn$glip; nglip<-normmfn$nglip
  # Output
  list(bbm=bbm,bbsd=bbsd,rhom=rhom,prmat=prmat,nmata=nmata,nmatb=nmatb,normm=normm,Z=Z,gausslip=abs((1/sqrt(2*pi*bbsd^2))*(-1/bbsd)*exp(-1/2)),glip=glip,nglip=nglip,Ll=Ll,Lu=Lu,Ul=Ul,Uu=Uu)}
bruteneval	<- function(m,s,q,t,x,y,Ll,Lu,Ul,Uu){ # m odd only # Evaluate ratio of alternating Cauchy sequence normalising constants
  vars2 <<- c(m,s,q,t,x,y,Ll,Lu,Ul,Uu)
  upper <- bruterhopi(m,s,q,t,x,y,Ll,Lu,Ul,Uu); lower <- bruterhopi(m+1,s,q,t,x,y,Ll,Lu,Ul,Uu)
  ratio <- min(lower$Z/upper$Z,upper$Z/lower$Z)
  list(upper=upper,lower=lower,ratio=ratio)}
brutemesh	<- function(upper,depth=1){ # Caculate local Lipschitz constants and bounding uniforms
  bm_upper <<- upper
  normm<-upper$normm; bbm<-upper$bbm; bbsd<-upper$bbsd; Z<-upper$Z; Ll<-upper$Ll; Lu<-upper$Lu; Ul<-upper$Ul; Uu<-upper$Uu
  meshlen<-(ceiling((Uu-Ll)/bbsd)+50)*depth; f1<-(Lu-Ll)/(Uu-Ll); f2<-(Ul-Lu)/(Uu-Ll); f3<-(Uu-Ul)/(Uu-Ll); mesh<-c(seq(Ll,Lu,length=ceiling(meshlen*f1)+1),seq(Lu,Ul,length=ceiling(meshlen*f2)+1),seq(Ul,Uu,length=ceiling(meshlen*f3)+1))
  lipV<-numeric(length(mesh)-1); for(i in 1:(length(mesh)-1)){ # Calculate local lipschitz constants
    Up<-mesh[i+1]; Lp<-mesh[i]; normmC<-normm[normm[,"L"]<=Lp,,drop=FALSE]; normmC<-normmC[normmC[,"U"]>=Up,,drop=FALSE]; lipsV<-numeric(dim(normmC)[1]); for(j in 1:dim(normmC)[1]){lipsV[j]<-brutelip(normmC[j,"mean"],normmC[j,"sd"],Lp,Up)*abs(normmC[j,"wei"])};
    lipV[i]<-sum(lipsV)}
  densV<-numeric(length(mesh)); ev <- numeric(length(mesh)); for(i in 1:(length(mesh))){ev[i]<-highevalab(normm,mesh[i]); if(ev[i] < 0){ev[i]<-highevalab2(normm,mesh[i])}; densV[i]<-ev[i]*dnorm(mesh[i],mean=bbm,sd=bbsd)} ; bm_densV <<- densV# Calculate density at mesh points
  bndV<-numeric(length(mesh)-1); bndZ<-0; for(i in 1:length(bndV)){bndV[i]<-max(densV[i],densV[i+1])+((mesh[i+1]-mesh[i])/2)*lipV[i]; bndZ<-bndZ+bndV[i]*(mesh[i+1]-mesh[i])} # Calculate bounding uniforms
  aprob<-Z/bndZ; pmat<-matrix(0,(length(mesh)-1),5); pmat[,1]<-mesh[1:(length(mesh)-1)]; pmat[,2]<-mesh[2:length(mesh)]; pmat[,3]<-mesh[2:length(mesh)]-mesh[1:(length(mesh)-1)]; pmat[,4]<-bndV; pmat[,5]<-pmat[,3]*pmat[,4]; pmat[,5]<-pmat[,5]/sum(pmat[,5]); colnames(pmat)<-c("s","t","t-s","bnd","prob")
  list(pmat=pmat,aprob=aprob)}
earpbrute 	<- function(s,q,t,x,y,Ll,Lu,Ul,Uu){ # Bound piecewise constant uniforms
  vars <<- c(s,q,t,x,y,Ll,Lu,Ul,Uu)
  m<--1; R1<-0; while(R1<=0.01){m<-m+4; altmat<-bruteneval(m,s,q,t,x,y,Ll,Lu,Ul,Uu); upper<-altmat$upper; R1<-altmat$ratio} # R1 denotes the ratio of the lower to upper normalising constant
  mesh<-0; R2<-0; while(R2<=0.01){mesh<-mesh+1; bndmat<-brutemesh(upper,mesh); pmat<-bndmat$pmat; R2<-bndmat$aprob; R3<-R1*R2} # R2 denotes the ratio of the upper to bounding normalising constant, R3 is the lower bound on the acceptance probability
  dind1I<-0; while(dind1I==0){
    band<-pmat[sample(1:dim(pmat)[1],1,prob=pmat[,"prob"]),,drop=FALSE]; draw<-runif(1,band[,"s"],band[,"t"]); u<-runif(1,0,band[,"bnd"])
    nev<-30+m; dind2I<-0; while(dind2I==0){cauc<-earh3C(nev,s,q,t,x,draw,y,Ll,Lu,Ul,Uu)*dnorm(draw,mean=upper$bbm,sd=upper$bbsd); if(u<=cauc[2]){dind1I<-dind2I<-1}else{if(u>=cauc[3]){dind2I<-1}else{nev<-nev+20}}}}
  list(draw=draw,altmat=altmat,upper=upper,bndmat=bndmat,pmat=pmat,R1=R1,R2=R2,R3=R3)}
earpbrutemix <- function(s,q,t,x,y,Ll,Lu,Ul,Uu){
  brute_list <<- list(s=s,q=q,t=t,x=x,y=y,Ll=Ll,Lu=Lu,Ul=Ul,Uu=Uu)
  cnt1I <- 0; cnt1M <- 100000000; cnt2M <- 25 # earpreg breakout count
  if((t-s) < 10^(-9)*s){
	  s <- mpfr(s,900); q <- mpfr(q,900); t <- mpfr(t,900)
  }
  bbm <- as.numeric(x + (q-s)*(y-x)/(t-s)); bbsd <- as.numeric(((t-q)*(q-s)/(t-s))^(0.5));P1<-as.numeric(-2/(q-s));P2<-as.numeric(-2/(t-q));Alu<-Uu-Ll;Auu<-Uu-Lu;All<-Ul-Ll;Aul<-Ul-Lu
  derived_list <<- list(bbm=bbm, bbsd=bbsd,P1=P1,P2=P2, All=All,Alu=Alu,Aul=Aul,Auu=Auu)

  mat <- matrix(c(-1,easiga(P1,x,Alu,Ll),easigb(P1,x,Alu,Ll),-1,easiga(P1,-x,Alu,-Uu),-easigb(P1,-x,Alu,-Uu),1,eaphia(P1,x,Alu,-1),eaphib(P1,x,Alu,-1),1,eaphia(P1,-x,Alu,-1),
                  -eaphib(P1,-x,Alu,-1),-1,easiga(P2,y,Alu,Ll),easigb(P2,y,Alu,Ll),
                  -1,easiga(P2,-y,Alu,-Uu),-easigb(P2,-y,Alu,-Uu),1,eaphia(P2,y,Alu,1),
                  eaphib(P2,y,Alu,1),1,eaphia(P2,-y,Alu,1),-eaphib(P2,-y,Alu,1),1,
                  easiga(P1,x,Alu,Ll)+easiga(P2,y,Alu,Ll),
                  easigb(P1,x,Alu,Ll)+easigb(P2,y,Alu,Ll),1,
                  easiga(P1,x,Alu,Ll)+easiga(P2,-y,Alu,-Uu),
                  easigb(P1,x,Alu,Ll)-easigb(P2,-y,Alu,-Uu),1,
                  easiga(P1,-x,Alu,-Uu)+easiga(P2,y,Alu,Ll),
                  -easigb(P1,-x,Alu,-Uu)+easigb(P2,y,Alu,Ll),1,
                  easiga(P1,-x,Alu,-Uu)+easiga(P2,-y,Alu,-Uu),
                  -easigb(P1,-x,Alu,-Uu)-easigb(P2,-y,Alu,-Uu),1,
                  easiga(P1,x,Auu,Lu),easigb(P1,x,Auu,Lu),1,easiga(P1,-x,Auu,-Uu),
                  -easigb(P1,-x,Auu,-Uu),1,easiga(P2,y,Auu,Lu),easigb(P2,y,Auu,Lu),1,
                  easiga(P2,-y,Auu,-Uu),-easigb(P2,-y,Auu,-Uu),-1,
                  easiga(P1,x,Auu,Lu)+easiga(P2,y,Auu,Lu),
                  easigb(P1,x,Auu,Lu)+easigb(P2,y,Auu,Lu),-1,
                  easiga(P1,x,Auu,Lu)+easiga(P2,-y,Auu,-Uu),
                  easigb(P1,x,Auu,Lu)-easigb(P2,-y,Auu,-Uu),1,
                  easiga(P1,x,Auu,Lu)+eaphia(P2,y,Auu,1),
                  easigb(P1,x,Auu,Lu)+eaphib(P2,y,Auu,1),1,
                  easiga(P1,x,Auu,Lu)+eaphia(P2,-y,Auu,1),
                  easigb(P1,x,Auu,Lu)-eaphib(P2,-y,Auu,1),-1,
                  easiga(P1,-x,Auu,-Uu)+easiga(P2,y,Auu,Lu),
                  -easigb(P1,-x,Auu,-Uu)+easigb(P2,y,Auu,Lu),-1,
                  easiga(P1,-x,Auu,-Uu)+easiga(P2,-y,Auu,-Uu),
                  -easigb(P1,-x,Auu,-Uu)-easigb(P2,-y,Auu,-Uu),
                  1,easiga(P1,-x,Auu,-Uu)+eaphia(P2,y,Auu,1),
                  -easigb(P1,-x,Auu,-Uu)+eaphib(P2,y,Auu,1),1,
                  easiga(P1,-x,Auu,-Uu)+eaphia(P2,-y,Auu,1),
                  -easigb(P1,-x,Auu,-Uu)-eaphib(P2,-y,Auu,1),1,
                  eaphia(P1,x,Auu,-1)+easiga(P2,y,Auu,Lu),
                  eaphib(P1,x,Auu,-1)+easigb(P2,y,Auu,Lu),1,
                  eaphia(P1,x,Auu,-1)+easiga(P2,-y,Auu,-Uu),
                  eaphib(P1,x,Auu,-1)-easigb(P2,-y,Auu,-Uu),-1,eaphia(P1,x,Auu,-1)+eaphia(P2,y,Auu,1),eaphib(P1,x,Auu,-1)+eaphib(P2,y,Auu,1),-1,eaphia(P1,x,Auu,-1)+eaphia(P2,-y,Auu,1),eaphib(P1,x,Auu,-1)-eaphib(P2,-y,Auu,1),1,eaphia(P1,-x,Auu,-1)+easiga(P2,y,Auu,Lu),-eaphib(P1,-x,Auu,-1)+easigb(P2,y,Auu,Lu),1,eaphia(P1,-x,Auu,-1)+easiga(P2,-y,Auu,-Uu),-eaphib(P1,-x,Auu,-1)-easigb(P2,-y,Auu,-Uu),-1,eaphia(P1,-x,Auu,-1)+eaphia(P2,y,Auu,1),-eaphib(P1,-x,Auu,-1)+eaphib(P2,y,Auu,1),-1,eaphia(P1,-x,Auu,-1)+eaphia(P2,-y,Auu,1),-eaphib(P1,-x,Auu,-1)-eaphib(P2,-y,Auu,1),1,easiga(P1,x,All,Ll),easigb(P1,x,All,Ll),1,easiga(P1,-x,All,-Ul),-easigb(P1,-x,All,-Ul),1,easiga(P2,y,All,Ll),easigb(P2,y,All,Ll),1,easiga(P2,-y,All,-Ul),-easigb(P2,-y,All,-Ul),-1,easiga(P1,x,All,Ll)+easiga(P2,y,All,Ll),easigb(P1,x,All,Ll)+easigb(P2,y,All,Ll),-1,easiga(P1,x,All,Ll)+easiga(P2,-y,All,-Ul),easigb(P1,x,All,Ll)-easigb(P2,-y,All,-Ul),1,easiga(P1,x,All,Ll)+eaphia(P2,y,All,1),easigb(P1,x,All,Ll)+eaphib(P2,y,All,1),1,easiga(P1,x,All,Ll)+eaphia(P2,-y,All,1),easigb(P1,x,All,Ll)-eaphib(P2,-y,All,1),-1,easiga(P1,-x,All,-Ul)+easiga(P2,y,All,Ll),-easigb(P1,-x,All,-Ul)+easigb(P2,y,All,Ll),-1,easiga(P1,-x,All,-Ul)+easiga(P2,-y,All,-Ul),-easigb(P1,-x,All,-Ul)-easigb(P2,-y,All,-Ul),1,easiga(P1,-x,All,-Ul)+eaphia(P2,y,All,1),-easigb(P1,-x,All,-Ul)+eaphib(P2,y,All,1),1,easiga(P1,-x,All,-Ul)+eaphia(P2,-y,All,1),-easigb(P1,-x,All,-Ul)-eaphib(P2,-y,All,1),1,eaphia(P1,x,All,-1)+easiga(P2,y,All,Ll),eaphib(P1,x,All,-1)+easigb(P2,y,All,Ll),1,eaphia(P1,x,All,-1)+easiga(P2,-y,All,-Ul),eaphib(P1,x,All,-1)-easigb(P2,-y,All,-Ul),-1,eaphia(P1,x,All,-1)+eaphia(P2,y,All,1),eaphib(P1,x,All,-1)+eaphib(P2,y,All,1),-1,eaphia(P1,x,All,-1)+eaphia(P2,-y,All,1),eaphib(P1,x,All,-1)-eaphib(P2,-y,All,1),1,eaphia(P1,-x,All,-1)+easiga(P2,y,All,Ll),-eaphib(P1,-x,All,-1)+easigb(P2,y,All,Ll),1,eaphia(P1,-x,All,-1)+easiga(P2,-y,All,-Ul),-eaphib(P1,-x,All,-1)-easigb(P2,-y,All,-Ul),-1,eaphia(P1,-x,All,-1)+eaphia(P2,y,All,1),-eaphib(P1,-x,All,-1)+eaphib(P2,y,All,1),-1,eaphia(P1,-x,All,-1)+eaphia(P2,-y,All,1),-eaphib(P1,-x,All,-1)-eaphib(P2,-y,All,1),-1,easiga(P1,x,Aul,Lu),easigb(P1,x,Aul,Lu),-1,easiga(P1,-x,Aul,-Ul),-easigb(P1,-x,Aul,-Ul),1,eaphia(P1,x,Aul,-1),eaphib(P1,x,Aul,-1),1,eaphia(P1,-x,Aul,-1),-eaphib(P1,-x,Aul,-1),-1,easiga(P2,y,Aul,Lu),easigb(P2,y,Aul,Lu),-1,easiga(P2,-y,Aul,-Ul),-easigb(P2,-y,Aul,-Ul),1,eaphia(P2,y,Aul,1),eaphib(P2,y,Aul,1),1,eaphia(P2,-y,Aul,1),-eaphib(P2,-y,Aul,1),1,easiga(P1,x,Aul,Lu)+easiga(P2,y,Aul,Lu),easigb(P1,x,Aul,Lu)+easigb(P2,y,Aul,Lu),1,easiga(P1,x,Aul,Lu)+easiga(P2,-y,Aul,-Ul),easigb(P1,x,Aul,Lu)-easigb(P2,-y,Aul,-Ul),1,easiga(P1,-x,Aul,-Ul)+easiga(P2,y,Aul,Lu),-easigb(P1,-x,Aul,-Ul)+easigb(P2,y,Aul,Lu),1,easiga(P1,-x,Aul,-Ul)+easiga(P2,-y,Aul,-Ul),-easigb(P1,-x,Aul,-Ul)-easigb(P2,-y,Aul,-Ul)),64,3,byrow=TRUE)
  mwt <- mat[,1]*exp(0.5*(bbsd*mat[,3])^2+bbm*mat[,3]+mat[,2]); mwt <- mwt[c(1:12,33:52,1:64,1:32)]; mmu <- bbm+(bbsd^2)*mat[,3]; mmu<-mmu[c(1:12,33:52,1:64,1:32)]
  umat1 <- earplays(mwt[1:32],mmu[1:32]); umat2<-earplays(mwt[33:96],mmu[33:96]); umat3<-earplays(mwt[97:128],mmu[97:128])
  nmat <- rbind(cbind(L=Ll,U=Lu,B=1,umat1),cbind(L=Lu,U=Ul,B=2,umat2),cbind(L=Ul,U=Uu,B=3,umat3)); nmat <- nmat[nmat[,4]!=0,];
  eb <- (nmat[,5]-bbm)/(bbsd^2); ea<-log(abs(nmat[,4]))-0.5*(bbsd*eb)^2-bbm*eb;esi<-sign(nmat[,4])
  Lp <- pnorm(nmat[,1],nmat[,5],sd=bbsd); Up<-pnorm(nmat[,2],nmat[,5],sd=bbsd); itwt<-nmat[,4]*(Up-Lp);
  nmat <- cbind(esi,ea,eb,wt=nmat[,4],itwt,L=nmat[,1],Lp,U=nmat[,2],Up,mmu=nmat[,5],B=nmat[,3]); pmat<-nmat[nmat[,1]==1, ,drop=FALSE]; pmat[,5] <- pmat[,5]/sum(pmat[,5])
  #cut to high precision if necessart
  if(nrow(pmat)==0){
    print("cut to high prec 1")
    dr<-earpbrute(s,q,t,x,y,Ll,Lu,Ul,Uu)$draw
  } else if(any(is.nan(pmat[,5]))){
    print("cut to high prec 2")
    dr<-earpbrute(s,q,t,x,y,Ll,Lu,Ul,Uu)$draw
  } else if(any(pmat[,5] < 0)){
    print("cut to high prec 3")
    dr<-earpbrute(s,q,t,x,y,Ll,Lu,Ul,Uu)$draw
  } else {
    counter <- 0
  dind1<-0; while(dind1==0){
    counter <- counter + 1
    dind2<-0;cnt2I<-0;while(dind2==0){
      cnt2I<-cnt2I+1
      sp  <- pmat[sample(1:(dim(pmat)[1]),1,replace=TRUE,prob=pmat[,5]), ,drop=FALSE]; bbsd <<- bbsd; dr	<- qnorm(runif(1,sp[7],sp[9]),sp[10],bbsd); pmt <- pmat[pmat[,11]==sp[11], ,drop=FALSE]; nmt <- nmat[nmat[,11]==sp[11],]
      if(is.vector(pmt)==1){
        pmt <- matrix(pmt,nrow=1)
      }
      if(is.vector(nmt)==1){
        nmt <- matrix(nmt,nrow=1)
      }
      u <- runif(1,0,sum(pmt[,1]*exp(pmt[,2]+pmt[,3]*dr))); if(u<=sum(nmt[,1]*exp(nmt[,2]+nmt[,3]*dr))){dind2<-1}
      if(dind2==0){if(cnt2I>=cnt2M){dind2<--1}}}
    if(dind2==1){dind3<-0; m<-3; while(dind3==0){cauc<-earh3C(m,s,q,t,x,dr,y,Ll,Lu,Ul,Uu);if(u<=cauc[2]){dind1<-dind3<-1}else{if(u>=cauc[3]){dind3<-1}else{m<-m+2}}}}
    if(dind1==0){cnt1I<-cnt1I+1; if(cnt1I>=cnt1M){dind1<-1; dr<-earpbrute(s,q,t,x,y,Ll,Lu,Ul,Uu)$draw}}}
  }
  # Index breakout counter and breakout with draw
  list(draw=dr)}
##### 4.2.5 - Selection of default double layer intermediary point function
earp		<- earpbrutemix
#### 4.3 - Simulation of Bisected Layer
eabl		<- function(s,q,t,x,w,y,Ll,Lu,Ul,Uu){
  Lmin<-min(x,w);Lmax<-max(x,w);Rmin<-min(w,y);Rmax<-max(w,y);if(w>=Ul){Ul<-w};if(w<=Lu){Lu<-w}
  Lmat <- matrix(c(Ll,Lu,Ul,Uu,Ll,Lu,Ul,Uu,Lu,Lmin,Ul,Uu,Ll,Lu,Ul,Uu,Ll,Lu,Ul,Uu,Lu,Lmin,Ul,Uu,Ll,Lu,Lmax,Ul,Ll,Lu,Lmax,Ul,Lu,Lmin,Lmax,Ul),9,4,byrow=TRUE); Rmat<-matrix(c(Ll,Lu,Ul,Uu,Lu,Rmin,Ul,Uu,Ll,Lu,Ul,Uu,Ll,Lu,Rmax,Ul,Lu,Rmin,Rmax,Ul,Ll,Lu,Rmax,Ul,Ll,Lu,Ul,Uu,Lu,Rmin,Ul,Uu,Ll,Lu,Ul,Uu),9,4,byrow=TRUE); LRmat<-cbind(Lmat,Rmat)
  bbrind<-deind<-0;bbrct<-1;m1<-3;u1<-runif(1,0,1)
  while(bbrind==0){
    while(deind==0){debd<-earhoC(m1,s,q,t,x,w,y,Ll,Lu,Ul,Uu);if(debd[2]<=0){m1<-m1+2}else{deind<-1}}
    bebe<-matrix(0,bbrct,4); for(i in 1:bbrct){bebe[i,1:2]<-eabetaC(m1,s,q,x,w,LRmat[i,1],LRmat[i,2],LRmat[i,3],LRmat[i,4]);bebe[i,3:4]<-eabetaC(m1,q,t,w,y,LRmat[i,5],LRmat[i,6],LRmat[i,7],LRmat[i,8])}
    bd<-c(sum(bebe[,1]*bebe[,3])/debd[1],sum(bebe[,2]*bebe[,4])/debd[2])
    if(u1<=bd[1]){bbrind<-1}else{if(u1>=bd[2]){bbrct<-bbrct+1;deind<-0}else{m1<-m1+2;deind<-0}}}
  list(lyr2=matrix(c(c(s,q,x,w,LRmat[bbrct,1:4]),c(q,t,w,y,LRmat[bbrct,5:8])),2,8,byrow=TRUE))}
#### 4.4 - Simulation of Refined Layer
earl		<- function(s,t,x,y,Ll,Lu,Ul,Uu){
  if(max((Uu-Ul),(Lu-Ll))>((t-s)/2)^(0.5)){while(max((Uu-Ul),(Lu-Ll))>((t-s)/2)^(0.5)){
    Ls<-Ll+(Lu-Ll)/2; Us<-Uu-(Uu-Ul)/2; mat<-matrix(c(Ll,Ls,Us,Uu,Ls,Lu,Us,Uu,Ll,Ls,Ul,Us,Ls,Lu,Ul,Us),4,4,byrow=TRUE)
    bbind<-deind<-0;bbct<-1;m1<-3;u1<-runif(1,0,1); while(bbind==0){
      while(deind==0){debd<-eabe3C(m1,s,t,x,y,Ll,Lu,Ul,Uu)[2:3];if(debd[2]<=0){m1<-m1+2}else{deind<-1}}
      be<-matrix(0,bbct,2); for(i in 1:bbct){be[i,]<-eabetaC(m1,s,t,x,y,mat[i,1],mat[i,2],mat[i,3],mat[i,4])}
      bd<-c(sum(be[,1])/debd[1],sum(be[,2])/debd[2])
      if(u1<=bd[1]){bbind<-1}else{if(u1>=bd[2]){bbct<-bbct+1;deind<-0}else{m1<-m1+2;deind<-0}}}
    Ll<-mat[bbct,1];Lu<-mat[bbct,2];Ul<-mat[bbct,3];Uu<-mat[bbct,4]}}
  list(ref=c(s,t,x,y,Ll,Lu,Ul,Uu))}

