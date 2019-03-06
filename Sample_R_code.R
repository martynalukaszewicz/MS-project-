## calc_coverage() function returns coverage for a chi-square statistic for product multinomial along with plots of sorted distances 
## of Pearson's chi-sqyuared statistic for varying sample sizes over 10^4 iterations plotted against quantiles of chi-squared(r-1) distribution.

## Sample R code with small example of determining sample as described in master's thesis project.
## expected probabilities for 6 sampling times:
n.0 <- c(rep(50,q),rep(250,q),rep(750,q),rep(1000,q))   ## choose starting n value, same for each t
l <- length(n.0)

prob.t1 <- c(0.895, 0.1, 0.005)
prob.t2 <- c(0.22, 0.62, 0.15 ,0.01)
prob.t3 <- c( 0.01, 0.275 ,0.6, 0.1065,0.005,0.0035)
prob.t4 <- c(0.004, 0.146 ,0.85)
prob.t5 <- c(0.005,0.57,0.425)
prob.t6 <- c(0.001, 0.999)   ## violation of minimum of expected count of 0.0035 for expected probability

q <- 6  ## 6 sampling times 

prob.adj <- vector("list",q)

prob.adj[[1]] <- prob.t1
prob.adj[[2]] <- prob.t2
prob.adj[[3]] <- prob.t3
prob.adj[[4]] <- prob.t4
prob.adj[[5]] <- prob.t5
prob.adj[[6]] <- prob.t6

k <- length(unlist(prob.adj))-q  ## degrees of freedom of chi-squared

alpha <- 0.05 ## choose significance level alpha, 1-alpha=confidence level
nsim <- 1000 ## choose number of simulations

calc_coverage <- function(q,prob.adj,k,n.0,l,alpha,nsim){

  prob.0 <- rep(prob.adj,l/q)
  nsim <- 100## choose number of simulations
  N.emp <- matrix(0,nsim,l/q) 
  p.log.p <- vector("list",l)
  p.emp <- vector("list",l)  ## empirical proportions of simulations
  chi.sq.int <- matrix(0,nsim,l)
  chi.sq <- matrix(0,nsim,l/q)
  g.sq.int <- matrix(0,nsim,l)
  g.sq <- matrix(0,nsim,l/q)
  
  
  for(jj in 1:nsim){
    for (i in 1:l){
      p.emp[[i]] <- t(rmultinom(nsim,size=n.0[i],prob=unlist(prob.0[[i]]))/n.0[i])
      chi.sq.int[jj,i] <- n.0[i]*sum((p.emp[[i]][jj,]-prob.0[[i]])^2/prob.0[[i]])
      
      for (m in 1:(l/q) ){
        chi.sq[,m] <- rowSums(chi.sq.int[,((m-1)*q+1):(m*q)])
        N.emp[jj,m] <- n.0[m*q]*k/chi.sq[jj,m]  
      } 
      p.log.p[[i]] <- p.emp[[i]]*log(p.emp[[i]]/matrix(prob.0[[i]],nsim,length(prob.0[[i]]),byrow = TRUE))
    }
  }
  
  for(jjj in 1:nsim){
    for (iii in 1:l){
      for (s in 1:ncol(unlist(p.log.p[[iii]]))){
        if(is.nan(p.log.p[[iii]][jjj,s])==TRUE){
          p.log.p[[iii]][jjj,s]=0                ## l'Hopital's Rule
        } 
      }
      g.sq.int[,iii] <- 2*n.0[iii]*rowSums(unlist(p.log.p[[iii]]))
      for (m in 1:(l/q) ){
        g.sq[,m] <- rowSums(g.sq.int[,((m-1)*q+1):(m*q)])
      }
    }
  }
  
  
  for (m in 1:(l/q)){
    qqplot(qchisq(ppoints(nsim),df=k), g.sq[,m],xlab = expression(paste(Chi^2, (k), ~"Quantile")),ylab=expression(paste(-2*log*Lambda)),cex.lab=1)
    qqline(distribution = function(p) qchisq(p,df=k),(g.sq[,m]))
    mtext(unique(n.0)[m])
  }
  
  
  chi.sq.N.emp.low <- matrix(0,nsim,l/q)
  chi.sq.N.emp.up <- matrix(0,nsim,l/q)
  count.N.emp <- matrix(0,nsim,l/q)
  count.N.hat <- matrix(0,nsim,l/q)
  V.inverse <- 1/(k-2)
  
  
  for (jj in 1:nsim){
    for (m in 1:(l/q)){
      chi.sq.N.emp.low[jj,m] <- qchisq(alpha/2,df=k)*V.inverse*N.emp[jj,m]  ## upper confidence interval for N.emp at alpha=0.05
      
      chi.sq.N.emp.up[jj,m] <- qchisq(alpha/2,df=k,lower.tail=FALSE)*V.inverse*N.emp[jj,m] ## lower confidence interval for N.emp at alpha=0.05
      
      if(n.0[m*q] > chi.sq.N.emp.low[jj,m]  & n.0[m*q] < chi.sq.N.emp.up[jj,m] )
      {
        count.N.emp[jj,m] <- 1
      }
      else{
        count.N.emp[jj,m]<- 0
      }
    }
  }
  
  for (jj in 1:nsim){
    for (m in 1:(l/q)){
      if(n.0[m*q] > (chi.sq.N.emp.low[jj,m]*(k-2)/k)  & n.0[m*q] < (chi.sq.N.emp.up[jj,m]*(k-2)/k))
      {
        count.N.hat[jj,m] <- 1
      }
      else{
        count.N.hat[jj,m] <- 0
      }
    }
  }
  
  
  ## table with resuls of actual coverage at alpha=0.05 of N.emp and N.hat=N.emp*(k-2)/k
  coverage <- rbind(colSums(count.N.emp)/nsim,colSums(count.N.hat)/nsim) 
  row.names(coverage) <- c("N.emp","N.emp(k-2)/k")
  return(coverage)
}

calc_coverage(q,prob.adj,k,n.0,l,alpha,nsim)

