###Functions for compadre
library(popbio)
library(popdemo)
library(Matrix)
library(MASS)
library(maps)
library(scales)
library(devtools)
library(ape)
library(caper)
library(phytools)


##------------------------------------------------------------------------
#mean centering function
mean_center <- function(x){ (x - mean(x))/sd(x)
}

##------------------------------------------------------------------------

lifeTimeRepEvents <- function(matU, matF, startLife = 1){

#Function to determine probability of reaching reproduction, age at maturity and reproductive lifespan (Code adapted from H. Caswell's matlab code):

uDim = dim(matU)[1]
surv = colSums(matU)
repLifeStages = colSums(matF)
repLifeStages[which(repLifeStages>0)] = 1

if(missing(matF) | missing(matU)){stop('matU or matF missing')}
if(sum(matF,na.rm=T)==0){stop('matF contains only 0 values')}

#Probability of survival to first reprod event
Uprime = matU
Uprime[,which(repLifeStages==1)] = 0
Mprime = matrix(0,2,uDim)
for (p in 1:uDim[1]) {
if (repLifeStages[p]==1) Mprime[2,p] = 1 else
Mprime[1,p] = 1-surv[p]
}

Bprime = Mprime%*%(ginv(diag(uDim)-Uprime))
pRep = Bprime[2,startLife]

out = data.frame(pRep = pRep)
#Age at first reproduction (La; Caswell 2001, p 124)
D = diag(c(Bprime[2,]))
Uprimecond = D%*%Uprime%*%ginv(D)
expTimeReprod = colSums(ginv(diag(uDim)-Uprimecond))

La = expTimeReprod[startLife]
out$La = La

#Mean life expectancy conditional on entering the life cycle in the first reproductive stage
firstRepLifeStage = min(which(repLifeStages==1))

N = solve(diag(uDim[1])-matU)

meanRepLifeExpectancy = colSums(N)[firstRepLifeStage]

out$meanRepLifeExpectancy = meanRepLifeExpectancy

#Life expectancy from mean maturity
remainingMatureLifeExpectancy = colSums(N)[startLife]-La

out$remainingMatureLifeExpectancy = remainingMatureLifeExpectancy

return(out)

}


##------------------------------------------------------------------------

#Calculates mean life expectancy
#From Owen Jones
meanLifeExpectancy <- function(matU = matU, startLife = 1){
  uDim=dim(matU)[1]
  N = solve(diag(uDim[startLife])-matU)
  eta = colSums(N)[startLife]
  return(eta)
}


##------------------------------------------------------------------------


makeLifeTable <- function(matU, matF = NULL, matC = NULL, startLife = 1, nSteps = 1000){
  
  matDim = ncol(matU)
  
  #Age-specific survivorship (lx) (See top function on page 120 in Caswell 2001):
  matUtemp = matU
  survivorship = array(NA, dim = c(nSteps, matDim))
  for (o in 1:nSteps){
    survivorship[o, ] = colSums(matUtemp %*% matU)
    matUtemp = matUtemp %*% matU
  }
  
  lx = survivorship[, startLife]
  lx = c(1, lx[1:(length(lx) - 1)])

  #Make room for dx and qx under assumption of 0.5 in age gap distributions
  
  #Start to assemble output object
  out = data.frame(x = 0:(length(lx)-1),lx = lx)
  
  if(!missing(matF)){
    if(sum(matF,na.rm=T)==0){
      warning("matF contains only 0 values")
    }
  #Age-specific fertility (mx, Caswell 2001, p. 120)
  ageFertility = array(0, dim = c(nSteps, matDim))
  fertMatrix = array(0, dim = c(nSteps, matDim))
  matUtemp2 = matU
  e = matrix(rep(1, matDim))
  for (q in 1:nSteps) {
    fertMatrix = matF %*% matUtemp2 * (as.numeric((ginv(diag(t(e) %*% matUtemp2)))))
    ageFertility[q, ] = colSums(fertMatrix)
    matUtemp2 = matUtemp2 %*% matU
  }  
  mx = ageFertility[, startLife]
  mx = c(0, mx[1:(length(mx) - 1)])
  out$mx = mx
  }
  
  if(!missing(matC)){

  #Age-specific clonality (cx)
  ageClonality = array(0, dim = c(nSteps, matDim))
  clonMatrix = array(0, dim = c(nSteps, matDim))
  matUtemp2 = matU
  e = matrix(rep(1, matDim))
  for (q in 1:nSteps) {
    clonMatrix = matC %*% matUtemp2 * (as.numeric((ginv(diag(t(e) %*% matUtemp2)))))
    ageClonality[q, ] = colSums(clonMatrix)
    matUtemp2 = matUtemp2 %*% matU
  }  
  cx = ageClonality[, startLife]
  cx = c(0, cx[1:(length(cx) - 1)])
  out$cx = cx
  }
  
  return(out)
  }


##------------------------------------------------------------------------
#calulculate expectional lifespan


exceptionalLife<-function(matU,startLife=1){
  popVector=rep(0,dim(matU)[1])
  popVector[startLife]=100
  lifespanLeftover=matrix(0,1000,1)
  for (n in 1:1000){
    lifespanLeftover[n]=sum(popVector)
    popVector=matU%*%popVector
  }
  Lexcept.95=min(which(lifespanLeftover<5))
  if(Lexcept.95==Inf) {Lexcept.95=999}
  Lexcept.99=min(which(lifespanLeftover<1))
  if(Lexcept.99==Inf) {Lexcept.99=999}
  
  return(c(Lexcept.95,Lexcept.99))
}


##------------------------------------------------------------------------
#calculate mean of a list of matrices

meanMatrix <- function(Matlist){
  
  mat_array <- do.call(cbind, Matlist)
  mat_array <- array(mat_array, dim=c(dim(Matlist[[1]]), length(Matlist)))
  mean_matrix <- apply(mat_array, c(1, 2), mean, na.rm = TRUE)
  
  return(mean_matrix)
}





##------------------------------------------------------------------------
####check each matrix to make sure its irriducable, primitive and ergodic

matrix_epi_check <- function(comadre_rows = c(), combMat = c() ){
  
  
  is_ergodic <- vector()
  is_primitive <- vector()
  is_irreducible <- vector()
  all_true <- vector()
  
  for(i in 1:length(comadre_rows)){
    tryCatch({
      is_ergodic[i] <- is.matrix_ergodic(combMat[[i]]$matA)
      
      is_primitive[i] <- is.matrix_primitive(combMat[[i]]$matA)
      
      is_irreducible[i] <- is.matrix_irreducible(combMat[[i]]$matA)
      
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    all_true[i] <- all(c(is_ergodic[i],is_primitive[i],is_irreducible[i]) == TRUE)
  } 
  
  test.frame <- data.frame(is_ergodic,is_primitive,is_irreducible)
  keep <- which(all_true == TRUE)
  discard_species <- which(all_true == FALSE)
  discard_species_names  <- unique(combined_data[discard_species,]$SpeciesAccepted)
  
  clean_species <- combined_data[keep,]$SpeciesAccepted
  
  return(list(is_epi = test.frame,
              all_true = all_true,
              clean_species = clean_species,
              clean_COMADRE = keep,
              discard_species_names = discard_species_names
  ))
}





##------------------------------------------------------------------------
#check if final colume is 0

is.matrix_post_rep <- function(A){
	
	N <- dim(A)[1]
	as.vector(A[N,N] > 0)
	
}


##------------------------------------------------------------------------
#calculate mean reporductive rate

meanRepo <- function(matA, matF){
   
    N <- dim(matF)[1]
	stable_stage <- eigen.analysis(matA)$stable.stage
  	repo_sum <- vector()
    for(j in 1:(N)){
      repo_sum[j] <- sum(matF[,j])
    }
    mean_repo_rate <- as.vector(repo_sum %*% stable_stage)
return(mean_repo_rate)
}


##------------------------------------------------------------------------
##plotting function


MultiDisPlot<-function(data.list, probs=c(95, 50), col="black", xlim="auto", ...) {
  #Sanitizing
  
  #require
  require(hdrcde)
  
  #data.list
  if(class(data.list) != "list") {
    stop("'data.list' must be a list of numerical values.")
  }
  
  #probs
  if(class(probs) != "numeric") {
    stop("'probs' must be a numeric.")
  }
  
  #col
  if(class(col) != "character") {
    stop("'col' must be a character string.")
  } 
  
  #Calculate the y limits (optional
  if(xlim == "auto") {
    xlim<-c(min(unlist(data.list), na.rm=TRUE) - 0.01*min(unlist(data.list), na.rm=TRUE), max(unlist(data.list), na.rm=TRUE) + 0.01*(max(unlist(data.list), na.rm=TRUE)))
  }
  
  #Calculating the hdr
  hdr.results<-lapply(data.list, hdr, prob=probs)
  
  #Empty plot
  plot(1,1, ylim= c(1,length(data.list)), xlim=xlim, col="white", ...)
  
  #Adding the lines
  for (j in 1:length(data.list)) {
    temp <- hdr.results[[j]]
    shift=0
    #border="black"
    
    for (k in 1:length(probs)) {
      #Plot the probabilities distribution
      temp2 <- temp$hdr[k, ]
      
      #Lines options
      lines(c(j+shift,j+shift)~ c(min(temp2[!is.na(temp2)]), max(temp2[!is.na(temp2)])), lwd=1+(k*2-2), lty=(length(probs)-(k-1)), col=col)  
      points(j+shift~temp$mode,pch=19, col=col)
    }
  }
  
}


##------------------------------------------------------------------------


#lx_spline
#fit a spline to predict lax values in order to remove time effect.

lx_spline <- function(lx, x, bins = 10){
  
  if(max(lx) > 1) stop("`lx` should be bounded between 0 and 1")
  if(sum(is.na(lx))>1) stop("There are missing values in `lx`")
  #if(sum(!diff(lx) <= 0)) stop("`lx` does not monotonically decline")
  if(sum(!diff(lx) <= 0)!=0)stop("`lx` does not monotonically decline")
  
  
  log_lx <- log10(lx*1000)
  
  if(any(is.infinite(log_lx)) == TRUE){
    log_lx[is.finite(log_lx) == FALSE] <- 0
    log_lx[which(log_lx <= 0)[1]:length(log_lx)] <- 0
  }
  
    
    lx_bins <- seq(0, max(x), max(x)/bins)
    
    lxlog2 <- approx(log_lx ~ x, xout = lx_bins)$y
    lx2 <- approx(lx ~ x, xout = lx_bins)$y
    
   
    lx2[lx2 > 1] <- 1
    lxlog2[lxlog2 > 3] <- 3
    lxlog2[lxlog2 < 0] <- 0
    
    if(any(lx2 <= 0) == TRUE){
      lx2[lx2 <= 0] <- 1e-09
      lx2[which(lx2 == 1e-09)[1]:length(lx2)] <- 1e-09
    }
    
    
    return(list(lx2,lxlog2))
}



##------------------------------------------------------------------------
#mean centering function
mul_resids <- function(mul_output, mul_data, Y_data_col){
  
  Expected_Y <- hdr(mul_output$`(Intercept)`)$mode + 
    (mul_data$data$mass_g)*hdr(mul_output$mass_g)$mode +
    (mul_data$data$matrix_size)*hdr(mul_output$matrix_size)$mode
  
  resids_mul <- mul_data$data[,Y_data_col]  - Expected_Y
  return(resids_mul)
  
}


