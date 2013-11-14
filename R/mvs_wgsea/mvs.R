##' Calculate a multivariate sample matrix using LD
##' 
##' As an alternative to case/control permutation, sampling multivariate normal 
##' distribution is computationally quicker as well as requiring only sample
##' population genotypes.
##' 
##' @param 
##' @param
##'
##' @return 
##' @author Olly Burren
##' @export
##' @seealso \code{\link{wilcoxon}}
##' @keywords htest
##' @examples
##'

library(snpStats)
library(corpcor)
library(mvtnorm)


##' For a given snpMatrix object generate a positive definite matrix suitable 
##' for generation of sigma for multivariate sampling.

mvs.sigma.ld<-function(gt){
	if(!is(gt,"SnpMatrix"))
		stop("Parameter is not a snpMatrix object!")
	#compute ld - note that for large amount of SNPs this might be computationally
	#expensive
	if((ncol(gt)-1)<=0)
		stop("Cannot calculate sigma for less than 2 snps")	
	ld <- ld(gt, stats=c("R.squared"), depth=ncol(gt)-1,symmetric=TRUE)
	## set missing values to 0
	ld[which(is.na(ld))]<-0;
	## here we attempt different values for diag in the
	## attempt to get a positive definite matrix
	
	## set diag  equal to 1 to have a shot at being positive definite matrix
	diag(ld)<-1
	if(!is.positive.definite(ld,,method="chol")){
		#this recurses through various values of diag if we exceed 1 then
		#we compute the closest matrix that is positive definite.
		ld<-attempt.pos.def(ld)
	}
	ld
}

##' For a given matrix recurse to find the closest positive definite matrix.

attempt.pos.def<-function(mat,diag.val=1.0001){
  print(paste("diag.val",diag.val))
	if(!is(mat,"Matrix"))
		stop("mat is not a Matrix!")
	if(diag.val >= 1.1){
	  print("Matrix is not positive definite. Finding closest approximation..")
		diag(mat)<-1
		return(as(make.positive.definite(mat),"Matrix"))
	}
	diag(mat)<-diag.val
	if(is.positive.definite(mat,,method="chol")==FALSE){
	  new.diag<-signif(1+((diag.val-trunc(diag.val))*10))
		mat<-attempt.pos.def(mat,new.diag)
	}else{
		return(mat)
	}
}

## this function takes a list of sigma matrices, number of 
## perms to perform and an optional character vector containing
## the snps to return, if the latter is ommitted get perms
## for all snps in sigma

compute.mvs.perms<-function(sigma,n.perms=1000,snp,names){
	do.call("rbind",lapply(seq_along(sigma),function(i){
	  mvs<-matrix()
	  x<-sigma[[i]]
	  flag=F;
	  ## regions with no snps in 
	  if(!is(x,"Matrix"))
	  	return()
	  ## region with 1 snp in just use normal random sampling
	  if(length(x)==1){
	  	mvs<-t(as.matrix(exp(-rexp(n.perms))))
	  ## region with more than one snp in compute mvs 
	  ## note we handle the case where sigma is still
	  ## not positve definite 
	  }else{
	  	withCallingHandlers(
       mvs<-mvs.perm(as.matrix(x),n.perms),
	  			warning=function(w) {
	  				print(paste("WARNING: mvs problem sigma region :",names(sigma)[i],"proably not positive definite trying again"))
	  				x<-attempt.pos.def(x)
	  				mvs.perm(as.matrix(x),n.perms)
	  				#attempt to make positive definite
	  				invokeRestart("muffleWarning")
	  			}
			)
    }
    
    rownames(mvs)<-rownames(x)
    ## in general our reference sigma will include
    ## many more SNPs than we are interested in, here
    ## we exclude snps that are not being tested
    if(length(snp.names) > 0){
    	index<-which(rownames(mvs) %in% snps.gr$name)
    	mvs[index,]
    }else{
    	mvs
    }
  }))
}

## compute permutations using mvs on a given correlation matrix sigma n times

mvs.perm<-function(sigma,n=1000){
	if(!is.matrix(sigma))
		stop("sigma parameter is not a matrix")		
	if(!is.positive.definite(sigma,,method="chol"))
		stop("sigma is not positive definite")
	## in original paper method="chol" was not defined so I assume used eigen default
	## this is slower than the choleski decomp ! Perhaps we should contact the author ?
	rd<-rmvnorm(n,mean=rep(0,ncol(sigma)),sigma=sigma,method="chol")
	t(rd)
}

