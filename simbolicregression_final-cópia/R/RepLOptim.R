RepLOptim <- function(parmean,parsd,fr,gr=NULL,inphess=NULL,lower=NULL,upper=NULL,
	rethess=FALSE,parmstder=FALSE,...,control=list()) 

#   RepLOptim -- Repeated local optimization 

#   Tries to minimize a function calling local optimizers several times from different random starting points
#   generated from multivariate normal distributions of independent variates. The standard deviations of the 
#   generating distributions are kept fixed, but their means are updated as better candidates for the global 
#   minimun are discovered

#   Arguments:

#   parmean   --  Vector of means for the parameter distribution generating starting points for the
#                 local optimizer. Also used as starting point of the first call to the optimizer.
#   parsd     --  Vector of standard deviations for the parameter distribution generating starting points 
#                 for the local optimizer. Also used as starting point of the first call to the optimizer.
#   fr        --  The function to be minimized. If method is neither "nlminb" or "L-BFGS-B", fr should
#                 accept a lbound and an ubound arguments for the parameter bounds, and should enforce
#                 these bounds before calling the local optimization routine.
#   gr        --  A function to return the gradient for the "nlminb", "BFGS", '"CG"'and '"L-BFGS-B"' methods.  
#                 If it is 'NULL', a finite-difference approximation will be used. For the '"SANN"' method 
#                 it specifies a function to generate a new candidate point.  If it is 'NULL' a default Gaussian
#                 Markov kernel is used.
#   inphess   --  A function to return the hessian for the "nlminb" method. Must return a square matrix of order 
#                 'length(parmean)' with the different hessian elements in its lower triangle. 
#                 It is ignored if method component of the control list is not set to its "nlminb" default. 
#   lower     --  Vector of parameter lower bounds. Set to -Inf (no bounds) by default
#   upper     --  Vector of parameter upper bounds. Set to Inf (no bounds) by default
#   rethess   --  Bolean flag indicating wether a numerically evaluated hessian matrix at the optimum 
#                 should be computed and returned. Not available for the nlminb method.
#   parmstder --  Bolean flag indicating wether parameter assymptotic standard errors based on the 
#                 inverse hessian approximation to the Fisher information matrix should be computed and returned.
#                 Only available if hessian is set to "true" and if a local miminum with a positive-definite hessian 
#                 was indeed found. For non-trivial problems of moderate/large dimensionality this requirement often 
#                 fais because of numerical problems.
#                 Note that the consistency of these estimates can only be assured when fr is based on a true 
#                 likelihood (as oposed to a quasi-likelihood) function
#   ...       --  Further arguments to be passed to fr, gr or to the local optimization routine selected

#  control   --  A list of control parameters that can supply any of the following components:
#    method  --     Local optimizer to be employed. Current alternatives are:
#                   "nlminb" (default) for the nlminb port routine, "nlm" for the nlm function and
#                   "Nelder-Mead", "L-BFGS-B",  "CG", "L-BFGS-B" and "SANN" for the corresponding methods of
#                   the optim function. 
#    maxrepet --    Maximum time of repetions of the same minimum objective value, before RepLOptim is stoped
#                   and the current best solution is returned. By default set to 5. 
#    maxnoimprov -- Maximum number of times the local optimizer is called without improvements in the minimum 
#                   objective value, before RepLOptim is stoped and the current best solution is returned. 
#                   By default set to 25. 
#    maxreplic  --  Maximum number of times the local optimizer is called and returns a valid solution before 
#                   RepLOptim is stoped and the current best solution is returned. By default set to 250. 
#    allrep     --  Total maximum number of replications (including those leading to non-valid solutions)
#                   performed. By default equals ten times the value of maxreplic. Ignored when objbnd is set to Inf. 
#    maxiter    --  Maximum number of iterations performed in each call to the local optimizer. By default set to 500
#                   execept with the "SANN" mehtod, when by default is set to 1500. 
#    maxeval    --  maximum number of function evaluations (nlminb method only) performed in each call to 
#                   the nlminb optimizer. By defaults set to 750.
#    RLOtol     --  The relative convergence tolerance of the local optimizer.  The local optimizer stops if it is 
#                   unable to reduce the value by a factor of ‘RLOtol *(abs(val) + reltol)’ at a step. 
#                   Ignored when method is set to "nlm". By default set to 1e-8.
#    HesEgtol   --  Numerical tolerance used to ensure that the hessian is non-singular. If the last eigenvalue 
#                   of the hessian is positive but the ratio between it and the first eigenvalue is below HesEgtol, 
#                   the hessian is considered to be semi-definite and the parameter assymptotic standard errors 
#                   are not computed. By default set to 1e-8.
#    objbnd    --   Upper bound for the objective. When different form "Inf" (default) only solutions leading to objective values #                   below objbnd are considered as valid.

#   Value:  A list with the following components
#
#         par            --   The best result found for the parameter vector.
#         val            --   The best value (minimum) found for the function fr.
#         iterations     --   Number the iterations performed by the local optimizer in the call
#                             that generated the best result.
#         vallist        --   A vector with the best values found for each starting point.
#         counts         --   Number of times the function fr was evaluated in the call that 
#                             generated the result returned.
#         convergence    --   Code with the convergence status returned by the local optimizer.
#         message        --   Message generated by the local optimizer.
#         hessian        --   Numerically evaluated hessian of fr at the result returned. 
#                             Only returned when the parameter hessian is set to true. 
#         hessegval      --   Eigenvalues of the hessian matrix. Used to confirm if a local minimum was indeed found.
#                             Only returned when the parameter hessian is set to true. 
#         stderrors      --   Assymptotic standard deviations of the parameters based on Fisher Information 
#                             matrix. Only returned when the parse parameter is set to true and the hessian is indeed
#                             positive definite.

{
	method_default <- "nlminb"
	maxrepet_default <- 5
	maxnoimprov_default <- 25
	maxreplic_default <- 250
	maxiter_default <- 500
	maxSANNiter_default <- 1500
	maxeval_default <- 750
	RLOtol_default <- 1e-8
	HesEgtol_default <- 1e-8
 
	if (!is.null(control$method)) method <- control$method else method <- method_default  
	if (!is.null(control$maxrepet)) maxrepet <- control$maxrepet else maxrepet <- maxrepet_default  
	if (!is.null(control$maxnoimprov)) maxnoimprov <- control$maxnoimprov else maxnoimprov <- maxnoimprov_default
	if (!is.null(control$maxreplic)) maxreplic <- control$maxreplic else maxreplic <- maxreplic_default
	if (!is.null(control$maxiter)) maxiter <- control$maxiter 
	else if (class(method)!="character") maxiter <- maxiter_default
		else if (method!="SANN") maxiter <- maxiter_default
			else maxiter <-  maxSANNiter_default
	if (!is.null(control$maxeval)) maxeval <- control$maxeval else maxeval <- maxeval_default
	allrep <- control$allrep
 	if (!is.null(control$RLOtol)) RLOtol <- control$RLOtol else RLOtol <- RLOtol_default  
 	if (!is.null(control$HesEgtol)) HesEgtol <- control$HesEgtol else HesEgtol <- HesEgtol_default  
 	if (!is.null(control$objbnd)) objbnd <- control$objbnd else objbnd <- -Inf 
 
	if (!is.null(control$parmean)) parmean <- control$parmean 
	if (!is.null(control$sdfactor)) 
		if (!is.null(control$parsd)) parsd <- control$sdfactor * control$parsd 
 		else parsd <- control$sdfactor * parsd 
 	else if (!is.null(control$parsd)) parsd <- control$parsd 
 	if (!is.null(control$inphess)) inphess <- control$inphess 
 	if (!is.null(control$lower)) lower <- control$lower 
 	if (!is.null(control$upper)) upper <- control$upper 
 	if (!is.null(control$rethess)) lower <- control$rethess 
 	if (!is.null(control$parmstder)) parmstder <- control$parmstder 
 	if (!is.null(control$EnfCnstrs)) EnfCnstrs <- control$EnfCnstrs 
  
	npar <- length(parmean)
	values <- NULL
	if (is.null(lower)) lower <- rep(-Inf,npar)
	if (is.null(upper)) upper <- rep(Inf,npar)
	if (is.finite(objbnd))  { if (is.null(allrep)) allrep <- 10*maxreplic }
	else allrep <- maxreplic

	bestres <- NULL	
	bestval <- Inf
	initpar <- bestpar <- parmean
	cnt <- repcnt <- noimpcnt <- 0 
	for (i in 1:maxreplic)  {
		if (cnt > allrep || repcnt >= maxrepet || noimpcnt >= maxnoimprov) break
		value <- Inf
		while (value >= objbnd && cnt < allrep && repcnt < maxrepet && noimpcnt < maxnoimprov)
		{
			if (is.function(method))
				tmpres <- method(initpar,gr=gr,lbound=lower,ubound=upper,hessian=rethess,...)
			else {
				if (method == "nlminb")
					tmpres <- nlminb(start=initpar,fr,gradient=gr,hessian=inphess,
						lower=lower,upper=upper,control=list(iter.max=maxiter,eval.max=maxeval),...)

				else if (method == "nlm") 
					tmpres <- nlm(fr,p=initpar,lbound=lower,ubound=upper,iterlim=maxiter,...)
				else if (method == "L-BFGS-B")
					tmpres <- optim(initpar,fr,gr=gr,method=method,lower=lower,upper=upper,
						control=list(maxit=maxiter),hessian=rethess,...)
				else if (method == "Nelder-Mead" || method == "BFGS" || method == "CG" || method == "SANN")
					tmpres <- optim(initpar,fr,gr=gr,method=method,control=list(maxit=maxiter),
						lbound=lower,ubound=upper,hessian=rethess,...)
			}
			if (!is.function(method))  {
				if (method == "nlminb") value <- tmpres$objective
				else if (method == "nlm") value <- tmpres$minimum
				else value<- tmpres$value
			}
			else value <- tmpres$value
			if (is.null(value) || is.na(value)) value <- objbnd
			cnt <- cnt+1
			values <- c(values,value)
			if (is.na(value)) { noimpcnt <- noimpcnt + 1 ; repcnt <- 0 }
			else {
				if (is.finite(bestval))  if (abs((value-bestval)/bestval) < RLOtol) repcnt <- repcnt + 1 
				else repcnt <- 0 
				if (value < bestval)  {
					bestval <- value
					if (is.function(method) || method != "nlm") bestpar <- tmpres$par
					else bestpar <- tmpres$estimate
					bestres <- tmpres
					noimpcnt <- 0 
				}
				else noimpcnt <- noimpcnt + 1
			}
			if (value >= objbnd && cnt < allrep) {
				u <- runif(n=npar)     # generate npar uniform random numbers
				initpar <- qnorm(u,mean=bestpar,sd=parsd) #  generate new parameters from a normal distribution
				lbndind <- initpar < lower   #  identify indices of parameters that fell below their lower bounds
				ubndind <- initpar > upper   #  identify indices of parameters that fell above their upper bounds
				initpar[lbndind] <- lower[lbndind] + u[lbndind] * (bestpar[lbndind]-lower[lbndind]) # and correct them
				initpar[ubndind] <- upper[ubndind] - u[ubndind] * (upper[ubndind]-bestpar[ubndind])
			}
		} 
	}

	if (!is.function(method) && method == "nlminb")  {
		iterations <- bestres$iterations
		counts <- bestres$evaluations
		hess <- NULL
		egval <- NULL
		parstd <- NULL
	} 
	else  {
		iterations <- NULL
		counts <- bestres$counts
		if (rethess==TRUE) {
			hess <- bestres$hessian
			egval <- eigen(hess,symmetric=TRUE,only.values=TRUE)$values
			if (parmstder==TRUE)
				if (egval[npar]/egval[1] < HesEgtol)
					parstd <- "Not computed because the hessian is not positive definite"
				else parstd <- sqrt(diag(solve(hess)))
			else parstd <- NULL
		}
		else {
			hess <- NULL
			egval <- NULL
			parstd <- NULL
		}
	} 
	if (!is.null(bestres))
		return( list(par=bestpar,val=bestval,iterations=iterations,vallist=values,counts=counts,
			convergence=bestres$convergence,message=bestres$message,hessian=hess,hessegval=egval,stderrors=parstd) )
	else
		return( list(par=NULL,val=Inf,iterations=NULL,vallist=NULL,counts=NULL,convergence=NULL,
			message="RepLOptim was unable to find any valid solution",hessian=NULL,hessegval=NULL,stderrors=NULL) )
}
