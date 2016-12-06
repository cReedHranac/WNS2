# sensitvityHaymanWNS.R

###############################################################

library(lhs)

nspaces=100 ## how many bins/intervals, 100 needed for PRCC

hypercube=randomLHS(n=nspaces, k=17) ## function with N columns
dimnames(hypercube)[[2]]=c(
	"alpha",
	'beta',
	'delta',
	'gamma',
	'epsilon',
	'tau',
	'mu',
	'phi',
	'psi',
	'omega',
	'kappa',
	"lambda",
	"sigma",
	"theta",
	"nu",
	"pi",
	"rho")  # named columns with parameter names

mins = with(c(as.list(batParams),as.list(growthParams)),
c(   			    ## set mins for each parameters-exclude those not to be varied if any
	RMR=0.9*RMR,
	TMRmin=0.9*TMRmin,
	Teu=0.9*Teu,
	Tlc=0.9*Tlc,
	Ttormin=0.9*Ttormin,
	Ceu=0.9*Ceu,
	Ct=0.9*Ct,
	S=0.9*S,
	ttormax=0.9*ttormax,
	tar=0.9*tar,
	teu=0.9*teu,
	mass=0.9*mass,
	beta1=0.9*beta1,
	beta2=0.9*beta2,
	beta3=0.9*beta3,
	mu1=0.9*mu1,
	mu2=0.9*mu2
))

maxs = with(c(as.list(batParams),as.list(growthParams)),
c( 				    ## set maxs for each parameters-exclude those not to be varied if any
	RMR=RMR,
	TMRmin=TMRmin,
	Teu=Teu,
	Tlc=Tlc,
	Ttormin=Ttormin,
	Ceu=Ceu,
	Ct=Ct,
	S=S,
	ttormax=ttormax,
	tar=tar,
	teu=teu,
	mass=mass,
	beta1=beta1,
	beta2=beta2,
	beta3=beta3,
	mu1=mu1,
	mu2=mu2
))

diffs=maxs-mins ## range of each variable

hypercubeadj = hypercube # create copy of hypercube samples to modify, hypercube adjusted; i.e. new matrix
for (i in 1:ncol(hypercube)){
	hypercubeadj[,i]=hypercube[,i]*diffs[i]+mins[i] # scale samples to difference and add minimum
}

# head(hypercubeadj)
# dim(hypercubeadj)

dimnames(hypercubeadj)[[2]]=c("RMR",
															'TMRmin',
															'Teu',
															'Tlc',
															'Ttormin',
															'Ceu',
															'Ct',
															'S',
															'ttormax',
															'tar',
															'teu',
															"mass",
															"beta1",
															"beta2",
															"beta3",
															"mu1",
															"mu2")
paramset<-hypercubeadj

######################################################################################
results<-matrix(NA,nrow=100,ncol=1)

dynamic <- TRUE

# for one parameter set....
for (i in 1:length(paramset[,1])){
	
	## range of environmental parameters for the model
	
# 	Ta = seq(from=0,to=19.4,by=0.2)  # range of ambient temperatures
# 	Hd = seq(from=88,to=99, by=0.2) # tamge of ambient rel. Humidity %s
# 	twinter = seq(from=1, to = 5041, by =60) # up to 7 months
	Ta = seq(from=0,to=19.4,by=2)  # range of ambient temperatures
	Hd = seq(from=88,to=99, by=2) # tamge of ambient rel. Humidity %s
	twinter = seq(from=1, to = 5041, by =60) # up to 7 months
	
	parametset.H<- expand.grid(Ta=Ta,Hd=Hd,twinter=twinter)
	
	## this is to run the model
	
	if (dynamic) {
		res_out <- matrix(0, nrow(parametset.H), 2)
		for (j in 1:nrow(parametset.H)) {
			res = dynamicEnergyPd(Ta=parametset.H$Ta[j], Hd=parametset.H$Hd[j], twinter=parametset.H$twinter[j], WNS=T, modelParams=paramset[i,])
			res_out[j,] <- res
			if (j %% 100 == 0) {
				cat("up to iteration", j, "of", nrow(parametset.H),"\n")
			}
		}
		res <- c(res_out[,1], res_out[,2])
	} else {
		res = energy.Pd(Ta=parametset.H$Ta, Hd=parametset.H$Hd, twinter=parametset.H$twinter, WNS=T, modelParams=paramset[i,])
	}
	
	## re-structure results for plotting/analysis
  
	res<-as.matrix(res)
	
	## NB 1st half res = grammes fat required
	##    2nd half res = Pd growth
	
	gfat<-res[1:(length(res)/2),]
	Pdgrowth<-res[((length(res)/2)+1):length(res),]
	all <- cbind(Pdgrowth,gfat,parametset.H)
	all.names <-c("Pdgrowth","Gfat","Ta","Humidity","twin")
	colnames(all) <-all.names
#	print(all)
  
	## select gfat above which survival does not occur to estimate time to reach this value
	#
	#  all.t <- all[all$Gfat >= 0.3*mass,] # can vary, but this is from Kunz et al.
	#
	## select lowest values of twin for each combination
	#
	#  res.min<-ddply(all.t, .(Ta, Humidity), subset, twin==min(twin))
	#  dim(res.min)
	
	pp.mortality<-function(months)
	{
		## estimate sd for pnorm from Kunz
		sd_ML<-0.4/sqrt(6) # SEM/sqrt(sample size)
		## convert months to hours (to match model results)
		hours <- months*30*24+1 ## convert months to hours
		res<-all[all$twin==hours,] # choose 6 mo data
		mortality<-pnorm(res$Gfat, mean=0.3*paramset[i,12],sd=sd_ML) # from Kunz
		
		## calculate probability  
		mortality<-signif(mortality,2) # round, otherwise some asymptotes ~ = 0/1
		pmodel <-ifelse(mortality==1,1,0) # total 1s from sample
		pmodel <-sum(pmodel)/length(mortality) #
		
		return(pmodel) # to check value
	}
	month.t <-c(6) # define number of months want values for
	resp<-lapply(month.t,FUN=pp.mortality)
	resp<-as.numeric(resp)
	
	## store results
	results[i,]<-resp
	cat("up to iteration", i, "of", nrow(paramset),"\n")
}

################################################################################
#####Functions to calculate and plot partial-rank correlation coefficients
#####(PRCCs) between parameter values and model output.
#####
#####Written by: Michael Buhnerkempe
#####      Date: Oct. 7, 2011
#####
##### Functions:
#####          prcc - calculate the partial rank correlation coefficient
#####     plot.prcc - plot a prcc object
#####
##### A brief example is presented at the end
################################################################################

################################################################################
## prcc - function to calculate partial rank correlation coefficients between
##           each of p parameters and k model outputs using n different
##           observations (number of parameter sets)
##
##    Arguments:
##          par.mat = n x p matrix containing the parameter values obtained
##                        from Latin Hypercube Sampling
##      model.output = n x k matrix containing the model outputs
##           routine = how should the PRCCs be calculated? One of:
##                        "blower" - calculated according to Appendix A in
##                                   Blower and Dowlatabadi (1994). DEFAULT.
##                        "regression" - calculated using a regression approach.
##                                   Here, the partial correlation coefficient
##                                   is defined as the correlation between the
##                                   residuals after regressing a model output
##                                   on all of the parameters except for the
##                                   parameter of interest and the residuals
##                                   after regressing the parameter of interest
##                                   on all of the other parameters. This can
##                                   be interpreted as the correlation between
##                                   of a parameter and the model output when
##                                   the effects of all other parameters have
##                                   been removed.
##         par.names = names of parameters
##      output.names = names of model outputs
##               ... = additional arguments to be passed to functions called
##                     within this function
##
##
##    Output attributes:
##      $par.matrix = original matrix of parameter set
##      $model.output = original model output
##      $(model output name) = prcc results for the named model output

prcc = function( par.mat, model.output, routine = "blower",
								 par.names = NA,output.names = NA, ...){
	
	#Make sure par.mat and model.output are matrices
	par.mat = as.matrix(par.mat)
	model.output = as.matrix(model.output)
	
	#How many parameter sets are there?
	n = length(par.mat[,1])
	
	#How many parameters are there?
	p = length(par.mat[1,])
	
	#How many model outputs are we calculating PRCCs for?
	k = length(model.output[1,])
	
	#Find the ranks of the parameter values and the model output
	par.rank = apply(par.mat,2,rank,...)
	output.rank = apply(model.output,2,rank,...)
	
	#What is the average rank?
	ave.rank = (1 + n)/2
	
	#Create a list object to store the PRCCs
	results = list()
	
	results$num.sets = n
	
	#Try to automatically get parameter and output names if they are not
	#given
	if( sum(is.na(par.names)) > 0){par.names=dimnames(par.mat)[[2]]}
	if( sum(is.na(output.names)) > 0){output.names=dimnames(model.output)[[2]]}
	
	########################################################################
	#Calculate the PRCCs using Appendix A from Blower and Dowlatabadi (1994)
	########################################################################
	if( routine == "blower" ){
		
		#Do the calculation for each model output
		for( i in 1:k ){
			
			work.mat = cbind(par.rank,output.rank[,i])
			
			C.temp = matrix(0,nrow=p+1,ncol=p+1)
			
			#Calculate the C matrix
			for( j in 1:(p+1) ){
				for( l in 1:(p+1) ){
					
					C.temp[j,l]=sum((work.mat[,j]-ave.rank)*(work.mat[,l]-ave.rank))/
						sqrt(sum((work.mat[,j]-ave.rank)^2)*
								 	sum((work.mat[,l]-ave.rank)^2))
				}
			}
			
			#Calculate the B matrix (qr helps with inversion)
			B.temp = solve(qr(C.temp))
			
			coeff.val = rep(0,p)
			
			#Calculate the PRCC
			for( j in 1:p ){
				
				coeff.val[j] = -B.temp[j,p+1]/sqrt(B.temp[j,j]*B.temp[p+1,p+1])
				
			}
			
			#Calculate the t-test statistics and p-values
			t.val = coeff.val*sqrt((n-2)/1-coeff.val)
			p.val = 2*pt(abs(t.val),df=(n-2),lower.tail=F)
			
			#Output the results
			results[[output.names[i]]] = data.frame(
				prcc = coeff.val,
				t.value = t.val,
				p.value = p.val,
				row.names = par.names)
		}
		
		return(results)
	}
	
	########################################################################
	#Calculate the PRCCs using regression methods
	########################################################################
	else if( routine == "regression" ){
		
		#Do the calculation for each model output
		for( i in 1:k ){
			
			coeff.val = rep(0,p)
			
			#Calculate the PRCC
			for( j in 1:p ){
				
				#Variation in output that can not be explained by all other predictors
				#(except the predictor of interest)
				fit.y = lm(output.rank[,i] ~ par.rank[,-j])
				
				#Variation in the predictor of interest that can not be explained
				#by the other predictors
				fit.x = lm(par.rank[,j] ~ par.rank[,-j])
				
				#PRCC is the correlation between the residuals of the two
				#regressions above
				coeff.val[j] = cor(fit.y$residuals,fit.x$residuals)
				
			}
			
			#Calculate the t-test statistics and p-values
			t.val = coeff.val*sqrt((n-2)/1-coeff.val)
			p.val = 2*pt(abs(t.val),df=(n-2),lower.tail=F)
			
			#Output the results
			results[[output.names[i]]] = data.frame(
				prcc = coeff.val,
				t.value = t.val,
				p.value = p.val,
				row.names = par.names)
		}
		
		return(results)
	}
	
	else{ return("Error: Calculation type is invalid. Must be either 'blower' or 'regression'") }
	
}


################################################################################
## plot.prcc - function to plot a prcc object
##
##    Arguments:
##          prcc.obj = a prcc object from the 'prcc' function
##             alpha = level of significance desired for cut-off lines
##               ... = additional arguments to be passed to functions called
##                     within this function
##
##
##    Output:
##      A plot that has a bar graph for each output variable giving the PRCCs.
##      Dashed red lines give the cutoff values for significance. It parameter
##      names are specified correctly, axis labels will be smart.

plot.prcc = function(prcc.obj,alpha=0.05,...){
	
	x11()
	par(mfrow=c(ceiling((length(prcc.obj)-1)/3),min(c(length(prcc.obj)-1,3))))
	
	for( i in 2:length(prcc.obj) ){
		
		names.list=dimnames(results[[i]])[[1]]
		
		#Bar graph with the prcc values. The function in names.arg converts character
		#strings containing the names of greek letters to actual greek letters
		barplot(prcc.obj[[i]][,"prcc"],
						names.arg= expression(  alpha,
																		beta,
																		delta,
																		gamma,
																		epsilon,
																		tau,
																		mu,
																		phi,
																		psi,
																		omega,
																		kappa,
																		lambda,
																		sigma,
																		theta,
																		nu,
																		pi,
																		rho),
						main=names(results)[i],ylab="PRCC",cex.lab=1.1,...)
		
		#Plot lines to show the cutoff values for the alpha level of significance
		t.cutoff=qt(alpha/2,df=prcc.obj$num.sets-2)
		sig.cutoffs=c( (-(t.cutoff^2)-sqrt((t.cutoff^4) + 4*(prcc.obj$num.sets-2)*(t.cutoff^2)))/(2*(prcc.obj$num.sets-2)),
									 (-(t.cutoff^2)+sqrt((t.cutoff^4) + 4*(prcc.obj$num.sets-2)*(t.cutoff^2)))/(2*(prcc.obj$num.sets-2)))
		abline(h=sig.cutoffs,lty=2,col="red")
		abline(h=0,lty=1,col="black")
	}
}

colnames(results)<-c("Survival")

PRCCresults=prcc(par.mat=hypercube,model.output=results ## results matrix here...
								 ,routine="blower" # NB removed par names so uses symbols, add [par.names="",]
								 ,output.names=c("Survival"))

plot.prcc(PRCCresults,ylim=c(-1,1),cex.sub=1.2)

## all works if enough paramets sets.... [otherwise get singular matrix warning..]
write.csv(PRCCresults,file="prccresults.csv")

d.plot.prcc = function(prcc.obj,alpha=0.05,...){
	x11()
	par(mfrow=c(ceiling((length(prcc.obj)-1)/3),min(c(length(prcc.obj)-1,3))))
	for( i in 2:length(prcc.obj) ){
		names.list=dimnames(results[[i]])[[1]]
		
		#Bar graph with the prcc values. The function in names.arg converts character
		#strings containing the names of greek letters to actual greek letters
		barplot(prcc.obj[[i]][,"prcc"],
						names.arg= expression("RMR",
																	'TMR'[min],
																	'T'[eu],
																	'T'[lc],
																	'T'[tor-min],
																	'C'[eu],
																	'C'[t],
																	'S',
																	't'[tor-max],
																	't'[ar],
																	't'[eu],
																	"mass",
																	beta[1],
																	beta[2],
																	beta[3],
																	mu[1],
																	mu[2]),
						main=names(results)[i],las=2,ylab="PRCC",cex.lab=1.1,...)
		#Plot lines to show the cutoff values for the alpha level of significance
		t.cutoff=qt(alpha/2,df=prcc.obj$num.sets-2)
		sig.cutoffs=c( (-(t.cutoff^2)-sqrt((t.cutoff^4) + 4*(prcc.obj$num.sets-2)*(t.cutoff^2)))/(2*(prcc.obj$num.sets-2)),
									 (-(t.cutoff^2)+sqrt((t.cutoff^4) + 4*(prcc.obj$num.sets-2)*(t.cutoff^2)))/(2*(prcc.obj$num.sets-2)))
		abline(h=sig.cutoffs,lty=2,col="red")
		abline(h=0,lty=1,col="black")
	}
}
d.plot.prcc(PRCCresults,ylim=c(-1,1),cex.sub=1.2)

parcoln<-expression("RMR",
           'TMR'[min],
           'T'[eu],
           'T'[lc],
           'T'[tor-min],
           'C'[eu],
           'C'[t],
           'S',
           't'[tor-max],
           't'[ar],
           't'[eu],
           "mass",
           beta[1],
           beta[2],
           beta[3],
           mu[1],
           mu[2])

pdf(paste("prcc_linear.pdf"),width=12, height=12)
par(mfrow=c(5,4))
for (i in 1:ncol(paramset))
 # plot(paramset[,i], results, main=colnames(paramset)[i])
  plot(paramset[,i], results, xlab=parcoln[i],ylab="PRCC",cex.lab=1.2)
dev.off()

pdf(paste("prcc_results.pdf"),width=10, height=10)
d.plot.prcc(PRCCresults,ylim=c(-1,1),cex.sub=1.2)
dev.off()

