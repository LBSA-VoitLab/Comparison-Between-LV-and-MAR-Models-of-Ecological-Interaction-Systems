### LV artificial data ###


rm(list = ls())


######################################################################################################################################
##### Smoother
#
### function to find a dataset apropriate smooth function 
# for theorical background see https://www.researchgate.net/publication/285256327_Parameter_estimation_in_canonical_biological_systems_models
#

#install.packages("splines")
library(splines)

#install.packages("fANCOVA")
library(fANCOVA)

# dataset = cbind(x1[,1],x1[,2],x2[,2]); colnames(dataset) = c('t','X1','X2')
# draw=TRUE
# data1_spline2 = 2
# smooth_type = 1 
# df = c(9,9)
# dataWeights = NULL
# splineMethod = "fmm"
# polyDegree012 = 1
# aicc1_gvc2 = 1
# span = NULL
# log_spline = FALSE

Smoother = function (
			dataset,
			draw = TRUE,
			data1_spline2 = 1, 
			smooth_type = 1,
			df = NULL,
			dataWeights = NULL, 
			splineMethod = "fmm",
			polyDegree012 = 1,
			aicc1_gvc2 = 1,
			span = NULL,
			log_spline = FALSE
			)           
	{
	### function to estimate parameters for Lotka Volterra system.
	# dataset # the data. First colunm is the timepoints, remaining colunms the dependent variables. All colunms should be named.
	# draw # TRUE draw the estimated slopes for every data point. 
	# data1_spline2 # 1 - data will be used to estimate parameters # 2 - spline samples will be used to estimate parameters 
	# smooth_type # 1 - smoothing spline that needs degrees of freedom to be defined. 2 - LOESS.
	# df # degrees of freedom for smoothing splines 2, must have 1 < df <= n. If NULL, df = n.
	# dataWeights # NULL - all datapoints have the same weights # 'heavyStart' - initial values 100 times heavier that other datapoints # it acepts handcrafted weights
	# splineMethod # method to estimate splines # "fmm" "periodic" "natural" # monotonic splines "hyman" "monoH.FC"
	# polyDegree012 # the degree of the local polynomials to be used. It can ben 0, 1 or 2.
	# aicc1_gvc2 # the criterion for automatic smoothing parameter selection: “aicc” denotes bias-corrected AIC criterion, “gcv” denotes generalized cross-validation. 
	# span # value from [0,1] to define the percentage of sample points in the moving interval
	# log_spline # TRUE uses log transform datapoints to create the splines. Splines will by in cartesian values.

	data_time = dataset[,1]														# create a time var to use in splines and slopes
	numVars = dim(dataset)[2]-1													# set the number of dependent variable	
	if (is.null(df)) {df=rep(dim(dataset)[1],numVars)}									# if df are not defined, go for a spline that passes in all points.	
	data_vars = matrix(rep(NA,dim(dataset)[1]*(dim(dataset)[2]-1)),ncol = numVars)				# unlisting the data
	for (iii in 1:(dim(dataset)[2]-1)) 
		{
		data_vars[,iii] = unlist(dataset[,iii+1])
		}
	colnames(data_vars) = colnames(dataset[,2:(dim(dataset)[2])])							# naming the data_vars colunms
	if (log_spline == TRUE) {data_vars = log(data_vars)}									# use log var values to calculate the splines if requested
	dataFrame = data.frame(t=dataset[,1],data_vars)							  			# create data frame for the data to use in loess  

	if ( data1_spline2 == 1 )
		{
		time_splines = data_time
		} else
		{
		time_splines = round(seq(head(data_time,1),tail(data_time,1),length.out = 10*length(data_time)),10)	# create a time sequence with a smaller resolution
		time_splines = unique(sort(c(data_time,time_splines)))
		}

	if ( !is.null(dataWeights) & length(dataWeights)==1)
		{ if ( dataWeights == 'heavyStart' ) { dataWeights = rep(1,dim(dataset)[1]);dataWeights[1]=100 } }	# set the weigths for when the initial values are 100 times heavier than the other data points

	# splines and LOESS (smoothing) 
    	splines_est = slopes_est = d2X_dt2_est = time_splines									# initiating the vars to store the results
	for (i in 1:numVars)                                                                            	# cycling thought the dependent vars
      	{
        	if (smooth_type == 1)
            	{
            	smoothSpline <- smooth.spline(data_time,data_vars[,i],w=dataWeights,df=df[i])                            				# smoothing with degrees of freedom (df)
            	smooth = predict(object = smoothSpline, x = time_splines, deriv = 0)                                 			# get the points of the fit of that linear model
            	f_of_x <- splinefun(smooth$x,smooth$y,method = splineMethod )  # "fmm" "periodic" "natural" "hyman" "monoH.FC"                                               				# create cubic spline for the linear model
            	}
	  	if (smooth_type == 2)
			{
  			loess1 = loess.as(dataFrame[,1], dataFrame[,i+1], 
						degree = polyDegree012, 
						criterion = c("aicc", "gcv")[aicc1_gvc2], 
						user.span = span, plot = F )										# create loess of the data points with span = .7
               	smooth = loess1$fit														# store the loess fit to build the spline
                	f_of_x <- splinefun(dataFrame[,1],smooth)                                     			# create cubic spline for a dependent variable
			print(summary(loess1))
			}

		if (log_spline == FALSE)
			{		
			splines_est = cbind( splines_est, f_of_x(time_splines) )											# store spline points for all dep. variables  	
			slopes_est = cbind( slopes_est, f_of_x(time_splines, deriv = 1) ) 							# store the slopes of the data points for all dep. variables	
			d2X_dt2_est = cbind( d2X_dt2_est, f_of_x(time_splines, deriv = 2) ) 							# store the 2nd derivative estimates
			} else
			{
			splines_est = cbind( splines_est, exp(f_of_x(time_splines)) )											# store spline points for all dep. variables  	
			slopes_est = cbind( slopes_est, f_of_x(time_splines, deriv = 1) * exp(f_of_x(time_splines)) ) 						# store the slopes of the data points for all dep. variables	
			d2X_dt2_est = cbind( d2X_dt2_est, (f_of_x(time_splines, deriv = 2) * exp(f_of_x(time_splines))^2 - (f_of_x(time_splines, deriv = 1) * exp(f_of_x(time_splines)))^2 ) / exp(f_of_x(time_splines)) ) 		
			}
		}

     if (draw == TRUE)
            {
		par(mfrow=c(round(sqrt(numVars),0)+1,round(sqrt(numVars),0)+1)) 
		for (i in 1:numVars)                                                                            			# cycling thought the dependent vars
	      	{	
            	plot(dataset[,1],dataset[,i+1],pch=20,col="grey",ylab=colnames(data_vars)[i])                               # plot the dependent variable data

            	slopeXsize = tail(time_splines,1)*.025                                             					# find a 2.5% time interval to draw the slopes
			if (data1_spline2 == 1)			      											# draw slopes
            		{segments(x0 = time_splines - slopeXsize, x1 = time_splines + slopeXsize, y0 = dataset[,i+1] - slopeXsize * slopes_est[,i+1], y1 = dataset[,i+1] + slopeXsize * slopes_est[,i+1],col='lightgreen')} else
				{segments(x0 = time_splines - slopeXsize, x1 = time_splines + slopeXsize, y0 = splines_est[,i+1] - slopeXsize * slopes_est[,i+1], y1 = splines_est[,i+1] + slopeXsize * slopes_est[,i+1],col='lightgreen')}

            	points(splines_est[,1],splines_est[,i+1],type='l',col="darkgreen",lwd=3)                      			# plot the spline	
			points(dataset[,1],dataset[,i+1],pch=20,col="grey")

			if ( round(sqrt(numVars),0)+1 > 4 )
				{
				if (i%%16 == 0)
					{
            			windows()                                                       						# creating a new plot window, one for every dependent variable
					par(mfrow=c(4,4))
					}
				}
 			}
		}

    	return( list (
			data = dataset,
			splines_est = splines_est,			
		   	slopes_est = slopes_est,
			d2X_dt2_est = d2X_dt2_est,
			df = df
			) )
	}

##### Smoother
######################################################################################################################################



######################################################################################################################################
##### LV_par_finder
#
### function to estimate parameters for Lotka Volterra system.
# for theorical background see https://www.researchgate.net/publication/285256327_Parameter_estimation_in_canonical_biological_systems_models
#

# rm(list = ls())

# smooth_out = smoother_out
# supplied_slopes = NULL
# alg1_lm2 = 1
# data_sample_alg = 'random_sample'
# data_sample_lm = NULL
# givenParNames = NULL

LV_pars_finder = function (
				smooth_out,
				supplied_slopes = NULL,
				alg1_lm2 = 2, 
				data_sample_alg = 'random_sample',
				data_sample_lm = NULL,
				givenParNames = NULL
				)           
	{
	### function to estimate parameters for Lotka Volterra system.
	# smooth_out # please enter ypur favorite spline info
	# supplied_slopes # sample of slopes supplied by the used. If not NULL the slopes will not be calculated.
	# alg1_lm2 # 1 - the function will use a system of equations solver to find the parameters, 2 - it will use linear regretion
	# data_sample_alg # the points of the sample to be used if the solver of linear system of equations is selected. 'random_sample' - draw a random sample from the data or spline
	# data_sample_lm # the points of the sample to be used if the solver of linear system of equations is selected. When NULL it will use all sample points or spline points.
	# givenParNames # NULL - canonical par names will be asign # if the parameters have different names then the cannonical ones, you can enter then where

	dataset = smooth_out$data													# extract data from smooth_out
	data_time = dataset[,1]														# create a time var to use in splines and slopes
	numVars = dim(dataset)[2]-1													# set the number of dependent variable	
	data_vars = dataset[,2:dim(dataset)[2]]											# extract dependent vars data from smooth_out
	dataFrame = data.frame(data_time,data_vars)							  			# create data frame for the data to use in loess  
	time_splines = smooth_out$splines_est[,1]											# extracting time for splines
	splines_est = smooth_out$splines_est											# extract splines ( S(t) ) 	
	slopes_est = smooth_out$slopes_est						 						# extract slopes ( dS/dt )	
	if ( length(data_time) == length(splines_est[,1]) ) {data1_spline2 = 1} else {data1_spline2 = 2}	# see if we are using the data or spline extended time
	if (!is.null(supplied_slopes)) {slopes_est = supplied_slopes}							# if we have slopes, use it

	leftSide = vector()                                                                             	# creating the left side vector that will house the slopes / var value

    	if (alg1_lm2 == 1)
        	{
		# data sample for algebraic solution
		if (is.null(data_sample_alg)) 
			{return("Unspecified data_sample_alg")} else                                  		# if sample not defined, I define it
			{
			if (data_sample_alg[1] == 'random_sample')								# if the user wants a random sample 
				{
				data_sample_alg = sort( sample(1:length(time_splines),numVars+1) )			# take a random sample from the data or splines
				} else 
				{
				if (data1_spline2 == 2)
					{
					data_sample_alg = which(time_splines %in% data_time[data_sample_alg]) 
					}
				}
			}

		solution = vector()													# vector to store the solution
    		if (data1_spline2 == 1)
			{
			print(cbind(
					sample = data_sample_alg,
					t=time_splines[data_sample_alg],
					data_vars[data_sample_alg,]
					))
 			rightSide = cbind(rep(1,numVars+1),data_vars[data_sample_alg,])                           # build a small matrix with the values of the dependent variables values
   			for (i in 1:numVars)                                                                      # cycling thought the dependent vars
      			{
				leftSide = slopes_est[data_sample_alg,i+1]/data_vars[data_sample_alg,i]  		# create the left side of the system of equations with slopes / vars values
				cat('Left side ',i,"\n",leftSide,"\n")
				solution = c(solution, solve(rightSide, leftSide)) 						# solve the system of equations to find the parameter values
				}                                          		
			}
    		if (data1_spline2 == 2)
			{
			print(cbind(
					sample = data_sample_alg,
					t=time_splines[data_sample_alg],
					splines_est[data_sample_alg,2:dim(splines_est)[2]]
					))
			rightSide = cbind(rep(1,numVars+1),splines_est[data_sample_alg,2:dim(splines_est)[2]])	# build a small matrix with the values of the dependent variables values
    			for (i in 1:numVars)                                                                      # cycling thought the dependent vars
      			{
				leftSide = slopes_est[data_sample_alg,i+1]/splines_est[data_sample_alg,i+1]		# create the left side of the system of equations with slopes / vars values
				cat('Left side ',i,"\n",leftSide,"\n")
				solution = c(solution, solve(rightSide, leftSide)) 						# solve the system of equations to find the parameter values
				}  
			}
		cat('Right side')
		print(rightSide)
		cat("Right side determinant is ",det(rightSide),"\n")									# check the determinant of the right side of the system if equations to see if it is solvable
		cat("Right side dimensions are ",dim(rightSide)," and the rank is ",qr(rightSide)$rank,"\n","\n")	# calculate the rank of the right side of the system if equations to see if it is solvable		
        	}
    	
    	if (alg1_lm2 == 2) 
		{
		solution = vector()                                                                         	# creating vector to store the linear regretion estimates
		if (data1_spline2 == 1) 												# if we are using data
			{
			if (is.null(data_sample_lm)) {data_sample_lm = 1:length(data_time)}				# if data_sample_lm is NULL use all points in the sample
			for (i in 1:numVars)                                                                      # cycling thought the dependent vars
      			{
				leftside = slopes_est[data_sample_lm,i+1]/data_vars[data_sample_lm,i]			# create the leftside of the system of equations with slopes and equation var - slope/eqVar
				rightside = data_vars[data_sample_lm,] 								# create the rigthside of the system of equations with all var values 
        			solution = c(solution,lm(leftside~rightside)$coef)         					# linear regretion estimates with spline points
				}
			}
		if (data1_spline2 == 2) 												# if we are using the spline values
			{
			if (is.null(data_sample_lm)) {data_sample_lm = 1:length(time_splines)}				# if data_sample_lm is NULL use all points in the spline
			for (i in 1:numVars)                                                                      # cycling thought the dependent vars
      			{
				leftside = slopes_est[data_sample_lm,i+1]/splines_est[data_sample_lm,i+1]		# create the leftside of the system of equations with slopes and equation splin - slope/eqSpline
				rightside = splines_est[data_sample_lm,2:dim(splines_est)[2]] 				# create the rigthside of the system of equations with the spline values for all vars 
        			solution = c(solution,lm(leftside~rightside)$coef)         					# linear regretion estimates with spline points
				}
			}
		}  

	if ( !is.null(givenParNames) & length(givenParNames) == length(solution) )					# are given names for the parameters available and as long as the solutions?
		{ names(solution)=givenParNames } else										# if yes, use them
		{
		canonicalParNames = vector() 												# if no, construct cannonical pars names
		for ( i in 1:numVars )
			{
			canonicalParNames = c(canonicalParNames,paste('a',i,sep=''))					# a's
			for ( ii in 1:numVars )
				{
				canonicalParNames = c(canonicalParNames,paste('b',i,ii,sep=''))				# b's	
				}
			}
		names(solution)=canonicalParNames											# use cannonical par names
		} 
	return(solution)															# return pars_est
    	}

##### LV_par_finder
######################################################################################################################################



######################################################################################################################################
##### AlgPoinFinder

AlgPointFind = function (
				smoother_out_APF,
				dif_equations = Equations,
				matrixEq = TRUE,
				randomSearchTries = NULL,
				divideByMean = FALSE
				)
	{

	##########################################################
	### Function to find the best combination of datapoints for the albegraic method
	# smoother_out_APF # smoother results
	# dif_equations # stuff needed to run the LV solver
	# matrixEq # FALSE - parameters not in matrix form # TRUE - parameters are in matrix form
	# randomSearchTries # NULL - will check all posibilities # number - will check the given number of possibilities
	# divideByMean # FALSE - do nothing # TRUE - divide the errors by the mean of the dep var. This will balance the SSEs of the different dep. vars. if their values are very different. 

   	nVars = dim(smoother_out_APF$data)[2]-1																		# getting the number of dependent variables							

	#install.packages("gtools")
	library(gtools)
	pointComb = combinations(n=length(smoother_out_APF$data[,1]), r=nVars+1, v=1:dim(smoother_out_APF$data)[1], set=TRUE, repeats.allowed=FALSE)  	# calculating the different combinations available for the point sample
	if ( !is.null(randomSearchTries) ) { pointComb = pointComb[sample(1:dim(pointComb)[1],size = randomSearchTries),] }					# if we are doing random search, choose the random combinations to try

	global_store = store1 = store2 = c(rep(NA,nVars+1),10^20)															# create the stores for the results	

	for (ii in 1:dim(pointComb)[1])																			# cycle all point combinations
		tryCatch({
			( parEst_algPointFind = LV_pars_finder(
							smooth_out = smoother_out_APF,
							alg1_lm2 = 1, 
							data_sample_alg = pointComb[ii,]
							) )
			if (matrixEq==TRUE)
				{
				# formating the parEst_algPointFind to matrix form - use if system is in matrix form
				estPars_mat_temp = matrix(parEst_algPointFind,nrow=(-1+sqrt(1+4*length(parEst_algPointFind)))/2,byrow=TRUE)
				estPars_mat_algPointFind = list(		
									a = estPars_mat_temp[,1],
									B = estPars_mat_temp[,2:dim(estPars_mat_temp)[2]]
									)
				parEst_algPointFind = estPars_mat_algPointFind
				}

			state_algPointFind = unlist(smoother_out_APF$splines_est[1,2:(nVars+1)])
			names(state_algPointFind) = colnames(smoother_out_APF$data[,2:(nVars+1)])

			# var estimates
			out_est = NULL																								# set out_est to NULL. For previous sucessfull runs are not used when solve LV fails 
			out_est = try( solveLV(times = smoother_out_APF$data[,1], initState = state_algPointFind, pars = parEst_algPointFind, equations = dif_equations),TRUE)	# try the numerical solver	
	
			if (class(out_est)!="try-error") 																					# if try does not detect an error (may create warnnings)	
				{if (dim(out_est)[1] == length(smoother_out_APF$data[,1]))																# if out_est is the same length as the data (important to calculate errors)	
					{
					if (sum(is.nan(out_est))==0)																				# if there is no NAs
						{
						if (nVars==1) 
							{varMeans = mean(smoother_out_APF$data[,2:(nVars+1)],na.rm=TRUE)} else 
							{varMeans = apply(X=smoother_out_APF$data[,2:(nVars+1)],MARGIN=2,FUN=function(data1){mean(data1,na.rm=TRUE)})}				# get the means of each variable
						varMaenMatrix = matrix(rep(1,dim(smoother_out_APF$data)[1]*nVars),ncol=nVars)											# create a unitary matrix with the same dim as the data
						if ( divideByMean == TRUE ) { for (iv in 1:dim(smoother_out_APF$data)[1]) {varMaenMatrix[iv,] = varMeans} }						# if we are dividing by the means we will populate the matrix with the means of the dep. vars. matrix with the var means repeated to divide the errors so that a high value var does not dominate the errors

						error1 = sum(((smoother_out_APF$data[,2:(nVars+1)]-out_est[,2:(nVars+1)])/varMaenMatrix)^2)								# calculate the error agianst the data
						if ( error1 < store1[nVars+2] ) {store1 = c(pointComb[ii,],error1)}												# if the latest error is the smallest store it as the best 
						error2 = sum(((smoother_out_APF$splines_est[which(smoother_out_APF$splines_est[,1] %in% smoother_out_APF$data[,1]),2:(nVars+1)]-out_est[,2:(nVars+1)])/varMaenMatrix)^2)	# calculate the error agianst the splines
						if ( error2 < store2[nVars+2] ) {store2 = c(pointComb[ii,],error2) }												# if the latest error is the smallest store it as the best 
						global_store = rbind(global_store,c(pointComb[ii,],error1))														# store the best errors in the global store
						}
					}
				}
		# print results and percentage of work done
		print( paste(
				paste(round(ii/dim(pointComb)[1]*100,3),'% complete || DF ='),
				paste(smoother_out_APF$df,collapse = " " ),
				paste(' || '),
				paste(round(store1,3),collapse = " " ),
 				paste(' || '),
				paste(round(store2,3),collapse = " " )
				))
		flush.console()
		})
	global_store ->> globalStore
	return( rbind(store1,store2) )
	}

##### AlgPoinFinder
######################################################################################################################################



#############################################################################
### fastMAR

fastMAR = function(
			data_vars, 
			log_transform = FALSE, 
			demeaned = FALSE, 
			abstol_ = 0.5, 
			maxit_ = 100,
			includeNoise = TRUE,
			estInit = FALSE,
			tSteps_ = NULL,
			nsim_ = 1
			) 
	{
    	### function to prepare data for MARSS ###
    	# data_vars # abundances of the dependent variables by column. the columns should by named.
    	# no data_time because I do not know if I can feed MARSS data with NAs 
	# log_transform # FALSE - use original values of the dep. variables. # TRUE use ln-values of the dep. variables.
	# demeaned # FALSE = use values as they are # TRUE = removes the mean before estimation, it returns the mean for the data estimates
	# abstol_ # The logLik.(iter-1)-logLik.(iter) convergence tolerance for the maximization routine. To meet convergence both the abstol and slope tests must be passed. 
	# maxit_ # Maximum number of iterations to be used in the maximization routine (if needed by method) (positive integer). 
	# noise # TRUE - include noise in the  time series estimation process # FALSE - do not include noise
	# estInit # FALSE - first values of the data set used as initial conditions # TRUE - initial conditions are estimated.


	if ( log_transform == TRUE ) {data_vars = log(data_vars)}

	n = dim(data_vars)[1]

	if (is.null(n)) 
		{
		numVar = 1 
		n = length(data_vars) 
		data_vars_toUSE = data_vars

		data_vars_means = mean(data_vars)
		data_vars_means_mat = rep(data_vars_means,n)
		data_vars_demeaned = data_vars - data_vars_means_mat
		if ( demeaned == TRUE ) {data_vars_toUSE = data_vars_demeaned} 

		### fitting the MAR model.
		model.gen = list(
	     		B = "unconstrained", 		# matrix(list("b11","b12","b21","b22"),nrow=2,ncol=2,byrow=T),  		## Specify the "Full model" matrix for the MAR community matrix
           		U = "unequal", 		# matrix(list("a1","a2"),nrow=2,ncol=1), 							## state eq. intercept
           		Q = "diagonal and unequal",	# matrix(list("Q",0,0,"Q"),nrow=2,ncol=2),   					## state eq. var-cov matrix
			Z = "identity",               									  	  			## conversion matrice from space to state
           		A = "zero",     			# matrix(list(0,0),nrow=2,ncol=1), 							## space eq. intercept
           		R = "zero"	                          											## space eq. var-cov matrix
			)

		if (estInit == FALSE) { model.gen = c(model.gen,list( x0 = matrix(list(data_vars_toUSE[1])[[1]],nrow=numVar,ncol=1) )) }

		MARvalEst1 = x_t = data_vars_toUSE[1]

		} else
		{
		numVar = dim(data_vars)[2]		
		data_vars_toUSE = data_vars

		data_vars_means = rep(NA,numVar)
		data_vars_means_mat = matrix(rep(NA,n*numVar),ncol=numVar)
		for (i in 1:numVar)
			{
			data_vars_means[i] = mean(data_vars[,i],na.rm=TRUE)
			data_vars_means_mat[,i] = rep(data_vars_means[i],n)
			}
		data_vars_demeaned = data_vars - data_vars_means_mat
		if ( demeaned == TRUE ) {data_vars_toUSE = data_vars_demeaned} 		

		data_vars_toUSE = t(data_vars_toUSE)

		### fitting the MAR model.
		model.gen = list(
	     		B = "unconstrained", 		# matrix(list("b11","b12","b21","b22"),nrow=2,ncol=2,byrow=T),  		## Specify the "Full model" matrix for the MAR community matrix
           		U = "unequal",		# matrix(list("a1","a2"),nrow=2,ncol=1), 							## state eq. intercept
           		Q = "diagonal and unequal",	# matrix(list("Q",0,0,"Q"),nrow=2,ncol=2),   					## state eq. var-cov matrix
			Z = "identity",               									  	  			## conversion matrice from space to state
           		A = "zero",     			# matrix(list(0,0),nrow=2,ncol=1), 							## space eq. intercept
           		R = "zero"	                          											## space eq. var-cov matrix
           		#x0 = matrix(list(data_vars_toUSE[,1])[[1]],nrow=numVar,ncol=1)							## mean of the distrib from which starting values are drawn
			)

		if (estInit == FALSE) { model.gen = c(model.gen,list( x0 = matrix(list(data_vars_toUSE[,1])[[1]],nrow=numVar,ncol=1) )) }

		MARvalEst1 = x_t = data_vars_toUSE[,1]
		}

	library(MARSS)
	library(MASS)
	MARSS_output1 = MARSS(data_vars_toUSE ,model=model.gen,control=list(conv.test.slope.tol=0.01,abstol=abstol_,maxit=maxit_,allow.degen=TRUE),form="marxss",silent=F) 	# abstol=0.0001,maxit=10000
	try(m.ci1 <- MARSSparamCIs(MARSS_output1,method="hessian",alpha=0.05,nboot=100),silent=T)
		
	MARparEst1 = cbind(matrix(MARSS_output1$par$B, ncol=numVar),MARSS_output1$par$U,MARSS_output1$par$Q)

	if (estInit == TRUE) { MARvalEst1 = x_t = as.vector(MARSS_output1$par$x0) }

	for (i in 2:n)
		{
		x_tplus1 = MARparEst1[,1:numVar] %*% x_t + MARparEst1[,numVar+1] + includeNoise*rnorm(n=numVar, mean = rep(0,numVar), sd = sqrt(MARparEst1[,(numVar+2)]))		#
		MARvalEst1 = rbind(MARvalEst1,t(x_tplus1))
		x_t = x_tplus1	
		}

	if ( is.null(tSteps_) ) { tSteps_ = n }
	MARSSsim1 = MARSSsimulate(object = MARSS_output1, tSteps = tSteps_, nsim = nsim_)$sim.data
	sim1dim = dim(MARSSsim1)

	remeaner = matrix(rep(0,numVar*n),ncol=numVar)
	remeaner_sim = array(rep(0,sim1dim[1]*sim1dim[2]*sim1dim[3]),dim=sim1dim)
	if ( demeaned == TRUE ) 
		{ 
		remeaner = data_vars_means_mat 
		for (iv in 1:nsim_) { remeaner_sim[,,iv] = data_vars_means_mat }
		}

	if ( log_transform == TRUE ) 
		{
		return( list(
			MARSSsim_data = exp(MARSSsim1 + remeaner_sim),
			MARSS_output = MARSS_output1,
			MARSSci_output = m.ci1,
			MAR_TS_Est = exp(MARvalEst1 + remeaner)
			))
		} else
		{	
		return( list(
			MARSSsim_data = MARSSsim1 + remeaner_sim,
			MARSS_output = MARSS_output1,
			MARSSci_output = m.ci1,
			MAR_TS_Est = MARvalEst1 + remeaner
			))
		}
	}


meanMaker = function(MARsim)
	{
	sim1 = MARsim$MARSSsim_data
	dim(sim1)

	if (dim(sim1)[1]==1)
		{
		meanStore = MARsim$MARSSsim_data[,,1]
		for (ii in 1:dim(sim1)[1])
			{
			for (i in 1:dim(sim1)[2])
				{
				meanStore[i] = mean(sim1[ii,i,])
				}
			}
		} else
		{
		meanStore = MARsim$MARSSsim_data[,,1]
		for (ii in 1:dim(sim1)[1])
			{
			for (i in 1:dim(sim1)[2])
				{
				meanStore[ii,i] = mean(sim1[ii,i,])
				}
			}
		}

	return(meanStore)
	}

### fastMAR
#############################################################################



###############################################################################                                                                  

fourColPlot = function(fakedata,out_est=NULL,mainLab = NULL,yLab = NULL)
	{
	# Function to create a four color plot
	varColor_light = c('lightblue','orange','lightgreen','grey')
	varColor = c('blue','orange3','darkgreen','black')
	if (!is.null(out_est)) { xlimMax = max(dim(out_est)[1],dim(fakedata)[1]) } else { xlimMax = dim(fakedata)[1] }
	plot(1,1,col='white',xlim=c(0,xlimMax),ylim=c(0,3),main=mainLab,xlab='t',ylab=yLab,cex.lab=1.5,cex.axis=1.5)
	for  ( i in 2:5 )
		{
		points(fakedata[,1],fakedata[,i],pch=1,col=varColor[i-1])
		if (!is.null(out_est)) { points(out_est[,1],out_est[,i],type='l',col=varColor_light[i-1],lwd=3) }
		}
	# legend(0,3,legend = c('data','estimate','x1','x2','x3','x4'),lty = c(0,1,0,0,0,0),pch = c(16,NA,16,16,16,16),col = c('lightgrey','lightgrey','lightblue','orange','lightgreen','grey'),bty = "n",cex=.6)
	# legend(0,3,legend = c('data','estimate','x1','x2','x3','x4'),lty = c(0,1,1,1,1,1),pch = c(16,NA,NA,NA,NA,NA),col = c('grey','grey','blue','orange3','green','black'),bty = "n",cex=.6)
	}

varColor_light = c('lightblue','orange','lightgreen','grey')
varColor = c('blue','orange3','darkgreen','black')
varNames = c('X1','X2','X3','X4')

# fourColPlot(fakedata = fakedata2, out_est = out_est2)
###########################################################################################################################################################



######################################################################################################################################
##### system of dif. eq.

#install.packages("deSolve")
library(deSolve)

Format_pars = function(truePars)
	{
	truePars_temp = matrix(truePars,nrow=(-1+sqrt(1+4*length(truePars)))/2,byrow=TRUE)
	truePars_mat = list(		
				a = truePars_temp[,1],
				B = truePars_temp[,2:dim(truePars_temp)[2]]
				)
	return(truePars_mat)
	}

Equations <- function(t, x, pars) 
        { 
        ### returns rate of change
        # t = model's time structure
        # initState = model initial state
        # pars = model parameters 

        with(as.list(c(x, pars)), 
            {
		
		eq = ( a * x + x * (B %*% x) )

		dx = eq

            return(list(c(dx),count=c(eq)))
            })
        }

solveLV <- function(times = t, initState, pars, equations = Equations) 
    {
     ### ode solves the model by integration.
     # pars = model parameters
     # equations = model equations
     # initState = model initial state
     # times = model time structure

    return(ode(times = times, y = initState, parms = pars, func = equations))
    }

##### system of dif. eq.
######################################################################################################################################



######################################################################################################################################
##### converging to a stable steady state

################################
### time, initState and truePars  

t = seq(1,500,1)

state1=c(
	x1 = 1.2,
	x2 = .3,
	x3 = 2,
	x4 = .001
	) 

truePars1 = c(		
   		a1 = 0.044, 
  		b11 = -0.08,
		b12 = 0.02,
		b13 = 0.08,	
		b14 = 0,
		a2 = 0.216,
		b21 = -.04,
		b22 = -.08, 
		b23 = .04,
		b24 = 0,
		a3 = 0.116,
		b31 = -.16,
		b32 = 0.16,
		b33 = -0.08,
		b34 = 0,
		# a3 = 0.276,
		# b31 = -0.39,
		# b32 = 0.241,
		# b33 = -0.83,
		# b34 = 0,
		a4 = .2,
		b41 = 0,
		b42 = 0,
		b43 = 0,
		b44 = -.1
            )

truePars1_mat = Format_pars(truePars = truePars1)

### time, initState and truePars 
################################

fakedata_STST_500 = solveLV(times = t, initState = state1, pars = truePars1_mat, equations = Equations)
# tail(out) # edit(edit) # View(out)
fakedata_STST = fakedata_STST_500[1:100,1:5]


#########################
### slope solution - stst

smoother_out = Smoother(
				dataset = fakedata_STST,
				draw = TRUE,
				data1_spline2 = 2, 
				smooth_type = 1,
				df = c(100,100,100,100),	
				log_spline = TRUE
				)
( est1 = LV_pars_finder(
				smooth_out = smoother_out,
				alg1_lm2 = 1, 
				data_sample_alg = c(5,10,20,30,50)		#c(3, 7, 10, 14, 22)		
				) )

plot(1,1,col='white')
legend(1,1,legend = c('data','splines','slopes'),pch = c(20,NA,NA),lty = c(0,1,1),col = c('grey','green','lightgreen'),bty = "n")
cbind(est1)
est1_mat = Format_pars(truePars = est1)
cbind(truePars1,est1,absDif = abs(truePars1-est1))
splineStart = smoother_out$splines_est[1,2:5];names(splineStart) = c('x1','x2','x3','x4')
out_est1 = solveLV(times = t, pars = est1_mat, initState = splineStart, equations = Equations)

par(mfrow = c(2,2))
for (i in 2:5)
	{
	plot(fakedata_STST[,1],fakedata_STST[,i],col=varColor[i-1],ylim=c(0,3),xlab='t',ylab=varNames[i-1],cex.lab=1.5,cex.axis=1.5)
	points(out_est1[,1],out_est1[,i],type='l',col=varColor_light[i-1],lwd=3)
	}

### slope solution - stst
#########################


#######################
### MAR artData1 - stst

MARest_orig =  fastMAR(
			data_vars = fakedata_STST[,2:5],
			log_transform = FALSE,
			demeaned = FALSE,
			abstol_= 0.01,
			maxit_ = 1000,
			includeNoise = TRUE,
			estInit = FALSE,
			tSteps_ = NULL,
			nsim_ = 100
			) 
MARest_orig$MARSSci_output

numVar=4
MARparEst1 = matrix(MARest_orig$MARSS_output$coef, ncol=numVar+2)
MARvalEst1 = x_t = fakedata_STST[1,2:5]
for (i in 2:500)
	{
	x_tplus1 = MARparEst1[,1:numVar] %*% x_t + MARparEst1[,numVar+1] + rnorm(n=numVar, mean = rep(0,numVar), sd = MARparEst1[,(numVar+2)])		#
	MARvalEst1 = rbind(MARvalEst1,t(x_tplus1))
	x_t = x_tplus1	
	}
MARest500_orig = MARvalEst1

MARest_log =  fastMAR(
			data_vars = fakedata_STST[,2:5],
			log_transform = TRUE,
			demeaned = FALSE,
			abstol_= 0.01,
			maxit_ = 1000,
			includeNoise = TRUE
			) 
MARest_log$MARSSci_output

numVar=4
MARparEst1 = matrix(MARest_log$MARSS_output$coef, ncol=numVar+2)
MARvalEst1 = x_t = log(fakedata_STST[1,2:5])
for (i in 2:500)
	{
	x_tplus1 = MARparEst1[,1:numVar] %*% x_t + MARparEst1[,numVar+1] + rnorm(n=numVar, mean = rep(0,numVar), sd = MARparEst1[,(numVar+2)])		#
	MARvalEst1 = rbind(MARvalEst1,t(x_tplus1))
	x_t = x_tplus1	
	}
MARest500_log = exp(MARvalEst1)

par(mfrow=c(2,2))
for ( i in 2:5 )
	{
	plot(fakedata_STST[,1],fakedata_STST[,i],col=varColor[i-1],xlim=c(0,500),ylim=c(0,3),xlab='t',ylab=varNames[i-1],cex.lab=1.5,cex.axis=1.5)
	points(MARest500_orig[,i-1] ,type='l',col=varColor_light[i-1],lwd=3)
	points(MARest500_log[,i-1] ,type='l',lty=2,col=varColor_light[i-1],lwd=3)
	#if (i==1) {legend(40,1.8,legend=c('original data','noisy data','MAR no transform','MAR log transform'),pch=c(NA,20,NA,NA),lty=c(1,NA,1,2),col=c('lightgrey',varColor[i],varColor_light[i],varColor_light[i]),bty = "n",cex=1,lwd=c(3,NA,1,1))}
	}

### MAR artData1 - stst
#######################

# SSE
SSE_lab = c(
		sum((fakedata_STST_500[,2:5] - out_est1[,2:5])^2),
		sum((fakedata_STST_500[,2:5] - MARest500_orig)^2),
		sum((fakedata_STST_500[,2:5] - MARest500_log)^2)
		)

# setwd("C:\\Users\\doliv\\OneDrive\\2_2019_America\\2020\\20200123_MAR\\Camparison_LV_MAR\\8_methods_of_ecology_and_evolution_2sub\\paper_figures")  
# tiff("fig5a.tiff", height = 8, width = 30, units = 'cm',compression = "lzw", res = 300)
par(mfrow=c(1,4))
fourColPlot(fakedata = fakedata_STST_500, out_est = NULL,mainLab='data',yLab='steady state')
fourColPlot(fakedata = fakedata_STST, out_est = out_est1[,1:5],mainLab='ALVI-MI',yLab = paste('SSE =',round(SSE_lab[1],6)))
fourColPlot(fakedata = fakedata_STST, out_est = cbind(1:dim(MARest500_orig)[1],MARest500_orig),mainLab='MAR',yLab = paste('SSE =',round(SSE_lab[2],6)))
fourColPlot(fakedata = fakedata_STST, out_est = cbind(1:dim(MARest500_log)[1],MARest500_log),mainLab='MAR logTrans',yLab = paste('SSE =',round(SSE_lab[3],6)))
# dev.off()


##### converging to a stable steady state
######################################################################################################################################



######################################################################################################################################
##### damped oscillations

################################
### time, initState and truePars  

t = seq(1,500,1)

state1=c(
	x1 = .3,
	x2 = .3,
	x3 = .4,
	x4 = .6
	) 

truePars1 = c(		
   		a1 = .3, #.235,
  		b11 = -.3,
		b12 = -.27,
		b13 = -.6,	
		b14 = -.045,
		a2 = .4,
		b21 = .2,
		b22 = -.4, 
		b23 = -.4,
		b24 = -.6,
		a3 = .7,  # .8
		b31 = -2.38,
		b32 = .35,
		b33 = -2.8,
		b34 = .35,
		a4 = .6,
		b41 = -.96,
		b42 = -.24,
		b43 = -.96,
		b44 = -.6
            )

truePars1_mat = Format_pars(truePars = truePars1)

### time, initState and truePars 
################################

fakedata_dampOsc_500 = solveLV(times = t, initState = state1, pars = truePars1_mat, equations = Equations)
# tail(out) # edit(edit) # View(out)
fakedata_dampOsc = fakedata_dampOsc_500[1:100,1:5]


#########################
### slope solution - dampOsc

smoother_out = Smoother(
				dataset = fakedata_dampOsc,
				draw = TRUE,
				data1_spline2 = 2, 
				smooth_type = 1,
				df = c(100,100,100,100),	
				log_spline = TRUE
				)
( est1 = LV_pars_finder(
				smooth_out = smoother_out,
				alg1_lm2 = 1, 
				data_sample_alg = c(5,10,20,30,50)		#c(3, 7, 10, 14, 22)		
				) )

plot(1,1,col='white')
legend(1,1,legend = c('data','splines','slopes'),pch = c(20,NA,NA),lty = c(0,1,1),col = c('grey','green','lightgreen'),bty = "n")
cbind(est1)
est1_mat = Format_pars(truePars = est1)
cbind(truePars1,est1,absDif = abs(truePars1-est1))
splineStart = smoother_out$splines_est[1,2:5];names(splineStart) = c('x1','x2','x3','x4')
out_est1 = solveLV(times = t, pars = est1_mat, initState = splineStart, equations = Equations)

par(mfrow = c(2,2))
for (i in 2:5)
	{
	plot(fakedata_dampOsc[,1],fakedata_dampOsc[,i],col=varColor[i-1],xlim=c(0,500),ylim=c(0,3),xlab='t',ylab=varNames[i-1],cex.lab=1.5,cex.axis=1.5)
	points(out_est1[,1],out_est1[,i],type='l',col=varColor_light[i-1],lwd=3)
	}

### slope solution - dampOsc
#########################


#######################
### MAR artData1 - dampOsc

MARest_orig =  fastMAR(
			data_vars = fakedata_dampOsc[,2:5],
			log_transform = FALSE,
			demeaned = FALSE,
			abstol_= 0.01,
			maxit_ = 1000,
			includeNoise = TRUE
			) 
MARest_orig$MARSSci_output

numVar=4
MARparEst1 = matrix(MARest_orig$MARSS_output$coef, ncol=numVar+2)
MARvalEst1 = x_t = fakedata_dampOsc[1,2:5]
for (i in 2:500)
	{
	x_tplus1 = MARparEst1[,1:numVar] %*% x_t + MARparEst1[,numVar+1] + rnorm(n=numVar, mean = rep(0,numVar), sd = MARparEst1[,(numVar+2)])		#
	MARvalEst1 = rbind(MARvalEst1,t(x_tplus1))
	x_t = x_tplus1	
	}
MARest500_orig = MARvalEst1

MARest_log =  fastMAR(
			data_vars = fakedata_dampOsc[,2:5],
			log_transform = TRUE,
			demeaned = FALSE,
			abstol_= 0.01,
			maxit_ = 1000,
			includeNoise = TRUE
			) 
MARest_log$MARSSci_output

numVar=4
MARparEst1 = matrix(MARest_log$MARSS_output$coef, ncol=numVar+2)
MARvalEst1 = x_t = log(fakedata_dampOsc[1,2:5])
for (i in 2:500)
	{
	x_tplus1 = MARparEst1[,1:numVar] %*% x_t + MARparEst1[,numVar+1] + rnorm(n=numVar, mean = rep(0,numVar), sd = MARparEst1[,(numVar+2)])		#
	MARvalEst1 = rbind(MARvalEst1,t(x_tplus1))
	x_t = x_tplus1	
	}
MARest500_log = exp(MARvalEst1)

par(mfrow=c(2,2))
for ( i in 2:5 )
	{
	plot(fakedata_dampOsc[,1],fakedata_dampOsc[,i],col=varColor[i-1],ylim=c(0,3),xlab='t',ylab=varNames[i-1],cex.lab=1.5,cex.axis=1.5)
	points(MARest500_orig[,i-1] ,type='l',col=varColor_light[i-1],lwd=3)
	points(MARest500_log[,i-1] ,type='l',lty=2,col=varColor_light[i-1],lwd=3)
	#if (i==1) {legend(40,1.8,legend=c('original data','noisy data','MAR no transform','MAR log transform'),pch=c(NA,20,NA,NA),lty=c(1,NA,1,2),col=c('lightgrey',varColor[i],varColor_light[i],varColor_light[i]),bty = "n",cex=1,lwd=c(3,NA,1,1))}
	}

### MAR artData1 - dampOsc
#######################

# SSE
SSE_lab = c(
		sum((fakedata_dampOsc_500[,2:5] - out_est1[,2:5])^2),
		sum((fakedata_dampOsc_500[,2:5] - MARest500_log)^2),
		sum((fakedata_dampOsc_500[,2:5] - MARest500_orig)^2)
		)

# setwd("C:\\Users\\doliv\\OneDrive\\2_2019_America\\2020\\20200123_MAR\\Camparison_LV_MAR\\8_methods_of_ecology_and_evolution_2sub\\paper_figures")  
# tiff("fig5b.tiff", height = 8, width = 30, units = 'cm',compression = "lzw", res = 300)
par(mfrow=c(1,4))
fourColPlot(fakedata = fakedata_dampOsc_500, out_est = NULL,mainLab='data',yLab='damped oscillations')
fourColPlot(fakedata = fakedata_dampOsc, out_est = out_est1[,1:5],mainLab='ALVI-MI',yLab = paste('SSE =',round(SSE_lab[1],6)))
fourColPlot(fakedata = fakedata_dampOsc, out_est = cbind(1:dim(MARest500_orig)[1],MARest500_orig),mainLab='MAR',yLab = paste('SSE =',round(SSE_lab[2],6)))
fourColPlot(fakedata = fakedata_dampOsc, out_est = cbind(1:dim(MARest500_log)[1],MARest500_log),mainLab='MAR logTrans',yLab = paste('SSE =',round(SSE_lab[3],6)))
# dev.off()

##### damped oscillations
######################################################################################################################################



######################################################################################################################################
##### initially erratic oscillations converging to a limit cycle

################################
### time, initState and truePars  

t = seq(1,500,1)

state1=c(
	x1 = .3,
	x2 = .3,
	x3 = .4,
	x4 = .6
	) 

truePars1 = c(		
   		a1 = 1, 
  		b11 = -1,
		b12 = -1.09,
		b13 = -1.52,	
		b14 = 0,
		a2 = .72,
		b21 = 0,
		b22 = -.72, 
		b23 = -0.3168,
		b24 = -0.9792,
		a3 = 1.53,  
		b31 = -3.672,
		b32 = 0,
		b33 = -1.53,
		b34 = -0.7191,
		a4 = 1.27,
		b41 = -1.5367,
		b42 = -0.6477,
		b43 = -0.4445,
		b44 = -1.27
            )

truePars1_mat = Format_pars(truePars = truePars1)

### time, initState and truePars 
################################

fakedata_estOsc1_500 = solveLV(times = t, initState = state1, pars = truePars1_mat, equations = Equations)
# tail(out) # edit(edit) # View(out)
fakedata_estOsc1 = fakedata_estOsc1_500[1:100,1:5]


#########################
### slope solution - estOsc1

smoother_out = Smoother(
				dataset = fakedata_estOsc1,
				draw = TRUE,
				data1_spline2 = 2, 
				smooth_type = 1,
				df = c(100,100,100,100),	
				log_spline = TRUE
				)
( est1 = LV_pars_finder(
				smooth_out = smoother_out,
				alg1_lm2 = 1, 
				data_sample_alg = c(5,10,20,30,50)		#c(3, 7, 10, 14, 22)		
				) )

plot(1,1,col='white')
legend(1,1,legend = c('data','splines','slopes'),pch = c(20,NA,NA),lty = c(0,1,1),col = c('grey','green','lightgreen'),bty = "n")
cbind(est1)
est1_mat = Format_pars(truePars = est1)
cbind(truePars1,est1,absDif = abs(truePars1-est1))
splineStart = smoother_out$splines_est[1,2:5];names(splineStart) = c('x1','x2','x3','x4')
out_est1 = solveLV(times = t, pars = est1_mat, initState = splineStart, equations = Equations)

par(mfrow = c(2,2))
for (i in 2:5)
	{
	plot(fakedata_estOsc1[,1],fakedata_estOsc1[,i],col=varColor[i-1],xlim=c(0,500),ylim=c(0,3),xlab='t',ylab=varNames[i-1],cex.lab=1.5,cex.axis=1.5)
	points(out_est1[,1],out_est1[,i],type='l',col=varColor_light[i-1],lwd=3)
	}

### slope solution - estOsc1
#########################


#######################
### MAR artData1 - estOsc1

MARest_orig =  fastMAR(
			data_vars = fakedata_estOsc1[,2:5],
			log_transform = FALSE,
			demeaned = FALSE,
			abstol_= 0.01,
			maxit_ = 1000,
			includeNoise = TRUE
			) 
MARest_orig$MARSSci_output

numVar=4
MARparEst1 = matrix(MARest_orig$MARSS_output$coef, ncol=numVar+2)
MARvalEst1 = x_t = fakedata_dampOsc[1,2:5]
for (i in 2:500)
	{
	x_tplus1 = MARparEst1[,1:numVar] %*% x_t + MARparEst1[,numVar+1] + rnorm(n=numVar, mean = rep(0,numVar), sd = MARparEst1[,(numVar+2)])		#
	MARvalEst1 = rbind(MARvalEst1,t(x_tplus1))
	x_t = x_tplus1	
	}
MARest500_orig = MARvalEst1

MARest_log =  fastMAR(
			data_vars = fakedata_estOsc1[,2:5],
			log_transform = TRUE,
			demeaned = FALSE,
			abstol_= 0.01,
			maxit_ = 1000,
			includeNoise = TRUE
			) 
MARest_log$MARSSci_output

numVar=4
MARparEst1 = matrix(MARest_log$MARSS_output$coef, ncol=numVar+2)
MARvalEst1 = x_t = log(fakedata_dampOsc[1,2:5])
for (i in 2:500)
	{
	x_tplus1 = MARparEst1[,1:numVar] %*% x_t + MARparEst1[,numVar+1] + rnorm(n=numVar, mean = rep(0,numVar), sd = MARparEst1[,(numVar+2)])		#
	MARvalEst1 = rbind(MARvalEst1,t(x_tplus1))
	x_t = x_tplus1	
	}
MARest500_log = exp(MARvalEst1)

par(mfrow=c(2,2))
for ( i in 2:5 )
	{
	plot(fakedata_estOsc1[,1],fakedata_estOsc1[,i],col=varColor[i-1],ylim=c(0,3),xlab='t',ylab=varNames[i-1],cex.lab=1.5,cex.axis=1.5)
	points(MARest500_orig[,i-1] ,type='l',col=varColor_light[i-1],lwd=3)
	points(MARest500_log[,i-1] ,type='l',lty=2,col=varColor_light[i-1],lwd=3)
	#if (i==1) {legend(40,1.8,legend=c('original data','noisy data','MAR no transform','MAR log transform'),pch=c(NA,20,NA,NA),lty=c(1,NA,1,2),col=c('lightgrey',varColor[i],varColor_light[i],varColor_light[i]),bty = "n",cex=1,lwd=c(3,NA,1,1))}
	}

### MAR artData1 - estOsc1
#######################

# SSE
SSE_lab = c(
		sum((fakedata_estOsc1_500[,2:5] - out_est1[,2:5])^2),
		sum((fakedata_estOsc1_500[,2:5] - MARest500_orig)^2),
		sum((fakedata_estOsc1_500[,2:5] - MARest500_log)^2)
		)

# setwd("C:\\Users\\doliv\\OneDrive\\2_2019_America\\2020\\20200123_MAR\\Camparison_LV_MAR\\8_methods_of_ecology_and_evolution_2sub\\paper_figures")  
# tiff("fig5c.tiff", height = 8, width = 30, units = 'cm',compression = "lzw", res = 300)
par(mfrow=c(1,4))
fourColPlot(fakedata = fakedata_estOsc1_500, out_est = NULL,mainLab='data',yLab='converging to a limit cycle')
fourColPlot(fakedata = fakedata_estOsc1, out_est = out_est1[,1:5],mainLab='ALVI-MI',yLab = paste('SSE =',round(SSE_lab[1],6)))
fourColPlot(fakedata = fakedata_estOsc1, out_est = cbind(1:dim(MARest500_orig)[1],MARest500_orig),mainLab='MAR',yLab = paste('SSE =',round(SSE_lab[2],6)))
fourColPlot(fakedata = fakedata_estOsc1, out_est = cbind(1:dim(MARest500_log)[1],MARest500_log),mainLab='MAR logTrans',yLab = paste('SSE =',round(SSE_lab[3],6)))
# dev.off()

##### initially erratic oscillations converging to a limit cycle
######################################################################################################################################



######################################################################################################################################
##### sustained oscillations

################################
### time, initState and truePars  

t = seq(1,500,1)

state1=c(
	x1 = .3,
	x2 = .3,
	x3 = .4,
	x4 = .6
	) 

truePars1 = c(		
   		a1 = .3, #.235,
  		b11 = -.3,
		b12 = -.27,
		b13 = -.6,	
		b14 = -.045,
		a2 = .4,
		b21 = .2,
		b22 = -.4, 
		b23 = -.4,
		b24 = -.6,
		a3 = .7,  # .8
		b31 = -2.38,
		b32 = .35,
		b33 = -2.45,
		b34 = .35,
		a4 = .6,
		b41 = -.96,
		b42 = -.24,
		b43 = -.96,
		b44 = -.3
            )

truePars1_mat = Format_pars(truePars = truePars1)

### time, initState and truePars 
################################

fakedata_estOsc2_500 = solveLV(times = t, initState = state1, pars = truePars1_mat, equations = Equations)
# tail(out) # edit(edit) # View(out)
fakedata_estOsc2 = fakedata_estOsc2_500[1:100,1:5]


#########################
### slope solution - estOsc2

smoother_out = Smoother(
				dataset = fakedata_estOsc2,
				draw = TRUE,
				data1_spline2 = 2, 
				smooth_type = 1,
				df = c(100,100,100,100),	
				log_spline = TRUE
				)
( est1 = LV_pars_finder(
				smooth_out = smoother_out,
				alg1_lm2 = 1, 
				data_sample_alg = c(5,10,20,30,50)		#c(3, 7, 10, 14, 22)		
				) )

plot(1,1,col='white')
legend(1,1,legend = c('data','splines','slopes'),pch = c(20,NA,NA),lty = c(0,1,1),col = c('grey','green','lightgreen'),bty = "n")
cbind(est1)
est1_mat = Format_pars(truePars = est1)
cbind(truePars1,est1,absDif = abs(truePars1-est1))
splineStart = smoother_out$splines_est[1,2:5];names(splineStart) = c('x1','x2','x3','x4')
out_est1 = solveLV(times = t, pars = est1_mat, initState = splineStart, equations = Equations)

par(mfrow = c(2,2))
for (i in 2:5)
	{
	plot(fakedata_estOsc2[,1],fakedata_estOsc2[,i],col=varColor[i-1],xlim=c(0,500),ylim=c(0,3),xlab='t',ylab=varNames[i-1],cex.lab=1.5,cex.axis=1.5)
	points(out_est1[,1],out_est1[,i],type='l',col=varColor_light[i-1],lwd=3)
	}

### slope solution - estOsc2
#########################


#######################
### MAR artData1 - estOsc2

MARest_orig =  fastMAR(
			data_vars = fakedata_estOsc2[,2:5],
			log_transform = FALSE,
			demeaned = FALSE,
			abstol_= 0.01,
			maxit_ = 1000,
			includeNoise = TRUE
			) 
MARest_orig$MARSSci_output

numVar=4
MARparEst1 = matrix(MARest_orig$MARSS_output$coef, ncol=numVar+2)
MARvalEst1 = x_t = fakedata_dampOsc[1,2:5]
for (i in 2:500)
	{
	x_tplus1 = MARparEst1[,1:numVar] %*% x_t + MARparEst1[,numVar+1] + rnorm(n=numVar, mean = rep(0,numVar), sd = MARparEst1[,(numVar+2)])		#
	MARvalEst1 = rbind(MARvalEst1,t(x_tplus1))
	x_t = x_tplus1	
	}
MARest500_orig = MARvalEst1

MARest_log =  fastMAR(
			data_vars = fakedata_estOsc2[,2:5],
			log_transform = TRUE,
			demeaned = FALSE,
			abstol_= 0.01,
			maxit_ = 1000,
			includeNoise = TRUE
			) 
MARest_log$MARSSci_output

numVar=4
MARparEst1 = matrix(MARest_log$MARSS_output$coef, ncol=numVar+2)
MARvalEst1 = x_t = log(fakedata_dampOsc[1,2:5])
for (i in 2:500)
	{
	x_tplus1 = MARparEst1[,1:numVar] %*% x_t + MARparEst1[,numVar+1] + rnorm(n=numVar, mean = rep(0,numVar), sd = MARparEst1[,(numVar+2)])		#
	MARvalEst1 = rbind(MARvalEst1,t(x_tplus1))
	x_t = x_tplus1	
	}
MARest500_log = exp(MARvalEst1)

par(mfrow=c(2,2))
for ( i in 2:5 )
	{
	plot(fakedata_estOsc1[,1],fakedata_estOsc1[,i],col=varColor[i-1],ylim=c(0,3),xlab='t',ylab=varNames[i-1],cex.lab=1.5,cex.axis=1.5)
	points(MARest500_orig[,i-1] ,type='l',col=varColor_light[i-1],lwd=3)
	points(MARest500_log[,i-1] ,type='l',lty=2,col=varColor_light[i-1],lwd=3)
	#if (i==1) {legend(40,1.8,legend=c('original data','noisy data','MAR no transform','MAR log transform'),pch=c(NA,20,NA,NA),lty=c(1,NA,1,2),col=c('lightgrey',varColor[i],varColor_light[i],varColor_light[i]),bty = "n",cex=1,lwd=c(3,NA,1,1))}
	}

### MAR artData1 - estOsc2
#######################

# SSE
SSE_lab = c(
		sum((fakedata_estOsc2_500[,2:5] - out_est1[,2:5])^2),
		sum((fakedata_estOsc2_500[,2:5] - MARest500_orig)^2),
		sum((fakedata_estOsc2_500[,2:5] - MARest500_log)^2)
		)

# setwd("C:\\Users\\doliv\\OneDrive\\2_2019_America\\2020\\20200123_MAR\\Camparison_LV_MAR\\8_methods_of_ecology_and_evolution_2sub\\paper_figures")  
# tiff("fig5d.tiff", height = 8, width = 30, units = 'cm',compression = "lzw", res = 300)
par(mfrow=c(1,4))
fourColPlot(fakedata = fakedata_estOsc2_500, out_est = NULL,mainLab='data',yLab='sustained oscillations')
fourColPlot(fakedata = fakedata_estOsc2, out_est = out_est1[,1:5],mainLab='ALVI-MI',yLab = paste('SSE =',round(SSE_lab[1],6)))
fourColPlot(fakedata = fakedata_estOsc2, out_est = cbind(1:dim(MARest500_orig)[1],MARest500_orig),mainLab='MAR',yLab = paste('SSE =',round(SSE_lab[2],6)))
fourColPlot(fakedata = fakedata_estOsc2, out_est = cbind(1:dim(MARest500_log)[1],MARest500_log),mainLab='MAR logTrans',yLab = paste('SSE =',round(SSE_lab[3],6)))
# dev.off()

##### sustained oscillations
######################################################################################################################################



######################################################################################################################################
##### deterministic chaos

################################
### time, initState and truePars  

t = seq(1,500,1)

state1=c(
	x1 = .3,
	x2 = .3,
	x3 = .4,
	x4 = .6
	) 

truePars1 = c(		
   		a1 = 1, 
  		b11 = -1,
		b12 = -1.09,
		b13 = -1.52,	
		b14 = 0,
		a2 = 0.72,
		b21 = 0,
		b22 = -0.72, 
		b23 = -0.3168,
		b24 = -0.9792,
		a3 = 1.53,
		b31 = -3.5649,
		b32 = 0,
		b33 = -1.53,
		b34 = -0.7191,
		a4 = 1.27,
		b41 = -1.5367,
		b42 = -0.6477,
		b43 = -0.4445,
		b44 = -1.27
            )

truePars1_mat = Format_pars(truePars = truePars1)

### time, initState and truePars 
################################

fakedata_chaos_500 = solveLV(times = t, initState = state1, pars = truePars1_mat, equations = Equations)
# tail(out) # edit(edit) # View(out)
fakedata_chaos = fakedata_chaos_500[1:100,1:5]


#########################
### slope solution - chaos

smoother_out = Smoother(
				dataset = fakedata_chaos,
				draw = TRUE,
				data1_spline2 = 2, 
				smooth_type = 1,
				df = c(100,100,100,100),	
				log_spline = TRUE
				)
( est1 = LV_pars_finder(
				smooth_out = smoother_out,
				alg1_lm2 = 1, 
				data_sample_alg = c(10,20,30,50,80)		#c(3, 7, 10, 14, 22)		
				) )

plot(1,1,col='white')
legend(1,1,legend = c('data','splines','slopes'),pch = c(20,NA,NA),lty = c(0,1,1),col = c('grey','green','lightgreen'),bty = "n")
cbind(est1)
est1_mat = Format_pars(truePars = est1)
cbind(truePars1,est1,absDif = abs(truePars1-est1))
splineStart = smoother_out$splines_est[1,2:5];names(splineStart) = c('x1','x2','x3','x4')
out_est1 = solveLV(times = t, pars = est1_mat, initState = splineStart, equations = Equations)

par(mfrow = c(2,2))
for (i in 2:5)
	{
	plot(fakedata_chaos[,1],fakedata_chaos[,i],col=varColor[i-1],xlim=c(0,500),ylim=c(0,3),xlab='t',ylab=varNames[i-1],cex.lab=1.5,cex.axis=1.5)
	points(out_est1[,1],out_est1[,i],type='l',col=varColor_light[i-1],lwd=3)
	}

### slope solution - chaos
#########################


#######################
### MAR artData1 - chaos

MARest_orig =  fastMAR(
			data_vars = fakedata_chaos[,2:5],
			log_transform = FALSE,
			demeaned = FALSE,
			abstol_= 0.01,
			maxit_ = 1000,
			includeNoise = TRUE
			) 
MARest_orig$MARSSci_output

numVar=4
MARparEst1 = matrix(MARest_orig$MARSS_output$coef, ncol=numVar+2)
MARvalEst1 = x_t = fakedata_dampOsc[1,2:5]
for (i in 2:500)
	{
	x_tplus1 = MARparEst1[,1:numVar] %*% x_t + MARparEst1[,numVar+1] + rnorm(n=numVar, mean = rep(0,numVar), sd = MARparEst1[,(numVar+2)])		#
	MARvalEst1 = rbind(MARvalEst1,t(x_tplus1))
	x_t = x_tplus1	
	}
MARest500_orig = MARvalEst1

MARest_log =  fastMAR(
			data_vars = fakedata_chaos[,2:5],
			log_transform = TRUE,
			demeaned = FALSE,
			abstol_= 0.01,
			maxit_ = 1000,
			includeNoise = TRUE
			) 
MARest_log$MARSSci_output

numVar=4
MARparEst1 = matrix(MARest_log$MARSS_output$coef, ncol=numVar+2)
MARvalEst1 = x_t = log(fakedata_dampOsc[1,2:5])
for (i in 2:500)
	{
	x_tplus1 = MARparEst1[,1:numVar] %*% x_t + MARparEst1[,numVar+1] + rnorm(n=numVar, mean = rep(0,numVar), sd = MARparEst1[,(numVar+2)])		#
	MARvalEst1 = rbind(MARvalEst1,t(x_tplus1))
	x_t = x_tplus1	
	}
MARest500_log = exp(MARvalEst1)

par(mfrow=c(2,2))
for ( i in 2:5 )
	{
	plot(fakedata_chaos[,1],fakedata_chaos[,i],col=varColor[i-1],ylim=c(0,3),xlab='t',ylab=varNames[i-1],cex.lab=1.5,cex.axis=1.5)
	points(MARest500_orig[,i-1] ,type='l',col=varColor_light[i-1],lwd=3)
	points(MARest500_log[,i-1] ,type='l',lty=2,col=varColor_light[i-1],lwd=3)
	#if (i==1) {legend(40,1.8,legend=c('original data','noisy data','MAR no transform','MAR log transform'),pch=c(NA,20,NA,NA),lty=c(1,NA,1,2),col=c('lightgrey',varColor[i],varColor_light[i],varColor_light[i]),bty = "n",cex=1,lwd=c(3,NA,1,1))}
	}

### MAR artData1 - chaos
#######################

# SSE
SSE_lab = c(
		sum((fakedata_chaos_500[,2:5] - out_est1[,2:5])^2),
		sum((fakedata_chaos_500[,2:5] - MARest500_orig)^2),
		sum((fakedata_chaos_500[,2:5] - MARest500_log)^2)
		)

# setwd("C:\\Users\\doliv\\OneDrive\\2_2019_America\\2020\\20200123_MAR\\Camparison_LV_MAR\\8_methods_of_ecology_and_evolution_2sub\\paper_figures")  
# tiff("fig5e.tiff", height = 8, width = 30, units = 'cm',compression = "lzw", res = 300)
par(mfrow=c(1,4))
fourColPlot(fakedata = fakedata_chaos_500, out_est = NULL,mainLab='data',yLab='chaos 1')
fourColPlot(fakedata = fakedata_chaos, out_est = out_est1[,1:5],mainLab='ALVI-MI',yLab = paste('SSE =',round(SSE_lab[1],6)))
fourColPlot(fakedata = fakedata_chaos, out_est = cbind(1:dim(MARest500_orig)[1],MARest500_orig),mainLab='MAR',yLab = paste('SSE =',round(SSE_lab[2],6)))
fourColPlot(fakedata = fakedata_chaos, out_est = cbind(1:dim(MARest500_log)[1],MARest500_log),mainLab='MAR logTrans',yLab = paste('SSE =',round(SSE_lab[3],6)))
# dev.off()

##### deterministic chaos
######################################################################################################################################



######################################################################################################################################
##### chaotic oscillations

################################
### time, initState and truePars  

t = seq(1,500,1)

state1=c(
	x1 = .3,
	x2 = .3,
	x3 = .4,
	x4 = .6
	) 

truePars1 = c(		
   		a1 = .3, #.235,
  		b11 = -.3,
		b12 = -.27,
		b13 = -.6,	
		b14 = -.045,
		a2 = .4,
		b21 = .2,
		b22 = -.4, 
		b23 = -.4,
		b24 = -.6,
		a3 = .8,  #
		b31 = -2.38,
		b32 = .35,
		b33 = -2.45,
		b34 = .35,
		a4 = .6,
		b41 = -.96,
		b42 = -.24,
		b43 = -.96,
		b44 = -.3
            )

truePars1_mat = Format_pars(truePars = truePars1)

### time, initState and truePars 
################################

fakedata_chaosOsc_500 = solveLV(times = t, initState = state1, pars = truePars1_mat, equations = Equations)
# tail(out) # edit(edit) # View(out)
fakedata_chaosOsc = fakedata_chaosOsc_500[1:100,1:5]


#########################
### slope solution - chaosOsc

smoother_out = Smoother(
				dataset = fakedata_chaosOsc,
				draw = TRUE,
				data1_spline2 = 2, 
				smooth_type = 1,
				df = c(100,100,100,100),	
				log_spline = TRUE
				)
( est1 = LV_pars_finder(
				smooth_out = smoother_out,
				alg1_lm2 = 1, 
				data_sample_alg = c(4,6,10,15,35)		#c(3, 7, 10, 14, 22)		
				) )

plot(1,1,col='white')
legend(1,1,legend = c('data','splines','slopes'),pch = c(20,NA,NA),lty = c(0,1,1),col = c('grey','green','lightgreen'),bty = "n")
cbind(est1)
est1_mat = Format_pars(truePars = est1)
cbind(truePars1,est1,absDif = abs(truePars1-est1))
splineStart = smoother_out$splines_est[1,2:5];names(splineStart) = c('x1','x2','x3','x4')
out_est1 = solveLV(times = t, pars = est1_mat, initState = splineStart, equations = Equations)

par(mfrow = c(2,2))
for (i in 2:5)
	{
	plot(fakedata_chaosOsc[,1],fakedata_chaosOsc[,i],col=varColor[i-1],xlim=c(0,500),ylim=c(0,3),xlab='t',ylab=varNames[i-1],cex.lab=1.5,cex.axis=1.5)
	points(out_est1[,1],out_est1[,i],type='l',col=varColor_light[i-1],lwd=3)
	}

### slope solution - chaosOsc
#########################


#######################
### MAR artData1 - chaosOsc

MARest_orig =  fastMAR(
			data_vars = fakedata_chaosOsc[,2:5],
			log_transform = FALSE,
			demeaned = FALSE,
			abstol_= 0.01,
			maxit_ = 1000,
			includeNoise = TRUE
			) 
MARest_orig$MARSSci_output

numVar=4
MARparEst1 = matrix(MARest_orig$MARSS_output$coef, ncol=numVar+2)
MARvalEst1 = x_t = fakedata_dampOsc[1,2:5]
for (i in 2:500)
	{
	x_tplus1 = MARparEst1[,1:numVar] %*% x_t + MARparEst1[,numVar+1] + rnorm(n=numVar, mean = rep(0,numVar), sd = MARparEst1[,(numVar+2)])		#
	MARvalEst1 = rbind(MARvalEst1,t(x_tplus1))
	x_t = x_tplus1	
	}
MARest500_orig = MARvalEst1

MARest_log =  fastMAR(
			data_vars = fakedata_chaosOsc[,2:5],
			log_transform = TRUE,
			demeaned = FALSE,
			abstol_= 0.01,
			maxit_ = 1000,
			includeNoise = TRUE
			) 
MARest_log$MARSSci_output

numVar=4
MARparEst1 = matrix(MARest_log$MARSS_output$coef, ncol=numVar+2)
MARvalEst1 = x_t = log(fakedata_dampOsc[1,2:5])
for (i in 2:500)
	{
	x_tplus1 = MARparEst1[,1:numVar] %*% x_t + MARparEst1[,numVar+1] + rnorm(n=numVar, mean = rep(0,numVar), sd = MARparEst1[,(numVar+2)])		#
	MARvalEst1 = rbind(MARvalEst1,t(x_tplus1))
	x_t = x_tplus1	
	}
MARest_log500 = exp(MARvalEst1)

par(mfrow=c(2,2))
for ( i in 2:5 )
	{
	plot(fakedata_chaosOsc[,1],fakedata_chaosOsc[,i],col=varColor[i-1],ylim=c(0,3),xlab='t',ylab=varNames[i-1],cex.lab=1.5,cex.axis=1.5)
	points(MARest500_orig[,i-1] ,type='l',col=varColor_light[i-1],lwd=3)
	points(MARest500_log[,i-1] ,type='l',lty=2,col=varColor_light[i-1],lwd=3)
	#if (i==1) {legend(40,1.8,legend=c('original data','noisy data','MAR no transform','MAR log transform'),pch=c(NA,20,NA,NA),lty=c(1,NA,1,2),col=c('lightgrey',varColor[i],varColor_light[i],varColor_light[i]),bty = "n",cex=1,lwd=c(3,NA,1,1))}
	}

### MAR artData1 - chaosOsc
#######################

# SSE
SSE_lab = c(
		sum((fakedata_chaosOsc_500[,2:5] - out_est1[,2:5])^2),
		sum((fakedata_chaosOsc_500[,2:5] - MARest500_orig)^2),
		sum((fakedata_chaosOsc_500[,2:5] - MARest500_log)^2)
		)

# setwd("C:\\Users\\doliv\\OneDrive\\2_2019_America\\2020\\20200123_MAR\\Camparison_LV_MAR\\8_methods_of_ecology_and_evolution_2sub\\paper_figures")  
# tiff("fig5f.tiff", height = 8, width = 30, units = 'cm',compression = "lzw", res = 300)
par(mfrow=c(1,4))
fourColPlot(fakedata = fakedata_chaosOsc_500, out_est = NULL,mainLab='data',yLab='chaos 2')
fourColPlot(fakedata = fakedata_chaosOsc, out_est = out_est1[,1:5],mainLab='ALVI-MI',yLab = paste('SSE =',round(SSE_lab[1],6)))
fourColPlot(fakedata = fakedata_chaosOsc, out_est = cbind(1:dim(MARest500_orig)[1],MARest500_orig),mainLab='MAR',yLab = paste('SSE =',round(SSE_lab[2],6)))
fourColPlot(fakedata = fakedata_chaosOsc, out_est = cbind(1:dim(MARest500_log)[1],MARest500_log),mainLab='MAR logTrans',yLab = paste('SSE =',round(SSE_lab[3],6)))
legend(0,3,legend = c('data','estimate',expression(italic('X'[1])),expression(italic('X'[2])),expression(italic('X'[3])),expression(italic('X'[4]))),lty = c(0,1,0,0,0,0),pch = c(16,NA,16,16,16,16),col = c('lightgrey','lightgrey','lightblue','orange','lightgreen','grey'),bty = "n",cex=1.5)
legend(0,3,legend = c('data','estimate',expression(italic('X'[1])),expression(italic('X'[2])),expression(italic('X'[3])),expression(italic('X'[4]))),lty = c(0,1,1,1,1,1),pch = c(16,NA,NA,NA,NA,NA),col = c('grey','grey','blue','orange3','green','black'),bty = "n",cex=1.5)
# dev.off()

##### chaotic oscillations
######################################################################################################################################


