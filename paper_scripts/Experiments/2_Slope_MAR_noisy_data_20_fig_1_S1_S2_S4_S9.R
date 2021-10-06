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
	# polyDegree012 # the degree of the local polynomials to be used. It can ben 0, 1 or 2. Only used on LOESS.
	# aicc1_gvc2 # the criterion for automatic smoothing parameter selection: “aicc” denotes bias-corrected AIC criterion, “gcv” denotes generalized cross-validation. 
	# span # value from [0,1] to define the percentage of sample points in the moving interval
	# log_spline # TRUE uses log transform datapoints to create the splines. Splines will by in cartesian values.
	
	data_time = dataset[,1]																				# create a time var to use in splines and slopes
	numVars = dim(dataset)[2]-1																			# set the number of dependent variable	
	data_vars = matrix(rep(NA,dim(dataset)[1]*(dim(dataset)[2]-1)),ncol = numVars)										# unlisting the data
	for (iii in 1:(dim(dataset)[2]-1)) 
		{
		data_vars[,iii] = unlist(dataset[,iii+1])
		}
	colnames(data_vars) = colnames(dataset[,2:(dim(dataset)[2])])													# naming the data_vars colunms
	if (log_spline == TRUE) {data_vars = log(data_vars)}															# use log var values to calculate the splines if requested
	dataFrame = data.frame(t=dataset[,1],data_vars)							  									# create data frame for the data to use in loess  

	if ( data1_spline2 == 1 )
		{
		time_splines = data_time
		} else
		{
		time_splines = round(seq(head(data_time,1),tail(data_time,1),length.out = 10*length(data_time)),10)						# create a time sequence with a smaller resolution
		time_splines = unique(sort(c(data_time,time_splines)))
		}

	if ( !is.null(dataWeights) & length(dataWeights)==1)
		{ if ( dataWeights == 'heavyStart' ) { dataWeights = rep(1,dim(dataset)[1]);dataWeights[1]=100 } }						# set the weigths for when the initial values are 100 times heavier than the other data points

	if (is.null(df)) {df_temp = rep(NA,numVars)}

	# splines and LOESS (smoothing) 
    	splines_est = slopes_est = d2X_dt2_est = time_splines															# initiating the vars to store the results
	for (i in 1:numVars)                                                                            							# cycling thought the dependent vars
      	{
        	if (smooth_type == 1)
            	{
			if (is.null(df))
				{
				smoothSpline <- smooth.spline(data_time,data_vars[,i],w=dataWeights)
				df_temp[i] = smoothSpline$df
				} else 								# smoothing with estimated degrees of freedom 
				{smoothSpline <- smooth.spline(data_time,data_vars[,i],w=dataWeights,df=df[i])}                             		# smoothing with degrees of freedom (df)
            	smooth = predict(object = smoothSpline, x = time_splines, deriv = 0)                                 					# get the points of the fit of that linear model
            	f_of_x <- splinefun(smooth$x,smooth$y,method = splineMethod )  # "fmm" "periodic" "natural" "hyman" "monoH.FC"                                               				# create cubic spline for the linear model
            	}
	  	if (smooth_type == 2)
			{
  			loess1 = loess.as(dataFrame[,1], dataFrame[,i+1], 
						degree = polyDegree012, 
						criterion = c("aicc", "gcv")[aicc1_gvc2], 
						user.span = span, plot = F )														# create loess of the data points with span = .7
               	smooth = loess1$fit																		# store the loess fit to build the spline
                	f_of_x <- splinefun(dataFrame[,1],smooth)                                     								# create cubic spline for a dependent variable
			print(summary(loess1))
			}

		if (log_spline == FALSE)
			{		
			splines_est = cbind( splines_est, f_of_x(time_splines) )												# store spline points for all dep. variables  	
			slopes_est = cbind( slopes_est, f_of_x(time_splines, deriv = 1) ) 										# store the slopes of the data points for all dep. variables	
			d2X_dt2_est = cbind( d2X_dt2_est, f_of_x(time_splines, deriv = 2) ) 										# store the 2nd derivative estimates
			} else
			{
			splines_est = cbind( splines_est, exp(f_of_x(time_splines)) )											# store spline points for all dep. variables  	
			slopes_est = cbind( slopes_est, f_of_x(time_splines, deriv = 1) * exp(f_of_x(time_splines)) ) 						# store the slopes of the data points for all dep. variables when y is in log form. dlog(y)/dy =1/y * dy/dt <=> y * dlog(y)/dy = dy/dt	
			d2X_dt2_est = cbind( d2X_dt2_est, f_of_x(time_splines, deriv = 2) * exp(f_of_x(time_splines)) + f_of_x(time_splines, deriv = 1)^2 * exp(f_of_x(time_splines)) )  	# store the 2nd derivative estimates when y is in log form. d^2(y)/(dt)^2 = d^2(log(y))/(dt)^2 * y + (d(log(y))/dt)^2 * y 
			}
		}
	if (is.null(df)) {df = df_temp}

     	if (draw == TRUE)
            {
		par(mfrow=c(round(sqrt(numVars),0)+1,round(sqrt(numVars),0)+1)) 
		for (i in 1:numVars)                                                                            						# cycling thought the dependent vars
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



######################################################################################################################################
##### Smoother_for_MAR

Smoother_for_MAR = function (
					dataset,
					df = NULL,
					dataWeights = NULL,
					splineMethod = "fmm",
					draw = TRUE,
					log_spline = FALSE
					)
	{
	### function to estimate parameters for MAR system. The domain of the spline will be the natural numbers. 
	# dataset # the data. First colunm is the timepoints, remaining colunms the dependent variables. All colunms should be named.
	# df # degrees of freedom for smoothing splines 2, must have 1 < df <= n. If NULL, df = n.
	# dataWeights # NULL - all datapoints have the same weights # 'heavyStart' - initial values 100 times heavier that other datapoints # it acepts handcrafted weights
	# splineMethod # method to estimate splines # "fmm" "periodic" "natural" # monotonic splines "hyman" "monoH.FC"
	# draw # TRUE draw the estimated slopes for every data point. 
	# log_spline # TRUE uses log transform datapoints to create the splines. Splines will by in cartesian values.

	#install.packages("splines")
	library(splines)

	data_time = dataset[,1]														# create a time var to use in splines and slopes
	numVars = dim(dataset)[2]-1	
	data_vars = matrix(rep(NA,dim(dataset)[1]*(dim(dataset)[2]-1)),ncol = numVars)				# unlisting the data
	for (iii in 1:(dim(dataset)[2]-1)) 
		{
		data_vars[,iii] = unlist(dataset[,iii+1])
		}
	colnames(data_vars) = colnames(dataset[,2:(dim(dataset)[2])])
	if (log_spline == TRUE) {data_vars = log(data_vars)}									# use log var values to calculate the splines if requested
	
	if ( !is.null(dataWeights) & length(dataWeights)==1)
		{ if ( dataWeights == 'heavyStart' ) { dataWeights = rep(1,dim(dataset)[1]);dataWeights[1]=100 } }	# set the weigths for when the initial values are 100 times heavier than the other data points

	time_splines = seq(data_time[1],tail(data_time,1),1)

	splines_est = time_splines	
	for (i in 1:numVars)                                                                            	# cycling thought the dependent vars
	     	{
		smoothSpline <- smooth.spline(data_time,data_vars[,i],w=dataWeights,df=df[i])                            				# smoothing with degrees of freedom (df)
		smooth = predict(object = smoothSpline, x = seq(data_time[1],tail(data_time,1),1), deriv = 0)                                 			# get the points of the fit of that linear model
		f_of_x <- splinefun(smooth$x,smooth$y,method = splineMethod )
		if (log_spline == FALSE)
			{ splines_est = cbind( splines_est, f_of_x(time_splines) ) } else { 			splines_est = cbind( splines_est, exp(f_of_x(time_splines)) ) }										# store spline points for all dep. variables  	
		}
	
	     if (draw == TRUE)
     	       	{
			par(mfrow=c(round(sqrt(numVars),0)+1,round(sqrt(numVars),0)+1)) 
			for (i in 1:numVars)                                                                            			# cycling thought the dependent vars
		      	{	
            		plot(dataset[,1],dataset[,i+1],pch=20,col="grey",ylab=colnames(data_vars)[i])                               # plot the dependent variable dat
            		points(splines_est[,1],splines_est[,i+1],type='l',col="darkgreen",lwd=3)                      			# plot the spline		

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
			df = df
			) )
	}	

##### Smoother_for_MAR
#############################################################################



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
	# tSteps_ # number of steps of the simulation # NULL = tSteps_ will be equal to the number of timepoints in the data sample
	# nsim_ # number of simulations with noise # NULL = nsim_ will be equal to 1


	if ( log_transform == TRUE ) {data_vars = log(data_vars)}

	nData = dim(data_vars)[1]

	if (is.null(nData)) 
		{
		numVar = 1 
		nData = length(data_vars) 
		data_vars_toUSE = data_vars

		if ( demeaned == TRUE ) 
			{
			data_vars_means = mean(data_vars,na.rm = TRUE)
			data_vars_means_mat = rep(data_vars_means,nData)
			data_vars_toUSE = data_vars_demeaned = data_vars - data_vars_means_mat
			} 

		### fitting the MAR model.
		model.gen = list(
	     		B = "unconstrained", 		# matrix(list("b11","b12","b21","b22"),nrow=2,ncol=2,byrow=T),  		## Specify the "Full model" matrix for the MAR community matrix
           		U = "unequal", 			# matrix(list("a1","a2"),nrow=2,ncol=1), 							## state eq. intercept
           		Q = "diagonal and unequal",	# matrix(list("Q",0,0,"Q"),nrow=2,ncol=2),   					## state eq. var-cov matrix
			Z = "identity",               									  	  			## conversion matrice from space to state
           		A = "zero",     			# matrix(list(0,0),nrow=2,ncol=1), 							## space eq. intercept
           		R = "zero",	                          											## space eq. var-cov matrix
			tinitx=1				# if 0 it assumes that series will start at 0, if 1 it starts at one. Laim
			)

		if (estInit == FALSE) 
			{ 
			model.gen = c(model.gen,list( x0 = matrix(list(data_vars_toUSE[1])[[1]],nrow=numVar,ncol=1) )) 
			model.gen$tinitx = 1
			}

		MARvalEst1 = x_t = data_vars_toUSE[1]

		} else
		{
		numVar = dim(data_vars)[2]		
		data_vars_toUSE = t(data_vars)

		if ( demeaned == TRUE ) 
			{
			mean_no_na = function(sample1) { mean(sample1,na.rm=TRUE) }
			data_vars_means = apply(X=data_vars, MARGIN=2, FUN=mean_no_na)
			data_vars_toUSE = apply(X=data_vars,MARGIN=1,FUN=function(line1){line1-data_vars_means} )
			data_vars_demeaned = t( data_vars_toUSE )
			} 		



		### fitting the MAR model.
		model.gen = list(
	     		B = "unconstrained", 		# matrix(list("b11","b12","b21","b22"),nrow=2,ncol=2,byrow=T),  		## Specify the "Full model" matrix for the MAR community matrix
           		U = "unequal",			# matrix(list("a1","a2"),nrow=2,ncol=1), 						## state eq. intercept
           		Q = "diagonal and unequal",	# matrix(list("Q",0,0,"Q"),nrow=2,ncol=2),   					## state eq. var-cov matrix
			Z = "identity",               									  	  		## conversion matrice from space to state
           		A = "zero",     			# matrix(list(0,0),nrow=2,ncol=1), 							## space eq. intercept
           		R = "zero",	                          											## space eq. var-cov matrix
			tinitx=0				# if 0 it assumes that series will start at 0, if 1 it starts at one. Laim
           		#x0 = matrix(list(data_vars_toUSE[,1])[[1]],nrow=numVar,ncol=1)							## mean of the distrib from which starting values are drawn
			)

		if (estInit == FALSE) 
			{ 
			model.gen = c(model.gen,list( x0 = matrix(list(data_vars_toUSE[,1])[[1]],nrow=numVar,ncol=1) )) 
			model.gen$tinitx = 1
			}

		MARvalEst1 = x_t = data_vars_toUSE[,1]
		}

	library(MARSS)
	library(MASS)
	MARSS_output1 = MARSS(y=data_vars_toUSE ,model=model.gen,control=list(conv.test.slope.tol=0.01,abstol=abstol_,maxit=maxit_,allow.degen=TRUE),form="marxss",silent=F) 	# abstol=0.0001,maxit=10000
	try(m.ci1 <- MARSSparamCIs(MARSS_output1,method="hessian",alpha=0.05,nboot=100),silent=T)
		
	MARparEst1 = cbind(matrix(MARSS_output1$par$B, ncol=numVar),MARSS_output1$par$U,MARSS_output1$par$Q)

	if (estInit == TRUE) { MARvalEst1 = x_t = as.vector(MARSS_output1$par$x0) }

	if ( is.null(tSteps_) ) { tSteps_ = nData }

	for (i in 2:tSteps_)
		{
		x_tplus1 = MARparEst1[,1:numVar] %*% x_t + MARparEst1[,numVar+1] + includeNoise*rnorm(n=numVar, mean = rep(0,numVar), sd = sqrt(MARparEst1[,(numVar+2)]))		#
		MARvalEst1 = rbind(MARvalEst1,t(x_tplus1))
		x_t = x_tplus1	
		}

	if ( is.null(nsim_) ) {nsim_=1}
	MARSSsim1 = MARSSsimulate(object = MARSS_output1, tSteps = tSteps_, nsim = nsim_)$sim.data
	sim1dim = dim(MARSSsim1)
	
	remeaner = matrix(rep(0,numVar*tSteps_),ncol=numVar)
	remeaner_sim = array(rep(0,sim1dim[1]*sim1dim[2]*sim1dim[3]),dim=sim1dim)
	if ( demeaned == TRUE ) 
		{ 
		for (ix in 1:tSteps_) {remeaner[ix,] = data_vars_means} 
		for (iv in 1:nsim_) { remeaner_sim[,,iv] = t(remeaner) }
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
### fastMAR
#############################################################################



#############################################################################
### MAR - meanMaker

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

### MAR - meanMaker
#############################################################################



###############################################################################                                                                  

fourColPlot = function(fakedata,out_est)
	{
	# Function to create a four color plot
	varColor_light = c('lightblue','orange','lightgreen','grey')
	varColor = c('blue','orange3','darkgreen','black')
	par(mfrow=c(1,1))
	plot(1,1,col='white',xlim=c(0,length(out_est[,1])),ylim=c(0,3),xlab='t',ylab=NA,cex.lab=1.5,cex.axis=1.5)
	for  ( i in 2:((dim(fakedata)[2]-1)/2+1) )
		{
		points(fakedata[,1],fakedata[,i],pch=1,col=varColor[i-1])
		points(out_est[,1],out_est[,i],type='l',col=varColor_light[i-1],lwd=3)
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
##### colors

colorPallet = c('black','grey','blue','darkgreen','green')

##### colors
######################################################################################################################################



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
##### time, initState and truePars 

t = seq(1,100,1)

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

##### time, initState and truePars
######################################################################################################################################



######################################################################################################################################
##### truedata

out1 = solveLV(times = t, initState = state1, pars = truePars1_mat, equations = Equations)
# tail(out) # edit(edit) # View(out)

par(mfrow=c(2,2))

plot(out1[,1],out1[,2],pch=20,col='grey',lwd=3,xlab='t',ylab='X1')
plot(out1[,1],out1[,3],pch=20,col='grey',lwd=3,xlab='t',ylab='X2')
plot(out1[,1],out1[,4],pch=20,col='grey',lwd=3,xlab='t',ylab='X3')
plot(out1[,1],out1[,5],pch=20,col='grey',lwd=3,xlab='t',ylab='X4')

legend(60,1.0,legend = c('true model'),pch = c(20),col = c('grey'),bty = "n")	# ,lty = c(0,1)

##### truedata
######################################################################################################################################



######################################################################################################################################
### noisy dataset (noisy dataset 3) - choose 40 points at random and add a random normal of mean 0 and variance mean(var)*.5
 
set.seed(100)
noisy_data3 = out1[c(1,sort(sample(2:100,39))),1:5]
for (i in 2:5)
	{
	for (ii in 1:dim(noisy_data3)[1])
		{
		temp1 = -1
		while (temp1 < 0)
			{ 
			temp1 = noisy_data3[ii,i] + rnorm(1,0,.2 * mean(out1[,i]))
			}
		noisy_data3[ii,i] = temp1
		}
	}

# setwd('C:\\Users\\doliv\\OneDrive\\2_2019_America\\2020\\20200123_MAR\\Camparison_LV_MAR\\5_methods_of_ecology_and_evolution\\paper_scripts')
# write.csv(noisy_data3,'noisy_data3_20.csv')
# noisy_data3 = read.csv('noisy_data3_20.csv')
# save(noisy_data3, file = "noisy_data3_20.Rdata")
 load(file = "noisy_data3_20.Rdata") # setting a seed is not enough to replicate the same ramdomness of data in different computers. As so, we saved our results and are being loaded now.

# setwd('C:\\Users\\doliv\\OneDrive\\2_2019_America\\2020\\20200123_MAR\\Camparison_LV_MAR\\5_methods_of_ecology_and_evolution\\paper_figures')
# tiff("Fig_1a.tiff", height = 20, width = 20, units = 'cm',compression = "lzw", res = 300)
# Function to create a four independent color plot
par(mfrow=c(2,2))
for  ( i in 2:((dim(out1)[2]-1)/2+1) )
	{
	plot(out1[,1],out1[,i],type='p',pch=20,col='lightgrey',xlim=c(0,length(out1[,1])),ylim=c(0,3),xlab='t',ylab=varNames[i-1],cex.lab=1.5,cex.axis=1.5,lwd=1)
	points(noisy_data3[,1],noisy_data3[,i],pch=20,col=varColor[i-1])
	}
	legend(45,1,legend = c('LV results','noise replicates'),lty = c(0,0),pch = c(20,20),col = c('darkgrey','black'),text.col=c('black','black'),bty = "n",cex=1)
# dev.off()

### noisy dataset (noisy dataset 3) - choose 40 points at random and add a random normal of mean 0 and variance mean(var)*.5
######################################################################################################################################



#############################################################################
##### clustered dataset (noisy data 4) - 11 points were chosen from the artificial data and for each point five observations were created multiplying it by a random normal of mean 1 and variance 0.5

set.seed(100)
data1 = out1
( data1_sparsed = data1[c(1,6,11,18,26,33,41,51,60,70,80),] )  #  #seq(0,101,10) #c(1,6,11,18,26,41,seq(51,101,10))

datapoint_n = 5
noisy_data4 = matrix(rep(NA,5),nrow=1)
for (i in 1:dim(data1_sparsed)[1])
	{
      for (ii in 1:datapoint_n)
		{
		noisy_data4 = rbind( noisy_data4, c(data1_sparsed[i,1], data1_sparsed[i,2:5] * rnorm(n = 4,mean = 1,sd = .2)) )
		}
	}
noisy_data4 = noisy_data4[2:dim(noisy_data4)[1],]

# setwd('C:\\Users\\doliv\\OneDrive\\2_2019_America\\2020\\20200123_MAR\\Camparison_LV_MAR\\5_methods_of_ecology_and_evolution\\paper_scripts')
# write.csv(noisy_data4,'noisy_data4_20.csv')
# noisy_data4 = read.csv('noisy_data4_20.csv')
# save(noisy_data4, file = "noisy_data4_20.Rdata")
 load(file = "noisy_data4_20.Rdata") # setting a seed is not enough to replicate the same ramdomness of data in different computers. As so, we saved our results and are being loaded now.

# tiff("Fig_1b.tiff", height = 20, width = 20, units = 'cm',compression = "lzw", res = 300)
par(mfrow=c(2,2))
plot(data1[,1],data1[,2],pch=20,col='lightgrey',ylim=c(0,3),xlab='t',ylab='X1',cex.lab=1.5,cex.axis=1.5,lwd=1)
# points(data1_sparsed[,1],data1_sparsed[,2],pch=20)
points(noisy_data4[,1],noisy_data4[,2],pch=20,col=varColor[1])
plot(data1[,1],data1[,3],pch=20,col='lightgrey',ylim=c(0,3),xlab='t',ylab='X2',cex.lab=1.5,cex.axis=1.5,lwd=1)
# points(data1_sparsed[,1],data1_sparsed[,3],pch=20)
points(noisy_data4[,1],noisy_data4[,3],pch=20,col=varColor[2])
plot(data1[,1],data1[,4],pch=20,col='lightgrey',ylim=c(0,3),xlab='t',ylab='X3',cex.lab=1.5,cex.axis=1.5,lwd=1)
# points(data1_sparsed[,1],data1_sparsed[,4],pch=20)
points(noisy_data4[,1],noisy_data4[,4],pch=20,col=varColor[3])
plot(data1[,1],data1[,5],pch=20,col='lightgrey',ylim=c(0,3),xlab='t',ylab='X4',cex.lab=1.5,cex.axis=1.5,lwd=1)
# points(data1_sparsed[,1],data1_sparsed[,5],pch=20)
points(noisy_data4[,1],noisy_data4[,5],pch=20,col=varColor[4])
legend(45,1,legend = c('LV results','noisy replicates'),pch = c(20,20),col = c('grey','black'),bty = "n")
# dev.off()

##### clustered dataset (noisy data 4) - 11 points were chosen from the artificial data and for each point five observations were created multiplying it by a random normal of mean 1 and variance 0.5
#############################################################################







#############################################################################
##### figure S1 - noisy and clustered data

# setwd('C:\\Users\\doliv\\OneDrive\\2_2019_America\\2020\\20200123_MAR\\Camparison_LV_MAR\\8_methods_of_ecology_and_evolution_2sub\\paper_scripts')
# noisy_data3 = read.csv('noisy_data3_20.csv')
 load(file = "noisy_data3_20.Rdata")
# noisy_data4 = read.csv('noisy_data4_20.csv')
 load(file = "noisy_data4_20.Rdata")

# setwd('C:\\Users\\doliv\\OneDrive\\2_2019_America\\2020\\20200123_MAR\\Camparison_LV_MAR\\8_methods_of_ecology_and_evolution_2sub\\paper_figures')
# tiff("Fig_S1.tiff", height = 40, width = 20, units = 'cm',compression = "lzw", res = 300)
# Function to create a four independent color plot
par(mfcol=c(4,2))
mainLab = c('Noisy Dataset',NA,NA,NA)
y_labs = c(expression(italic('X'[1])),expression(italic('X'[2])),expression(italic('X'[3])),expression(italic('X'[4])))
for ( i in 2:5 )
	{
	plot(out1[,1],out1[,i],type='p',pch=20,col=colorPallet[2],xlim=c(0,length(out1[,1])),ylim=c(0,3),xlab='',main='',ylab='',cex.lab=1.5,cex.axis=1.5,lwd=1)
	points(noisy_data3[,1],noisy_data3[,i],col=colorPallet[1])	#,pch=20
	mtext(expression(italic('t')),side=1,line=2.1,cex=1.5)
	mtext(y_labs[i-1],side=2,line=2.1,cex=1.5)
	mtext(mainLab[i-1],side=3,line=1,cex=1.5)
	}
legend(45,1,legend = c('LV results','noisy dataset'),pch = c(20,1),col = c(colorPallet[2],colorPallet[1]),bty = "n",cex=1.3)
mainLab = c('Replicate Dataset',NA,NA,NA)
for (i in 2:5)
	{
	plot(data1[,1],data1[,i],pch=20,col=colorPallet[2],ylim=c(0,3),main='',xlab='',ylab='',cex.lab=1.5,cex.axis=1.5,lwd=1)
	points(noisy_data4[,1],noisy_data4[,i],col=colorPallet[1])	#,pch=20
	mtext(expression(italic('t')),side=1,line=2.1,cex=1.5)
	mtext(y_labs[i-1],side=2,line=2.1,cex=1.5)
	mtext(mainLab[i-1],side=3,line=1,cex=1.5)
	}
legend(45,1,legend = c('LV results','replicate dataset'),pch = c(20,1),col = c(colorPallet[2],colorPallet[1]),bty = "n",cex=1.3)

# dev.off()

##### figure S1 - noisy and clustered data 20%
#############################################################################







#############################################################################
##### figure S2 - smoothing on X1 

##########
### splines on X1

x = seq(head(out1[,1],1),tail(out1[,1],1),.1)
# setwd('C:\\Users\\doliv\\OneDrive\\2_2019_America\\2020\\20200123_MAR\\Camparison_LV_MAR\\8_methods_of_ecology_and_evolution_2sub\\paper_figures')
# tiff("Fig_S2.tiff", height = 40, width = 20, units = 'cm',compression = "lzw", res = 300)
par(mfcol=c(4,2))
DFlist = c(5,10,15,35)
mainLab = c('Splines',NA,NA,NA)
for (i in 1:4)
	{
	smoothSpline <- smooth.spline(noisy_data3[,1],noisy_data3[,2],df=DFlist[i])                         		# smoothing with degrees of freedom (df)
	smooth = predict(object = smoothSpline, x = x, deriv = 0)                                 			# get the points of the fit of that linear model
	f_of_x <- splinefun(smooth$x,smooth$y)
	plot(out1[,1],out1[,2],pch=20,col=colorPallet[2],ylim=c(0,3),main='',xlab='',ylab='',cex.lab=1.5,cex.axis=1.5,lwd=1)
	points(noisy_data3[,1],pch=1,noisy_data3[,2],col=colorPallet[1]) 
	legend(30,1,legend = c('LV results','noisy data',paste(DFlist[i],'DF-spline',sep="")),lty = c(0,0,1),pch = c(20,1,NA),col = c(colorPallet[2],colorPallet[1],'darkgreen'),text.col=c('black','black','black'),bty = "n",cex=1,lwd=c(1,1,3))
	points(x,f_of_x(x) ,type='l',lty=1,col='darkgreen',lwd=3)
	mtext(expression(italic('t')),side=1,line=2.1,cex=1.5)
	mtext(rep(y_labs[1],4),side=2,line=2.1,cex=1.5)
	mtext(mainLab[i],side=3,line=1,cex=1.5)
	}

### splines on X1
##########


##########
### LOESS on X1
spanList = c(0.9,.5,.2,.1)
mainLab = c('LOESS',NA,NA,NA) 
for (i in 1:4)
	{
	( est_data1_sparse_noisy = Smoother(
				dataset = noisy_data3,
				draw = FALSE,
				data1_spline2 = 2, 
				smooth_type = 2,
				splineMethod = "fmm",
				polyDegree012 = 1,
				aicc1_gvc2 = 1,
				span = spanList[i],	
				log_spline = FALSE
				) )

	plot(out1[,1],out1[,2],pch=20,col=colorPallet[2],ylim=c(0,3),main='',xlab='',ylab='',cex.lab=1.5,cex.axis=1.5,lwd=1)
	points(noisy_data3[,1],pch=1,noisy_data3[,2],col=colorPallet[1]) 
	legend(30,1,legend = c('LV results','noisy data',paste('LOESS span =',spanList[i])),lty = c(0,0,1),pch = c(20,1,NA),col = c(colorPallet[2],colorPallet[1],'darkgreen'),text.col=c('black','black','black'),bty = "n",cex=1,lwd=c(1,1,3))
	points(est_data1_sparse_noisy$splines_est[,1],est_data1_sparse_noisy$splines_est[,2],type='l',col='darkgreen',lwd=3)
	mtext(expression(italic('t')),side=1,line=2.1,cex=1.5)
	mtext(rep(y_labs[1],4),side=2,line=2.1,cex=1.5)
	mtext(mainLab[i],side=3,line=1,cex=1.5)
	}
# dev.off()

### LOESS on X1
##########

##### figure S2 - smoothing on X1 
#############################################################################







#############################################################################
##### find alg1 point sub sample for noisy_data3 and estimates

smoother_out = Smoother(
				dataset = noisy_data3,
				draw = TRUE,
				data1_spline2 = 2, 
				smooth_type = 1,
				df = c(5,8,11,5),	
				log_spline = TRUE
				)
# uncomment the next lines to find a subsample for this dataset
#AlgPointFind(
#		smoother_out_APF = smoother_out,
#		dif_equations = Equations,
#		matrixEq = TRUE
#		#randomSearchTries = 100
#		)

# "100 % complete || DF = 5 8 11 5  ||  4 5 6 7 33 17.981  ||  3 6 14 23 33 0.799"
# "100 % complete || DF = 5 6 10 5  ||  1 10 26 31 40 18.143  ||  3 7 14 23 33 0.62"
# "100 % complete || DF = 11 11 11 11  ||  4 6 10 15 24 1.357  ||  4 6 10 15 24 0.265"
# "100 % complete || DF = 3 3 11 5     ||  6 7 21 22 34 21.129  ||  3 10 19 34 40 0.528"

smoother_out_alg1 = smoother_out = Smoother(
				dataset = noisy_data3,
				draw = TRUE,
				data1_spline2 = 2, 
				smooth_type = 1,
				df = c(5,8,11,5),	
				log_spline = TRUE
				)
( est1 = LV_pars_finder(
				smooth_out = smoother_out,
				alg1_lm2 = 1, 
				data_sample_alg = c(4,5,6,7,33)#c(3, 6, 14, 21, 27)		#c(3, 7, 10, 14, 22)		
				) )

plot(1,1,col='white')
legend(1,1,legend = c('data','splines','slopes'),pch = c(20,NA,NA),lty = c(0,1,1),col = c('grey','green','lightgreen'),bty = "n")
cbind(est1)
est1_mat = Format_pars(truePars = est1)
cbind(truePars1,est1,absDif = abs(truePars1-est1))
splineStart = smoother_out$splines_est[1,2:5];names(splineStart) = c('x1','x2','x3','x4')
out_est_noisydata3_alg1 = out_est1 = solveLV(times = t, pars = est1_mat, initState = splineStart, equations = Equations)
#noisy_data3[1,2:5]

par(mfrow = c(2,2))
for (i in 2:5)
	{
	plot(out1[,1],out1[,i],pch=20,col='lightgrey',ylim=c(0,3),xlab='t',ylab=varNames[i-1],cex.lab=1.5,cex.axis=1.5,lwd=1)
	points(noisy_data3[,1],noisy_data3[,i],pch=20,col=varColor[i-1]) 
	points(smoother_out$splines_est[,1],smoother_out$splines_est[,i],type='l',col=varColor[i-1])
	points(out_est_noisydata3_alg1[,1],out_est_noisydata3_alg1[,i],type='l',col=varColor[i-1],lwd=3)
	if (i == 2) {legend(35,1.2,legend = c('data','sparse noisy data','smooth','estimates'),lty = c(0,0,1,1),pch = c(20,20,NA,NA),col = c('darkgrey',varColor[i-1],varColor[i-1],varColor[i-1]),text.col=c('black','black','black'),bty = "n",cex=1,lwd=c(NA,NA,1,3))}
	}

# SSE original data againts alg solution (this is in table one)
sum( ( out1[,1:5] - out_est_noisydata3_alg1[,1:5] )^2 )

# SSE noisy data againts alg solution (this is to justify why alg pars est are not similar to true pars)
sum( ( noisy_data3[,1:5] - out_est_noisydata3_alg1[noisy_data3[,1],1:5] )^2 )

# SSE noisy data againts original data (this is to justify why alg pars est are not similar to true pars)
sum( ( noisy_data3[,1:5] - out1[noisy_data3[,1],1:5] )^2 )

##### find alg1 point sub sample for noisy_data3 and estimates
############################################################################################################################


############################################################################################################################
##### MAR for noisy data 3

# Completing the data to have the same time structure as the original data
TS_noisy_data = cbind(1:100,matrix(rep(NA,100*4),ncol=4))
for (i in 1:100)
	{
	if (i%in%noisy_data3[,1])
		{
		TS_noisy_data[i,] = noisy_data3[noisy_data3[,1]==i,] 
		} 
	}

# MAR estimates without data log transformation
MARest1 =  fastMAR(
			data_vars = TS_noisy_data[,2:5],
			log_transform=FALSE,
			demeaned=FALSE,
			abstol_=0.01,
			maxit_ = 3000,
			includeNoise = FALSE,
			estInit = FALSE,
			tSteps_ = NULL,
			nsim_ = 100
			) 
MARest1$MARSSci_output

# MAR estimates with data log transformation
MARest2 =  fastMAR(
			data_vars = TS_noisy_data[,2:5],
			log_transform=TRUE,
			demeaned=FALSE,
			abstol_=0.01,
			maxit_ = 3000,
			includeNoise = FALSE,
			estInit = FALSE,
			tSteps_ = NULL,
			nsim_ = 100
			) 
MARest2$MARSSci_output

# plots
par(mfrow=c(2,2))
for ( i in 1:(dim(TS_noisy_data)[2]-1) )
	{
	plot(out1[,1],out1[,i+1],type='l',col='lightgrey',lwd=3,ylim=c(0,3),ylab = varNames[i])

	# for (ii in 1:100) { points(MARest1$MARSSsim_data[i,,ii],type='l',col='grey') }
	points(meanMaker(MARest1)[i,],type='l',col = 'purple')

	# for (ii in 1:100) { points(MARest2$MARSSsim_data[i,,ii],type='l',col='grey') }
	points(meanMaker(MARest2)[i,],type='l',col = 'purple',lty=2)


	points(out1[,1],out1[,i+1],type='l',col='lightgrey',lwd=3,ylim=c(0,3),ylab = varNames[i])
	points(TS_noisy_data[,i+1],pch=20,col=varColor[i],xlab = 'n')
	points(MARest1$MAR_TS_Est[,i],col='green',type='l',lty=1)
	points(MARest2$MAR_TS_Est[,i],col='green',type='l',lty=2)
	if (i==1) {legend(40,1.8,legend=c('original data','noisy data','MAR no transform','MAR log transform'),pch=c(NA,20,NA,NA),lty=c(1,NA,1,2),col=c('lightgrey',varColor[i],varColor_light[i],varColor_light[i]),bty = "n",cex=1,lwd=c(3,NA,1,1))}

	# expected states
	lines(MARest1$MARSS_output$states[i,])
	#lines(MARest1$MARSS_output$states[i,] - 2*MARest1$MARSS_output$states.se[i,],lty=3)
	#lines(MARest1$MARSS_output$states[i,] + 2*MARest1$MARSS_output$states.se[i,],lty=3)

	# expected log states
	#lines(exp(MARest2$MARSS_output$states[i,]))
	#lines(exp(MARest2$MARSS_output$states[i,] - 2*MARest1$MARSS_output$states.se[i,]),lty=3)
	#lines(exp(MARest2$MARSS_output$states[i,] + 2*MARest1$MARSS_output$states.se[i,]),lty=3)
	}

# put time back in the data
MARest1_orig = cbind(1:dim(MARest1$MAR_TS_Est)[1],MARest1$MAR_TS_Est)
MARest1_log = cbind(1:dim(MARest2$MAR_TS_Est)[1],MARest2$MAR_TS_Est)

# errors
sum( ( out1[,2:5] - MARest1_orig[,2:5] )^2 )
sum( ( out1[,2:5] - MARest1_log[,2:5] )^2 )

sum( ( out1[,2:5] - t(MARest1$MARSS_output$states) )^2 )
sum( ( out1[,2:5] - t(exp(MARest2$MARSS_output$states)) )^2 )

##### MAR for noisy data 3
###############################################################################################################


############################################################################################################################
##### MAR with smooth data for noisy data 3

# smoothing the data for MAR
smoother_for_MAR_out = Smoother_for_MAR(	
						dataset = noisy_data3,
						df = c(5,8,11,5),
						dataWeights = NULL,
						splineMethod = "fmm",
						draw = TRUE,
						log_spline = TRUE
						)
smoother_for_MAR_out
smooth_TS = smoother_for_MAR_out$splines_est

# Completing the smooth data to have the same time structure as the original data
TS_smooth_data = cbind(1:100,matrix(rep(NA,100*4),ncol=4))
for (i in 1:100)
	{
	if (i%in%smooth_TS[,1])
		{
		TS_smooth_data[i,] = smooth_TS[smooth_TS[,1]==i,] 
		} 
	}

# MAR estimates without data log transformation, with spline smooth data  
MARest1_smoothData =  fastMAR(
			data_vars = TS_smooth_data[,2:5],
			log_transform=FALSE,
			demeaned=TRUE,
			abstol_=0.01,
			maxit_ = 1000,
			includeNoise = FALSE,
			estInit = FALSE,
			tSteps_ = NULL,
			nsim_ = 100
			) 
MARest1_smoothData$MARSSci_output

# MAR estimates with data log transformation, with spline smooth data  
MARest2_smoothData =  fastMAR(
			data_vars = TS_smooth_data[,2:5],
			log_transform=TRUE,
			demeaned=TRUE,
			abstol_=0.01,
			maxit_ = 1000,
			includeNoise = FALSE,
			estInit = FALSE,
			tSteps_ = NULL,
			nsim_ = 100
			) 
MARest2_smoothData$MARSSci_output

# plots
par(mfrow=c(2,2))
for ( i in 1:(dim(TS_noisy_data)[2]-1) )
	{
	plot(out1[,1],out1[,i+1],type='l',col='lightgrey',lwd=3,ylim=c(0,3),ylab = varNames[i])		# original data

	# for (ii in 1:100) { points(MARest1$MARSSsim_data[i,,ii],type='l',col='grey') }			# MAR estimates without data log transformation, with spline smooth data, with noise
	points(meanMaker(MARest1_smoothData)[i,],type='l',col = 'purple')

	# for (ii in 1:100) { points(MARest2$MARSSsim_data[i,,ii],type='l',col='grey') }
	points(meanMaker(MARest2_smoothData)[i,],type='l',col = 'purple',lty=2)					# MAR estimates with data log transformation, with spline smooth data, with noise

	points(TS_smooth_data[,i+1],pch=20,col=varColor[i],xlab = 'n')						# treated data
	points(MARest1_smoothData$MAR_TS_Est[,i],col='green',type='l',lty=1)						# MAR estimates without data log transformation, with spline smooth data
	points(MARest2_smoothData$MAR_TS_Est[,i],col='green',type='l',lty=2)						# MAR estimates with data log transformation, with spline smooth data
	if (i==1) {legend(40,1.8,legend=c('original data','noisy data','MAR no transform','MAR log transform'),pch=c(NA,20,NA,NA),lty=c(1,NA,1,2),col=c('lightgrey',varColor[i],varColor_light[i],varColor_light[i]),bty = "n",cex=1,lwd=c(3,NA,1,1))}

	# expected states
	lines(MARest1_smoothData$MARSS_output$states[i,]+mean(TS_smooth_data[,i+1],na.rm = TRUE))
	#lines(MARest1_smoothData$MARSS_output$states[i,]+mean(TS_smooth_data[,i+1],na.rm = TRUE) - 2*MARest1_smoothData$MARSS_output$states.se[i,],lty=3)
	#lines(MARest1_smoothData$MARSS_output$states[i,]+mean(TS_smooth_data[,i+1],na.rm = TRUE) + 2*MARest1_smoothData$MARSS_output$states.se[i,],lty=3)

	# expected log states
	#lines(exp(MARest2_smoothData$MARSS_output$states[i,]+mean(log(TS_smooth_data[,i+1]),na.rm = TRUE)))
	#lines(exp(MARest2_smoothData$MARSS_output$states[i,]+mean(log(TS_smooth_data[,i+1]),na.rm = TRUE)) - 2*MARest2_smoothData$MARSS_output$states.se[i,],lty=3)
	#lines(exp(MARest2_smoothData$MARSS_output$states[i,]+mean(log(TS_smooth_data[,i+1]),na.rm = TRUE)) + 2*MARest2_smoothData$MARSS_output$states.se[i,],lty=3)
	}

# put time back in the data
MARest1_smoothData_orig = cbind(1:dim(MARest1_smoothData$MAR_TS_Est)[1],MARest1_smoothData$MAR_TS_Est)
MARest1_smoothData_log = cbind(1:dim(MARest2_smoothData$MAR_TS_Est)[1],MARest2_smoothData$MAR_TS_Est)



# errors
sum( ( out1[,2:5] - MARest1_smoothData_orig[,2:5] )^2 )
sum( ( out1[,2:5] - MARest1_smoothData_log[,2:5] )^2 )

MARest1_smoothData_remeaned = MARest1_smoothData$MARSS_output$states
for (i in 1:4)
	{
	MARest1_smoothData_remeaned[i,] = MARest1_smoothData$MARSS_output$states[i,]+mean(TS_smooth_data[,i+1],na.rm = TRUE)
	}
sum( ( out1[,2:5] - t(MARest1_smoothData_remeaned) )^2 )

MARest2_smoothData_remeaned = MARest2_smoothData$MARSS_output$states
for (i in 1:4)
	{
	MARest2_smoothData_remeaned[i,] = exp(MARest2_smoothData$MARSS_output$states[i,]+mean(log(TS_smooth_data[,i+1]),na.rm = TRUE))
	}
sum( ( out1[,2:5] - t(MARest2_smoothData_remeaned) )^2 )

##### MAR with smooth data for noisy data 3
############################################################################################################################





#################################################################################################
##### noisy data4 - alg1 solution

### first attempt - just put it in

smoother_out = Smoother(
				dataset = noisy_data4,
				draw = TRUE,
				data1_spline2 = 2, 
				smooth_type = 1,
				df = c(5,6,8,8),	# c(11,11,11,11) for lm2 # c(8,8,8,8) for alg1	
				log_spline = TRUE
				)

( est1 = LV_pars_finder(
				smooth_out = smoother_out,
				alg1_lm2 = 1, 
				data_sample_alg = 5*c(3,4,5,6,9)	# data_sample_lm = NULL for lm2	# data_sample_alg = 5*c(2,4,5,6,9) for alg1
				) )

plot(1,1,col='white')
legend(1,1,legend = c('data','splines','slopes'),pch = c(20,NA,NA),lty = c(0,1,1),col = c('grey','green','lightgreen'),bty = "n")
cbind(est1)
est1_mat = Format_pars(truePars = est1)
cbind(truePars1,est1,absDif = abs(truePars1-est1))
splineStart = smoother_out$splines_est[1,2:5];names(splineStart) = c('x1','x2','x3','x4')
out_est_noisy_data4_alg = out_est1 = solveLV(times = t, pars = est1_mat, initState = splineStart, equations = Equations)
#noisy_data3[1,2:5]

# setwd("C:\\Users\\doliv\\OneDrive\\2_2019_America\\2020\\20200123_MAR\\Camparison_LV_MAR\\3_paper_new_functions\\paper_figures")     # P51
# tiff("figS2.tiff", height = 20, width = 20, units = 'cm',compression = "lzw", res = 300)
par(mfrow = c(2,2))
for (i in 2:5)
	{
	plot(out1[,1],out1[,i],pch=20,col='lightgrey',ylim=c(0,3),xlab='t',ylab=varNames[i-1],cex.lab=1.5,cex.axis=1.5,lwd=1)
	points(noisy_data4[,1],noisy_data4[,i],pch=20,col=varColor[i-1]) 
	points(smoother_out$splines_est[,1],smoother_out$splines_est[,i],type='l',col=varColor[i-1])
	points(out_est_noisy_data4_alg[,1],out_est_noisy_data4_alg[,i],type='l',col=varColor[i-1],lwd=3)
	if (i == 2) {legend(35,1.2,legend = c('data','sparse noisy data','smooth','estimates'),lty = c(0,0,1,1),pch = c(20,20,NA,NA),col = c('darkgrey',varColor[i-1],varColor[i-1],varColor[i-1]),text.col=c('black','black','black'),bty = "n",cex=1,lwd=c(NA,NA,1,3))}
	}
# dev.off()

sum( ( out1[,2:5] - out_est_noisy_data4_alg[,2:5] )^2 )

##### noisy data4 - alg1 solution
#################################################################################################


#################################################################################################
##### MAR for noisy data 4

noisy_data4_for_MAR = matrix(rep(NA,(dim(noisy_data4)[1]/5)*dim(noisy_data4)[2]),ncol = dim(noisy_data4)[2])
noisy_data4_for_MAR[,1] = unique(noisy_data4[,1])
for (i in 0:(dim(noisy_data4_for_MAR)[1]-1))
	{
	noisy_data4_for_MAR[i+1,2:5] = apply(X = noisy_data4[(5*i+1):(5*i+5),2:5],MARGIN = 2, FUN='mean')
	}

TS_noisy_data4 = cbind(1:100,matrix(rep(NA,100*4),ncol=4))
for (i in 1:100)
	{
	if (i%in%noisy_data4_for_MAR[,1])
		{
		TS_noisy_data4[i,] = noisy_data4_for_MAR[noisy_data4_for_MAR[,1]==i,] 
		} 
	}

# setwd('C:\\Users\\doliv\\OneDrive\\2_2019_America\\2020\\20200123_MAR\\Camparison_LV_MAR\\3_paper_new_functions\\paper_scripts')
# write.csv(TS_noisy_data4,'TS_noisy_data4.csv')
# TS_noisy_data4 = read.csv('TS_noisy_data4.csv')
# save(TS_noisy_data4, file = "TS_noisy_data4.Rdata")
# load(file = "TS_noisy_data4.Rdata")

MARest3 =  fastMAR(
			data_vars = TS_noisy_data4[,2:5],
			log_transform = FALSE,
			demeaned = TRUE,
			abstol_= 0.01,
			maxit_ = 1000,
			includeNoise = FALSE,			
			estInit = FALSE,
			tSteps_ = NULL,
			nsim_ = 100
			) 
MARest3$MARSSci_output

# SSE
sum( ( out1[,2:5] - MARest3$MAR_TS_Est )^2 )


MARest4 =  fastMAR(
			data_vars = TS_noisy_data4[,2:5],
			log_transform = TRUE,
			demeaned = TRUE,
			abstol_= 0.01,
			maxit_ = 1000,
			includeNoise = FALSE,			
			estInit = FALSE,
			tSteps_ = NULL,
			nsim_ = 100
			) 
MARest4$MARSSci_output

# SSE
sum( ( out1[,2:5] - MARest4$MAR_TS_Est )^2 )


par(mfrow=c(2,2))
for ( i in 1:(dim(TS_noisy_data4)[2]-1) )
	{
	plot(out1[,1],out1[,i+1],type='l',col='lightgrey',lwd=3,ylim=c(0,3),ylab = varNames[i])
	points(noisy_data4[,1],noisy_data4[,i+1],pch=20,col=varColor[i],xlab = 'n')
	points(MARest3$MAR_TS_Est[,i],col=varColor_light[i],type='l',lty=1)
	points(MARest4$MAR_TS_Est[,i],col=varColor_light[i],type='l',lty=2)
	if (i==1) {legend(40,1.8,legend=c('original data','noisy data','MAR no transform','MAR log transform'),pch=c(NA,20,NA,NA),lty=c(1,NA,1,2),col=c('lightgrey',varColor[i],varColor_light[i],varColor_light[i]),bty = "n",cex=1,lwd=c(3,NA,1,1))}

	# for (ii in 1:100) { points(MARest3$MARSSsim_data[i,,ii],type='l',col='grey') }
	points(meanMaker(MARest3)[i,],type='l',col = 'purple')

	# for (ii in 1:100) { points(MARest4$MARSSsim_data[i,,ii],type='l',col='grey') }
	points(meanMaker(MARest4)[i,],type='l',col = 'purple',lty=2)

	# expected states
	lines(MARest3$MARSS_output$states[i,] + mean(TS_noisy_data4[,i+1],na.rm = TRUE))
	#lines(MARest3$MARSS_output$states[i,]+ mean(TS_noisy_data4[,i+1],na.rm = TRUE) - 2*MARest3$MARSS_output$states.se[i,],lty=3)
	#lines(MARest3$MARSS_output$states[i,]+ mean(TS_noisy_data4[,i+1],na.rm = TRUE) + 2*MARest3$MARSS_output$states.se[i,],lty=3)

	# expected log states
	#lines(exp(MARest4$MARSS_output$states[i,] + mean(log(TS_noisy_data4[,i+1]),na.rm = TRUE)))
	#lines(exp(MARest4$MARSS_output$states[i,] + mean(log(TS_noisy_data4[,i+1]),na.rm = TRUE) - 2*MARest4$MARSS_output$states.se[i,]),lty=3)
	#lines(exp(MARest4$MARSS_output$states[i,] + mean(log(TS_noisy_data4[,i+1]),na.rm = TRUE) + 2*MARest4$MARSS_output$states.se[i,]),lty=3)
	}

MARest3_remeaned = MARest3$MARSS_output$states
for (i in 1:4)
	{
	MARest3_remeaned[i,] = MARest3$MARSS_output$states[i,] + mean(TS_noisy_data4[,i+1],na.rm = TRUE)
	}
sum( ( out1[,2:5] - t(MARest3_remeaned) )^2 )

MARest4_remeaned = MARest4$MARSS_output$states
for (i in 1:4)
	{
	MARest4_remeaned[i,] = exp(MARest4$MARSS_output$states[i,] + mean(log(TS_noisy_data4[,i+1]),na.rm = TRUE))
	}
sum( ( out1[,2:5] - t(MARest4_remeaned) )^2 )

##### MAR for noisy data 4
#######################################################################################3


############################################################################################################################
##### MAR with smooth data for noisy data 4

# smoothing the data for MAR
smoother_for_MAR_out = Smoother_for_MAR(	
						dataset = noisy_data4,
						df = c(5,8,11,5),
						dataWeights = NULL,
						splineMethod = "fmm",
						draw = TRUE,
						log_spline = TRUE
						)
smoother_for_MAR_out
smooth_TS = smoother_for_MAR_out$splines_est

# Completing the smooth data to have the same time structure as the original data
TS_smooth_data = cbind(1:100,matrix(rep(NA,100*4),ncol=4))
for (i in 1:100)
	{
	if (i%in%smooth_TS[,1])
		{
		TS_smooth_data[i,] = smooth_TS[smooth_TS[,1]==i,] 
		} 
	}

# MAR estimates without data log transformation, with spline smooth data  
MARest3_smoothData =  fastMAR(
			data_vars = TS_smooth_data[,2:5],
			log_transform=FALSE,
			demeaned=TRUE,
			abstol_=0.01,
			maxit_ = 1000,
			includeNoise = FALSE,
			estInit = FALSE,
			tSteps_ = NULL,
			nsim_ = 100
			) 
MARest3_smoothData$MARSSci_output

# MAR estimates with data log transformation, with spline smooth data  
MARest4_smoothData =  fastMAR(
			data_vars = TS_smooth_data[,2:5],
			log_transform=TRUE,
			demeaned=TRUE,
			abstol_=0.01,
			maxit_ = 1000,
			includeNoise = FALSE,
			estInit = FALSE,
			tSteps_ = NULL,
			nsim_ = 100
			) 
MARest4_smoothData$MARSSci_output

# plots
par(mfrow=c(2,2))
for ( i in 1:(dim(TS_noisy_data)[2]-1) )
	{
	plot(out1[,1],out1[,i+1],type='l',col='lightgrey',lwd=3,ylim=c(0,3),ylab = varNames[i])		# original data

	# for (ii in 1:100) { points(MARest1$MARSSsim_data[i,,ii],type='l',col='grey') }			# MAR estimates without data log transformation, with spline smooth data, with noise
	points(meanMaker(MARest1_smoothData)[i,],type='l',col = 'purple')

	# for (ii in 1:100) { points(MARest2$MARSSsim_data[i,,ii],type='l',col='grey') }
	points(meanMaker(MARest2_smoothData)[i,],type='l',col = 'purple',lty=2)					# MAR estimates with data log transformation, with spline smooth data, with noise

	points(TS_smooth_data[,i+1],pch=20,col=varColor[i],xlab = 'n')						# treated data
	points(MARest3_smoothData$MAR_TS_Est[,i],col='green',type='l',lty=1)						# MAR estimates without data log transformation, with spline smooth data
	points(MARest4_smoothData$MAR_TS_Est[,i],col='green',type='l',lty=2)						# MAR estimates with data log transformation, with spline smooth data
	if (i==1) {legend(40,1.8,legend=c('original data','noisy data','MAR no transform','MAR log transform'),pch=c(NA,20,NA,NA),lty=c(1,NA,1,2),col=c('lightgrey',varColor[i],varColor_light[i],varColor_light[i]),bty = "n",cex=1,lwd=c(3,NA,1,1))}

	# expected states
	lines(MARest3_smoothData$MARSS_output$states[i,]+mean(TS_smooth_data[,i+1],na.rm = TRUE))
	#lines(MARest3_smoothData$MARSS_output$states[i,]+mean(TS_smooth_data[,i+1],na.rm = TRUE) - 2*MARest1_smoothData$MARSS_output$states.se[i,],lty=3)
	#lines(MARest3_smoothData$MARSS_output$states[i,]+mean(TS_smooth_data[,i+1],na.rm = TRUE) + 2*MARest1_smoothData$MARSS_output$states.se[i,],lty=3)

	# expected log states
	#lines(exp(MARest4_smoothData$MARSS_output$states[i,]+mean(log(TS_smooth_data[,i+1]),na.rm = TRUE)))
	#lines(exp(MARest4_smoothData$MARSS_output$states[i,]+mean(log(TS_smooth_data[,i+1]),na.rm = TRUE)) - 2*MARest2_smoothData$MARSS_output$states.se[i,],lty=3)
	#lines(exp(MARest4_smoothData$MARSS_output$states[i,]+mean(log(TS_smooth_data[,i+1]),na.rm = TRUE)) + 2*MARest2_smoothData$MARSS_output$states.se[i,],lty=3)
	}

# put time back in the data
MARest3_smoothData_orig = cbind(1:dim(MARest3_smoothData$MAR_TS_Est)[1],MARest3_smoothData$MAR_TS_Est)
MARest4_smoothData_log = cbind(1:dim(MARest4_smoothData$MAR_TS_Est)[1],MARest4_smoothData$MAR_TS_Est)

# errors
sum( ( out1[,2:5] - MARest3_smoothData_orig[,2:5] )^2 )
sum( ( out1[,2:5] - MARest4_smoothData_log[,2:5] )^2 )

MARest3_smoothData_remeaned = MARest3_smoothData$MARSS_output$states
for (i in 1:4)
	{
	MARest3_smoothData_remeaned[i,] = MARest3_smoothData$MARSS_output$states[i,]+mean(TS_smooth_data[,i+1],na.rm = TRUE)
	}
sum( ( out1[,2:5] - t(MARest3_smoothData_remeaned) )^2 )

MARest4_smoothData_remeaned = MARest4_smoothData$MARSS_output$states
for (i in 1:4)
	{
	MARest4_smoothData_remeaned[i,] = exp(MARest4_smoothData$MARSS_output$states[i,]+mean(log(TS_smooth_data[,i+1]),na.rm = TRUE))
	}
sum( ( out1[,2:5] - t(MARest4_smoothData_remeaned) )^2 )

##### MAR with smooth data for noisy data 4
############################################################################################################################





#############################################################################
##### figure 1 - noisy and clustered data 20% - alg1

# setwd('C:\\Users\\doliv\\OneDrive\\2_2019_America\\2020\\20200123_MAR\\Camparison_LV_MAR\\6_methods_of_ecology_and_evolution_new_images\\paper_scripts')
# noisy_data3 = read.csv('noisy_data3.csv')
# load(file = "noisy_data3.Rdata")
# noisy_data4 = read.csv('noisy_data4.csv')
# load(file = "noisy_data4.Rdata")

# setwd('C:\\Users\\doliv\\OneDrive\\2_2019_America\\2020\\20200123_MAR\\Camparison_LV_MAR\\8_methods_of_ecology_and_evolution_2sub\\paper_figures')
# tiff("Fig_1.tiff", height = 40, width = 20, units = 'cm',compression = "lzw", res = 300)
# Function to create a four independent color plot
par(mfcol=c(4,2))
mainLab = c('Noisy Dataset',NA,NA,NA)
for (i in 2:5)
	{
	plot(out1[,1],out1[,i],pch=20,col=colorPallet[2],ylim=c(0,3),main='',xlab='',ylab='',cex.lab=1.5,cex.axis=1.5,lwd=1)
	points(noisy_data3[,1],noisy_data3[,i],pch=1,col=colorPallet[1]) 
	#points(smoother_out_alg1$splines_est[,1],smoother_out_alg1$splines_est[,i],type='l',col=varColor[i-1])
	points(MARest1$MAR_TS_Est[,i-1],col=colorPallet[4],type='l',lty=1,lwd=3)
	points(MARest2$MAR_TS_Est[,i-1],col=colorPallet[5],type='l',lty=1,lwd=3)
	points(MARest1_smoothData_orig[,i],col='orange',type='l',lty=1,lwd=3)
	points(MARest1_smoothData_log[,i],col='khaki',type='l',lty=1,lwd=3)	
	points(out_est_noisydata3_alg1[,1],out_est_noisydata3_alg1[,i],type='l',col=colorPallet[3],lwd=3)
	mtext(expression(italic('t')),side=1,line=2.1,cex=1.5)
	mtext(y_labs[i-1],side=2,line=2.1,cex=1.5)
	mtext(mainLab[i-1],side=3,line=1,cex=1.5)
	if (i == 2) {legend(35,1.25,legend = c('LV results','noisy dataset','ALVI-MI method','MAR no transform','MAR log transform','MAR no transform with smoothing','MAR log transform with smoothing'),lty = c(0,0,1,1,1,1,1),pch = c(20,1,NA,NA,NA,NA,NA),col = c(colorPallet[2],colorPallet[1],colorPallet[3],colorPallet[4],colorPallet[5],'orange','khaki'),text.col=c('black','black','black','black','black','black','black'),bty = "n",cex=1,lwd=c(NA,NA,3,3,3,3,3))}
	}

MARest1_smoothData_orig

mainLab = c('Replicate Dataset',NA,NA,NA)
for (i in 2:5)
	{
	plot(out1[,1],out1[,i],pch=20,col=colorPallet[2],ylim=c(0,3),main='',xlab='',ylab='',cex.lab=1.5,cex.axis=1.5,lwd=1)
	points(noisy_data4[,1],noisy_data4[,i],pch=1,col=colorPallet[1]) 
	#points(smoother_out$splines_est[,1],smoother_out$splines_est[,i],type='l',col=varColor[i-1])
	points(MARest3$MAR_TS_Est[,i-1],col=colorPallet[4],type='l',lty=1,lwd=3)
	points(MARest4$MAR_TS_Est[,i-1],col=colorPallet[5],type='l',lty=1,lwd=3)
	points(MARest3_smoothData_orig[,i],col='orange',type='l',lty=1,lwd=3)
	points(MARest4_smoothData_log[,i],col='khaki',type='l',lty=1,lwd=3)
	points(out_est_noisy_data4_alg[,1],out_est_noisy_data4_alg[,i],type='l',col=colorPallet[3],lwd=3)
	mtext(expression(italic('t')),side=1,line=2.1,cex=1.5)
	mtext(y_labs[i-1],side=2,line=2.1,cex=1.5)
	mtext(mainLab[i-1],side=3,line=1,cex=1.5)
	if (i == 2) {legend(35,1.25,legend = c('LV results','replicate dataset','ALVI-MI method','MAR no transform','MAR log transform','MAR no transform with smoothing','MAR log transform with smoothing'),lty = c(0,0,1,1,1,1,1),pch = c(20,1,NA,NA,NA,NA,NA),col = c(colorPallet[2],colorPallet[1],colorPallet[3],colorPallet[4],colorPallet[5],'orange','khaki'),text.col=c('black','black','black','black','black','black','black'),bty = "n",cex=1,lwd=c(NA,NA,3,3,3,3,3))}
	}
# dev.off()

##### figure 1 - noisy and clustered data 20% - alg1
#############################################################################







#############################################################################
##### lm2 solution for noisy data 3

smoother_out = Smoother(
				dataset = noisy_data3,
				draw = TRUE,
				data1_spline2 = 2, 
				smooth_type = 1,
				df = c(6,11,11,11),	#c(6,11,15,15),	
				log_spline = TRUE
				)
( est1 = LV_pars_finder(
				smooth_out = smoother_out,
				alg1_lm2 = 2, 
				data_sample_lm = NULL				# 1:15
				) )

plot(1,1,col='white')
legend(1,1,legend = c('data','splines','slopes'),pch = c(20,NA,NA),lty = c(0,1,1),col = c('grey','green','lightgreen'),bty = "n")
cbind(est1)
est1_mat = Format_pars(truePars = est1)
cbind(truePars1,est1,absDif = abs(truePars1-est1))
splineStart = smoother_out$splines_est[1,2:5];names(splineStart) = c('x1','x2','x3','x4')
out_est_noisy_data3_lm = solveLV(times = t, pars = est1_mat, initState = splineStart, equations = Equations)
#noisy_data3[1,2:5]

par(mfrow = c(2,2))
for (i in 2:5)
	{
	plot(out1[,1],out1[,i],pch=20,col='lightgrey',ylim=c(0,3),xlab='t',ylab=varNames[i-1],cex.lab=1.5,cex.axis=1.5,lwd=1)
	points(noisy_data3[,1],noisy_data3[,i],pch=20,col=varColor[i-1]) 
	points(smoother_out$splines_est[,1],smoother_out$splines_est[,i],type='l',col=varColor[i-1])
	points(out_est_noisy_data3_lm[,1],out_est_noisy_data3_lm[,i],type='l',col=varColor[i-1],lwd=3)
	if (i == 2) {legend(35,1.2,legend = c('data','sparse noisy data','smooth','estimates'),lty = c(0,0,1,1),pch = c(20,20,NA,NA),col = c('darkgrey',varColor[i-1],varColor[i-1],varColor[i-1]),text.col=c('black','black','black'),bty = "n",cex=1,lwd=c(NA,NA,1,3))}
	}

# SSE
sum( ( out1[,2:5] - out_est_noisy_data3_lm[,2:5] )^2 )

##### lm2 solution for noisy data 3
#############################################################################



#############################################################################
##### lm2 solution for noisy data 4

### first attempt - just put it in

smoother_out = Smoother(
				dataset = noisy_data4,
				draw = TRUE,
				data1_spline2 = 2, 
				smooth_type = 1,
				df = c(8,8,8,8),	# c(11,11,11,11) for lm2 # c(8,8,8,8) for alg1	
				log_spline = TRUE
				)
( est1 = LV_pars_finder(
				smooth_out = smoother_out,
				alg1_lm2 = 2, 
				data_sample_alg = 5*c(2,4,5,6,9)	# data_sample_lm = NULL for lm2	# data_sample_alg = 5*c(2,4,5,6,9) for alg1
				) )

plot(1,1,col='white')
legend(1,1,legend = c('data','splines','slopes'),pch = c(20,NA,NA),lty = c(0,1,1),col = c('grey','green','lightgreen'),bty = "n")
cbind(est1)
est1_mat = Format_pars(truePars = est1)
cbind(truePars1,est1,absDif = abs(truePars1-est1))
splineStart = smoother_out$splines_est[1,2:5];names(splineStart) = c('x1','x2','x3','x4')
out_est_noisy_data4_lm = solveLV(times = t, pars = est1_mat, initState = splineStart, equations = Equations)
#noisy_data3[1,2:5]

par(mfrow = c(2,2))
for (i in 2:5)
	{
	plot(out1[,1],out1[,i],pch=20,col='lightgrey',ylim=c(0,3),xlab='t',ylab=varNames[i-1],cex.lab=1.5,cex.axis=1.5,lwd=1)
	points(noisy_data4[,1],noisy_data4[,i],pch=20,col=varColor[i-1]) 
	points(smoother_out$splines_est[,1],smoother_out$splines_est[,i],type='l',col=varColor[i-1])
	points(out_est_noisy_data4_lm[,1],out_est_noisy_data4_lm[,i],type='l',col=varColor[i-1],lwd=3)
	if (i == 2) {legend(35,1.2,legend = c('data','sparse noisy data','smooth','estimates'),lty = c(0,0,1,1),pch = c(20,20,NA,NA),col = c('darkgrey',varColor[i-1],varColor[i-1],varColor[i-1]),text.col=c('black','black','black'),bty = "n",cex=1,lwd=c(NA,NA,1,3))}
	}

# SSE
sum( ( out1[,2:5] - out_est_noisy_data4_lm[,2:5] )^2 )


##### lm2 solution for noisy data 4
#############################################################################



#############################################################################
##### figure S4 - noisy and clustered data 20% - lin reg

# setwd('C:\\Users\\doliv\\OneDrive\\2_2019_America\\2020\\20200123_MAR\\Camparison_LV_MAR\\6_methods_of_ecology_and_evolution_new_images\\paper_scripts')
# noisy_data3 = read.csv('noisy_data3.csv')
# load(file = "noisy_data3.Rdata")
# noisy_data4 = read.csv('noisy_data4.csv')
# load(file = "noisy_data4.Rdata")

# setwd('C:\\Users\\doliv\\OneDrive\\2_2019_America\\2020\\20200123_MAR\\Camparison_LV_MAR\\8_methods_of_ecology_and_evolution_2sub\\paper_figures')
# tiff("Fig_S4.tiff", height = 40, width = 20, units = 'cm',compression = "lzw", res = 300)
# Function to create a four independent color plot
par(mfcol=c(4,2))
mainLab = c('Noisy Dataset',NA,NA,NA)
for (i in 2:5)
	{
	plot(out1[,1],out1[,i],pch=20,col=colorPallet[2],ylim=c(0,3),main='',xlab='',ylab='',cex.lab=1.5,cex.axis=1.5,lwd=1)
	points(noisy_data3[,1],noisy_data3[,i],pch=1,col=colorPallet[1]) 
	#points(smoother_out_alg1$splines_est[,1],smoother_out_alg1$splines_est[,i],type='l',col=varColor[i-1])
	points(MARest1$MAR_TS_Est[,i-1],col=colorPallet[4],type='l',lty=1,lwd=3)
	points(MARest2$MAR_TS_Est[,i-1],col=colorPallet[5],type='l',lty=1,lwd=3)
	points(MARest1_smoothData_orig[,i],col='orange',type='l',lty=1,lwd=3)
	points(MARest1_smoothData_log[,i],col='khaki',type='l',lty=1,lwd=3)
	points(out_est_noisy_data3_lm[,1],out_est_noisy_data3_lm[,i],type='l',col=colorPallet[3],lwd=3)
	mtext(expression(italic('t')),side=1,line=2.1,cex=1.5)
	mtext(y_labs[i-1],side=2,line=2.1,cex=1.5)
	mtext(mainLab[i-1],side=3,line=1,cex=1.5)
	if (i == 2) {legend(35,1.25,legend = c('LV results','noisy dataset','ALVI-LR method','MAR no transform','MAR log transform','MAR no transform with smoothing','MAR log transform with smoothing'),lty = c(0,0,1,1,1,1,1),pch = c(20,1,NA,NA,NA,NA,NA),col = c(colorPallet[2],colorPallet[1],colorPallet[3],colorPallet[4],colorPallet[5],'orange','khaki'),text.col=c('black','black','black','black','black','black','black'),bty = "n",cex=1,lwd=c(NA,NA,3,3,3,3,3))}
	}
mainLab = c('Replicate Dataset',NA,NA,NA)
for (i in 2:5)
	{
	plot(out1[,1],out1[,i],pch=20,col=colorPallet[2],ylim=c(0,3),main=,xlab='',ylab='',cex.lab=1.5,cex.axis=1.5,lwd=1)
	points(noisy_data4[,1],noisy_data4[,i],pch=1,col=colorPallet[1]) 
	#points(smoother_out$splines_est[,1],smoother_out$splines_est[,i],type='l',col=varColor[i-1])
	points(MARest3$MAR_TS_Est[,i-1],col=colorPallet[4],type='l',lty=1,lwd=3)
	points(MARest4$MAR_TS_Est[,i-1],col=colorPallet[5],type='l',lty=1,lwd=3)
	points(MARest3_smoothData_orig[,i],col='orange',type='l',lty=1,lwd=3)
	points(MARest4_smoothData_log[,i],col='khaki',type='l',lty=1,lwd=3)
	points(out_est_noisy_data4_lm[,1],out_est_noisy_data4_lm[,i],type='l',col=colorPallet[3],lwd=3)
	mtext(expression(italic('t')),side=1,line=2.1,cex=1.5)
	mtext(y_labs[i-1],side=2,line=2.1,cex=1.5)
	mtext(mainLab[i-1],side=3,line=1,cex=1.5)
	if (i == 2) {legend(35,1.25,legend = c('LV results','replicate dataset','ALVI-LR method','MAR no transform','MAR log transform','MAR no transform with smoothing','MAR log transform with smoothing'),lty = c(0,0,1,1,1,1,1),pch = c(20,1,NA,NA,NA,NA,NA),col = c(colorPallet[2],colorPallet[1],colorPallet[3],colorPallet[4],colorPallet[5],'orange','khaki'),text.col=c('black','black','black','black','black','black','black'),bty = "n",cex=1,lwd=c(NA,NA,3,3,3,3,3))}
	}
# dev.off()

##### figure S4 - noisy and clustered data 20% - lin reg
#############################################################################















##########################################################################################
##### par estimation - alg for original data - fig S9

smoother_out_original = Smoother(
				dataset = out1[,1:5],
				draw = TRUE,
				data1_spline2 = 2, 
				smooth_type = 1,
				df = c(100,100,100,100),	
				log_spline = FALSE
				)



( est_alg = LV_pars_finder(
				smooth_out = smoother_out_original,
				alg1_lm2 = 1, 
				data_sample_alg = c(5,15,25,35,45)	
				) )

plot(1,1,col='white')
legend(1,1,legend = c('data','splines','slopes'),pch = c(20,NA,NA),lty = c(0,1,1),col = c('grey','green','lightgreen'),bty = "n")
cbind(est_alg)
est1_mat_alg = Format_pars(truePars = est_alg)
cbind(truePars1,est_alg,absDif = abs(truePars1-est_alg))
splineStart = smoother_out_original$splines_est[1,2:5];names(splineStart) = c('x1','x2','x3','x4')
out_est_alg = solveLV(times = t, pars = est1_mat_alg, initState = splineStart, equations = Equations)


( est_lm = LV_pars_finder(
				smooth_out = smoother_out_original,
				alg1_lm2 = 2, 
				data_sample_alg = c(5,15,25,35,45)	
				) )

plot(1,1,col='white')
legend(1,1,legend = c('data','splines','slopes'),pch = c(20,NA,NA),lty = c(0,1,1),col = c('grey','green','lightgreen'),bty = "n")
cbind(est_lm)
est1_mat_lm = Format_pars(truePars = est_lm)
cbind(truePars1,est_lm,absDif = abs(truePars1-est_lm))
splineStart = smoother_out_original$splines_est[1,2:5];names(splineStart) = c('x1','x2','x3','x4')
out_est_lm = solveLV(times = t, pars = est1_mat_lm, initState = splineStart, equations = Equations)


# setwd('C:\\Users\\doliv\\OneDrive\\2_2019_America\\2020\\20200123_MAR\\Camparison_LV_MAR\\8_methods_of_ecology_and_evolution_2sub\\paper_figures')
# tiff("Fig_S9.tiff", height = 40, width = 20, units = 'cm',compression = "lzw", res = 300)
# Function to create a four independent color plot
par(mfcol=c(4,2))
mainLab = c('ALVI-MI',NA,NA,NA)
for (i in 2:5)
	{
	plot(out1[,1],out1[,i],col=colorPallet[1],ylim=c(0,3),main='',xlab='',ylab='',cex.lab=1.5,cex.axis=1.5,lwd=1)
	points(out_est_alg[,1],out_est_alg[,i],type='l',col=colorPallet[3],lwd=2)
	mtext(expression(italic('t')),side=1,line=2.1,cex=1.5)
	mtext(y_labs[i-1],side=2,line=2.1,cex=1.5)
	mtext(mainLab[i-1],side=3,line=1,cex=1.5)
	if (i == 2) {legend(20,1,legend = c('LV results','ALVI-MI method'),lty = c(0,1),pch = c(1,NA),col = c(colorPallet[1],colorPallet[3]),text.col=c('black','black'),bty = "n",cex=1.5,lwd=c(NA,2))}
	}
mainLab = c('ALVI_LR',NA,NA,NA)
for (i in 2:5)
	{
	plot(out1[,1],out1[,i],col=colorPallet[1],ylim=c(0,3),main='',xlab='',ylab='',cex.lab=1.5,cex.axis=1.5,lwd=1)
	points(out_est_lm[,1],out_est_lm[,i],type='l',col=colorPallet[3],lwd=2)
	mtext(expression(italic('t')),side=1,line=2.1,cex=1.5)
	mtext(y_labs[i-1],side=2,line=2.1,cex=1.5)
	mtext(mainLab[i-1],side=3,line=1,cex=1.5)
	if (i == 2) {legend(20,1,legend = c('LV results','ALVI-LR method'),lty = c(0,1),pch = c(1,NA),col = c(colorPallet[1],colorPallet[3]),text.col=c('black','black'),bty = "n",cex=1.5,lwd=c(NA,2))}
	}
# dev.off()


# SSE original data againts alg solution 
sum( ( out1[,1:5] - out_est_alg[,1:5] )^2 )

# SSE original data againts lm solution 
sum( ( out1[,1:5] - out_est_lm[,1:5] )^2 )

##### par estimation - alg for original data - fig S9
############################################################################################################################







