

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
		add_later = vector()
		for (i in 1:length(data_time))
			{
			for (ii in 1:length(time_splines))
				{
				if (all.equal(data_time[i],time_splines[ii]) == TRUE)
					{
					time_splines[ii] = data_time[i]
					} else 
					{
					add_later = c(add_later,data_time[i])
					}
				}
			}
		time_splines = unique(sort(c(add_later,time_splines)))
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

Pars = truePars1 = c(		
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
for (i in 2:5)
	{
	plot(out1[,1],out1[,i],pch=20,col='grey',lwd=3,xlab='t',ylab=c('t','X1','X2','X3','X4')[i])
	}
legend(60,1,legend = c('true model'),pch = c(20),col = c('grey'),bty = "n")	# ,lty = c(0,1)

##### truedata
######################################################################################################################################



######################################################################################################################################
##### Var colors and names
varColor_light = c('lightblue','orange','lightgreen','grey')
varColor = c('blue','orange3','darkgreen','black')
varNames = c('X1','X2','X3','X4')
colorPallet = c('black','grey','blue','darkgreen','green')
##### Var colors and names
######################################################################################################################################







######################################################################################################################################
##### noisyData5

set.seed(1)
rnorm(1)
t = seq(1,100,.1)

d_output1 = c(t=min(t),state1) # initial state and starting responce matrix
for ( i in seq(min(t)+(t[2]-t[1]),max(t),t[2]-t[1]) )	# cycle to all timepoints
	{
	# retriving the variable values from the responce matrix
	if ( is.null(dim(d_output1))==TRUE )	# If is the first step
		{
		X1 = d_output1[2]
		X2 = d_output1[3]
		X3 = d_output1[4]
		X4 = d_output1[5]
		} else 					# if not first step
		{
		X1 = d_output1[dim(d_output1)[1],2]
		X2 = d_output1[dim(d_output1)[1],3]
		X3 = d_output1[dim(d_output1)[1],4]
		X4 = d_output1[dim(d_output1)[1],5]
		}

	Xs = with(as.list(c(c(X1,X2,X3,X4), Pars)), 
            {

		# discrete equations
		X1_temp = ( X1 * (a1 + b11*X1 + b12*X2 + b13*X3 + b14*X4 ) * (t[2]-t[1]) + X1 ) * rnorm(1,1,.005)
		X2_temp = ( X2 * (a2 + b21*X1 + b22*X2 + b23*X3 + b24*X4 ) * (t[2]-t[1]) + X2 ) * rnorm(1,1,.005)
		X3_temp = ( X3 * (a3 + b31*X1 + b32*X2 + b33*X3 + b34*X4 ) * (t[2]-t[1]) + X3 ) * rnorm(1,1,.005)
		X4_temp = ( X4 * (a4 + b41*X1 + b42*X2 + b43*X3 + b44*X4 ) * (t[2]-t[1]) + X4 ) * rnorm(1,1,.005)

            return(list(c(X1_temp,X2_temp,X3_temp,X4_temp)))
            })
	# updating the responce matrix
	d_output1 = rbind(d_output1,c(t=i,X1=Xs[[1]][1],X2=Xs[[1]][2],X3=Xs[[1]][3],X4=Xs[[1]][4]))
	}

noisyData5 = d_output1[c(1,sort(sample(seq(11,991,10),39))),1:5]
rownames(noisyData5) = NULL
# saveRDS(d_output1, file = "d_output1.RData")
# saveRDS(noisyData5, file = "noisyData5.RData")
# setwd('C://Users//dolivenca3//OneDrive//2_2019_America//2020//20200123_MAR//Camparison_LV_MAR//11_Frontiers//V3//paper_scripts//Fig_1_S1_S2_S3_S5_S6')
d_output1 = readRDS(file = "d_output1.RData")
noisyData5 = readRDS(file = "noisyData5.RData")

# visualization
par(mfrow=c(2,2))
for (i in 2:5)
	{
	plot(out1[,1],out1[,i],pch=20,col='grey',lwd=3,xlab='t',ylab=c('t','X1','X2','X3','X4')[i])
	lines(d_output1[,1],d_output1[,i])
	points(noisyData5[,1],noisyData5[,i],pch=16,col='lightblue')
	}
legend(45,1,legend = c('dif. eq. model','disc. model','noisyData5'),pch = c(20,20,20),col = c('grey','black','lightblue'),bty = "n")	# ,lty = c(0,1)

##### noisyData5
######################################################################################################################################



#################################################################################################
##### noisyData5 - alg1 solution

smoother_out = Smoother(
				dataset = noisyData5,
				draw = TRUE,
				data1_spline2 = 2, 
				smooth_type = 1,
				df = c(10,8,11,15), 
				log_spline = TRUE
				)

# uncomment the next lines to find a subsample for this dataset
#AlgPointFind(
#		smoother_out_APF = smoother_out,
#		dif_equations = Equations,
#		matrixEq = TRUE,
#		randomSearchTries = 100
#		) 
# "100 % complete || DF = 6 6 10 6  ||  3 9 13 15 17 1.031  ||  2 7 11 12 34 0.347"
# "100 % complete || DF = 8 8 11 15  ||  3 8 11 13 36 0.691  ||  2 7 11 12 36 0.361" 
# "100 % complete || DF = 11 8 11 15  ||  3 7 10 12 25 0.63  ||  3 8 10 12 25 0.452"
# "100 % complete || DF = 10 8 11 15  ||  3 8 10 14 25 0.629  ||  3 7 10 12 25 0.435"# "100 % complete || DF = 10 8 10 15  ||  4 6 11 14 21 0.965  ||  3 4 8 11 21 0.698"
# "100 % complete || DF = 15 8 11 15  ||  3 10 18 19 38 0.73  ||  2 5 8 12 36 0.552"
# saveRDS(globalStore, file = "globalStore_noisyData5_10_8_11_15") 
# globalStore_noisyData5 = readRDS(file = "globalStore_noisyData5_10_8_11_15")

( estND5 = LV_pars_finder(
				smooth_out = smoother_out,
				alg1_lm2 = 1, 
				data_sample_alg = c(3,8,10,14,25) # c(3, 9, 13, 15, 17)	# data_sample_lm = NULL for lm2	# data_sample_alg = 5*c(2,4,5,6,9) for alg1
				) )
cbind(estND5)
estND5_mat = Format_pars(truePars = estND5)
t = seq(1,100,1)
cbind(truePars1,estND5,absDif = abs(truePars1-estND5))
splineStart_ND5 = smoother_out$splines_est[1,2:5];names(splineStart_ND5) = c('x1','x2','x3','x4')
out_est_noisyData5_alg = solveLV(times = t, pars = estND5_mat, initState = splineStart_ND5, equations = Equations)

par(mfrow=c(2,2))
for (i in 2:5)
	{
	plot(out1[,1],out1[,i],pch=20,col='grey',lwd=3,xlab='t',ylab=c('t','X1','X2','X3','X4')[i])
	lines(d_output1[,1],d_output1[,i])
	points(noisyData5[,1],noisyData5[,i],pch=16,col='lightblue')
	points(smoother_out$splines_est[,1],smoother_out$splines_est[,i],type='l',col='lightgreen')
	points(out_est_noisyData5_alg[,1],out_est_noisyData5_alg[,i],type='l',col='blue',lwd=3)
	}
legend(45,1,legend = c('dif. eq. model','disc. model','noisyData5'),pch = c(20,20,20),col = c('grey','black','blue'),bty = "n")	# ,lty = c(0,1)


# SSE original data againts alg solution (this is in table one)
sum( ( out1[,1:5] - out_est_noisyData5_alg[,1:5])^2 )

# SSE noisy data againts alg solution (this is to justify why alg pars est are not similar to true pars)
sum( (noisyData5[,1:5] - out_est_noisyData5_alg[which(round(out_est_noisyData5_alg[,1],1) %in% round(noisyData5[,1],1)),1:5] )^2 )
# SSE spline against alg solution 
sum( (smoother_out$splines_est[which(round(smoother_out$splines_est[,1],5) %in% round(noisyData5[,1],1)),1:5] - out_est_noisyData5_alg[which(round(out_est_noisyData5_alg[,1],1) %in% round(noisyData5[,1],1)),1:5] )^2 )

# SSE of original parameter against parameter estimates
sum( ( Pars - estND5 )^2 )

##### noisyData5 - alg1 solution
#################################################################################################



############################################################################################################################
##### MAR for noisyData5

# Completing the data to have the same time structure as the original data
t=seq(1,100,1)
TS_noisy_data = cbind(t,matrix(rep(NA,length(t)*4),ncol=4))
TS_noisy_data[which(round(t,1) %in% round(noisyData5[,1],1)),] = noisyData5

# MAR estimates without data log transformation
MARest1 =  fastMAR(
			data_vars = TS_noisy_data[,2:5],
			log_transform=FALSE,
			demeaned=FALSE,
			abstol_=0.01,
			maxit_ = 2000,
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
	}

# put time back in the data
MARest1_orig = cbind(1:dim(MARest1$MAR_TS_Est)[1],MARest1$MAR_TS_Est)
MARest1_log = cbind(1:dim(MARest2$MAR_TS_Est)[1],MARest2$MAR_TS_Est)

# errors
sum( ( out1[,1:5] - MARest1_orig[,1:5] )^2 )
sum( ( out1[,1:5] - MARest1_log[,1:5] )^2 )

##### MAR for noisyData5
###############################################################################################################



############################################################################################################################
##### MAR with smooth data for noisyData5

# smoothing the data for MAR
smoother_for_MAR_out = Smoother_for_MAR(	
						dataset = noisyData5,
						df = c(8,8,11,15),
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
	}

# put time back in the data
MARest1_smoothData_orig = cbind(1:dim(MARest1_smoothData$MAR_TS_Est)[1],MARest1_smoothData$MAR_TS_Est)
MARest1_smoothData_log = cbind(1:dim(MARest2_smoothData$MAR_TS_Est)[1],MARest2_smoothData$MAR_TS_Est)

# errors
sum( ( out1[,2:5] - MARest1_smoothData_orig[,2:5] )^2 )
sum( ( out1[,2:5] - MARest1_smoothData_log[,2:5] )^2 )

##### MAR with smooth data for noisyData5
############################################################################################################################







######################################################################################################################################
##### repData5 replicateData5

set.seed(1)
rnorm(1)

Fig2_data = function(t=seq(1,100,.1),init_state)
	{
	d_output1 = c(t=min(t),init_state) # initial state and starting responce matrix
	for ( i in seq(min(t)+(t[2]-t[1]),max(t),t[2]-t[1]) )	# cycle to all timepoints
		{
		# retriving the variable values from the responce matrix
		if ( is.null(dim(d_output1))==TRUE )	# If is the first step
			{
			X1 = d_output1[2]
			X2 = d_output1[3]
			X3 = d_output1[4]
			X4 = d_output1[5]	
			} else 					# if not first step
			{
			X1 = d_output1[dim(d_output1)[1],2]
			X2 = d_output1[dim(d_output1)[1],3]
			X3 = d_output1[dim(d_output1)[1],4]
			X4 = d_output1[dim(d_output1)[1],5]
			}
	
		Xs = with(as.list(c(c(X1,X2,X3,X4), Pars)), 
      	      {
	
			# discrete equations
			X1_temp = ( X1 * (a1 + b11*X1 + b12*X2 + b13*X3 + b14*X4 ) * (t[2]-t[1]) + X1 ) * rnorm(1,1,.005)
			X2_temp = ( X2 * (a2 + b21*X1 + b22*X2 + b23*X3 + b24*X4 ) * (t[2]-t[1]) + X2 ) * rnorm(1,1,.005)
			X3_temp = ( X3 * (a3 + b31*X1 + b32*X2 + b33*X3 + b34*X4 ) * (t[2]-t[1]) + X3 ) * rnorm(1,1,.005)
			X4_temp = ( X4 * (a4 + b41*X1 + b42*X2 + b43*X3 + b44*X4 ) * (t[2]-t[1]) + X4 ) * rnorm(1,1,.005)	

	            return(list(c(X1_temp,X2_temp,X3_temp,X4_temp)))
	            })
		# updating the responce matrix
		d_output1 = rbind(d_output1,c(t=i,X1=Xs[[1]][1],X2=Xs[[1]][2],X3=Xs[[1]][3],X4=Xs[[1]][4]))
		}
	rownames(d_output1) = NULL
	return(d_output1)
	}

d_output1 = Fig2_data(t=seq(1,100,.1),init_state=state1)
d_output2 = Fig2_data(t=seq(1,100,.1),init_state=state1)
d_output3 = Fig2_data(t=seq(1,100,.1),init_state=state1)
d_output4 = Fig2_data(t=seq(1,100,.1),init_state=state1)
d_output5 = Fig2_data(t=seq(1,100,.1),init_state=state1)

repData5 = d_output1[1,]
repData5_mean = d_output1[1,]
for (i in c(1,51,101,171,251,341,401,501,591,691,791))
	{
	repData5 = rbind(repData5,d_output1[i,1:5])
	repData5 = rbind(repData5,d_output2[i,1:5])
	repData5 = rbind(repData5,d_output3[i,1:5])
	repData5 = rbind(repData5,d_output4[i,1:5])
	repData5 = rbind(repData5,d_output5[i,1:5])

 	temp1 = colMeans(t(data.frame(d_output1[i,1:5],
			d_output2[i,1:5],
			d_output3[i,1:5],
			d_output4[i,1:5],
			d_output5[i,1:5]
			)))
			
	repData5_mean = rbind(repData5_mean,temp1)
	}
rownames(repData5) = NULL
repData5 = repData5[2:dim(repData5)[1],]
rownames(repData5_mean) = NULL
repData5_mean = repData5_mean[2:dim(repData5_mean)[1],]

# saveRDS(list(d_output1,d_output2,d_output3,d_output4,d_output5), file = "List_output_ND6.RData")
# saveRDS(repData5, file = "repData5.RData")
# saveRDS(repData5_mean, file = "repData5_mean.RData")
# setwd("C://Users//dolivenca3//OneDrive//2_2019_America//2020//20200123_MAR//Camparison_LV_MAR//11_Frontiers//V3//paper_scripts//Fig_1_S1_S2_S3_S4_S5_S9")
List_output_rep5Data = readRDS(file = "List_output_rep5Data.RData")
repData5 = readRDS(file = "repData5.RData")
repData5_mean = readRDS(file = "repData5_mean.RData")

# visualization
par(mfrow=c(2,2))
for (i in 2:5)
	{
	plot(out1[,1],out1[,i],pch=20,col='grey',lwd=3,xlab='t',ylab=c('t','X1','X2','X3','X4')[i])
	lines(List_output_rep5Data[[1]][,1],List_output_rep5Data[[1]][,i])
	lines(List_output_rep5Data[[2]][,1],List_output_rep5Data[[2]][,i])
	lines(List_output_rep5Data[[3]][,1],List_output_rep5Data[[3]][,i])
	lines(List_output_rep5Data[[4]][,1],List_output_rep5Data[[4]][,i])
	lines(List_output_rep5Data[[5]][,1],List_output_rep5Data[[5]][,i])
	points(repData5[,1],repData5[,i],pch=16,col='lightblue')
	}
legend(45,1,legend = c('dif. eq. model','disc. model','noisyData5'),pch = c(20,20,20),col = c('grey','black','lightblue'),bty = "n")	# ,lty = c(0,1)

##### repData5
######################################################################################################################################



#################################################################################################
##### repData5 - alg1 solution

smoother_out = Smoother(
				dataset = repData5_mean,
				draw = TRUE,
				data1_spline2 = 2, 
				smooth_type = 1,
				df = c(8,10,10,9),	
				log_spline = TRUE
				)

# uncomment the next lines to find a subsample for this dataset
#AlgPointFind(
#		smoother_out_APF = smoother_out,
#		dif_equations = Equations,
#		matrixEq = TRUE,
#		randomSearchTries = 100
#		) 
# "100 % complete || DF = 11 11 11 11  ||  2 3 4 6 9 0.097  ||  2 3 4 6 9 0.097"
# "100 % complete || DF = 8 10 8 9  ||  2 5 6 7 10 0.058  ||  2 4 6 7 10 0.024"
# "100 % complete || DF = 8 10 9 9  ||  2 3 5 6 10 0.039  ||  2 3 4 6 7 0.027"
# "100 % complete || DF = 8 10 10 9  ||  2 3 5 6 10 0.034  ||  2 3 4 6 7 0.03"
# saveRDS(globalStore, file = "globalStore_repData5_8_10_10_9") 
# globalStore_repData5 = readRDS(file = "globalStore_repData5_8_10_10_9")

( est_repData5 = LV_pars_finder(
				smooth_out = smoother_out,
				alg1_lm2 = 1, 
				data_sample_alg = c(2,3,4,6,7)	# data_sample_lm = NULL for lm2	# data_sample_alg = 5*c(2,4,5,6,9) for alg1
				) )
cbind(est_repData5)
est_repData5_mat = Format_pars(truePars = est_repData5)
t = seq(1,100,1)
cbind(truePars1,est_repData5,absDif = abs(truePars1-est_repData5))
splineStart_repData5 = smoother_out$splines_est[1,2:5];names(splineStart_repData5) = c('x1','x2','x3','x4')
out_est_repData5_alg = solveLV(times = t, pars = est_repData5_mat, initState = splineStart_repData5, equations = Equations)

par(mfrow=c(2,2))
for (i in 2:5)
	{
	plot(out1[,1],out1[,i],pch=20,col='grey',lwd=3,xlab='t',ylab=c('t','X1','X2','X3','X4')[i])
	points(repData5[,1],repData5[,i],pch=16,col='lightblue')
	points(smoother_out$splines_est[,1],smoother_out$splines_est[,i],type='l',col='lightgreen')
	points(out_est_repData5_alg[,1],out_est_repData5_alg[,i],type='l',col='blue',lwd=3)
	}
legend(45,1,legend = c('dif. eq. model','disc. model','noisyData5'),pch = c(20,20,20),col = c('grey','black','blue'),bty = "n")	# ,lty = c(0,1)


# SSE original data againts alg solution (this is in table one)
sum( ( out1[,1:5] - out_est_repData5_alg[,1:5] )^2 )

# SSE spline against alg solution 
sum( (smoother_out$splines_est[which(round(smoother_out$splines_est[,1],5) %in% round(out_est_repData5_alg[,1],1)),1:5] - out_est_repData5_alg[which(round(out_est_repData5_alg[,1],1) %in% round(repData5[,1],1)),1:5] )^2 )

# SSE of original parameter against parameter estimates
sum( ( Pars - est_repData5 )^2 )

##### repData5 - alg1 solution
#################################################################################################



############################################################################################################################
##### MAR for repData5

# Completing the data to have the same time structure as the original data
t = seq(1,100,1)
TS_noisy_data = cbind(t,matrix(rep(NA,length(t)*4),ncol=4))
TS_noisy_data[which(round(t,1) %in% round(repData5[,1],1)),] = repData5_mean

# MAR estimates without data log transformation
MARest3 =  fastMAR(
			data_vars = TS_noisy_data[,2:5],
			log_transform=FALSE,
			demeaned=FALSE,
			abstol_=0.01,
			maxit_ = 1500,
			includeNoise = FALSE,
			estInit = FALSE,
			tSteps_ = NULL,
			nsim_ = 100
			) 
MARest3$MARSSci_output

# MAR estimates with data log transformation
MARest4 =  fastMAR(
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
MARest4$MARSSci_output

# plots
par(mfrow=c(2,2))
for ( i in 1:(dim(TS_noisy_data)[2]-1) )
	{
	plot(out1[,1],out1[,i+1],type='l',col='lightgrey',lwd=3,ylim=c(0,3),ylab = varNames[i])

	# for (ii in 1:100) { points(MARest3$MARSSsim_data[i,,ii],type='l',col='grey') }
	points(meanMaker(MARest3)[i,],type='l',col = 'purple')

	# for (ii in 1:100) { points(MARest4$MARSSsim_data[i,,ii],type='l',col='grey') }
	points(meanMaker(MARest4)[i,],type='l',col = 'purple',lty=2)


	points(out1[,1],out1[,i+1],type='l',col='lightgrey',lwd=3,ylim=c(0,3),ylab = varNames[i])
	points(TS_noisy_data[,i+1],pch=20,col=varColor[i],xlab = 'n')
	points(MARest3$MAR_TS_Est[,i],col='green',type='l',lty=1)
	points(MARest4$MAR_TS_Est[,i],col='green',type='l',lty=2)
	if (i==1) {legend(40,1.8,legend=c('original data','noisy data','MAR no transform','MAR log transform'),pch=c(NA,20,NA,NA),lty=c(1,NA,1,2),col=c('lightgrey',varColor[i],varColor_light[i],varColor_light[i]),bty = "n",cex=1,lwd=c(3,NA,1,1))}
	}

# put time back in the data
MARest2_orig = cbind(1:dim(MARest3$MAR_TS_Est)[1],MARest3$MAR_TS_Est)
MARest2_log = cbind(1:dim(MARest4$MAR_TS_Est)[1],MARest4$MAR_TS_Est)

# errors
sum( ( out1[,2:5] - MARest2_orig[,2:5] )^2 )
sum( ( out1[,2:5] - MARest2_log[,2:5] )^2 )

##### MAR for repData5
###############################################################################################################



############################################################################################################################
##### MAR with smooth data for repData5

# smoothing the data for MAR
smoother_for_MAR_out = Smoother_for_MAR(	
						dataset = repData5,
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

	# for (ii in 1:100) { points(MARest3$MARSSsim_data[i,,ii],type='l',col='grey') }			# MAR estimates without data log transformation, with spline smooth data, with noise
	points(meanMaker(MARest3_smoothData)[i,],type='l',col = 'purple')

	# for (ii in 1:100) { points(MARest4$MARSSsim_data[i,,ii],type='l',col='grey') }
	points(meanMaker(MARest4_smoothData)[i,],type='l',col = 'purple',lty=2)					# MAR estimates with data log transformation, with spline smooth data, with noise

	points(TS_smooth_data[,i+1],pch=20,col=varColor[i],xlab = 'n')						# treated data
	points(MARest3_smoothData$MAR_TS_Est[,i],col='green',type='l',lty=1)						# MAR estimates without data log transformation, with spline smooth data
	points(MARest4_smoothData$MAR_TS_Est[,i],col='green',type='l',lty=2)						# MAR estimates with data log transformation, with spline smooth data
	if (i==1) {legend(40,1.8,legend=c('original data','noisy data','MAR no transform','MAR log transform'),pch=c(NA,20,NA,NA),lty=c(1,NA,1,2),col=c('lightgrey',varColor[i],varColor_light[i],varColor_light[i]),bty = "n",cex=1,lwd=c(3,NA,1,1))}
	}

# put time back in the data
MARest2_smoothData_orig = cbind(1:dim(MARest3_smoothData$MAR_TS_Est)[1],MARest3_smoothData$MAR_TS_Est)
MARest2_smoothData_log = cbind(1:dim(MARest4_smoothData$MAR_TS_Est)[1],MARest4_smoothData$MAR_TS_Est)

# errors
sum( ( out1[,2:5] - MARest2_smoothData_orig[,2:5] )^2 )
sum( ( out1[,2:5] - MARest2_smoothData_log[,2:5] )^2 )

##### MAR with smooth data for repData5
############################################################################################################################



#############################################################################
##### figure 1 - noisy and clustered data 20% - alg1

# data																					
# setwd("C://Users//dolivenca3//OneDrive//2_2019_America//2020//20200123_MAR//Camparison_LV_MAR//11_Frontiers//V3//paper_scripts//Fig_1_S1_S2_S3_S5_S6")
noisyData5 = readRDS(file = "noisyData5.RData")
repData5 = readRDS(file = "repData5.RData")
repData5_mean = readRDS(file = "repData5_mean.RData")

# setwd("C:\\Users\\dolivenca3\\OneDrive\\2_2019_America\\2020\\20200123_MAR\\Camparison_LV_MAR\\11_Frontiers\\V3\\Figures")
# tiff("Fig_1_.tiff", height = 40, width = 20, units = 'cm',compression = "lzw", res = 300)
# Function to create a four independent color plot
par(mfcol=c(4,2))
mainLab = c('Noisy Dataset',NA,NA,NA)
y_labs = c(expression(italic('X'[1])),expression(italic('X'[2])),expression(italic('X'[3])),expression(italic('X'[4])))
for (i in 2:5)
	{
	plot(out1[,1],out1[,i],pch=20,col=colorPallet[2],ylim=c(0,3),main='',xlab='',ylab='',cex.lab=1.5,cex.axis=1.5,lwd=1)
	points(noisyData5[,1],noisyData5[,i],pch=1,col=colorPallet[1]) 
	#points(smoother_out_alg1$splines_est[,1],smoother_out_alg1$splines_est[,i],type='l',col=varColor[i-1])
	points(MARest1_orig[,i],col=colorPallet[4],type='l',lty=1,lwd=3)
	points(MARest1_log[,i],col=colorPallet[5],type='l',lty=1,lwd=3)
	points(MARest1_smoothData_orig[,i],col='orange',type='l',lty=1,lwd=3)
	points(MARest1_smoothData_log[,i],col='khaki',type='l',lty=1,lwd=3)	
	points(out_est_noisyData5_alg[,1],out_est_noisyData5_alg[,i],type='l',col=colorPallet[3],lwd=3)
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
	points(repData5[,1],repData5[,i],pch=1,col=colorPallet[1]) 
	#points(smoother_out$splines_est[,1],smoother_out$splines_est[,i],type='l',col=varColor[i-1])
	points(MARest2_orig[,i],col=colorPallet[4],type='l',lty=1,lwd=3)
	points(MARest2_log[,i],col=colorPallet[5],type='l',lty=1,lwd=3)
	points(MARest2_smoothData_orig[,i],col='orange',type='l',lty=1,lwd=3)
	points(MARest2_smoothData_log[,i],col='khaki',type='l',lty=1,lwd=3)
	points(out_est_repData5_alg[,1],out_est_repData5_alg[,i],type='l',col=colorPallet[3],lwd=3)
	mtext(expression(italic('t')),side=1,line=2.1,cex=1.5)
	mtext(y_labs[i-1],side=2,line=2.1,cex=1.5)
	mtext(mainLab[i-1],side=3,line=1,cex=1.5)
	if (i == 2) {legend(35,1.25,legend = c('LV results','replicate dataset','ALVI-MI method','MAR no transform','MAR log transform','MAR no transform with smoothing','MAR log transform with smoothing'),lty = c(0,0,1,1,1,1,1),pch = c(20,1,NA,NA,NA,NA,NA),col = c(colorPallet[2],colorPallet[1],colorPallet[3],colorPallet[4],colorPallet[5],'orange','khaki'),text.col=c('black','black','black','black','black','black','black'),bty = "n",cex=1,lwd=c(NA,NA,3,3,3,3,3))}
	}
# dev.off()


x=1:5
error = cbind(
		rbind(sum( ( out1[,x] - out_est_noisyData5_alg[,x])^2 ),
			sum( ( out1[,x] - MARest1_orig[,x] )^2 ),
			sum( ( out1[,x] - MARest1_log[,x] )^2 ),
			sum( ( out1[,x] - MARest1_smoothData_orig[,x] )^2 ),
			sum( ( out1[,x] - MARest1_smoothData_log[,x] )^2 )),
		rbind(sum( ( out1[,x] - out_est_repData5_alg[,x])^2 ),
			sum( ( out1[,x] - MARest2_orig[,x] )^2 ),
			sum( ( out1[,x] - MARest2_log[,x] )^2 ),
			sum( ( out1[,x] - MARest2_smoothData_orig[,x] )^2 ),
			sum( ( out1[,x] - MARest2_smoothData_log[,x] )^2 ))
		)
rownames(error) = c('ALVI-MI method','MAR no transform','MAR log transform','MAR no transform with smoothing','MAR log transform with smoothing')
colnames(error) = c('noisyData5','repData5')
error

error_detailed = cbind(
		rbind(sum( ( out1[,2] - out_est_noisyData5_alg[,2])^2 ),
			sum( ( out1[,2] - MARest1_orig[,2] )^2 ),
			sum( ( out1[,2] - MARest1_log[,2] )^2 ),
			sum( ( out1[,2] - MARest1_smoothData_orig[,2] )^2 ),
			sum( ( out1[,2] - MARest1_smoothData_log[,2] )^2 )),
		rbind(sum( ( out1[,3] - out_est_noisyData5_alg[,3])^2 ),
			sum( ( out1[,3] - MARest1_orig[,3] )^2 ),
			sum( ( out1[,3] - MARest1_log[,3] )^2 ),
			sum( ( out1[,3] - MARest1_smoothData_orig[,3] )^2 ),
			sum( ( out1[,3] - MARest1_smoothData_log[,3] )^2 )),
		rbind(sum( ( out1[,4] - out_est_noisyData5_alg[,4])^2 ),
			sum( ( out1[,4] - MARest1_orig[,4] )^2 ),
			sum( ( out1[,4] - MARest1_log[,4] )^2 ),
			sum( ( out1[,4] - MARest1_smoothData_orig[,4] )^2 ),
			sum( ( out1[,4] - MARest1_smoothData_log[,4] )^2 )),
		rbind(sum( ( out1[,5] - out_est_noisyData5_alg[,5])^2 ),
			sum( ( out1[,5] - MARest1_orig[,5] )^2 ),
			sum( ( out1[,5] - MARest1_log[,5] )^2 ),
			sum( ( out1[,5] - MARest1_smoothData_orig[,5] )^2 ),
			sum( ( out1[,5] - MARest1_smoothData_log[,5] )^2 )),

		rbind(sum( ( out1[,2] - out_est_repData5_alg[,2])^2 ),
			sum( ( out1[,2] - MARest2_orig[,2] )^2 ),
			sum( ( out1[,2] - MARest2_log[,2] )^2 ),
			sum( ( out1[,2] - MARest2_smoothData_orig[,2] )^2 ),
			sum( ( out1[,2] - MARest2_smoothData_log[,2] )^2 )),
		rbind(sum( ( out1[,3] - out_est_repData5_alg[,3])^2 ),
			sum( ( out1[,3] - MARest2_orig[,3] )^2 ),
			sum( ( out1[,3] - MARest2_log[,3] )^2 ),
			sum( ( out1[,3] - MARest2_smoothData_orig[,3] )^2 ),
			sum( ( out1[,3] - MARest2_smoothData_log[,3] )^2 )),
		rbind(sum( ( out1[,4] - out_est_repData5_alg[,4])^2 ),
			sum( ( out1[,4] - MARest2_orig[,4] )^2 ),
			sum( ( out1[,4] - MARest2_log[,4] )^2 ),
			sum( ( out1[,4] - MARest2_smoothData_orig[,4] )^2 ),
			sum( ( out1[,4] - MARest2_smoothData_log[,4] )^2 )),
		rbind(sum( ( out1[,5] - out_est_repData5_alg[,5])^2 ),
			sum( ( out1[,5] - MARest2_orig[,5] )^2 ),
			sum( ( out1[,5] - MARest2_log[,5] )^2 ),
			sum( ( out1[,5] - MARest2_smoothData_orig[,5] )^2 ),
			sum( ( out1[,5] - MARest2_smoothData_log[,5] )^2 ))
		)
rownames(error_detailed) = c('ALVI-MI','MAR','MAR log transform','MAR with smoothing','MAR log transform with smoothing')
colnames(error_detailed) = c('noisyData5 - X1','noisyData5 - X2','noisyData5 - X3','noisyData5 - X4','repData5 - X1','repData5 - X2','repData5 - X3','repData5 - X4')
error_detailed

##### figure 1 - noisy and clustered data 20% - alg1
#############################################################################



#############################################################################
##### figure S1 - noisy and clustered data

# data
# setwd("C://Users//dolivenca3//OneDrive//2_2019_America//2020//20200123_MAR//Camparison_LV_MAR//11_Frontiers//V3//paper_scripts//Fig_1_S1_S2_S3_S5_S6")
noisyData5 = readRDS(file = "noisyData5.RData")
repData5 = readRDS(file = "repData5.RData")
repData5_mean = readRDS(file = "repData5_mean.RData")

# setwd("C:\\Users\\dolivenca3\\OneDrive\\2_2019_America\\2020\\20200123_MAR\\Camparison_LV_MAR\\11_Frontiers\\V3\\Figures")
# tiff("Fig_S1.tiff", height = 40, width = 20, units = 'cm',compression = "lzw", res = 300)
# Function to create a four independent color plot
par(mfcol=c(4,2))
mainLab = c('Noisy Dataset',NA,NA,NA)
y_labs = c(expression(italic('X'[1])),expression(italic('X'[2])),expression(italic('X'[3])),expression(italic('X'[4])))
for ( i in 2:5 )
	{
	plot(out1[,1],out1[,i],type='p',pch=20,col=colorPallet[2],xlim=c(0,length(out1[,1])),ylim=c(0,3),xlab='',main='',ylab='',cex.lab=1.5,cex.axis=1.5,lwd=1)
	points(noisyData5[,1],noisyData5[,i],col=colorPallet[1])	#,pch=20
	mtext(expression(italic('t')),side=1,line=2.1,cex=1.5)
	mtext(y_labs[i-1],side=2,line=2.1,cex=1.5)
	mtext(mainLab[i-1],side=3,line=1,cex=1.5)
	}
legend(45,1,legend = c('LV results','noisy dataset'),pch = c(20,1),col = c(colorPallet[2],colorPallet[1]),bty = "n",cex=1.3)
mainLab = c('Replicate Dataset',NA,NA,NA)
for (i in 2:5)
	{
	plot(out1[,1],out1[,i],pch=20,col=colorPallet[2],ylim=c(0,3),main='',xlab='',ylab='',cex.lab=1.5,cex.axis=1.5,lwd=1)
	points(repData5[,1],repData5[,i],col=colorPallet[1])	#,pch=20
	mtext(expression(italic('t')),side=1,line=2.1,cex=1.5)
	mtext(y_labs[i-1],side=2,line=2.1,cex=1.5)
	mtext(mainLab[i-1],side=3,line=1,cex=1.5)
	}
legend(45,1,legend = c('LV results','replicate dataset'),pch = c(20,1),col = c(colorPallet[2],colorPallet[1]),bty = "n",cex=1.3)

# dev.off()

##### figure S1 - noisy and clustered data 20%
#############################################################################



#############################################################################
##### figure S3 - smoothing on X1 

##########
### splines on X1

x = seq(head(out1[,1],1),tail(out1[,1],1),.1)
# setwd("C:\\Users\\dolivenca3\\OneDrive\\2_2019_America\\2020\\20200123_MAR\\Camparison_LV_MAR\\11_Frontiers\\V3\\Figures")
# tiff("Fig_S3.tiff", height = 40, width = 20, units = 'cm',compression = "lzw", res = 300)
par(mfcol=c(4,2))
DFlist = c(5,10,15,35)
mainLab = c('Splines',NA,NA,NA)
for (i in 1:4)
	{
	smoothSpline <- smooth.spline(noisyData5[,1],noisyData5[,2],df=DFlist[i])                         		# smoothing with degrees of freedom (df)
	smooth = predict(object = smoothSpline, x = x, deriv = 0)                                 			# get the points of the fit of that linear model
	f_of_x <- splinefun(smooth$x,smooth$y)
	plot(out1[,1],out1[,2],pch=20,col=colorPallet[2],ylim=c(0,3),main='',xlab='',ylab='',cex.lab=1.5,cex.axis=1.5,lwd=1)
	points(noisyData5[,1],pch=1,noisyData5[,2],col=colorPallet[1]) 
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
				dataset = noisyData5,
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
	points(noisyData5[,1],pch=1,noisyData5[,2],col=colorPallet[1]) 
	legend(30,1,legend = c('LV results','noisy data',paste('LOESS span =',spanList[i])),lty = c(0,0,1),pch = c(20,1,NA),col = c(colorPallet[2],colorPallet[1],'darkgreen'),text.col=c('black','black','black'),bty = "n",cex=1,lwd=c(1,1,3))
	points(est_data1_sparse_noisy$splines_est[,1],est_data1_sparse_noisy$splines_est[,2],type='l',col='darkgreen',lwd=3)
	mtext(expression(italic('t')),side=1,line=2.1,cex=1.5)
	mtext(rep(y_labs[1],4),side=2,line=2.1,cex=1.5)
	mtext(mainLab[i],side=3,line=1,cex=1.5)
	}
# dev.off()

### LOESS on X1
##########

##### figure S3 - smoothing on X1 
############################################################################################################



############################################################################################################
##### figure S5 - noisy and clustered data - lin reg

# data
# setwd("C://Users//dolivenca3//OneDrive//2_2019_America//2020//20200123_MAR//Camparison_LV_MAR//11_Frontiers//V3//paper_scripts//Fig_1_S1_S2_S3_S5_S6")
noisyData5 = readRDS(file = "noisyData5.RData")
repData5 = readRDS(file = "repData5.RData")
repData5_mean = readRDS(file = "repData5_mean.RData")



##### noisyData5 - lm2 solution

smoother_out = Smoother(
				dataset = noisyData5,
				draw = TRUE,
				data1_spline2 = 2, 
				smooth_type = 1,
				df = c(10,8,11,15), 
				log_spline = TRUE
				)

( estND5_lr = LV_pars_finder(
				smooth_out = smoother_out,
				alg1_lm2 = 2
				) )
cbind(estND5_lr)
estND5_lr_mat = Format_pars(truePars = estND5_lr)
t = seq(1,100,1)
cbind(truePars1,estND5_lr,absDif = abs(truePars1-estND5_lr))
splineStart_ND5_lr = smoother_out$splines_est[1,2:5];names(splineStart_ND5_lr) = c('x1','x2','x3','x4')
out_est_noisyData5_lr = solveLV(times = t, pars = estND5_lr_mat, initState = splineStart_ND5_lr, equations = Equations)

par(mfrow=c(2,2))
for (i in 2:5)
	{
	plot(out1[,1],out1[,i],pch=20,col='grey',lwd=3,xlab='t',ylab=c('t','X1','X2','X3','X4')[i])
	lines(d_output1[,1],d_output1[,i])
	points(noisyData5[,1],noisyData5[,i],pch=16,col='lightblue')
	points(smoother_out$splines_est[,1],smoother_out$splines_est[,i],type='l',col='lightgreen')
	points(out_est_noisyData5_lr[,1],out_est_noisyData5_lr[,i],type='l',col='blue',lwd=3)
	}
legend(45,1,legend = c('dif. eq. model','disc. model','noisyData5'),pch = c(20,20,20),col = c('grey','black','blue'),bty = "n")	# ,lty = c(0,1)


# SSE original data againts lr solution (this is in table one)
sum( ( out1[,1:5] - out_est_noisyData5_lr[,1:5])^2 )

# SSE noisy data againts lr solution (this is to justify why alg pars est are not similar to true pars)
sum( (noisyData5[,1:5] - out_est_noisyData5_lr[which(round(out_est_noisyData5_lr[,1],1) %in% round(noisyData5[,1],1)),1:5] )^2 )
# SSE spline against lr solution 
sum( (smoother_out$splines_est[which(round(smoother_out$splines_est[,1],5) %in% round(noisyData5[,1],1)),1:5] - out_est_noisyData5_lr[which(round(out_est_noisyData5_lr[,1],1) %in% round(noisyData5[,1],1)),1:5] )^2 )

# SSE of original parameter against parameter estimates
sum( ( Pars - estND5_lr )^2 )

##### noisyData5 - lm2 solution



##### repData5 - lm2 solution

smoother_out = Smoother(
				dataset = repData5_mean,
				draw = TRUE,
				data1_spline2 = 2, 
				smooth_type = 1,
				df = c(8,10,10,9),	
				log_spline = TRUE
				)

( est_repData5_lr = LV_pars_finder(
				smooth_out = smoother_out,
				alg1_lm2 = 2
				) )
cbind(est_repData5_lr)
est_repData5_lr_mat = Format_pars(truePars = est_repData5_lr)
t = seq(1,100,1)
cbind(truePars1,est_repData5_lr,absDif = abs(truePars1-est_repData5_lr))
splineStart_ND6_lr = smoother_out$splines_est[1,2:5];names(splineStart_ND6_lr) = c('x1','x2','x3','x4')
out_est_repData5_lr = solveLV(times = t, pars = est_repData5_lr_mat, initState = splineStart_ND6_lr, equations = Equations)

par(mfrow=c(2,2))
for (i in 2:5)
	{
	plot(out1[,1],out1[,i],pch=20,col='grey',lwd=3,xlab='t',ylab=c('t','X1','X2','X3','X4')[i])
	lines(d_output1[,1],d_output1[,i])
	points(repData5[,1],repData5[,i],pch=16,col='lightblue')
	points(smoother_out$splines_est[,1],smoother_out$splines_est[,i],type='l',col='lightgreen')
	points(out_est_repData5_lr[,1],out_est_repData5_lr[,i],type='l',col='blue',lwd=3)
	}
legend(45,1,legend = c('dif. eq. model','disc. model','noisyData5'),pch = c(20,20,20),col = c('grey','black','blue'),bty = "n")	# ,lty = c(0,1)


# SSE original data againts alg solution (this is in table one)
sum( ( out1[,1:5] - out_est_repData5_lr[,1:5] )^2 )

# SSE spline against alg solution 
sum( (smoother_out$splines_est[which(round(smoother_out$splines_est[,1],5) %in% round(out_est_repData5_lr[,1],1)),1:5] - out_est_repData5_lr[which(round(out_est_repData5_lr[,1],1) %in% round(repData5[,1],1)),1:5] )^2 )

# SSE of original parameter against parameter estimates
sum( ( Pars - est_repData5_lr )^2 )

##### repData5 - lm2 solution



# setwd("C:\\Users\\dolivenca3\\OneDrive\\2_2019_America\\2020\\20200123_MAR\\Camparison_LV_MAR\\11_Frontiers\\V3\\Figures")
# tiff("Fig_S5.tiff", height = 40, width = 20, units = 'cm',compression = "lzw", res = 300)
# Function to create a four independent color plot
par(mfcol=c(4,2))
mainLab = c('Noisy Dataset',NA,NA,NA)
for (i in 2:5)
	{
	plot(out1[,1],out1[,i],pch=20,col=colorPallet[2],ylim=c(0,3),main='',xlab='',ylab='',cex.lab=1.5,cex.axis=1.5,lwd=1)
	points(noisyData5[,1],noisyData5[,i],pch=1,col=colorPallet[1]) 
	#points(smoother_out_alg1$splines_est[,1],smoother_out_alg1$splines_est[,i],type='l',col=varColor[i-1])
	points(MARest1_orig[,i],col=colorPallet[4],type='l',lty=1,lwd=3)
	points(MARest1_log[,i],col=colorPallet[5],type='l',lty=1,lwd=3)
	points(MARest1_smoothData_orig[,i],col='orange',type='l',lty=1,lwd=3)
	points(MARest1_smoothData_log[,i],col='khaki',type='l',lty=1,lwd=3)
	points(out_est_noisyData5_lr[,1],out_est_noisyData5_lr[,i],type='l',col=colorPallet[3],lwd=3)
	mtext(expression(italic('t')),side=1,line=2.1,cex=1.5)
	mtext(y_labs[i-1],side=2,line=2.1,cex=1.5)
	mtext(mainLab[i-1],side=3,line=1,cex=1.5)
	if (i == 2) {legend(35,1.25,legend = c('LV results','noisy dataset','ALVI-LR method','MAR no transform','MAR log transform','MAR no transform with smoothing','MAR log transform with smoothing'),lty = c(0,0,1,1,1,1,1),pch = c(20,1,NA,NA,NA,NA,NA),col = c(colorPallet[2],colorPallet[1],colorPallet[3],colorPallet[4],colorPallet[5],'orange','khaki'),text.col=c('black','black','black','black','black','black','black'),bty = "n",cex=1,lwd=c(NA,NA,3,3,3,3,3))}
	}
mainLab = c('Replicate Dataset',NA,NA,NA)
for (i in 2:5)
	{
	plot(out1[,1],out1[,i],pch=20,col=colorPallet[2],ylim=c(0,3),main=,xlab='',ylab='',cex.lab=1.5,cex.axis=1.5,lwd=1)
	points(repData5[,1],repData5[,i],pch=1,col=colorPallet[1]) 
	#points(smoother_out$splines_est[,1],smoother_out$splines_est[,i],type='l',col=varColor[i-1])
	points(MARest2_orig[,i],col=colorPallet[4],type='l',lty=1,lwd=3)
	points(MARest2_log[,i],col=colorPallet[5],type='l',lty=1,lwd=3)
	points(MARest2_smoothData_orig[,i],col='orange',type='l',lty=1,lwd=3)
	points(MARest2_smoothData_log[,i],col='khaki',type='l',lty=1,lwd=3)
	points(out_est_repData5_lr[,1],out_est_repData5_lr[,i],type='l',col=colorPallet[3],lwd=3)
	mtext(expression(italic('t')),side=1,line=2.1,cex=1.5)
	mtext(y_labs[i-1],side=2,line=2.1,cex=1.5)
	mtext(mainLab[i-1],side=3,line=1,cex=1.5)
	if (i == 2) {legend(35,1.25,legend = c('LV results','replicate dataset','ALVI-LR method','MAR no transform','MAR log transform','MAR no transform with smoothing','MAR log transform with smoothing'),lty = c(0,0,1,1,1,1,1),pch = c(20,1,NA,NA,NA,NA,NA),col = c(colorPallet[2],colorPallet[1],colorPallet[3],colorPallet[4],colorPallet[5],'orange','khaki'),text.col=c('black','black','black','black','black','black','black'),bty = "n",cex=1,lwd=c(NA,NA,3,3,3,3,3))}
	}
# dev.off()

##### figure S5 - noisy and clustered data - lin reg
##########################################################################################



##########################################################################################
##### par estimation - alg for original data - fig S2

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


# setwd("C:\\Users\\dolivenca3\\OneDrive\\2_2019_America\\2020\\20200123_MAR\\Camparison_LV_MAR\\11_Frontiers\\V3\\Figures")
# tiff("Fig_S2.tiff", height = 40, width = 20, units = 'cm',compression = "lzw", res = 300)
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

##### par estimation - alg for original data - fig S2
############################################################################################################################



######################################################################################################################################
##### Fig S6 - different levels of process noise

set.seed(1)
rnorm(1)
t = seq(1,100,.1)

Process_noise_data = function(sd_ = 0.005)
	{
	d_output1 = c(t=min(t),state1) # initial state and starting responce matrix
	for ( i in seq(min(t)+(t[2]-t[1]),max(t),t[2]-t[1]) )	# cycle to all timepoints
		{
		# retriving the variable values from the responce matrix
		if ( is.null(dim(d_output1))==TRUE )	# If is the first step
			{
			X1 = d_output1[2]
			X2 = d_output1[3]	
			X3 = d_output1[4]
			X4 = d_output1[5]
			} else 					# if not first step
			{
			X1 = d_output1[dim(d_output1)[1],2]
			X2 = d_output1[dim(d_output1)[1],3]
			X3 = d_output1[dim(d_output1)[1],4]
			X4 = d_output1[dim(d_output1)[1],5]
			}
	
		Xs = with(as.list(c(c(X1,X2,X3,X4), Pars)), 
	            {		
			# discrete equations
			X1_temp = ( X1 * (a1 + b11*X1 + b12*X2 + b13*X3 + b14*X4 ) * (t[2]-t[1]) + X1 ) * rnorm(1,1,sd_)
			X2_temp = ( X2 * (a2 + b21*X1 + b22*X2 + b23*X3 + b24*X4 ) * (t[2]-t[1]) + X2 ) * rnorm(1,1,sd_)
			X3_temp = ( X3 * (a3 + b31*X1 + b32*X2 + b33*X3 + b34*X4 ) * (t[2]-t[1]) + X3 ) * rnorm(1,1,sd_)
			X4_temp = ( X4 * (a4 + b41*X1 + b42*X2 + b43*X3 + b44*X4 ) * (t[2]-t[1]) + X4 ) * rnorm(1,1,sd_)		

	            return(list(c(X1_temp,X2_temp,X3_temp,X4_temp)))
	            })
		# updating the responce matrix
		d_output1 = rbind(d_output1,c(t=i,X1=Xs[[1]][1],X2=Xs[[1]][2],X3=Xs[[1]][3],X4=Xs[[1]][4]))
		}
	
	noisyData = d_output1[c(1,sort(sample(seq(11,991,10),39))),1:5]
	rownames(noisyData) = NULL
	return(noisyData)
	}
Proc_data_001 = Process_noise_data(sd=0.001)
Proc_data_005 = Process_noise_data(sd=0.005)
Proc_data_01 = Process_noise_data(sd=0.01)
Proc_data_05 = Process_noise_data(sd=0.05)

# visualization
# setwd("C:\\Users\\dolivenca3\\OneDrive\\2_2019_America\\2020\\20200123_MAR\\Camparison_LV_MAR\\11_Frontiers\\V3\\Figures"
# tiff("Fig_new_S6.tiff", height = 40, width = 10, units = 'cm',compression = "lzw", res = 300)
par(mfrow=c(4,1))
for (i in 2:5)
	{
	plot(out1[,1],out1[,i],pch=20,col='grey',lwd=3,xlab='',ylab='',cex.lab=1.5,cex.axis=1.5)
	lines(Proc_data_001[,1],Proc_data_001[,i],col='blue',lwd=3)
	lines(Proc_data_005[,1],Proc_data_005[,i],col='green',lwd=3)
	lines(Proc_data_01[,1],Proc_data_01[,i],col='orange',lwd=3)
	lines(Proc_data_05[,1],Proc_data_05[,i],col='cyan',lwd=3)
	mtext(expression(italic('t')),side=1,line=2.1,cex=1.5)
	mtext(y_labs[i-1],side=2,line=2.1,cex=1.5)
	}
legend(60,1.3,legend = c('SD','0','0.001','0.005','0.01','0.05'),pch=c(NA,16,NA,NA,NA,NA),lty = c(NA,NA,1,1,1,1),lwd=c(NA,NA,3,3,3,3),col = c('black','grey','blue','green','orange','cyan'),bty = "n",cex=2)	# ,lty = c(0,1)
# dev.off()

##### Fig S6 - different levels of process noise
######################################################################################################################################

