### MAR - multivariate auto regression ###


rm(list=ls())


varColor_light = c('lightblue','orange','lightgreen','grey')
varColor = c('blue','orange3','darkgreen','black')
varNames = c('X1','X2','X3','X4')

### population component control ###
# initial conditions
N1_init = 5
N2_init = 10
N3_init = 15
N4_init = 20

# growth rate when pop = 0
r1 = 1.02106	# runif(1,.5,2)
r2 = 2 		# runif(1,.5,2)
r3 = 0.740671	# runif(1,.5,2)
r4 = 2		# runif(1,.5,2)

# intra and inter spicies interaction component
alpha11 = 0.45			# runif(1,.01,1.9)
alpha12 = 0.5			# runif(1,.01,.5)
alpha13 = -.1			# runif(1,.01,1.9)
alpha14 = -.3			# runif(1,.01,.5)
alpha21 = -0.5		 	# runif(1,.01,.5)
alpha22 = 0.544635		# runif(1,.01,1.9)
alpha23 = .3			# runif(1,.01,.5)
alpha24 = .1			# runif(1,.01,1.9)
alpha31 = -.2			# runif(1,.01,.5)
alpha32 = -0.3			# runif(1,.01,1.9)
alpha33 = 0.9			# -abs(runif(1,.01,.5))
alpha34 = .3			# runif(1,.01,1.9)
alpha41 = .3			# runif(1,.01,.5)
alpha42 = -.1 			# runif(1,.01,1.9)
alpha43 = -.3			# runif(1,.01,.5)
alpha44 = .3			#-abs(runif(1,.01,1.9))

### random component control ###
rand_on = 0
rand_mean = 0
rand_sd = .1			# randomness : .01 small / .1 equivalent to reality / .5 big 

### interesting parameters ###
# r1 = 0.5162979; r2 = 0.5182981; alpha11 = 1.155862; alpha12 = 0.05086701; alpha21 = 0.3608059; alpha22 = 0.1497929
# r1 = 2.3; r2 = 1.6; alpha11 = .9; alpha12 = .1; alpha21 = .5; alpha22 = 0.1981737
# r1 = 1.02106; r2 = 1.551703; alpha11 = 1.615045; alpha12 = 0.3177388; alpha21 = 0.2395321; alpha22 = 1.455365



### MAR / GOMPERTZ in log and cartesian space in matrix form ###

#set.seed(100)

Ni = c(N1_init,N2_init,N3_init,N4_init)

r = c(r1,r2,r3,r4)

Id = diag(length(Ni))

A = t(matrix(c(alpha11,alpha12,alpha13,alpha14,
		alpha21,alpha22,alpha23,alpha24,
		alpha31,alpha32,alpha33,alpha34,
		alpha41,alpha42,alpha43,alpha44),
		nrow=4))

ni_store = log(Ni)
for (i in 1:30)
	{
	n_t = log(Ni)
	n_tplus1 = r + (A) %*% n_t + rand_on * rnorm(n=2, mean = rand_mean, sd = rand_sd)
	ni_store = rbind(ni_store,t(n_tplus1))

	Ni = exp(n_tplus1)	
	}
MAR_artData1 = exp(ni_store)

par(mfrow=c(2,1))
plot(0,0,col='white',xlim=c(0,dim(ni_store)[1]+.1*dim(ni_store)[1]),ylim=c(0,5),main = 'log')
for (ii in 1:length(Ni))
	{
	points(log(MAR_artData1[,ii]),type='l',col=varColor[ii],pch='.')
	legend(dim(MAR_artData1)[1],log(MAR_artData1[dim(MAR_artData1)[1],ii]),legend=varNames[ii],col=varColor[ii],cex=.5,bty = "n")
	}

plot(10,10,col='white',xlim=c(0,dim(ni_store)[1]+.1*dim(ni_store)[1]),ylim=c(-10,80),main = 'cart')
for (ii in 1:length(Ni))
	{
	points(MAR_artData1[,ii],type='l',col=varColor[ii],pch='.')
	legend(dim(MAR_artData1)[1],MAR_artData1[dim(MAR_artData1)[1],ii],legend=varNames[ii],text.col=varColor[ii],cex=.5,bty = "n")
	}


#######################
### MAR_artData1_noisy

set.seed(100)

MAR_artData1_noisy = MAR_artData1
for (i in 1:4)
	{
	MAR_artData1_noisy = MAR_artData1 + rnorm(dim(MAR_artData1)[1],0,.2*mean(MAR_artData1[,i]))
	}

MAR_artData1_noisy[MAR_artData1_noisy < 0] = .005
MAR_artData1_noisy

for (ii in 1:4) { points (MAR_artData1_noisy[,ii],col=varColor[ii]) }

# setwd('C:\\Users\\doliv\\OneDrive\\2_2019_America\\2020\\20200123_MAR\\Camparison_LV_MAR\\8_methods_of_ecology_and_evolution_2sub\\paper_scripts_GitHub')
# write.csv(MAR_artData1_noisy,'MAR_artData1_noisy.csv')
# MAR_artData1_noisy = read.csv('MAR_artData1_noisy.csv')
# save(MAR_artData1_noisy, file = "MAR_artData1_noisy.Rdata")
 load(file = "MAR_artData1_noisy.Rdata")

cbind(apply(MAR_artData1,2,mean)*.2)

### MARartData_noise
#######################


#######################
### MARartData1

par(mfrow=c(2,2))
for (ii in 1:length(Ni))
	{
	plot(MAR_artData1[,ii],type = 'l',lwd=3,col='lightgrey',xlim=c(0,dim(ni_store)[1]+.1*dim(ni_store)[1]),ylim=c(0,80),xlab='t',ylab=varNames[ii])
	points(MAR_artData1_noisy[,ii],pch=20,col=varColor[ii])
	if (ii==1) { legend(15,80,legend=c('MAR results','noisy replicates'),lty=c(1,NA),lwd=c(3,NA),pch=c(NA,20),col=c('lightgrey',varColor[ii]),bty = "n") }
	}

### MARartData1
#######################



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



################################
### time, initState and truePars 

t_MAR1 = seq(1,dim(MAR_artData1_noisy)[1],1)			# time
MAR_artData1_noisy_t = cbind(t_MAR1,MAR_artData1_noisy)	# data
state_fig1 = MAR_artData1_noisy[1,]					# state
names(state_fig1) = c('x1','x2','x3','x4')			# varNames
#truePars_fig1 = c(		
#			a1 = 1.259,
#			b11 = -0.005
#		      )
#truePars_mat_fig1 = Format_pars(truePars = truePars_fig1)	# pars

### time, initState and truePars 
################################



#########################
### slope estimates

smoother_out = Smoother(
				dataset = MAR_artData1_noisy_t,
				draw = TRUE,
				data1_spline2 = 2, 
				smooth_type = 1,
				df = c(15,15,15,15),	
				log_spline = TRUE
				)
plot(1,1,col='white')
legend(1,1,legend = c('data','splines','slopes'),pch = c(20,NA,NA),lty = c(0,1,1),col = c('grey','green','lightgreen'),bty = "n")

AlgPointFind(
		smoother_out_APF = smoother_out,
		dif_equations = Equations,
		matrixEq = TRUE
		)
# 10% "100 % complete || DF = 15 15 15 15  ||  2 3 8 14 27 1190.243  ||  2 3 8 14 27  516.843"
# 20% "100 % complete || DF = 15 15 15 15  ||  1 3 5 15 27 2022.96   ||  2 6 15 18 26 750.132"



( est1 = LV_pars_finder(
				smooth_out = smoother_out,
				alg1_lm2 = 1, 
				data_sample_alg = c(2, 6, 15, 18, 26)		#c(3, 7, 10, 14, 22)		
				) )
cbind(est1)
est1_mat = Format_pars(truePars = est1)
splineStart = smoother_out$splines_est[1,2:5];names(splineStart) = c('x1','x2','x3','x4')
out_est1 = solveLV(times = t_MAR1, pars = est1_mat, initState = splineStart, equations = Equations)

par(mfrow=c(2,2),mar=c(5, 6, 4, 2) + 0.1)
for (iii in 1:4)
	{
	plot(MAR_artData1[,iii],type = 'l',lwd=3,col='lightgrey',xlim=c(0,dim(ni_store)[1]+.1*dim(ni_store)[1]),ylim=c(-10,80),xlab='t',ylab=varNames[iii])
	points(MAR_artData1_noisy[,iii],pch=20,col=varColor[iii])
	if (iii==1) { legend(15,80,legend=c('data','noisy data'),lty=c(1,NA),lwd=c(3,NA),pch=c(NA,20),col=c('lightgrey',varColor[ii]),bty = "n") }
	points(out_est1[,1],out_est1[,iii+1],type='l',col='cyan',lwd=3)
	}

sum( (MAR_artData1-out_est1[,2:5])^2 )

### slope estimates
#########################


#######################
### MAR estimates

MARest_orig1 =  fastMAR(
			data_vars = MAR_artData1_noisy,
			log_transform = FALSE,
			demeaned = FALSE,
			abstol_= 0.01,
			maxit_ = 500,
			includeNoise = FALSE,
			estInit = FALSE,
			tSteps_ = NULL,
			nsim_ = 100
			) 
MARest_orig1$MARSSci_output

MARest_log1 =  fastMAR(
			data_vars = MAR_artData1_noisy,
			log_transform = TRUE,
			demeaned = FALSE,
			abstol_= 0.01,
			maxit_ = 10000,
			includeNoise = FALSE,
			estInit = FALSE,
			tSteps_ = NULL,
			nsim_ = 100
			) 
MARest_log1$MARSSci_output

par(mfrow=c(2,2),mar=c(5, 6, 4, 2) + 0.1)
for (iii in 1:4)
	{
	plot(MAR_artData1[,iii],type = 'l',lwd=3,col='lightgrey',xlim=c(0,dim(ni_store)[1]+.1*dim(ni_store)[1]),ylim=c(-10,80),xlab='t',ylab=varNames[ii])
	points(MAR_artData1_noisy[,iii],pch=20,col=varColor[iii])
	if (iii==1) { legend(15,80,legend=c('data','noisy data'),lty=c(1,NA),lwd=c(3,NA),pch=c(NA,20),col=c('lightgrey',varColor[ii]),bty = "n") }
	points(MARest_orig1$MAR_TS_Est[,iii] ,type='l',col='green',lwd=1)
	points(MARest_log1$MAR_TS_Est[,iii] ,type='l',col='green',lty=2,lwd=1)
	points(meanMaker(MARest_orig1)[iii,],type='l',col='purple',lty=1)
	points(meanMaker(MARest_log1)[iii,],type='l',col='purple',lty=2)
	}

sum( (MAR_artData1-MARest_orig1$MAR_TS_Est)^2 )
sum( (MAR_artData1-MARest_log1$MAR_TS_Est)^2 )

### MAR estimates
#######################



################################################
##### MAR with smooth data 

# smoothing the data for MAR
smoother_for_MAR_out = Smoother_for_MAR(	
						dataset = MAR_artData1_noisy_t,
						df = c(15,15,15,15),
						dataWeights = NULL,
						splineMethod = "fmm",
						draw = TRUE,
						log_spline = TRUE
						)
smoother_for_MAR_out
TS_smooth_data = smooth_TS = smoother_for_MAR_out$splines_est

# MAR estimates without data log transformation, with spline smooth data  
MARest_orig1_smoothData =  fastMAR(
			data_vars = TS_smooth_data[,2:5],
			log_transform=FALSE,
			demeaned=FALSE,
			abstol_=0.01,
			maxit_ = 1000,
			includeNoise = FALSE,
			estInit = FALSE,
			tSteps_ = NULL,
			nsim_ = 100
			) 
MARest_orig1_smoothData$MARSSci_output

# MAR estimates with data log transformation, with spline smooth data  
MARest_log1_smoothData =  fastMAR(
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
MARest_log1_smoothData$MARSSci_output

sum( (MAR_artData1-MARest_orig1_smoothData$MAR_TS_Est)^2 )
sum( (MAR_artData1-MARest_log1_smoothData$MAR_TS_Est)^2 )

##### MAR with smooth data 
################################################



#####################
##### fig2a 1st try

par(mfrow=c(2,2),mar=c(5, 6, 4, 2) + 0.1)
for (iii in 1:4)
	{
	plot(MAR_artData1[,iii],type = 'l',lwd=3,col='lightgrey',xlim=c(0,dim(ni_store)[1]+.1*dim(ni_store)[1]),ylim=c(-10,80),xlab='t',ylab=varNames[iii])
	points(MAR_artData1_noisy[,iii],pch=20,col=varColor[iii])
	if (iii==1) { legend(15,80,legend=c('data','noisy data','slope est','MAR','MAR logTrans'),lty=c(1,NA,1,1,2),lwd=c(3,NA,3,1,1),pch=c(NA,20,NA,NA,NA),col=c('lightgrey',varColor[ii],'cyan','green','green'),bty = "n") }
	points(MARest_orig1$MAR_TS_Est[,iii] ,type='l',col='green',lwd=1)
	points(MARest_log1$MAR_TS_Est[,iii] ,type='l',col='green',lty=2,lwd=1)
	points(MARest_orig1_smoothData$MAR_TS_Est[,iii] ,type='l',col='orange',lwd=1)
	points(MARest_log1_smoothData$MAR_TS_Est[,iii] ,type='l',col='orange',lty=2,lwd=1)
	points(out_est1[,1],out_est1[,iii+1],type='l',col='cyan',lwd=3)
	}

##### fig2a 1st try
#####################








##########################################################################
##### more dynamics

### pseudo sustained oxcilations
r1 = 1.02106	# runif(1,.5,2)
r2 = 2 		# runif(1,.5,2)
r3 = 0.740671	# runif(1,.5,2)
r4 = 2		# runif(1,.5,2)

# intra and inter spicies interaction component
alpha11 = 0.45			# runif(1,.01,1.9)
alpha12 = 0.7			# runif(1,.01,.5)
alpha13 = -.1			# runif(1,.01,1.9)
alpha14 = -.3			# runif(1,.01,.5)
alpha21 = -0.7		 	# runif(1,.01,.5)
alpha22 = 0.544635		# runif(1,.01,1.9)
alpha23 = .3			# runif(1,.01,.5)
alpha24 = .1			# runif(1,.01,1.9)
alpha31 = -.2			# runif(1,.01,.5)
alpha32 = -.3			# runif(1,.01,1.9)
alpha33 = 0.9			# -abs(runif(1,.01,.5))
alpha34 = .3			# runif(1,.01,1.9)
alpha41 = .3			# runif(1,.01,.5)
alpha42 = -.1 			# runif(1,.01,1.9)
alpha43 = -.3			# runif(1,.01,.5)
alpha44 = .3			#-abs(runif(1,.01,1.9))



### sustained oxilations
# growth rate when pop = 0
r1 = 1.02106	# runif(1,.5,2)
r2 = 2 		# runif(1,.5,2)
r3 = 2.740671	# runif(1,.5,2)
r4 = 2.74067		# runif(1,.5,2)

# intra and inter spicies interaction component
alpha11 = -0.45			# runif(1,.01,1.9)
alpha12 = .5			# runif(1,.01,.5)
alpha13 = 0			# runif(1,.01,1.9)
alpha14 = 0			# runif(1,.01,.5)
alpha21 = .5		 	# runif(1,.01,.5)
alpha22 = -0.544635		# runif(1,.01,1.9)
alpha23 = 0			# runif(1,.01,.5)
alpha24 = 0			# runif(1,.01,1.9)
alpha31 = -.5			# runif(1,.01,.5)
alpha32 = 0			# runif(1,.01,1.9)
alpha33 = -0.9			# -abs(runif(1,.01,.5))
alpha34 = .1			# runif(1,.01,1.9)
alpha41 = -.5			# runif(1,.01,.5)
alpha42 = 0			# runif(1,.01,1.9)
alpha43 = .1			# runif(1,.01,.5)
alpha44 = -.9			#-abs(runif(1,.01,1.9))



### slow dampning oxcilations
# growth rate when pop = 0
r1 = 1.02106	# runif(1,.5,2)
r2 = 2 		# runif(1,.5,2)
r3 = 2.740671	# runif(1,.5,2)
r4 = 2.74067		# runif(1,.5,2)

# intra and inter spicies interaction component
alpha11 = -0.45			# runif(1,.01,1.9)
alpha12 = .5			# runif(1,.01,.5)
alpha13 = -.3			# runif(1,.01,1.9)
alpha14 = .5			# runif(1,.01,.5)
alpha21 = .3		 	# runif(1,.01,.5)
alpha22 = -0.544635		# runif(1,.01,1.9)
alpha23 = -.1			# runif(1,.01,.5)
alpha24 = -.1			# runif(1,.01,1.9)
alpha31 = -.22			# runif(1,.01,.5)
alpha32 = -.1			# runif(1,.01,1.9)
alpha33 = -0.9			# -abs(runif(1,.01,.5))
alpha34 = .1			# runif(1,.01,1.9)
alpha41 = -.5			# runif(1,.01,.5)
alpha42 = .1			# runif(1,.01,1.9)
alpha43 = .1			# runif(1,.01,.5)
alpha44 = -.1			#-abs(runif(1,.01,1.9))



### Another dampned oxcilations
# growth rate when pop = 0
r1 = 5	# runif(1,.5,2)
r2 = 5 		# runif(1,.5,2)
r3 = 2.740671	# runif(1,.5,2)
r4 = 0.5		# runif(1,.5,2)

# intra and inter spicies interaction component
alpha11 = -0.045			# runif(1,.01,1.9)
alpha12 = -0.05			# runif(1,.01,.5)
alpha13 = 0			# runif(1,.01,1.9)
alpha14 = -.7			# runif(1,.01,.5)
alpha21 = 0.05		 	# runif(1,.01,.5)
alpha22 = -0.01		# runif(1,.01,1.9)
alpha23 = -.5			# runif(1,.01,.5)
alpha24 = 0			# runif(1,.01,1.9)
alpha31 = 0			# runif(1,.01,.5)
alpha32 = 1			# runif(1,.01,1.9)
alpha33 = -.02			# -abs(runif(1,.01,.5))
alpha34 = -.5			# runif(1,.01,1.9)
alpha41 = .7			# runif(1,.01,.5)
alpha42 = 0			# runif(1,.01,1.9)
alpha43 = .5			# runif(1,.01,.5)
alpha44 = -.02			#-abs(runif(1,.01,1.9))


# another pseudo stable oscilation
# growth rate when pop = 0
r1 = 1.02106	# runif(1,.5,2)
r2 = 2 		# runif(1,.5,2)
r3 = 2.740671	# runif(1,.5,2)
r4 = 2.74067		# runif(1,.5,2)

# intra and inter spicies interaction component
alpha11 = -0.45			# runif(1,.01,1.9)
alpha12 = .5			# runif(1,.01,.5)
alpha13 = -.3			# runif(1,.01,1.9)
alpha14 = .5			# runif(1,.01,.5)
alpha21 = .3		 	# runif(1,.01,.5)
alpha22 = -0.544635		# runif(1,.01,1.9)
alpha23 = -.1			# runif(1,.01,.5)
alpha24 = -.1			# runif(1,.01,1.9)
alpha31 = -.22			# runif(1,.01,.5)
alpha32 = -.1			# runif(1,.01,1.9)
alpha33 = -0.9			# -abs(runif(1,.01,.5))
alpha34 = .1			# runif(1,.01,1.9)
alpha41 = -.5			# runif(1,.01,.5)
alpha42 = .1			# runif(1,.01,1.9)
alpha43 = .1			# runif(1,.01,.5)
alpha44 = -.1			#-abs(runif(1,.01,1.9))








#############################################################################
##### clustered dataset (noisy data 4) - 11 points were chosen from the artificial data and for each point five observations were created multiplying it by a random normal of mean 1 and variance 0.5

MAR_artData1_withTime = cbind(1:dim(MAR_artData1)[1],MAR_artData1)

set.seed(100)
( MAR_artData1_forClust = MAR_artData1_withTime[c(1,3,5,7,9,11,13,15,18,22,26,30),] )  #  #seq(0,101,10) #c(1,6,11,18,26,41,seq(51,101,10))

datapoint_n = 5
MAR_artData1_clustered = matrix(rep(NA,5),nrow=1)
for (i in 1:dim(MAR_artData1_forClust)[1])
	{
      for (ii in 1:datapoint_n)
		{
		MAR_artData1_clustered = rbind( MAR_artData1_clustered, c(MAR_artData1_forClust[i,1], MAR_artData1_forClust[i,2:5] * rnorm(n = 4,mean = 1,sd = .2)) )
		}
	}
MAR_artData1_clustered = MAR_artData1_clustered[2:dim(MAR_artData1_clustered)[1],]

# setwd('C:\\Users\\doliv\\OneDrive\\2_2019_America\\2020\\20200123_MAR\\Camparison_LV_MAR\\8_methods_of_ecology_and_evolution_2sub\\paper_scripts_GitHub')
# write.csv(MAR_artData1_clustered,'MAR_artData1_clustered.csv')
# MAR_artData1_clustered = read.csv('MAR_artData1_clustered.csv')
# save(MAR_artData1_clustered, file = "MAR_artData1_clustered.Rdata")
 load(file = "MAR_artData1_clustered.Rdata")



# tiff("Fig_1b.tiff", height = 20, width = 20, units = 'cm',compression = "lzw", res = 300)
par(mfrow=c(2,2))
for (ii in 1:length(Ni))
	{
	plot(MAR_artData1[,ii],type = 'l',lwd=3,col='lightgrey',xlim=c(0,dim(ni_store)[1]+.1*dim(ni_store)[1]),ylim=c(0,80),xlab='t',ylab=varNames[ii])
	points(MAR_artData1_clustered[,1],MAR_artData1_clustered[,ii+1],pch=20,col=varColor[ii])
	if (ii==1) { legend(15,80,legend=c('MAR results','replicates data'),lty=c(1,NA),lwd=c(3,NA),pch=c(NA,20),col=c('lightgrey',varColor[ii]),bty = "n") }
	}
# dev.off()

##### clustered dataset (noisy data 4) - 11 points were chosen from the artificial data and for each point five observations were created multiplying it by a random normal of mean 1 and variance 0.5
#############################################################################



#########################
### slope estimates

MAR_artData1_clustered_means = MAR_artData1_forClust[1:11,]
for (i in 1:(dim(MAR_artData1_forClust)[1]-1))
	{
	MAR_artData1_clustered_means[i,] = apply(MAR_artData1_clustered[(i*5-4):(i*5),],2,mean)
	}

smoother_out = Smoother(
				dataset = MAR_artData1_clustered_means,
				draw = TRUE,
				data1_spline2 = 2, 
				smooth_type = 1,
				df = c(11,11,11,11),
				log_spline = TRUE
				)
plot(1,1,col='white')
legend(1,1,legend = c('data','splines','slopes'),pch = c(20,NA,NA),lty = c(0,1,1),col = c('grey','green','lightgreen'),bty = "n")

# remove the # comment to do a point subsample 
AlgPointFind(
		smoother_out_APF = smoother_out,
		dif_equations = Equations,
		matrixEq = TRUE
		)
#  "100 % complete || DF = 11 11 11 11  ||  1 2 6 7 8 581.513  ||  1 2 6 7 8 581.513"

( est1 = LV_pars_finder(
				smooth_out = smoother_out,
				alg1_lm2 = 1, 
				data_sample_alg = c(1,2,6,7,8) # c(5,6,7,8,9)		#c(3, 7, 10, 14, 22)		
				) )
cbind(est1)
est1_mat = Format_pars(truePars = est1)
splineStart = smoother_out$splines_est[1,2:5];names(splineStart) = c('x1','x2','x3','x4')
out_est2 = solveLV(times = t_MAR1, pars = est1_mat, initState = splineStart, equations = Equations)

par(mfrow=c(2,2),mar=c(5, 6, 4, 2) + 0.1)
for (iii in 1:4)
	{
	plot(MAR_artData1[,iii],type = 'l',lwd=3,col='lightgrey',xlim=c(0,dim(ni_store)[1]+.1*dim(ni_store)[1]),ylim=c(-10,80),xlab='t',ylab=varNames[iii])
	points(MAR_artData1_clustered[,1],MAR_artData1_clustered[,iii+1],pch=20,col=varColor[iii])
	if (iii==1) { legend(15,80,legend=c('data','noisy data'),lty=c(1,NA),lwd=c(3,NA),pch=c(NA,20),col=c('lightgrey',varColor[ii]),bty = "n") }
	points(out_est2[,1],out_est2[,iii+1],type='l',col='cyan',lwd=3)
	}

# SSE
sum( ( MAR_artData1 - out_est2[,2:5] )^2 )

### slope estimates
#########################


#######################
### MAR estimates

MAR_artData1_clustered_means

TS_MAR_artData1_clustered_means = cbind(1:31,matrix(rep(NA,31*4),ncol=4))
for (i in 1:100)
	{
	if (i%in%MAR_artData1_clustered_means[,1])
		{
		TS_MAR_artData1_clustered_means[i,] = MAR_artData1_clustered_means[MAR_artData1_clustered_means[,1]==i,] 
		} 
	}

# setwd('C:\\Users\\doliv\\OneDrive\\2_2019_America\\2020\\20200123_MAR\\Camparison_LV_MAR\\3_paper_new_functions\\paper_scripts')
# write.csv(TS_MAR_artData1_clustered_means,'TS_MAR_artData1_clustered_means.csv')
# TS_MAR_artData1_clustered_means = read.csv('TS_MAR_artData1_clustered_means.csv')
# save(TS_MAR_artData1_clustered_means, file = "TS_MAR_artData1_clustered_means.Rdata")
# load(file = "TS_MAR_artData1_clustered_means")


MARest_orig2 =  fastMAR(
			data_vars = TS_MAR_artData1_clustered_means[,2:5],
			log_transform = FALSE,
			demeaned = TRUE,
			abstol_= 0.01,
			maxit_ = 1000,
			includeNoise = FALSE,
			estInit = FALSE,
			tSteps_ = NULL,
			nsim_ = 100
			) 
MARest_orig2$MARSSci_output

# SSE
sum( ( MAR_artData1 - MARest_orig2$MAR_TS_Est )^2 )


MARest_log2 =  fastMAR(
			data_vars = TS_MAR_artData1_clustered_means[,2:5],
			log_transform = TRUE,
			demeaned = TRUE,
			abstol_= 0.01,
			maxit_ = 1000,
			includeNoise = FALSE,
			estInit = FALSE,
			tSteps_ = NULL,
			nsim_ = 100
			) 
MARest_log2$MARSSci_output

# SSE
sum( ( MAR_artData1 - MARest_log2$MAR_TS_Est )^2 )


par(mfrow=c(2,2),mar=c(5, 6, 4, 2) + 0.1)
for (iii in 1:4)
	{
	plot(MAR_artData1[,iii],type = 'l',lwd=3,col='lightgrey',xlim=c(0,dim(ni_store)[1]+.1*dim(ni_store)[1]),ylim=c(-10,80),xlab='t',ylab=varNames[ii])
	points(MAR_artData1_clustered[,1],MAR_artData1_clustered[,iii+1],pch=20,col=varColor[iii])
	if (iii==1) { legend(15,80,legend=c('data','noisy data'),lty=c(1,NA),lwd=c(3,NA),pch=c(NA,20),col=c('lightgrey',varColor[ii]),bty = "n") }
	points(MARest_orig2$MAR_TS_Est[,iii] ,type='l',col='green',lwd=1)
	points(MARest_log2$MAR_TS_Est[,iii] ,type='l',col='green',lty=2,lwd=1)
	points(meanMaker(MARest_orig2)[iii,],type='l',col='purple',lty=1)
	points(meanMaker(MARest_log2)[iii,],type='l',col='purple',lty=2)
	}

### MAR estimates
#######################



################################################
##### MAR with smooth data 

# smoothing the data for MAR
smoother_for_MAR_out = Smoother_for_MAR(	
						dataset = MAR_artData1_clustered,
						df = c(12,12,12,12),
						dataWeights = NULL,
						splineMethod = "fmm",
						draw = TRUE,
						log_spline = TRUE
						)
smoother_for_MAR_out
smooth_TS = smoother_for_MAR_out$splines_est

# Completing the smooth data to have the same time structure as the original data
TS_smooth_data = cbind(1:31,matrix(rep(NA,31*4),ncol=4))
for (i in 1:31)
	{
	if (i%in%smooth_TS[,1])
		{
		TS_smooth_data[i,] = smooth_TS[smooth_TS[,1]==i,] 
		} 
	}

# MAR estimates without data log transformation, with spline smooth data  
MARest_orig2_smoothData =  fastMAR(
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
MARest_orig2_smoothData$MARSSci_output

# MAR estimates with data log transformation, with spline smooth data  
MARest_log2_smoothData =  fastMAR(
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
MARest_log2_smoothData$MARSSci_output

sum( (MAR_artData1-MARest_orig2_smoothData$MAR_TS_Est)^2 )
sum( (MAR_artData1-MARest_log2_smoothData$MAR_TS_Est)^2 )

##### MAR with smooth data 
################################################



######################################################################################################################################
##### colors

colorPallet = c('black','grey','blue','darkgreen','green')

##### colors
######################################################################################################################################



#############################################################################
##### figure 2 - noisy and clustered MAR data 20% - alg1

# setwd("C:\\Users\\doliv\\OneDrive\\2_2019_America\\2020\\20200123_MAR\\Camparison_LV_MAR\\3_paper_new_functions\\paper_scripts")     # P51
# MAR_artData1_noisy = read.csv('MAR_artData1_noisy.csv')
# load(file = "MAR_artData1_noisy.Rdata")
# MAR_artData1_clustered = read.csv('MAR_artData1_clustered.csv')
# load(file = "MAR_artData1_clustered.Rdata")

# setwd('C:\\Users\\doliv\\OneDrive\\2_2019_America\\2020\\20200123_MAR\\Camparison_LV_MAR\\8_methods_of_ecology_and_evolution_2sub\\paper_figures')
# tiff("Fig_2.tiff", height = 40, width = 20, units = 'cm',compression = "lzw", res = 300)
# Function to create a four independent color plot
y_labs = c(expression(italic('X'[1])),expression(italic('X'[2])),expression(italic('X'[3])),expression(italic('X'[4])))
par(mfcol=c(4,2))
mainLab = c('Noisy MAR Dataset',NA,NA,NA)
for (i in 1:4)
	{
	plot(MAR_artData1[,i],type = 'l',lwd=3,col=colorPallet[2],xlim=c(0,dim(ni_store)[1]+.1*dim(ni_store)[1]),ylim=c(-10,80),main='',xlab='',ylab='',cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
	points(MAR_artData1_noisy[,i],pch=1,col=colorPallet[1])
	if (i==1) { legend(10,80,legend=c('MAR results','noisy MAR','ALVI-MI method','MAR no transform','MAR log transform','MAR no transform with smoothness','MAR log transform with smoothness'),lty=c(1,NA,1,1,1,1,1),lwd=c(3,NA,3,3,3,3,3),pch=c(NA,1,NA,NA,NA,NA,NA),col=c(colorPallet[2],colorPallet[1],colorPallet[3],colorPallet[4],colorPallet[5],'orange','khaki'),bty = "n") }
	points(MARest_orig1$MAR_TS_Est[,i] ,type='l',col=colorPallet[4],lwd=3)
	points(MARest_log1$MAR_TS_Est[,i] ,type='l',col=colorPallet[5],lwd=3)
	points(MARest_orig1_smoothData$MAR_TS_Est[,i] ,type='l',col='orange',lwd=3)
	points(MARest_log1_smoothData$MAR_TS_Est[,i] ,type='l',col='khaki',lwd=3)
	points(out_est1[,1],out_est1[,i+1],type='l',col=colorPallet[3],lwd=3)
	mtext(expression(italic('t')),side=1,line=2.1,cex=1.5)
	mtext(y_labs[i],side=2,line=2.1,cex=1.5)
	mtext(mainLab[i],side=3,line=1,cex=1.5)
	}
mainLab = c('Replicate MAR Dataset',NA,NA,NA)
for (i in 1:4)
	{
	plot(MAR_artData1[,i],type = 'l',lwd=3,col=colorPallet[2],xlim=c(0,dim(ni_store)[1]+.1*dim(ni_store)[1]),ylim=c(-10,80),main='',xlab='',ylab='',cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
	points(MAR_artData1_clustered[,1],MAR_artData1_clustered[,i+1],pch=1,col=colorPallet[1])
	if (i==1) { legend(10,80,legend=c('MAR results','replicate MAR','ALVI-MI method','MAR no transform','MAR log transform','MAR no transform with smoothness','MAR log transform with smoothness'),lty=c(1,NA,1,1,1,1,1),lwd=c(3,NA,3,3,3,3,3),pch=c(NA,1,NA,NA,NA,NA,NA),col=c(colorPallet[2],colorPallet[1],colorPallet[3],colorPallet[4],colorPallet[5],'orange','khaki'),bty = "n") }
	points(MARest_orig2$MAR_TS_Est[,i] ,type='l',col=colorPallet[4],lwd=3)
	points(MARest_log2$MAR_TS_Est[,i] ,type='l',col=colorPallet[5],lwd=3)
	points(MARest_orig2_smoothData$MAR_TS_Est[,i] ,type='l',col='orange',lwd=3)
	points(MARest_log2_smoothData$MAR_TS_Est[,i] ,type='l',col='khaki',lwd=3)
	points(out_est2[,1],out_est2[,i+1],type='l',col=colorPallet[3],lwd=3)
	mtext(expression(italic('t')),side=1,line=2.1,cex=1.5)
	mtext(y_labs[i],side=2,line=2.1,cex=1.5)
	mtext(mainLab[i],side=3,line=1,cex=1.5)
	}
# dev.off()

##### figure 2 - noisy and clustered MAR data 20% - alg1
#############################################################################




#############################################################################
### MAR noisy slope estimates ALVI - LR

smoother_out = Smoother(
				dataset = MAR_artData1_noisy_t,
				draw = TRUE,
				data1_spline2 = 2, 
				smooth_type = 1,
				df = c(14,15,11,2),	
				log_spline = TRUE
				)
plot(1,1,col='white')
legend(1,1,legend = c('data','splines','slopes'),pch = c(20,NA,NA),lty = c(0,1,1),col = c('grey','green','lightgreen'),bty = "n")


( est1 = LV_pars_finder(
				smooth_out = smoother_out,
				alg1_lm2 = 2
				) )
cbind(est1)
est1_mat = Format_pars(truePars = est1)
splineStart = smoother_out$splines_est[1,2:5];names(splineStart) = c('x1','x2','x3','x4')
out_est1 = solveLV(times = t_MAR1, pars = est1_mat, initState = splineStart, equations = Equations)

par(mfrow=c(2,2),mar=c(5, 6, 4, 2) + 0.1)
for (iii in 1:4)
	{
	plot(MAR_artData1[,iii],type = 'l',lwd=3,col='lightgrey',xlim=c(0,dim(ni_store)[1]+.1*dim(ni_store)[1]),ylim=c(-10,80),xlab='t',ylab=varNames[iii])
	points(MAR_artData1_noisy[,iii],pch=20,col=varColor[iii])
	if (iii==1) { legend(15,80,legend=c('data','noisy data'),lty=c(1,NA),lwd=c(3,NA),pch=c(NA,20),col=c('lightgrey',varColor[ii]),bty = "n") }
	points(out_est1[,1],out_est1[,iii+1],type='l',col='cyan',lwd=3)
	}

sum( (MAR_artData1-out_est1[,2:5])^2 )

### MAR noisy slope estimates ALVI - LR
#############################################################################


#############################################################################
### MAR clustered slope estimates ALVI - LR

MAR_artData1_clustered_means = MAR_artData1_forClust[1:11,]
for (i in 1:(dim(MAR_artData1_forClust)[1]-1))
	{
	MAR_artData1_clustered_means[i,] = apply(MAR_artData1_clustered[(i*5-4):(i*5),],2,mean)
	}

smoother_out = Smoother(
				dataset = MAR_artData1_clustered_means,
				draw = TRUE,
				data1_spline2 = 2, 
				smooth_type = 1,
				df = c(11,11,11,5),
				log_spline = TRUE
				)
plot(1,1,col='white')
legend(1,1,legend = c('data','splines','slopes'),pch = c(20,NA,NA),lty = c(0,1,1),col = c('grey','green','lightgreen'),bty = "n")

( est1 = LV_pars_finder(
				smooth_out = smoother_out,
				alg1_lm2 = 2 
				) )
cbind(est1)
est1_mat = Format_pars(truePars = est1)
splineStart = smoother_out$splines_est[1,2:5];names(splineStart) = c('x1','x2','x3','x4')
out_est2 = solveLV(times = t_MAR1, pars = est1_mat, initState = splineStart, equations = Equations)

par(mfrow=c(2,2),mar=c(5, 6, 4, 2) + 0.1)
for (iii in 1:4)
	{
	plot(MAR_artData1[,iii],type = 'l',lwd=3,col='lightgrey',xlim=c(0,dim(ni_store)[1]+.1*dim(ni_store)[1]),ylim=c(-10,80),xlab='t',ylab=varNames[iii])
	points(MAR_artData1_clustered[,1],MAR_artData1_clustered[,iii+1],pch=20,col=varColor[iii])
	if (iii==1) { legend(15,80,legend=c('data','noisy data'),lty=c(1,NA),lwd=c(3,NA),pch=c(NA,20),col=c('lightgrey',varColor[ii]),bty = "n") }
	points(out_est2[,1],out_est2[,iii+1],type='l',col='cyan',lwd=3)
	}

# SSE
sum( ( MAR_artData1 - out_est2[,2:5] )^2 )

### MAR clustered slope estimates ALVI - LR
#############################################################################



