# MAR real data


rm(list=ls())


library(MARSS)

#####################
##### greywhales

#######
### MAR
data(graywhales) 
# This data is missing values. No 1986, 1988, 1989, 1990, 1991, 1994, 1996. 
# Fixing that.
graywhales_corrected = rbind(graywhales[1:34,],
					c(1986,NA),
					graywhales[35,],
					cbind(c(1988:1991),rep(NA,4)),
					graywhales[36:37,],
					c(1994,NA),
					graywhales[38,],
					c(1996,NA),
					graywhales[39,]
					)
years = graywhales_corrected[,1]
loggraywhales = log(graywhales_corrected[,2])
kem = MARSS(graywhales_corrected[,2])
kem_log = MARSS(loggraywhales)

par(mfrow=c(2,2))
# kem_log log
plot(graywhales_corrected[,1],log(graywhales_corrected[,2]),col='red',xlab='Years',ylab='log( Gray Whales )')
points(graywhales_corrected[,1],kem_log$states, type = "l", lwd = 2)
lines(graywhales_corrected[,1],kem_log$states - 2*kem_log$states.se,lty=3)
lines(graywhales_corrected[,1],kem_log$states + 2*kem_log$states.se,lty=3)

plot(1,1,col='white')

# kem_log cart
plot(graywhales_corrected[,1],graywhales_corrected[,2],col='red',xlab='Years',ylab='Gray Whales')
points(graywhales_corrected[,1],exp(kem_log$states), type = "l", lwd = 2)
lines(graywhales_corrected[,1],exp(kem_log$states - 2*kem_log$states.se),lty=3)
lines(graywhales_corrected[,1],exp(kem_log$states + 2*kem_log$states.se),lty=3)

# kem
plot(graywhales_corrected[,1],graywhales_corrected[,2],col='red',xlab='Years',ylab='Gray Whales')
points(graywhales_corrected[,1],kem$states, type = "l", lwd = 2)
lines(graywhales_corrected[,1],kem$states - 2*kem$states.se,lty=3)
lines(graywhales_corrected[,1],kem$states + 2*kem$states.se,lty=3)
### MAR
#######

##### greywhales
#####################



#####################
### Ives Lake Data

#######
### MAR 

data(ivesDataByWeek)
ivesDataByWeek

plank.dat = t(ivesDataByWeek[,c("Large Phyto","Small Phyto","Daphnia","Non-daphnia")])	# getting the data and log it
d.plank.dat = (plank.dat - apply(plank.dat, 1, mean, na.rm = TRUE))					# centering the data at zero by subtracting the mean
varNames = c("Large Phyto","Small Phyto","Daphnia","Non-daphnia")

Z = factor(rownames(d.plank.dat))
U = "zero"
A = "zero"
B = "unconstrained"
Q = matrix(list(0), 4, 4)
diag(Q) = c("Phyto", "Phyto", "Zoo", "Zoo")
R = diag(c(0.04, 0.04, 0.16, 0.16))
plank.model = list(Z = Z, U = U, Q = Q, R = R, B = B, A = A, tinitx=1)

kem.plank = MARSS(d.plank.dat,model = plank.model)
B.est = matrix(kem.plank$par$B, 4, 4)
rownames(B.est) = colnames(B.est) = c("LP", "SP", "D", "ND")
print(B.est, digits = 2)


# cart
par(mfrow = c(2,2))
for (i in 1:dim(plank.dat)[1])
	{
	plot(ivesDataByWeek[,2],plank.dat[i,],col='red',xlab='Years',ylab=varNames[i] )
	points(ivesDataByWeek[,2],kem.plank$states[i,]+mean(plank.dat[i,],na.rm=TRUE), type = "l", lwd = 2)
	lines(ivesDataByWeek[,2],kem.plank$states[i,]+mean(plank.dat[i,],na.rm=TRUE) - 2*kem.plank$states.se[i,],lty=3)
	lines(ivesDataByWeek[,2],kem.plank$states[i,]+mean(plank.dat[i,],na.rm=TRUE) + 2*kem.plank$states.se[i,],lty=3)
	}

### MAR 
#######


#######
### MAR - logTrans

data(ivesDataByWeek)
ivesDataByWeek

plank.dat_log = t(log(ivesDataByWeek[,c("Large Phyto","Small Phyto","Daphnia","Non-daphnia")]))	# getting the data and log it
d.plank.dat_log = (plank.dat_log - apply(plank.dat_log, 1, mean, na.rm = TRUE))					# centering the data at zero by subtracting the mean
varNames = c("Large Phyto","Small Phyto","Daphnia","Non-daphnia")

Z = factor(rownames(d.plank.dat_log))
U = "zero"
A = "zero"
B = "unconstrained"
Q = matrix(list(0), 4, 4)
diag(Q) = c("Phyto", "Phyto", "Zoo", "Zoo")
R = diag(c(0.04, 0.04, 0.16, 0.16))
plank.model = list(Z = Z, U = U, Q = Q, R = R, B = B, A = A, tinitx=1)

kem.plank_log = MARSS(d.plank.dat_log,model = plank.model)
B.est = matrix(kem.plank_log$par$B, 4, 4)
rownames(B.est) = colnames(B.est) = c("LP", "SP", "D", "ND")
print(B.est, digits = 2)

# log zero centered
par(mfrow = c(2,2))
for (i in 1:dim(plank.dat_log)[1])
	{
	plot(ivesDataByWeek[,2],d.plank.dat_log[i,],col='red',xlab='Years',ylab=paste('log(',varNames[i],')') )
	points(ivesDataByWeek[,2],kem.plank_log$states[i,], type = "l", lwd = 2)
	lines(ivesDataByWeek[,2],kem.plank_log$states[i,] - 2*kem.plank_log$states.se[i,],lty=3)
	lines(ivesDataByWeek[,2],kem.plank_log$states[i,] + 2*kem.plank_log$states.se[i,],lty=3)
	}

# cart zero centered
par(mfrow = c(2,2))
for (i in 1:dim(plank.dat_log)[1])
	{
	plot(ivesDataByWeek[,2],exp(d.plank.dat_log[i,]),col='red',xlab='Years',ylab=varNames[i])
	points(ivesDataByWeek[,2],exp(kem.plank_log$states[i,]), type = "l", lwd = 2)
	lines(ivesDataByWeek[,2],exp(kem.plank_log$states[i,] - 2*kem.plank_log$states.se[i,]),lty=3)
	lines(ivesDataByWeek[,2],exp(kem.plank_log$states[i,] + 2*kem.plank_log$states.se[i,]),lty=3)
	}

# log 
par(mfrow = c(2,2))
for (i in 1:dim(plank.dat_log)[1])
	{
	plot(ivesDataByWeek[,2],plank.dat_log[i,],col='red',xlab='Years',ylab=paste('log(',varNames[i],')') )
	points(ivesDataByWeek[,2],kem.plank_log$states[i,]+mean(plank.dat_log[i,],na.rm=TRUE), type = "l", lwd = 2)
	}

# cart 
par(mfrow = c(2,2))
for (i in 1:dim(plank.dat_log)[1])
	{
	plot(ivesDataByWeek[,2],exp(plank.dat_log[i,]),col='red',xlab='Years',ylab=paste('log(',varNames[i],')') )
	points(ivesDataByWeek[,2],exp(kem.plank_log$states[i,]+mean(plank.dat_log[i,],na.rm=TRUE)), type = "l", lwd = 2)
	
	lines(ivesDataByWeek[,2],exp(kem.plank_log$states[i,]+mean(plank.dat_log[i,],na.rm=TRUE) - 2*kem.plank_log$states.se[i,]),lty=3)
	lines(ivesDataByWeek[,2],exp(kem.plank_log$states[i,]+mean(plank.dat_log[i,],na.rm=TRUE) + 2*kem.plank_log$states.se[i,]),lty=3)
	}

### MAR - logTrans
#######

### Ives Lake Data
#####################



#############################################################################



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



######################################################################################################################################
##### colors

colorPallet = c('black','grey','blue','darkgreen','green')

##### colors
######################################################################################################################################



######################################################################################################################################
##### graywhales

################################
### time, initState and truePars 

data(graywhales)
graywhales
# This data is missing values. No 1986, 1988, 1989, 1990, 1991, 1994, 1996. 
# Fixing that.
graywhales_corrected = rbind(graywhales[1:34,],
					c(1986,NA),
					graywhales[35,],
					cbind(c(1988:1991),rep(NA,4)),
					graywhales[36:37,],
					c(1994,NA),
					graywhales[38,],
					c(1996,NA),
					graywhales[39,]
					)

# Insert the number of time units that the simulation will run
time_graywhales_corrected = graywhales_corrected[,1]
state_graywhales_corrected = graywhales_corrected[1,2]
truePars_graywhales_corrected = c(		
						a1 = 1.259,
						b11 = -0.005
		      			)
truePars_mat_graywhales_corrected = Format_pars(truePars = truePars_graywhales_corrected)

### time, initState and truePars
################################

out1_greywhales = solveLV(times = time_graywhales_corrected, 
				initState = state_graywhales_corrected, 
				pars = truePars_mat_graywhales_corrected, 
				equations = Equations)
# head(out1) #tail(out1,10) # edit(out1) # View(out1)

par(mfrow=c(1,1))
plot(graywhales_corrected[,1],graywhales_corrected[,2],col='red',xlab='Years',ylab='Gray Whales')
points(out1_greywhales[,1],out1_greywhales[,2], type = "l", lwd = 2)
legend(1985,7000,legend = c("data","Slope est"),pch = c(1,NA),lty = c(NA,1),col = c("red","black"),bty = "n") # 

#########################
### slope estimates

graywhales_rm.na = graywhales_corrected[is.na(graywhales_corrected[,2])==FALSE,]

smoother_out = Smoother(
				dataset = graywhales_rm.na,
				draw = TRUE,
				data1_spline2 = 2, 
				smooth_type = 1,
				df = 3,
				log_spline = TRUE
				)
plot(1,1,col='white')
legend(1,1,legend = c('data','splines','slopes'),pch = c(20,NA,NA),lty = c(0,1,1),col = c('grey','green','lightgreen'),bty = "n")

AlgPointFind(
		smoother_out_APF = smoother_out,
		dif_equations = Equations,
		matrixEq = TRUE
		)
# 	"100 % complete || DF = 15  ||  2 14 170511189.416  ||  2 14 114963636.776"
# log "100 % complete || DF = 15  ||  1  2 168599305.82   ||  1 18 114535854.508"
# log "100 % complete || DF =  3  ||  4  5 153703188.988  ||  2 20   2090624.691"

( est1 = LV_pars_finder(
				smooth_out = smoother_out,
				alg1_lm2 = 1, 
				data_sample_alg = c(4,5)		#c(3, 7, 10, 14, 22)		
				) )
cbind(est1)
est1_mat = Format_pars(truePars = est1)
splineStart = smoother_out$splines_est[1,2:dim(smoother_out$data)[2]];names(splineStart) = 'x1'
out_est1 = solveLV(times = time_graywhales_corrected, pars = est1_mat, initState = splineStart, equations = Equations)

### slope estimates
#########################


#########################
### MAR estimates

MARest_orig =  fastMAR(
			data_vars = graywhales_corrected[,2],
			log_transform = FALSE,
			demeaned = FALSE,
			abstol_= 0.01,
			maxit_ = 100,
			includeNoise = FALSE,
			estInit = FALSE,
			tSteps_ = NULL,
			nsim_ = NULL
			) 
MARest_orig$MARSSci_output

MARest_log =  fastMAR(
			data_vars = graywhales_corrected[,2],
			log_transform = TRUE,
			demeaned = FALSE,
			abstol_= 0.01,
			maxit_ = 5000,
			includeNoise = FALSE,
			estInit = FALSE,
			tSteps_ = NULL,
			nsim_ = NULL
			) 
MARest_log$MARSSci_output


plot(graywhales_corrected[,1],graywhales_corrected[,2],xlab='t',cex.axis=1.5)
points(graywhales_corrected[,1],MARest_orig$MAR_TS_Est ,type='l',lwd=3)
points(graywhales_corrected[,1],MARest_log$MAR_TS_Est ,type='l',lty=2,lwd=3)
#if (i==1) {legend(1952,25000,legend=c('original data','noisy data','MAR no transform','MAR log transform'),pch=c(NA,20,NA,NA),lty=c(1,NA,1,2),col=c('lightgrey',varColor[i],varColor_light[i],varColor_light[i]),bty = "n",cex=1,lwd=c(3,NA,1,1))}

points( graywhales_corrected[,1],meanMaker(MARest_orig),type='l',col='purple',lty=1)
points( graywhales_corrected[,1],meanMaker(MARest_log),type='l',col='purple',lty=2)

### MAR estimates
#########################



#######################
### MAR estimates with smoothing

# smoothing the data for MAR
smoother_for_MAR_out = Smoother_for_MAR(	
						dataset = graywhales_rm.na,
						df = c(8),
						dataWeights = NULL,
						splineMethod = "fmm",
						draw = TRUE,
						log_spline = FALSE
						)
smoother_for_MAR_out
TS_smooth_data = smooth_TS = smoother_for_MAR_out$splines_est

# MAR estimates without data log transformation, with spline smooth data  
MARest1_smoothData =  fastMAR(
			data_vars = TS_smooth_data[,2],
			log_transform=FALSE,
			demeaned=FALSE,
			abstol_=0.01,
			maxit_ = 1000,
			includeNoise = FALSE,
			estInit = FALSE,
			tSteps_ = NULL,
			nsim_ = NULL
			) 
MARest1_smoothData$MARSSci_output


# MAR estimates with data log transformation, with spline smooth data  
MARest2_smoothData =  fastMAR(
			data_vars = TS_smooth_data[,2],
			log_transform=TRUE,
			demeaned=FALSE,
			abstol_=0.01,
			maxit_ = 1000,
			includeNoise = FALSE,
			estInit = FALSE,
			tSteps_ = NULL,
			nsim_ = NULL
			) 
MARest2_smoothData$MARSSci_output

### MAR estimates with smoothing
#######################



# setwd("C:\\Users\\doliv\\OneDrive\\2_2019_America\\2020\\20200123_MAR\\Camparison_LV_MAR\\8_methods_of_ecology_and_evolution_2sub\\paper_figures")     # P51
# tiff("fig4a.tiff", height = 15, width = 15, units = 'cm',compression = "lzw", res = 300)
par(mfrow=c(1,1),mar=c(5, 6, 4, 2) + 0.1)
plot(graywhales_rm.na,pch=1,col=colorPallet[1],xlab='years',ylab='number of graywhales',cex.lab=2,cex.axis=1.5)
points(graywhales_corrected[,1],MARest_orig$MAR_TS_Est, type = "l", lwd = 3,col=colorPallet[4])
points(graywhales_corrected[,1],MARest_log$MAR_TS_Est, type = "l", lwd = 3,col=colorPallet[5])
points(graywhales_corrected[,1],MARest1_smoothData$MAR_TS_Est, type = "l", lwd = 3,col='orange')
points(graywhales_corrected[,1],MARest2_smoothData$MAR_TS_Est, type = "l", lwd = 3,col='khaki')
points(out_est1[,1],out_est1[,2],type='l',col=colorPallet[3],lwd=3)
legend(1970,10000,legend=c('ALVI-MI method','MAR no transform','MAR log transform','MAR no transform with smoothing','MAR log transform with smoothing','data'),lty=c(1,1,1,1,1,NA),lwd=c(3,3,3,3,3,1),pch=c(NA,NA,NA,NA,NA,1),col=c(colorPallet[3],colorPallet[4],colorPallet[5],'orange','khaki',colorPallet[1]),bty = "n",cex=.8)
# dev.off()

cbind(graywhales_corrected,t(kem$states))
sum( (graywhales_rm.na[,2] - out_est1[out_est1[,1]%in%graywhales_rm.na[,1],2])^2 )
sum( (graywhales_rm.na[,2] - MARest_orig$MAR_TS_Est[!is.na(graywhales_corrected[,2])])^2 )
sum( (graywhales_rm.na[,2] - MARest_log$MAR_TS_Est[!is.na(graywhales_corrected[,2])])^2 )
sum( (graywhales_rm.na[,2] - MARest1_smoothData$MAR_TS_Est[TS_smooth_data[,1]%in%graywhales_rm.na[,1] ])^2 )
sum( (graywhales_rm.na[,2] - MARest2_smoothData$MAR_TS_Est[TS_smooth_data[,1]%in%graywhales_rm.na[,1] ])^2 )


##### graywhales
######################################################################################################################################



######################################################################################################################################
##### ivesDataByWeek

################################
### time, initState and truePars 

data(ivesDataByWeek)
ivesDataByWeek
varNames = c("Large Phyto","Small Phyto","Daphnia","Non-daphnia")

data_plank_NAs = ivesDataByWeek[,c("Study week","Large Phyto","Small Phyto","Daphnia","Non-daphnia")]
data_plank = data_plank_NAs[1,]
for (i in 2:dim(data_plank_NAs)[1])
	{
	if ( sum(is.na(data_plank_NAs[i,])) == 0 )
		{
		data_plank = rbind(data_plank,data_plank_NAs[i,])
		} else
		{
		#data_plank = rbind(data_plank,c(data_plank_NAs[i,1],.001,.001,.001,.001))
		}
	}

time_plank = 27:295
state_plank = data_plank[1,2:dim(data_plank)[2]]
truePars_plank = c(		
			a1 = 0.344,
			b11 = 0,
			b12 = -0.059,
			b13 = 0,
			b14 = 0,

			a2 = -0.236,
			b21 = 0.0005,
			b22 = 0,
			b23 = 0,
			b24 = 0,

			a3 = 0.344,
			b31 = 0,
			b32 = -0.059,
			b33 = 0,
			b34 = 0,

			a4 = -0.236,
			b41 = 0.0005,
			b42 = 0,
			b43 = 0,
			b44 = 0
      		)
truePars_mat_plank = Format_pars(truePars = truePars_plank)

### time, initState and truePars 
################################

out1_plank = solveLV(times = time_plank, 
				initState = state_plank, 
				pars = truePars_mat_plank, 
				equations = Equations)
# head(out1) #tail(out1,10) # edit(out1) # View(out1)

par(mfrow = c(2,2))
for (i in 1:length(state_plank))
	{
	plot(data_plank[,1],data_plank[,i+1],col='red',xlab='Years',ylab=paste('log(',varNames[i],')') )
	points(out1_plank[,1],out1_plank[,i+1], type = "l", lwd = 2)
	}







#########################
### slope estimates

smoother_out = Smoother(
				dataset = data_plank,
				draw = TRUE,
				data1_spline2 = 2, 
				smooth_type = 1,
				polyDegree012 = 1,
				df = c(5,5,5,5), #c(50,50,50,50),
				log_spline = TRUE
				)
plot(1,1,col='white')
legend(1,1,legend = c('data','splines','slopes'),pch = c(20,NA,NA),lty = c(0,1,1),col = c('grey','green','lightgreen'),bty = "n")

AlgPointFind(
		smoother_out_APF = smoother_out,
		dif_equations = Equations,
		matrixEq = TRUE
		)
# 	"100 % complete || DF = 15  ||  2 14 170511189.416  ||  2 14 114963636.776"
# log "100 % complete || DF = 15  ||  1  2 168599305.82   ||  1 18 114535854.508"
# log "100 % complete || DF =  3  ||  4  5 153703188.988  ||  2 20   2090624.691"

( est1 = LV_pars_finder(
				smooth_out = smoother_out,
				alg1_lm2 = 1, 
				data_sample_alg = c(3,26,35,49,70)		#c(3, 7, 10, 14, 22)	#c(3,26,35,49,70)		
				) )
cbind(est1)
est1_mat = Format_pars(truePars = est1)
splineStart = smoother_out$splines_est[1,2:dim(smoother_out$data)[2]];names(splineStart) = c('X1','X2','X3','X4')
out_est1 = solveLV(times = seq(27,295,1), pars = est1_mat, initState = splineStart, equations = Equations)

par(mfrow=c(2,2),mar=c(5, 6, 4, 2) + 0.1)
for (iii in 1:4)
	{
	plot(data_plank[,1],data_plank[,iii+1],pch=20,lwd=3,col='grey',xlab='t',ylab=names(data_plank)[iii+1])
	points(out_est1[,1],out_est1[,iii+1] ,type='l',col='cyan',lwd=3)
	}

### slope estimates
#########################


#########################
### MAR estimates

MARest_orig =  fastMAR(
			data_vars = data_plank_NAs[,2:5],
			log_transform = FALSE,
			demeaned = FALSE,
			abstol_= 0.01,
			maxit_ = 200,
			includeNoise = FALSE,
			estInit = FALSE,
			tSteps_ = NULL,
			nsim_ = 100
			) 
MARest_orig$MARSSci_output

MARest_log =  fastMAR(
			data_vars = data_plank_NAs[,2:5],
			log_transform = TRUE,
			demeaned = FALSE,
			abstol_= 0.01,
			maxit_ = 700,
			includeNoise = FALSE,
			estInit = TRUE,
			tSteps_ = NULL,
			nsim_ = 100
			) 
MARest_log$MARSSci_output


par(mfrow=c(2,2),mar=c(5, 6, 4, 2) + 0.1)
for (iii in 1:4)
	{
	plot(data_plank[,1],data_plank[,iii+1],pch=20,lwd=3,col='grey',xlab='t',ylab=colnames(data_plank)[iii+1],cex.lab=2,cex.axis=1.5)
	points(data_plank_NAs[,1],MARest_orig$MAR_TS_Est[,iii] ,type='l',col='green',lwd=3)
	points(data_plank_NAs[,1],MARest_log$MAR_TS_Est[,iii] ,type='l',col='green',lty=2,lwd=3)
	points(data_plank_NAs[,1],meanMaker(MARest_orig)[iii,],type='l',col='purple',lty=1)
	points(data_plank_NAs[,1],meanMaker(MARest_log)[iii,],type='l',col='purple',lty=2)
	}

### MAR estimates
#########################


# setwd("C:\\Users\\doliv\\OneDrive\\2_2019_America\\2020\\20200123_MAR\\Camparison_LV_MAR\\5_methods_of_ecology_and_evolution\\paper_figures")     # P51
# tiff("fig5b.tiff", height = 20, width = 20, units = 'cm',compression = "lzw", res = 300)
par(mfrow=c(2,2),mar=c(5, 6, 4, 2) + 0.1)
for (iii in 1:4)
	{
	plot(data_plank[,1],data_plank[,iii+1],pch=20,lwd=3,col='grey',xlab='t',ylab=colnames(data_plank)[iii+1],cex.lab=2,cex.axis=1.5)
	points(out_est1[,1],out_est1[,iii+1] ,type='l',col='cyan',lwd=3)
	points(data_plank_NAs[,1],MARest_orig$MAR_TS_Est[,iii] ,type='l',col='green',lwd=3)
	points(data_plank_NAs[,1],MARest_log$MAR_TS_Est[,iii] ,type='l',col='green',lty=2,lwd=3)
	if (iii == 4) { legend(50,1200,legend=c('data','slope est','MAR','MAR logTrans'),lty=c(NA,1,1,2),lwd=c(3,3,2,2),pch=c(20,NA,NA,NA),col=c('lightgrey','cyan','green','green'),bty = "n",cex=1) }
	}
# dev.off()


sum( (data_plank[,2:5] - out_est1[out_est1[,1]%in%data_plank[,1],2:5])^2 )
sum( (data_plank[,2:5] - MARest_orig$MAR_TS_Est[data_plank_NAs[,1]%in%data_plank[,1],])^2 )
sum( (data_plank[,2:5] - MARest_log$MAR_TS_Est[data_plank_NAs[,1]%in%data_plank[,1],])^2 )

##### ivesDataByWeek
######################################################################################################################################



######################################################################################################################################
##### wolves moose

#######
### MAR 

# isleRoyal = read.csv('C:\\Users\\dolivenca3\\OneDrive\\2_2019_America\\2020\\20200123_MAR\\Camparison_LV_MAR\\8_methods_of_ecology_and_evolution_2sub\\paper_scripts_GitHub\\Holme_isle_royale_data\\isleRoyale.csv')
isleRoyal = read.csv('C:\\Users\\doliv\\OneDrive\\2_2019_America\\2020\\20200123_MAR\\Camparison_LV_MAR\\8_methods_of_ecology_and_evolution_2sub\\paper_scripts_GitHub\\Holme_isle_royale_data\\isleRoyale.csv')
yr1960to2011 <- isleRoyal[, "year"] >= 1960 & isleRoyal[, "year"] <= 2011

royale.dat <- log(t(isleRoyal[yr1960to2011,2:3]))
x <- isleRoyal[, "year"]
y <- log(isleRoyal[, 2:3])

graphics::matplot(x, y,
	ylab = "Log count", xlab = "Year", type = "l",
	lwd = 3, bty = "L", col = "black"
	)
points(isleRoyal[2:53,1],royale.dat[1,])
points(isleRoyal[2:53,1],royale.dat[2,],pch=2)
legend("topleft", c("Wolf", "Moose"), lty = c(1, 2),pch = c(1,2), bty = "n")


# if missing values are in the data, they should be NAs
z.royale.dat <- zscore(royale.dat)
royale_stats = rbind(apply(royale.dat,MARGIN=1,FUN=mean),apply(royale.dat,MARGIN=1,FUN=sd))
colnames(royale_stats) = rownames(royale.dat); rownames(royale_stats) = c('mean','sd')
royale_stats

par(mfrow=c(1,1))
royale.model.2 <- list(
	Z = "identity", B = "unconstrained",
	Q = "diagonal and unequal", R = "zero", U = "zero"
	)
kem.2 <- MARSS(z.royale.dat, model = royale.model.2)
try(kem.2.ci1 <- MARSSparamCIs(kem.2,method="hessian",alpha=0.05,nboot=100),silent=T)

numVar=2
n = 100
MARparEst1 = cbind(matrix(kem.2$par$B, ncol=numVar),c(0,0),kem.2$par$Q)


MARvalEst1 = x_t = z.royale.dat[,1]
	for (i in 2:n)
		{
		x_tplus1 = MARparEst1[,1:numVar] %*% x_t + MARparEst1[,numVar+1] + 0*rnorm(n=numVar, mean = rep(0,numVar), sd = sqrt(MARparEst1[,(numVar+2)]))		#
		MARvalEst1 = rbind(MARvalEst1,t(x_tplus1))
		x_t = x_tplus1	
		}
MARvalEst1



par(mfrow=c(2,1))
ylimlim=rbind(c(NA,NA),c(0,100),c(0,3000))
for (i in 2:3)
	{
	plot(isleRoyal[,1],isleRoyal[,i],pch=i-1,col='grey',ylim=ylimlim[i,])
	points(isleRoyal[yr1960to2011,1],isleRoyal[yr1960to2011,i],pch=i-1)

	points(1960+0:(dim(MARvalEst1)[1]-1),exp(MARvalEst1[,i-1]*royale_stats[2,i-1]+royale_stats[1,i-1]),type='l',col='green')

	points(1960+0:(dim(MARvalEst1)[1]-1),exp(MARvalEst1[,i-1]*royale_stats[2,i-1]+royale_stats[1,i-1]-2*sqrt(MARparEst1[i-1,4])),type='l',col='green',lty=2)
	points(1960+0:(dim(MARvalEst1)[1]-1),exp(MARvalEst1[,i-1]*royale_stats[2,i-1]+royale_stats[1,i-1]+2*sqrt(MARparEst1[i-1,4])),type='l',col='green',lty=2)
	}


##### MAR
###########################################



#########################
### slope estimates

isleRoyal[yr1960to2011,]

smoother_out = Smoother(
				dataset = isleRoyal[yr1960to2011,],
				draw = TRUE,
				data1_spline2 = 2, 
				smooth_type = 1,
				df = c(15,15),
				log_spline = TRUE
				)
plot(1,1,col='white')
legend(1,1,legend = c('data','splines','slopes'),pch = c(20,NA,NA),lty = c(0,1,1),col = c('grey','green','lightgreen'),bty = "n")

# uncommewnt the next lines to find the point subsample
#AlgPointFind(
#		smoother_out_APF = smoother_out,
#		dif_equations = Equations,
#		matrixEq = TRUE
#		)
# "100 % complete || DF = 15 15  ||  33 36 39 5543504.488  ||  33 36 39 4039227.828"

( est1 = LV_pars_finder(
				smooth_out = smoother_out,
				alg1_lm2 = 1, 
				data_sample_alg = c(33, 36, 39)		#c(3, 7, 10, 14, 22)		
				) )
cbind(est1)
est1_mat = Format_pars(truePars = est1)
splineStart = smoother_out$splines_est[1,2:dim(smoother_out$data)[2]];names(splineStart) = c('x1','x2')
out_est1 = solveLV(times = isleRoyal[,1], pars = est1_mat, initState = splineStart, equations = Equations)

par(mfrow=c(2,1))
ylimlim=rbind(c(NA,NA),c(0,100),c(0,3000))
for (i in 2:3)
	{
	plot(isleRoyal[,1],isleRoyal[,i],pch=i-1,col='grey',ylim=ylimlim[i,])
	points(isleRoyal[yr1960to2011,1],isleRoyal[yr1960to2011,i],pch=i-1)

	points(1960+0:(dim(MARvalEst1)[1]-1),exp(MARvalEst1[,i-1]*royale_stats[2,i-1]+royale_stats[1,i-1]),type='l',col='green')

	points(1960+0:(dim(MARvalEst1)[1]-1),exp(MARvalEst1[,i-1]*royale_stats[2,i-1]+royale_stats[1,i-1]-2*sqrt(MARparEst1[i-1,4])),type='l',col='green',lty=2)
	points(1960+0:(dim(MARvalEst1)[1]-1),exp(MARvalEst1[,i-1]*royale_stats[2,i-1]+royale_stats[1,i-1]+2*sqrt(MARparEst1[i-1,4])),type='l',col='green',lty=2)

	points(out_est1[,1],out_est1[,i],type='l',col='cyan',lwd=3)
	}

### slope estimates
#########################


#########################
### MAR estimates

MARest_orig =  fastMAR(
			data_vars = isleRoyal[yr1960to2011,2:3],
			log_transform = FALSE,
			demeaned = FALSE,
			abstol_= 0.01,
			maxit_ = 500,
			includeNoise = FALSE,
			estInit = FALSE,
			tSteps_ = 100,
			nsim_ = 100
			) 
MARest_orig$MARSSci_output

MARest_log =  fastMAR(
			data_vars = isleRoyal[yr1960to2011,2:3],
			log_transform = TRUE,
			demeaned = TRUE,
			abstol_= 0.01,
			maxit_ = 5000,
			includeNoise = FALSE,
			estInit = FALSE,
			tSteps_ = 100,
			nsim_ = 100
			) 
MARest_log$MARSSci_output


par(mfrow=c(2,1))
ylimlim=rbind(c(NA,NA),c(0,100),c(0,3000))
for (i in 2:3)
	{
	plot(isleRoyal[,1],isleRoyal[,i],pch=i-1,col='grey',ylim=ylimlim[i,])
	points(isleRoyal[yr1960to2011,1],isleRoyal[yr1960to2011,i],pch=i-1)

	points(1960+0:(dim(MARvalEst1)[1]-1),exp(MARvalEst1[,i-1]*royale_stats[2,i-1]+royale_stats[1,i-1]),type='l',col='green')

	points(1960+0:(dim(MARvalEst1)[1]-1),exp(MARvalEst1[,i-1]*royale_stats[2,i-1]+royale_stats[1,i-1]-2*sqrt(MARparEst1[i-1,4])),type='l',col='green',lty=2)
	points(1960+0:(dim(MARvalEst1)[1]-1),exp(MARvalEst1[,i-1]*royale_stats[2,i-1]+royale_stats[1,i-1]+2*sqrt(MARparEst1[i-1,4])),type='l',col='green',lty=2)

	points(out_est1[,1],out_est1[,i],type='l',col='cyan',lwd=3)

	points(isleRoyal[,1],MARest_orig$MAR_TS_Est[1:61,i-1] ,type='l',lwd=3)
	points(isleRoyal[,1],MARest_log$MAR_TS_Est[1:61,i-1] ,type='l',lty=2,lwd=3)
	}


#if (i==1) {legend(1952,25000,legend=c('original data','noisy data','MAR no transform','MAR log transform'),pch=c(NA,20,NA,NA),lty=c(1,NA,1,2),col=c('lightgrey',varColor[i],varColor_light[i],varColor_light[i]),bty = "n",cex=1,lwd=c(3,NA,1,1))}

### MAR estimates
#########################


#######################
### MAR estimates with smoothing

# smoothing the data for MAR
smoother_for_MAR_out = Smoother_for_MAR(	
						dataset = isleRoyal[yr1960to2011,],
						df = c(15,15),
						dataWeights = NULL,
						splineMethod = "fmm",
						draw = TRUE,
						log_spline = FALSE
						)
smoother_for_MAR_out
TS_smooth_data = smooth_TS = smoother_for_MAR_out$splines_est

# MAR estimates without data log transformation, with spline smooth data  
MARest1_smoothData =  fastMAR(
			data_vars = TS_smooth_data[,2:3],
			log_transform=FALSE,
			demeaned=FALSE,
			abstol_=0.01,
			maxit_ = 1000,
			includeNoise = FALSE,
			estInit = FALSE,
			tSteps_ = 100,
			nsim_ = 100
			) 
MARest1_smoothData$MARSSci_output


# MAR estimates with data log transformation, with spline smooth data  
MARest2_smoothData =  fastMAR(
			data_vars = TS_smooth_data[,2:3],
			log_transform=TRUE,
			demeaned=TRUE,
			abstol_=0.01,
			maxit_ = 1000,
			includeNoise = FALSE,
			estInit = FALSE,
			tSteps_ = 100,
			nsim_ = 100
			) 
MARest2_smoothData$MARSSci_output

### MAR estimates with smoothing
#######################

# setwd("C:\\Users\\doliv\\OneDrive\\2_2019_America\\2020\\20200123_MAR\\Camparison_LV_MAR\\8_methods_of_ecology_and_evolution_2sub\\paper_figures")     # P51
# tiff("fig4b.tiff", height = 30, width = 15, units = 'cm',compression = "lzw", res = 300)
par(mfrow=c(2,1),mar=c(5, 6, 4, 2) + 0.1)
ylimlim=rbind(c(NA,NA),c(0,100),c(0,3000))
ylablab=c(NA,'wolves','moose')
for (i in 2:3)
	{
	plot(isleRoyal[,1],isleRoyal[,i],col='grey',ylim=ylimlim[i,],xlab='years',ylab=ylablab[i],cex.lab=2,cex.axis=2)
	points(isleRoyal[yr1960to2011,1],isleRoyal[yr1960to2011,i])

	points(isleRoyal[,1],MARest_orig$MAR_TS_Est[1:61,i-1], type = "l", lwd = 3,col='darkgreen')
	points(isleRoyal[,1],MARest_log$MAR_TS_Est[1:61,i-1], type = "l", lwd = 3,col='green')
	
	#points(1960+0:(dim(MARvalEst1)[1]-1),exp(MARvalEst1[,i-1]*royale_stats[2,i-1]+royale_stats[1,i-1]),type='l',col='green',lwd=3)
	#points(1960+0:(dim(MARvalEst1)[1]-1),exp(MARvalEst1[,i-1]*royale_stats[2,i-1]+royale_stats[1,i-1]-2*sqrt(MARparEst1[i-1,4])),type='l',col='green',lty=2)
	#points(1960+0:(dim(MARvalEst1)[1]-1),exp(MARvalEst1[,i-1]*royale_stats[2,i-1]+royale_stats[1,i-1]+2*sqrt(MARparEst1[i-1,4])),type='l',col='green',lty=2)

	points(isleRoyal[,1],MARest1_smoothData$MAR_TS_Est[1:61,i-1], type = "l", lwd = 3,col='orange')
	points(isleRoyal[,1],MARest2_smoothData$MAR_TS_Est[1:61,i-1], type = "l", lwd = 3,col='khaki')

	points(out_est1[,1],out_est1[,i],type='l',col='blue',lwd=3)

	if (i==2) { legend(1960,100,legend=c('data','ALVI-MI method','MAR no transformation','MAR log transformation','MAR no transformation with smoothing','MAR log transformation with smoothing'),lty=c(NA,1,1,1,1,1),lwd=c(1,3,3,3,3,3),pch=c(1,NA,NA,NA,NA,NA),col=c('black','blue','darkgreen','green','orange','khaki'),bty = "n",cex=1) }
	}

# dev.off()


cbind(graywhales_corrected,t(kem$states))
sum( (isleRoyal - out_est1[,1:3])^2 )
sum( (isleRoyal[,2:3] - MARest_orig$MAR_TS_Est[1:61,])^2 )
sum( (isleRoyal[,2:3] - MARest_log$MAR_TS_Est[1:61,])^2 )
sum( (isleRoyal[,2:3] - MARest1_smoothData$MAR_TS_Est[1:61,])^2 )
sum( (isleRoyal[,2:3] - MARest2_smoothData$MAR_TS_Est[1:61,])^2 )




sum( (isleRoyal[,2:3] - cbind(exp(MARvalEst1[,1]*royale_stats[2,1]+royale_stats[1,1]),exp(MARvalEst1[,2]*royale_stats[2,2]+royale_stats[1,2])))^2 )

##### wolves moose
######################################################################################################################################


