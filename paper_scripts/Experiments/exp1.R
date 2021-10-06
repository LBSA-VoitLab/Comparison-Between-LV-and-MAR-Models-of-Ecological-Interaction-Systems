# Please uncomment the next few lines to find the best subsample for this dataset





smoother_out = Smoother(
				dataset = noisy_data3,
				draw = TRUE,
				data1_spline2 = 2, 
				smooth_type = 1,
				df = c(5,8,11,5),	
				log_spline = TRUE
				)

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


############################################################################################################################
##### MAR for noisy data 3

TS_noisy_data = cbind(1:100,matrix(rep(NA,100*4),ncol=4))
for (i in 1:100)
	{
	if (i%in%smooth_TS[,1])
		{
		TS_noisy_data[i,] = smooth_TS[smooth_TS[,1]==i,] 
		} 
	}

MARest1_splines =  fastMAR(
			data_vars = TS_noisy_data[,2:5],
			log_transform=FALSE,
			demeaned=TRUE,
			abstol_=0.5,
			maxit_ = 3000,
			includeNoise = FALSE,
			estInit = FALSE,
			tSteps_ = NULL,
			nsim_ = 100
			) 
MARest1_splines$MARSSci_output

MARest2_splines =  fastMAR(
			data_vars = TS_noisy_data[,2:5],
			log_transform=TRUE,
			demeaned=TRUE,
			abstol_=0.5,
			maxit_ = 3000,
			includeNoise = FALSE,
			estInit = FALSE,
			tSteps_ = NULL,
			nsim_ = 100
			) 
MARest2_splines$MARSSci_output

par(mfrow=c(2,2))
for ( i in 1:(dim(TS_noisy_data)[2]-1) )
	{
	plot(out1[,1],out1[,i+1],type='l',col='lightgrey',lwd=3,ylim=c(0,3),ylab = varNames[i])

	# for (ii in 1:100) { points(MARest1$MARSSsim_data[i,,ii],type='l',col='grey') }
	points(meanMaker(MARest1_splines)[i,],type='l',col = 'purple')

	# for (ii in 1:100) { points(MARest2$MARSSsim_data[i,,ii],type='l',col='grey') }
	points(meanMaker(MARest2_splines)[i,],type='l',col = 'purple',lty=2)


	points(out1[,1],out1[,i+1],type='l',col='lightgrey',lwd=3,ylim=c(0,3),ylab = varNames[i])
	points(TS_noisy_data[,i+1],pch=20,col=varColor[i],xlab = 'n')
	points(MARest1_splines$MAR_TS_Est[,i],col='green',type='l',lty=1)
	points(MARest2_splines$MAR_TS_Est[,i],col='green',type='l',lty=2)
	if (i==1) {legend(40,1.8,legend=c('original data','noisy data','MAR no transform','MAR log transform'),pch=c(NA,20,NA,NA),lty=c(1,NA,1,2),col=c('lightgrey',varColor[i],varColor_light[i],varColor_light[i]),bty = "n",cex=1,lwd=c(3,NA,1,1))}
	}

MARest1_orig = cbind(1:dim(MARest1_splines$MAR_TS_Est)[1],MARest1_splines$MAR_TS_Est)
MARest1_log = cbind(1:dim(MARest2_splines$MAR_TS_Est)[1],MARest2_splines$MAR_TS_Est)

sum( ( out1[,2:5] - MARest1_orig[,2:5] )^2 )
sum( ( out1[,2:5] - MARest1_log[,2:5] )^2 )

##### MAR for noisy data 3
###############################################################################################################

