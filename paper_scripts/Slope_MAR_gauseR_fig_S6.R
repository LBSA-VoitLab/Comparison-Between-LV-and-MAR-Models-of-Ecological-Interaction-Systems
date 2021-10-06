

rm(list = ls())


#install.packages("gauseR")
library(gauseR)

# data
# check http://127.0.0.1:13667/library/gauseR/html/00Index.html
# or hit index in the bottom of ?gauseR




##############################################################################
### data - fig 1

# load Gause competition data
data(gause_1934_science_f02_03)
# extract monoculture data for P.c.
Pcmono<-gause_1934_science_f02_03[gause_1934_science_f02_03$Treatment=="Pc",]

# calculate lag and per-capita growth
lagged_data_Pc <- get_lag(x=Pcmono$Volume_Species1,
                          time = Pcmono$Day)
Pcmono$dNNdt_Pc <- percap_growth(x=lagged_data_Pc$x, laggedx=lagged_data_Pc$laggedx,
                                 dt=lagged_data_Pc$dt)

# fit linear model to get dNNdt ~ r + s*N
mod_Pc<-lm(dNNdt_Pc~Volume_Species1, Pcmono)
rsn_pars<-coef(mod_Pc)

# transform into logistic growth parameters
logistic_pars<-c(r=unname(rsn_pars["(Intercept)"]),
     K=unname(-rsn_pars["(Intercept)"]/rsn_pars["Volume_Species1"]))

#fit with nls, using linear model estimates as starting values for parameters
nls_mod<-nls(Volume_Species1~get_logistic(time = Day, N0, r, K),
             data=Pcmono,
             start=c(N0=unname(Pcmono$Volume_Species1[which.min(Pcmono$Day)]),
             r=unname(logistic_pars["r"]), K=unname(logistic_pars["K"])))
summary(nls_mod)

# plot results
plot(Volume_Species1~Day, Pcmono, type="b", ylab="P. caudatum Volume")
timesq<-seq(0, 30, length=100)
Ntest<-get_logistic(time = timesq, N0=coef(nls_mod)["N0"], r=coef(nls_mod)["r"],
       K=coef(nls_mod)["K"])

lines(timesq, Ntest, col="red")



fig1data = Pcmono[,c(3,5)]
plot(fig1data[,1],fig1data[,2],xlab = "Days",ylab = "Paramecium caudatum std. volume")

### data - fig 1
##############################################################################



##############################################################################
### data - fig 2

###############################
### Figure 2 original code http://127.0.0.1:16583/library/gauseR/html/lv_interaction.html

# load data from competition experiment
data(gause_1934_science_f02_03)

# subset data to include just mixtures
mixturedata<-gause_1934_science_f02_03[gause_1934_science_f02_03$Treatment=="Mixture",]

# get time-lagged observations for each species
Pc_lagged<-get_lag(x = mixturedata$Volume_Species1, time = mixturedata$Day)
Pa_lagged<-get_lag(x = mixturedata$Volume_Species2, time = mixturedata$Day)

# calculate per-capita growth rates
Pc_dNNdt<-percap_growth(x = Pc_lagged$x, laggedx = Pc_lagged$laggedx, dt = Pc_lagged$dt)
Pa_dNNdt<-percap_growth(x = Pa_lagged$x, laggedx = Pa_lagged$laggedx, dt = Pa_lagged$dt)

# fit linear models to dNNdt, based on average
# abundances between current and lagged time steps
Pc_mod_dat<-data.frame(Pc_dNNdt=Pc_dNNdt, Pc=Pc_lagged$laggedx, Pa=Pa_lagged$laggedx)
mod_comp_Pc<-lm(Pc_dNNdt~Pc+Pa, data=Pc_mod_dat)

Pa_mod_dat<-data.frame(Pa_dNNdt=Pa_dNNdt, Pa=Pa_lagged$laggedx, Pc=Pc_lagged$laggedx)
mod_comp_Pa<-lm(Pa_dNNdt~Pa+Pc, data=Pa_mod_dat)

# model summaries
summary(mod_comp_Pc)
summary(mod_comp_Pa)

# extract parameters
# note - linear regressions give us dynamics in the form:
# dni/nidt ~ (Intercept) + (n1_slope) * n1 + (n2_slope) n2
# and thus:
# dni/dt = n1*((Intercept) + (n1_slope) * n1 + (n2_slope) n2)

# growth rates
r1 <- unname(coef(mod_comp_Pc)["(Intercept)"])
r2 <- unname(coef(mod_comp_Pa)["(Intercept)"])

# self-limitation
a11 <- unname(coef(mod_comp_Pc)["Pc"])
a22 <- unname(coef(mod_comp_Pa)["Pa"])

# effect of Pa on Pc
a12 <- unname(coef(mod_comp_Pc)["Pa"])
# effect of Pc on Pa
a21 <- unname(coef(mod_comp_Pa)["Pc"])

# run ODE:
# make paramter vector:
parms <- c(r1, r2, a11, a12, a21, a22)
initialN <- c(1, 1)
out <- deSolve::ode(y=initialN, times=1:25, func=lv_interaction, parms=parms)
matplot(out[,1], out[,-1], type="l",
   xlab="time", ylab="N", col=c("black","red"), lty=c(1,3), lwd=2, ylim=c(0, 150))
legend("topleft", c("Pc", "Pa"), col=c(1,2), lwd=2, lty=c(1,3))

# now, plot in points from data
points(mixturedata$Day, mixturedata$Volume_Species1, col=1)
points(mixturedata$Day, mixturedata$Volume_Species2, col=2)

### Figure 2 original code http://127.0.0.1:16583/library/gauseR/html/lv_interaction.html
###############################

fig2data<-gause_1934_science_f02_03[gause_1934_science_f02_03$Treatment=="Mixture",]
points(fig2data[,3],fig2data[,5],ylim=c(0,150),pch=6)
points(fig2data[,3],fig2data[,7],col='red',pch=2)

### data - fig 2
##############################################################################



##############################################################################
### data - fig 3

###############################
### Figure 3 original code http://127.0.0.1:16583/library/gauseR/html/lv_optim.html

# load data from competition experiment
data(gause_1934_book_f32)

# keep all data - no separate treatments exist for this experiment
predatorpreydata<-gause_1934_book_f32

# get time-lagged observations for each species
prey_lagged<-get_lag(x = predatorpreydata$Individuals_Prey, time = predatorpreydata$Day)
predator_lagged<-get_lag(x = predatorpreydata$Individuals_Predator, time = predatorpreydata$Day)

# calculate per-capita growth rates
prey_dNNdt<-percap_growth(x = prey_lagged$x, laggedx = prey_lagged$laggedx, dt = prey_lagged$dt)
predator_dNNdt<-percap_growth(x = predator_lagged$x,
      laggedx = predator_lagged$laggedx, dt = predator_lagged$dt)

# fit linear models to dNNdt, based on average
# abundances between current and lagged time steps
prey_mod_dat<-data.frame(prey_dNNdt=prey_dNNdt, prey=prey_lagged$laggedx,
      predator=predator_lagged$laggedx)
mod_prey<-lm(prey_dNNdt~prey+predator, data=prey_mod_dat)

predator_mod_dat<-data.frame(predator_dNNdt=predator_dNNdt,
      predator=predator_lagged$laggedx, prey=prey_lagged$laggedx)
mod_predator<-lm(predator_dNNdt~predator+prey, data=predator_mod_dat)

# model summaries
summary(mod_prey)
summary(mod_predator)

# extract parameters
# growth rates
r1 <- unname(coef(mod_prey)["(Intercept)"])
r2 <- unname(coef(mod_predator)["(Intercept)"])

# self-limitation
a11 <- unname(coef(mod_prey)["prey"])
a22 <- unname(coef(mod_predator)["predator"])

# effect of Pa on Pc
a12 <- unname(coef(mod_prey)["predator"])
# effect of Pc on Pa
a21 <- unname(coef(mod_predator)["prey"])

# run ODE:
# make paramet

### Figure 3 original code http://127.0.0.1:16583/library/gauseR/html/lv_optim.html
###############################


gause_1934_book_f32 		# Didinium/Paramecium predator/prey experiment
fig3data = gause_1934_book_f32[,c(3,5,7)] 

plot(fig3data[,1],fig3data[,2],xlab='Days',ylab='dot Prey | sqr Predator')
points(fig3data[,1],fig3data[,3],col='blue',pch=0)

### data - fig 3
##############################################################################



##############################################################################
### data - fig 4

mclaren_1994_f03	# Wolf, Moose, and Fir dynamics from Isle Royale
fig4data = cbind(Year = mclaren_1994_f03[1:34,3],Wolves = mclaren_1994_f03[1:34,6],Moose = mclaren_1994_f03[37:70,6],Fir_trees = mclaren_1994_f03[107:140,c(5)])  #73:106

par(mfrow = c(2,2))
plot(fig4data[,1],fig4data[,2],col='red',xlab='Year',ylab='Wolves')
plot(fig4data[,1],fig4data[,3],pch=2,xlab='Year',ylab='Moose')
plot(fig4data[,1],fig4data[,4],col='blue',pch=0,xlab='Year',ylab='Fir_tree_rings (tree growth)')

### data - fig 4
##############################################################################



##############################################################################
### data - fig 5

huffaker_1963
fig5data = cbind(huffaker_1963[1:60,4:5],huffaker_1963[61:120,5])
colnames(fig5data) = c('Weeks','Eotetranychus_sexmaculatus','Typhlodromus_occidentalis')
fig5data[1,3] = fig5data[2,3] = fig5data[3,3]
par(mfrow = c(2,1))
plot(fig5data[,1],fig5data[,2],xlab='Weeks',ylab='Eotetranychus sexmaculatus')
plot(fig5data[,1],fig5data[,3],col='blue',pch=0,xlab='Weeks',ylab='Typhlodromus occidentalis')

### data - fig 5
##############################################################################



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
						varMeans = apply(smoother_out_APF$data[,2:(nVars+1)],2,mean)													# get the means of each variable
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
##### fig1

################################
### time, initState and truePars  

fig1data				# data
cbind(timesq, Ntest)		# Muhlbauer solution
t_fig1 = seq(2,24,1)		# time
state_fig1 = c(1.595)		# state
names(state_fig1) = c('x1')	# varNames
truePars_fig1 = c(		
			a1 = 1.259,
			b11 = -0.005
		      )
truePars_mat_fig1 = Format_pars(truePars = truePars_fig1)	# pars

### time, initState and truePars 
################################


#########################
### slope estimates

smoother_out = Smoother(
				dataset = fig1data,
				draw = TRUE,
				data1_spline2 = 2, 
				smooth_type = 1,
				df = 8,	
				log_spline = FALSE
				)
plot(1,1,col='white')
legend(1,1,legend = c('data','splines','slopes'),pch = c(20,NA,NA),lty = c(0,1,1),col = c('grey','green','lightgreen'),bty = "n")

( est1 = LV_pars_finder(
				smooth_out = smoother_out,
				alg1_lm2 = 2
				) )
cbind(truePars_fig1,est1,absDif = abs(truePars_fig1-est1))
est1_mat = Format_pars(truePars = est1)
splineStart = smoother_out$splines_est[1,2];names(splineStart) = 'x1'
out_est1 = solveLV(times = t_fig1, pars = est1_mat, initState = splineStart, equations = Equations)

par(mfrow=c(1,1),mar=c(5, 6, 4, 2) + 0.1)
plot(fig1data[,1],fig1data[,2],xlab = "Days",ylab = "Paramecium caudatum std. volume",cex.lab=2,cex.axis=2) #
lines(timesq, Ntest)
points(smoother_out$splines_est[,1],smoother_out$splines_est[,2],type='l',col='green')
points(out_est1[,1],out_est1[,2],type='l',col='cyan',lwd=3)

### slope estimates
#########################

#######################
### MAR estimates

MARest_orig =  fastMAR(
			data_vars = fig1data[,2],
			log_transform = FALSE,
			demeaned = TRUE,
			abstol_= 0.01,
			maxit_ = 100,
			includeNoise = FALSE,
			estInit = FALSE,
			tSteps_ = NULL,
			nsim_ = 100
			) 
MARest_orig$MARSSci_output


MARest_log =  fastMAR(
			data_vars = fig1data[,2],
			log_transform = TRUE,
			demeaned = TRUE,
			abstol_= 0.01,
			maxit_ = 200,
			includeNoise = FALSE,
			estInit = FALSE,
			tSteps_ = NULL,
			nsim_ = 100
			) 
MARest_log$MARSSci_output


plot(fig1data[,1],fig1data[,2],xlab='t',cex.axis=1.5)
points(2:24,MARest_orig$MAR_TS_Est ,type='l',lwd=3)
points(2:24,MARest_log$MAR_TS_Est ,type='l',lty=2,lwd=3)
#if (i==1) {legend(40,1.8,legend=c('original data','noisy data','MAR no transform','MAR log transform'),pch=c(NA,20,NA,NA),lty=c(1,NA,1,2),col=c('lightgrey',varColor[i],varColor_light[i],varColor_light[i]),bty = "n",cex=1,lwd=c(3,NA,1,1))}

points( 2:24,meanMaker(MARest_orig),type='l',col='purple',lty=1)
points( 2:24,meanMaker(MARest_log),type='l',col='purple',lty=2)

### MAR estimates
#######################


#######################
### MAR estimates with smoothing

# smoothing the data for MAR
smoother_for_MAR_out = Smoother_for_MAR(	
						dataset = fig1data,
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
			data_vars = TS_smooth_data[,2],
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

### MAR estimates with smoothing
#######################


# setwd("C:\\Users\\doliv\\OneDrive\\2_2019_America\\2020\\20200123_MAR\\Camparison_LV_MAR\\8_methods_of_ecology_and_evolution_2sub\\paper_figures")     # P51
# tiff("FigS6a.tiff", height = 7.5, width = 7.5, units = 'cm',compression = "lzw", res = 300)
par(mfrow=c(1,1))		#,mar=c(5, 6, 4, 2) + 0.1
plot(fig1data[,1], fig1data[,2], col=colorPallet[1], main='', xlab="Days", ylab="", cex.main=.8, cex.lab=1.5, cex.axis=1.5) #
lines(timesq, Ntest,col=colorPallet[2],lwd=3)
points(out_est1[,1],out_est1[,2],type='l',col=colorPallet[3],lwd=3)
points(2:(length(MARest_orig$MAR_TS_Est)+1),MARest_orig$MAR_TS_Est,type='l',col=colorPallet[4],lwd=3)
points(2:(length(MARest_orig$MAR_TS_Est)+1),MARest_log$MAR_TS_Est,type='l',col=colorPallet[5],lwd=3)
points(2:(length(MARest_orig$MAR_TS_Est)+1),MARest1_smoothData$MAR_TS_Est,type='l',col='orange',lwd=3)
points(2:(length(MARest_orig$MAR_TS_Est)+1),MARest2_smoothData$MAR_TS_Est,type='l',col='khaki',lwd=3)
points(out_est1[,1],out_est1[,2],type='l',col=colorPallet[3],lwd=3)
mtext(expression(paste(italic('P. caudatum'))),side=3,line=1,cex=1.5)
mtext("std. volume",side=2,line=2.5,cex=1.5)
#legend(13,75,legend = c("data",'Muhlbauer est','Slope est','MAR','MAR logTrans'),lty = c(0,1,1,1,1),pch = c(1,NA,NA,NA,NA),col = colorPallet,text.col=c('black','black','black','black','black'),bty = "n",cex=1,lwd=c(1,3,3,3,3))
# dev.off()


Muhlbauer = cbind(timesq, Ntest)
cbind(1:length(Muhlbauer[,1]),round(Muhlbauer[,1],0))
Muhlbauer = Muhlbauer[c(3,6,10,13,16,20,23,26,30,33,36,39,43,46,49,53,56,59,63,66,69,72,76,79,82,86,89,92,96,99),]
Muhlbauer[,1] = round(Muhlbauer[,1],0)
Muhlbauer = Muhlbauer[Muhlbauer[,1]%in%fig1data[,1],]
sum( (fig1data[,2]-Muhlbauer[,2])^2 )

sum( (fig1data[,2]-out_est1[,2])^2 )
sum( (fig1data[,2]-MARest_orig$MAR_TS_Est)^2 )
sum( (fig1data[,2]-MARest_log$MAR_TS_Est)^2 )
sum( (fig1data[,2]-MARest1_smoothData$MAR_TS_Est)^2 )
sum( (fig1data[,2]-MARest2_smoothData$MAR_TS_Est)^2 )

##### fig1
######################################################################################################################################



######################################################################################################################################
##### fig2

################################
### time, initState and truePars  

t_fig2 = seq(2,24,1)
fig2data = fig2data[,c(3,5,7)]
state_fig2 = c(1.595,1.632)
names(state_fig2) = c('x1','x2')
truePars_fig2 = c(		
			a1 = 1.259,
			b11 = -0.005,
			b12 = -0.008,

			a2 = 1.026,
			b21 = -0.002,
			b22 = -0.007
		      )
truePars_mat_fig2 = Format_pars(truePars = truePars_fig2)

### time, initState and truePars 
################################


#########################
### slope estimates

smoother_out = Smoother(
				dataset = fig2data,
				draw = TRUE,
				data1_spline2 = 2, 
				smooth_type = 1,
				df = c(10,7),	
				log_spline = TRUE
				)
plot(1,1,col='white')
legend(1,1,legend = c('data','splines','slopes'),pch = c(20,NA,NA),lty = c(0,1,1),col = c('grey','green','lightgreen'),bty = "n")

( est1 = LV_pars_finder(
				smooth_out = smoother_out,
				alg1_lm2 = 2
				) )
cbind(truePars_fig2,est1,absDif = abs(truePars_fig2-est1))
est1_mat = Format_pars(truePars = est1)
splineStart = smoother_out$splines_est[1,2:3];names(splineStart) = c("Pc", "Pa")
out_est1 = solveLV(times = t_fig2, pars = est1_mat, initState = splineStart, equations = Equations)

par(mfrow=c(1,1),mar=c(5, 6, 4, 2) + 0.1)
matplot(out[,1], out[,-1], type="l",
   xlab="Days", ylab="Std. Volume", col=c("black","red"), lty=c(1,3), lwd=2, ylim=c(0, 150),cex.lab=2,cex.axis=2)
legend("topleft", c("Pc", "Pa"), col=c(1,2),pch=c(1,2), lwd=2, lty=c(1,3),bty = "n",cex=1.5)
points(mixturedata$Day, mixturedata$Volume_Species1, col=1)
points(mixturedata$Day, mixturedata$Volume_Species2, col=2, pch=2)

for (v in 2:3)
	{
	points(smoother_out$splines_est[,1],smoother_out$splines_est[,v],type='l',col='green')
	points(fig2data[,1],fig2data[,v],ylim=c(0,150))
	points(out_est1[,1],out_est1[,v],type='l',col='cyan',lwd=3)
	}

### slope estimates
#########################

#######################
### MAR estimates

MARest_orig =  fastMAR(
			data_vars = fig2data[,2:3],
			log_transform = FALSE,
			demeaned = FALSE,
			abstol_= 0.01,
			maxit_ = 500,
			includeNoise = FALSE,
			estInit = FALSE,
			tSteps_ = NULL,
			nsim_ = 100
			) 
MARest_orig$MARSSci_output

MARest_log =  fastMAR(
			data_vars = fig2data[,2:3],
			log_transform = TRUE,
			demeaned = FALSE,
			abstol_= 0.01,
			maxit_ = 500,
			includeNoise = FALSE,
			estInit = FALSE,
			tSteps_ = NULL,
			nsim_ = 100
			) 
MARest_log$MARSSci_output

par(mfrow=c(1,1),mar=c(5, 6, 4, 2) + 0.1)
matplot(out[,1], out[,-1], type="l",
   xlab="Days", ylab="Std. Volume", col=c("black","red"), lty=c(1,3), lwd=2, ylim=c(0, 150),cex.lab=2,cex.axis=2)
legend("topleft", c("Pc", "Pa"), col=c(1,2),pch=c(1,2), lwd=2, lty=c(1,3),bty = "n",cex=1.5)
points(mixturedata$Day, mixturedata$Volume_Species1, col=1)
points(mixturedata$Day, mixturedata$Volume_Species2, col=2, pch=2)

for (v in 2:3)
	{
	points(2:24,MARest_orig$MAR_TS_Est[,v-1] ,type='l',col='green',lwd=3)
	points(2:24,MARest_log$MAR_TS_Est[,v-1] ,type='l',col='green',lty=2,lwd=3)
	points(2:24,meanMaker(MARest_orig)[v-1,],type='l',col='purple',lty=1)
	points(2:24,meanMaker(MARest_log)[v-1,],type='l',col='purple',lty=1)
	}

### MAR estimates
#######################


#########################
### MAR estimates with smoothing

# smoothing the data for MAR
smoother_for_MAR_out = Smoother_for_MAR(	
						dataset = fig2data,
						df = c(10,7),
						dataWeights = NULL,
						splineMethod = "fmm",
						draw = TRUE,
						log_spline = TRUE
						)
smoother_for_MAR_out
TS_smooth_data = smooth_TS = smoother_for_MAR_out$splines_est

# MAR estimates without data log transformation, with spline smooth data  
MARest1_smoothData =  fastMAR(
			data_vars = TS_smooth_data[,2:3],
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
			data_vars = TS_smooth_data[,2:3],
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

### MAR estimates with smoothing
#########################


# tiff("FigS6b.tiff", height = 15, width = 7.5, units = 'cm',compression = "lzw", res = 300)
par(mfrow=c(2,1))		 #,mar=c(5, 6, 4, 2) + 0.1  \n
main_lab = c( expression(paste(italic('P. caudatum'))),expression(paste(italic('P. aurelia'))) )
data_fig2 = mixturedata[,c(3,5,7)]
for (v in 2:3)	 
	{
	plot(data_fig2[,1], data_fig2[,v], col=colorPallet[1], main='', xlab="Days", ylab="", cex.main=.8 , cex.lab=1.5, cex.axis=1.5 )
	points(out[,1], out[,v], type='l', col=colorPallet[2], lwd=3)
	points(2:24,MARest_orig$MAR_TS_Est[,v-1], type='l', col=colorPallet[4], lwd=3)
	points(2:24,MARest_log$MAR_TS_Est[,v-1], type='l', col=colorPallet[5], lwd=3)
	points(2:24,MARest1_smoothData$MAR_TS_Est[,v-1],type='l',col='orange',lwd=3)
	points(2:24,MARest2_smoothData$MAR_TS_Est[,v-1],type='l',col='khaki',lwd=3)
	points(out_est1[,1], out_est1[,v], type='l', col=colorPallet[3], lwd=3)
	mtext(main_lab[v-1],side=3,line=1,cex=1.5)
	mtext("Std. Volume",side=2,line=2.5,cex=1.5)
	}
#legend(17,110, c('data','Muhlbauer est','Slope est','MAR','MAR logTrans'), col=colorPallet,pch=c(1,NA,NA,NA,NA), lwd=c(NA,3,3,3,3), lty=c(NA,1,1,1,1),bty = "n",cex=1)
# dev.off()

Muhlbauer2 = out
sum( (fig2data[,2:3]-Muhlbauer2[,2:3])^2 )

sum( (fig2data[,2:3]-out_est1[,2:3])^2 )
sum( (fig2data[,2:3]-MARest_orig$MAR_TS_Est)^2 )
sum( (fig2data[,2:3]-MARest_log$MAR_TS_Est)^2 )
sum( (fig2data[,2:3]-MARest1_smoothData$MAR_TS_Est)^2 )
sum( (fig2data[,2:3]-MARest2_smoothData$MAR_TS_Est)^2 )

##### fig2
######################################################################################################################################



######################################################################################################################################
##### fig3

################################
### time, initState and truePars  

t_fig3 = seq(1,17,1)
fig3data

state_fig3 = c(1.655,0.117)
names(state_fig3) = c('x1','x2')
truePars_fig3 = c(		
			a1 = 1.099,
			b11 = -0.013,
			b12 = -0.078,

			a2 = -0.89,
			b21 = 0.084,
			b22 = -0.002
		      )
truePars_mat_fig3 = Format_pars(truePars = truePars_fig3)

fig3data_nozero = fig3data
fig3data_nozero[fig3data_nozero==0]=10^-5

out_fig3data = solveLV(times = seq(1,17,.1), initState = state_fig3, pars = truePars_mat_fig3, equations = Equations)

### time, initState and truePars 
################################


#########################
### slope estimates

smoother_out = Smoother(
				dataset = fig3data_nozero,
				draw = TRUE,
				data1_spline2 = 2, 
				smooth_type = 1,
				df = c(14,10),	
				log_spline = TRUE
				)
plot(1,1,col='white')
legend(1,1,legend = c('data','splines','slopes'),pch = c(20,NA,NA),lty = c(0,1,1),col = c('grey','green','lightgreen'),bty = "n")

( est1 = LV_pars_finder(
				smooth_out = smoother_out,
				alg1_lm2 = 2
				) )
cbind(truePars_fig3,est1,absDif = abs(truePars_fig3-est1))
est1_mat = Format_pars(truePars = est1)
splineStart = smoother_out$splines_est[1,2:3];names(splineStart) = c('Individuals_Prey', 'Individuals_Predator')
out_est1 = solveLV(times = seq(1,17,.1), pars = est1_mat, initState = splineStart, equations = Equations)

par(mfrow=c(1,1),mar=c(5, 6, 4, 2) + 0.1)
plot(fig3data[,1],fig3data[,2],ylim=c(0,60),xlab='Days',ylab="Individuals",cex.lab=2,cex.axis=2)
points(fig3data[,1],fig3data[,3],col='blue',pch=0,ylim=c(0,200),xlab='Days',ylab="Individuals",cex.lab=2,cex.axis=2)
points(out_fig3data[,1],out_fig3data[,2],type='l',lwd=3)
points(out_fig3data[,1],out_fig3data[,3],type='l',col='blue',lty=3,lwd=3)
legend(7,60,legend = c("P. caudatum (prey)","D. nasutum (predator)","Slope est","MAR est"),lty = c(1,3,1,1),pch=c(1,0,NA,NA),col = c('black','blue','cyan','green'),bty = "n",cex=1.5,lwd=c(1,1,3,3)) 

for (v in 2:3)
	{
	points(smoother_out$splines_est[,1],smoother_out$splines_est[,v],type='l',col='green')
	points(out_est1[,1],out_est1[,v],type='l',col='cyan',lwd=3)
	}

### slope estimates
#########################

#######################
### MAR estimates



MARest_orig =  fastMAR(
			data_vars = fig3data_nozero[,2:3],
			log_transform = FALSE,
			demeaned = FALSE,
			abstol_= 0.01,
			maxit_ = 100,
			includeNoise = FALSE,
			estInit = FALSE,
			tSteps_ = NULL,
			nsim_ = 100
			) 
MARest_orig$MARSSci_output

MARest_log =  fastMAR(
			data_vars = fig3data_nozero[,2:3],
			log_transform = TRUE,
			demeaned = FALSE,
			abstol_= 0.01,
			maxit_ = 100,
			includeNoise = FALSE,
			estInit = FALSE,
			tSteps_ = NULL,
			nsim_ = 100
			) 
MARest_log$MARSSci_output

par(mfrow=c(1,1),mar=c(5, 6, 4, 2) + 0.1)
plot(fig3data[,1],fig3data[,2],ylim=c(0,60),xlab='Days',ylab="Individuals",cex.lab=2,cex.axis=2)
points(fig3data[,1],fig3data[,3],col='blue',pch=0,ylim=c(0,200),xlab='Days',ylab="Individuals",cex.lab=2,cex.axis=2)
points(out_fig3data[,1],out_fig3data[,2],type='l',lwd=3)
points(out_fig3data[,1],out_fig3data[,3],type='l',col='blue',lty=3,lwd=3)
legend(7,60,legend = c("P. caudatum (prey)","D. nasutum (predator)","Slope est","MAR est"),lty = c(1,3,1,1),pch=c(1,0,NA,NA),col = c('black','blue','cyan','green'),bty = "n",cex=1.5,lwd=c(1,1,3,3)) 

for (v in 2:3)
	{
	points(MARest_orig$MAR_TS_Est[,v-1] ,type='l',col='green',lwd=3)
	points(MARest_log$MAR_TS_Est[,v-1] ,type='l',col='green',lty=2,lwd=3)
	points(meanMaker(MARest_orig)[v-1,],type='l',col='purple',lty=1)
	points(meanMaker(MARest_log)[v-1,],type='l',col='purple',lty=2)
	}

### MAR estimates
#######################


#########################
### MAR estimates with smoothing

# smoothing the data for MAR
smoother_for_MAR_out = Smoother_for_MAR(	
						dataset = fig3data_nozero,
						df = c(14,10),
						dataWeights = NULL,
						splineMethod = "fmm",
						draw = TRUE,
						log_spline = TRUE
						)
smoother_for_MAR_out
TS_smooth_data = smooth_TS = smoother_for_MAR_out$splines_est

# MAR estimates without data log transformation, with spline smooth data  
MARest1_smoothData =  fastMAR(
			data_vars = TS_smooth_data[,2:3],
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
			data_vars = TS_smooth_data[,2:3],
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

### MAR estimates with smoothing
#########################

# tiff("FigS6c.tiff", height = 15, width = 7.5, units = 'cm',compression = "lzw", res = 300)
par(mfrow=c(2,1))		# ,mar=c(5, 6, 4, 2) + 0.1
main_lab = c( expression(paste(italic('P. caudatum'), " (Prey)")),expression(paste(italic('D. nasutum'), " (Predator)")) )
for (v in 2:3)
	{
	plot(fig3data[,1], fig3data[,v], ylim=c(0,60), col=colorPallet[1], main='', xlab='Days', ylab="", cex.main=.8, cex.lab=1.5, cex.axis=1.5)
	points(out_fig3data[,1],out_fig3data[,2],type='l',col=colorPallet[2],lwd=3)
	points(out_est1[,1],out_est1[,v],type='l',col=colorPallet[3],lwd=3)
	points(MARest_orig$MAR_TS_Est[,v-1] ,type='l',col=colorPallet[4],lwd=3)
	points(MARest_log$MAR_TS_Est[,v-1] ,type='l',col=colorPallet[5],lwd=3)
	points(MARest1_smoothData$MAR_TS_Est[,v-1],type='l',col='orange',lwd=3)
	points(MARest2_smoothData$MAR_TS_Est[,v-1],type='l',col='khaki',lwd=3)
	points(out_est1[,1],out_est1[,v],type='l',col=colorPallet[3],lwd=3)
	mtext(main_lab[v-1],side=3,line=1,cex=1.5)
	mtext("Individuals",side=2,line=2.5,cex=1.5)
	#if (v==2) { legend(7,65, c('data','Muhlbauer est','Slope est','MAR','MAR logTrans'), col=col=colorPallet,pch=c(1,NA,NA,NA,NA), lwd=c(NA,3,3,3,3), lty=c(NA,1,1,1,1),bty = "n",cex=1) }
	}
# dev.off()


Muhlbauer = out_fig3data[,1:3]
Muhlbauer = Muhlbauer[Muhlbauer[,1]%in%fig3data[,1],]
sum( (fig3data[,2:3]-Muhlbauer[,2:3])^2 )
 
sum( (fig3data[,2:3]-out_est1[out_est1[,1]%in%t_fig3,2:3])^2 )
sum( (fig3data[,2:3]-MARest_orig$MAR_TS_Est)^2 )
sum( (fig3data[,2:3]-MARest_log$MAR_TS_Est)^2 )
sum( (fig3data[,2:3]-MARest1_smoothData$MAR_TS_Est)^2 )
sum( (fig3data[,2:3]-MARest2_smoothData$MAR_TS_Est)^2 )

##### fig3
######################################################################################################################################



######################################################################################################################################
##### fig4

################################
### time, initState and truePars  


t_fig4 = seq(1959,1992,1)
fig4data
state_fig4 = c(24.312,699.9,0.326)
names(state_fig4) = c('x1','x2','x3')
truePars_fig4 = c(		
			a1 = 0.01,
			b11 = -0.003,
			b12 = 0.00004,
			b13 = 0,
	
			a2 = 2.021,
			b21 = -0.088,
			b22 = 0,
			b23 = 0.002,
	
			a3 = 0.238,
			b31 = 0,
			b32 = -0.0002,
			b33 = -0.139
		      )
truePars_mat = truePars_mat_fig4 = Format_pars(truePars = truePars_fig4)

out_fig4data = solveLV(times = t_fig4, initState = state_fig4, pars = truePars_mat_fig4, equations = Equations)

### time, initState and truePars 
################################


#########################
### slope estimates

smoother_out = Smoother(
				dataset = fig4data,
				draw = TRUE,
				data1_spline2 = 2, 
				smooth_type = 1,
				df = c(28,24,28),	
				log_spline = TRUE
				)
plot(1,1,col='white')
legend(1,1,legend = c('data','splines','slopes'),pch = c(20,NA,NA),lty = c(0,1,1),col = c('grey','green','lightgreen'),bty = "n")

( est1 = LV_pars_finder(
				smooth_out = smoother_out,
				alg1_lm2 = 2
				) )
cbind(truePars_fig4,est1,absDif = abs(truePars_fig4-est1))
est1_mat = Format_pars(truePars = est1)
splineStart = smoother_out$splines_est[1,2:4];names(splineStart) = c("x1", "x2","x3")
out_est1 = solveLV(times = t_fig4, pars = est1_mat, initState = splineStart, equations = Equations)

par(mfrow=c(2,2))		#,mar=c(5, 6, 4, 2) + 0.1   \n
plot(fig4data[,1],fig4data[,2],col='red',xlab='Year',ylab='Wolves (Individuals)',cex.lab=1.5,cex.axis=1.5)
points(out_fig4data[,1],out_fig4data[,2],type='l',col='red',lwd=2)
points(out_est1[,1],out_est1[,2],col='cyan',type='l',lwd=3)
points(smoother_out$splines_est[,1],smoother_out$splines_est[,2],type='l',col='green')

plot(fig4data[,1],fig4data[,3],pch=2,xlab='Year',ylab='Moose (Individuals)',cex.lab=1.5,cex.axis=1.5)
points(out_fig4data[,1],out_fig4data[,3],type='l',lty=3,lwd=3)
points(out_est1[,1],out_est1[,3],type='l',col='cyan',lty=3,lwd=3)
points(smoother_out$splines_est[,1],smoother_out$splines_est[,3],type='l',col='green')

plot(fig4data[,1],fig4data[,4],col='blue',pch=0,xlab='Year',ylab='Fir tree rings (tree growth)',cex.lab=1.5,cex.axis=1.5)
points(out_fig4data[,1],out_fig4data[,4],type='l',col='blue',lty=3,lwd=3)
points(out_est1[,1],out_est1[,4],type='l',col='cyan',lty=3,lwd=3)
points(smoother_out$splines_est[,1],smoother_out$splines_est[,4],type='l',col='green')

### slope estimates
#########################

#######################
### MAR estimates

MARest_orig =  fastMAR(
			data_vars = fig4data[,2:4],
			log_transform = FALSE,
			demeaned = TRUE,
			abstol_= 0.01,
			maxit_ = 500,
			includeNoise = FALSE,
			estInit = FALSE,
			tSteps_ = NULL,
			nsim_ = 100
			) 
MARest_orig$MARSSci_output

MARest_log =  fastMAR(
			data_vars = fig4data[,2:4],
			log_transform = TRUE,
			demeaned = FALSE,
			abstol_= 0.01,
			maxit_ = 5000,
			includeNoise = FALSE,
			estInit = FALSE,
			tSteps_ = NULL,
			nsim_ = 100
			) 
MARest_log$MARSSci_output

par(mfrow=c(2,2))		#,mar=c(5, 6, 4, 2) + 0.1   \n
plot(fig4data[,1],fig4data[,2],col='red',xlab='Year',ylab='Wolves (Individuals)',cex.lab=1.5,cex.axis=1.5)
points(out_fig4data[,1],out_fig4data[,2],type='l',col='red',lwd=2)
points(fig4data[,1],MARest_orig$MAR_TS_Est[,1],col='green',type='l',lwd=1)
points(fig4data[,1],MARest_log$MAR_TS_Est[,1],col='green',type='l',lty=2,lwd=1)
	points(fig4data[,1],meanMaker(MARest_orig)[1,],type='l',col='purple',lty=1)
	points(fig4data[,1],meanMaker(MARest_log)[1,],type='l',col='purple',lty=2)

plot(fig4data[,1],fig4data[,3],col='black',xlab='Year',ylab='Wolves (Individuals)',cex.lab=1.5,cex.axis=1.5)
points(out_fig4data[,1],out_fig4data[,3],type='l',col='black',lwd=2)
points(fig4data[,1],MARest_orig$MAR_TS_Est[,2],col='green',type='l',lwd=1)
points(fig4data[,1],MARest_log$MAR_TS_Est[,2],col='green',type='l',lty=2,lwd=1)
	points(fig4data[,1],meanMaker(MARest_orig)[2,],type='l',col='purple',lty=1)
	points(fig4data[,1],meanMaker(MARest_log)[2,],type='l',col='purple',lty=2)

plot(fig4data[,1],fig4data[,4],col='blue',xlab='Year',ylab='Wolves (Individuals)',cex.lab=1.5,cex.axis=1.5)
points(out_fig4data[,1],out_fig4data[,4],type='l',col='blue',lwd=2)
points(fig4data[,1],MARest_orig$MAR_TS_Est[,3],col='green',type='l',lwd=1)
points(fig4data[,1],MARest_log$MAR_TS_Est[,3],col='green',type='l',lty=2,lwd=1)
	points(fig4data[,1],meanMaker(MARest_orig)[3,],type='l',col='purple',lty=1)
	points(fig4data[,1],meanMaker(MARest_log)[3,],type='l',col='purple',lty=2)

### MAR estimates
#######################


#########################
### MAR estimates with smoothing

# smoothing the data for MAR
smoother_for_MAR_out = Smoother_for_MAR(	
						dataset = fig4data,
						df = c(8,8,8),
						dataWeights = NULL,
						splineMethod = "fmm",
						draw = TRUE,
						log_spline = TRUE
						)
smoother_for_MAR_out
TS_smooth_data = smooth_TS = smoother_for_MAR_out$splines_est

# MAR estimates without data log transformation, with spline smooth data  
MARest1_smoothData =  fastMAR(
			data_vars = TS_smooth_data[,2:4],
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
			data_vars = TS_smooth_data[,2:4],
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

### MAR estimates with smoothing
#########################

# tiff("FigS6d1.tiff", height = 15, width = 7.5, units = 'cm',compression = "lzw", res = 300)
par(mfrow=c(2,1))		#,mar=c(5, 6, 4, 2) + 0.1   

plot(fig4data[,1], fig4data[,2], col=colorPallet[1], main='', xlab='Year', ylab='Individuals', cex.main=.8, cex.lab=1.5, cex.axis=1.5)
points(out_fig4data[,1],out_fig4data[,2],type='l',col=colorPallet[2],lwd=3)
points(out_est1[,1],out_est1[,2],col=colorPallet[3],type='l',lwd=3)
points(fig4data[,1],MARest_orig$MAR_TS_Est[,1],col=colorPallet[4],type='l',lwd=3)
points(fig4data[,1],MARest_log$MAR_TS_Est[,1],col=colorPallet[5],type='l',lwd=3)
points(1959:1992,MARest1_smoothData$MAR_TS_Est[,1],type='l',col='orange',lwd=3)
points(1959:1992,MARest2_smoothData$MAR_TS_Est[,1],type='l',col='khaki',lwd=3)
points(out_est1[,1],out_est1[,2],col=colorPallet[3],type='l',lwd=3)
mtext('Wolves',side=3,line=1,cex=1.5)

plot(fig4data[,1], fig4data[,3], col=colorPallet[1], main='', xlab='Year', ylab='Individuals', cex.main=.8, cex.lab=1.5, cex.axis=1.5)
points(out_fig4data[,1],out_fig4data[,3],type='l',col=colorPallet[2],lwd=3)
points(out_est1[,1],out_est1[,3],type='l',col=colorPallet[3],lwd=3)
points(fig4data[,1],MARest_orig$MAR_TS_Est[,2],col=colorPallet[4],type='l',lwd=3)
points(fig4data[,1],MARest_log$MAR_TS_Est[,2],col=colorPallet[5],type='l',lwd=3)
points(1959:1992,MARest1_smoothData$MAR_TS_Est[,2],type='l',col='orange',lwd=3)
points(1959:1992,MARest2_smoothData$MAR_TS_Est[,2],type='l',col='khaki',lwd=3)
points(out_est1[,1],out_est1[,3],type='l',col=colorPallet[3],lwd=3)
mtext('Moose',side=3,line=1,cex=1.5)
# dev.off()

# tiff("FigS6d2.tiff", height = 15, width = 7.5, units = 'cm',compression = "lzw", res = 300)
par(mfrow=c(2,1))		#,mar=c(5, 6, 4, 2) + 0.1   
plot(fig4data[,1],fig4data[,4],col=colorPallet[1],main='', xlab='Year', ylab='tree growth', cex.main=.8, cex.lab=1.5, cex.axis=1.5)
points(out_fig4data[,1],out_fig4data[,4],type='l',col=colorPallet[2],lwd=3)
points(out_est1[,1],out_est1[,4],type='l',col=colorPallet[3],lwd=3)
points(fig4data[,1],MARest_orig$MAR_TS_Est[,3],col=colorPallet[4],type='l',lwd=3)
points(fig4data[,1],MARest_log$MAR_TS_Est[,3],col=colorPallet[5],type='l',lwd=3)
points(1959:1992,MARest1_smoothData$MAR_TS_Est[,3],type='l',col='orange',lwd=3)
points(1959:1992,MARest2_smoothData$MAR_TS_Est[,3],type='l',col='khaki',lwd=3)
points(out_est1[,1],out_est1[,4],type='l',col=colorPallet[3],lwd=3)
mtext('Fir tree rings',side=3,line=1,cex=1.5)

plot(1,1,col='white',xlab='',ylab='',axes = FALSE)
legend(.6,1.3,legend = c('data','Mühlbauer est','ALVI-LR method','MAR no transform','MAR log transform','MAR no transform with smoothing','MAR log transform with smoothing'),lty = c(NA,1,1,1,1,1,1),lwd=c(NA,3,3,3,3,3,3),pch = c(1,NA,NA,NA,NA,NA,NA),col = c(colorPallet,'orange','khaki'),bty = "n",cex=.5) # ,
# dev.off()

Muhlbauer = out_fig4data[,1:4]
sum( (fig4data[,2:4]-Muhlbauer[,2:4])^2 )

sum( (fig4data[,2:4]-out_est1[,2:4])^2 )
sum( (fig4data[,2:4]-MARest_orig$MAR_TS_Est)^2 )
sum( (fig4data[,2:4]-MARest_log$MAR_TS_Est)^2 )
sum( (fig4data[,2:4]-MARest1_smoothData$MAR_TS_Est)^2 )
sum( (fig4data[,2:4]-MARest2_smoothData$MAR_TS_Est)^2 )

##### fig4
######################################################################################################################################



######################################################################################################################################
##### fig5

################################
### time, initState and truePars  

t_fig5 = seq(1,60,1)
fig5data
state_fig5 = c(64.029,5.23)
names(state_fig5) = c('x1','x2')
truePars_fig5 = c(		
			a1 = 0.187,
			b11 = 0,
			b12 = -0.028,

			a2 = -0.377,
			b21 = 0.0012,
			b22 = -0.024
      		)
truePars_mat_fig5 = Format_pars(truePars = truePars_fig5)

out_fig5data = solveLV(times = seq(1,60,.1), initState = state_fig5, pars = truePars_mat_fig5, equations = Equations)

### time, initState and truePars 
################################


#########################
### slope estimates

smoother_out = Smoother(
				dataset = fig5data,
				draw = TRUE,
				data1_spline2 = 2, 
				smooth_type = 1,
				df = c(15,20),	
				log_spline = TRUE
				)
plot(1,1,col='white')
legend(1,1,legend = c('data','splines','slopes'),pch = c(20,NA,NA),lty = c(0,1,1),col = c('grey','green','lightgreen'),bty = "n")

( est1 = LV_pars_finder(
				smooth_out = smoother_out,
				alg1_lm2 = 2
				) )
cbind(truePars_fig5,est1,absDif = abs(truePars_fig5-est1))
est1_mat = Format_pars(truePars = est1)
splineStart = smoother_out$splines_est[1,2:3];names(splineStart) = c("x1", "x2")
out_est1 = solveLV(times = seq(1,60,.1), pars = est1_mat, initState = splineStart, equations = Equations)

par(mfrow=c(2,1))		#,mar=c(5, 6, 4, 2) + 0.1   \n
plot(fig5data[,1],fig5data[,2],ylim=c(0,1300),main = "Eotetranychus sexmaculatus (Prey)",xlab='Week',ylab="Individuals",cex.lab=1.5,cex.axis=1.5)
points(out_fig5data[,1],out_fig5data[,2],type='l',lwd=3)
points(smoother_out$splines_est[,1],smoother_out$splines_est[,2],type='l',col='green')
points(out_est1[,1],out_est1[,2],type='l',col='cyan',lwd=3)
legend(47,1400,legend = c('data','Slope est','spline'), pch = c(1,NA,NA), lty = c(1,1,1),col = c('black','cyan','green'), bty = "n") 

plot(fig5data[,1],fig5data[,3],col='blue',ylim=c(0,30),main = "Eotetranychus sexmaculatus (Prey)",xlab='Week',ylab="Individuals",cex.lab=1.5,cex.axis=1.5)
points(out_fig5data[,1],out_fig5data[,3],col='blue',type='l',lwd=3)
points(smoother_out$splines_est[,1],smoother_out$splines_est[,3],type='l',col='green')
points(out_est1[,1],out_est1[,3],type='l',col='cyan',lwd=3)
legend(47,30,legend = c('data','Slope est','spline'), pch = c(1,NA,NA), lty = c(1,1,1),col = c('blue','cyan','green'), bty = "n") 

### slope estimates
#########################

#######################
### MAR estimates

MARest_orig =  fastMAR(
			data_vars = fig5data[,2:3],
			log_transform = FALSE,
			demeaned = TRUE,
			abstol_= 0.01,
			maxit_ = 100,
			includeNoise = FALSE,
			estInit = FALSE,
			tSteps_ = NULL,
			nsim_ = 100
			) 
MARest_orig$MARSSci_output

MARest_log =  fastMAR(
			data_vars = fig5data[,2:3],
			log_transform = TRUE,
			demeaned = TRUE,
			abstol_= 0.01,
			maxit_ = 600,
			includeNoise = FALSE,
			estInit = FALSE,
			tSteps_ = NULL,
			nsim_ = 100
			) 
MARest_log$MARSSci_output

par(mfrow=c(2,1))		#,mar=c(5, 6, 4, 2) + 0.1   \n
plot(fig5data[,1],fig5data[,2],ylim=c(0,1300),main = "Eotetranychus sexmaculatus (Prey)",xlab='Week',ylab="Individuals",cex.lab=1.5,cex.axis=1.5)
points(out_fig5data[,1],out_fig5data[,2],type='l',lwd=3)
points(MARest_orig$MAR_TS_Est[,1] ,type='l',col='green',lwd=3)
points(MARest_log$MAR_TS_Est[,1] ,type='l',col='green',lty=2,lwd=3)
legend(47,1400,legend = c('data','Slope est','spline'), pch = c(1,NA,NA), lty = c(1,1,1),col = c('black','cyan','green'), bty = "n") 
	points(meanMaker(MARest_orig)[1,],type='l',col='purple',lty=1)
	points(meanMaker(MARest_log)[1,],type='l',col='purple',lty=2)

plot(fig5data[,1],fig5data[,3],col='blue',ylim=c(0,30),main = "Eotetranychus sexmaculatus (Prey)",xlab='Week',ylab="Individuals",cex.lab=1.5,cex.axis=1.5)
points(out_fig5data[,1],out_fig5data[,3],col='blue',type='l',lwd=3)
points(MARest_orig$MAR_TS_Est[,2] ,type='l',col='green',lwd=3)
points(MARest_log$MAR_TS_Est[,2] ,type='l',col='green',lty=2,lwd=3)
legend(47,30,legend = c('data','Slope est','spline'), pch = c(1,NA,NA), lty = c(1,1,1),col = c('blue','cyan','green'), bty = "n") 
	points(meanMaker(MARest_orig)[2,],type='l',col='purple',lty=1)
	points(meanMaker(MARest_log)[2,],type='l',col='purple',lty=2)

### MAR estimates
#######################


#########################
### MAR estimates with smoothing

# smoothing the data for MAR
smoother_for_MAR_out = Smoother_for_MAR(	
						dataset = fig5data,
						df = c(15,20),
						dataWeights = NULL,
						splineMethod = "fmm",
						draw = TRUE,
						log_spline = TRUE
						)
smoother_for_MAR_out
TS_smooth_data = smooth_TS = smoother_for_MAR_out$splines_est

# MAR estimates without data log transformation, with spline smooth data  
MARest1_smoothData =  fastMAR(
			data_vars = TS_smooth_data[,2:3],
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
			data_vars = TS_smooth_data[,2:3],
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

### MAR estimates with smoothing
#########################


# setwd("C:\\Users\\doliv\\OneDrive\\2_2019_America\\2020\\20200123_MAR\\Camparison_LV_MAR\\8_methods_of_ecology_and_evolution_2sub\\paper_figures")     # P51
# tiff("FigS6e.tiff", height = 15, width = 7.5, units = 'cm',compression = "lzw", res = 300)
par(mfrow=c(2,1))		#,mar=c(5, 6, 4, 2) + 0.1   \n
plot(fig5data[,1],fig5data[,2],ylim=c(0,1300),col=colorPallet[1],main = expression(paste(italic('E. sexmaculatus'))),xlab='Week',ylab="",cex.main=1.5,cex.lab=1.5,cex.axis=1.5)
points(out_fig5data[,1],out_fig5data[,2],type='l',col=colorPallet[2],lwd=3)
points(out_est1[,1],out_est1[,2],type='l',col=colorPallet[3],lwd=3)
points(MARest_orig$MAR_TS_Est[,1] ,type='l',col=colorPallet[4],lwd=3)
points(MARest_log$MAR_TS_Est[,1] ,type='l',col=colorPallet[5],lwd=3)
points(MARest1_smoothData$MAR_TS_Est[,1],type='l',col='orange',lwd=3)
points(MARest2_smoothData$MAR_TS_Est[,1],type='l',col='khaki',lwd=3)
points(out_est1[,1],out_est1[,2],type='l',col=colorPallet[3],lwd=3)
mtext("(Prey)",side=3,line=.4,cex=1.5)
mtext("Individuals",side=2,line=2.5,cex=1.5)
#legend(0,1400,legend = c('data','Muhlbauer est','Slope est','MAR','MAR logTrans'), pch = c(1,NA,NA,NA,NA), lty = c(NA,1,1,1,1),col=colorPallet, bty = "n", lwd=c(NA,3,3,3,3),cex=.7) 

plot(fig5data[,1], fig5data[,3],col=colorPallet[1],ylim=c(0,30),main = expression(paste(italic("T. occidentalis"))), xlab='Week', ylab="", cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
points(out_fig5data[,1],out_fig5data[,3],type='l',col=colorPallet[2],lwd=3)
points(out_est1[,1],out_est1[,3],type='l',col=colorPallet[3],lwd=3)
points(MARest_orig$MAR_TS_Est[,2] ,type='l',col=colorPallet[4],lwd=3)
points(MARest_log$MAR_TS_Est[,2] ,type='l',col=colorPallet[5],lwd=3)
points(MARest1_smoothData$MAR_TS_Est[,2],type='l',col='orange',lwd=3)
points(MARest2_smoothData$MAR_TS_Est[,2],type='l',col='khaki',lwd=3)
points(out_est1[,1],out_est1[,3],type='l',col=colorPallet[3],lwd=3)
mtext("(Predator)",side=3,line=.4,cex=1.5)
mtext("Individuals",side=2,line=2.5,cex=1.5)
#legend(0,30,legend = c('data','Muhlbauer est','Slope est','MAR','MAR logTrans'), pch = c(1,NA,NA,NA,NA), lty = c(NA,1,1,1,1),col=colorPallet, bty = "n", lwd=c(NA,3,3,3,3),cex=.7) 
# dev.off()

Muhlbauer = out_fig5data[,1:3]
Muhlbauer = Muhlbauer[Muhlbauer[,1]%in%fig5data[,1],]
sum( (fig5data[,2:3]-Muhlbauer[,2:3])^2 )


sum( (fig5data[,2:3]-out_est1[out_est1[,1]%in%fig5data[,1],2:3])^2 )
sum( (fig5data[,2:3]-MARest_orig$MAR_TS_Est)^2 )
sum( (fig5data[,2:3]-MARest_log$MAR_TS_Est)^2 )
sum( (fig5data[,2:3]-MARest1_smoothData$MAR_TS_Est)^2 )
sum( (fig5data[,2:3]-MARest2_smoothData$MAR_TS_Est)^2 )

##### fig5
######################################################################################################################################











######################################################################################################################################
##### fig4 - minus 1981, 1982

################################
### time, initState and truePars  

t_fig4 = seq(1959,1992,1)
fig4data
t_fig4_1 = seq(1959,1980,1)
fig4data_1 = fig4data[1:22,]
t_fig4_2 = seq(1983,1992,1)
fig4data_2 = fig4data[25:34,]

state_fig4 = c(24.312,699.9,0.326)
names(state_fig4) = c('x1','x2','x3')
truePars_fig4 = c(		
			a1 = 0.01,
			b11 = -0.003,
			b12 = 0.00004,
			b13 = 0,
	
			a2 = 2.021,
			b21 = -0.088,
			b22 = 0,
			b23 = 0.002,
	
			a3 = 0.238,
			b31 = 0,
			b32 = -0.0002,
			b33 = -0.139
		      )
truePars_mat = truePars_mat_fig4 = Format_pars(truePars = truePars_fig4)

out_fig4data = solveLV(times = t_fig4, initState = state_fig4, pars = truePars_mat_fig4, equations = Equations)

### time, initState and truePars 
################################


#########################
### slope estimates

smoother_out = Smoother(
				dataset = fig4data,
				draw = TRUE,
				data1_spline2 = 2, 
				smooth_type = 1,
				df = c(8,8,8),	
				log_spline = TRUE
				)
plot(1,1,col='white')
legend(1,1,legend = c('data','splines','slopes'),pch = c(20,NA,NA),lty = c(0,1,1),col = c('grey','green','lightgreen'),bty = "n")

AlgPointFind(
		smoother_out_APF = smoother_out,
		dif_equations = Equations,
		matrixEq = TRUE
		)
#      "100 % complete || DF = 10 8 10  ||  5 26 28 31 531508.248  ||  5 26 28 31 207607.651"
# log  "100 % complete || DF = 10 8 10  ||  4 27 28 29 477711.377  ||  4 27 28 29 239168.672"    # BETTER SSE

( est1 = LV_pars_finder(
				smooth_out = smoother_out,
				alg1_lm2 = 1, 
				data_sample_alg = c(15, 20, 21, 24)		#c(3, 7, 10, 14, 22)		
				) )
cbind(truePars_fig4,est1,absDif = abs(truePars_fig4-est1))
est1_mat = Format_pars(truePars = est1)
splineStart = smoother_out$splines_est[1,2:4];names(splineStart) = c("x1", "x2","x3")
out_est1 = solveLV(times = t_fig4, pars = est1_mat, initState = splineStart, equations = Equations)

par(mfrow=c(2,2))		#,mar=c(5, 6, 4, 2) + 0.1   \n
plot(fig4data[,1],fig4data[,2],col='red',xlab='Year',ylab='Wolves (Individuals)',cex.lab=1.5,cex.axis=1.5)
points(out_fig4data[,1],out_fig4data[,2],type='l',col='red',lwd=2)
points(out_est1[,1],out_est1[,2],col='cyan',type='l',lwd=3)
points(smoother_out$splines_est[,1],smoother_out$splines_est[,2],type='l',col='green')

plot(fig4data[,1],fig4data[,3],pch=2,xlab='Year',ylab='Moose (Individuals)',cex.lab=1.5,cex.axis=1.5)
points(out_fig4data[,1],out_fig4data[,3],type='l',lty=3,lwd=3)
points(out_est1[,1],out_est1[,3],type='l',col='cyan',lty=3,lwd=3)
points(smoother_out$splines_est[,1],smoother_out$splines_est[,3],type='l',col='green')

plot(fig4data[,1],fig4data[,4],col='blue',pch=0,xlab='Year',ylab='Fir tree rings (tree growth)',cex.lab=1.5,cex.axis=1.5)
points(out_fig4data[,1],out_fig4data[,4],type='l',col='blue',lty=3,lwd=3)
points(out_est1[,1],out_est1[,4],type='l',col='cyan',lty=3,lwd=3)
points(smoother_out$splines_est[,1],smoother_out$splines_est[,4],type='l',col='green')

### slope estimates
#########################


#########################
### slope estimates fig4data_1

smoother_out = Smoother(
				dataset = fig4data_1,
				draw = TRUE,
				data1_spline2 = 2, 
				smooth_type = 1,
				df = c(6,8,8),	
				log_spline = TRUE
				)
plot(1,1,col='white')
legend(1,1,legend = c('data','splines','slopes'),pch = c(20,NA,NA),lty = c(0,1,1),col = c('grey','green','lightgreen'),bty = "n")

AlgPointFind(
		smoother_out_APF = smoother_out,
		dif_equations = Equations,
		matrixEq = TRUE
		)
# log  "100 % complete || DF = 8 8 8  ||  3 7 9 19 182655.037  ||  3 7 9 19 151585.223"    # BETTER SSE
# log  "100 % complete || DF = 6 8 8  ||  1 11 13 21 167470.462  ||  1 11 13 21 140513.878"

( est1 = LV_pars_finder(
				smooth_out = smoother_out,
				alg1_lm2 = 1, 
				data_sample_alg = c(1, 11, 13, 21)		#c(3, 7, 10, 14, 22)		
				) )
cbind(truePars_fig4,est1,absDif = abs(truePars_fig4-est1))
est1_mat = Format_pars(truePars = est1)
splineStart = smoother_out$splines_est[1,2:4];names(splineStart) = c("x1", "x2","x3")
out_est1_1 = solveLV(times = t_fig4_1, pars = est1_mat, initState = splineStart, equations = Equations)

par(mfrow=c(2,2))		#,mar=c(5, 6, 4, 2) + 0.1   \n
plot(fig4data[,1],fig4data[,2],col='red',xlab='Year',ylab='Wolves (Individuals)',cex.lab=1.5,cex.axis=1.5)
points(out_fig4data[,1],out_fig4data[,2],type='l',col='red',lwd=2)
points(out_est1_1[,1],out_est1_1[,2],col='cyan',type='l',lwd=3)
points(smoother_out$splines_est[,1],smoother_out$splines_est[,2],type='l',col='green')

plot(fig4data[,1],fig4data[,3],pch=2,xlab='Year',ylab='Moose (Individuals)',cex.lab=1.5,cex.axis=1.5)
points(out_fig4data[,1],out_fig4data[,3],type='l',lty=3,lwd=3)
points(out_est1_1[,1],out_est1_1[,3],type='l',col='cyan',lty=3,lwd=3)
points(smoother_out$splines_est[,1],smoother_out$splines_est[,3],type='l',col='green')

plot(fig4data[,1],fig4data[,4],col='blue',pch=0,xlab='Year',ylab='Fir tree rings (tree growth)',cex.lab=1.5,cex.axis=1.5)
points(out_fig4data[,1],out_fig4data[,4],type='l',col='blue',lty=3,lwd=3)
points(out_est1_1[,1],out_est1_1[,4],type='l',col='cyan',lty=3,lwd=3)
points(smoother_out$splines_est[,1],smoother_out$splines_est[,4],type='l',col='green')

### slope estimates fig4data_1
#########################


#########################
### slope estimates fig4data_2

smoother_out = Smoother(
				dataset = fig4data_2,
				draw = TRUE,
				data1_spline2 = 2, 
				smooth_type = 1,
				df = c(8,8,8),	
				log_spline = TRUE
				)
plot(1,1,col='white')
legend(1,1,legend = c('data','splines','slopes'),pch = c(20,NA,NA),lty = c(0,1,1),col = c('grey','green','lightgreen'),bty = "n")

AlgPointFind(
		smoother_out_APF = smoother_out,
		dif_equations = Equations,
		matrixEq = TRUE
		)
# log   "100 % complete || DF = 6 6 6  ||  2 3 6 10 233394.599  ||  2 3 4 10 120040.643"   # BETTER SSE
# log   "100 % complete || DF = 8 8 8  ||  2 3 4  9 257219.425  ||  2 3 4  9 210523.642"
# log	  "100 % complete || DF = 10 10 10  ||  1 2 4 6 283651.391  ||  1 2 4 6 283651.391"
# log	  "100 % complete || DF = 6 6 5  ||  1 2 4  9 227689.001  ||  1 2 4  9 99427.541"

( est1 = LV_pars_finder(
				smooth_out = smoother_out,
				alg1_lm2 = 1, 
				data_sample_alg = c(2, 3, 4, 9)		#c(3, 7, 10, 14, 22)		
				) )
cbind(truePars_fig4,est1,absDif = abs(truePars_fig4-est1))
est1_mat = Format_pars(truePars = est1)
splineStart = smoother_out$splines_est[1,2:4];names(splineStart) = c("x1", "x2","x3")
out_est1_2 = solveLV(times = t_fig4_2, pars = est1_mat, initState = splineStart, equations = Equations)

par(mfrow=c(2,2))		#,mar=c(5, 6, 4, 2) + 0.1   \n
plot(fig4data[,1],fig4data[,2],col='red',xlab='Year',ylab='Wolves (Individuals)',cex.lab=1.5,cex.axis=1.5)
points(out_fig4data[,1],out_fig4data[,2],type='l',col='red',lwd=2)
points(out_est1_2[,1],out_est1_2[,2],col='cyan',type='l',lwd=3)
points(smoother_out$splines_est[,1],smoother_out$splines_est[,2],type='l',col='green')

plot(fig4data[,1],fig4data[,3],pch=2,xlab='Year',ylab='Moose (Individuals)',cex.lab=1.5,cex.axis=1.5)
points(out_fig4data[,1],out_fig4data[,3],type='l',lty=3,lwd=3)
points(out_est1_2[,1],out_est1_2[,3],type='l',col='cyan',lty=3,lwd=3)
points(smoother_out$splines_est[,1],smoother_out$splines_est[,3],type='l',col='green')

plot(fig4data[,1],fig4data[,4],col='blue',pch=0,xlab='Year',ylab='Fir tree rings (tree growth)',cex.lab=1.5,cex.axis=1.5)
points(out_fig4data[,1],out_fig4data[,4],type='l',col='blue',lty=3,lwd=3)
points(out_est1_2[,1],out_est1_2[,4],type='l',col='cyan',lty=3,lwd=3)
points(smoother_out$splines_est[,1],smoother_out$splines_est[,4],type='l',col='green')

### slope estimates fig4data_2
#########################



#########################
### MAR estimates

MARest_orig =  fastMAR(
			data_vars = fig4data[,2:4],
			log_transform = FALSE,
			demeaned = FALSE,
			abstol_= 0.01,
			maxit_ = 500,
			includeNoise = FALSE,
			estInit = TRUE,
			tSteps_ = NULL,
			nsim_ = 100
			) 
MARest_orig$MARSSci_output

MARest_log =  fastMAR(
			data_vars = fig4data[,2:4],
			log_transform = TRUE,
			demeaned = FALSE,
			abstol_= 0.01,
			maxit_ = 5000,
			includeNoise = FALSE,
			estInit = FALSE,
			tSteps_ = NULL,
			nsim_ = 100
			) 
MARest_log$MARSSci_output

par(mfrow=c(2,2))		#,mar=c(5, 6, 4, 2) + 0.1   \n
plot(fig4data[,1],fig4data[,2],col='red',xlab='Year',ylab='Wolves (Individuals)',cex.lab=1.5,cex.axis=1.5)
points(out_fig4data[,1],out_fig4data[,2],type='l',col='red',lwd=2)
points(fig4data[,1],MARest_orig$MAR_TS_Est[,1],col='green',type='l',lwd=1)
points(fig4data[,1],MARest_log$MAR_TS_Est[,1],col='green',type='l',lty=2,lwd=1)
	points(fig4data[,1],meanMaker(MARest_orig)[1,],type='l',col='purple',lty=1)
	points(fig4data[,1],meanMaker(MARest_log)[1,],type='l',col='purple',lty=2)

plot(fig4data[,1],fig4data[,3],col='black',xlab='Year',ylab='Wolves (Individuals)',cex.lab=1.5,cex.axis=1.5)
points(out_fig4data[,1],out_fig4data[,3],type='l',col='black',lwd=2)
points(fig4data[,1],MARest_orig$MAR_TS_Est[,2],col='green',type='l',lwd=1)
points(fig4data[,1],MARest_log$MAR_TS_Est[,2],col='green',type='l',lty=2,lwd=1)
	points(fig4data[,1],meanMaker(MARest_orig)[2,],type='l',col='purple',lty=1)
	points(fig4data[,1],meanMaker(MARest_log)[2,],type='l',col='purple',lty=2)

plot(fig4data[,1],fig4data[,4],col='blue',xlab='Year',ylab='Wolves (Individuals)',cex.lab=1.5,cex.axis=1.5)
points(out_fig4data[,1],out_fig4data[,4],type='l',col='blue',lwd=2)
points(fig4data[,1],MARest_orig$MAR_TS_Est[,3],col='green',type='l',lwd=1)
points(fig4data[,1],MARest_log$MAR_TS_Est[,3],col='green',type='l',lty=2,lwd=1)
	points(fig4data[,1],meanMaker(MARest_orig)[3,],type='l',col='purple',lty=1)
	points(fig4data[,1],meanMaker(MARest_log)[3,],type='l',col='purple',lty=2)

### MAR estimates
#########################

# setwd("C:\\Users\\doliv\\OneDrive\\2_2019_America\\2020\\20200123_MAR\\Camparison_LV_MAR\\6_methods_of_ecology_and_evolution_new_images\\paper_figures")     # P51
# tiff("FigS7.tiff", height = 15, width = 15, units = 'cm',compression = "lzw", res = 300)
par(mfrow=c(2,2))		#,mar=c(5, 6, 4, 2) + 0.1   \n
plot(fig4data[,1], fig4data[,2], col=colorPallet[1], main='Wolves', xlab='Year', ylab='Individuals', cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
points(out_fig4data[,1],out_fig4data[,2],type='l',col=colorPallet[2],lwd=3)
points(out_est1[,1],out_est1[,2],col=colorPallet[3],type='l',lwd=3)
points(out_est1_1[,1],out_est1_1[,2],col='orange',type='l',lwd=3)
points(out_est1_2[,1],out_est1_2[,2],col='orange',type='l',lwd=3)
points(fig4data[,1],MARest_orig$MAR_TS_Est[,1],col=colorPallet[4],type='l',lwd=3)
points(fig4data[,1],MARest_log$MAR_TS_Est[,1],col=colorPallet[5],type='l',lwd=3)

plot(fig4data[,1], fig4data[,3], col=colorPallet[1], main='Moose', xlab='Year', ylab='Individuals', cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
points(out_fig4data[,1],out_fig4data[,3],type='l',col=colorPallet[2],lwd=3)
points(out_est1[,1],out_est1[,3],col=colorPallet[3],type='l',lwd=3)
points(out_est1_1[,1],out_est1_1[,3],col='orange',type='l',lwd=3)
points(out_est1_2[,1],out_est1_2[,3],col='orange',type='l',lwd=3)
points(fig4data[,1],MARest_orig$MAR_TS_Est[,2],col=colorPallet[4],type='l',lwd=3)
points(fig4data[,1],MARest_log$MAR_TS_Est[,2],col=colorPallet[5],type='l',lwd=3)

plot(fig4data[,1],fig4data[,4],col=colorPallet[1],main='Fir tree rings', xlab='Year', ylab='tree growth', cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
points(out_fig4data[,1],out_fig4data[,4],type='l',col=colorPallet[2],lwd=3)
points(out_est1[,1],out_est1[,4],col=colorPallet[3],type='l',lwd=3)
points(out_est1_1[,1],out_est1_1[,4],col='orange',type='l',lwd=3)
points(out_est1_2[,1],out_est1_2[,4],col='orange',type='l',lwd=3)
points(fig4data[,1],MARest_orig$MAR_TS_Est[,3],col=colorPallet[4],type='l',lwd=3)
points(fig4data[,1],MARest_log$MAR_TS_Est[,3],col=colorPallet[5],type='l',lwd=3)

plot(1,1,col='white',xlab='',ylab='',axes = FALSE)
legend(.6,1.3,legend = c('data','Muhlbauer est','Slope est - all years','Slope est - minus 1981,1982','MAR','MAR logTrans'),lty = c(NA,1,1,1,1,1),lwd=c(NA,3,3,3,3,3),pch = c(1,NA,NA,NA,NA,NA),col = c(colorPallet[1:3],'orange',colorPallet[4:5]),bty = "n",cex=1) # ,
# dev.off()

Muhlbauer = out_fig4data[,1:4]
sum( (fig4data[,2:4]-Muhlbauer[,2:4])^2 )

sum( (fig4data[,2:4]-out_est1[,2:4])^2 )
sum( (fig4data[,2:4]-MARest_orig$MAR_TS_Est)^2 )
sum( (fig4data[,2:4]-MARest_log$MAR_TS_Est)^2 )

##### fig4 - minus 1981, 1982
######################################################################################################################################

