
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
	#lines(MARest1$MARSS_output$states[i,])
	#lines(MARest1$MARSS_output$states[i,] - 2*MARest1$MARSS_output$states.se[i,],lty=3)
	#lines(MARest1$MARSS_output$states[i,] + 2*MARest1$MARSS_output$states.se[i,],lty=3)
	#lines(sim.data$sim.data[i,,1])

	}

plot(MARest1$MARSS_output)

 MARest1$MARSS_output$states

sim.data=MARSSsimulate(MARest1$MARSS_output)
sim.data=MARSSsimulate(MARSS_output1, nsim = 1, tSteps = 100)




	# expected states
	#lines(MARest1_smoothData$MARSS_output$states[i,]+mean(TS_smooth_data[,i+1],na.rm = TRUE))
	#lines(MARest1_smoothData$MARSS_output$states[i,]+mean(TS_smooth_data[,i+1],na.rm = TRUE) - 2*MARest1_smoothData$MARSS_output$states.se[i,],lty=3)
	#lines(MARest1_smoothData$MARSS_output$states[i,]+mean(TS_smooth_data[,i+1],na.rm = TRUE) + 2*MARest1_smoothData$MARSS_output$states.se[i,],lty=3)

	# expected log states
	lines(exp(MARest2_smoothData$MARSS_output$states[i,]+mean(log(TS_smooth_data[,i+1]),na.rm = TRUE)))
	lines(exp(MARest2_smoothData$MARSS_output$states[i,]+mean(log(TS_smooth_data[,i+1]),na.rm = TRUE)) - 2*MARest2_smoothData$MARSS_output$states.se[i,],lty=3)
	lines(exp(MARest2_smoothData$MARSS_output$states[i,]+mean(log(TS_smooth_data[,i+1]),na.rm = TRUE)) + 2*MARest2_smoothData$MARSS_output$states.se[i,],lty=3)






plot(out1[,1],log(out1[,i+1]),type='l',col='lightgrey',lwd=3,ylab = varNames[i])

points(log(TS_smooth_data[,i+1]),pch=20,col=varColor[i],xlab = 'n')						# treated data

lines(MARest2_smoothData$MARSS_output$states[i,]+mean(log(TS_smooth_data[,i+1]),na.rm = TRUE))


