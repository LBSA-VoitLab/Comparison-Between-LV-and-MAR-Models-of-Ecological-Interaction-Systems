### t.test ###

### data ###
alg = c(
1.160,
0.226,
1221.929,
385.279,

2515.934,
4781.476,
2190.5862,
1118350,
2251804,
153703251,
5556585
)

LR = c(
1.076,
0.349,
1360.438,
754.249,

2588.491,
6203.687,
24182.900,
17700169,
12114811,
168134153,
32326913
)

MAR = c(
3.674,
1.513,
636.272,
1144.325,

18159.950,
13016.262,
4618.7096,
2199174,
4245588,
217478452,
10267681
)

MAR_logTrans = c(
5.335,
1.395,
381.430,
454.437,

51738.263,
22162.242,
10146.0397,
1777559,
4450871,
166186993,
10947294
)

MAR_smooth = c(
2.638,
1.394,
830.475,
1070.461,

17876.615,
9292.995,
4779.6096,
9904281,
5596187,
143728190,
11136306
)

MAR_logTrans_smooth = c(
5.017,
1.481,
1141.293,
583.539,

16398.792,
12446.589,
8795.03,
3973876,
4554089,
212424910,
14342449
)


cbind(alg,LR,MAR)
cbind(LR,MAR,LR-MAR)
cbind(alg,MAR,alg-MAR)


hist(MAR)
plot(alg[1:6],MAR[1:6])
boxplot(cbind(alg[1:6],MAR[1:6]))
plot(LR[1:6],MAR[1:6])
boxplot(cbind(LR[1:6],MAR[1:6]))
plot(LR[1:6],alg[1:6])

?wilcox.test

wilcox.test(LR,MAR,mu=0,alternative = "less",paired=TRUE,conf.level=.95,exact=TRUE)				# H0: LR is equal or higher than MAR	
wilcox.test(alg,MAR,mu=0,alternative = "less",paired=TRUE,conf.level=.95,exact=TRUE)			# H0: alg is iqual or higher than MAR
wilcox.test(alg,MAR_logTrans,mu=0,alternative = "less",paired=TRUE,conf.level=.95,exact=TRUE)		# H0: alg is iqual or higher than MAR_logTrans
wilcox.test(alg,MAR_smooth,mu=0,alternative = "less",paired=TRUE,conf.level=.95,exact=TRUE)		# H0: alg is iqual or higher than MAR_smooth
wilcox.test(alg,MAR_logTrans_smooth,mu=0,alternative = "less",paired=TRUE,conf.level=.95,exact=TRUE)	# H0: alg is iqual or higher than MAR_logTrans_smooth
wilcox.test(alg,LR,mu=0,alternative = "less",paired=TRUE,conf.level=.95,exact=TRUE)				# H0: alg is equal or higher than LR	
wilcox.test(MAR,MAR_logTrans,mu=0,alternative = "less",paired=TRUE,conf.level=.95,exact=TRUE)		# H0: MAR is iqual or higher than MAR_logTrans
wilcox.test(MAR_logTrans,MAR,mu=0,alternative = "less",paired=TRUE,conf.level=.95,exact=TRUE)		# H0: MAR_logTrans is iqual or higher than MAR
wilcox.test(MAR_smooth,MAR,mu=0,alternative = "less",paired=TRUE,conf.level=.95,exact=TRUE)		# H0: MAR is iqual or higher than MAR_smooth
wilcox.test(MAR_logTrans_smooth,MAR_logTrans,mu=0,alternative = "less",paired=TRUE,conf.level=.95,exact=TRUE)	# H0: MAR_logTransform is iqual or higher than MAR_logTrans_smooth




wilcox.test(LR,alg,mu=0,alternative = "less",paired=TRUE,conf.level=.95,exact=TRUE)




##### Pars absDif

parLR = c(
0.0143,
0.0070,
0.0002,
0.2083,
1.0698,
0.1918,
0.8240,
0.1879
)

parALG = c(
0.0139,
0.0056,
0.0004,
0.1553,
0.0982,
0.6080,
0.3621,
0.9779
)

wilcox.test(parLR,parALG,mu=0,alternative = "less",paired=TRUE,conf.level=.95,exact=TRUE)		# H0: parLR is equal or higher than parALG	




##### noisy vs replicates

noisy_point_fits = c(
1.076,
1.160,
3.674,
5.335,
2.638,
5.017,
1360.438,
1221.929,	
636.272,	
381.430,	
830.475,
1141.293
)

reps_point_fits = c(
0.349,
0.226,
1.513,
1.395,
1.394,
1.481,
754.249,
385.279,
1144.325,
454.437,
1070.461,
583.539
)

wilcox.test(reps_point_fits,noisy_point_fits,mu=0,alternative = "less",paired=TRUE,conf.level=.95,exact=TRUE)		# H0: reps is equal or higher than noisy	



noisy_pars_errors = c(
0.0143,
0.01392
)

reps_par_errors = c(
0.0070,
0.0056
)

wilcox.test(reps_par_errors,noisy_pars_errors,mu=0,alternative = "less",paired=TRUE,conf.level=.95,exact=TRUE)		# H0: reps is equal or higher than noisy	

