### t.test ###

### data ###
alg = c(
4.362,
2.953,
1141.60,
1462.294,
2515.934,
4781.476,
2190.5862,
1118350,
2251804,
153703251,
5556585
)

LR = c(
9.009,
11.410,
2727.131,
2517.640,
2588.491,
6203.687,
24182.900,
17700169,
12114811,
168134153,
32326913
)

MAR = c(
37.314,
6.247,
908.288,
933.052,
18159.950,
13016.262,
4618.7096,
2199174,
4245588,
217478452,
10267681
)

MAR_logTrans = c(
53.311,
5.439,
37446.550,
177.967,
51738.263,
22162.242,
10146.0397,
1777559,
4450871,
166186993,
10947294
)

MAR_smooth = c(
7.364,
13.874,
1190.739,
1773.447,
17876.615,
9292.995,
4779.6096,
9904281,
5596187,
143728190,
11136306
)

MAR_logTrans_smooth = c(
8.260,
42.831,
1604.136,
1481.969,
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
1.97538,
0.24117,
1.35692,
0.88350,
0.00021,
0.20832,
1.06975,
0.19183,
0.82403,
0.18794
)

parALG = c(
2.68019,
0.32537,
46.60082,
0.66838,
0.00038,
0.15534,
0.09820,
0.60795,
0.36212,
0.97794
)

wilcox.test(parLR,parALG,mu=0,alternative = "less",paired=TRUE,conf.level=.95,exact=TRUE)		# H0: parLR is equal or higher than parALG	




##### noisy vs replicates

noisy_point_fits = c(
9.009,
4.362,
35.657,
55.187,
6885.041,	
1141.595,	
1764.455,	
6885.041
)

reps_point_fits = c(
11.410,
2.953,	
6.947,
6.053,
2075.618,	
1340.290,	
2088.561,	
3264.232
)

wilcox.test(reps_point_fits,noisy_point_fits,mu=0,alternative = "less",paired=TRUE,conf.level=.95,exact=TRUE)		# H0: reps is equal or higher than noisy	



noisy_pars_errors = c(
1.97538,
2.68019,
1.35692,
46.60082)

reps_par_errors = c(
0.24117,
0.32537,
0.88350,
0.66838)

wilcox.test(reps_par_errors,noisy_pars_errors,mu=0,alternative = "less",paired=TRUE,conf.level=.95,exact=TRUE)		# H0: reps is equal or higher than noisy	

