


main_lab = c( expression(paste(italic('Paramecium caudatum'), " (Prey)")),expression(paste(italic('Didinium nasutum'), " (Predator)")) )


y_labs = c(expression(italic('X'[1])),expression(italic('X'[2])),expression(italic('X'[3])),expression(italic('X'[4])))
plot(1,1,ylab=y_lab)


mtext("Awesome Y variable", side=2, line=2.0, cex=3)



	mtext(expression(italic('t')),side=1,line=2.1,cex=1.5)
	mtext(y_labs[i-1],side=2,line=2.1,cex=1.5)
	mtext(mainLab[i-1],side=3,line=1,cex=1.5)