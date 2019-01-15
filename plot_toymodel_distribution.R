bin = 100000
 x = seq(0.0,1.0,length=bin)
 peak1 = rep(0.5,bin*0.4)
 vallay = rep(0.000001, bin*0.5)
 peak2_height = 8.0 - 0.000001 * bin * 0.5
 peak2 = rep(peak2_height, bin*0.1)
 y = c(peak1, vallay, peak2)
 y[1] = 0
 y[bin] = 0
#  plot(x, y, type="l", xlab="Parameter Value", ylab="Probability Density")
 
 # plot
 library(ggplot2)
 library(ggsci)
 data = data.frame(ParameterValue=x, ProbabilityDensity=y)
 g = ggplot(data, aes(x = ParameterValue, y = ProbabilityDensity))
 g = g + geom_line(size=1)
 plot(g)
 ggsave(file = "toydata_dist.png", plot=g)