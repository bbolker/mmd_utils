## from https://stat.ethz.ch/pipermail/r-sig-mixed-models/2011q2/013602.html
require(gamm4); require(mgcv); require(lme4)
set.seed(0)
dat <- gamSim(6,n=400,scale=.2,dist="poisson")
b<-gamm4(y~s(x0), random = ~ (1|fac), family=poisson,
         data=dat)
r1<-as.vector(residuals(b$gam))
r2<-as.vector(b$gam$residuals)
r3<-as.vector(b$mer at resid)
r4<-as.vector(log(dat$y)-predict(b$gam))

residsDF<-data.frame(
	values = c(r1,r2,r3,r4),
	group = rep(c("r1","r2","r3","r4"), each = nrow(dat))
)

with(residsDF, boxplot(values ~ group))

assign("differences", with(residsDF,
	data.frame(
		d1=mean(r1-r2),
		d2=mean(r1-r3),
		d3=mean(r1-r4),
		d4=mean(r2-r3),
		d5=mean(r2-r4),
		d6=mean(r3-r4)
	)))
differences
  #so r2 and r3 are equivalent
