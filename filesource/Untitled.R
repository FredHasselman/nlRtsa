```{r set-options, echo=FALSE,include=FALSE}
require(knitr)
require(formatR)
require(dplyr)
options(width=200)
knitr::opts_chunk$set(cache=FALSE,prompt=FALSE,comment=">",message=FALSE,echo=TRUE,warning=FALSE,tidy=TRUE,strip.white=TRUE,size="small", fig.align = "center",fig.show='hold')

```



{D_{PSD} ≈ 3/2 + ((14/33)*tanh(sap*log(1+sqrt(2))))}

# Create a plot with marginal distributions.
library(ggplot2)
library(gridExtra)
library(scales)
df <- data.frame(x = rnorm(n = 100), y = rnorm(n = 100), group = factor(sample(x=c(0,1), size = 100, replace = TRUE)))
scatterP <- ggplot(df, aes(x = x, y =y, colour = group)) + geom_point() + gg.theme()
xDense <- ggplot(df, aes(x = x, fill = group)) + geom_density(aes(y= ..count..),trim=F,alpha=.5) + gg.theme("noax") + theme(legend.position = "none")
yDense <- ggplot(df, aes(x = y, fill = group)) + geom_density(aes(y= ..count..),trim=F,alpha=.5) + coord_flip() + gg.theme("noax") + theme(legend.position = "none")


grid.arrange(xDense, gg.plotHolder(), scatterP, yDense, ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))

plot(lmSimulate(lmModel("ppl")))

as.numeric(timeLag(beamchaos, method="mutual", plot=TRUE))


beam.d1 <- infoDim(beamchaos, dim=10)

## print a summary of the results
print(beam.d1)

## plot the information dimension curves without
## regression lines
plot(beam.d1, fit=FALSE, legend=FALSE)

## plot an extended data analysis plot
eda.plot(beam.d1)

in.SIGN()



## create test series
set.seed(100)
x <- rnorm(1024)
walk <- cumsum(x)

## calculate the Hurst coefficient of a random
## walk series using various techniques
methods <- c("aggabs","aggvar","diffvar","higuchi")
z <- lapply(methods, function(method, walk){
  hurstBlock(ifelse1(method=="higuchi",diff(walk),walk), method=method)
},walk=walk )
names(z) <- methods

## plot results
old.plt <- splitplot(2,2,1)
for (i in 1:4){
  if (i > 1)
    splitplot(2,2,i)
  plot(z[[i]], key=FALSE)
  mtext(paste(attr(z[[i]],"stat.name"), round(as.numeric(z[[i]]),3), sep=", H="),
        line=0.5, adj=1)
}
par(old.plt)

## calculate various SDF estimates for the
## sunspots series. remove mean component for a
## better comparison.

require(ifultools)
data <- as.numeric(sunspots)
methods <- c("direct","wosa","multitaper",
             "lag window")

S <- lapply(methods, function(x, data) SDF(data, method=x), data)
x <- attr(S[[1]], "frequency")[-1]
y <- lapply(S,function(x) decibel(as.vector(x)[-1]))
names(y) <- methods

## create a stack plot of the data
stackPlot(x, y, col=1:4)


#
# PSDslope <- function(y = ts(rnorm(n = 100), frequency = 1), fs = frequency(y), nfft = 1024, fitRange = round(.1*nfft), plot = FALSE){
#   require(oce)
#
#   if(!is.ts(y)) ts(y, frequency = fs)
#
#   perioGram <- pwelch(x = y, nfft = nfft, fs = fs, plot = FALSE)
#   psd <- data.frame(perioGram$freq,perioGram$spec)
#   psd[1,]=NA
#   slope <- lm(log10(psd[1:fitRange,2])~log10(psd[1:fitRange,1]))
#
#   return(slope)
# }

PSDslope()

in.IT("animation")
saveGIF({
  ani.options(nmax = 30)
  brownian.motion(pch = 21, cex = 5, col = "red", bg = "yellow")
}, interval = 0.05, movie.name = "bm_demo.gif", ani.width = 600, ani.height = 600)


prm <- round(c(seq(-2,4,by = 0.05), seq(3.9,-2,by=-0.05)), digits = 2)
saveGIF({
  ani.options(nmax = length(prm))
  for(r in prm){
    tsdisplay(growth.ac(r=r, type = "logistic"), lag.max = 100, plot.type = "scatter", main = paste0('Logistic Growth\n r = ',r))
  }
}, interval = 0.1, movie.name = "logist.gif", ani.width = 600, ani.height = 600)

saveGIF({
  brownian.motion(n=1, pch = 21, cex = 5, col = "red", bg = "yellow")
}, movie.name = "brownian_motion.gif", interval = 0.1, nmax = 100, ani.width = 600,
ani.height = 600)



saveHTML({
  par(mar = c(3, 3, 1, 0.5), mgp = c(2, 0.5, 0), tcl = -0.3,
      cex.axis = 0.8, cex.lab = 0.8, cex.main = 1)
  ani.options(interval = 0.05, nmax = ifelse(interactive(),
                                             150, 10))
  brownian.motion(pch = 21, cex = 5, col = "red", bg = "yellow")
}, description = c("Random walk on the 2D plane: for each point",
                   "(x, y), x = x + rnorm(1) and y = y + rnorm(1)."),
title = "Demonstration of Brownian Motion")

ani.options(oopt)




Fs <- 1000
t <- seq(0, 0.296, 1/Fs)
x <- cos(2 * pi * t * 200) + rnorm(n=length(t))
xts <- ts(x, frequency=Fs)
s <- spectrum(xts, spans=c(3,2), main="random + 200 Hz", log='no')
w <- pwelch(xts, nfft = 64, plot=FALSE)
lines(w$freq, w$spec, col="red")

x <- FNS(beamchaos, dim=10, tlag=10, olag=15)

## print the results
print(x)



require(fractal)
old.plt <- par("plt")
models <- c("ppl","fdp","fgn","dfbm")
for (i in seq(along=models)){
  splitplot(2,2,i)
  plot(lmSDF(lmModel(models[i])),
       reference.grid=FALSE, log.axes="xy")
}
par(plt=old.plt)

s <- fd.psd(rnorm(1024))
plot.loglog()
