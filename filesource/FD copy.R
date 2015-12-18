# R version of some analyses reported in:
# Hasselman, F.. (2013). When the blind curve is finite: dimension estimation and model inference based on empirical waveforms. Frontiers in
#  Physiology, 4(April), 10-13.

BCdata <- read.delim("blindcurve.dat", header=FALSE)
df.y  <- read.delim("planarext.dat", header=FALSE)
df.y  <- as.data.frame(y[,2:length(y)])
c  <- ncol(y)

FDdata <- data.frame(cbind(ts=(1:c), PSDsap=NA, PSDfd=NA, DFAsap=NA, DFAfd=NA, SDAsap=NA, SDAfd=NA))

# Conversion formulas for self-affinity parameter estimates (sap) to Dimension (fd) suggested in Hasselman (2013)
# PSD slope (if signal is continuous (sampled) consider Wijnants et al. (2013) log-log fitting procedure)
psd2fd <- function(sap){return(fd <- 3/2 + ((14/33)*tanh(sap*log(1+sqrt(2)))) )}
# DFA slope (this is H in DFA)
dfa2fd <- function(sap){return(fd <- 2-(tanh(log(3)*sap)) )}
# SDA slope (simple 1-sap, but note that for some signals different PSD slope values project to 1 SDA slope)
sda2fd <- function(sap){return(fd <- 1-sap)}

cnt=0
for(k in 1:c){
  # time series for PSD and SDA are detrended and normalized using SD based on N
  N             <- length(y[,k])
  ys            <- as.vector(detrend(y[,k]))
  ys            <- (ys - mean(ys)) / (sd(ys)*sqrt((N-1)/N))

  out           <- Dpsd(ys)
  FDdata[k,1:3] <- c(k,out,psd2fd(out))
  rm(out)

  ifelse(k==12, dmethod<-"poly1", dmethod<-"bridge")
  out           <- DFA(y[,k], detrend=dmethod, sum.order=1, scale.max=trunc(length(y[,k])/4), scale.min=4, scale.ratio=2^(1/4), verbose=FALSE)
  FDdata[k,4:5] <- c(attributes(out)$logfit[]$coefficients['x'],dfa2fd(attributes(out)$logfit[]$coefficients['x']))
  rm(out)

  out           <- dispersion(ys)
  maxbin        <- length(out$sd)-2
  lsfit         <- lm(log(out$sd[2:maxbin]) ~ log(out$scale[2:maxbin]))
  FDdata[k,6:7] <- c(coef(lsfit)[2],sda2fd(coef(lsfit)[2]))
  rm(out,lsfit,ys)
}

# Check differences to Blind Curve data from Matlab
Diff<-FDdata-BCdata
FDiff<-as.data.frame(rbind(cbind(1:12,Diff$PSDfd,seq(1,1,length.out=12)),cbind(1:12,Diff$DFAfd,seq(2,2,length.out=12)),cbind(1:12,Diff$SDAfd,seq(3,3,length.out=12))))
names(FDiff)<-c("ts","difference","analysis")
FDiff$analysis<-factor(FDiff$analysis)
levels(FDiff$analysis)[1] <- "PSD"
levels(FDiff$analysis)[2] <- "DFA"
levels(FDiff$analysis)[3] <- "SDA"

ggplot(FDiff) +  geom_hline(xintercept = 0) +
  geom_dotplot(aes(x=factor(ts), y=difference, fill=factor(analysis)), binaxis = "y", stackdir = "center",position="jitter")+
   theme_bw(base_size=12, base_family="Arial") +
   scale_y_continuous("FD difference: Rplus - Matlab",breaks=seq(-0.45,0.15,by=0.05),limits=c(-0.45,0.15),expand=c(0,0)) +
   xlab("Time Series used in Hasselman (2013)") +
   guides(fill=guide_legend(title="Analysis")) +
   theme(axis.line    = element_line(colour="black"),
         axis.title.x = element_text(vjust=-1.1),
         plot.margin  = unit(c(1,1,1,1), "lines"),
         panel.grid.major = element_line(colour = "grey"))


# stackPlot(x=2:12,y=cbind(Diff[2:12,3],Diff[2:12,5],Diff[2:12,7]),main="Difference between informed FD estimates obtained from Rplus & Matlab", xlab="Time Series",ylab=c("iFD: PSD","iFD: DFA","iFD: SDA"),ylim=c(-.08,.13),lwd=c(3,3,3))

x <- FNN(beamchaos, tlag=10, olag=15 )

## print the results
print(x)

## plot the results
plot(x)


plot(ecgrr)

DFA.walk <- DFA(rnorm(1024), detrend="poly1", sum.order=1)

## print the results
print(DFA.walk)

## plot a summary of the results
eda.plot(DFA.walk)


z <- dispersion(rnorm(1024))
plot(log(z$scale),log(z$sd))
