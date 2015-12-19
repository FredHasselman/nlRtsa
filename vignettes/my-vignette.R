## ---- eval=FALSE---------------------------------------------------------
#  require(devtools)
#  install_github("FredHasselman/nlRtsa")
#  library(nlRtsa)

## ------------------------------------------------------------------------
in.IT(c("signal","pracma"))

## ------------------------------------------------------------------------
# Sawtooth
x <- seq(-3.2,3.2, length.out = 256)
y <- 2*sin(10*x) - 1*sin(20*x) + (2/3)*sin(30*x) - (1/2)*sin(40*x) + (2/5)*sin(50*x) - (1/4)*sin(60*x)

# Plot the sawtooth wave as constructed by the Fourier series above
plot(x,y, xlab ='Time (a.u.)', ylab = 'Variable (a.u.)', main ='Sawtooth wave', type = "l")

# Perform a Fast Fourier Transform and calculate the Power and Frequency
Y <- fft(y)
Pyy <- Y*Conj(Y)/256
f <- 1000/256*(0:127)

# Plot the power spectrum of the sawtooth wave
plot(f[1:50],Pyy[1:50], type="b",xlab='Frequency (a.u.)', ylab ='Power (a.u.)', pch=21, bg='grey60', main = 'Power Spectrum')

## ------------------------------------------------------------------------
# A time vector
t <- pracma::linspace(x1 = 0, x2 = 50, n = 256)
# There are three sine components
x <- sin(2*pi*t/.1) + sin(2*pi*t/.3) + sin(2*pi*t/.5)
# Add random noise!
y <- x + 1*randn(size(t))

# Plote the noise.
plot(t, y, type = "l", xlab = 'Time (a.u.)', ylab = 'Variable (a.u.)', main = 'A very noisy signal')

# Get the frequency domain
Y <- fft(y)
Pyy <- Y*Conj(Y)/256
f <- 1000/256*(0:127)

# Plot the power spectrum of this noise
plot(f[1:50],Pyy[1:50], type="b",xlab='Frequency (a.u.)', ylab='Power (a.u.)', pch=21, bg='grey60', main = 'Power Spectrum')

