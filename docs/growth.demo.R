#  -----------------------------------------------------------------------
#                            Growth Data Analyses
#  -----------------------------------------------------------------------

#  -----------------------------------------------------------------------
#
#                          Overview of the analyses
#
#  These analyses are intended to illustrate the analysis of nonperiod data
#  where a spline basis is the logical choice.  These analyses complement
#  the daily weather data in that sense.
#
#  The growth data have the additional feature of being essentially
#  monotonic or, to say the same thing in another way, have an essentially
#  positive first derivative or velocity.  This requires monotone smoothing.
#  Moreover, most of the interpretability of the growth data comes from
#  inspecting the acceleration of the height curves, so that great emphasis
#  is placed here on getting a good sensible and stable acceleration
#  estimate.
#
#  Finally, a large prortion of the variation in the growth curve data is 
#  due to phase variation, mainly through the variation in the timing of the
#  pubertal growth spurt.  Registration therefore plays a major role and is
#  especially illustrated here.
#
#  Most of the analyses are carried out on the Berkeley growth data, which
#  have the advantage of being freely distributable, whereas as more recent
#  and larger data bases require special permission from the agencies that
#  are responsible for them.  Not much is lost, however, since the quality
#  of the Berkeley data are quite comparable to those of other datasets.
#  The primary analyses are the monotone smoothing of the data.  The right
#  smoothing level is taken as known, and was determined by other analyses
#  in the Matlab language.  The monotone smoothing function used here
#  requires the use of low-level code in C and C++, but even with that help,
#  computation times are substantially longer than in Matlab.
#  Following monotone smoothing, the growth data are registered, an
#  essential step because of the large variation in the timing of the
#  pubertal growth spurt.  The pubertal growth spurts are aligned using
#  landmark registration, and the land-mark registered curves are then
#  registered using continuous registration.
#  The final analysis is of a set of data on a single boy where the
#  measurements are taken every three days or so, rather than twice a year.
#  These data show that growth is rather more complex than the traditional
#  data could have revealed.
#  -----------------------------------------------------------------------

#  -----------------------------------------------------------------------
#                           Berkeley Growth Data
#  -----------------------------------------------------------------------

#  Last modified 2008.06.21;  previously modified 21 March 2006
###
###
### 0.  Access the data (available in the 'fda' package)
###
###

attach(growth)
(nage <- length(age))
(ncasem <- ncol(hgtm))
(ncasef <- ncol(hgtf))

(ageRng <- range(age))
agefine <- seq(ageRng[1],ageRng[2],length=101)
###
###
### 1.  Smooth the data (ignore monotonicity) --------------
###
###
#  This smooth uses the usual smoothing methods to smooth the data,
#  but is not guaranteed to produce a monotone fit.  This may not
#  matter much for the estimate of the height function, but it can
#  have much more serious consequences for the velocity and
#  accelerations.  See the monotone smoothing method below for a
#  better solution, but one with a much heavier calculation overhead.

#  -----------  Create fd objects   ----------------------------
#  A B-spline basis with knots at age values and order 6 is used

# A single call to smooth.basisPar would give us a cubic spline.  
# However, to get a smooth image of acceleration,
# we need a quintic spline (degree 5, order 6) 

# .... 


hgtm = growth$hgtm
hgtf = growth$hgtf
age = growth$age
rng = range(age)

#from growthsetup.R
#growthdata <- list(hgtm  = hgtm,  hgtf = hgtf, age = age)
#save(growthdata, file = "growthdata")


knots  <- growth$age
norder <- 6
nbasis <- length(knots) + norder - 2
hgtbasis <- create.bspline.basis(range(knots), nbasis, norder, knots)

#  --- Smooth these objects, penalizing the 4th derivative  --
#  This gives a smoother estimate of the acceleration functions

Lfdobj <- 4
lambda <- 1e-2
growfdPar <- fdPar(hgtbasis, Lfdobj, lambda)

# Need 'hgtm', 'hgtf', e.g., from attach(growth)
hgtmfd <- smooth.basis(growth$age, growth$hgtm, growfdPar)$fd
hgtffd <- smooth.basis(growth$age, growth$hgtf, growfdPar)$fd

#  plot data and smooth, residuals, velocity, and acceleration

#  Males:

hgtmfit <- eval.fd(age,     hgtmfd)
hgtmhat <- eval.fd(agefine, hgtmfd)
velmhat <- eval.fd(agefine, hgtmfd, 1)
accmhat <- eval.fd(agefine, hgtmfd, 2)

par(mfrow=c(2,2),pty="s",ask=TRUE)
children <- 1:ncasem
for (i in children) {
    plot(age, hgtm[,i], ylim=c(60,200),
         xlab="Years", ylab="", main=paste("Height for male",i))
    lines(agefine, hgtmhat[,i], col=2)
    resi <- hgtm[,i] - hgtmfit[,i]
    ind  <- resi >= -.7 & resi <= .7
    plot(age[ind], resi[ind], type="b", ylim=c(-.7,.7),
         xlab="Years", ylab="", main="Residuals")
    abline(h=0, lty=2)
    ind <- velmhat[,i] >= 0 & velmhat[,i] <= 20
    plot(agefine[ind], velmhat[ind,i], type="l", ylim=c(0,20),
         xlab="Years", ylab="", main="Velocity")
    abline(h=0, lty=2)
    ind <- accmhat[,i] >= -6 & accmhat[,i] <= 6
    plot(agefine[ind], accmhat[ind,i], type="l", ylim=c(-6,6),
         xlab="Years", ylab="", main="Acceleration")
    abline(h=0, lty=2)
}    

# Females:

hgtffit <- eval.fd(age,     hgtffd)
hgtfhat <- eval.fd(agefine, hgtffd)
velfhat <- eval.fd(agefine, hgtffd, 1)
accfhat <- eval.fd(agefine, hgtffd, 2)

par(mfrow=c(2,2),pty="s",ask=TRUE)
children <- 1:ncasef
for (i in children) {
    plot(age, hgtf[,i], ylim=c(60,200),
         xlab="Years", ylab="", main=paste("Height for female",i))
    lines(agefine, hgtfhat[,i], col=2)
    resi <- hgtf[,i] - hgtffit[,i]
    ind  <- resi >= -.7 & resi <= .7
    plot(age[ind], resi[ind], type="b", ylim=c(-.7,.7),
         xlab="Years", ylab="", main="Residuals")
    abline(h=0, lty=2)
    ind <- velfhat[,i] >= 0 & velfhat[,i] <= 20
    plot(agefine[ind], velfhat[ind,i], type="l", ylim=c(0,20),
         xlab="Years", ylab="", main="Velocity")
    abline(h=0, lty=2)
    ind <- accfhat[,i] >= -6 & accfhat[,i] <= 6
    plot(agefine[ind], accfhat[ind,i], type="l", ylim=c(-6,6),
         xlab="Years", ylab="", main="Acceleration")
    abline(h=0, lty=2)
}
###
###
### 2.  Smooth the data monotonically  
###
###

#  These analyses use a function written entirely in S-PLUS called
#  smooth.monotone that fits the data with a function of the form
#                   f(x) = b_0 + b_1 D^{-1} exp W(x)
#     where  W  is a function defined over the same range as X,
#                 W + ln b_1 = log Df and w = D W = D^2f/Df.
#  The constant term b_0 in turn can be a linear combinations of covariates:
#                         b_0 = zmat * c.
#  The fitting criterion is penalized mean squared error:
#    PENSSE(lambda) = \sum [y_i - f(x_i)]^2 +
#                     \lambda * \int [L W(x)]^2 dx
#  where L is a linear differential operator defined in argument Lfdobj.
#  The function W(x) is expanded by the basis in functional data object
#  Because the fit must be calculated iteratively, and because S-PLUS
#  is so slow with loopy calculations, these fits are VERY slow.  But
#  they are best quality fits that I and my colleagues, notably
#  R. D. Bock, have been able to achieve to date.
#  The Matlab version of this function is much faster.

#  ------  First set up a basis for monotone smooth   --------

#  We use b-spline basis functions of order 6
#  Knots are positioned at the ages of observation.

norder <- 6
nbasis <- nage + norder - 2
wbasis <- create.bspline.basis(rng, nbasis, norder, age)

#  starting values for coefficient

cvec0 <- matrix(0,nbasis,1)
Wfd0  <- fd(cvec0, wbasis)

Lfdobj    <- 3          #  penalize curvature of acceleration
lambda    <- 10^(-0.5)  #  smoothing parameter
growfdPar <- fdPar(Wfd0, Lfdobj, lambda)

#  Set up design matrix and wgt vector

zmat  <- matrix(1,nage,1)
wgt   <- rep(1,nage)

#  ---------------------  Now smooth the data  --------------------

# Males:

cvecm <- matrix(0, nbasis, ncasem)
betam <- matrix(0, 2,      ncasem)
RMSEm <- matrix(0, 1,      ncasem)

attach(growth)
children <- 1:ncasem
for (icase in children) {
   hgt     <- hgtm[,icase]
   smoothList <-
		smooth.monotone(age, hgt, growfdPar, wgt, zmat,
		                conv=0.001, dbglev=0)
   Wfd     <- smoothList$Wfdobj
   beta    <- smoothList$beta
   Flist   <- smoothList$Flist
   iternum <- smoothList$iternum
   cvecm[,icase] <- Wfd$coefs
   betam[,icase] <- beta
   hgthat <- beta[1] + beta[2]*monfn(age, Wfd)
   RMSE   <- sqrt(mean((hgt - hgthat)^2*wgt)/mean(wgt))
   RMSEm[icase] <- RMSE
   cat(c(icase, iternum),paste("  ",round(Flist$f,4),
       "  ",round(RMSE, 4),"\n"))
}

# Females:

cvecf <- matrix(0, nbasis, ncasef)
betaf <- matrix(0, 2, ncasef)
RMSEf <- matrix(0, 1, ncasef)

children <- 1:ncasef
for (icase in children) {
   hgt    <- hgtf[,icase]
   smoothList <-
		smooth.monotone(age, hgt, growfdPar, wgt, zmat,
		                 conv=0.001, dbglev=0)
   Wfd     <- smoothList$Wfd
   beta    <- smoothList$beta
   Flist   <- smoothList$Flist
   iternum <- smoothList$iternum
   cvecf[,icase] <- Wfd$coefs
   betaf[,icase] <- beta
   hgthat <- beta[1] + beta[2]*monfn(age, Wfd)
   RMSE   <- sqrt(mean((hgt - hgthat)^2*wgt)/mean(wgt))
   RMSEf[icase] <- RMSE
   cat(c(icase, iternum),paste("  ",round(Flist$f,4),
       "  ",round(RMSE, 4),"\n"))
}

#  -------------  plot the results  --------------------

#  Males:

par(mfrow=c(2,2),pty="s",ask=TRUE)
children <- 1:ncasem
for (i in children) {
    Wfd  <- fd(cvecm[,i],wbasis)
    beta <- betam[,i]
    hgtmfit <- beta[1] + beta[2]*monfn(age, Wfd)
    hgtmhat <- beta[1] + beta[2]*monfn(agefine, Wfd)
    velmhat <- beta[2]*eval.monfd(agefine, Wfd, 1)
    accmhat <- beta[2]*eval.monfd(agefine, Wfd, 2)
    plot(age, hgtm[,i], ylim=c(60,200),
         xlab="Years", ylab="", main=paste("Height for male",i))
    lines(agefine, hgtmhat, col=2)
    resi <- hgtm[,i] - hgtmfit
    ind  <- resi >= -.7 & resi <= .7
    plot(age[ind], resi[ind], type="b", ylim=c(-.7,.7),
         xlab="Years", ylab="", main="Residuals")
    abline(h=0, lty=2)
    ind <- velmhat >= 0 & velmhat <= 20
    plot(agefine[ind], velmhat[ind], type="l", ylim=c(0,20),
         xlab="Years", ylab="", main="Velocity")
    ind <- accmhat >= -6 & accmhat <= 6
    plot(agefine[ind], accmhat[ind], type="l", ylim=c(-6,6),
         xlab="Years", ylab="", main="Acceleration")
    abline(h=0, lty=2)
}

#  Females:

par(mfrow=c(2,2),pty="s",ask=TRUE)
children <- 1:ncasef
for (i in children) {
    Wfd  <- fd(cvecf[,i],wbasis)
    beta <- betaf[,i]
    hgtffit <- beta[1] + beta[2]*monfn(age, Wfd)
    hgtfhat <- beta[1] + beta[2]*monfn(agefine, Wfd)
    velfhat <- beta[2]*eval.monfd(agefine, Wfd, 1)
    accfhat <- beta[2]*eval.monfd(agefine, Wfd, 2)
    plot(age, hgtf[,i], ylim=c(60,200), 
         xlab="Years", ylab="", main=paste("Height for female",i))
    lines(agefine, hgtfhat, col=2)
    resi <- hgtf[,i] - hgtffit
    ind  <- resi >= -.7 & resi <= .7
    plot(age[ind], resi[ind], type="b", ylim=c(-.7,.7),
         xlab="Years", ylab="", main="Residuals")
    abline(h=0, lty=2)
    ind <- velfhat >= 0 & velfhat <= 20
    plot(agefine[ind], velfhat[ind], type="l", ylim=c(0,20),
         xlab="Years", ylab="", main="Velocity")
    ind <- accfhat >= -6 & accfhat <= 6
    plot(agefine[ind], accfhat[ind], type="l", ylim=c(-6,6),
         xlab="Years", ylab="", main="Acceleration")
    abline(h=0, lty=2)
}


#registration from growthreg

ncasem <- 39
ncasef <- 54

#  set up the basis for function Wfd

rng     <- c(1,18)
nbasisw <- 15
norder  <- 5
basisw  <- create.bspline.basis(rng, nbasisw, norder)

# set up the mean velocity curve as the preliminary target for
#  registration

hgtfmeanfd <- mean.fd(hgtffd)
y0fd <- deriv.fd(hgtfmeanfd,  1)

#  curves to be registered

yfd  <- deriv.fd(hgtffd, 1)

#  set up functional parameter object for function Wfd

coef0  <- matrix(0,nbasisw,ncasef)
Wfd0   <- fd(coef0, basisw)
lambda <- 10
WfdPar <- fdPar(Wfd0, 2, lambda)

#  register the data.  It might be a good idea to disable
#  buffered output in the Misc menu for the R Console in order
#  to track progress of this fairly slow process.

reglist <- register.fd(y0fd, yfd, WfdPar)

yregfd  <- reglist$regfd  #  registered curves
Wfd     <- reglist$Wfd    #  functions defining warping functions

#  evaluate the registered curves and warping functions

agefine <- seq(1, 18, len=101)
ymat    <- eval.fd(agefine, yfd)
y0vec   <- eval.fd(agefine, y0fd)
yregmat <- eval.fd(agefine, yregfd)
warpmat <- eval.monfd(agefine, Wfd)
warpmat <- 1 + 17*warpmat/(matrix(1,101,1)%*%warpmat[101,])

#  plot the results for each girl:
#    blue:  unregistered curve
#    red:   target curve
#    green: registered curve

par(mfrow=c(1,2),pty="s",ask=T)
for (i in 1:ncasef) {
  plot (agefine, ymat[,i], type="l", ylim=c(0,20), col=4,
        xlab="Year", ylab="Velocity", main=paste("Case",i))
  lines(agefine, y0vec, lty=2, col=2)
  lines(agefine, yregmat[,i],  col=3)
  plot (agefine, warpmat[,i], type="l",
        xlab="Clock year", ylab="Biological Year")
  abline(0,1,lty=2)
}

#  Comments:  we see that not all curves are properly registered.
#     Curves 7, 11, 13 and 25, to mention a few, are so far from
#     the target that the registration is unsuccessful.  This
#     argues for a preliminary landmark registration of the 
#     velocity curves prior to the continuous registration 
#     process.  However, we will see some improvement below.  

#  compute the new mean curve as a target

y0fd2 <- mean.fd(yregfd)

#  plot the unregistered mean and the registered mean

par(mfrow=c(1,1),pty="s",ask=F)
plot(y0fd2, col=4, xlab="Year", ylab="Mean Velocity")
lines(y0fd, col=3)
legend(10,15, c("Registered", "Unregistered"), lty=c(1,1), col=c(4,3))

#  Comment:  The new mean has a sharper peak at the pubertal
#      growth spurt, which is what we wanted to achieve.

#  define the registered curves and the new curves to be registered

yfd2 <- yregfd

#  register the curves again, this time to a better target

reglist2 <- register.fd(y0fd2, yfd2, WfdPar)

yregfd2  <- reglist2$regfd  #  registered curves
Wfd2     <- reglist2$Wfd    #  functions defining warping functions

y0vec2   <- eval.fd(agefine, y0fd2)
yregmat2 <- eval.fd(agefine, yregfd2)
warpmat2 <- eval.monfd(agefine, Wfd2)
warpmat2 <- 1 + 17*warpmat2/(matrix(1,101,1)%*%warpmat2[101,])

#  plot the results for each girl:
#    blue:  unregistered curve
#    red:   target curve
#    green: registered curve

par(mfrow=c(1,2),pty="s",ask=T)
for (i in 1:ncasef) {
  #  plot velocity curves
  plot (agefine, ymat[,i], type="l", ylim=c(0,20), col=4,
        xlab="Year", ylab="Velocity", main=paste("Case",i))
  lines(agefine, y0vec2, lty=2, col=2)
  lines(agefine, yregmat[,i],   col=3, lty=3)
  lines(agefine, yregmat2[,i],  col=3)
  #  plot warping functions
  plot (agefine, warpmat2[,i], type="l",
        xlab="Clock year", ylab="Biological Year")
  abline(0,1,lty=2)
}

#  compute the new mean curve as a target

y0fd3 <- mean.fd(yregfd2)

#  plot the unregistered mean and the registered mean

par(mfrow=c(1,1),pty="s",ask=F)
plot(y0fd3, col=4, xlab="Year", ylab="Mean Velocity")
lines(y0fd2, col=3)
lines(y0fd, col=3, lty=3)
legend(10,15, c("Registered twice", "Registered once", "Unregistered"), 
       lty=c(1,1,3), col=c(4,3,3))

#  Comment:  The second round of registered made hardly any
#    difference for either the individual curves or the mean curve.


