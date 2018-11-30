# R prototype for the 1D kinematic wave: 1D domain

# goals:
# - establish formulation: based on Chow nonlinear explicit scheme
#   That is a backward diference in both space and time, evaluated at (i+1,t+1)
# - Nonlinear scheme. Evaluated at (i+1, j+1) with
# - Both spatial & temporal as backward differences
# - Directly uses A_t + Q_x = q, so no need to substitute Q in the PDE

# - test difussion
#   con: there is what seems a reasonable amount of numerical difussion. Accept that this account for real difussion.
# - test re-initialitation : sink & sources

                                        # input
intt <- 24*3600    # [s] total integration time
dto  <- 15*60      # [s] outer timestep [15' intervals]
nto  <- intt / dto # [-] number of outer timesteps: number of times that the routing function will be called

nx <-  200
dx <-   50  # m

slo0  <- 0.05     # [0/1] bottom slope
nman  <- 0.035    # [-]   manning coefficient

a       <- (nman * dx^(2./3.) / sqrt(slo0))^0.6     # alpha
b       <- 0.6                                      # beta
epsilon <- 1.0E-10
maxiter <- 100

dt   <- 60        # [s] inner timestep [evaluate sensitivity to this]
# rating curve: plot(0:200,a*(0:200)^b/dx,type='l',col='navyblue', xlab='Q [m3/s]', ylab='h [m]', main='overland flow rating curve', las=1)
nt <- nto*dto/dt
ts <- 1:nt * dt

timeo <- 1:nto * dto                                      # time outputs [at outer times]
rain  <- 60/1000/3600 * exp(-(2*3600-timeo)^2/3600^2)     # [m/s]  60mm/h peak
qin   <- 50           * exp(-(4*3600-timeo)^2/1800^2)     # [m3/s] upstream boundary conditions. 50 m3/s peak
qin   <- 50           * exp(-(4*3600-ts)^2/1800^2)        # [m3/s] upstream boundary conditions. 50 m3/s peak
ql    <- rep(1/100,nx)                                    # [m2/s] lateral inflow constant 1 m3/s / 100 m  inflow along the domain
ql[]  <- 0                                                # disregard ql here
ckts  <- matrix(NA, nrow=nt, ncol=2)                      # output maximum celerity time series
ckts[,1] <- ts
# ---
calcAlpha <- function(slo0, h, dx, n) {
    # obtain alpha in the rating relation for a kinematic wave [A = alpha * Q ^ beta]
    # assuming a rectangular channel

    # slope: downstream slope [0/1]
    # h : water depth [m]
    # dx: cell resolution
    # n : manning coefficient
    # note that for non-chanelized flow, P <- dx

    P <- dx + 2 * h                   # rectangular channel assumption
    (n * P^(2./3.) / sqrt(slo0))^0.6
  }

Q <- matrix(NA,nrow=nx,ncol=nto)  # state vector evolution [1 column/outer timestep]

ckf <- sqrt(slo0)/nman * 5./3.    # Courant number constant part
#t0 = 0.0
it <- 0
for (ito in 1:nto) { # outer
#for (ito in 1:5) { # outer
    cat('outer timestep:',ito,'\n')
  # sink / sources
  #tnex <- t0 + dto
  nti <- dto / dt

  if (ito == 1) {
    Qb <- rep(0, nx) # initial conditions (Q[,k], being k the time index in the DF scheme; i.e. the state at the previously computed time)
    Qa <- rep(NA,nx)                                       # [,k+1]
  }

  for (k in 1:nti) { # inner routing timestep            #       n: j+1, p: j
    it <- it + 1                                         # internal incremental timestep
    # check Courant condition
    h  <- a * Qb ^ b / dx                                # water depth state vector
    ck <- ckf * max(h)^(2./3.)                           # Courant number
    cat('dt:',dt,'| dx/ck:',dx/ck,'|ck:',ck,'\n')
    ckts[it,2] <- ck

    # routing
    for (i in 0:(nx-1)) {
      #qln <- ql[i]             # q_{i+1}^{k+1}   [L2T-1]
      #qlp <- ql[i]             # q_{i+1}^{k}     [L2T-1]

      # (9.6.14) Chow
      # (9.6.7) linear scheme for 1st estimate
        if (i == 0) {
        Qa[i+1] <- (dt/dx*qin[it] + a * b * Qb[i+1]   * ((Qb[i+1]+qin[it])/2)^(b-1)  + dt * ql[i+1]) /
            (dt/dx + a*b*((Qb[i+1]+qin[it])/2)^(b-1))
        f <- 1000
        Co <- dt / dx * qin[it] + a*Qb[i+1]^b + dt * ql[i+1]                    # (9.6.14)       \frac{\delta t}{\delta x}Q_{i+1}^{k+1}
      } else {
        Qa[i+1] <- (dt/dx*Qa[i]    + a * b * Qb[i+1]   * ((Qb[i+1]+Qa[i])/2)^(b-1)     + dt * ql[i+1]) /
            (dt/dx + a*b*((Qb[i+1]+Qa[i])/2)^(b-1))
        f <- 1000
        Co <- dt / dx * Qa[i] + a*Qb[i+1]^b + dt * ql[i+1]                    # (9.6.14)       \frac{\delta t}{\delta x}Q_{i+1}^{k+1}
      }

      while (abs(f) > epsilon) {
        #cat('f:',f,'\n')
        f <- dt / dx * Qa[i+1] + a*Qa[i+1]^b - Co
        g <- dt / dx + a*b*Qa[i+1]^(b-1)
        Qa[i+1] <- Qa[i+1]  - f / g
        if (is.na(Qa[i+1])) {
          Qa[i+1] <- 0
          f <- 0
        }
      }
      if (Qa[i+1] < 0) { Qa[i+1] <- 0 }
    } # end space loop
    Qb <- Qa
  } # end inner routing timestep
  Q[,ito] <- Qa
} # end outer time loop

## plots

# maximum celerity vs dt
plot(ckts[,1],dx/ckts[,2],type='l',col='navyblue', ylab='Courant number [dx/ck]')
abline(h=dt, col='green4')
lines(ckts[,1],dx/ckts[,2],type='l',col='green2',lty=2)

# flow along the channel
X11(); plot(Q[,1],type='l',ylim=c(0,100))
for(ito in 1:nto) {
    lines(Q[,ito],col='green4',lwd=2)
    abline(h=50)

    Sys.sleep(0.1)
   # lines(Q[,it],col='white')
}

# sequence of hydrographs at several locations in the channel
X11();plot(ts, qin, type='l')      # inflow hydrograph
      for (i in c(1:nx)) {
          lines(timeo, Q[i,],col='navyblue',lty=1)   # outflow hydrograph

          Sys.sleep(0.1)
         lines(timeo, Q[i,],col='white')   # outflow hydrograph

      }

# integral of flow along time for each pixel in the domain
X11();plot(rowSums(Q)*dto, type='l')
      points(1,sum(qin)*dto, col='darkred')
      #lines(rowSums(Q)*dto, col='blue')
ovfs <- a * Q[,nto]^b * dx # [m3] final storage
ovcs <- cumsum(ovfs)        # cummulative final storage from upstream to downstream
lines(rowSums(Q)*dto+ovcs, lty=2) # cummulated storage plus outflow per pixel. The difference between the total inflow (the darkred point and this line
                                  # is the dispersion. As seen, mass is not conserved but the loss is just about 0.01% as maximum)

X11();plot(rowSums(Q)/max(rowSums(Q)), type='l')



