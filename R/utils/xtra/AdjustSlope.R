AdjustSlope <- function(x) {
    # x :: input slope [0/1]
    xsig <- ifelse(x >= 0, 1, -1)
    x <- abs(x)
    xbig <- x >= 0.04
    x[xbig] <- 0.05247 + 0.06363 * x[xbig] - 0.182 * exp (-62.38 * x[xbig])
    xsig * x
}

# test
if (1 > 2) {
  x <- seq(-0.5,0.5,by=0.001)
  plot(x,AdjustSlope(x),type='l',col='blue')
  abline(0,1,col='grey')
,}


# initial condition
nx <-  10
ny <-  20
dx <- dy <- 500

h <- matrix(0,ny,nx)

KinematicWave <- function(dx,dy) {

                                        # [i,j] == ['','right']
                                        # [,j-1] = ['','left']
  for (j in 1:nx) {
    for (i in 1:ny) {
      if (wave=='kine') {
        levleft  <- topo[i,j-1]               # for kinematic this is not needed here every timestep
        levright <- topo[i,j]
      } else if(wave='diff') {
        levleft  <- topo[i,j-1] + h[i,j-1]
        levright <- topo[i,j]   + h[i,j]
      }
      # slope
      slope <- AdjustSlope((levleft-lefright) / dx)
      # hydraulic radius assumed as h


    } # for j
  }   # for i
                                       # slope




}
