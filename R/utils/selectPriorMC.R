selectPriorMC <- function(MCwpar, glike, m, stat='E', decreasing=TRUE) {
  # +++ purpose +++
  # Select a subset of parameter from a previous simulation,
  # given its corresponding performance statistics
  # The function is currently restricted to a single criterion (i.e. one statistics in 'glike')
  #
  # Record of revisions:
  #   Date          Programmer         Description of change
  #   ----          ----------         ---------------------
  #   15/10/2013    J. Garcia-Pintado  Original code
  #
  # list(len=nhrus), INTENT(IN) :: MCwpar         Each element is a [npar,mi] matrix. NA for non-randomized parameters
  # list(len=nhrus), INTENT(IN) :: glike          Each element is a [nstat,mi] matrix
  # INTEGER, INTENT(IN)         :: m              Current simulation ensemble size
  # CHARACTER, INTENT(IN)       :: stats          Name of the statistics to be used as criterion for parameter selection
  # LOGICAL, INTENT(IN)         :: decreasing     If TRUE, higher values of 'stat' indicate better performance
  # INTEGER                     :: mi             Previous simulation ensemble size
  # INTEGER                     :: npar           Number of parameter, including initial conditions
  # INTEGER                     :: nstat          Number of statistics

  nhrus <- length(MCwpar)
  if (length(glike) != nhrus)
    stop('selectPriorMC :: length(glike) != nhrus')

  npar <- nrow(MCwpar[[1]])
  mi   <- ncol(MCwpar[[1]])

  if (m > mi)
    stop('selectPriorMC :: m > mi')

  if (mi == m)
    return(MCwpar)

  #nstats <- nrow(glike[[1]])
  if (mi != ncol(glike[[1]]))
    stop('selectPriorMC :: glike ensemble size != mi')

  opar <- vector('list', nhrus)

  for (iw in 1:nhrus) {
    stats <- glike[[iw]][stat,]
    imd   <- order(stats, decreasing=decreasing)[1:m]
    opar[[iw]] <- MCwpar[[iw]][,imd]
  }
  return(opar)
}
