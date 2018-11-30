uzfrac <- function(uzrat) {
    uzfrac <- uzrat * 0.0
    b2 <- uzrat < 2.0          # below 2
    k1 <- 3.0 - uzrat[b2]
    uzfrac[b2] <- 1.0 - (uzrat[b2]*0.5)*((1.0/(4.0 - uzrat[b2]))^k1)
    k2 <- 2.0 * uzrat[!b2] - 3.0
    uzfrac[!b2] <- (1.0/(1.0 + k2))^k2
    return(uzfrac)
}

uzrat <- seq(from=0,to=10,by=0.1)
plot(uzrat, uzfrac(uzrat), type='l', col='blue')     # very similar to a Gaussian, but with a mor compact support, being vey close to 0 for uzrat=3
                                                     # and uzfrac=0.5 for uzrat=2
