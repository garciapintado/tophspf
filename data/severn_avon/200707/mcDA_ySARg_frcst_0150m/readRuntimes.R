yt <- RF$y$uzs$s$data$time                   # [POSIXct] SAR overpass times
yt <-  yt[yt > RF$staT & yt < RF$endT]

RF$runtimes <- c(RF$staT, yt)                # POSIXct

rm(yt)
