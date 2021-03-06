plot.VAR.shocks <- function(shocks, series_name=NA, start_plot=NA, type="const", quarterly=FALSE, bar_width=10, right_margin=10, bar_col=NA, Ylab=NA, Xlab=NA, leg=NA) 
{
  # Number of lines of margin,  c(bottom, left, top, right)
  par(mar = c(5, 4, 4, right_margin) + 0.1)
  
  #select variable to plot
  if(quarterly){
    series <- lapply(shocks[[series_name]], mon2qtr)
  } else {
    series <- shocks[[series_name]]
  }
  
  #select type to plot
  if(is.na(start_plot)) {
    series.type <- series[[type]]
    start_plot <- tsp(series[[type]])[1]
  } else {
    series.type <- window(series[[type]], start=start_plot)
  }
  index.comp  <- seq(tsp(series.type)[1], tsp(series.type)[2], 1/tsp(series.type)[3])
  components  <- t(series.type)
  
  neg <- components*(components < 0)
  pos <- components*(components > 0)
  pos.cum <- array(NA, dim = dim(components))
  neg.cum <- array(NA, dim = dim(components))
  pos.cum[1, ] <- pos[1, ]
  neg.cum[1, ] <- neg[1, ]
  for(i in 2:nrow(components))
  {
    pos.cum[i, ] <- pos.cum[i-1, ] + pos[i, ]
    neg.cum[i, ] <- neg.cum[i-1, ] + neg[i, ]
  }
  
  pos.cum[pos.cum == 0] <- NA 
  neg.cum[neg.cum == 0] <- NA 
  
  plot.limits <- c(min(apply(neg, 2, sum)), max(apply(pos, 2, sum)))
  plot(index.comp, pos.cum[nrow(components), ], ty = "h", lend = 1,  lwd = bar_width, 
       col = bar_col[length(bar_col)], ylim = plot.limits, xlab = NA, ylab = NA, bty="L")
  lines(index.comp, neg.cum[nrow(components), ], ty = "h", lend = 1, lwd = bar_width, 
        col = bar_col[length(bar_col)])
  for(i in (nrow(components)-1):1) {
    lines(index.comp, pos.cum[i, ], ty = "h", lend = 1, lwd = bar_width, col = bar_col[i])
    lines(index.comp, neg.cum[i, ], ty = "h", lend = 1, lwd = bar_width, col = bar_col[i])
  }
  
  legend_names <- gsub(series_name, "other", rownames(components))
  
  if(is.na(leg)) {
    if(type=="omit") {
      leg <- c(tsp(series.type)[2]+0.5,
               mean(series[["series.omit"]])+1*sd(series[["series.omit"]]))
    } else {
      leg <- c(tsp(series.type)[2]+0.5,
               mean(series[["series"]])+1*sd(series[["series"]]))
    }
    
  }
  
  if(type == "omit") {
    lines(index.comp, window(series[["series.omit"]], start=start_plot), lwd = 3)
    par(xpd=TRUE)
    legend(x=leg[1], y=leg[2], legend=c(legend_names, series_name), 
           col = c(bar_col,"black"), bg = "white", lty = 1, lwd = 4, cex = 0.75, bty="n")   
    par(xpd=FALSE)
    abline(h = 0)
  } else {
    lines(index.comp, window(series[["series"]], start=start_plot), lwd = 3)
    par(xpd=TRUE)
    legend(x=leg[1], y=leg[2], legend=c(legend_names, series_name), 
           col = c(bar_col,"black"), bg = "white", lty = 1, lwd = 4, cex = 0.75, bty="n")  
    par(xpd=FALSE)
    abline(h = 0)
  }
  
  mtext(Xlab, side=1, line=3, cex=1.5, font=2)
  mtext(Ylab, side=2, line=3, cex=1.25, font=2)
  
  # Return to default
  par(mar = c(5, 4, 4, 2) + 0.1)
}