plot.MFVAR.nowcast.news <- function(nowcast_news)
{
  pos <- nowcast_news$smoother_update
  neg <- nowcast_news$smoother_update
  pos[which(nowcast_news$smoother_update<0)] <- 0
  neg[which(nowcast_news$smoother_update>0)] <- 0
  
  colours <- brewer.pal(12, "Paired")
  
  par(mar = c(5, 4, 4, 10) + 0.1)
  
  plot(NA, ty = "h", lend = 1, ylim = c(sum(neg), sum(pos)), xlim = rep(ncast[1]+ncast[2]/4,2), 
       xlab = NA, ylab = "Percent (annualised)", xaxt="n", bty="L")
  lines(rep(ncast[1]+ncast[2]/4,12), rev(cumsum(pos)), ty = "h", lend = 1, lwd = 100, col = rev(colours))
  lines(rep(ncast[1]+ncast[2]/4,12), rev(cumsum(neg)), ty = "h", lend = 1, lwd = 100, col = rev(colours))
  points(ncast[1]+ncast[2]/4, sum(nowcast_news$smoother_update), pch=16, cex=2)
  #text(ncast[1]+ncast[2]/4, sum(nowcast_news$smoother_update)+1, round(sum(nowcast_news$smoother_update), 1))
  abline(h=0)
  #grid(col="dark grey")
  
  legend_names <- gsub("gdp", "other", colnames(gdp_model$y$pctile_50))
  par(xpd=TRUE)
  legend(ncast[1]+500, sum(nowcast_news$smoother_update), legend=c(legend_names, "gdp"), 
         col = c(colours,"black"), bg = "white", lty = 1, lwd = 4, cex = 1, bty="n")
  par(xpd=FALSE)
  mtext(text="Vintage decomposition", side=1, line=1)
  
  par(mar = c(5, 4, 4, 2) + 0.1)
}