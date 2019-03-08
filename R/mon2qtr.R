mon2qtr <-
function(monthly, aggregation="mean")
{
    month.start = round(tsp(monthly)[1]%%1*12)+1
    year.start = round(tsp(monthly)[1])
    month.end = round(tsp(monthly)[2]%%1*12)+1
    year.end = round(tsp(monthly)[2])

    quarter.start <- NULL
    if(any(month.start == c(2,5,8,11))) {
        message("Dropping two months from start.")
        quarter.start <- month.start + 2
    } else if(any(month.start == c(3,6,9,12))) {
        message("Dropping one month from start.")
        quarter.start <- month.start + 1
    } else {
        quarter.start <- month.start
    }

    quarter.end <- NULL
    if(any(month.end == c(1,4,7,10))) {
        message("Dropping one month from end.")
        quarter.end <- month.end - 1
    } else if(any(month.end == c(2,5,8,11))) {
        message("Dropping two months from end.")
        quarter.end <- month.end - 2
    } else {
        quarter.end <- month.end
    }

    monthly <- window(monthly, 
        start=c(year.start, quarter.start), 
        end=c(year.end, quarter.end))
    
    agg.transform <- function(x,agg="mean") {
        y <- switch(agg,
            "mean" = apply(x,2,mean),
            "sum"  = apply(x,2,sum),
            "last" = apply(x,2,function(z) z[length(z)]))
        return(y)
    }

    quarterly <- NULL
    for(i in 1:NCOL(monthly)) {
        monthly.reshape <- matrix(as.matrix(monthly)[,i],3,NROW(monthly)/3)
        quarterly <- cbind(quarterly, agg.transform(monthly.reshape))
    }

    quarterly <- ts(quarterly,freq=4,start=c(year.start,(quarter.start-1)/3+1))
    colnames(quarterly) <- colnames(monthly)
    
    return(quarterly)
}
