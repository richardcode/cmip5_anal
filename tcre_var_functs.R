#Create the budget distribution calculation function
calc_budget_dist <- function(warm_thresh,temps,emms,years) {
    periods <- c()
    budgets <- c()
    end_points <- c()
    for (i in c(1:length(years))){
        start_temp <- temps[i]
        start_emms <- emms[i]
        if (is.na(start_temp)!=TRUE & is.na(start_emms)!=TRUE){
            diff_temps <- temps - start_temp
            
            exceed_year <- years[ (diff_temps>=warm_thresh & years>years[i]) ][1]
            if (is.na(exceed_year)!=TRUE) {
                #Linearly interpolate exactly when warming threshold crossed
                t_exceed <- approx(x=c(diff_temps[years==(exceed_year-1)],diff_temps[years==exceed_year]),y=c(exceed_year-1,exceed_year),xout=warm_thresh)$y
                #Calculate the correct cumulative emissions up to this point
                d_cumems <- sum(emms[(years>=(years[i]+1) & (years<=(exceed_year-1)))]) + 0.5*emms[years==years[i]] + (t_exceed-(exceed_year-1))*emms[years==(exceed_year-1)]
                period <- t_exceed - years[i]
                
                periods <- c(periods,period)
                budgets <- c(budgets,d_cumems)
                end_points <- c(end_points,exceed_year)
                
            }
        }
    }
    return(list(periods,budgets,end_points))
}


calc_budget_dist_allmods <- function(warm_thresh,df,temps_n,emms_n,years) {
out <- list(c(),c(),c())
for (j in c(1:nrow(df[df$Variable==temps_n,]))){
    out_e <- calc_budget_dist(warm_thresh,as.numeric(df[df$Variable==temps_n,][j,-c(1:5)]),as.numeric(df[df$Variable==emms_n,][j,-c(1:5)]),years)
    if (length(out_e[[2]]>0)){
        out[[1]]<- c(out[[1]],out_e[[1]])
        out[[2]]<- c(out[[2]],out_e[[2]])
        out[[3]]<- c(out[[3]],out_e[[3]])
    }
}
return(out)
}

calc_budget_dist_obsens <- function(warm_thresh,temps_e,emms_e,years) {
    out <- list(c(),c(),c())
    for (j in c(1:nrow(temps_e))){
        out_e <- calc_budget_dist(warm_thresh,temps_e[j,],emms_e[j,],years)
    if (length(out_e[[2]]>0)){
        out[[1]]<- c(out[[1]],out_e[[1]])
        out[[2]]<- c(out[[2]],out_e[[2]])
        out[[3]]<- c(out[[3]],out_e[[3]])
    }
    }
    return(out)
}
