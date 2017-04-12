source("tcre_var_functs.R")

#Load in the data from the CMIP5 ensembles
cmip5_data <- read.csv('./Data/ESM_cmip5_tempems.csv')
#Get models that are both in RCP2.6 and 4xCO2

rcp26_data <- cmip5_data[cmip5_data$Scenario=='RCP26',-c(6:(2004-1850+6))]


warm_threshs <- c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)
for (warm_thresh in warm_threshs) {
    


out_rcp26_pos <- list(c(),c())
out_rcp26 <- list(c(),c())
for (j in c(1:nrow(rcp26_data[rcp26_data$Variable=='Temperature|rel to 1861-80',]))){
    if (sum(as.numeric(rcp26_data[rcp26_data$Variable=='Total anthropogenic carbon flux',][j,-c(1:5)]) < 0.0 , na.rm=TRUE)==0){
        log_pos_emms <- c(1:length(as.numeric(rcp26_data[rcp26_data$Variable=='Total anthropogenic carbon flux',][j,-c(1:5)])))
    } else {
        log_pos_emms <- c(1:(which((as.numeric(rcp26_data[rcp26_data$Variable=='Total anthropogenic carbon flux',][j,-c(1:5)]) < 0.0) )[1]-1))
    }
    
    temps_j <- as.numeric(rcp26_data[rcp26_data$Variable=='Temperature|rel to 1861-80',][j,-c(1:5)])
    temps_j_p <- as.numeric(rcp26_data[rcp26_data$Variable=='Temperature|rel to 1861-80',][j,-c(1:5)])[log_pos_emms]
    emms_j <- as.numeric(rcp26_data[rcp26_data$Variable=='Total anthropogenic carbon flux',][j,-c(1:5)])
    emms_j_p <- as.numeric(rcp26_data[rcp26_data$Variable=='Total anthropogenic carbon flux',][j,-c(1:5)])[log_pos_emms]
    out_rcp26_e <- calc_budget_dist(warm_thresh,temps_j,emms_j,c(2005:2159))
    out_rcp26_e_p <- calc_budget_dist(warm_thresh,temps_j_p,emms_j_p,c(2005:2159)[log_pos_emms])
    out_rcp26[[1]]<- c(out_rcp26[[1]],out_rcp26_e[[1]])
    out_rcp26[[2]]<- c(out_rcp26[[2]],out_rcp26_e[[2]])
    out_rcp26_pos[[1]]<- c(out_rcp26_pos[[1]],out_rcp26_e_p[[1]])
    out_rcp26_pos[[2]]<- c(out_rcp26_pos[[2]],out_rcp26_e_p[[2]])
}

df_all <- data.frame(matrix(nrow=0,ncol=3))
df_all <- rbind(df_all,matrix(c(rep('RCP2.6',length(out_rcp26[[2]])),out_rcp26[[1]],out_rcp26[[2]]/warm_thresh),ncol=3))
df_all <- rbind(df_all,matrix(c(rep('RCP2.6- Positive emissions',length(out_rcp26_pos[[2]])),out_rcp26_pos[[1]],out_rcp26_pos[[2]]/warm_thresh),ncol=3))


df_all[,2] <- as.numeric(as.vector(df_all[,2]))
df_all[,3] <- as.numeric(as.vector(df_all[,3]))

colnames(df_all) <- c('Scenario','Duration','Budget')

p <- ggplot(df_all, aes(Budget, colour = Scenario)) + geom_density(alpha=0.1) + labs(x = "Budget (GtC/°C)")  + ggtitle(paste('Warming: ',as.character(warm_thresh),'°C',sep=''))
ggsave(paste('Figures/rcp26_posneg_',as.character(warm_thresh),'.png',sep=''),plot=p,dpi=300,width=5,height=5)

}
