source("tcre_var_functs.R")

#Load in the data from the CMIP5 ensembles
cmip5_data <- read.csv('./Data/ESM_cmip5_tempems.csv')

esms <- c('Bern3D','DCESS','GENIE','IGSM','UVi')


emics_data <- cmip5_data[cmip5_data$Model %in% emics,]
esms_data <- cmip5_data[!(cmip5_data$Model %in% emics),]

emics_rcp26_data <-  emics_data[(emics_data$Scenario=='RCP26'),-c(6:(2004-1850+6))]
emics_rcp45_data <-  emics_data[(emics_data$Scenario=='RCP45'),-c(6:(2004-1850+6))]
emics_rcp6_data <-  emics_data[(emics_data$Scenario=='RCP6'),-c(6:(2004-1850+6))]
emics_rcp85_data <-  emics_data[(emics_data$Scenario=='RCP85'),-c(6:(2004-1850+6))]
esms_rcp26_data <-  esms_data[(esms_data$Scenario=='RCP26'),-c(6:(2004-1850+6))]
esms_rcp45_data <-  esms_data[(esms_data$Scenario=='RCP45'),-c(6:(2004-1850+6))]
esms_rcp6_data <-  esms_data[(esms_data$Scenario=='RCP6'),-c(6:(2004-1850+6))]
esms_rcp85_data <-  esms_data[(esms_data$Scenario=='RCP85'),-c(6:(2004-1850+6))]

library(ggplot2)
cbPalette <- c("#101FDD", "#0DCAEA", "#DDB510", "#EA0D0D")
warm_threshs <- c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)
for (warm_thresh in warm_threshs) {
    out_rcp26_emics <-calc_budget_dist_allmods(warm_thresh,emics_rcp26_data,'Temperature|rel to 1861-80','Total anthropogenic carbon flux',c(2005:2159))
    out_rcp45_emics <-calc_budget_dist_allmods(warm_thresh,emics_rcp45_data,'Temperature|rel to 1861-80','Total anthropogenic carbon flux',c(2005:2159))
    out_rcp6_emics <-calc_budget_dist_allmods(warm_thresh,emics_rcp6_data,'Temperature|rel to 1861-80','Total anthropogenic carbon flux',c(2005:2159))
    out_rcp85_emics <-calc_budget_dist_allmods(warm_thresh,emics_rcp85_data,'Temperature|rel to 1861-80','Total anthropogenic carbon flux',c(2005:2159))
    out_rcp26_esms <-calc_budget_dist_allmods(warm_thresh,esms_rcp26_data,'Temperature|rel to 1861-80','Total anthropogenic carbon flux',c(2005:2159))
    out_rcp45_esms <-calc_budget_dist_allmods(warm_thresh,esms_rcp45_data,'Temperature|rel to 1861-80','Total anthropogenic carbon flux',c(2005:2159))
    out_rcp6_esms <-calc_budget_dist_allmods(warm_thresh,esms_rcp6_data,'Temperature|rel to 1861-80','Total anthropogenic carbon flux',c(2005:2159))
    out_rcp85_esms <-calc_budget_dist_allmods(warm_thresh,esms_rcp85_data,'Temperature|rel to 1861-80','Total anthropogenic carbon flux',c(2005:2159))



df_all <- data.frame(matrix(nrow=0,ncol=4))
df_all <- rbind(df_all,matrix(c(rep('EMICs',length(out_rcp26_emics[[2]])),rep('RCP2.6',length(out_rcp26_emics[[2]])),out_rcp26_emics[[1]],out_rcp26_emics[[2]]/warm_thresh),ncol=4))
df_all <- rbind(df_all,matrix(c(rep('EMICs',length(out_rcp45_emics[[2]])),rep('RCP4.5',length(out_rcp45_emics[[2]])),out_rcp45_emics[[1]],out_rcp45_emics[[2]]/warm_thresh),ncol=4))
df_all <- rbind(df_all,matrix(c(rep('EMICs',length(out_rcp6_emics[[2]])),rep('RCP6.0',length(out_rcp6_emics[[2]])),out_rcp6_emics[[1]],out_rcp6_emics[[2]]/warm_thresh),ncol=4))
df_all <- rbind(df_all,matrix(c(rep('EMICs',length(out_rcp85_emics[[2]])),rep('RCP8.5',length(out_rcp85_emics[[2]])),out_rcp85_emics[[1]],out_rcp85_emics[[2]]/warm_thresh),ncol=4))
df_all <- rbind(df_all,matrix(c(rep('ESMs',length(out_rcp26_esms[[2]])),rep('RCP2.6',length(out_rcp26_esms[[2]])),out_rcp26_esms[[1]],out_rcp26_esms[[2]]/warm_thresh),ncol=4))
df_all <- rbind(df_all,matrix(c(rep('ESMs',length(out_rcp45_esms[[2]])),rep('RCP4.5',length(out_rcp45_esms[[2]])),out_rcp45_esms[[1]],out_rcp45_esms[[2]]/warm_thresh),ncol=4))
df_all <- rbind(df_all,matrix(c(rep('ESMs',length(out_rcp6_esms[[2]])),rep('RCP6.0',length(out_rcp6_esms[[2]])),out_rcp6_esms[[1]],out_rcp6_esms[[2]]/warm_thresh),ncol=4))
df_all <- rbind(df_all,matrix(c(rep('ESMs',length(out_rcp85_esms[[2]])),rep('RCP8.5',length(out_rcp85_esms[[2]])),out_rcp85_esms[[1]],out_rcp85_esms[[2]]/warm_thresh),ncol=4))


df_all[,3] <- as.numeric(as.vector(df_all[,3]))
df_all[,4] <- as.numeric(as.vector(df_all[,4]))

colnames(df_all) <- c('Ensemble','Scenario','Duration','Budget')




p <- ggplot(df_all, aes(Budget, colour = Scenario,linetype=Ensemble)) + geom_density(alpha=0.1) + labs(x = "Budget (GtC/°C)") + scale_colour_manual(values=cbPalette) + ggtitle(paste('Warming: ',as.character(warm_thresh),'°C',sep=''))
ggsave(paste('Figures/esmsemics_rcps_',as.character(warm_thresh),'.png',sep=''),plot=p,dpi=300,width=5,height=5)

}




