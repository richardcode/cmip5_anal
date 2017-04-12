source("tcre_var_functs.R")

#Load in the data from the CMIP5 ensembles
cmip5_data <- read.csv('./Data/ESM_cmip5_tempems.csv')
#Get models that are both in RCP2.6 and 4xCO2

step_data <- cmip5_data[cmip5_data$Scenario=='abrupt4xCO2',]
rcp26_data <- cmip5_data[cmip5_data$Scenario=='RCP26',-c(6:(2004-1850+6))]

#Get overlapping set of models
com_mod <- intersect(unique(step_data$Model),unique(rcp26_data$Model))
step_data <- step_data[step_data$Model %in% com_mod,-c(6,7)]
rcp26_data <- rcp26_data[rcp26_data$Model %in% com_mod,]

#Select on emissions
plot(c(1:308),step_data[step_data$Variable=='Total anthropogenic carbon flux',-c(1:5)][1,],
col = NA,
#     type = "p",
#     xaxt = "n", yaxt = "n",
    ylim = c(-1,2),
xlab = expression("Year after doubling"),
ylab = 'Emissions (GtC/yr)',
cex = 2.0,
)
grid()

for (i in c(1:nrow(step_data[step_data$Variable=='Total anthropogenic carbon flux',-c(1:5)]))) {
    umid_data <- step_data[step_data$Variable=='Total anthropogenic carbon flux',-c(1:5)][i,]
    lines(c(1:308),umid_data[1,],col='grey',lwd=0.4)
}




library(ggplot2)
warm_threshs <- c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)
for (warm_thresh in warm_threshs) {
    out_rcp26 <-calc_budget_dist_allmods(warm_thresh,rcp26_data,'Temperature|rel to 1861-80','Total anthropogenic carbon flux',c(2005:2159))
    out_step <-calc_budget_dist_allmods(warm_thresh,step_data,'Temperature|anom for piControl','Total anthropogenic carbon flux',c(1:308))

df_all <- data.frame(matrix(nrow=0,ncol=4))
df_all <- rbind(df_all,matrix(c(rep('Common',length(out_rcp26[[2]])),rep('RCP2.6',length(out_rcp26[[2]])),out_rcp26[[1]],out_rcp26[[2]]/warm_thresh),ncol=4))
df_all <- rbind(df_all,matrix(c(rep('Common',length(out_step[[2]])),rep('abpupt4xCO2',length(out_step[[2]])),out_step[[1]],out_step[[2]]/warm_thresh),ncol=4))


df_all[,3] <- as.numeric(as.vector(df_all[,3]))
df_all[,4] <- as.numeric(as.vector(df_all[,4]))

colnames(df_all) <- c('Ensemble','Scenario','Duration','Budget')

p <- ggplot(df_all, aes(Budget, colour = Scenario)) + geom_density(alpha=0.1) + labs(x = "Budget (GtC/°C)")  + ggtitle(paste('Warming: ',as.character(warm_thresh),'°C',sep=''))