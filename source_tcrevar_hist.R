source("tcre_var_functs.R")

#Load in the data from the CMIP5 ensembles
cmip5_data <- read.csv('./Data/ESM_cmip5_tempems.csv')
#Get models that are both in RCP2.6 and 4xCO2


hist_data <-  cmip5_data[(cmip5_data$Scenario %in% c('RCP26','RCP45','RCP6','RCP85')),c(1:(2004-1850+6))]


#Get the HadCRUT ensemble
hadcrut_dir <- '/Users/richardmillar/Documents/Thesis_Mac_work/HadCRUT/'
hc_files <- list.files(hadcrut_dir)
hc_ens <- matrix(nrow=166,ncol=100)
hc_stds <- matrix(nrow=166,ncol=100)
for (f in c(1:100)) {
    data <- read.table(paste(hadcrut_dir,hc_files[f],sep=''))
    hc_ens[,f] <- data[,2]
    hc_stds[,f] <- data[,3]
}

hc_ens <- t(hc_ens)

#Load the observational data
hc_data <- read.table('./Data/HadCRUT.4.4.0.0.annual_ns_avg.txt')
hc_years <- hc_data[,1]
hc_temps <- hc_data[,2]
base_start <- 1861
base_end <- 1880
hc_temps <- hc_temps - mean(hc_temps[hc_years>=base_start & hc_years<=base_end])
#Remove incomplete 2016 temp data
hc_temps <- hc_temps[-length(hc_years)]
hc_years <- hc_years[-length(hc_years)]
gcp_data <- read.csv('./Data/gcbdata_2016.csv')
gcp_data$tot_ems <- gcp_data$fossil.fuel.and.cement.emissions + gcp_data$land.use.change.emissions
#Cut the emissions data to be after 1850 only
gcp_data<- gcp_data[gcp_data$Year>=1850,]

n_samp <- 100
#1 sigma uncertainty in fossil fuel emissions (+/- 5%)
ffi_ems_sf <- rnorm(n_samp,mean=1.0,sd=0.05)
#1 sigma uncertainty in land use emissions (+/-0.5GtC/yr 1959-2015, +/-33% pre 1959)
lu_ems_sf <- rnorm(n_samp,mean=0.0,sd=1.0)

ffi_ens <- outer(ffi_ems_sf,gcp_data$fossil.fuel.and.cement.emissions)
luc_ens <- mat.or.vec(n_samp, ncol(ffi_ens))
luc_ens[,gcp_data$Year<1959] <- outer(1.0+0.33*lu_ems_sf,gcp_data$land.use.change.emissions[gcp_data$Year<1959])
luc_ens[,gcp_data$Year>=1959] <- matrix(rep(0.5*lu_ems_sf,times=sum(gcp_data$Year>=1959)),nrow=n_samp) + matrix(rep(gcp_data$land.use.change.emissions[gcp_data$Year>=1959],n_samp),nrow=n_samp,byrow=TRUE)

gcp_ens <- ffi_ens + luc_ens


warm_threshs <- c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)
for (warm_thresh in warm_threshs) {

out_hist <-calc_budget_dist_allmods(warm_thresh,hist_data,'Temperature|rel to 1861-80','Total anthropogenic carbon flux',c(1850:2004))
out_obs <- calc_budget_dist_obsens(warm_thresh,hc_ens[,hc_years<=2004],gcp_ens[,hc_years<=2004],hc_years[hc_years<=2004])

#Replace the CMIP5 emissions by the GCP emissions
out_hist_obsems <- list(c(),c())
for (j in c(1:nrow(hist_data[hist_data$Variable=='Temperature|rel to 1861-80',]))){
    out_hist_e <- calc_budget_dist(warm_thresh,as.numeric(hist_data[hist_data$Variable=='Temperature|rel to 1861-80',][j,-c(1:5)]),gcp_ens[j,hc_years<=2004],c(1850:2004))
    out_hist_obsems[[1]]<- c(out_hist_obsems[[1]],out_hist_e[[1]])
    out_hist_obsems[[2]]<- c(out_hist_obsems[[2]],out_hist_e[[2]])
}

#Replace the CMIP5 temps by the HC temps
out_hist_obstemp <- list(c(),c())
for (j in c(1:nrow(hist_data[hist_data$Variable=='Total anthropogenic carbon flux',]))){
    out_hist_e <- calc_budget_dist(warm_thresh,hc_ens[j,hc_years<=2004],as.numeric(hist_data[hist_data$Variable=='Total anthropogenic carbon flux',][j,-c(1:5)]),c(1850:2004))
    out_hist_obstemp[[1]]<- c(out_hist_obstemp[[1]],out_hist_e[[1]])
    out_hist_obstemp[[2]]<- c(out_hist_obstemp[[2]],out_hist_e[[2]])
}



df_all <- data.frame(matrix(nrow=0,ncol=3))
df_all <- rbind(df_all,matrix(c(rep('CMIP5',length(out_hist[[2]])),out_hist[[1]],out_hist[[2]]/warm_thresh),ncol=3))
df_all <- rbind(df_all,matrix(c(rep('CMIP5-HadCRUT',length(out_hist_obstemp[[2]])),out_hist_obstemp[[1]],out_hist_obstemp[[2]]/warm_thresh),ncol=3))
df_all <- rbind(df_all,matrix(c(rep('CMIP5-GCP',length(out_hist_obsems[[2]])),out_hist_obsems[[1]],out_hist_obsems[[2]]/warm_thresh),ncol=3))
df_all <- rbind(df_all,matrix(c(rep('Observations',length(out_obs[[2]])),out_obs[[1]],out_obs[[2]]/warm_thresh),ncol=3))




df_all[,2] <- as.numeric(as.vector(df_all[,2]))
df_all[,3] <- as.numeric(as.vector(df_all[,3]))

colnames(df_all) <- c('Combination','Duration','Budget')

p <- ggplot(df_all, aes(Budget, colour = Combination)) + geom_density(alpha=0.1) + labs(x = "Budget (GtC/°C)")  + ggtitle(paste('Warming: ',as.character(warm_thresh),'°C',sep='')) + geom_hline(yintercept=0, colour="white", size=1)+geom_hline(yintercept=0, colour="brown", size=2,linetype='dashed')
ggsave(paste('Figures/varsource_',as.character(warm_thresh),'.png',sep=''),plot=p,dpi=300,width=5,height=5)

}
























