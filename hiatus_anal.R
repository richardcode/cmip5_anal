source("tcre_var_functs.R")

#Load in the data from the CMIP5 ensembles
cmip5_data <- read.csv('./Data/ESM_cmip5_tempems.csv')

cmip5_hist_data <-  cmip5_data[(cmip5_data$Scenario %in% c('RCP26','RCP45','RCP6','RCP85')),c(1:(2004-1850+6))]

cmip5_hist_data_pre <- cmip5_hist_data[,c(1:(1997-1850+6))]
cmip5_hist_data_pos <- cmip5_hist_data[,c(c(1:5),c((1998-1850+6):(2004-1850+6)))]


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



#Cut the observational ensemble pre and post hiatus
gcp_ens_pre <- gcp_ens[,hc_years<1998]
gcp_ens_pos <- gcp_ens[,hc_years>=1998]
hc_ens_pre <- hc_ens[,hc_years<1998]
hc_ens_pos <- hc_ens[,hc_years>=1998]

warm_threshs <- c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)
for (warm_thresh in warm_threshs) {
    out_obs_pre <- calc_budget_dist_obsens(warm_thresh,hc_ens_pre,gcp_ens_pre,hc_years[hc_years<1998])
    out_obs_pos <- calc_budget_dist_obsens(warm_thresh,hc_ens_pos,gcp_ens_pos,hc_years[hc_years>=1998])
    out_cmip5_pre <-calc_budget_dist_allmods(warm_thresh,cmip5_hist_data_pre,'Temperature|rel to 1861-80','Total anthropogenic carbon flux',c(1850:1997))
    out_cmip5_pos <-calc_budget_dist_allmods(warm_thresh,cmip5_hist_data_pos,'Temperature|rel to 1861-80','Total anthropogenic carbon flux',c(1998:2004))
    
df_all <- data.frame(matrix(nrow=0,ncol=4))
df_all <- rbind(df_all,matrix(c(rep('Pre-Hiatus',length(out_obs_pre[[2]])),rep('Observations',length(out_obs_pre[[2]])),out_obs_pre[[1]],out_obs_pre[[2]]/warm_thresh),ncol=4))
df_all <- rbind(df_all,matrix(c(rep('Post-Hiatus',length(out_obs_pos[[2]])),rep('Observations',length(out_obs_pos[[2]])),out_obs_pos[[1]],out_obs_pos[[2]]/warm_thresh),ncol=4))
df_all <- rbind(df_all,matrix(c(rep('Pre-Hiatus',length(out_cmip5_pre[[2]])),rep('CMIP5',length(out_cmip5_pre[[2]])),out_cmip5_pre[[1]],out_cmip5_pre[[2]]/warm_thresh),ncol=4))
df_all <- rbind(df_all,matrix(c(rep('Post-Hiatus',length(out_cmip5_pos[[2]])),rep('CMIP5',length(out_cmip5_pos[[2]])),out_cmip5_pos[[1]],out_cmip5_pos[[2]]/warm_thresh),ncol=4))



df_all[,3] <- as.numeric(as.vector(df_all[,3]))
df_all[,4] <- as.numeric(as.vector(df_all[,4]))

colnames(df_all) <- c('Period','Type','Duration','Budget')

p <- ggplot(df_all, aes(Budget, colour = Period,linetype=Type)) + geom_density(alpha=0.1) + labs(x = "Budget (GtC/°C)") + scale_colour_manual(values=cbPalette) + ggtitle(paste('Warming: ',as.character(warm_thresh),'°C',sep=''))
ggsave(paste('Figures/preposthiatus_',as.character(warm_thresh),'.png',sep=''),plot=p,dpi=300,width=5,height=5)

}















