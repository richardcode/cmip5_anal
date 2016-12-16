#Create the budget distribution calculation function
calc_budget_dist <- function(warm_thresh,temps,emms,years) {
    periods <- c()
    budgets <- c()
    for (i in c(1:length(years))){
        start_temp <- temps[i]
        if (is.na(start_temp)!=TRUE){
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
                
            }
        }
        }
    return(list(periods,budgets))
}

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


#Load in the data from the CMIP5 ensembles
cmip5_data <- read.csv('./Data/ESM_cmip5_tempems.csv')
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

#1 sigma uncertainty in fossil fuel emissions (+/- 5%)
#1 sigma uncertainty in land use emissions (+/-0.5GtC/yr)

#Loop through historical data to find the distribution of budgets for specified warming
warm_thresh <- 0.5

out_hc <- calc_budget_dist(warm_thresh,hc_temps,gcp_data$tot_ems,hc_years)

#Load the HadCRUT-CW data
hccw_data <- read.table('./Data/had4_krig_annual_v2_0_0.txt')
hccw_temps <- hccw_data[,2]
hccw_temps <- hccw_temps - mean(hccw_temps[hc_years>=base_start & hc_years<=base_end])

out_hccw <- calc_budget_dist(warm_thresh,hccw_temps,gcp_data$tot_ems,hc_years)

#Load the NOAA temperature data
noaa_data <- read.csv('./Data/noaatemp_1880-2015.csv',skip=3)
colnames(noaa_data) <- c('Years','Temperature')

out_noaa <- calc_budget_dist(warm_thresh,noaa_data$Temperature,gcp_data$tot_ems[hc_years>=noaa_data$Years[1]],noaa_data$Years)

#Load the NASA temperature data
giss_data <- read.csv('./Data/gisstemp_loti.csv',skip=1,header=TRUE,na.strings="***")

out_giss <- calc_budget_dist(warm_thresh,giss_data$J.D,gcp_data$tot_ems[hc_years>=giss_data$Year[1]],giss_data$Year)

#Loop over the HadCRUT ensemble to get distribution of budgets for same level of warming
out_hc_ens <- list(c(),c())
for (j in c(1:ncol(hc_ens))) {
    out_hc_ens_m <- calc_budget_dist(warm_thresh,hc_ens[,j],gcp_data$tot_ems,hc_years)
    out_hc_ens[[1]] <- c(out_hc_ens[[1]],out_hc_ens_m[[1]])
    out_hc_ens[[2]] <- c(out_hc_ens[[2]],out_hc_ens_m[[2]])
}


#Load the model data
cmip5_data <- read.csv('./Data/ESM_cmip5_tempems.csv',header=TRUE)
pct_data <- cmip5_data[cmip5_data$Scenario='1pctCO2',]
hist_data <-  cmip5_data[(cmip5_data$Scenario %in% c('RCP26','RCP45','RCP6','RCP85')),c(1:(2004-1850+6))]
rcp26_data <-  cmip5_data[(cmip5_data$Scenario=='RCP26'),-c(6:(2004-1850+6))]
rcp45_data <-  cmip5_data[(cmip5_data$Scenario=='RCP45'),-c(6:(2004-1850+6))]
rcp6_data <-  cmip5_data[(cmip5_data$Scenario=='RCP6'),-c(6:(2004-1850+6))]
rcp85_data <-  cmip5_data[(cmip5_data$Scenario=='RCP85'),-c(6:(2004-1850+6))]

out_hist <- list(c(),c())
for (j in c(1:nrow(hist_data[hist_data$Variable=='Temperature|rel to 1861-80',]))){
    out_hist_e <- calc_budget_dist(warm_thresh,as.numeric(hist_data[hist_data$Variable=='Temperature|rel to 1861-80',][j,-c(1:5)]),as.numeric(hist_data[hist_data$Variable=='Total anthropogenic carbon flux',][j,-c(1:5)]),c(1850:2004))
    out_hist[[1]]<- c(out_hist[[1]],out_hist_e[[1]])
    out_hist[[2]]<- c(out_hist[[2]],out_hist_e[[2]])
}

out_rcp26 <- list(c(),c())
for (j in c(1:nrow(rcp26_data[rcp26_data$Variable=='Temperature|rel to 1861-80',]))){
    out_rcp26_e <- calc_budget_dist(warm_thresh,as.numeric(rcp26_data[rcp26_data$Variable=='Temperature|rel to 1861-80',][j,-c(1:5)]),as.numeric(rcp26_data[rcp26_data$Variable=='Total anthropogenic carbon flux',][j,-c(1:5)]),c(2005:2159))
    out_rcp26[[1]]<- c(out_rcp26[[1]],out_rcp26_e[[1]])
    out_rcp26[[2]]<- c(out_rcp26[[2]],out_rcp26_e[[2]])
}

out_rcp45 <- list(c(),c())
for (j in c(1:nrow(rcp45_data[rcp45_data$Variable=='Temperature|rel to 1861-80',]))){
    out_rcp45_e <- calc_budget_dist(warm_thresh,as.numeric(rcp45_data[rcp45_data$Variable=='Temperature|rel to 1861-80',][j,-c(1:5)]),as.numeric(rcp45_data[rcp45_data$Variable=='Total anthropogenic carbon flux',][j,-c(1:5)]),c(2005:2159))
    out_rcp45[[1]]<- c(out_rcp45[[1]],out_rcp45_e[[1]])
    out_rcp45[[2]]<- c(out_rcp45[[2]],out_rcp45_e[[2]])
}

out_rcp6 <- list(c(),c())
for (j in c(1:nrow(rcp6_data[rcp6_data$Variable=='Temperature|rel to 1861-80',]))){
    out_rcp6_e <- calc_budget_dist(warm_thresh,as.numeric(rcp6_data[rcp6_data$Variable=='Temperature|rel to 1861-80',][j,-c(1:5)]),as.numeric(rcp6_data[rcp6_data$Variable=='Total anthropogenic carbon flux',][j,-c(1:5)]),c(2005:2159))
    out_rcp6[[1]]<- c(out_rcp6[[1]],out_rcp6_e[[1]])
    out_rcp6[[2]]<- c(out_rcp6[[2]],out_rcp6_e[[2]])
}

out_rcp85 <- list(c(),c())
for (j in c(1:nrow(rcp85_data[rcp85_data$Variable=='Temperature|rel to 1861-80',]))){
    out_rcp85_e <- calc_budget_dist(warm_thresh,as.numeric(rcp85_data[rcp85_data$Variable=='Temperature|rel to 1861-80',][j,-c(1:5)]),as.numeric(rcp85_data[rcp85_data$Variable=='Total anthropogenic carbon flux',][j,-c(1:5)]),c(2005:2159))
    out_rcp85[[1]]<- c(out_rcp85[[1]],out_rcp85_e[[1]])
    out_rcp85[[2]]<- c(out_rcp85[[2]],out_rcp85_e[[2]])
}


df_all <- data.frame(matrix(nrow=0,ncol=3))
#df_all <- rbind(df_all,matrix(c(rep('rcp26',length(out_rcp26[[2]])),out_rcp26[[1]],out_rcp26[[2]]),ncol=3))
#df_all <- rbind(df_all,matrix(c(rep('rcp45',length(out_rcp45[[2]])),out_rcp45[[1]],out_rcp45[[2]]),ncol=3))
#df_all <- rbind(df_all,matrix(c(rep('rcp6',length(out_rcp6[[2]])),out_rcp6[[1]],out_rcp6[[2]]),ncol=3))
#df_all <- rbind(df_all,matrix(c(rep('rcp85',length(out_rcp85[[2]])),out_rcp85[[1]],out_rcp85[[2]]),ncol=3))
df_all <- rbind(df_all,matrix(c(rep('hist',length(out_hist[[2]])),out_hist[[1]],out_hist[[2]]),ncol=3))
df_all <- rbind(df_all,matrix(c(rep('hc',length(out_hc[[2]])),out_hc[[1]],out_hc[[2]]),ncol=3))
df_all <- rbind(df_all,matrix(c(rep('hccw',length(out_hccw[[2]])),out_hccw[[1]],out_hccw[[2]]),ncol=3))
df_all <- rbind(df_all,matrix(c(rep('hc_ens',length(out_hc_ens[[2]])),out_hc_ens[[1]],out_hc_ens[[2]]),ncol=3))
df_all <- rbind(df_all,matrix(c(rep('giss',length(out_giss[[2]])),out_giss[[1]],out_giss[[2]]),ncol=3))
df_all <- rbind(df_all,matrix(c(rep('noaa',length(out_noaa[[2]])),out_noaa[[1]],out_noaa[[2]]),ncol=3))
df_all[,2] <- as.numeric(as.vector(df_all[,2]))
df_all[,3] <- as.numeric(as.vector(df_all[,3]))

colnames(df_all) <- c('dist','duration','budget')

library(ggplot2)
ggplot(df_all, aes(budget, colour = dist)) + geom_density(alpha = 0.0)











