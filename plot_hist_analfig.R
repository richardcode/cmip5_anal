base_start <- 1861
base_end <- 1880

#Load the observational data
hc_data <- read.table('/Users/rm604/Data/HadCRUT4/global/HadCRUT.4.5.0.0.annual_ns_avg.txt')
hc_years <- hc_data[,1]
hc_temps <- hc_data[,2]

hc_temps <- hc_temps - mean(hc_temps[hc_years>=base_start & hc_years<=base_end])
#Remove incomplete 2016 temp data
hc_temps <- hc_temps[-c(length(hc_years),length(hc_years)-1)]
hc_years <- hc_years[-c(length(hc_years),length(hc_years)-1)]
gcp_data <- read.csv('./Data/gcbdata_2016.csv')
gcp_data$tot_ems <- gcp_data$fossil.fuel.and.cement.emissions + gcp_data$land.use.change.emissions
#Cut the emissions data to be after 1850 only
gcp_data<- gcp_data[gcp_data$Year>=1850,]


#Get the HadCRUT ensemble
hadcrut_dir <- '/Users/rm604/Data/HadCRUT4/global/ensemble/'
hc_files <- list.files(hadcrut_dir)
hc_ens <- matrix(nrow=166,ncol=100)
hc_stds <- matrix(nrow=166,ncol=100)
for (f in c(1:100)) {
    data <- read.table(paste(hadcrut_dir,hc_files[f],sep=''))
    data <- data[-c(length(hc_years),length(hc_years)-1),]
    hc_ens[,f] <- data[,2]
    hc_stds[,f] <- data[,3]
}
hc_ens <- t(hc_ens)
base_mean <- apply(hc_ens[,hc_years>=base_start & hc_years<=base_end],1,mean)
hc_ens <- sweep(hc_ens,1,base_mean,"-")

#Load in the data from the CMIP5 ensembles
cmip5_data <- read.csv('./Data/ESM_cmip5_tempems.csv')

emics <- c('Bern3D','DCESS','GENIE','IGSM','UVi')


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

gcp_data_c <- cumsum(gcp_data$tot_ems)
gcp_data_c <- gcp_data_c - gcp_data_c[gcp_data$Year==1869]
gcp_ens_c <- t(apply(gcp_ens, 1, cumsum))
gcp_ens_c <- sweep(gcp_ens_c,1,gcp_ens_c[,gcp_data$Year==1869],"-")


#Load the HadCRUT-CW data
hccw_data <- read.table('./Data/had4_krig_annual_v2_0_0.txt')
hccw_temps <- hccw_data[,2]
hccw_temps <- hccw_temps - mean(hccw_temps[hc_years>=base_start & hc_years<=base_end])

#Load the NOAA temperature data
noaa_data <- read.csv('./Data/noaatemp_1880-2015.csv',skip=3)
colnames(noaa_data) <- c('Years','Temperature')
noaa_data$Temperature <- noaa_data$Temperature -   mean(noaa_data$Temperature[noaa_data$Years>=base_start & noaa_data$Years<=base_end])

#Load the NASA temperature data
giss_data <- read.csv('./Data/gisstemp_loti.csv',skip=1,header=TRUE,na.strings="***")
giss_data$J.D <- giss_data$J.D - mean(giss_data$J.D[giss_data$Year>=base_start & giss_data$Year<=base_end])

#Load the attributable warming data
awi_data <- read.csv('./Data/AWI_Tmean_Quant.csv',header=TRUE)
#Mean over the months in a year
awi_data$DATE <- floor(awi_data$DATE)
awi_data_ann <- data.frame(matrix(nrow=0,ncol=ncol(awi_data)))
for (i in unique(awi_data$DATE)){
    awi_data_ann <- rbind(awi_data_ann,colMeans(awi_data[awi_data$DATE==i,]))
}
colnames(awi_data_ann) <- colnames(awi_data)
awi_data <- awi_data_ann[c(-nrow(awi_data_ann)),]
awi_data$T_ANT_FIT <- awi_data$T_ANT_FIT - mean(awi_data$T_ANT_FIT[awi_data$DATE>=base_start & awi_data$DATE<=base_end])

#Load the attributable warming data with CW masking technique
awicw_data <- read.csv('./Data/AWI_CW_Tmean_Quant.csv',header=TRUE)
#Mean over the months in a year
awicw_data$DATE <- floor(awicw_data$DATE)
awicw_data_ann <- data.frame(matrix(nrow=0,ncol=ncol(awicw_data)))
for (i in unique(awicw_data$DATE)){
    awicw_data_ann <- rbind(awicw_data_ann,colMeans(awicw_data[awicw_data$DATE==i,]))
}
colnames(awicw_data_ann) <- colnames(awicw_data)
awicw_data <- awicw_data_ann[c(-nrow(awicw_data_ann)),]
awicw_data$T_ANT_FIT <- awicw_data$T_ANT_FIT - mean(awicw_data$T_ANT_FIT[awicw_data$DATE>=base_start & awicw_data$DATE<=base_end])


#Do estimates of TCRE using different period
periods_start <- c(1970,1980,1990,2000,2006)
periods_end <- c(1979,1989,1999,2009,2015)
periods_ave <- (periods_start + periods_end)*0.5

period_tcre_out <- matrix(nrow=length(periods_start),ncol=4)
period_tcre_ens <- matrix(nrow=length(periods_start),ncol=4)
for (i in c(1:length(periods_start))){
    
    c_ems <- mean(gcp_data_c[hc_years>=periods_start[i] & hc_years<=periods_end[i]])
    
    tcre_hc <- mean(hc_temps[hc_years>=periods_start[i] & hc_years<=periods_end[i]]) / c_ems*1000.
    tcre_cw <- mean(hccw_temps[hc_years>=periods_start[i] & hc_years<=periods_end[i]]) / c_ems*1000.
    tcre_awi <- mean(awi_data$T_ANT_FIT[awi_data$DATE>=periods_start[i] & awi_data$DATE<=periods_end[i]]) / c_ems*1000.
    tcre_awicw <- mean(awicw_data$T_ANT_FIT[awicw_data$DATE>=periods_start[i] & awicw_data$DATE<=periods_end[i]]) / c_ems *1000.
    
    period_tcre_out[i,] <- c(tcre_hc,tcre_cw,tcre_awi,tcre_awicw)
    
    #Calculate percentiles from across the ensemble
    c_ems_e <- apply(gcp_ens_c[,hc_years>=periods_start[i] & hc_years<=periods_end[i]],1,mean)
    tcre_hc_e <- apply(hc_ens[,hc_years>=periods_start[i] & hc_years<=periods_end[i]],1,mean) / c_ems_e * 1000.0
    
    period_tcre_ens[i,] <- c(quantile(tcre_hc_e,0.05),quantile(tcre_hc_e,0.5),quantile(tcre_hc_e,0.95),mean(tcre_hc_e))
}


#Make plot of the historical temperatures and cumulative emissions
library(RColorBrewer)
cols <- brewer.pal(9,"Set1")

#png(file="./Data/hist_fig.png",width=600,height=500,res=200)
pdf(file="./Figures/hist_fig.pdf")
op <- par(mar=c(5.0, 6.0, 2.0, 3.0) + 0.1)
par(mfrow=c(2,2))

xmin <- -30.
xmax = 700.0
ymax = 1.2
ymin = -0.5
plot(1,type='n',xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab=expression('Cumulative CO'[2]*' emissions since 1870 (GtC)'), ylab=expression('Warming relative to 1861-1880 (K)'),xaxs='i',yaxs='i')

for (i in c(1:nrow(hc_ens))){
    par(new=TRUE)
    par(xpd=FALSE)
    lines(gcp_ens_c[i,],hc_ens[i,],type='l',col='grey',lwd=0.3)

}
par(new=TRUE)
par(xpd=FALSE)
lines(gcp_data_c,hc_temps,type='l',col='black',lwd=1.0)
par(new=TRUE)
par(xpd=FALSE)
lines(gcp_data_c,hccw_temps,type='l',col=cols[1],lwd=1.0)
#par(new=TRUE)
#par(xpd=FALSE)
#lines(gcp_data_c[hc_years>=noaa_data$Years[1]],noaa_data$Temperature,type='l',col=cols[2],lwd=1.0)
#par(new=TRUE)
#par(xpd=FALSE)
#lines(gcp_data_c[hc_years>=giss_data$Year[1]],giss_data$J.D,type='l',col=cols[3],lwd=1.0)
par(new=TRUE)
par(xpd=FALSE)
lines(gcp_data_c[hc_years>=awi_data$DATE[1]],awi_data$T_ANT_FIT,type='l',col=cols[4],lwd=1.0)
par(new=TRUE)
par(xpd=FALSE)
lines(gcp_data_c[hc_years>=awicw_data$DATE[1]],awicw_data$T_ANT_FIT,type='l',col=cols[2],lwd=1.0)
grid()
legend('bottomright',legend=c('HadCRUT4','HadCRUT4-CW','AWI-HadCRUT','AWI-HadCRUT4-CW','HadCRUT4/GCP ensemble'),col=c('black',cols[1],cols[4],cols[2],'grey'),lty=1,ncol=1,cex=0.65)
mtext('a)',at=-300,side=1,cex=1.0)


#Plot the TCRE ranges
xmin <- 1970
xmax = 2015
ymax = 2.2
ymin = 0.5
plot(1,type='n',xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab='',ylab=expression('All-forcing estimated TCRE (K)'),xaxs='i',yaxs='i')
for (i in c(1:length(periods_ave))){
    lines(rep(periods_ave[i],2),c(period_tcre_ens[i,1],period_tcre_ens[i,3]),lwd=3,col='grey')
}
points(periods_ave,period_tcre_out[,1],col='black',pch=18)
points(periods_ave,period_tcre_out[,2],col=cols[1],pch=18)
points(periods_ave,period_tcre_out[,3],col=cols[4],pch=18)
points(periods_ave,period_tcre_out[,4],col=cols[2],pch=18)
grid()
text(1976,1.0,'1970-79',cex=1.0,srt=90)
text(1986,1.3,'1980-89',cex=1.0,srt=90)
text(1996,1.45,'1990-99',cex=1.0,srt=90)
text(2006,1.6,'2000-09',cex=1.0,srt=90)
text(2012,1.6,'2006-15',cex=1.0,srt=90)

mtext('b)',at=1950,side=1,cex=1.0)

#Load the CO2 attributable warming ensemble
attco2_data <- read.csv('./Data/attribwarm_co2_hadcrut4.txt',header=FALSE)
attwarm_years <- c(1750:2015)
attco2_med <- as.vector(apply(attco2_data,2,quantile,0.5))
attco2_max <- as.vector(apply(attco2_data,2,quantile,0.95))
attco2_min <- as.vector(apply(attco2_data,2,quantile,0.05))

xmin <- -30.
xmax = 700.0
ymax = 1.0
ymin = -0.1
plot(1,type='n',xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab=expression('Cumulative CO'[2]*' emissions since 1870 (GtC)'), ylab='',#ylab=expression('CO'[2]*" -induced \n attributable warming relative to 1861-1880 (K)"),
xaxs='i',yaxs='i')
polygon(c(gcp_data_c,rev(gcp_data_c)),c(attco2_max[attwarm_years>=hc_years[1]],rev(attco2_min[attwarm_years>=hc_years[1]])),col='grey',border=NA)
lines(gcp_data_c,attco2_med[attwarm_years>=hc_years[1]],col='black',lwd=1.0)
grid()
mtext(side = 2, text = c(expression('CO'[2]*'-induced attrbutable warming'), "relative to 1861-1880 (K)"),
line = c(2.7, 1.9),cex=0.8)
mtext('c)',at=-300,side=1,cex=1.0)

#Calculate the period TCREs for attributable CO2 warming
period_tcre_co2 <- matrix(nrow=length(periods_start),ncol=4)
period_tcre_co2_m <- c()
for (i in c(1:length(periods_start))){
    c_ems <- mean(gcp_data_c[hc_years>=periods_start[i] & hc_years<=periods_end[i]])
    tcre_mean <- mean(attco2_med[attwarm_years>=periods_start[i] & attwarm_years<=periods_end[i]]) / c_ems*1000.
    period_tcre_co2_m <- c(period_tcre_co2_m,tcre_mean)
    coll <- c()
    for (j in c(1:nrow(attco2_data))){
        tcre_e <- mean(as.numeric(attco2_data[j,attwarm_years>=periods_start[i] & attwarm_years<=periods_end[i]])) / mean(gcp_ens_c[sample(1:100,1,replace=T),hc_years>=periods_start[i] & hc_years<=periods_end[i]])*1000.
        coll <- c(coll,tcre_e)
    }
    period_tcre_co2[i,] <- c(quantile(coll,0.05),quantile(coll,0.5),quantile(coll,0.95),mean(coll))
    
}

#Plot the Co2 attributable warming TCRE ranges
xmin <- 1970
xmax = 2015
ymax = 2.0
ymin = 0.5
plot(1,type='n',xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab='',ylab=expression('CO'[2]*'-only estimated TCRE (K)'),xaxs='i',yaxs='i')
for (i in c(1:length(periods_ave))){
    lines(rep(periods_ave[i],2),c(period_tcre_co2[i,1],period_tcre_co2[i,3]),lwd=3,col='grey')
}
points(periods_ave,period_tcre_co2_m,col='black',pch=18)
text(1976,1.0,'1970-79',cex=1.0,srt=90)
text(1986,1.0,'1980-89',cex=1.0,srt=90)
text(1996,1.0,'1990-99',cex=1.0,srt=90)
text(2006,1.0,'2000-09',cex=1.0,srt=90)
text(2012,1.0,'2006-15',cex=1.0,srt=90)

mtext('d)',at=1950,side=1,cex=1.0)

grid()
dev.off()





