source("tcre_var_functs.R")

#Load in the data from the CMIP5 ensembles
cmip5_data <- read.csv('./Data/ESM_cmip5_tempems.csv')

hist_data <-  cmip5_data[(cmip5_data$Scenario %in% c('RCP85')),c(1:(2015-1850+6))]
opc_data <-  cmip5_data[(cmip5_data$Scenario %in% c('1pctCO2')),]

hmods <- unique(hist_data$Model)
omods <- unique(opc_data$Model)

u_mods <- intersect(omods,hmods)

hist_data <- hist_data[hist_data$Model %in% u_mods,]
opc_data <- opc_data[opc_data$Model %in% u_mods,]

co2_warm <- data.frame(matrix(nrow=0,ncol=ncol(hist_data)))
#Calculate the fraction of historical warming attributable to non-CO2 forcing at any point
for (i in c(1:length(unique(hist_data$Model)))){
    co2_warm_row <- hist_data[(hist_data$Model==u_mods[i] & hist_data$Variable=='Temperature|rel to 1861-80'),]
    #co2_warm_row[,-c(1:5)] <- rollapply(as.numeric(co2_warm_row[,-c(1:5)]),9,mean,fill=NA,align="center")
    opc_data_i <- opc_data[opc_data$Model==u_mods[i],]
    opc_data_i[opc_data_i$Variable=='Temperature|anom from piControl',-c(1:5)] <- rollapply(as.numeric(opc_data_i[opc_data_i$Variable=='Temperature|anom from piControl',-c(1:5)]),9,mean,fill=NA,align="center")
    for (j in c(6:(2015-1850+6))){
        #Get the linearly interpolated 1pct warming at the current cumulative emissions
        year_j <- j-6+1850
        if (year_j>1850) {
            cum_e_j <- hist_data[(hist_data$Model==u_mods[i] & hist_data$Variable=='Total cumulative CO2 emissions|since start of 1870'),j]
            t_i <- which.max(opc_data_i[opc_data_i$Variable=='Total cumulative CO2 emissions',-c(1:5)]>cum_e_j)
            res <- try(approx(opc_data_i[opc_data_i$Variable=='Total cumulative CO2 emissions',c(6:(t_i+5))],opc_data_i[opc_data_i$Variable=='Temperature|anom from piControl',c(6:(t_i+5))],cum_e_j))
            if (class(res) != "try-error"){
                cont <- approx(opc_data_i[opc_data_i$Variable=='Total cumulative CO2 emissions',c(6:(t_i+5))],opc_data_i[opc_data_i$Variable=='Temperature|anom from piControl',c(6:(t_i+5))],cum_e_j)$y
                co2_warm_row[,j] <-   cont
            }
            else {
                co2_warm_row[,j]<-NA
            }
        } else {
            co2_warm_row[,j]<-NA
        }
        }
    co2_warm <- rbind(co2_warm,co2_warm_row)
        
    }


#Make plot of the historical temperatures and cumulative emissions
library(RColorBrewer)
cols <- brewer.pal(length(unique(hist_data$Model)),"Paired")

xmin <- -30.
xmax = 650.0
ymax = 1.7
ymin = -1.
plot(1,type='n',xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab=expression('Cumulative CO'[2]*' emissions (GtC)'), ylab=expression('Warming relative to 1861-1870 (K)'),xaxs='i',yaxs='i')
for (i in c(1:length(unique(hist_data$Model)))){
    
    temps <- hist_data[hist_data$Model==u_mods[i] & hist_data$Variable=='Temperature|rel to 1861-80',-c(1:5)]
    cum_emms <- hist_data[hist_data$Model==u_mods[i] & hist_data$Variable=='Total cumulative CO2 emissions|since start of 1870',-c(1:5)]
    co2_temps <- co2_warm[co2_warm$Model==u_mods[i] & co2_warm$Variable=='Temperature|rel to 1861-80',-c(1:5)]
    #temps <- opc_data[opc_data$Model==u_mods[i] & opc_data$Variable=='Temperature|anom from piControl',-c(1:5)]
    #cum_emms <- opc_data[opc_data$Model==u_mods[i] & opc_data$Variable=='Total cumulative CO2 emissions',-c(1:5)]
    
    par(new=TRUE)
    par(xpd=FALSE)
    lines(cum_emms,temps,type='l',col=cols[i],lwd=0.8)
    par(new=TRUE)
    par(xpd=FALSE)
    lines(cum_emms,co2_temps,type='l',col=cols[i],lwd=0.8,lty=2)
}
legend('bottomright',legend=unique(hist_data$Model),col=cols,lty=1,ncol=2,cex=0.7)
grid()

#Make distribution of non-CO2 warming in 2015
distrib <- hist_data[hist_data$Variable=='Temperature|rel to 1861-80',(2015-1850+6)] - co2_warm[,(2015-1850+6)]





