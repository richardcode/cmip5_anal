import numpy as np
from numpy.random import normal
import sys
sys.path.insert(0,'/Users/rm604/Documents/FAIR/')
from fair_main_rev import fair_scm
from pandas import DataFrame
from statsmodels.api import OLS
from scipy.odr import *
import statsmodels.tools.tools


def tls_lin_mod(B,x):

  return B[0]*x[0] + B[1]*x[1] +B[2]*x[2]






#forc_file = './Data/ipcc_tod_rftimeseries_extended.txt'
forc_file = './Data/Annualforcings_Mar2014_GHGrevised.txt'

data = np.genfromtxt(forc_file,skip_header=4)
#Remove 2016
data = data[:-1,:]

co2_rf_be = data[:,1]
years = data[:,0]


#Append the extra data from the NOAA (https://www.esrl.noaa.gov/gmd/aggi/aggi.html)
#years = np.append(years,[2012.,2013.,2014.,2015.])
#co2_rf_be = np.append(co2_rf_be,[1.845,1.882,1.908,1.939])

#Create 200 member ensemble of CO2 ERF
co2_rf_2011 = np.array([0.8*co2_rf_be[years==2011],co2_rf_be[years==2011],1.2*co2_rf_be[years==2011]])
co2_rf_sig = (co2_rf_2011[2] - co2_rf_2011[0]) / (2*1.654)

sf_dist = normal(loc=co2_rf_2011[1], scale=co2_rf_sig, size=200) / co2_rf_2011[1]
co2_rf_dist = sf_dist[:,np.newaxis] * co2_rf_be[np.newaxis,:]

#Create 200 member ensemble of total WMGHG
wmghg_rf_be = data[:,1] + data[:,2]
sig_comb = np.sqrt(((2.16 - 1.44)/ (2*1.654))**2 + ((1.2 - 0.8) / (2*1.654))**2)
sf_dist = normal(loc=wmghg_rf_be[years==2011], scale=sig_comb, size=200) / wmghg_rf_be[years==2011]
wmghg_rf_dist = sf_dist[:,np.newaxis] * wmghg_rf_be[np.newaxis,:]



#Load the distribution of rf_anthro
rf_anthro_file = '/Users/rm604/Documents/TCRE1p5/cmip5_anal/Data/rf_anthro_1750_2015.csv'
rf_anthro_e = np.transpose(np.genfromtxt(rf_anthro_file,delimiter=','))

rf_non_wmghg = rf_anthro_e - wmghg_rf_dist

#Load the natural forcing
rf_nat_file = '/Users/rm604/Documents/TCRE1p5/cmip5_anal/Data/rf_nat_1750_2015.csv'
rf_nat_e = np.transpose(np.genfromtxt(rf_nat_file,delimiter=','))

#Load the HadCRUT4 data
d_file = '/Users/rm604/Data/HadCRUT4/global/HadCRUT.4.5.0.0.annual_ns_avg.txt'
hc_data = np.genfromtxt(d_file)
hc_years = (hc_data[:,0])
hc_data = (hc_data[:,1] - np.mean(hc_data[:,1][np.logical_and(hc_years>1860,hc_years<=1880)]))[:-2]
hc_years = hc_years[:-2]


f2x = 3.74 + 0.2*3.74*(sf_dist*co2_rf_2011[1] - co2_rf_2011[1]) / 1.65
#Setup default input parameters and tile to the number of ensemble members and put in f_2x values
def_params = np.array([1.75,2.5,4.1,239.0,0.2173,0.2240,0.2824,0.2763,1000000,394.4,36.54,4.304,100.0,35.0,0.02,4.5,3.74,278.0,2.123,97.0])
e_params = np.tile(def_params,(rf_anthro_e.shape[0],1))
e_params[:,16] = f2x

#t_co2_d = fair_scm(co2_rf_be,mode='forcing_driven')
#t_nco2_d = fair_scm(rf_non_co2[5],mode='forcing_driven')
#t_nat_d = fair_scm(rf_nat_e[0],mode='forcing_driven')

t_co2 = fair_scm(co2_rf_dist,mode='forcing_driven',input_params=e_params,comb='forcparams')
t_co2 = t_co2 - (np.mean(t_co2[:,np.logical_and(years>1860,years<=1880)],axis=-1))[:,np.newaxis]
t_wmghg = fair_scm(wmghg_rf_dist,mode='forcing_driven',input_params=e_params,comb='forcparams')
t_wmghg = t_wmghg - (np.mean(t_wmghg[:,np.logical_and(years>1860,years<=1880)],axis=-1))[:,np.newaxis]
t_nwmghg = fair_scm(rf_non_wmghg,mode='forcing_driven',input_params=e_params,comb='forcparams')
t_nwmghg = t_nwmghg - (np.mean(t_nwmghg[:,np.logical_and(years>1860,years<=1880)],axis=-1))[:,np.newaxis]
t_nat = fair_scm(rf_nat_e,mode='forcing_driven',input_params=e_params,comb='forcparams')
t_nat = t_nat - (np.mean(t_nat[:,np.logical_and(years>1860,years<=1880)],axis=-1))[:,np.newaxis]
t_anthro = fair_scm(rf_anthro_e,mode='forcing_driven',input_params=e_params,comb='forcparams')
t_anthro = t_anthro - (np.mean(t_anthro[:,np.logical_and(years>1860,years<=1880)],axis=-1))[:,np.newaxis]


#Do the attributable warming regression
w_co2=[]
w_nwmghg=[]
w_wmghg=[]
w_nat=[]
for i in range(0,t_co2.shape[0]):
  #TLS regression
  #x = np.array([ (t_wmghg)[i,years>=hc_years[0]],(t_nwmghg)[i,years>=hc_years[0]],t_nat[i,years>=hc_years[0]] ])
  #y = hc_data
  #linear = Model(tls_lin_mod)
  #mydata = Data(x, y)
  #myodr = ODR(mydata, linear, beta0=[1., 1.,1.])
  #myoutput = myodr.run()
  #w_co2.append(myoutput.beta[0]*t_co2[i])
  #w_nwmghg.append(myoutput.beta[1]*t_nwmghg[i])
  #w_wmghg.append(myoutput.beta[0]*t_wmghg[i])
  #w_nat.append(myoutput.beta[2]*t_nat[i])


  #OLS regression
  y = hc_data
  x = DataFrame({'x1': (t_wmghg)[i,years>=hc_years[0]], 'x2': (t_nwmghg)[i,years>=hc_years[0]], 'x3':t_nat[i,years>=hc_years[0]]})
  x = statsmodels.tools.tools.add_constant(x)
  model = OLS(y, x)
  result = model.fit()
  sf = result.params
  w_co2.append(sf['x1']*t_co2[i])
  w_nwmghg.append(sf['x2']*t_nwmghg[i])
  w_wmghg.append(sf['x1']*t_wmghg[i])
  w_nat.append(sf['x3']*t_nat[i])


w_co2=np.array(w_co2)
w_wmghg=np.array(w_wmghg)
w_nwmghg = np.array(w_nwmghg)
w_nat=np.array(w_nat)


w_co2_2015 = w_co2[:,years==2015]

#Write out the attributable CO2-induced warming timeseries
f_out_name = './Data/attribwarm_co2_hadcrut4.txt'
f_out=open(f_out_name,'w')
for i in range(0,w_co2.shape[0]):
  f_out.write(','.join([str(f) for f in w_co2[i].tolist()])+'\n')
f_out.close()










