import numpy as np
from scipy.io import netcdf

def prod_ESMs_cmip5_tempems():

  #Process the 1%/yr data (OPPY) from netcdf to csv format
  nc_file = 'Data/knutti_data.nc'
  df = netcdf.netcdf_file(nc_file,'r')
  mods_num = df.variables['model'].data
  mods_name = df.Models.split(' ')
  mods_dict = {mods_num[i]:mods_name[i] for i in range(0,len(mods_name))}

  years = df.variables['Year'].data
  wanted_vars = ['Temperature_abrupt','Emissions_abrupt','Temperature','Emissions']
  exp_name = ['1pctCO2','abrupt4xCO2']
  #var_name = ['tas','fCO2antt']
  var_name = ['Temperature|anom for piControl','Total anthropogenic carbon flux']
  units = ['K','PgC/yr']
  mv = df.missing_value

  #Build the scenario/variable dictionary
  scen_dict = {'Temperature_abrupt':[exp_name[1],var_name[0],units[0]],
               'Emissions_abrupt':[exp_name[1],var_name[1],units[1]],
               'Temperature':[exp_name[0],var_name[0],units[0]],
               'Emissions':[exp_name[0],var_name[1],units[1]]}

  #Prep for CSV output
  out_data =[]
  for i in range(0,len(mods_name)):
    for j in range(0,len(wanted_vars)):
      nc_data =  np.copy(df.variables[wanted_vars[j]].data[i])
      nc_data[nc_data==mv] = 'nan'
      d_row = [mods_name[i]] + scen_dict[wanted_vars[j]]+['0'] + list(nc_data)
      out_data.append(d_row)
      if j in [1,3]:
        d_row = [mods_name[i]]+[scen_dict[wanted_vars[j]][0]]+['Total cumulative CO2 emissions|since first year','PgC','0'] + list(np.cumsum(nc_data))
        out_data.append(d_row)

  #Load the csv files for the RCPs
  rcps = ['26','45','6','85']
  rcp_years = np.arange(1850,2100)
  dt_rcp = np.dtype({'names':(['Model']+[str(f) for f in rcp_years]),'formats':([np.object] + len(rcp_years)*[np.float])})
  padding = (len(years) - len(rcp_years)) * [np.nan]

  for sen in rcps:
    t_data = np.genfromtxt('Data/AR5WGISPM10_rawmodeldata_FINAL_rcp'+sen+'temp.csv',skiprows=2,delimiter=',',dtype=dt_rcp)
    cum_e_data = np.genfromtxt('Data/AR5WGISPM10_rawmodeldata_FINAL_rcp'+sen+'cumc.csv',skiprows=2,delimiter=',',dtype=dt_rcp)

    #Make the arrays into a useful format
    for i in range(0,t_data.shape[0]):
      mod = t_data['Model'][i]
      if '.dat' in mod:
        mod = mod.strip('_'+mod.split('_')[-1])

      t_row = [mod,'RCP'+sen,'Temperature|rel to 1861-80','K',1850] + list((t_data[i]).tolist()[1:]) +padding
      cumc = list((cum_e_data[i]).tolist()[1:])
      ann_ems = [cumc[x+1]-cumc[x] for x in range(0,len(cumc)-1)]
      e_row = [mod,'RCP'+sen,'Total anthropogenic carbon flux','PgC/yr',1850] + ann_ems + padding
      cumc_row = [mod,'RCP'+sen,'Total cumulative CO2 emissions|since start of 1870','PgC',1850] + cumc + padding

      out_data.append(t_row)
      out_data.append(e_row)
      out_data.append(cumc_row)


  #Write the combined data to a CSV file
  out_file = 'Data/ESM_cmip5_tempems.csv'
  f_out = open(out_file,'w')
  f_out.write('Model,Scenario,Variable,Unit,First_Year,'+','.join([str(f) for f in years]))

  for i in range(0,len(out_data)):
    f_out.write(','.join([str(f) for f in out_data[i]])+'\n')

  f_out.close()

  return

