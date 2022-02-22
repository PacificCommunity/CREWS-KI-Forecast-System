# -*- coding: utf-8 -*-
"""
Created on Thu May  6 17:47:05 2021

@author: judithg
"""
import os
import datetime as dt
import netCDF4 as nc
import numpy as np
import numpy.matlib


def download_CNEMS(now):  
    print('Downloading Sea Level Anomaly data from http://nrt.cmems-du.eu')
    folder_tmp ='../tmp/' 
    now = now - dt.timedelta(2)# NCEP needs at least 3 hours to upload their forecast from UTC00
    # now = now.replace(hour=0,minute=0,second=0,microsecond=0)
    then = now+dt.timedelta(days=9.5,hours=1)
    
    
    os.system('python -m motuclient --motu http://nrt.cmems-du.eu/motu-web/Motu\
                              --service-id GLOBAL_ANALYSIS_FORECAST_PHY_001_024-TDS --product-id global-analysis-forecast-phy-001-024-hourly-t-u-v-ssh\
                              --longitude-min 172.97 --longitude-max 173.03 --latitude-min 1.47 --latitude-max 1.53\
                              --date-min '"'%(ini)s'"' --date-max '"'%(end)s'"'\
                              --depth-min 0.493 --depth-max 0.4942 --variable thetao --variable uo --variable vo --variable zos\
                              --out-dir %(folder)s --out-name sla.nc\
                              --user aespejo --pwd Copernicus580000!' %{"ini":now.strftime('%Y-%m-%d %H:%M:%S') ,"end":then.strftime('%Y-%m-%d %H:%M:%S'),"folder":folder_tmp})
                              
                              
    nc_fname = folder_tmp  + 'sla.nc'
    ncc = nc.Dataset(nc_fname)    
    time_var = ncc.variables['time']  
    time_or=nc.num2date(time_var[:],time_var.units, time_var.calendar)
    time_str=[time_or[i].strftime('%Y-%b-%d %H:%M') for i in range(0,len(time_or))]
    Time=[dt.datetime.strptime(time_str[i],'%Y-%b-%d %H:%M') for i in range(0,len(time_str))]
    # hourssince=np.array(ncc['time'])
    # time0 = dt.datetime(1950,1,1)#  hours since 1950-01-01
    # timer = [time0 + dt.timedelta(minutes=hourssince[i]*60-30,seconds=0) for i in range(0,len(hourssince))]
    # date_time = pd.to_datetime(timer)
    sla = (np.array(ncc['zos'])- 0.5442)*100# This offset comes from comparing DUACS altimetry with the model, it has been aplied to transform Sea Surface Height to Sea Level Anomaly (height from the MSL) 
    sla=sla.squeeze()
    
    
    
    
    fn = folder_tmp  + 'sla_hourly.nc'
    try:
        os.remove(fn)
    except:
        print('The system cannot find the file specified')
    ds = nc.Dataset(fn, 'w', format='NETCDF4')
    time = ds.createDimension('time',None)
    times = ds.createVariable('time', 'f8',('time',))
    times.units='hours since 1950-01-01 00:00:00'
    times.calendar='gregorian'
    SLA = ds.createVariable('SLA','f4', dimensions=('time',))
    SLA.units = 'cm'
    SLA[:]=sla
    times[:]=[nc.date2num(x,units=times.units,calendar=times.calendar)-0.5 for x in Time]
    ds.close()
    
    
    
    # # plot time series
    # now = dt.datetime.utcnow()
    # #-- create plot with tidal displacements, high and low tides and dates
    # plt.figure(figsize=(15,5)) 
    # fig,ax1 = plt.subplots(num=1)
    # ax1.plot(Time,sla,'k')
    # # ax1.plot(date_time[htindex],TIDE.data[htindex],'r*')
    # # ax1.plot(date_time[ltindex],TIDE.data[ltindex],'b*')
    # ax1.plot(numpy.matlib.repmat(now,5,1),np.linspace(-100,100,5),'r')
    # plt.ylim((-15,15))
    # for h in range(0,len(Time),24):
    #     ax1.axvline(Time[h],color='gray',lw=0.5,ls='dashed',dashes=(11,5))
    # ylbl = ax1.axes.get_yticks()
    # for h in range(0,len(ylbl)):
    #     ax1.axhline(ylbl[h],color='gray',lw=0.5,ls='dashed',dashes=(11,5))
    # # ax1.set_xlim(0,7*24)
    # ax1.set_ylabel('{0} SLA [cm]'.format('tpxo8'))
    # # ax1.set_title(u'{0:0.6f}\u00b0N {1:0.6f}\u00b0W'.format(LAT,LON))
    # fig.subplots_adjust(left=0.10,right=0.98,bottom=0.30,top=0.7)
    # plt.show()
    # #plt.close()

    print('sea level anomaly stored as ' + fn)