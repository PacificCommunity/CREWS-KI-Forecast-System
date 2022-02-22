# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 08:58:04 2021

@author: antonioh
"""

# import os
import numpy as np
import numpy.matlib
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import datetime as dt
import netCDF4 as nc
import scipy.io as sio
from scipy.interpolate import interp1d


def readvars_nc(nc_fname,varname):
    ncc = nc.Dataset(nc_fname,'r') 

    if varname == 'time':
        time_var = ncc.variables[varname] 
        time_or = nc.num2date(ncc[varname],time_var.units, time_var.calendar)
        time_str=[time_or[i].strftime('%Y-%m-%d %H:%M') for i in range(0,len(time_or))]
        var = [dt.datetime.strptime(time_str[i],'%Y-%m-%d %H:%M') for i in range(0,len(time_str))]
    else:
        var = np.array(ncc[varname])
        
    ncc.close()
    return(var)


def rbfphi_gaussian(r, const):
    return np.exp(-0.5*r*r/(const*const))

def RBF_Interpolation(rbf_constant, rbf_coeff, nodes, x):

    phi = rbfphi_gaussian   # gaussian RBFs
    rbf_coeff = rbf_coeff.flatten()

    dim, n = nodes.shape
    dim_p, n_p = x.shape

    f = np.zeros(n_p)
    r = np.zeros(n)

    for i in range(n_p):
        r = np.linalg.norm(
            np.repeat([x[:,i]], n, axis=0)-nodes.T,
            axis=1
        )
        s = rbf_coeff[n] + np.sum(rbf_coeff[:n] * phi(r, rbf_constant))

        # linear part
        for k in range(dim):
            s = s + rbf_coeff[k+n+1] * x[k,i]

        f[i] = s

    return(f)

def plot_risk_maps(time_w,Risk_Cat,thresholds,fig_folder):
    
    # if tiff==1:
        ##-geotiff image
        # gtf_img = "../extras/Domain_big.tif"
        # src = rasterio.open(gtf_img)
    
    riskcolors=np.array([[0.0745, 0.6235, 1.0000], [0.8706, 0.8706, 0.3725],[1, 0.4118, 0.1608]])
    newcmp = ListedColormap(riskcolors)   
    plt.ioff()
    # for ts in range(24,len(time_w)):
       
    risk = np.amax(Risk_Cat,0)

    
    fig, ax = plt.subplots(figsize=(10,10))
    ax.set_aspect('equal')
    tcf=ax.scatter(thresholds[:,0],thresholds[:,1],10,risk,cmap=newcmp,zorder=1,vmin=-0.5, vmax=2.5)
    # if tiff==1:
    #     show(src.read(), transform=src.transform)
    ax.set_xlim(172.8,173.2)
    ax.set_ylim(1.3,1.65)
    xlbl = ax.axes.get_xticks()
    for h in range(0,len(xlbl)):
        ax.axvline(xlbl[h],color='gray',lw=0.5,ls='dashed',dashes=(11,5))
    ylbl = ax.axes.get_yticks()
    for h in range(0,len(ylbl)):
        ax.axhline(ylbl[h],color='gray',lw=0.5,ls='dashed',dashes=(11,5))
    # ax.plot(points_reef[:,0],points_reef[:,1],'*k')
    cbar = fig.colorbar(tcf, ax=ax)
    cbar.solids.set_edgecolor("face")
    cbar. set_ticks([0,1,2])
    cbar. set_ticklabels(["No flood", "Minor flood", "Moderate flood"])
    ax.set_xlabel('Longitude (degrees)')
    ax.set_ylabel('Latitude (degrees)')
    ax.set_title(' Max flood risk from '+time_w[48].strftime('%Y%m%d%H%M')+' to '+time_w[-1].strftime('%Y%m%d%H%M'))
    # fileprint=fig_folder + 'Tarawa_fld_' + time_w[ts].strftime('%Y%m%d%H%M') 
    fileprint=fig_folder + 'Tarawa_Max_fld_risk'
    plt.savefig(fileprint)
    # plt.show()
    plt.close(fig)

    return()

def risk_time_series(now,TWL,time_w,time_tide_min,sla,tide_min,thresholds,fig_folder):
    
    plt.ioff()
    for npo in range(len(TWL[1,:])):
        
        xor = [nc.date2num(i,'minutes since 1950-01-01 00:00:00') for i in time_w]
        xnew = [nc.date2num(i,'minutes since 1950-01-01 00:00:00') for i in time_tide_min]
        # xor = np.linspace(0, len(time_w)-1, len(time_w), endpoint=True).T
        # xnew = np.linspace(0, len(time_w)-1, len(time_w)*60, endpoint=True).T
        f = interp1d(xor,TWL[:,npo],kind='cubic')
        twl_min=f(xnew)
        f = interp1d(xor,sla,kind='cubic')
        sla_min=f(xnew)
        # plot time series
        
        #-- create plot with tidal displacements, high and low tides and dates
        
 
        #-- differentiate to calculate high and low tides
        diff = np.zeros_like(time_tide_min, dtype=np.float64)
        #-- forward differentiation for starting point
        diff[0] = tide_min[1] - tide_min[0]
        #-- backward differentiation for end point
        diff[-1] = tide_min[-1] - tide_min[-2]
        #-- centered differentiation for all others
        diff[1:-1] = (tide_min[2:] - tide_min[0:-2])/2.0
        htindex, = np.nonzero((np.sign(diff[0:-1]) >= 0) & (np.sign(diff[1:]) < 0))
        # ltindex, = np.nonzero((np.sign(diff[0:-1]) <= 0) & (np.sign(diff[1:]) > 0))
        
        
        
        
        plt.figure(figsize=(25,10)) 
        # fig,ax1 = plt.subplots(num=1)
        plt.fill_between(time_tide_min,twl_min,0,color='paleturquoise',label='Total water level')
        plt.fill_between(time_tide_min,tide_min+sla_min,0,color='c',label='Ocean water level')
        plt.plot(time_tide_min,twl_min,color='paleturquoise',linewidth=0.5)
        plt.plot(time_tide_min,tide_min,'--k',linewidth=0.5,label='Astronomical tide')
        plt.plot(np.matlib.repmat(now,5,1),np.linspace(-1,thresholds[npo,4]+0.6,5),color= 'gray',linewidth=1)
        plt.plot(time_w,np.matlib.repmat(thresholds[npo,4],len(time_w),1),color=[1, 0.4118, 0.1608])
        plt.text(time_tide_min[-4000], thresholds[npo,4]+0.02, 'Moderate flood threshold',color=[1, 0.4118, 0.1608],size=15)
        plt.plot(time_w,np.matlib.repmat(thresholds[npo,2],len(time_w),1),color=[0.8706, 0.8706, 0.3725])
        plt.text(time_tide_min[-3800], thresholds[npo,2]+0.02, 'Minor flood threshold',color=[0.8706, 0.8706, 0.3725],size=15)
        #ax1.plot(time_tide_min[htindex],twl_min[htindex],"v",color= 'gray',markersize=5)
        for h in range(2,len(htindex)-1):
            text=time_tide_min[htindex[h]].strftime("%H:%M")
            plt.plot(time_tide_min[htindex[h]],twl_min[htindex[h]],"v",color= 'gray',markersize=5)
            plt.text(time_tide_min[htindex[h]]-dt.timedelta(hours=3),twl_min[htindex[h]]+0.08,text,color= 'gray',size=15) 
        plt.plot(now,thresholds[npo,4]+0.6,"v",color= 'gray',markersize=25)
        plt.text(now-dt.timedelta(hours=3),thresholds[npo,4]+0.6+0.09, 'Now',color='gray',size=15)
        plt.ylim((-1,thresholds[npo,4]+0.6))       
        plt.xlim((time_tide_min[60*24],time_tide_min[-1]))
        plt.xlabel('UTC Time',fontsize=15)
        plt.ylabel('Water level [m]',fontsize=15)
        plt.yticks(fontsize=15)
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        for h in range(0,len(time_w),24):
            plt.axvline(time_w[h].replace(hour=0),color='gray',lw=0.5,ls='dashed',dashes=(11,5))
        ylbl = plt.yticks()[0]
        for h in range(0,len(ylbl)):
            plt.axhline(ylbl[h],color='gray',lw=0.5,ls='dashed',dashes=(11,5))
        plt.legend(loc='upper center', ncol = 3,fontsize=20)
        # plt.show()
        fileprint=fig_folder + 'TSeries_' + str(npo) 
        plt.savefig(fileprint)
        plt.close()
    return()

  
def TWL_points2nc(result_folder,Risk_Cat,TWL,time_w,time_tide_min,sla,tide_min,thresholds):

    x=thresholds[:,0]
    y=thresholds[:,1]
    thresholdss = thresholds[:,[2,4]]
    time_10min = [time_tide_min[i] for i in range(1440, len(time_tide_min), 10)]# vhange the starting time into range to the lenght of the spin-up  in minutes
        
    xorhour = [nc.date2num(i,'minutes since 1950-01-01 00:00:00') for i in time_w]
    xormin = [nc.date2num(i,'minutes since 1950-01-01 00:00:00') for i in time_tide_min]
    xor10min = [nc.date2num(i,'minutes since 1950-01-01 00:00:00') for i in time_10min]
    
    
    f = interp1d(xormin,tide_min)
    tide_10min=f(xor10min)
    f = interp1d(xorhour,sla)
    sla_10min=f(xor10min)   
    
    riskmax=np.max(Risk_Cat,axis=0)

    TWL_10min = np.zeros((len(time_10min),len(x)))       
    for npo in range(len(x)):          
        f = interp1d(xorhour,TWL[:,npo],kind='cubic')
        TWL_10min[:,npo] =f(xor10min)

    
    
    fn = result_folder  + 'Risk_results.nc'          
    ds = nc.Dataset(fn, 'w', format='NETCDF4')
    time = ds.createDimension('time', None)
    index = ds.createDimension('index', len(riskmax))
    thres = ds.createDimension('thres', None)
    index = np.arange(0, len(riskmax))
    times = ds.createVariable('time', 'f8', ('time',))
    times.units='hours since 1950-01-01 00:00:00'
    times.calendar='gregorian'
    times[:] = [nc.date2num(x,units=times.units,calendar=times.calendar) for x in time_10min]
    lonnc = ds.createVariable('lon', 'f4', ('index',))
    lonnc.units ='degrees_east'
    lonnc[:] = x
    latnc = ds.createVariable('lat', 'f4', ('index',))
    latnc.units ='degrees_north'
    latnc[:] = y
    TWLnc= ds.createVariable('TWL', 'f8', ('index','time'))
    TWLnc.units = 'm'
    TWLnc[:,:] =  TWL_10min.T
    SLAnc= ds.createVariable('SLA', 'f8', ('time'))
    SLAnc.units = 'm'
    SLAnc[:] =  sla_10min
    Tidenc= ds.createVariable('Tide', 'f8', ('time'))
    Tidenc.units = 'm'
    Tidenc[:] =  tide_10min
    Riskmaxnc= ds.createVariable('Riskmax', 'f4', ('index'))
    Riskmaxnc.units = 'No risk=0, Minor=1, Moderate=2'
    Riskmaxnc[:] =  riskmax
    Tresholdsnc= ds.createVariable('Thresholds', 'f8', ('index','thres'))
    Tresholdsnc.units = 'm'
    Tresholdsnc[:,:] =  thresholdss
    ds.close()

   
###############################################################################
def make_flood(now,fmap,fseries):
    
    out_name='../runs/' + now.strftime("%Y") + now.strftime("%m") + now.strftime("%d") +  now.strftime("%H")  +'/'
    result_folder = out_name + 'results/'
    folder_tmp ='../tmp/'
    
    nc_fname = folder_tmp  + 'sla_hourly.nc'
    time_sla = readvars_nc(nc_fname,'time')
    sla = readvars_nc(nc_fname,'SLA')/100
    
    nc_fname = folder_tmp  + 'tide_hourly.nc'
    time_tide = readvars_nc(nc_fname,'time')
    tide = readvars_nc(nc_fname,'tide')/100
    
    nc_fname = folder_tmp  + 'tide_minute.nc'
    time_tide_min = readvars_nc(nc_fname,'time')
    tide_min = readvars_nc(nc_fname,'tide')/100
    
    
    nc_fname = folder_tmp  + 'wind_and_waves.nc'
    time_w = readvars_nc(nc_fname,'time')
    Hs = readvars_nc(nc_fname,'Hs')
    Tm = readvars_nc(nc_fname,'Tm')
    Tp = readvars_nc(nc_fname,'Tp')
    holo_o=Hs/(1.56*Tp**2)# Hs-wave lenght relation in deep waters
    Wx12 = readvars_nc(nc_fname,'Windx')[327]# 327 is the point used to predict water levels inside the lagoon, we call it , Point 12
    Wy12 = readvars_nc(nc_fname,'Windy')[327]
    Hs12 = Hs[327]
    Tm12 = Tm[327]
    Dir12 = readvars_nc(nc_fname,'Dir')[327]
    ncc = nc.Dataset(nc_fname,'r') 
    point_type = ncc.point_type
    ncc.close()
    
    
    ##- load lagoon pooints
    f = open('../extras/Lagoon_profiles_xyxy.txt', 'r') 
    lago_prof = np.genfromtxt(f, delimiter='  ')
    f.close()
    npointslago = len(lago_prof)
    points_lago = lago_prof[:,2:4]
    points_lago_pos=np.arange(len(lago_prof))
    
    
    ##- load rbr coefficients obtained from ADCIRC+SWAN database
    fl_name = '../extras/Lagoon_RBF_coefficients.mat'
    infos = sio.whosmat(fl_name)
    names, sizes, types  = zip(*infos)
    MAT = sio.loadmat(fl_name)
    maximos_surge_wave = MAT['maximos_surge_wave']
    minimos_surge_wave = MAT['minimos_surge_wave']
    rbf_surge_coeficients= MAT['rbf_surge_coeficients']
    
    
    ##- load forereef pooints
    f = open('../extras/Forereef_profiles_vf.txt', 'r') 
    reef_prof = np.genfromtxt(f, delimiter='  ')
    f.close()
    npointsreef = len(reef_prof)
    points_reef = reef_prof[:,2:4]
    points_reef_pos=np.arange(len(lago_prof),len(lago_prof)+len(reef_prof))
    
    dep_reef = reef_prof[:,5]
    
    ##- load rbr coefficients obtained from BEWARE database
    fl_name = '../extras/Forereef_RBF_coeff_min_max.mat'
    infos = sio.whosmat(fl_name)
    names, sizes, types  = zip(*infos)
    MAT = sio.loadmat(fl_name)
    maximos_reef = MAT['maximos']
    minimos_reef = MAT['minimos']
    rbf_reef = MAT['RBF_coeficients']
    prof_reef = MAT['Prof']
    
    
    
    ##- load lagoon pooints
    f = open('../extras/Tarawa_inundation_thresholds.dat', 'r') 
    thresholds = np.genfromtxt(f, delimiter='  ')
    f.close()
    
    
    
    
    ##############################################################################
    
    Dir12=270-Dir12
    Hsx12=Hs12*np.cos(np.radians(Dir12))
    Hsy12=Hs12*np.sin(np.radians(Dir12))
    
    
    datos = np.vstack((Hsx12,Hsy12,Tm12,Wx12,Wy12))
    datos = datos.T
    datos_n=np.empty((len(datos[:,0]),len(datos[0,:])))
    for ii in range(len(datos[0,:])):
        datos_n[:,ii]=(datos[:,ii]-minimos_surge_wave[ii])/(maximos_surge_wave[ii]-minimos_surge_wave[ii])
    
    
    levels_lagoon=np.empty((len(time_w),len(lago_prof)))
    
    for i in range(len(lago_prof)):
            
        rbf_constant=rbf_surge_coeficients[0][i][0][0][0][1]
        rbf_nodes=rbf_surge_coeficients[0][i][0][0][0][4]
        rbf_coeff=rbf_surge_coeficients[0][i][0][0][0][7]  
        su = RBF_Interpolation(rbf_constant,rbf_coeff,rbf_nodes, datos_n.T) 
        su[su<0] = 0
        levels_lagoon[:,i]=su
    
   
    ru_lagoon=1.1*0.35*np.matlib.repmat(lago_prof[:,4],len(time_w),1)*(Hs[points_lago_pos].T*1.56*Tp[points_lago_pos].T**2)**(1/2)
    TWL_lagoon=levels_lagoon+ru_lagoon+np.tile(sla,(len(lago_prof),1)).T+np.tile(tide,(len(lago_prof),1)).T
    
    
    levels_forereef=np.empty((len(time_w),len(prof_reef)))
    
    for i in range(len(prof_reef)):
            
        dep=np.matlib.repmat(prof_reef[i,5],1,len(tide))+tide+sla
        dep = np.squeeze(dep.T)
        #-Min max depth in BEWARE
        dep[dep<-1]=-1
        dep[dep>3]=3
        Hsp=Hs[i+len(lago_prof)]
        holop=holo_o[i+len(lago_prof)]
        holop[holop<0.0050]=0.0050
        holop[holop>0.0500]=0.0500
        
        datos = np.vstack((Hsp,holop,dep))
        datos = datos.T
        
        datos_n=np.empty((len(datos[:,0]),len(datos[0,:])));
        for ii in range(len(datos[0,:])):
            datos_n[:,ii]=(datos[:,ii]-minimos_reef[ii])/(maximos_reef[ii]-minimos_reef[ii])
        
    
        rbf_constant=rbf_reef[0][i][0][0][0][1]
        rbf_nodes=rbf_reef[0][i][0][0][0][4]
        rbf_coeff=rbf_reef[0][i][0][0][0][7]
    
         
        ru = RBF_Interpolation(rbf_constant,rbf_coeff,rbf_nodes, datos_n.T)
        ru[ru<0] = 0
        

        plt.figure(figsize=(25,10)) 
        plt.plot(ru)
        plt.plot(Hs[i+len(lago_prof)])
      
        Betaf=0.05 # Asume intermediate beach slope 1/20
        L0=Hsp*1/holop
        L0[L0>975]=975# 25s Tp maximun wave length
    
        ru2=0.73*Betaf*(Hsp*L0)**(1/2)## runup reflective 
        #Runup aproximation for the waves below 1 m
        ru[Hsp<1]=ru2[Hsp<1]
        levels_forereef[:,i]=ru
    
    
    TWL_forerref=levels_forereef+np.tile(sla,(len(prof_reef),1)).T+np.tile(tide,(len(prof_reef),1)).T;
    
    TWL = np.concatenate((TWL_lagoon,TWL_forerref),axis=1)
    
    
    Risk_Cat=np.zeros(TWL.shape)
    for ts in range(len(TWL[1,:])):
        Risk_Cat[(TWL[:,ts]>=thresholds[ts,4]),ts]=2
        Risk_Cat[(TWL[:,ts]>=thresholds[ts,2]) & (TWL[:,ts]<=thresholds[ts,4]),ts]=1                      
      
    Risk_Cat=Risk_Cat[47:,:]
    
    TWL_points2nc(result_folder,Risk_Cat,TWL,time_w,time_tide_min,sla,tide_min,thresholds)
    
    if fmap==1:
        plot_risk_maps(time_w,Risk_Cat,thresholds,result_folder)
        
    if fseries==1:
        risk_time_series(now,TWL,time_w,time_tide_min,sla,tide_min,thresholds,result_folder)
    
