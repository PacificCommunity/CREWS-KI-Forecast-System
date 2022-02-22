import os
from pathlib import Path
import shutil
import pandas as pd
import datetime
import xarray as xr
import numpy as np

def ingest2GUI(now):
# function to ingest results into GUI

    #Paths
    webSrvPath = 'C:/apache-tomcat-9.0.48/webapps/ROOT/ds_kiri/'
    dirpath='D:/CREWS_KI/operational_v1/runs/'
    realTimePath = 'C:/Users/test/Documents/Latest/'
    ingesterPath = 'C:/ingest/'
    baseUrl = 'http://202.0.156.66:8080/ds_kiri/inundation/'
    
    #Move nc files
    run_path = dirpath + now.strftime("%Y") + now.strftime("%m") + now.strftime("%d") + now.strftime("%H")
    dirname = Path(run_path)
    file_list = ['Gilbert_Group.nc', 'Risk_results.nc', 'Tarawa.nc']

    for x in file_list:
        print('copying '+x)
        shutil.copy(str(dirname)+'/results/'+x, 'C:/Users/test/Documents/Latest/'+x)

    #Move Inundation Files
    for i in range(0, 327):
        inunfilename="TSeries_"+str(i)+".png"
        print('copying '+inunfilename)
        shutil.copy(str(dirname)+'/results/'+inunfilename, 'C:/apache-tomcat-9.0.48/webapps/ROOT/ds_kiri/inundation/'+inunfilename)

    #Process inundation
    ##STEP 1
    print('Generating Risk')
    fn = "%sRisk_results.nc"%(realTimePath)
    ds = xr.open_dataset(fn)
    df = ds.to_dataframe()
    df.to_csv('%srisk.csv'%(ingesterPath),index=True, header=True)
    
    #STEP 2
    df = pd.read_csv('%srisk.csv'%(ingesterPath), delimiter = ",")
    os.remove('%srisk.csv'%(ingesterPath)) 
    df['time'] =  pd.to_datetime(df['time'], format='%Y-%m-%d %H:%M:%S')
    df['time'] = df['time'].apply(lambda dt: datetime.datetime(dt.year, dt.month, dt.day, \
    dt.hour,10*round((float(dt.minute) + float(dt.second)/60) / 10)))
    df_unique = df.drop_duplicates(subset = ["index"])
    df_final=df_unique.drop(['thres', 'time', 'TWL', 'SLA','Tide','Thresholds'], axis = 1)
    df_final['Primary image'] = np.nan

    df_final_tidy = df_final.rename(columns = {'Riskmax': 'Coastal Inundation risk Levels'}, inplace = False)

    for idx, row in df_final_tidy.iterrows():
        x = int(row['index'])
        df_final_tidy.loc[idx,'Primary image'] = '%sTSeries_%s.png'%(baseUrl,str(x))
        y = int(row['Coastal Inundation risk Levels'])
        riskVal = ""
        if y == 0:
            riskVal= "Low Risk"
        elif y ==1:
            riskVal= "Moderate Risk"
        elif y==2:
            riskVal= "High Risk"
        df_final_tidy.loc[idx,'Coastal Inundation risk Levels'] = riskVal

    df_final_tidy.to_csv('%stations.csv'%(ingesterPath),index=False, header=True)
    shutil.move('%stations.csv'%(ingesterPath), '%sstations.csv'%(webSrvPath))

    print("Ingestion Successful")    
    return()    

