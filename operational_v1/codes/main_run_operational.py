# -*- coding: utf-8 -*-
"""
Created on Sat Jun 12 18:16:18 2021

@author: antonioh
"""
import datetime as dt
from datetime import timedelta
import os
import shutil
import requests
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry
from bs4 import BeautifulSoup
import step1_download_NCEP as step1
import step2_download_CMEMS as step2
import step3_gen_tide_TPOX8 as step3
import step4_make_wave_forcing as step4
import step5_make_wind_forcing as step5
import step6_parallelize_run as step6
import step7_postprocess_output as step7 
import step8_make_flood_risk as step8
import step9_archive_output as step9
import step10_ingest2GUI as step10


def list_available_runs(url):
    session = requests.Session()
    retry = Retry(connect=5, backoff_factor=1)
    adapter = HTTPAdapter(max_retries=retry)
    session.mount(url, adapter)

    try:
        req = session.get(url).text
        session.close() 
        soup = BeautifulSoup(req, 'html.parser')     
        x = (soup.find_all('a'))
        runs = []
        for i in x:
            file_name = i.extract().get_text()
            runs.append(int(file_name))
    except:
        
        runs = []
        print('Keep working on making the dowinloading process more robut')
        
      
    return(runs)



def delete_ndaysbefore(now,ndays):
    folder_name='../runs/'
    flist = os.listdir(folder_name)
    for d in flist:
        rundate=dt.datetime.strptime(d+'0000',"%Y%m%d%H%M%S")
        if rundate<now-timedelta(ndays):
            shutil.rmtree('../runs/'+ d)
    return()

def delete_nmonthsbefore(now,nmonths):
    folder_name='../archives/'
    flist = os.listdir(folder_name)
    for d in flist:
        nb = d.split("_")[0]
        rundate=dt.datetime.strptime(nb+'0000',"%Y%m%d%H%M%S")
        if rundate<now-timedelta(nmonths * 30):
            shutil.rmtree('../archives/'+ d)
    return()
            
    

###############################################################################

#- Find the llatest available run in nomads.ncep.noaa.gov
now = dt.datetime.utcnow()
url = 'https://nomads.ncep.noaa.gov/cgi-bin/filter_gfswave.pl?dir=%2Fgfs.' + now.strftime("%Y") + now.strftime("%m") + now.strftime("%d")
runs = list_available_runs(url)
if len(runs)==0:
    now = dt.datetime.utcnow()-timedelta(1)
    url = 'https://nomads.ncep.noaa.gov/cgi-bin/filter_gfswave.pl?dir=%2Fgfs.' + now.strftime("%Y") + now.strftime("%m") + now.strftime("%d")
    runs = list_available_runs(url)

#- Define the run to be used
runs=sorted(runs)
now = now.replace(hour=runs[-1],minute=0,second=0,microsecond=0)


# delete previous runs older than 14 days
try:
    delete_ndaysbefore(now,14)
except Exception as e:
    print(e)
    
    

step1.download_NCEP(now)
step2.download_CNEMS(now)
step3.gen_tide(now)
step4.make_waves(now)
step5.make_winds(now)
step6.par_run(now)
step7.postprocess_SWAN(now,1)
step8.make_flood(now,1,1)

# delete previous archived swan outputs older than 3 months
try:
    delete_nmonthsbefore(now,3)
except Exception as e:
    print(e)
    
step9.archive_output(now)
step10.ingest2GUI(now)

