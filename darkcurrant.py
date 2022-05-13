#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 08:19:23 2021

@author: noboru
"""

from findramp import findramp 
from findread import findread
from astropy.stats import sigma_clipped_stats as sc
from astropy.io import fits
from os.path import join,isfile
from hxrg.type.ramp import hxramp 
from hxrg.type.read import hxread 
from tqdm import tqdm
from matplotlib import pyplot as plt
import seaborn as sns
sns.set_theme()
import numpy as np
path = '/nirps_raw/characterization' #place where to store the superbias file and the graphs
gpm = np.ones((4096,4096),dtype=bool)
#if gpm.fits is there take it
if isfile(join(path,'gpm.fits')):
    print('using gpm.fits')
    gpm = np.asarray(fits.getdata(join(path,'gpm.fits')),dtype=bool)
def repport(dc,show=True,gain=1.27):
    dc = dc*gain
    stats = sc(dc,sigma=5)
    print("Dark currant: %f +/- %f"%(stats[1],stats[2]))
    dc = dc[(dc>np.percentile(dc,0.13)) & (dc<np.percentile(dc,99.87))]
    count, bins_count = np.histogram(dc, bins=200)
    pdf = count / sum(count)
    cdf = np.cumsum(pdf)
    f,ax = plt.subplots(2,1)
    ax[0].hist(dc,bins=200)
    ax[1].plot(bins_count[1:], cdf)
    ax[1].set(xlabel='Dark current (e$^-$/second)',title='Cummulative Sum')
    ax[0].set(xlabel='Dark currant (e$^-$/second)',ylabel='Count',title='Dark current distribution')
    plt.tight_layout()
    f.savefig(join(path,'darkcurrant_dist.png'))
    if show:
        plt.show()
def makedarkcurrent(uid,sideonly=True,toponly=False):
    fr = findread()
    ls = fr.findreads(uid)
    read = hxread()
    R = hxramp(sideonly=sideonly,toponly=toponly)
    print("Making uid: %s"%uid)
    for f in tqdm(ls):
        read(f)
        R<<read
    R.fit()
    return R.a[gpm]
def darkcurrant(uids):
    dc = [];
    for uid in uids:
        fr = findramp()
        rfile = fr.findramp(uid)
        slope = fits.open(rfile)['slope'].data
        dc.append(sc(slope[gpm],sigma=5))
    return dc
def darkcurrantUdeM(ls):
    dc = []; 
    for rfile in ls:   
        slope = fits.getdata(rfile)[0]
        dc.append(sc(slope[gpm],sigma=5))
    return dc

if '__main__' in __name__:
    
    lsuids =    ['20210406204047',
             '20210406205100',
             '20210406210113',
             '20210406211126',
             '20210406212139',
             '20210406213152',
             '20210406214205',
             '20210406215218',
             '20210406220231',
             '20210406221245',
             '20210406222258',
             '20210406223311',
             '20210406224324',
             '20210406225337',
             '20210406230350',
             '20210406231403',
             '20210406232416',
             '20210406233429',
             '20210406234442',
             '20210406235455',
             '20210407000508',
             '20210407001521',
             '20210407002534',
             '20210407003547',
             '20210407004600',
             '20210407005613',
             '20210407010626',
             '20210407011639',
             '20210407012652',
             '20210407013706',
             '20210407014719',
             '20210407015732',
             '20210407020745',
             '20210407021758',
             '20210407022811',
             '20210407023824',
             '20210407024837',
             '20210407025850',
             '20210407030903',
             '20210407031916',
             '20210407032929',
             '20210407033942',
             '20210407034955',
             '20210407040008',
             '20210407041021',
             '20210407042034',
             '20210407043047',
             '20210407044100',
             '20210407045113',
             '20210407050127',
             '20210407051140',
             '20210407052153',
             '20210407053206',
             '20210407054219',
             '20210407055232',
             '20210407060245',
             '20210407061258',
             '20210407062311',
             '20210407063324',
             '20210407064337']#data for NIRPS
    '''
    ls_udem = ['/nasmtl01/HxRGdata/grosdisk/grosdisk/data/20181113223015/NIRPS_UdeM_DARK_01.fits',
               '/nasmtl01/HxRGdata/grosdisk/grosdisk/data/20181113223015/NIRPS_UdeM_DARK_02.fits',
               '/nasmtl01/HxRGdata/grosdisk/grosdisk/data/20181113223015/NIRPS_UdeM_DARK_03.fits',
               '/nasmtl01/HxRGdata/grosdisk/grosdisk/data/20181113223015/NIRPS_UdeM_DARK_04.fits',
               '/nasmtl01/HxRGdata/grosdisk/grosdisk/data/20181113223015/NIRPS_UdeM_DARK_05.fits',
               '/nasmtl01/HxRGdata/grosdisk/grosdisk/data/20181113223015/NIRPS_UdeM_DARK_06.fits',
               '/nasmtl01/HxRGdata/grosdisk/grosdisk/data/20181113223015/NIRPS_UdeM_DARK_07.fits',
               '/nasmtl01/HxRGdata/grosdisk/grosdisk/data/20181113223015/NIRPS_UdeM_DARK_08.fits',
               '/nasmtl01/HxRGdata/grosdisk/grosdisk/data/20181113223015/NIRPS_UdeM_DARK_09.fits',
               '/nasmtl01/HxRGdata/grosdisk/grosdisk/data/20181113223015/NIRPS_UdeM_DARK_10.fits',
               '/nasmtl01/HxRGdata/grosdisk/grosdisk/data/20181113223015/NIRPS_UdeM_DARK_11.fits',
               '/nasmtl01/HxRGdata/grosdisk/grosdisk/data/20181113223015/NIRPS_UdeM_DARK_12.fits',
               '/nasmtl01/HxRGdata/grosdisk/grosdisk/data/20181113223015/NIRPS_UdeM_DARK_13.fits',
               '/nasmtl01/HxRGdata/grosdisk/grosdisk/data/20181113223015/NIRPS_UdeM_DARK_14.fits',
               '/nasmtl01/HxRGdata/grosdisk/grosdisk/data/20181113223015/NIRPS_UdeM_DARK_15.fits',
               '/nasmtl01/HxRGdata/grosdisk/grosdisk/data/20181113223015/NIRPS_UdeM_DARK_16.fits'
               ]
    dc = darkcurrant(lsuids)
    dc_mean,dc_med,dc_std = zip(*dc)
    darkc_copl,_,darkerr_copl = sc(dc_mean)
    with open('dark_copl.txt','w') as f:
        for i in range(len(dc_mean)):
            f.write("%f,%f,%f\n"%(dc_mean[i],dc_med[i],dc_std[i]))
    
    dc = darkcurrantUdeM(ls_udem)
    dc_mean,dc_med,dc_std = zip(*dc)
    with open('dark_udem.txt','w') as f:
        for i in range(len(dc_mean)):
            f.write("%f,%f,%f\n"%(dc_mean[i],dc_med[i],dc_std[i]))
    darkc_udem,_,darkerr_udem = sc(dc_mean)
    print("NIRPS [COPL]: %.4f+/-%.4f adu/s"%(darkc_copl,darkerr_copl))
    print("NIRPS [UdeM]: %.4f+/-%.4f adu/s"%(darkc_udem,darkerr_udem))
    '''
    # from datetime import datetime
    # gain = 1.27
    # lfile = '/nirps_raw/characterization/dc-%s.csv'%(datetime.now().strftime("%Y%m%d%H%M"))
    # with open(lfile,'a') as f:
    #     f.write("uid,mean,median,std\n")
    # for uid in lsuids:
    #     dc = makedarkcurrent(uid)
    #     stats = sc(dc*gain,sigma=5)
    #     with open(lfile,'a') as f:
    #         f.write("%s,%f,%f,%f\n"%(uid,stats[0],stats[1],stats[2]))
    #     print("%s,%f,%f,%f\n"%(uid,stats[0],stats[1],stats[2]))
    #repport(dc,show=True)
    uid = lsuids[10]
    print(uid)

    dc = makedarkcurrent(uid,sideonly=True,toponly=False)
    repport(dc)
    
    
