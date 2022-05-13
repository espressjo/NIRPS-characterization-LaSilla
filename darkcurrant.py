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


path = '/nirps_raw/nirps/characterization'#place where to store the superbias file and the graphs
p_data = "/nirps_raw/nirps/reads"
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
def makedarkcurrent(fname,sideonly=True,toponly=False):
    if '/' not in fname:
        fname = join(p_data,fname)
    
    read = hxread()
    R = hxramp(sideonly=sideonly,toponly=toponly)
    cube = fits.getdata(fname)
    for im in tqdm(cube):
        read(np.asarray(im,dtype=float))
        read.H = fits.getheader(fname)
        R<<read
    R.fit()
    return R.a[gpm]


if '__main__' in __name__:
    darks = ["NIRPS_2022-05-05T15_31_01_590.fits",
            "NIRPS_2022-05-05T15_40_24_489.fits",
            "NIRPS_2022-05-05T15_49_47_388.fits",
            "NIRPS_2022-05-05T15_59_10_285.fits"]

    '''
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
    from datetime import datetime
    gain = 1.27
    lfile = join(path,'dc-%s.csv'%(datetime.now().strftime("%Y%m%d%H%M")))
    with open(lfile,'a') as f:
        f.write("uid,mean,median,std\n")
    for fname in darks:
        dc = makedarkcurrent(fname)
        stats = sc(dc*gain,sigma=5)
        with open(lfile,'a') as f:
            f.write("%s,%f,%f,%f\n"%(fname,stats[0],stats[1],stats[2]))
        print("%s,%f,%f,%f\n"%(fname,stats[0],stats[1],stats[2]))
    repport(dc,show=True)


'''
    dc = makedarkcurrent(uid,sideonly=True,toponly=False)
    repport(dc)
'''  
    
