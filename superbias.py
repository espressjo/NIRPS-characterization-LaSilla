#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 15:40:33 2021

@author: espressjo
"""
import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt
import seaborn as sns
from os.path import basename,isfile,isdir,join
sns.set_theme()

ls = ['/nirps_raw/r20210403/NIRPS_HA_DARK093_0001.fits',
      '/nirps_raw/r20210403/NIRPS_HA_DARK093_0002.fits',
      '/nirps_raw/r20210403/NIRPS_HA_DARK093_0003.fits',
      '/nirps_raw/r20210403/NIRPS_HA_DARK093_0004.fits',
      '/nirps_raw/r20210403/NIRPS_HA_DARK093_0005.fits',
      '/nirps_raw/r20210403/NIRPS_HA_DARK093_0006.fits',
      '/nirps_raw/r20210403/NIRPS_HA_DARK093_0007.fits',
      '/nirps_raw/r20210403/NIRPS_HA_DARK093_0008.fits',
      '/nirps_raw/r20210403/NIRPS_HA_DARK093_0009.fits',
      '/nirps_raw/r20210403/NIRPS_HA_DARK093_0010.fits',
      '/nirps_raw/r20210403/NIRPS_HA_DARK093_0011.fits',
      '/nirps_raw/r20210403/NIRPS_HA_DARK093_0012.fits',
      '/nirps_raw/r20210403/NIRPS_HA_DARK093_0013.fits',
      '/nirps_raw/r20210403/NIRPS_HA_DARK093_0014.fits',
      '/nirps_raw/r20210403/NIRPS_HA_DARK093_0015.fits',
      '/nirps_raw/r20210403/NIRPS_HA_DARK093_0016.fits',
      '/nirps_raw/r20210403/NIRPS_HA_DARK093_0017.fits',
      '/nirps_raw/r20210403/NIRPS_HA_DARK093_0018.fits',
      '/nirps_raw/r20210403/NIRPS_HA_DARK093_0019.fits',
      '/nirps_raw/r20210403/NIRPS_HA_DARK093_0020.fits'] #list of dark ramp file

path = '/nirps_raw/characterization' #place where to store the superbias file and the graphs

if not isdir(path):
    from os import mkdir
    print("%s created."%path)
    mkdir(path)

def side(im):
    '''
    Return the median of side ref. pixels of an image

    '''
    l = np.append(im[:,0:4].ravel(),im[:,-4:].ravel())
    return np.median(l)
if '__main__' in __name__:
    from datetime import datetime
    fname = datetime.now().strftime("_%Y%m%d")
    time = datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
    uid = datetime.now().strftime("%y%m%d%H%M%S")
    data = [fits.open(fname)[2].data for fname in ls] #get all data
    bias_level = [side(im) for im in data]
    for i in range(len(data)):
        data[i]-=bias_level[i]
    
    #side ref. pixels bias level
    f,ax = plt.subplots(figsize=(14,14))
    name_x = [basename(f).replace('NIRPS_','').replace('.fits','') for f in ls]
    ax.set_xticklabels(name_x,rotation=-45)
    ax.plot(name_x,bias_level,'d',markersize=3)
    ax.set(xlabel='File name',ylabel='Side ref. pixels bias level (ADU)')
    plt.tight_layout()
    
    f.savefig(join(path,'bias%s.png'%fname))
    print('bias%s.png saved.'%fname)
    #show cross section of all data
    f,ax = plt.subplots(figsize=(14,14))
    y2 = [np.median(im,axis=0) for im in data]
    #ax.set_xticklabels(name_x,rotation=-45)
    for i in range(len(name_x)):
        ax.plot(y2[i],label=name_x[i])
        
    ax.set(xlabel='Pixels',ylabel='Bias cross section')
    ax.legend()
    plt.tight_layout()
    print('cross%s.png saved.'%fname)
    f.savefig(join(path,'cross%s.png'%fname))

    superbias = np.median(data,axis=0)+np.median(bias_level)

    
    hdu = fits.PrimaryHDU(data = superbias)#np.uint16(superbias))
    names = [basename(f) for f in ls]
    for name in names:
        hdu.header.add_comment(name)
    hdu.header['DATE'] = (time,"Creation date") 
    hdu.header['UNIQUEID'] = (uid, "Unique id of this file")
    if isfile(join(path,"superbias%s.fits"%fname)):
        if 'yes' in input('%s exist, do you want to overwrite it?[yes/no]'%("superbias%s.fits"%fname)):
            hdu.writeto(join(path,"superbias%s.fits"%fname),overwrite=True)
    else:
        hdu.writeto(join(path,"superbias%s.fits"%fname))
    plt.show()
    