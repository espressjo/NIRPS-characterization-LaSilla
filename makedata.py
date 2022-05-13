#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 08:24:58 2021

@author: noboru
"""
from findread import findread
import numpy as np
from astropy.io import fits
from os.path import join
from tqdm import tqdm
from hxrg.fits2ramp.refpixcorr import refpxcorr as rpc
from datetime import datetime

lsuid = [ '20210403171440',
          '20210403175352',
          '20210403183253',
          '20210403191153',
          '20210403195054',
          '20210403203006',
          '20210403210907',
          '20210403214818',
          '20210403222719',
          '20210403230620']#list of flat (200 reads separated by darks) (unique ID)


path = "/nirps_raw/characterization"
sbias_file = "superbias_20210416.fits"

#load the bias file
hdu = fits.open(join(path,sbias_file))
superbias = np.asarray(hdu[0].data,dtype=float)
superbias_uid = hdu[0].header['UNIQUEID']
hdu.close()

def _with_bias():
    fr = findread()
    print("Making data")
    for uid in (lsuid):
        ls = fr.findreads(uid)
        cube = np.zeros((len(ls),4096,4096))
        print("working on %s."%uid)
        for i in tqdm(range(len(ls))):
            cube[i,:,:] = rpc(np.asarray(fits.getdata(ls[i]),dtype=float)-superbias).refpxcorrtop()
        hdu = fits.PrimaryHDU(data = cube)
        time = datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
        ID = datetime.now().strftime("%y%m%d%H%M%S")
        hdu.header['DATE'] = (time,"Creation date")
        hdu.header['UNIQUEID'] = (ID,"Unique id of this file")
        hdu.header['RUID'] = (uid,"Unique id of the ramp used to built this file")
        hdu.header['SBIASID'] = (superbias_uid,"Unique id of the super bias")
        hdu.header['SBIASF'] = (sbias_file,"Super bias filename")
        
        #save the files
        hdu.writeto(join(path,uid+".bias.fits"),overwrite=True)
        print("%s saved."%join(path,uid+".fits"))
def _no_bias():
    fr = findread()
    print("Making data")
    for uid in (lsuid):
        ls = fr.findreads(uid)
        bias = np.asarray(fits.getdata(ls[0]),dtype=float)
        ls = ls[1:]
        cube = np.zeros((len(ls),4096,4096))
        print("working on %s."%uid)
        for i in tqdm(range(len(ls))):
            cube[i,:,:] = rpc(np.asarray(fits.getdata(ls[i]),dtype=float)-bias).refpxcorrtop()
        hdu = fits.PrimaryHDU(data = cube)
        time = datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
        ID = datetime.now().strftime("%y%m%d%H%M%S")
        hdu.header['DATE'] = (time,"Creation date")
        hdu.header['UNIQUEID'] = (ID,"Unique id of this file")
        hdu.header['RUID'] = (uid,"Unique id of the ramp used to built this file")

        
        #save the files
        hdu.writeto(join(path,uid+".fits"),overwrite=True)
        print("%s saved."%join(path,uid+".fits"))
if '__main__' in __name__:
    from sys import argv
    if '-selfbias' in argv:
        _no_bias()
    else:
        _with_bias()
        
        
