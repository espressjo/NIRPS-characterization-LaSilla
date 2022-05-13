#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 09:43:00 2021

@author: noboru
"""
from os.path import isfile,join
from findread import findread
from findramp import findramp
from hxrg.fits2ramp.refpixcorr import refpxcorr as rpc
from astropy.io import fits
import numpy as np
from datetime import datetime
from astropy.stats import sigma_clipped_stats as sc
from scipy.ndimage import median_filter
from hxrg.type.read import hxread
from hxrg.type.ramp import hxramp
from scipy.ndimage import zoom
from tqdm import tqdm
class badpixelmap():
    def __init__(self):
        self.path = '/nirps_raw/nirps/characterization'
        self.pathtmp = '/nirps_raw/nirps/.tmp'
        self.nonlin = 'nonlin.fits'
        self.time_fmt = "%Y-%m-%dT%H:%M:%S"
        self.ls_dark = [];
        self.ls_flat = [];
        self.p_data = "/nirps_raw/nirps/reads"
        
        self.darks = [ "NIRPS_2022-05-05T15_31_01_590.fits",
                      "NIRPS_2022-05-05T15_40_24_489.fits",
                      "NIRPS_2022-05-05T15_49_47_388.fits",
                      "NIRPS_2022-05-05T15_59_10_285.fits"]
              
        self.flats= ["NIRPS_2022-05-05T16_09_12_198.fits",
                     "NIRPS_2022-05-05T16_30_17_336.fits",
                     "NIRPS_2022-05-05T16_51_22_473.fits",
                     "NIRPS_2022-05-05T17_12_27_612.fits",
                     "NIRPS_2022-05-05T17_33_32_752.fits",
                     "NIRPS_2022-05-05T18_15_43_030.fits",
                     "NIRPS_2022-05-05T18_36_48_170.fits",
                     "NIRPS_2022-05-05T19_18_58_450.fits",
                     "NIRPS_2022-05-05T19_40_03_591.fits",
                     "NIRPS_2022-05-05T20_01_08_734.fits",
                     "NIRPS_2022-05-05T20_22_13_876.fits",
                     "NIRPS_2022-05-05T20_43_19_017.fits",
                     "NIRPS_2022-05-05T21_04_18_580.fits",
                     "NIRPS_2022-05-05T21_25_29_302.fits",
                     "NIRPS_2022-05-05T22_07_45_164.fits",
                     "NIRPS_2022-05-05T22_49_55_447.fits"]
        
        if not isfile(join(self.path,self.nonlin)):
            print("no nonlin.fits found.")
            self.nonlin = ''
    def _make_cds_data(self,f,n=2):
        '''
        For NIRPS we want to avoid using the bottom ref. pix. and 
        we do not want to use the odd-even option.
        Parameters
        ----------
        f: file (cube of data)
        n: if you want to create someting else like read5-read1, n=5
        -------
        TYPE
            CDS np.array of float with corrected ref.px. using side and top.
        '''
        arr = fits.getdata(join(self.p_data,f))[:n,:,:]
        return rpc(np.asarray(arr[n-1],dtype=float)-np.asarray(arr[0],dtype=float)).refpxcorrtop(False)


    def hotpixels(self):
        fr = findread()
        cube = [];
        hp = np.zeros((4096,4096),dtype=int)
        for f in self.darks:
            cube.append(self._make_cds_data(f))
        
        dark = np.median(cube,axis=0)
        mean,_,std = sc(dark,sigma=5)
        hp[dark>(mean+5*std)]=1
        hdu = fits.PrimaryHDU(data=hp)
        hdu.header['DATE'] = (datetime.now().strftime(self.time_fmt),"Creation date")
        hdu.header['UNIQUEID'] = (datetime.now().strftime("%y%m%d%H%M%S"),"Unique identification number")
        hdu.header.add_comment("Built using:")
        
        for f in self.darks:
            if "HIERARCH ESO DET RAMP ID" in fits.getheader(join(self.p_data,f)):               
                hdu.header.add_comment(fits.getheader(join(self.p_data,f))["HIERARCH ESO DET RAMP ID"])

        hdu.writeto(join(self.pathtmp,"hotpixels.fits"),overwrite=True)
        print("hotpixels.fits write successfully.")
        return hp
    def badpixels(self):

        #check slope. We need to build ramp file
        c = fits.getdata(join(self.p_data,self.flats[0]))
        c_H = fits.getheader(join(self.p_data,self.flats[0]))
        read = hxread()
        ramp = hxramp()
        print("Making a ramp file with the 1st flat")
        for im in tqdm(c):
            read(np.asarray(im,dtype=float))
            read.H = c_H
            ramp << read
        ramp.fit()
        slope = ramp.a
        #fname_ramp = fr.findramp(self.ls_flat[0])
        #slope = fits.open(fname_ramp)['slope'].data
        mean,med,std = sc(slope[1800:2200,1800:2200],sigma=5)
        n = int(40000./(mean*5.5733))
        print("40,000ADU ~ n=%d"%n)
        cube_bpm = [];
        for f in self.flats:
            cube_bpm.append(self._make_cds_data(f,n=n))
        im = np.median(cube_bpm,axis=0)
        
        pix_bin = 64 # pixel size of the binning box
        nbin = 4096//pix_bin
        percentile = 50
        box = np.zeros([nbin,nbin])
        for i in (range(nbin)):
            for j in range(nbin):
                # -1 sigma of distrubution. Should be good to remove illuminated pixels
                box[i,j] = np.nanpercentile(im[i*pix_bin:i*pix_bin+pix_bin,j*pix_bin:j*pix_bin+pix_bin],percentile)

        lowf = zoom(box,pix_bin)
        im/=lowf

        mask = np.zeros_like(im)
        mask[im<0.85] = 1
        mask[im>1.15] = 1
        mask[:4,:]=0
        mask[-4:,:] = 0
        mask[:,:4] = 0
        mask[:,-4:] = 0
        
        
        hdu = fits.PrimaryHDU(data=mask)
        hdu.header['DATE'] = (datetime.now().strftime(self.time_fmt),"Creation date")
        hdu.header['UNIQUEID'] = (datetime.now().strftime("%y%m%d%H%M%S"),"Unique identification number")
        hdu.header.add_comment("Built using:")
        for f in self.flats:
            if "HIERARCH ESO DET RAMP ID" in fits.getheader(join(self.p_data,f)):               
                hdu.header.add_comment(fits.getheader(join(self.p_data,f))["HIERARCH ESO DET RAMP ID"])
        hdu.writeto(join(self.pathtmp,"badpixels.fits"),overwrite=True)
    def nonlinbpm(self):
        if self.nonlin=='':
            return
        nlbpm = fits.getdata(join(self.path,self.nonlin))[2]
        hdu = fits.PrimaryHDU(data=nlbpm)
        hdu.header['DATE'] = (datetime.now().strftime(self.time_fmt),"Creation date")
        hdu.header['UNIQUEID'] = (datetime.now().strftime("%y%m%d%H%M%S"),"Unique identification number")
        hdu.header.add_comment("Built using:")
        hdu.header.add_comment(fits.getheader(join(self.path,self.nonlin))['UNIQUEID'])
        hdu.writeto(join(self.pathtmp,'nl-bpm.fits'))
    def gpm(self):
        gpm = np.ones((4096,4096),dtype=int)
        uid = [];
        if isfile(join(self.pathtmp,'hotpixels.fits')):
            print("Using hotpixels")
            hp = fits.getdata(join(self.pathtmp,'hotpixels.fits'))
            uid.append(fits.getheader(join(self.pathtmp,'hotpixels.fits'))['UNIQUEID'])
            gpm[hp==1]=0
        if isfile(join(self.pathtmp,'nl-bpm.fits')):
            print("Using nl-bpm")
            nl = fits.getdata(join(self.pathtmp,'nl-bpm.fits'))
            uid.append(fits.getheader(join(self.pathtmp,'nl-bpm.fits'))['UNIQUEID'])
            gpm[nl==1]=0
        if isfile(join(self.pathtmp,'badpixels.fits')):
            print("Using badpixels")
            bp = fits.getdata(join(self.pathtmp,'badpixels.fits'))
            uid.append(fits.getheader(join(self.pathtmp,'badpixels.fits'))['UNIQUEID'])
            gpm[bp==1]=0
        hdu = fits.PrimaryHDU(data=gpm)
        hdu.header['DATE'] = (datetime.now().strftime(self.time_fmt),"Creation date")
        hdu.header['UNIQUEID'] = (datetime.now().strftime("%y%m%d%H%M%S"),"Unique identification number")
        hdu.header.add_comment("Built using:")
        for ui in uid:
            hdu.header.add_comment(ui)
        hdu.writeto(join(self.path,'gpm.fits'),overwrite=True)
if __name__=='__main__':
    
       bpm = badpixelmap()     
       print("Computing hot pixels")
       #bpm.hotpixels()
       print("Computing bad pixels")
       bpm.badpixels()
       print("Creating nl bpm")
       #bpm.nonlinbpm()
       
       bpm.gpm()
       
       
       
    