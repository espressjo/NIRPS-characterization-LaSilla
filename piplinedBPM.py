#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 11:00:56 2021

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
from scipy.ndimage import zoom
class pipelinedBPM():
    def __init__(self):
        self.path = '/nirps_raw/nirps/characterization/check-gpm'
        self.pathtmp = '/nirps_raw/nirps/.tmp/check-gpm'
        self.nonlin = 'nonlin.fits'
        self.time_fmt = "%Y-%m-%dT%H:%M:%S"
        self.ls_dark = [];
        self.ls_flat = [];
        self.p_data = "/nirps_raw/nirps/read"
        self.darks = [ "NIRPS_2022-05-05T15_31_01_590.fits",
                      "NIRPS_2022-05-05T15_40_24_489.fits",
                      "NIRPS_2022-05-05T15_49_47_388.fits",
                      "NIRPS_2022-05-05T15_59_10_285.fits"]
       
        
        self.flats= [  "NIRPS_2022-05-05T16_09_12_198.fits",
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
    def _make_cds_data(self,f):
        '''
        For NIRPS we want to avoid using the bottom ref. pix. and 
        we do not want to use the odd-even option.
        Parameters
        ----------
        f: file (cube of data)
        -------
        TYPE
            CDS np.array of float with corrected ref.px. using side and top.
        '''
        arr = fits.getdata(join(self.p_data,f))[:2,:,:]
        return rpc(np.asarray(arr[1],dtype=float)-np.asarray(arr[0],dtype=float)).refpxcorrtop(False)
    def _make_cds(self,fr1,fr2):
        '''
        For NIRPS we want to avoid using the bottom ref. pix. and 
        we do not want to use the odd-even option.
        Parameters
        ----------
        fr1 : STR
            File name of read 1
        fr2 : STR
            File name of read N
        Returns
        -------
        TYPE
            CDS np.array of float with corrected ref.px. using side and top.
        '''
        return rpc(np.asarray(fits.getdata(fr2),dtype=float)-np.asarray(fits.getdata(fr1),dtype=float)).refpxcorrtop(False)
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
        for uid in self.ls_dark:
            hdu.header.add_comment(uid)
        hdu.writeto(join(self.pathtmp,"hotpixels.fits"),overwrite=True)
        print("hotpixels.fits write successfully.")
        return hp
    def badpixels_new(self):
        fr = findramp()
        freads = findread()
        #check slope on 1st ramp
       
        cube_bpm = [];
        for uid in self.ls_flat:
            fname_ramp = fr.findramp(uid)
            slope = fits.open(fname_ramp)['slope'].data
            mean,med,std = sc(slope[1800:2200,1800:2200],sigma=5)
            n = int(40000./(mean*5.5733))
            
            ls = freads.findreads(uid)
            if len(ls)-1 <n:
                print('Skipping %s since f<40,000ADU.'%uid)   
            else:
                cube_bpm.append(self._make_cds(ls[0], ls[n]))
                
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
        mask[im<0.7] = 1
        mask[im>1.3] = 1
        mask[:4,:]=0
        mask[-4:,:] = 0
        mask[:,:4] = 0
        mask[:,-4:] = 0

        hdu = fits.PrimaryHDU(data=mask)
        hdu.header['DATE'] = (datetime.now().strftime(self.time_fmt),"Creation date")
        hdu.header['UNIQUEID'] = (datetime.now().strftime("%y%m%d%H%M%S"),"Unique identification number")
        hdu.header.add_comment("Built using:")
        for uid in self.ls_flat:
            hdu.header.add_comment(uid)
        hdu.writeto(join(self.pathtmp,"badpixels.fits"),overwrite=True)
    def badpixels(self):
        fr = findramp()
        freads = findread()
        #check slope on 1st ramp
       
        cube_bpm = [];
        for uid in self.ls_flat:
            fname_ramp = fr.findramp(uid)
            slope = fits.open(fname_ramp)['slope'].data
            mean,med,std = sc(slope[1800:2200,1800:2200],sigma=5)
            n = int(40000./(mean*5.5733))
            
            ls = freads.findreads(uid)
            if len(ls)-1 <n:
                print('Skipping %s since f<40,000ADU.')   
            else:
                cube_bpm.append(self._make_cds(ls[0], ls[n]))
        im = np.median(cube_bpm,axis=0)
        filt =    median_filter(im,size=(20,20))
        filt/=np.median(filt)
        im/=filt
        mean,_,std = sc(im,sigma=5)
        bpm = np.zeros((4096,4096))
        bpm[im<(mean-5*std)]=1
        bpm[:4,:]=0
        bpm[-4:,:] = 0
        bpm[:,:4] = 0
        bpm[:,-4:] = 0
        
        _p = [128*i for i in range(32) if all([i!=0,i%2==0])]
        for p in _p:
            bpm[:,p-1]=0 
            bpm[:,p] = 0
            mn,_,st = sc(im[4:-4,p-1:p+1])
            _s1 = bpm[:,p-1:p+1]
            _s2 = im[:,p-1:p+1]
            _s1[_s2<(mn-5*st)]=1
            bpm[:,p-1:p+1] = _s1
        
        hdu = fits.PrimaryHDU(data=bpm)
        hdu.header['DATE'] = (datetime.now().strftime(self.time_fmt),"Creation date")
        hdu.header['UNIQUEID'] = (datetime.now().strftime("%y%m%d%H%M%S"),"Unique identification number")
        hdu.header.add_comment("Built using:")
        for uid in self.ls_flat:
            hdu.header.add_comment(uid)
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
        hdu.writeto(join(self.pathtmp,'nl-bpm.fits'),overwrite=True)
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
    def change_value(self,array,value=256):
        arr = np.zeros(array.shape)
        arr[array==1]=value
        return arr
        
    def make_eso_cmpl_file(self,withnl = False,ext=''):
        hdu = fits.open(join(self.pathtmp,"badpixels.fits"))
        if withnl:
            nl = fits.getdata(join(self.path,'nonlin.fits'))[2,:,:]
            hdu[0].data[nl==1]=1
        
        H = hdu[0].header
        H['HIERARCH ESO INS MODE'] = 'HA'
        H['HIERARCH ESO PRO CATG'] = 'BAD_PIXEL_MASK'
        H['HIERARCH ESO DET WIN1 BINX'] = 1 
        H['HIERARCH ESO DET WIN1 BINY'] = 1 
        H['HIERARCH ESO DET WIN1 NX'] = 4088 
        H['HIERARCH ESO DET WIN1 NY'] = 4088 
        H['MJD-OBS'] = 59215#59335
        H['HIERARCH ESO QC EXT0 ROX0 ROY0 CONAD'] = 1/1.27
        #bad pixel map en HA
        Phdu = fits.PrimaryHDU(header=H)
        #imext = fits.ImageHDU(self.change_value(np.rot90(hdu[0].data[4:-4,4:-4],2),value=8192))
        imext = fits.ImageHDU(self.change_value(hdu[0].data[4:-4,4:-4],value=8192))
        hdul = fits.HDUList([Phdu,imext])
        if withnl:
            hdul.writeto(join(self.path,"badpixels_HA.nl.fits%s"%ext),overwrite=True)
        else:
            hdul.writeto(join(self.path,"badpixels_HA.fits%s"%ext),overwrite=True)
        #bad pixel map en HE
        H['HIERARCH ESO INS MODE'] = 'HE'
        Phdu = fits.PrimaryHDU(header=H)
        #imext = fits.ImageHDU(self.change_value(np.rot90(hdu[0].data[4:-4,4:-4],2),value=8192))
        imext = fits.ImageHDU(self.change_value(hdu[0].data[4:-4,4:-4],value=8192))
        hdul = fits.HDUList([Phdu,imext])
        hdul.writeto(join(self.path,"badpixels_HE.fits%s"%ext),overwrite=True)
        
        
        
        #hot pixels
        hdu = fits.open(join(self.pathtmp,"hotpixels.fits"))
        H = hdu[0].header
        H['HIERARCH ESO INS MODE'] = 'HA'
        H['HIERARCH ESO PRO CATG'] = 'HOT_PIXEL_MASK'
        H['HIERARCH ESO DET WIN1 BINX'] = 1 
        H['HIERARCH ESO DET WIN1 BINY'] = 1 
        H['HIERARCH ESO DET WIN1 NX'] = 4088 
        H['HIERARCH ESO DET WIN1 NY'] = 4088 
        H['MJD-OBS'] = 59215#59335
        
        H['HIERARCH ESO DET OUT1 GAIN'] = 1.27
        Phdu = fits.PrimaryHDU(header=H)
        #imext = fits.ImageHDU(self.change_value(np.rot90(hdu[0].data[4:-4,4:-4],2),value=256))
        imext = fits.ImageHDU(self.change_value(hdu[0].data[4:-4,4:-4],value=256))
        hdul = fits.HDUList([Phdu,imext])
        hdul.writeto(join(self.path,"hotpixels_HA.fits%s"%ext),overwrite=True)
        #hot pixel map en HE
        H['HIERARCH ESO INS MODE'] = 'HE'
        Phdu = fits.PrimaryHDU(header=H)
        
        #imext = fits.ImageHDU(self.change_value(np.rot90(hdu[0].data[4:-4,4:-4],2),value=256))
        imext = fits.ImageHDU(self.change_value(hdu[0].data[4:-4,4:-4],value=256))
        hdul = fits.HDUList([Phdu,imext])
        hdul.writeto(join(self.path,"hotpixels_HE.fits%s"%ext),overwrite=True)
        print("done")
        
if __name__=='__main__':
    
    bpm = pipelinedBPM()   
    bpm.ls_dark = bpm.ls_dark_april
    bpm.ls_flat = bpm.ls_flat_april
    print("Computing hot pixels")
    bpm.hotpixels()
    print("Computing bad pixels")
    bpm.badpixels_new()
    #print("Creating nl bpm")
    #bpm.nonlinbpm()
    bpm.make_eso_cmpl_file(withnl=False,ext='.april')
    bpm.make_eso_cmpl_file(withnl=True,ext='.april')
    #bpm.gpm()
    #now for jan2april
    bpm = pipelinedBPM()   
    bpm.ls_dark = bpm.ls_dark_jan2april
    bpm.ls_flat = bpm.ls_flat_jan2april
    print("Computing hot pixels")
    bpm.hotpixels()
    print("Computing bad pixels")
    bpm.badpixels_new()
    #print("Creating nl bpm")
    #bpm.nonlinbpm()
    bpm.make_eso_cmpl_file(withnl=False,ext='.jan2april')
    bpm.make_eso_cmpl_file(withnl=True,ext='.jan2april')
    