#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  4 11:53:16 2021

@author: noboru
"""
from astropy.io import fits
import numpy as np
from scipy.ndimage import zoom,median_filter
from tqdm import tqdm

def combine(ls):
    '''
    You do not need this bit of code. It is a simple median combine of all the images.

    '''
    if ls[0].isdigit():
        #this must be a list of unique ID
        from findramp import findramp
        fr = findramp()
        ls = [fr.findramp(f) for f in ls]        
    print("Combining images...")
    cube = [fits.open(f)[1].data for f in tqdm(ls)]
    return np.median(np.asarray(cube),axis=0)    

def hotpixels(im,pix_bin = 64,low=-0.3,hi=4):

    #:::::: Preprocessing ::::::
    '''
    This is the preprocessing that is already done before running the bad pixel routine in the pipline
    '''
    for i in range(4096):
        im[i,:]-=np.nanmedian(im[i,:])
     # first pixel of each amplifier
    first_col_x = (np.append(np.arange(16)*256,np.arange(16)*256-1+256))
    first_col_x.sort()
    # median-filter the first and last ref pixel, which trace the
    # behavior of the first-col per amp.
    amp0 = (median_filter(im[:,0],7) + median_filter(im[:,-1],7))/2
    for i in first_col_x: # subtract first column behavior
        im[:,i] -= amp0
    
    for i in range(4096):
        im[:,i]-=np.nanmedian(im[:,i])
    #:::::: Preprocessing END ::::::
    
    #::::::: Smoothing the image #:::::::
    '''
    this is very similar to the smoothing median filter described in the HDRL manual
    '''
    nbin = 4096//pix_bin
    percentile = 50  
    box = np.zeros([nbin,nbin])
    for i in (range(nbin)):
        for j in range(nbin):
            # -1 sigma of distrubution. Should be good to remove illuminated pixels
            box[i,j] = np.nanpercentile(im[i*pix_bin:i*pix_bin+pix_bin,j*pix_bin:j*pix_bin+pix_bin],percentile)    
    lowf = zoom(box,pix_bin)
    im-=lowf
    #::::::: Smoothing the image end :::::::    
        
    #hard coded limit
    '''
        I believe there is a way to set specific limits in the HDRL
    '''
    mask = np.zeros_like(im)
    mask[im>hi] = 256# should be aroung 4. We consider a pixel good as long as its hotnest is much lower than the typicall target flux
    mask[im<low] = 256#should be around -0.3    
    #lets flag all the reference pixels as good pixels
    mask[:4,:]=0
    mask[-4:,:] = 0
    mask[:,:4] = 0
    mask[:,-4:] = 0
   
    #::::::::::::::::::::::::::::::::::::
    #make everything ESO ready
    phdu = fits.PrimaryHDU()
    image = fits.ImageHDU(data=mask[4:-4,4:-4])
    phdu.header['HIERARCH ESO INS MODE'] = 'HA'
    phdu.header['HIERARCH ESO PRO CATG'] = 'HOT_PIXEL_MASK'
    phdu.header['HIERARCH ESO DET WIN1 BINX'] = 1 
    phdu.header['HIERARCH ESO DET WIN1 BINY'] = 1 
    phdu.header['HIERARCH ESO DET WIN1 NX'] = 4088 
    phdu.header['HIERARCH ESO DET WIN1 NY'] = 4088 
    phdu.header['MJD-OBS'] = 59215#59335
    phdu.header['HIERARCH ESO DET OUT1 GAIN'] = 1.27
    hdul = fits.HDUList([phdu,image])
    hdul.writeto("hotpixels_HA.fits",overwrite=True)
    hdul[1].header['HIERARCH ESO INS MODE'] = 'HE'
    hdul.writeto("hotpixels_HE.fits",overwrite=True)
    #::::::::::::::::::::::::::::::::::::
    return
def badpixels(im,pix_bin=64,low=0.7,hi=1.3):
    #:::::: Preprocessing ::::::
    '''
    This is the preprocessing that is already done before running the bad pixel routine in the pipline
    '''
    # first pixel of each amplifier
    first_col_x = (np.append(np.arange(16)*256,np.arange(16)*256-1+256))
    first_col_x.sort()
    # median-filter the first and last ref pixel, which trace the
    # behavior of the first-col per amp.
    amp0 = (median_filter(im[:,0],7) + median_filter(im[:,-1],7))/2
    for i in first_col_x: # subtract first column behavior
        im[:,i] -= amp0
    #:::::: Preprocessing END ::::::
        
    #::::::: Smoothing the image #:::::::
    '''
    this is very similar to the smoothing median filter described in the HDRL manual
    '''
    nbin = 4096//pix_bin
    percentile = 50
    box = np.zeros([nbin,nbin])
    for i in (range(nbin)):
        for j in range(nbin):
            # -1 sigma of distrubution. Should be good to remove illuminated pixels
            box[i,j] = np.nanpercentile(im[i*pix_bin:i*pix_bin+pix_bin,j*pix_bin:j*pix_bin+pix_bin],percentile)

    lowf = zoom(box,pix_bin)
    im/=lowf #however we divide the filter to the image. Is this possible in HDRL ?    
    #We want to see the illumination of a pixel compaired with its neighbours to set hard limits
    #::::::: Smoothing the image END :::::::        
        
    #::::::::::::: create the mask ::::::::::::::
    mask = np.zeros_like(im)
    mask[im<low] = 8192
    mask[im>hi] = 8192
    mask[:4,:]=0
    mask[-4:,:] = 0
    mask[:,:4] = 0
    mask[:,-4:] = 0  
    
    #::::::::::::::::::::::::::::::::::::
        
    #:::::::: make everything ESO ready ::::::::::
    phdu = fits.PrimaryHDU()
    image = fits.ImageHDU(data=mask[4:-4,4:-4])
    hdul = fits.HDUList([phdu,image])
        
    hdul[0].header['HIERARCH ESO INS MODE'] = 'HA'
    hdul[0].header['HIERARCH ESO PRO CATG'] = 'BAD_PIXEL_MASK'
    hdul[0].header['HIERARCH ESO DET WIN1 BINX'] = 1 
    hdul[0].header['HIERARCH ESO DET WIN1 BINY'] = 1 
    hdul[0].header['HIERARCH ESO DET WIN1 NX'] = 4088 
    hdul[0].header['HIERARCH ESO DET WIN1 NY'] = 4088 
    hdul[0].header['MJD-OBS'] = 59215#59335
    hdul[0].header['HIERARCH ESO QC EXT0 ROX0 ROY0 CONAD'] = 1/1.27

    hdul.writeto("badpixels_HA.fits",overwrite=True)
    #bad pixel map en HE
    hdul[0].header['HIERARCH ESO INS MODE'] = 'HE'
    hdul.writeto("badpixels_HE.fits",overwrite=True)
    #::::::::::::::::::::::::::::::::::::
    return
    
if '__main__' in __name__:
    #ls_dark and ls_flat can be a list of unique ID or a list files e.i., /home/jonathan/NIRPS_GEN_LED_001.fits
    ls_dark = [ '20210403135159',
                '20210403140207',
                '20210403141214',
                '20210403142222',
                '20210403143229',
                '20210403144237',
                '20210403145244',
                '20210403150252',
                '20210403151259',
                '20210403152307',
                '20210403153314',
                '20210403154322',
                '20210403155329',
                '20210403160336',
                '20210403161344',
                '20210403162351',
                '20210403163359',
                '20210403164406',
                '20210403165414',
                '20210403170421']#unique ID
    ls_flat = [ '20210403171440',
                '20210403175352',
                '20210403183253',
                '20210403191153',
                '20210403195054',
                '20210403203006',
                '20210403210907',
                '20210403214818',
                '20210403222719',
                '20210403230620']#unique ID for flat
    
    dark = combine(ls_dark,name='dark.fits')
    flat = combine(ls_flat,name='flat.fits')
    print("Creating hot pixel map")
    hotpixels(dark)
    print("Creating bad pixel map")
    badpixels(flat)

    
 