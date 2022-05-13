#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 13 10:59:20 2021

@author: noboru
"""
from astropy.io import fits
class fix_nonlin():
    def __init__(self,fname):
        hdu = fits.open(fname)
        self.hdul = self.fix_ext(hdu)
        self.H = self.fix_header(hdu)
        if 'y' in input("Do you want to update the file?[yes/no]"):
            self.hdul[0].header = self.H
            #set ref. px. of C2 and C3 to zero
            self.hdul[0].data = self.set_ref_px_to_zero(self.hdul[0].data)
            self.hdul[1].data = self.set_ref_px_to_zero(self.hdul[1].data)
            self.hdul.writeto(fname,overwrite=True)
        else:
            if 'y' in input("Do you want to rename the file?[yes/no]"):
                new_name = input("Enter new name?:\n")
                from os import getcwd
                p = getcwd()
                from os.path import join
                self.hdul[0].header = self.H
                #set ref. px. of C2 and C3 to zero
                self.hdul[0].data = self.set_ref_px_to_zero(self.hdul[0].data)
                self.hdul[1].data = self.set_ref_px_to_zero(self.hdul[1].data)
                self.hdul.writeto(join(p,new_name),overwrite=True)
    def set_ref_px_to_zero(self,im):
        im[:4,:]=0
        im[-4:,:]=0
        im[:,:4]=0
        im[:,-4:]=0
        return im
    def fix_ext(self,hdu):
        if len(hdu)==1:
            if hdu[0].data.shape ==(3,4096,4096):
                c3 = hdu[0].data[0]
                c2= hdu[0].data[1]
                err = hdu[0].data[2]
                return fits.HDUList([fits.PrimaryHDU(data=c3),fits.ImageHDU(data=c2),fits.ImageHDU(data=err)])
            else:
                print("Cannot understand the file format")
        elif len(hdu)==3:
            return hdu
    def fix_header(self,hdu):
        if 'UID' in hdu[0].header:
            return hdu[0].header
        if 'UNIQUEID' in hdu[0].header:
            hdu[0].header['UID'] = hdu[0].header['UNIQUEID']
            del hdu[0].header['UNIQUEID']
        else:
            from datetime import datetime
            hdu[0].header['UID'] = datetime.now().strftime("%y%m%d%H%M%S")
        return hdu[0].header
if '__main__' in __name__:
    from sys import argv
    if len(argv)!=2:
        print("fixnonfit.py name-of-file-to-fix.fits")
        exit(0)
    fname = argv[1]
    fix_nonlin(fname)