#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 13 10:59:20 2021

@author: noboru
"""
from astropy.io import fits
class fix_bias():
    def __init__(self,fname):
        self.hdu = fits.open(fname)
        
        self.H = self.fix_header(self.hdu)
        if 'y' in input("Do you want to update the file?[yes/no]"):
            self.hdu.header = self.H
            self.hdu.writeto(fname,overwrite=True)
        else:
            if 'y' in input("Do you want to rename the file?[yes/no]"):
                new_name = input("Enter new name?:\n")
                from os import getcwd
                p = getcwd()
                from os.path import join
                self.hdu.header = self.H
                self.hdu.writeto(join(p,new_name),overwrite=True)

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
    fix_bias(fname)