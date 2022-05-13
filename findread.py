#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 09:55:06 2021

The only meaningfull method is find_reads

Given a unique ID number (str) it will return
a list of reads.


@author: noboru
"""
from datetime import datetime,timedelta
from os.path import join,isfile,isdir
from os import listdir
from natsort import natsorted
from astropy.io import fits

class findread():
    def __init__(self):
        self.datapath = '/nirps_raw/read'
    def is_file(self,p):
        return isfile(join(p,'NIRPS_R01_R01.fits'))
    def reduce_list(self,ls,uid):
        start,stop = self.init_dt(uid)
        nls = [];
        for i in ls:
            if i[-1]=='/':
                i = i[:-1]
            date = int(i.split('/')[-1])
            if date>=start and date<=stop:
                nls.append(i)
        return nls
    def make_list(self,ramp,p):
        ls = [join(p,f) for f in listdir(p) if 'NIRPS_R%2.2d_R'%ramp in f]
        ls = natsorted(ls)
        return ls
    def init_dt(self,uid,ramp=False):
        y = int(uid[:4])
        m = int(uid[4:6])
        d = int(uid[6:8])
        H = int(uid[8:10])
        M = int(uid[10:12])
        S = int(uid[12:14])
        dt = datetime(y,m,d,H,M,S)
        if ramp:
            delta = timedelta(days=1)
            dt1 = dt-delta
            dt2 = dt+delta
            return [dt1.strftime("r%Y%m%d"),dt.strftime("r%Y%m%d"),dt2.strftime("r%Y%m%d")]
        delta = timedelta(hours=1)
        dt1 = dt-delta
        dt2 = dt+delta
        start = dt1.strftime("%Y%m%d%H%M%S")
        stop = dt2.strftime("%Y%m%d%H%M%S")
        return int(start),int(stop)
    def findreads(self,ID):
    
        ls = [join(self.datapath,f) for f in listdir(self.datapath) if all([isdir(join(self.datapath,f)),self.is_file(join(self.datapath,f))]) ]
    
        ls = self.reduce_list(ls, ID)
        for l in ls:
            i = 1
            while(isfile(join(l,'NIRPS_R%2.2d_R01.fits'%i))):
                if ID in fits.getheader(join(l,'NIRPS_R%2.2d_R01.fits'%i))['SEQID']:
                    return self.make_list(i,l)
                i+=1
        return []
    
if '__main__' in __name__:
    fr = findread()
    ls = fr.findreads("20210403135159")
    for l in ls:
        print(l)