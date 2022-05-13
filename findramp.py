#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 10:15:54 2021

@author: noboru
"""
from os.path import isdir,join
from datetime import datetime,timedelta
from os import listdir
from astropy.io import fits

class findramp():
    def __init__(self):
        self.datapath = '/nirps_raw/'
    def _kinit(self,uid):
        y = int(uid[:4])
        m = int(uid[4:6])
        d = int(uid[6:8])
        H = int(uid[8:10])
        M = int(uid[10:12])
        S = int(uid[12:14])
        dt = datetime(y,m,d,H,M,S)
        delta = timedelta(days=1)
        dt1 = dt-delta
        dt2 = dt+delta
        ls = [];
        if isdir(join(self.datapath,dt1.strftime("r%Y%m%d"))):
            ls.append(join(self.datapath,dt1.strftime("r%Y%m%d")))
        if isdir(join(self.datapath,dt.strftime("r%Y%m%d"))):
            ls.append(join(self.datapath,dt.strftime("r%Y%m%d")))
        if isdir(join(self.datapath,dt2.strftime("r%Y%m%d"))):
            ls.append(join(self.datapath,dt2.strftime("r%Y%m%d")))
        return ls
        
    def findramp(self,uid):
        ls = self._kinit(uid)
        for d in ls:
            ls = [join(d,f) for f in listdir(d) if 'fits' in f]
            for fname in ls:
                try:
                    if uid==fits.getheader(fname)["UNIQUEID"]:
                        return fname
                except :
                    pass