#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 15:40:33 2021

@author: espressjo

Since we don't want to use the bottom reference pixels, 
we need to recreate the ramp file to compute the superbias
"""
import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt
import seaborn as sns
from os.path import basename,isfile,isdir,join
from findread import findread
from findramp import findramp

from hxrg.type.read import hxread
from hxrg.type.ramp import hxramp
from tqdm import tqdm

class superbias():
    def __init__(self):
        self.path = '/nirps_raw/characterization' #place where to store the superbias file and the graphs
        self.tmppath = '/nirps_raw/characterization/.tmp'#tmp path 
        
        
    def makeramp(self,uid):
        #given an uid, it will return a hxramp object
        fr = findread()
        framp = findramp()
        rampfile = framp.findramp(uid)
        ls = fr.findreads(uid)
        read = hxread()
        Ramp = hxramp(toponly = True)
        Ramp.H = fits.getheader(rampfile)
        print("Working on %s"%uid)
        for f in tqdm(ls):
            read(f)
            Ramp<<read
        Ramp.fit()
        return Ramp
    def extractbias(self,lsuids):
        biases = [self.makeramp(uid).b for uid in lsuids]
        return np.asarray(biases)
            
sns.set_theme()
lsuids = ['20210403135159',
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
          '20210403170421']

def side(im):
    '''
    Return the median of side ref. pixels of an image

    '''
    l = np.append(im[:,0:4].ravel(),im[:,-4:].ravel())
    return np.median(l)

if '__main__' in __name__:
    sb = superbias()
    from datetime import datetime
    fname = datetime.now().strftime("_%Y%m%d")
    time = datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
    uid = datetime.now().strftime("%y%m%d%H%M%S")
    
    if isfile(join(sb.tmppath,"tmp.bias.fits")):
        if 'yes' in input("tmp.bias.fits found, Do you want to upload \nthe file or re-create the biases?[yes/no] "):
            data = fits.getdata(join(sb.tmppath,"tmp.bias.fits"))
        else:
            data = sb.extractbias(lsuids)
    else:
        data = sb.extractbias(lsuids)
    
    hdu = fits.PrimaryHDU(data = data)
    hdu.header.add_comment("Built using uids:")
    for uid in lsuids:
        hdu.header.add_comment(uid)
    hdu.writeto(join(sb.tmppath,"tmp.bias.fits"),overwrite=True)
    
    bias_level = [side(im) for im in data]
    for i in range(len(data)):
        data[i]-=bias_level[i]
    
    #side ref. pixels bias level
    f,ax = plt.subplots(figsize=(14,14))
    name_x = [uid for uid in lsuids]
    ax.set_xticklabels(name_x,rotation=-45)
    ax.plot(name_x,bias_level,'d',markersize=3)
    ax.set(xlabel='File name',ylabel='Side ref. pixels bias level (ADU)')
    plt.tight_layout()
    
    f.savefig(join(sb.path,'bias%s.png'%fname))
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
    f.savefig(join(sb.path,'cross%s.png'%fname))

    superbias = np.median(data,axis=0)+np.median(bias_level)

    
    hdu = fits.PrimaryHDU(data = superbias)#np.uint16(superbias))
    names = [f for f in lsuids]
    for name in names:
        hdu.header.add_comment(name)
    hdu.header['DATE'] = (time,"Creation date") 
    hdu.header['UNIQUEID'] = (uid, "Unique id of this file")
    if isfile(join(sb.path,"superbias%s.fits"%fname)):
        if 'yes' in input('%s exist, do you want to overwrite it?[yes/no]'%("superbias%s.fits"%fname)):
            hdu.writeto(join(sb.path,"superbias%s.fits"%fname),overwrite=True)
    else:
        hdu.writeto(join(sb.path,"superbias%s.fits"%fname))
    plt.show()
    