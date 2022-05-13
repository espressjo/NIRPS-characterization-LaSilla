#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  7 16:08:29 2021

@author: noboru
"""

#capacitive coupling
from hxrg.type.ramp import hxramp
from numpy import zeros,argwhere
from astropy.io import fits
from os.path import join
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
cpath = '/nirps_raw/characterization'
class ccoupling():
    def __init__(self):
        self.data = zeros((4096,4096))
        self.candidate = zeros((4096,4096))
    def get_data(self,uid):
        from findramp import findramp as fr
        fr = fr()
        f = fr.findramp(uid)
        R = hxramp()
        R(f)
        R.unfold()
        self.data = R.a
        return 
    def cap_coupling(self):
        
        #mask = zeros((4096,4096),dtype=bool)
        #mask[(self.data >10) & (self.data<12)]=True
        w = argwhere((self.data >10) & (self.data<12))
        w = [c for c in w if all([c[0]>200,c[0]<3800,c[1]>200,c[1]<3800])]
        for c in w:
            print(c[0],c[1])
        
        cube =([np.asarray(self.data[coord[0]-10:coord[0]+11,coord[1]-10:coord[1]+11]) for coord in w ])
        cube = [im for im in cube if all([len(im[im>2])<2,len(im[im>13])==0])]
        print("There is %d hot pixels in this ramp."%len(cube))
        cube = [ (im/np.nanmedian(im))-1 for im in cube]
        cube = [im/np.sum(im) for im in cube]
        return np.nanmedian(cube,axis=0)        
        #
        
        #hdu = fits.PrimaryHDU(data=np.asarray(cube))
        #hdu.writeto("/home/jonathan/c11.fits",overwrite=True)
        
        #hdu = fits.PrimaryHDU(data=np.asarray(np.nanmedian(cube,axis=0)))
        #hdu.writeto("/home/jonathan/c12.fits",overwrite=True)
    def plot(self,array,color='y',alpha=0.9):
        """
        **Description**:
            This function generate a 3D plot from a 5x5 array. The output of stackHotPixel should be used as input of plotIPC.
        **color**:
            color of the plot
        **alpha**:
            transparancy of the bar plot
            """
        from matplotlib import pyplot as plt
        import seaborn as sns
        sns.set_theme()
        A = np.array(array)*100
        A[A<.1]=0
        bkp = A[2,2]
        A[2,2]=6
        fig = plt.figure(figsize=(12,12))
        ax = fig.add_subplot(111,projection='3d')
        x,y = np.meshgrid(range(5),range(5))
        z = np.zeros(len(A.ravel()))
        dx = np.ones(len(x.ravel()))-0.1
        dy = np.ones(len(y.ravel()))-0.1
        ax.bar3d(x.ravel(),y.ravel(),z,dx,dy,A.ravel(),alpha=0.7)
        ax.set_zlim([0,8])
        for i in range(5):
            for ii in range(5):
                if np.any([ii!=2,i!=2]):
                    if A[ii,i]<0.001:
                        ax.text(i+0.4,ii+0.3,A[ii,i],"%d%%"%float(A[ii,i]),c='w',fontsize=8,alpha=0.9)
                    else:
                        ax.text(i+0.3,ii+0.3,A[ii,i],"%2.1f%%"%float(A[ii,i]),c='w',fontsize=8,alpha=0.9)
                else:
                    ax.text(i+0.35,ii+0.35,A[ii,i],"%d%%"%float(bkp),c='w',fontsize=8)
        ax.set_xlim([0,5])
        ax.set_ylim([0,5])
        ax.set_xlabel('X-axis')
        ax.set_ylabel('Y-axis')
        ax.grid(False)
        ax.zaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        ax.xaxis.set_ticklabels([])
        ax.xaxis.set_ticks([],[])
        ax.yaxis.set_ticks([],[])
        ax.zaxis.set_ticks([],[])
        for line in ax.zaxis.get_ticklines():
            line.set_visible(False)
        plt.tight_layout()
        plt.show()
        return fig
        #ax.get_zaxis().set_visible(False)
        
        #ax.set_axis_off()       
        
        
if '__main__' in __name__:
    lsuids =    ['20210406204047',
         '20210406205100',
         '20210406210113',
         '20210406211126',
         '20210406212139',
         '20210406213152',
         '20210406214205',
         '20210406215218',
         '20210406220231',
         '20210406221245',
         '20210406222258',
         '20210406223311',
         '20210406224324',
         '20210406225337',
         '20210406230350',
         '20210406231403',
         '20210406232416',
         '20210406233429',
         '20210406234442',
         '20210406235455',
         '20210407000508',
         '20210407001521',
         '20210407002534',
         '20210407003547',
         '20210407004600',
         '20210407005613',
         '20210407010626',
         '20210407011639',
         '20210407012652',
         '20210407013706',
         '20210407014719',
         '20210407015732',
         '20210407020745',
         '20210407021758',
         '20210407022811',
         '20210407023824',
         '20210407024837',
         '20210407025850',
         '20210407030903',
         '20210407031916',
         '20210407032929',
         '20210407033942',
         '20210407034955',
         '20210407040008',
         '20210407041021',
         '20210407042034',
         '20210407043047',
         '20210407044100',
         '20210407045113',
         '20210407050127',
         '20210407051140',
         '20210407052153',
         '20210407053206',
         '20210407054219',
         '20210407055232',
         '20210407060245',
         '20210407061258',
         '20210407062311',
         '20210407063324',
         '20210407064337']#data for NIRPS
    cmaps = [];
    from tqdm import tqdm
    for uid in tqdm(lsuids):
        
        coupling = ccoupling()
        coupling.get_data(uid)
        cmaps.append(coupling.cap_coupling())
    capmap = np.nanmedian(cmaps,axis=0)
    hdu = fits.PrimaryHDU(data=capmap)
    hdu.header.add_comment("built using: ")
    from datetime import datetime
    fmt = "%Y-%m-%dT%H:%M:%S"
    fmtuid = "%Y%m%d%H%M%S"
    hdu.header['DATE'] = datetime.now().strftime(fmt)
    hdu.header['UID'] = datetime.now().strftime(fmtuid)
    
    for uid in lsuids:
        hdu.header.add_comment(uid)
    hdu.writeto(join(cpath,'capacitive_coupling_model.fits'),overwrite=True)
        
    _c = ((1.0-capmap[10,10]) +1.)**2 -1.
    print("Gain will be off by %.2f %%"%(_c*100))
    fig = coupling.plot(capmap[8:13,8:13])
    fig.savefig(join(cpath,'capacitive_coupling.png'))
    