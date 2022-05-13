#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 13 09:05:58 2021

@author: noboru
"""
from findread import findread as fr
from hxrg.type.read import hxread 
from hxrg.type.ramp import hxramp
from astropy.stats import sigma_clipped_stats as sc
import numpy as np
from os.path import join,isfile
from tqdm import tqdm
from os.path import basename
from astropy.io import fits

class readoutnoise():
    def __init__(self):

        self.path = '/nirps_raw/nirps/characterization'
        self.pathtmp = '/nirps_raw/nirps/characterization/.tmp'
        self.p_data = "/nirps_raw/nirps/reads"
        self.ro_time = 5.24288
        self.N = 1#number of read
        self.gain = 1.27
        self.gpm = np.asarray(fits.getdata(join(self.path,'gpm.fits')),dtype=bool) 
        #We don't want to use the reference pixels to calculate the gain
        self.gpm[:4,:]=False
        self.gpm[-4:,:]=False
        self.gpm[:,:4]=False
        self.gpm[:,-4:]=False
        self.toponly = True#use top ref.px only
    def _amp_ro(self,amp,mask):
        return sc(amp[mask],sigma=5)[2]
    def compute_ro(self,cds):
        return np.asarray([self._amp_ro(cds[:,i*128:(i+1)*128],self.gpm[:,i*128:(i+1)*128]) for i in range(32)])
            
    def make(self,fname):
        if '/' not in fname:
            fname = join(self.p_data,fname)
        cube = fits.getdata(fname)
        H = fits.getheader(fname)
        r = hxread()
        R = hxramp(toponly=self.toponly)
        r(np.asarray(cube[0],dtype=float))
        r.H = H
        R<<r
        per_amp_ro = [];
        
        for i in tqdm(range(cube.shape[0]-1)):
            r(cube[i])
            r.H = H
            R<<r
            R.fit()
            per_amp_ro.append(self.compute_ro(R.a*self.N*self.ro_time))
            self.N+=1
        if self.toponly:
            np.save(join(self.pathtmp,"%s.readout.to"%basename(fname)),np.asarray(per_amp_ro))
        else:
            np.save(join(self.pathtmp,"%s.readout"%basename(fname)),np.asarray(per_amp_ro))
        return np.asarray(per_amp_ro)
    def check_cds_noise(self,fnames):
        stats = [];
        for fname in fnames:
            if self.toponly:
                stats.append(sc(np.load(join(self.pathtmp,"%s.readout.to.npy"%fname))[0,:]))
            else:
                stats.append(sc(np.load(join(self.pathtmp,"%s.readout.npy"%fname))[0,:]))
        mn,_,_ = zip(*stats)
        print(mn)
        _m,_,_s = sc(mn)
        print("CDS readout noise is %.2f +/- %.2f"%(_m*self.gain,_s*self.gain))
    def _single(self,fname):
        if self.toponly:
            if not isfile(join(self.pathtmp,'%s.readout.to.npy'%fname)):
                print("%s not found"%fname)
                return [],[]
            arr = np.load(join(self.pathtmp,'%s.readout.to.npy'%fname))
        else:
            if not isfile(join(self.pathtmp,'%s.readout.npy'%fname)):
                print("%s not found"%fname)
                return [],[]
            arr = np.load(join(self.pathtmp,'%s.readout.npy'%fname))    
        stats = [sc(r) for r in arr]
        S,_,err = zip(*stats)
        return np.asarray(S)*self.gain,np.asarray(err)*self.gain
    def make_graph(self,fnames):
        from matplotlib import pyplot as plt
        import seaborn as sns
        sns.set_theme()
        fig,ax = plt.subplots(2,1,figsize=(12,9))
        self.toponly = True
        for fname in fnames:
            y,yerr = self._single(fname)
            if len(y)<1:
                continue
            x = np.arange(2,len(y)+2,1)
            #ax.errorbar(x,y,yerr=yerr,fmt='')
            ax[0].plot(x,y,markersize=2,label=fname.replace(".fits","").replace("NIRPS_",""))
        self.toponly = False
        for fname in fnames:
            y,yerr = self._single(fname)
            if len(y)<1:
                continue
            x = np.arange(2,len(y)+2,1)
            #ax.errorbar(x,y,yerr=yerr,fmt='')
            ax[1].plot(x,y,markersize=2,label=fname.replace(".fits","").replace("NIRPS_",""))
        ax[0].set_xlim([2,len(y)])
        ax[0].set(title='Readout noise of %d ramps [top only]'%(len(fnames)),ylabel='Readout noise (electron)')
        ax[1].set_xlim([2,len(y)])
        ax[0].set_ylim(ax[1].get_ylim())
        ax[1].set(title='Readout noise of %d ramps'%(len(fnames)),xlabel='Read number',ylabel='Readout noise (electron)')
        plt.tight_layout()
        fig.savefig(join(self.path,"ro_noise_all.png"))
        plt.show()
    def plot(self,uids):
        from matplotlib import pyplot as plt
        import seaborn as sns
        sns.set_theme()
        fig,ax = plt.subplots()
        for uid in uids:
            y,yerr = self._single(uid)
            if len(y)<1:
                continue
            x = np.arange(2,len(y)+2,1)
            #ax.errorbar(x,y,yerr=yerr,fmt='')
            ax.plot(x,y,'d',markersize=2,label=uid)
        more = ""
        if self.toponly:
            more = "[top only]"
        ax.set_xlim([2,len(y)])
        ax.set(title='Readout noise of %d ramps %s'%(len(uids),more),xlabel='Read number',ylabel='Readout noise (electron)')
        #ax.legend()
        #if self.toponly:
        #fig.savefig(join(self.path,"all_ro_noise.png"))
        plt.show()
    def plot_comparison(self,fname,noshow=False):
        from matplotlib import pyplot as plt
        import seaborn as sns
        sns.set_theme()
        fig,ax = plt.subplots()
        self.toponly=True
        yto,yerrto = self._single(fname)
        self.toponly=False
        y,yerr = self._single(fname)
        if len(y)<1 or len(yto)<1:
                print("something is wrong!")
                exit(0)
        x = np.arange(2,len(y)+2,1)
        ax.errorbar(x,y,yerr=yerr,fmt='',c='r')
        ax.errorbar(x,yto,yerr=yerrto,c='b',fmt='')
        ax.plot(x,y,'d',c='r',markersize=2,label="top/bottom ref. px.")
        ax.plot(x,yto,'d',c='b',markersize=2,label="top ref. px.")
        ax.set(title='Readout noise (file: %s)'%fname,xlabel='Read number',ylabel='Readout noise (electron)')
        ax.set_xlim([2,len(y)])
        ax.legend()
        fig.savefig(join(self.path,"%s.comp.png"%(fname)))
        if noshow:
            return
        
        plt.show()
    def plot_rapport(self,fname,noshow=False):
        from matplotlib import pyplot as plt
        import seaborn as sns
        sns.set_theme()
        fig,ax = plt.subplots()
        self.toponly=True
        yto,yerrto = self._single(fname)

        if len(yto)<1:
                print("something is wrong!")
                exit(0)
        x = np.arange(2,len(yto)+2,1)
        a = yto[1]/(3**(-1.0/2.0))
        y2 = [a*i**(-1.0/2.0) for i in x]
        ax.plot(x,y2,'--',markersize=2,label="N$^{-1/2}$ decay")
        ax.errorbar(x,yto,yerr=yerrto,c='b',fmt='')
        ax.plot(x,yto,'d',c='b',markersize=2,label="NIRPS [COPL]")
        ax.set(title='Readout noise (file: %s)'%fname,xlabel='Read number',ylabel='Readout noise (electron)')
        ax.set_xlim([2,len(yto)])
        ax.legend()
        fig.savefig(join(self.path,"%s.rapport.png"%(fname)))
        if noshow:
            return
        
        plt.show()
if '__main__' in __name__:
    
    darks = ["NIRPS_2022-05-05T15_31_01_590.fits",
            "NIRPS_2022-05-05T15_40_24_489.fits",
            "NIRPS_2022-05-05T15_49_47_388.fits",
            "NIRPS_2022-05-05T15_59_10_285.fits"]

    '''
    for fname in darks:
        print("Working on %s"%fname)
        ro_noise = readoutnoise()
        ro_noise.toponly = False
        ro_array = ro_noise.make(fname)
    for fname in darks:
        print("Working on %s"%fname)
        ro_noise = readoutnoise()
        ro_noise.toponly = True#default setup
        ro_array = ro_noise.make(fname)
    '''
    #to plot comparaison between top and top/bottom ref. px methods
    ro_noise = readoutnoise()  
    
    #ro_noise.plot_comparison(lsuids[0],noshow=True)
    #ro_noise.plot_comparison(lsuids[1])
    ro_noise.make_graph(darks)
    ro_noise.plot_rapport(darks[0])
    ro_noise.check_cds_noise(darks)
    #ro_noise.plot(lsuids)
    #to plot the results
    #ro_noise = readoutnoise()   
    #ro_noise.plot(lsuids)
        
        
                
        