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
from astropy.io import fits

class readoutnoise():
    def __init__(self):
        self.path = '/nirps_raw/characterization'
        self.pathtmp = '/nirps_raw/characterization/.tmp'
        self.ro_time = 5.24288
        self.N = 1#number of read
        self.gain = 1.33
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
            
    def make(self,uid):
        ls = fr.findreads(uid)
        r = hxread()
        R = hxramp(toponly=self.toponly)
        r(ls[0])
        R<<r
        per_amp_ro = [];
        
        for f in tqdm(ls[1:]):
            r(f)
            R<<r
            R.fit()
            per_amp_ro.append(self.compute_ro(R.a*self.N*self.ro_time))
            self.N+=1
        if self.toponly:
            np.save(join(self.pathtmp,"%s.readout.to"%uid),np.asarray(per_amp_ro))
        else:
            np.save(join(self.pathtmp,"%s.readout"%uid),np.asarray(per_amp_ro))
        return np.asarray(per_amp_ro)
    def check_cds_noise(self,lsuid):
        stats = [];
        for uid in lsuid:
            if self.toponly:
                stats.append(sc(np.load(join(self.pathtmp,"%s.readout.to.npy"%uid))[0,:]))
            else:
                stats.append(sc(np.load(join(self.pathtmp,"%s.readout.npy"%uid))[0,:]))
        mn,_,_ = zip(*stats)
        print(mn)
        _m,_,_s = sc(mn)
        print("CDS readout noise is %.2f +/- %.2f"%(_m*self.gain,_s*self.gain))
    def _single(self,uid):
        if self.toponly:
            if not isfile(join(self.pathtmp,'%s.readout.to.npy'%uid)):
                print("%s not found"%uid)
                return [],[]
            arr = np.load(join(self.pathtmp,'%s.readout.to.npy'%uid))
        else:
            if not isfile(join(self.pathtmp,'%s.readout.npy'%uid)):
                print("%s not found"%uid)
                return [],[]
            arr = np.load(join(self.pathtmp,'%s.readout.npy'%uid))    
        stats = [sc(r) for r in arr]
        S,_,err = zip(*stats)
        return np.asarray(S)*self.gain,np.asarray(err)*self.gain
    def make_graph(self,uids):
        from matplotlib import pyplot as plt
        import seaborn as sns
        sns.set_theme()
        fig,ax = plt.subplots(2,1,figsize=(12,9))
        self.toponly = True
        for uid in uids:
            y,yerr = self._single(uid)
            if len(y)<1:
                continue
            x = np.arange(2,len(y)+2,1)
            #ax.errorbar(x,y,yerr=yerr,fmt='')
            ax[0].plot(x,y,markersize=2,label=uid)
        self.toponly = False
        for uid in uids:
            y,yerr = self._single(uid)
            if len(y)<1:
                continue
            x = np.arange(2,len(y)+2,1)
            #ax.errorbar(x,y,yerr=yerr,fmt='')
            ax[1].plot(x,y,markersize=2,label=uid)
        ax[0].set_xlim([2,len(y)])
        ax[0].set(title='Readout noise of %d ramps [top only]'%(len(uids)),ylabel='Readout noise (electron)')
        ax[1].set_xlim([2,len(y)])
        ax[0].set_ylim(ax[1].get_ylim())
        ax[1].set(title='Readout noise of %d ramps'%(len(uids)),xlabel='Read number',ylabel='Readout noise (electron)')
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
    def plot_comparison(self,uid,noshow=False):
        from matplotlib import pyplot as plt
        import seaborn as sns
        sns.set_theme()
        fig,ax = plt.subplots()
        self.toponly=True
        yto,yerrto = self._single(uid)
        self.toponly=False
        y,yerr = self._single(uid)
        if len(y)<1 or len(yto)<1:
                print("something is wrong!")
                exit(0)
        x = np.arange(2,len(y)+2,1)
        ax.errorbar(x,y,yerr=yerr,fmt='',c='r')
        ax.errorbar(x,yto,yerr=yerrto,c='b',fmt='')
        ax.plot(x,y,'d',c='r',markersize=2,label="top/bottom ref. px.")
        ax.plot(x,yto,'d',c='b',markersize=2,label="top ref. px.")
        ax.set(title='Readout noise (uid: %s)'%uid,xlabel='Read number',ylabel='Readout noise (electron)')
        ax.set_xlim([2,len(y)])
        ax.legend()
        fig.savefig(join(self.path,"%s.comp.png"%(uid)))
        if noshow:
            return
        
        plt.show()
    def plot_rapport(self,uid,noshow=False):
        from matplotlib import pyplot as plt
        import seaborn as sns
        sns.set_theme()
        fig,ax = plt.subplots()
        self.toponly=True
        yto,yerrto = self._single(uid)

        if len(yto)<1:
                print("something is wrong!")
                exit(0)
        x = np.arange(2,len(yto)+2,1)
        a = yto[1]/(3**(-1.0/2.0))
        y2 = [a*i**(-1.0/2.0) for i in x]
        ax.plot(x,y2,'--',markersize=2,label="N$^{-1/2}$ decay")
        ax.errorbar(x,yto,yerr=yerrto,c='b',fmt='')
        ax.plot(x,yto,'d',c='b',markersize=2,label="NIRPS [COPL]")
        ax.set(title='Readout noise (uid: %s)'%uid,xlabel='Read number',ylabel='Readout noise (electron)')
        ax.set_xlim([2,len(yto)])
        ax.legend()
        fig.savefig(join(self.path,"%s.rapport.png"%(uid)))
        if noshow:
            return
        
        plt.show()
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
    fr = fr()
    '''
    for uid in lsuids:
        print("Working on %s"%uid)
        ro_noise = readoutnoise()
        
        ro_array = ro_noise.make(uid)
    '''
    #to plot comparaison between top and top/bottom ref. px methods
    ro_noise = readoutnoise()  
    
    #ro_noise.plot_comparison(lsuids[0],noshow=True)
    #ro_noise.plot_comparison(lsuids[1])
    #ro_noise.make_graph(lsuids)
    ro_noise.plot_rapport(lsuids[0])
    ro_noise.check_cds_noise(lsuids)
    #ro_noise.plot(lsuids)
    #to plot the results
    #ro_noise = readoutnoise()   
    #ro_noise.plot(lsuids)
        
        
                
        