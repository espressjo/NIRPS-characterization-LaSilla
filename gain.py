#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 14:19:23 2021

@author: espressjo
"""
from os.path import join,isfile
from astropy.io import fits
import numpy as np
from scipy.optimize import curve_fit
from tqdm import tqdm
from astropy.stats import sigma_clipped_stats as sc
from matplotlib import pyplot as plt
import seaborn as sns
sns.set_theme()

class gain():
    def __init__(self):
        self.path = '/nirps_raw/characterization'
        self.pathtmp = '/nirps_raw/characterization/.tmp'
        self.gpm = np.asarray(fits.getdata(join(self.path,'gpm.fits')),dtype=bool) 
        #We don't want to use the reference pixels to calculate the gain
        self.gpm[:4,:]=False
        self.gpm[-4:,:]=False
        self.gpm[:,:4]=False
        self.gpm[:,-4:]=False
        self.c3 = fits.getdata(join(self.path,'nonlin.fits'))[0]
        self.c2 = fits.getdata(join(self.path,'nonlin.fits'))[1]
        self.size = (64,64)
        self.uid1 = ''
        self.uid2 = ''
        
    def var_err(self,f,var):
        def unc(f):
            return np.sqrt(16**2+f)
        return var*2*unc(f)/f
    def makeXY(self,uid1,uid2,remakefile=False):
        self.uid1 = uid1
        self.uid2 = uid2
        if isfile(join(self.pathtmp,"%sx%s.gain.diff.fits"%(uid1,uid2))) and not remakefile:
            self.diff = fits.getdata(join(self.pathtmp,"%sx%s.gain.diff.fits"%(uid1,uid2)))
            self.S = fits.getdata(join(self.pathtmp,"%sx%s.gain.s.fits"%(uid1,uid2)))
            print("%sx%s.gain.s.fits loaded."%(uid1,uid2))
            print("%sx%s.gain.diff.fits loaded."%(uid1,uid2))
            return
        if isfile(join(self.pathtmp,"%sx%s.gain.diff.fits"%(uid2,uid1))) and not remakefile:
            self.diff = fits.getdata(join(self.pathtmp,"%sx%s.gain.diff.fits"%(uid2,uid1)))
            self.S = fits.getdata(join(self.pathtmp,"%sx%s.gain.s.fits"%(uid2,uid1)))
            print("%sx%s.gain.s.fits loaded."%(uid2,uid1))
            print("%sx%s.gain.diff.fits loaded."%(uid2,uid1))
            return
        
        print("Creating files.")
        cube1 = fits.getdata(join(self.path,uid1+'.fits'))
        cube2 = fits.getdata(join(self.path,uid2+'.fits'))
        
        self.diff = np.zeros((len(cube1),4096,4096))
        self.S = np.zeros((len(cube1),4096,4096))
        for i in tqdm(range(len(cube1))):
            self.diff[i,:,:] = self._corr_nl(cube1[i,:,:])-self._corr_nl(cube2[i,:,:])
            self.S[i,:,:] = 0.5*(self._corr_nl(cube1[i,:,:])+self._corr_nl(cube2[i,:,:]))
        hdu = fits.PrimaryHDU(data=self.diff)
        hdu.writeto(join(self.pathtmp,"%sx%s.gain.diff.fits"%(uid1,uid2)),overwrite=True)
        hdu = fits.PrimaryHDU(data=self.S)
        hdu.writeto(join(self.pathtmp,"%sx%s.gain.s.fits"%(uid1,uid2)),overwrite=True)
    def _sub3d(self,im,x,y):
        return im[:,y-int(self.size[1]/2):y+int(self.size[1]/2),x-int(self.size[0]/2):x+int(self.size[0]/2)]
    def _sub(self,im,x,y):
        return im[y-int(self.size[1]/2):y+int(self.size[1]/2),x-int(self.size[0]/2):x+int(self.size[0]/2)]
    def _sub3d_corner(self,im,x,y):
        return im[:,y:y+self.size[1],x:x+self.size[0]]
    def _sub_corner(self,im,x,y):
        return im[y:y+self.size[1],x:x+self.size[0]]
    def makesvar_center(self,x,y):
        '''
        x and y is at the center of a box of size self.size

        Parameters
        ----------
        x : TYPE
            DESCRIPTION.
        y : TYPE
            DESCRIPTION.

        Returns
        -------
        TYPE
            DESCRIPTION.

        '''
        gpm = self._sub(self.gpm,x,y)
        s = self._sub3d(self.S,x,y)
        diff = self._sub3d(self.diff,x,y)
        S = [sc(im[gpm],maxiters=4)[0] for im in s]
        var = [sc(im[gpm],maxiters=4)[2]**2/2.0 for im in diff]
        return np.asarray(S),np.asarray(var)
    def makesvar_corner(self,x,y):
        '''
        x and y is at the upper corner (0,0) of a box of size self.size

        Parameters
        ----------
        x : TYPE
            DESCRIPTION.
        y : TYPE
            DESCRIPTION.

        Returns
        -------
        TYPE
            DESCRIPTION.

        '''
        gpm = self._sub_corner(self.gpm,x,y)
        s = self._sub3d_corner(self.S,x,y)
        diff = self._sub3d(self.diff,x,y)
        S = [sc(im[gpm],maxiters=4)[0] for im in s]
        var = [sc(im[gpm],maxiters=4)[2]**2/2.0 for im in diff]
        return np.asarray(S),np.asarray(var)
    def _corr_nl(self,im):
        return im+self.c3*im**3+self.c2*im**2 
     
    def fit(self,s,var,sat=10000):
        def lin(t,a,b):
            return t*a+b
        var = var[s<sat]
        s = s[s<sat]
        var_unc = self.var_err(s,var)
        beta,cov = curve_fit(lin,s,var,sigma=var_unc,p0=[0.7,300],maxfev=10000)
        #print(beta,np.diag(cov))
        gain = 1.0/beta[0]
        gain_err,ro_err = np.sqrt(np.diag(cov))
        gain_err = gain*(gain_err/beta[0])
        ro = np.sqrt(beta[1])
        ro_err = ro*(ro_err/beta[1])
        
        return gain,gain_err,ro,ro_err
    def run(self):
        from multiprocessing import Process
        
        ps = [ Process(target=self.p_run,args=(amp,)) for amp in range(1,33)]
        for i in range(32):
            print("Launching amp #%d"%(i+1))
            ps[i].start()
        for i in range(32):
            ps[i].join()
    def p_run(self,amp):
        _map = np.zeros((4,4096,128))
        if 128%self.size[0]!=0:
            print("size [x] must be a multiple of 128")
            exit(0)
        if 4096%self.size[0]!=0:
            print("size [y] must be a multiple of 4096")
            exit(0)
        xsteps = int(128/self.size[0])
        ysteps = int(4096/self.size[1])
        
        for y in tqdm(range(ysteps)):
            for x in range(xsteps):
                try:
                    s,var = self.makesvar_corner(x*self.size[0]+(amp-1)*128,y*self.size[1])
                    gain,gainerr,ro,roerr = self.fit(s,var)
                except :
                    
                    gain=0
                    gainerr=0 
                    ro=0
                    roerr=0
                
                _map[0,y*self.size[1]:(y+1)*self.size[1],x*self.size[0]:(x+1)*self.size[0]] = gain
                _map[1,y*self.size[1]:(y+1)*self.size[1],x*self.size[0]:(x+1)*self.size[0]] = gainerr
                _map[2,y*self.size[1]:(y+1)*self.size[1],x*self.size[0]:(x+1)*self.size[0]] = ro
                _map[3,y*self.size[1]:(y+1)*self.size[1],x*self.size[0]:(x+1)*self.size[0]] = roerr
               
        hdu = fits.PrimaryHDU()
        hdu.data = _map
        hdu.writeto(join(self.pathtmp,'%sx%s.amp.%d'%(self.uid1,self.uid2,amp)),overwrite=True)
    def recombine(self):
        from os import listdir
        from os.path import basename
        from scanf import scanf
        
        ls = [join(self.pathtmp,f) for f in listdir(self.pathtmp) if '%sx%s.amp'%(self.uid1,self.uid2) in f]

        _map = np.zeros((4,4096,4096))
        for f in ls:
            amp = scanf("%s.amp.%d",basename(f))[1]
            print("Parsing amp #%.2d"%amp)
            _map[:,:,(amp-1)*128:(amp)*128] = fits.getdata(f)

        hdu = fits.PrimaryHDU()
        hdu.data = _map
        from datetime import datetime
        time = datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
        uid = datetime.now().strftime("%y%m%d%H%M%S")
        hdu.header['DATE'] = (time,"Creation date")
        hdu.header['UNIQUEID'] = (uid,"Unique identification number for this file.")
        hdu.header['RUID1'] = (self.uid1,"Built from this ramp uid.")
        hdu.header['RUID2'] = (self.uid2,"Built from this ramp uid.")
        hdu.writeto(join(self.pathtmp,'%sx%s.gain.fits'%(self.uid1,self.uid2)),overwrite=True) 
    def master_map(self):
        from scanf import scanf
        from os import listdir
        from os.path import basename
        ls = [join(self.pathtmp,f) for f in listdir(self.pathtmp) if ".gain.fits" in f]
        uid_pairs = [scanf("%sx%s.gain.fits",basename(f)) for f in ls]
        print("Will use uid pair:")
        for uid in uid_pairs:
            print("%s X %s"%(uid[0],uid[1]))
        cube = [fits.getdata(f)[0] for f in ls]
        
        gg= np.median(cube,axis=0)
        hdu = fits.PrimaryHDU(data=gg)
        from datetime import datetime
        time = datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
        uid = datetime.now().strftime("%y%m%d%H%M%S")
        hdu.header['DATE'] = (time,"Creation date")
        hdu.header['UNIQUEID'] = (uid,"Unique identification number for this file.")
        hdu.header.add_comment("Built using:")
        for uid in uid_pairs:
            hdu.header.add_comment("%s,%s"%(uid[0],uid[1]))
        hdu.writeto(join(self.path,"master_gain.fits"))
        print("master_gain.fits created")
        fig,ax = plt.subplots()
        im = ax.imshow(gg,cmap='YlGn',vmin=.7,vmax=1.9)
        plt.colorbar(im)
        
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        g,_,g_err = sc(gg)
        ax.set(title='Gain %2.2f+/-%2.2f adu/e$^-$'%(g,g_err))
        plt.tight_layout()
        fig.savefig(join(self.path,'gain.png'))
        sns.set_theme()
        f1,ax1 = plt.subplots()
        gg_flat = gg.ravel()
        gg_flat = gg_flat[(gg_flat>.9) & (gg_flat<1.9)]
        g,_,g_err = sc(gg_flat)
        ax1.hist(gg_flat,bins=50)
        ax1.set(title='Median, %.2f+/-%.2f'%(g,g_err),xlabel='Gain',ylabel='Occurance')
        plt.tight_layout()
        f1.savefig(join(self.path,'gain.dist.png'))
        plt.show()        
if '__main__' in __name__:
    path = '/nirps_raw/characterization'
    
    lsuid = [ '20210403171440',
              '20210403175352',
              '20210403183253',
              '20210403191153',
              '20210403195054',
              '20210403203006',
              '20210403210907',
              '20210403214818',
              '20210403222719',
              '20210403230620']#list of flat (200 reads separated by darks) (unique ID)
    #ls = [join(path,uid+'.fits') for uid in lsuid] #make a list with all the corrected files
    
    ll = [('20210403171440','20210403175352'),
          ('20210403183253','20210403191153'),
          ('20210403195054','20210403203006'),
          ('20210403210907','20210403214818'),
          ('20210403222719','20210403230620')]
    
    #from multiprocessing import Pool
    for i in range(len(ll)):
        g = gain()
        l = ll[i]
        g.makeXY(l[0],l[1], remakefile=True)
        g.run()
        # for amp in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32]:
        #     print("working on ramp # %d"%amp)
        #     g.p_run(amp)

        g.recombine()
    g = gain()
    g.master_map()
    
    
    
    
    '''
    lo = [];
    hi = [];
    for l in ll:
        g = gain()
        g.size = (32,32)
        g.makeXY(l[0],l[1])
        s,var = g.makesvar(2000,2000)
        lo.append(g.fit(s,var,sat=8000))
        hi.append(g.fit(s,var,sat=25000))

    _gain,gain_err,ro,ro_err = zip(*lo)
    print('Stats (lo): %.2f+/-%.2f adu/e-, ro: %.2f+/-%.2f (ADU)'%(sc(_gain)[0],sc(_gain)[2],sc(ro)[0],sc(ro)[2])) 
    _gain,gain_err,ro,ro_err = zip(*hi)
    print('Stats (hi): %.2f+/-%.2f adu/e-, ro: %.2f+/-%.2f (ADU)'%(sc(_gain)[0],sc(_gain)[2],sc(ro)[0],sc(ro)[2])) 
    g = gain()
    g.size = (256,256)
    g.makeXY('20210403195054','20210403203006')
    f,ax = plt.subplots()
    ax.plot(s,var,'s',markersize=3)
    gain,gain_err,ro,ro_err = g.fit(s,var,sat=8000)
    ax.plot(s,1.0/gain*s+ro**2,label='gain$_{lo}$: %.2f+/-%.2f adu/e$^-$, ro: %.2f+/-%.2f (ADU)'%(gain,gain_err,ro,ro_err))
    gain,gain_err,ro,ro_err = g.fit(s,var,sat=25000)
    ax.plot(s,1.0/gain*s+ro**2,label='gain$_{hi}$: %.2f+/-%.2f adu/e$^-$, ro: %.2f+/-%.2f (ADU)'%(gain,gain_err,ro,ro_err))
    
    ax.set(xlabel='Signal (ADU)',ylabel='Variance $\sigma^2$',title="Pixels: %dx%d"%(2000,2000))
    ax.legend()
    plt.tight_layout()
    #f.savefig('/home/jonathan/gain_%dx%128.png'%(2000,2000))
    plt.show()    
    
    
    
    g.size = (64,64)
    s,var = g.makesvar(2000,2000)
            
    f,ax = plt.subplots()
    ax.plot(s,var,'s',markersize=3)
    gain,gain_err,ro,ro_err = g.fit(s,var,sat=8000)
    ax.plot(s,1.0/gain*s+ro**2,label='gain$_{lo}$: %.2f+/-%.2f adu/e$^-$, ro: %.2f+/-%.2f (ADU)'%(gain,gain_err,ro,ro_err))
    gain,gain_err,ro,ro_err = g.fit(s,var,sat=25000)
    ax.plot(s,1.0/gain*s+ro**2,label='gain$_{hi}$: %.2f+/-%.2f adu/e$^-$, ro: %.2f+/-%.2f (ADU)'%(gain,gain_err,ro,ro_err))
    
    ax.set(xlabel='Signal (ADU)',ylabel='Variance $\sigma^2$',title="Pixels: %dx%d"%(2000,2000))
    ax.legend()
    plt.tight_layout()
    #f.savefig('/home/jonathan/gain_%dx%128.png'%(2000,2000))
    plt.show()
    '''
    