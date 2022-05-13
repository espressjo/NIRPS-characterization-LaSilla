#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 08:31:18 2021

@author: noboru
"""
import numpy as np
from tqdm import tqdm
from astropy.stats import sigma_clipped_stats as sc
from datetime import datetime
from os.path import join,basename,isdir
from astropy.io import fits
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from os import getcwd

def plot():
    import seaborn as sns
    sns.set_theme()
    f = plt.figure(figsize=(12,6))#,facecolor='beige')
    ax1 = plt.subplot2grid((3, 3), (0, 0), rowspan=2,colspan=3)
    ax2 = plt.subplot2grid((3, 3), (2, 0), rowspan=1,colspan=3)
    #ax1.grid()
    #ax2.grid()
    #ax1.set_facecolor('lightcyan')
    #ax2.set_facecolor('lightcyan')
    ax = [ax1,ax2]
    return f,ax
def flin(t,a,b):
    return a*t+b
def f3rd(t,a3,a2,a1,a0):
    return a3*t**3+a2*t**2+a1*t+a0
class nl():
    def __init__(self,f=''):
        '''
            ls: list of reads
        '''
        self.tmpdir = '/var/tmp'
        self.sample_time = 5.5733#5.24288
        if f!='':
            print("Uploading %s."%(f))
            self.cube = fits.getdata(f)
            self.nb_reads = self.cube.shape[0]
            self.time = (np.asarray(range(self.nb_reads))+1)*self.sample_time
            self.uid = basename(f).replace('.fits','')
        self.order=3
        self.a3 = np.zeros((4096,4096))
        self.a2 = np.zeros((4096,4096))
        self.res = np.zeros((4096,4096))
        self.chi2 = np.zeros((4096,4096))
        
        self.ite = 4#iteration to solve for a3 and a2
        self.saturation = 50000# saturation up to which to solve for NL
        
        self.lin_estimate = [100,0]
        self.poly_estimate = [-7e-13,2e-08,1e-04,-2]
        self.fname = ''
        self.rangex = [4,4092]
        self.rangey = [4,4092]
    def set_range_x(self,start,stop):
        self.rangex = [start,stop] 
    def set_range_y(self,start,stop):
        self.rangey = [start,stop]
    def set_name(self,name):
        self.fname = name
    def fit(self,time,f_meas):
        a2 = 0 # 2nd order coefficient
        a3 = 0
        time = time[f_meas<self.saturation]
        f = f_meas[f_meas<self.saturation]        
        for ite in range(self.ite):
            y2 = f + a2*f**2 + a3*f**3
            correction = np.polyfit(y2, y2-np.polyval(np.polyfit(time, y2, 1), time)  , 3)
            a2 -= correction[1]
            a3 -= correction[0]
        return a3,a2

    def unc(self,f):
        return np.sqrt(17**2+np.where(f<0,0,f))
    
    def get_chi2(self,x,y,beta):
        '''
        returned residuals of linear fit and reduced chi squared

        Parameters
        ----------
        x : TYPE
            DESCRIPTION.
        y : TYPE
            DESCRIPTION.
        y_ideal : TYPE
            DESCRIPTION.
        y_ideal_unc : TYPE
            DESCRIPTION.
        beta : TYPE
            DESCRIPTION.

        Returns
        -------
        residuals_std : TYPE
            DESCRIPTION.
        reduced_chi2 : TYPE
            DESCRIPTION.

        '''
        #nonlinearity = (y_ideal/y)
        #nonlinearity_unc = np.sqrt(nonlinearity**2*(self.unc(y)**2/y**2+y_ideal_unc**2/y_ideal**2))
       
        residuals = (y-flin(x,*beta))/y
        residuals_std = sc(residuals)[2]
        chi2 = np.sum(((flin(x,*beta)-y)/self.unc(y))**2)
        reduced_chi2 = chi2/(len(x)-len(beta))
        return residuals_std,reduced_chi2
    def fit_cf(self,x,f_meas):
        '''
        fit pixel using curve_fit and error estimate

        Parameters
        ----------
        time : TYPE
            DESCRIPTION.
        f_meas : TYPE
            DESCRIPTION.

        Returns
        -------
        a3,a2,res,reduced chi^2

        '''
        a2 = 0 # 2nd order coefficient, curvefit
        a3 = 0 #curvefit
        time = x[f_meas<self.saturation]
        f = f_meas[f_meas<self.saturation]
        for ite in range(self.ite):
            y2 = f + a2*f**2 + a3*f**3
            beta,lin_cov = curve_fit(flin,time,y2,sigma=self.unc(y2),p0=self.lin_estimate,maxfev=10000)
            signal_fit_cov = lin_cov[1,0] 
            err = np.sqrt(np.diag(lin_cov))
            ideal_fit_unc = np.sqrt(err[0]**2+(err[1]*time)**2+ 2*signal_fit_cov**2)
            err_2fit = ideal_fit_unc+self.unc(y2)            
            correction,c_cov = curve_fit(f3rd,y2, y2-flin(time,*beta),  p0=self.poly_estimate,sigma=err_2fit ,maxfev=10000)

            a2 -= correction[1]
            a3 -= correction[0]       
        beta,lin_cov = curve_fit(flin,time,y2,sigma=self.unc(y2),p0=self.lin_estimate,maxfev=10000)
        chi2_res,red_chi2 = self.get_chi2(time,y2,beta)
        return a3,a2,chi2_res,red_chi2
    def fit_graph(self,f_meas):
        '''
        Show an example of fit for one pixel using curve_fit and polyfit

        Parameters
        ----------
        f_meas : TYPE
            DESCRIPTION.

        Returns
        -------
        TYPE
            DESCRIPTION.

        '''
        a2 = 0 # 2nd order coefficient, curvefit
        a3 = 0 #curvefit
        _a2 = 0 # polyfit
        _a3 = 0 #polyfit
        x = 5.5733*np.arange(1,len(f_meas)+1,1)
        time = x[f_meas<self.saturation]
        f = f_meas[f_meas<self.saturation]    
        fig2,ax2 = plot()
        ax2[0].plot(time,f,'d',markersize=3,label='$f_{original}$')
        fig,ax = plot()
        #def ffit(t,a,b):
        #    return a*t+b
        def fffit(t,a3,a2,a1,a0):
            return a3*t**3+a2*t**2+a1*t+a0
        #y2 -> use curvefit
        #_y2 -> use polyfit
        for ite in range(self.ite):
            y2 = f + a2*f**2 + a3*f**3
            _y2 = f + _a2*f**2 + _a3*f**3
            #fit to find the S_ideal 
            _beta,_lin_cov = np.polyfit(time, _y2, 1,cov=True)
            beta,lin_cov = curve_fit(flin,time,y2,sigma=self.unc(y2),p0=[100,0],maxfev=10000)
            
            signal_fit_cov = lin_cov[1,0] 
            #_signal_fit_cov = _lin_cov[1,0]
            #_err = np.sqrt(np.diag(_lin_cov))
            err = np.sqrt(np.diag(lin_cov))
            ideal_fit_unc = np.sqrt(err[0]**2+(err[1]*time)**2+ 2*signal_fit_cov**2)
            #_ideal_fit_unc = np.sqrt(_err[0]**2+(_err[1]*time)**2+ 2*_signal_fit_cov**2)
            err_2fit = ideal_fit_unc+self.unc(y2)
            #_err_2fit = _ideal_fit_unc+self.unc(_y2)
            
            _correction,_c_cov = np.polyfit(_y2, _y2-flin(time,*_beta)  , 3,cov=True)
            correction,c_cov = curve_fit(fffit,y2, y2-flin(time,*beta),  p0=[-7e-13,2e-08,1e-04,-2],sigma=err_2fit ,maxfev=10000)

            a2 -= correction[1]
            a3 -= correction[0]
            _a2 -= _correction[1]
            _a3 -= _correction[0]

            ax[1].plot(time,y2-np.polyval(np.polyfit(time,y2,1),time),'.',label = 'iteration {0}'.format(ite))
            ax[0].plot(time,y2,'.',label = 'iteration {0}'.format(ite),alpha = 0.5)
        beta,lin_cov = curve_fit(flin,time,y2,sigma=self.unc(y2),p0=[100,0],maxfev=10000)
        chi2_res,red_chi2 = self.get_chi2(time,y2,beta)
        print(chi2_res,red_chi2)
        res = 1-(flin(time,*beta)/y2)
        _res = 1-(flin(time,*_beta)/_y2)
        #graph comparatif curve_fit vs polyfit        
        props = dict(boxstyle='round', facecolor='lightgreen', alpha=0.5)
        beta_,lin_cov_ = curve_fit(flin,time,y2,sigma=self.unc(y2),p0=[100,0],maxfev=10000)
        _beta_,_lin_cov_ = curve_fit(flin,time,_y2,sigma=self.unc(_y2),p0=[100,0],maxfev=10000)
        _err_ = np.sqrt(np.diag(_lin_cov_))
        err_ = np.sqrt(np.diag(lin_cov_))
        
        ax2[0].text(.8, 0.3, "Curvefit\nSlope: %.2f+/-%.2f\nb: %.2f+/-%.2f\nPolyfit\nSlope: %.2f+/-%.2f\nb: %.2f+/-%.2f"%(beta_[0],err_[0],beta_[1],err_[1],_beta_[0],_err_[0],_beta_[1],_err_[1] ),
                    transform=ax2[0].transAxes, fontsize=8,verticalalignment='top', bbox=props)
        ax2[0].plot(time,_y2,'d',markersize=3,label='$f_{polyfit}$')
        ax2[0].plot(time,y2,'d',markersize=3,label='$f_{curvefit}$')
        ax2[0].plot(time,flin(time,*beta),label='$f_{ideal,curvefit}$')
        ax2[0].plot(time,flin(time,*_beta),label='$f_{ideal,polyfit}$')
        ax2[1].plot(time,res,'d',markersize=3,label='curvefit')  
        ax2[1].plot(time,_res,'d',markersize=3,label='polyfit') 
        
        ax2[1].text(.8, 0.2, "res. (curvefit): %f\nres. (polyfit): %f"%(sc(res)[2],sc(_res)[2]),
                    transform=ax2[1].transAxes, fontsize=8,verticalalignment='top', bbox=props)
        ax2[0].legend()
        ax2[1].legend()
        ax2[1].set(xlabel = 'Time (second)', ylabel = 'Residuals')
        ax2[0].set(xlabel = 'Time (second)', ylabel = 'flux (DN)')
        
        print('res: %f'%(sc(res)[2]))
        #self.get_chi2(time,)
        ax[1].legend()
        ax[0].legend()
        ax[1].set(xlabel = 'time', ylabel = 'flux - lin model')
        ax[0].set(xlabel = 'time', ylabel = 'flux')
        
        return a3,a2
    def fit_graph_rapport(self,f_meas,xpos=2000,ypos=2000):
        '''
        Show an example of fit for one pixel using curve_fit and polyfit

        Parameters
        ----------
        f_meas : TYPE
            DESCRIPTION.

        Returns
        -------
        TYPE
            DESCRIPTION.

        '''
        a2 = 0 # 2nd order coefficient, curvefit
        a3 = 0 #curvefit
        _a2 = 0 # polyfit
        _a3 = 0 #polyfit
        x = 5.5733*np.arange(1,len(f_meas)+1,1)
        time = x[f_meas<self.saturation]
        f = f_meas[f_meas<self.saturation]    
        fig2,ax2 = plot()
        ax2[0].plot(time,f,'d',markersize=3,label='$f_{original}$')
        fig,ax = plot()
        #def ffit(t,a,b):
        #    return a*t+b
        def fffit(t,a3,a2,a1,a0):
            return a3*t**3+a2*t**2+a1*t+a0
        #y2 -> use curvefit
        #_y2 -> use polyfit
        for ite in range(self.ite):
            y2 = f + a2*f**2 + a3*f**3
            _y2 = f + _a2*f**2 + _a3*f**3
            #fit to find the S_ideal 
            _beta,_lin_cov = np.polyfit(time, _y2, 1,cov=True)
            beta,lin_cov = curve_fit(flin,time,y2,sigma=self.unc(y2),p0=[100,0],maxfev=10000)
            
            signal_fit_cov = lin_cov[1,0] 
            #_signal_fit_cov = _lin_cov[1,0]
            #_err = np.sqrt(np.diag(_lin_cov))
            err = np.sqrt(np.diag(lin_cov))
            ideal_fit_unc = np.sqrt(err[0]**2+(err[1]*time)**2+ 2*signal_fit_cov**2)
            #_ideal_fit_unc = np.sqrt(_err[0]**2+(_err[1]*time)**2+ 2*_signal_fit_cov**2)
            err_2fit = ideal_fit_unc+self.unc(y2)
            #_err_2fit = _ideal_fit_unc+self.unc(_y2)
            
            _correction,_c_cov = np.polyfit(_y2, _y2-flin(time,*_beta)  , 3,cov=True)
            correction,c_cov = curve_fit(fffit,y2, y2-flin(time,*beta),  p0=[-7e-13,2e-08,1e-04,-2],sigma=err_2fit ,maxfev=10000)

            a2 -= correction[1]
            a3 -= correction[0]
            _a2 -= _correction[1]
            _a3 -= _correction[0]

            ax[1].plot(time,y2-np.polyval(np.polyfit(time,y2,1),time),'.',label = 'iteration {0}'.format(ite))
            ax[0].plot(time,y2,'.',label = 'iteration {0}'.format(ite),alpha = 0.5)
        beta,lin_cov = curve_fit(flin,time,y2,sigma=self.unc(y2),p0=[100,0],maxfev=10000)
        chi2_res,red_chi2 = self.get_chi2(time,y2,beta)
        print(chi2_res,red_chi2)
        res = 1-(flin(time,*beta)/y2)
        _res = 1-(flin(time,*_beta)/_y2)
        #graph comparatif curve_fit vs polyfit        
        props = dict(boxstyle='round', facecolor='lightgreen', alpha=0.5)
        beta_,lin_cov_ = curve_fit(flin,time,y2,sigma=self.unc(y2),p0=[100,0],maxfev=10000)
        _beta_,_lin_cov_ = curve_fit(flin,time,_y2,sigma=self.unc(_y2),p0=[100,0],maxfev=10000)
        _err_ = np.sqrt(np.diag(_lin_cov_))
        err_ = np.sqrt(np.diag(lin_cov_))
        
        ax2[0].text(.8, 0.3, "Slope: %.2f+/-%.2f\nb: %.2f+/-%.2f"%(beta_[0],err_[0],beta_[1],err_[1]),
                    transform=ax2[0].transAxes, fontsize=8,verticalalignment='top', bbox=props)
        ax2[0].plot(time,y2,'d',markersize=3,label='$f_{corr}$')
        ax2[0].plot(time,flin(time,*beta),label='$f_{ideal}$')
        ax2[1].plot(time,res,'d',markersize=3,label='Residuals')  
        ax2[1].set_ylim([-.01,.01])
        ax2[1].text(.8, 0.2, "res. : %f "%(sc(res)[2]),
                    transform=ax2[1].transAxes, fontsize=8,verticalalignment='top', bbox=props)
        ax2[0].set(title='Pixel X:%d Y:%d'%(xpos,ypos))
        ax2[0].legend()
        ax2[1].legend()
        ax2[1].set(xlabel = 'Time (second)', ylabel = 'Residuals')
        ax2[0].set(xlabel = 'Time (second)', ylabel = 'flux (DN)')
        
        print('res: %f'%(sc(res)[2]))
        #self.get_chi2(time,)
        ax[1].legend()
        ax[0].legend()
        ax[1].set(xlabel = 'time', ylabel = 'flux - lin model')
        ax[0].set(xlabel = 'time', ylabel = 'flux')
        plt.show()
        return a3,a2
    def fit_a2(self,time,f_meas):
        a2 = 0 # 2nd order coefficient
        time = time[f_meas<self.saturation]
        f = f_meas[f_meas<self.saturation]        
        for ite in range(self.ite):
            y2 = f + a2*f**2
            correction = np.polyfit(y2, y2-np.polyval(np.polyfit(time, y2, 1), time)  , 2)
            a2 -= correction[0]
            
        return 0,a2
    def make_coeff_map_cf(self,p):
        '''
        Create a NL map using curve_fit.

        Parameters
        ----------
        p : STR
            Path where to save the file.

        Returns
        -------
        None.

        '''
        print("Solving for a3 and a2")
       
        for x in tqdm(range(self.rangex[0],self.rangex[1])):#original 4,4092
            for y in range(self.rangey[0],self.rangey[1]):#original 4,4092
                try:
                    a3,a2,res,chi2 = self.fit_cf(self.time,self.cube[:,y,x])
                except:
                    a3 = 0
                    a2 = 0
                    res = 0
                    chi2 = 0
                self.a3[y,x] = a3
                self.a2[y,x] = a2
                self.res[y,x] = res
                self.chi2[y,x] = chi2
        dt = datetime.now()
        print('creating FITS file')
        creation_time = dt.strftime("%Y-%m-%dT%H:%M:%S")
        hdu = fits.PrimaryHDU()
        hdu.data = np.asarray([self.a3,self.a2,self.res,self.chi2])
        hdu.header['DATE'] = creation_time
        if self.fname=='':
            self.fname = "nl_cf_artigau_%s.fits"%(dt.strftime("%Y%m%d%H%M%S"))
            
        hdu.writeto(join(p,self.fname),overwrite=True)
        print('Done.')
    def make_coeff_map(self,order=3):
        print("Solving for a3 and a2")
        if order!=3 and order!=2:
            print("only order 2 and 3 are supported")
            return
        if order!=3:
            self.order=2
            fit = self.fit_a2
        else:
            fit = self.fit
        for x in tqdm(range(self.rangex[0],self.rangex[1])):#original 4,4092
            for y in range(self.rangey[0],self.rangey[1]):#original 4,4092
                try:
                    a3,a2 = fit(self.time,self.cube[:,y,x])
                except:
                    a3 = 0
                    a2 = 0
                self.a3[y,x] = a3
                self.a2[y,x] = a2
                
        print('Done.')
        
    def save_coeff(self,p):
        dt = datetime.now()
        creation_time = dt.strftime("%Y-%m-%dT%H:%M:%S")
        if self.order==3:
            primary_hdu = fits.PrimaryHDU(self.a3)
            image_hdu = fits.ImageHDU(self.a2)
        else:
            primary_hdu = fits.PrimaryHDU(self.a2)
        primary_hdu.header['DATE'] = creation_time
        #for n in self.ls:
        #    primary_hdu.header.add_comment(n)
        if self.order==3:    
            hdul = fits.HDUList([primary_hdu, image_hdu])
        else:
            hdul = fits.HDUList(primary_hdu)
        if self.fname=='':
            self.fname = "nl_%s.fits.subset"%(dt.strftime("%Y%m%d%H%M%S"))
            
        hdul.writeto(join(p,self.fname),overwrite=True)
        hdul.close()
    def subsection_fit(self):
        '''
        Fit the coefficient over a small subsection to check if everything is alright

        '''
        
        coeff = np.zeros((4,100,100))
        
        
        for y in range(2000,2100):
            for x in range(2000,2100):
                try:
                    a3,a2,res,chi2 = self.fit_cf(self.time,self.cube[:,y,x])
                except :
                    a3=0
                    a2=0 
                    res=0
                    chi2=0
                
                coeff[0,y-2000,x-2000] = a3
                coeff[1,y-2000,x-2000] = a2
                coeff[2,y-2000,x-2000] = res
                coeff[3,y-2000,x-2000] = chi2
               
        hdu = fits.PrimaryHDU()
        hdu.data = coeff
        hdu.writeto('fit.check.fits',overwrite=True)
    def p_run(self,amp):
        '''
        amp [1,32]

        Parameters
        ----------
        data : TYPE
            DESCRIPTION.
        amp : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        
        coeff = np.zeros((4,4096,128))
        #print("Amp: %d"%amp)
        for y in (range(4096)):
            for x in range(128*(amp-1),(amp)*128):
                try:
                    a3,a2,res,chi2 = self.fit_cf(self.time,self.cube[:,y,x])
                except :
                    a3=0
                    a2=0 
                    res=0
                    chi2=0
                
                coeff[0,y,x-((amp-1)*128)] = a3
                coeff[1,y,x-((amp-1)*128)] = a2
                coeff[2,y,x-((amp-1)*128)] = res
                coeff[3,y,x-((amp-1)*128)] = chi2
               
        hdu = fits.PrimaryHDU()
        hdu.data = coeff
        hdu.writeto(join(self.tmpdir,'%s.amp.%d'%(self.uid,amp)),overwrite=True)
    def run(self):
        from multiprocessing import Process
        
        ps = [ Process(target=self.p_run,args=(amp,)) for amp in range(1,33)]
        for i in range(32):
            print("Launching amp #%d"%(i+1))
            ps[i].start()
        for i in range(32):
            ps[i].join()
            
        
        
    def recombine(self,path=None):
        from os import listdir
        from os import getcwd
        from scanf import scanf
        
        ls = [join(self.tmpdir,f) for f in listdir(self.tmpdir) if any(['%s.amp'%self.uid in f,'%s.bias.amp'%self.uid in f])]

        coeff = np.zeros((4,4096,4096))
        for f in ls:
            amp = scanf("%s.amp.%d",basename(f))[1]
            print("Parsing amp #%.2d"%amp)
            coeff[:,:,(amp-1)*128:(amp)*128] = fits.getdata(f)

        hdu = fits.PrimaryHDU()
        hdu.data = coeff
        from datetime import datetime
        time = datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
        uid = datetime.now().strftime("%y%m%d%H%M%S")
        hdu.header['DATE'] = (time,"Creation date")
        hdu.header['UNIQUEID'] = (uid,"Unique identification number for this file.")
        hdu.header['RUID'] = (self.uid,"Built from this ramp uid.")
        hdu.writeto(join(path,'%s.nl.fits'%self.uid),overwrite=True) 
    def make_master_map(self,ls,fname,path=''):
        '''
        Will create a master nonlinearity map

        Parameters
        ----------
        ls : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        lsuid = [basename(f).replace('.nl.fits','') for f in ls]
        if path=='':
            path = getcwd()
        from astropy.stats import median_absolute_deviation as mad
        from datetime import datetime
        nb_ext = fits.getdata(ls[0]).shape[0]

        nl = np.zeros((3,4096,4096))

        res = np.asarray([fits.getdata(f)[2,:,:] for f in ls])
        chi2 = np.asarray([fits.getdata(f)[3,:,:] for f in ls])
        md = np.median(res,axis=0)
        msd = mad(res,axis=0)
        top = md+5*msd 
        mask_outliers = np.asarray([ np.where(im>top,True,False) for im in res],dtype=bool)
        mask_chi2 = np.asarray([ np.where(im>=1.1,True,False) for im in chi2],dtype=bool)
        bc2 = np.all(mask_chi2,axis=0)
        nl[2] = np.asarray(bc2,dtype=int)
        mask3d = mask_outliers | mask_chi2
        for i in range(2):
            ext = np.asarray([fits.getdata(f)[i] for f in ls])   
            a = np.ma.array(ext,mask=mask3d)
            mean_ext = a.mean(axis=0)
            nl[i,:,:] = mean_ext.data
        def make_file(nl,name,ls=[],path=''):
            hdu = fits.PrimaryHDU()
            #import pdb 
            #pdb.set_trace()
            hdu.data = nl
            time = datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S")
            uid = datetime.utcnow().strftime("%y%m%d%H%M%S")
            hdu.header['DATE'] = (time,"Creation time")
            hdu.header['UNIQUEID'] = (uid,"Unique identification number.")
            hdu.header.add_comment("Built using:")
            for u in lsuid:
                hdu.header.add_comment(u)
            if len(ls)>0:
                for f in ls:
                    hdu.header.add_comment(f)
            if path=='':
                path = getcwd()
            
            hdu.writeto(join(path,name),overwrite=True)
        make_file(nl,fname,path=path)
    # def make_flat(self,uid,cpath='/nirps_raw/characterization',sbias='superbias_20210416.fits'):
    #     from os.path import isfile
    #     if not isfile(join(cpath,'nonlin.fits')) :
    #         print("nonlin.fits not found")
    #         exit(0)
    #     if not isfile(join(cpath,sbias)) and sbias!='':
    #         print("%s not found"%sbias)
    #         exit(0)
    #     c3,c2,_ = fits.getdata(join(cpath,'nonlin.fits'))
    #     if sbias!='':
    #         bias = fits.getdata(join(cpath,sbias))
        
    def make_report(self,nonlinf='nonlin.fits',cpath='/nirps_raw/characterization',sbias='superbias_20210416.fits'):
        nl = fits.getdata(join(cpath,nonlinf))
        c3 = nl[0]
        c2 = nl[1]
        bpm = nl[2]
        if sbias=='':
            bias = 0
        else:
            bias = sc(fits.getdata(join(cpath,sbias)),sigma=5)[0]
        
        statsC3 = sc(c3[bpm==0],sigma=5)
        statsC2 = sc(c2[bpm==0],sigma=5)
        statsbpm = len(bpm[bpm==1].ravel())
        print("C3: %.2E +/- %.2E"%(statsC3[0],statsC3[2]))
        print("C2: %.2E +/- %.2E"%(statsC2[0],statsC2[2]))
        print("%f percent of pixels are BAD."%( (float(statsbpm)/4096**2)*100 ) )
        from matplotlib import pyplot as plt
        import seaborn as sns
        sns.set_theme()
        f,ax = plt.subplots(1,2)
        distc3 = c3[(c3>np.percentile(c3,0.15)) & (c3<np.percentile(c3,99.5))]
        distc2 = c2[(c2>np.percentile(c2,0.3)) & (c2<np.percentile(c2,99.9))]
        ax[0].set_yscale('log')
        ax[1].set_yscale('log')
        ax[0].hist(distc3,bins=200)
        ax[0].set(xlabel='Non-linearity value',ylabel='Count',title='C$_3$')
        ax[1].hist(distc2,bins=200)
        ax[1].set(xlabel='Non-linearity value',ylabel='Count',title='C$_2$')
        f.suptitle('Non-linearity Coefficient distribution')
        
        plt.tight_layout()
        f.savefig(join(cpath,'nl-hist.png'))
        print("%s saved successfully."%join(cpath,'nl-hist.png'))
        #check the 3% and 5 % NL
        flux = np.arange(1,65000,1)
        perc = (((flux+statsC3[0]*flux**3+statsC2[0]*flux**2)/flux)-1)*100
        def closest(arr,val):
            return np.argmin(abs(arr-val))
        perc3 = flux[closest(perc,3)]
        perc5 = flux[closest(perc,5)]
        print("3 percent non-linear @ %f, 5 percent non-linear @ %f"%(perc3,perc5))
        
        plt.show()
if '__main__' in __name__:
    
        
    
    #from multiprocessing import Pool
    
    path = '/nirps_raw/characterization'
    path_tmp = '/nirps_raw/characterization/.tmp'
    '''
    if not isdir(path_tmp):
        from os import mkdir
        mkdir(path_tmp)
    bias_subtracted = True
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
    if bias_subtracted:
        ext = ".bias.fits" #
    else:
        ext = ".fits"
    ls = [join(path,uid+ext) for uid in lsuid] #make a list with all the corrected files
    
    for fname in ls:
        print("Working on %s"%(basename(fname)))
        nla = nl(fname)
        nla.tmpdir = path_tmp

        nla.run()    
        nla.recombine(path=path_tmp)
    
    lsnl = [join(path_tmp,uid+".bias.nl.fits") for uid in lsuid]
    nla.make_master_map(lsnl,'nonlin.fits',path=path)
    '''
    n = nl()
    #with bias subtracted
    f_meas = fits.getdata('/nirps_raw/characterization/20210403171440.bias.fits')[:,2000,2000]
    #n.fit_graph_rapport(f_meas,xpos=2000,ypos=2000)
    
    
    n.make_report()
    
    
    