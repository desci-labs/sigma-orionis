#!/usr/bin/env

"""

 run visbin and plot
 execfile('dovisbin.py')

"""


# ===================== Import Packages ====================

import os,pdb,glob,sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import units as u
    
# ===================== Define Functions ===================



def get_ctr(fp,dx,dy):

    ### get source center coords (deg)
    hdr = fits.getheader(fp+'.fits')
    ra_obs = hdr['OBSRA']
    de_obs = hdr['OBSDEC']
    ra_deg = ra_obs+dx/3600.
    de_deg = de_obs+dy/3600.

    ### change dec to sexagesimal
    ideg = np.fix(de_deg)
    imin = np.fix(np.abs(de_deg-ideg)*60.)
    xsec = ((np.abs(de_deg-ideg)*60.) - imin)*60.
    de_ctr = str(int(ideg))+'d'+str(int(imin))+'m'+str(xsec)[0:7]

    ### change ra to sexagesimal
    ihr = np.fix(ra_deg/15.)
    xmin = np.abs(ra_deg/15.*60. - ihr*60.)
    imin = np.fix(xmin)
    xsec = (xmin-imin)*60.0
    ra_ctr = str(int(ihr))+'h'+str(int(imin))+'m'+str(xsec)[0:7]

    ### join in correct format
    c_ptc = 'ICRS '+ ra_ctr + ' ' + de_ctr

    return c_ptc

        

def run_visbin (fp,i,PA,c_ptc):

    visbin(vis = fp+'.vis',
           incl = i,
           posang = PA,
           binwidth = 8e4,
           nbins = 30,
           bmin = -1,
           binfreq = False,
           fit_phasecenter = False,
           phasecenter = c_ptc,
           outfile = fp.split('/')[0],
           overwrite = True)

    os.system('rm -r '+fp.split('/')[0]+'/'+fp.split('/')[0]+'.visbin')
    os.system('mv '+fp.split('/')[0]+'.visbin '+fp.split('/')[0]+'/')

    
    
def plot_fit(fp):

    uvfit = np.genfromtxt(fp.split('/')[0]+'/'+fp.split('/')[0]+'.visbin/'+fp.split('/')[0]+'.visbin.binned',
                          dtype={'names': ('dist','Re','e_Re','Im','e_Im'),'formats': (float,float,float,float,float)})
    
    Dist,Re,e_Re,Im,e_Im = uvfit['dist'][1:]/1e3,uvfit['Re'][1:],uvfit['e_Re'][1:],uvfit['Im'][1:],uvfit['e_Im'][1:] 

    mpl.rc('xtick',labelsize=10) 
    mpl.rc('ytick',labelsize=10)
    mpl.rc('xtick.major',size=5,pad=7,width=.5)
    mpl.rc('ytick.major',size=5,pad=7,width=.5)
    mpl.rc('axes',linewidth = .5)
    mpl.rc('lines',markersize=5)

    fig = plt.figure(figsize = (8,6))
    ax = fig.add_subplot(111)
    ax.set_xlabel(r'UV Distance'+' ['+r'k$\lambda$'+']',fontsize=10)
    ax.set_ylabel(r'Re(v) [Jy]',fontsize=10)

    ax.plot(Dist,Re,'-',lw=1,color='gray')
    ax.errorbar(Dist,Re,yerr=e_Re,fmt='o',ms=8,mew=1,mec='gray',mfc='dodgerblue',ecolor='gray')
    
    ax.plot(Dist,Im,'-',lw=1,color='gray')
    ax.errorbar(Dist,Im,yerr=e_Im,fmt='o',ms=4,mew=0.5,mec='gray',mfc='lightgray',ecolor='gray')

    ax.axhline(y=0.0,color='darkgray',linestyle=':')
    # ax.text(np.max(Dist)-100.,np.max(Re),'i = '+str(i[0])+r'$\pm$'+str(e_i[0])+' deg',fontsize=8,color='Black')
    # ax.text(np.max(Dist)-100.,np.max(Re)-0.03,'PA = '+str(PA[0])+r'$\pm$'+str(e_PA[0])+' deg',fontsize=8,color='Black')

    fig.savefig(fp.split('/')[0]+'/'+fp.split('/')[0]+'.visbin/plot_visbin.eps', bbox_inches='tight',format='eps', dpi=100)
    plt.close('all')


    
# ========================== Code ==========================

fn = '86'
R  =  0.62
PA = 14
dx = -0.2223
dy = -1.5663

fp = 'S'+str(fn)+'/S'+str(int(fn))+'_cont'
c_ptc = get_ctr(fp,dx,dy)
i = np.arccos(R)*(180./np.pi)

run_visbin(fp,i,PA,c_ptc) 
plot_fit(fp)


