#!/usr/bin/env

'''

Get spectrum from ALMA cube
Assumes no emission in first channel

'''


# ======================== Import Packages ==========================

import sys,os,pdb,glob
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.table import Table,vstack
from astropy import constants as const
from astropy import units as u
from scipy.optimize import curve_fit
import changecoords
import get_lsrk


# ========================== Define Fuctions ==========================



def get_cube(fn):

    '''
    GET IMAGE CUBE AND HEADER
    
    '''
    
    cube, hdr = fits.getdata(fn,header=True)
    cube = np.squeeze(cube)
    cube[np.isnan(cube)] = 0

    # plt.imshow(cube[3,:,:],origin='lower')
    # plt.show()
    # pdb.set_trace()
        
    return cube,hdr



def get_grid(cube,hdr,xoff,yoff):

    '''
    MAKE RADIUS GRID (IN ARCSEC) CENTERED ON SOURCE
    
    '''
    
    nx,ny = np.shape(cube)[1],np.shape(cube)[2]              # size of input cube
    dx,dy = hdr['CDELT1']*3600.,hdr['CDELT2']*3600.          # pixel scale (as/pix; -ve ra!)
    cx,cy = nx/2,ny/2                                        # cube center (python pix)
    
    xr = (np.arange(-1*dx*(cx),dx*(nx-cx),dx)-xoff)          # x-dist from source center (as)
    yr = (np.arange(-1*dy*(cy),dy*(ny-cy),dy)-yoff).T        # y-dist from source center (as)
    xx,yy  = np.meshgrid(xr,yr)   
    rr = np.sqrt(xx**2+yy**2)

    # ctr = np.where(rr == np.min(rr))
    # plt.imshow(np.sum(cube,0),origin='lower')
    # plt.axvline(x=float(ctr[1]),color='black',linestyle=':')
    # plt.axhline(y=float(ctr[0]),color='black',linestyle=':')
    # plt.show()
    # pdb.set_trace()

    return rr



def get_rms(cube,rr,sr):

    '''
    GET RMS WITHIN APERTURE
    
    '''

    ind = np.where( rr <= sr )
    rms = []
    for i in np.arange(np.shape(cube)[0]):
        rms = np.append(rms,np.nanstd(cube[i,:,:][ind]))
    rms = np.mean(rms[1:])

    return rms



def get_beamsize(hdr,sr):

    '''
    CALCULATE NUMBER OF PIXELS PER BEAM
    AND NUMBER OF BEAMS IN SOURCE APERTURE 
    USING HEADER INFORMATION
    
    '''

    ### get info from header
    ps = abs(hdr['CDELT1'])
    bmaj = hdr['BMAJ']
    bmin = hdr['BMIN']

    ### calculate values
    pix_per_beam = np.pi/4.0/np.log(2.0) * (bmaj/ps) * (bmin/ps)
    pix_per_aper = np.pi*(sr/(abs(hdr['CDELT1'])*3600.))**2
    nbeam = pix_per_aper/pix_per_beam

    return pix_per_beam,nbeam



def get_flx(cube,rr,sr,bp):

    '''
    GET FLUX DENSITY WITHIN APERTURE
    (MATCHES FLUX DENSITY IN CASA STATISTICS)
    
    '''
    
    nvel = cube.shape[0]
    flx = np.zeros(nvel)
    
    for k in range(0,nvel):
        
        slc = cube[k,:,:]
        # w = np.logical_and( rr < sr , slc > -nlim*rms )
        ind = np.where(rr < sr)
        flx[k]=np.nansum(slc[ind])/bp
               
    return flx


    
def get_vel(hdr):

    '''
    GET VELOCITY CHANNELS (KM/S)
    USING HEADER FREQUENCY INFORMATION
    
    '''

    ### GET HEADER INFO
    rf_hz = hdr['RESTFRQ']         # REST FREQ
    fo_hz = hdr['CRVAL3']          # REF FREQ
    df_hz = hdr['CDELT3']          # FREQ SPACING
    cf_px = hdr['CRPIX3'] - 1      # REF PIXEL
    nf    = hdr['NAXIS3']          # DIMENSIONS

    ### CALC CHANNEL FREQ
    freq  = (np.arange(nf)-cf_px)*df_hz+fo_hz

    ### CALC CHANNEL VELOCITIES    
    vel   = np.zeros(nf)
    for q in range(0,nf):
		vel[q] = float( "{0:0.3f}".format((-1*const.c.value/rf_hz*(freq[q]-rf_hz))/1e3) )

    return vel



def gauss(x,*p):

    A,mu,sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))



def gauss_dofit(xval,yval,xfit,*p_guess):

    popt, pcov = curve_fit(gauss,xval,yval,p0=p_guess)
    yfit = gauss(xfit,*popt) 

    return popt,pcov,yfit



def plot_spec(name,flx,vel,rms,nb,sr,rv):

    
    flx = 1e3*flx[1:]
    vel = vel[1:]
    rms = 1e3*rms

    outfile = 'CO_spec_'+str(sr)+'/'+data['name'][i]+'_'+linename[n]
    # outfile = data['name'][i]+'_'+linename[n]
    np.savetxt(outfile+'.txt',np.c_[vel,flx],fmt='%10.2f %10.4f')

    fig = plt.figure(figsize = (8,6))
    ax = fig.add_subplot(111)
    ax.set_xlabel("Velocity (km/s)",fontsize=15) 
    ax.set_ylabel("Flux (mJy/beam)",fontsize=15)
    ax.tick_params(axis='both',which='major',labelsize=12)
    
    ax.plot(vel,flx,drawstyle='steps-mid',color='blue',lw=2)
    ax.axhline(y=0,color='gray',linestyle=':')
    ax.axhline(y=rms*3.0,color='orange',lw=3,alpha=0.3)
    ax.axhline(y=rms*4.0,color='red',lw=3,alpha=0.3)
    
    if rv != -99.0:
        ax.axvline(x=rv,color='gray',lw=3,alpha=0.3)
    else:
        ax.axvline(x=13,color='gray',lw=1,linestyle=':')

    fig.savefig(outfile+'.pdf',bbox_inches='tight',dpi=100,alpha=True,rasterized=True)
    plt.close('all')
    

    # ind = np.where(flx > rms_test*nl)
    # if len(ind[0]) > 0:
    #     print "v-range: "+str(vel[ind][0])+'-'+str(vel[ind][-1])+'km/s'
    # if len(ind[0]) == 0:
    #     ind = np.where( (vel>=2.0) & (vel<=6.0))
    #     print 'noisy spectrum, using default 2-6 km/s v-range'
    # ax.axvline(x=vel[ind][0],color='red',linestyle=':')
    # ax.axvline(x=vel[ind][-1],color='red',linestyle=':')
    # p_guess = [35.0, 4.5, 1.0]
    # xfit = np.arange(vel[0],vel[-1],0.1)
    # popt,pcov,yfit = gauss_dofit(vel,flx,xfit,*p_guess)
    # ax.plot(xfit,yfit,color='green')
    # ax.text(-1,popt[0]*0.9,"{0:.1f}".format(popt[1])+' km/s',color='green',size=10)
    # print "v-center: "+"{0:.1f}".format(popt[1])+'km/s'
            
    
def rv_convert(rv,ra,de):

    if rv != -99.0:
        corr = get_lsrk.get_corr(ra,de)
        rv_lsrk = rv + corr
    else:
        print "unknown RV"
        pdb.set_trace()
    
    return rv_lsrk


 
# ============================= Code ==================================


linename = ['12CO','13CO','C18O']
data = Table.read('/Users/mansdell/Dropbox/thesis/sori/catalog/sori_cat.fits',format='fits')
data.sort('name')

data = data[np.where(data['name'] == 'S55')]
skip = np.array(['S29','S56','S74','S75','S76','S88'])
skip = np.array([])


for i,val in enumerate(data['name']):

    if val in skip:
        print "\n>>> "+data['name'][i]+" is known gas detection, skipping..."
        continue
    else:        
        print "\n>>> Getting spectrum of "+data['name'][i]+'\n'

    field = int(data['name'][i][1:])
    dx,dy = data[i]['dx'],data[i]['dy']
    
    sr = 0.3
    nl = 3.0
    
    for n,nval in enumerate(linename):

        fn = data['name'][i]+'/S'+str(field)+'_'+linename[n]+'.fits'
        cube,hdr = get_cube(fn)

        rr = get_grid(cube,hdr,dx,dy)
        rms = get_rms(cube,rr,sr)
        bp,bn = get_beamsize(hdr,sr)

        print "    " + nval+'\n'
        print "      [" + str(dx) +','+ str(dy) + "]"
        print "      {0:0.1f} mJy/beam RMS".format(rms*1e3)
        print "      {0:0.1f} pixels per beam".format(bp)
        print "      {0:0.1f} beams in aperture\n".format(bn)

        flx = get_flx(cube,rr,sr,bp)
        vel = get_vel(hdr)

        if data['RV'][i] != -99.0:
            rv_lsrk = rv_convert(data['RV'][i],data['ra'][i],data['de'][i])
        else:
            rv_lsrk = -99.0
            
        plot_spec(data['name'][i],flx,vel,rms,bn,sr,rv_lsrk)

    
