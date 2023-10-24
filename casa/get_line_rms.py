#!/usr/bin/env

import numpy as np
import sys,os,pdb,glob
from astropy.table import Table
from astropy.io import fits
from astropy import constants as const



def get_cube(filename):

    '''
    GET IMAGE CUBE AND HEADER
    CONVERT IMAGE TO MJY
    
    '''
    
    cube, hdr = fits.getdata(filename,header=True)
    cube = np.squeeze(cube)*1e3
    cube[np.isnan(cube)] = 0

    return cube,hdr



def get_rms(cube,hdr,dxdy=[0,0],annulus=[4,9]):

    '''
    GET RMS WITHIN ANNULUS
    ASSUMES FIRST CHANNEL IS NOISY
    ASSUMES NO SOURCE EMISSION IN SECOND CHANNEL
    
    '''

    ### MAKE RADIAL GRID CENTERED ON SOURCE
    nx,ny = np.shape(cube)[1],np.shape(cube)[2]                  # SIZE OF INPUT CUBE
    dx,dy = hdr['CDELT1']*3600.,hdr['CDELT2']*3600.              # PIXEL SCALE (as/pix)
    cx,cy = nx/2,ny/2                                            # CUBE CENTER (pix)   
    xr = (np.arange(-1*dx*(cx),dx*(nx-cx),dx)-dxdy[0])           # X-DIST FROM SOURCE (as)
    yr = (np.arange(-1*dy*(cy),dy*(ny-cy),dy)-dxdy[1]).T         # Y-DIST FROM SOURCE (as)
    xx,yy = np.meshgrid(xr,yr)
    rr = np.sqrt(xx**2+yy**2)

    ### MEASURE RMS IN EACH CHANNEL
    rms = np.zeros(cube.shape[0])
    for i,val in enumerate(np.arange(0,cube.shape[0])):
        chan = cube[val,:,:]
        ind = np.where( (rr <= annulus[1]) & (rr >= annulus[0]) )    
        rms[i] = np.nanstd(chan[ind])

    ### IGNORE FIRST CHANNEL
    rms_avg = round(np.mean(rms[1:]),4)
    rms_std  = round(np.std(rms[1:]),4)
    
    return rms_avg,rms_std



def get_vel(hdr):

    '''
    GET VELOCITY CHANNELS OF CUBE (KM/S)
    FROM HEADER FREQUENCY INFORMATION
    USING "RADIO DEFINITION" OF CONVERTING dF TO dV
    GOES TO PRECISION OF 0.001KM/S!

    INPUTS:
       HDR = FITS HEADER OF CUBE

    OUTPUTS:
       VEL  = VELOCITIES OF EACH CHANNEL
       FREQ = FREQUENCIES OF EACH CHANNEL
    
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
        vel[q] = round((-1*const.c.value/rf_hz*(freq[q]-rf_hz))/1e3,4)

    return vel,freq



line = 'C18O'

data = Table.read('/Users/meganansdell/Dropbox/thesis/sori/alma/cont_results.txt',format='ascii.ipac')
dirs = np.array(data['name'])
rms_all = np.zeros(len(dirs))
for i,val in enumerate(dirs):

    field = int(dirs[i][1:])
    cubefile = dirs[i]+'/S'+str(field)+'_'+line+'.fits'

    cube,hdr = get_cube(cubefile)

    vel,freq = get_vel(hdr)
    rms,rms_std = get_rms(cube,hdr,dxdy=[float(data['dx'][i]),float(data['dy'][i])],annulus=[4,9])

    rms_all[i] = rms
    print(val+' '+str(rms)+' '+str(rms_std))


