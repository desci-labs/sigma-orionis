#!/usr/bin/env

import numpy as np
import sys,os,pdb,glob
from astropy.table import Table

linename = ['12CO','13CO','C18O']
restfreq = ['230.538GHz','220.39868GHz','219.56035GHz']

width    = '1.0km/s'
start    = '-11.0km/s'
nchan    = 42

fitspw   = '0,3,5,8,10,13,15,18'
linespw  = '1,2,4,6,7,9,11,12,14,16,17,19'
cell     = '0.03arcsec'
imsize   = [640,640]
robust   = 0.5

skip = np.array(['S29','S56','S74','S75','S76','S88'])
contdata = Table.read('/Users/mansdell/Dropbox/thesis/sori/alma/cont_results.txt',format='ascii.ipac')
ind_nd = np.where( contdata['F_cont']/contdata['E_cont'] >= 4.0 )
data = contdata[ind_nd]
dirs = np.array(data['name'])

dirs = np.array(['S55'])
linename = ['13CO']
restfreq = ['220.39868GHz']

for i,val in enumerate(dirs):

    if val in skip:
        print "\n>>> "+dirs[i]+" is known gas detection, skipping..."
        continue
    
    else:
        print "\n>>>Extracting "+dirs[i]

    field = int(dirs[i][1:])
    linevis = dirs[i]+'/S'+str(field)+'.vis.contsub'
    
    for n,val in enumerate(linename):

        print "\n    Analyzing "+linename[n]
        
        lineimg = dirs[i]+'/S'+str(field)+'_'+linename[n]
        linefreq = restfreq[n]

        maskfile = dirs[i]+'/mask_cont.crtf'
        if (os.path.isfile(maskfile) == False):
            print "    Missing mask file!"
            pdb.set_trace()

        if (os.path.isfile(lineimg+'.fits') == True):
            print "    FITS file exists, skipping"
            continue
        
        else:
            for ext in ['.flux','.image','.mask','.model','.psf','.residual']:
                if (os.path.isdir(lineimg+ext) == True): rmtables(lineimg+ext)

        print "      clearing visibility file..."
        clearcal(vis=linevis)
        delmod(vis=linevis)
            
        print "      performing clean..."
        clean(vis           = linevis,
              imagename     = lineimg,
              mode          = 'velocity',
              start         = start,
              width         = width,
              nchan         = nchan,
              # mask          = maskfile,
              outframe      = 'lsrk',
              veltype       = 'radio',
              restfreq      = linefreq,
              niter         = 300,
              threshold     = '0.0mJy',
              interactive   = False,
              imsize        = imsize,
              cell          = cell,
              weighting     ='briggs',
              robust        = robust,
              imagermode    = 'csclean')

        exportfits(imagename=lineimg+'.image',fitsimage=lineimg+'.fits')
        
        for ext in ['.flux','.image','.mask','.model','.psf','.residual']:
            if (os.path.isdir(lineimg+ext) == True): rmtables(lineimg+ext)
                
