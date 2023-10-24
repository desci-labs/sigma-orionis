#!/usr/bin/env

import numpy as np
import sys,os,pdb,glob
from astropy.table import Table

linename = ['12CO','13CO','C18O']
boxwidth = 150.
pscale = 0.03

skip = np.array(['S29','S56','S74','S75','S76','S88'])
contdata = Table.read('/Users/mansdell/Dropbox/thesis/sori/alma/cont_results.txt',format='ascii.ipac')

# ind_nd = np.where( contdata['F_cont']/contdata['E_cont'] < 4.0 )
# data = contdata[ind_nd]
# data = data[np.where(data['name'] == 'S86')]
data = contdata

dirs = np.array(data['name'])
for i,val in enumerate(dirs):

    if val in skip:
        print "\n>>> "+dirs[i]+" is known gas detection, skipping..."
        continue
    else:
        print "\n>>> Making maps for "+dirs[i]+'\n'

    field = int(dirs[i][1:])
    xc,yc = 320. - data['dx'][i]/pscale, 320. + data['dy'][i]/pscale
    print "    Source location = {0:0.1f}".format(xc)+' '+"{0:0.1f}".format(yc)+'\n'

    for n,val in enumerate(linename):

        
        lineimg = dirs[i]+'/S'+str(field)+'_'+linename[n]

        if (os.path.isfile(lineimg+'_mom0_g2.fits') == True):
            print "    FITS file exists, skipping"
            continue
        else:
            immoments(imagename  = lineimg+'.fits',
                      outfile    = lineimg+'_mom0',
                      moments    = [0],
                      includepix = [-10.0,1000.0],
                      chans      = ('range=[29km/s,33km/s]'))
            exportfits(imagename=lineimg+'_mom0',fitsimage=lineimg+'_mom0_g2.fits')
            os.system('rm -rf '+lineimg+'_mom0')
            box = rg.box([xc-boxwidth,yc-boxwidth],[xc+boxwidth,yc+boxwidth])
            ia.fromimage(outfile=lineimg+'_mom0_g2_crop.image',infile=lineimg+'_mom0_g2.fits',region=box )
            ia.close() 
            exportfits(imagename = lineimg+'_mom0_g2_crop.image',fitsimage = lineimg+'_mom0_g2_crop.fits')
            os.system('rm -rf '+lineimg+'_mom0_g2_crop.image')


            
        if (os.path.isfile(lineimg+'_mom0_g1.fits') == True):
            print "    FITS file exists, skipping"
            continue
        else:
            immoments(imagename  = lineimg+'.fits',
                      outfile    = lineimg+'_mom0',
                      moments    = [0],
                      includepix = [-10.0,1000.0],
                      chans      = ('range=[22km/s,26km/s]'))
            exportfits(imagename=lineimg+'_mom0',fitsimage=lineimg+'_mom0_g1.fits')
            os.system('rm -rf '+lineimg+'_mom0')
            box = rg.box([xc-boxwidth,yc-boxwidth],[xc+boxwidth,yc+boxwidth])
            ia.fromimage(outfile=lineimg+'_mom0_g1_crop.image',infile=lineimg+'_mom0_g1.fits',region=box )
            ia.close() 
            exportfits(imagename = lineimg+'_mom0_g1_crop.image',fitsimage = lineimg+'_mom0_g1_crop.fits')
            os.system('rm -rf '+lineimg+'_mom0_g1_crop.image')

        print "    "+linename[n] +" done"

