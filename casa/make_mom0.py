#!/usr/bin/env

import numpy as np
import sys,os,pdb,glob
from astropy.table import Table
import get_lsrk

linename = ['12CO','13CO','C18O']
boxwidth = 150.
pscale = 0.03

skip = np.array(['S29','S56','S74','S75','S76','S88'])
contdata = Table.read('/Users/mansdell/Dropbox/thesis/sori/catalog/sori_cat.fits',format='fits')

# ind_nd = np.where( contdata['F_cont']/contdata['E_cont'] < 4.0 )
# data = contdata[ind_nd]
# data = data[np.where(data['name'] == 'S86')]
data = contdata
data.sort('name')

dirs = np.array(data['name'])
dirs = np.array(['S55'])
linename = ['13CO']

for i,val in enumerate(dirs):

    if val in skip:
        print "\n>>> "+dirs[i]+" is known gas detection, skipping..."
        continue
    else:
        print "\n>>> Making maps for "+dirs[i]+'\n'

    field = int(dirs[i][1:])

    ind = np.where(data['name'] == dirs[i])
    if len(ind[0]) == 0:
        print "cannot find field in sigOri catalog"
        pdb.set_trace()
    
    xc,yc = 320. - float(data['dx'][ind]/pscale), 320. + float(data['dy'][ind]/pscale)
    print "    Source location = {0:0.2f}".format(xc)+' '+"{0:0.2f}".format(yc)

    for n,val in enumerate(linename):

        
        lineimg = dirs[i]+'/S'+str(field)+'_'+linename[n]

        if (os.path.isfile(lineimg+'_mom0.fits') == True):
            
            print "    FITS file exists, skipping"
            continue
        
        else:

            
            if data['RV'][ind] != -99.0:
                corr = get_lsrk.get_corr(float(data['ra'][ind]),float(data['de'][ind]))
                rv_lsrk = float(data['RV'][ind]) + corr
                vmin,vmax = rv_lsrk-1,rv_lsrk+1
            else:
                vmin,vmax = 12,14
            if n==0: print "    Velocity range = {0:0.2f}".format(vmin)+' '+"{0:0.2f}".format(vmax)+'\n'
                      
            immoments(imagename  = lineimg+'.fits',
                      outfile    = lineimg+'_mom0',
                      moments    = [0],
                      includepix = [-10.0,1000.0],
                      chans      = ('range=['+str(vmin)+'km/s,'+str(vmax)+'km/s]'))
            
            exportfits(imagename=lineimg+'_mom0',fitsimage=lineimg+'_mom0.fits')
            os.system('rm -rf '+lineimg+'_mom0')
            
            box = rg.box([xc-boxwidth,yc-boxwidth],[xc+boxwidth,yc+boxwidth])
            ia.fromimage(outfile=lineimg+'_mom0_crop.image',infile=lineimg+'_mom0.fits',region=box )
            ia.close() 
            exportfits(imagename = lineimg+'_mom0_crop.image',fitsimage = lineimg+'_mom0_crop.fits')
            os.system('rm -rf '+lineimg+'_mom0_crop.image')

        print "    "+linename[n] +" done"

