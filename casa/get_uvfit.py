#!/usr/bin/env

import numpy as np
import sys,os,pdb,glob
from astropy.table import Table


data = Table.read('/Users/meganansdell/Dropbox/thesis/sori/alma/cont_results.txt',format='ascii.ipac')
dirs = np.array(data['name'])

for i,val in enumerate(dirs):

    
    print "\n\n>>>Extracting "+dirs[i]

    field = int(dirs[i][1:])
    contvis = dirs[i]+'/S'+str(field)+'_cont.vis'

    flx = 1e-3*data['F_cont'][i]
    dx  = data['dx'][i]
    dy  = data['dy'][i]

    if data['F_cont'][i]/data['E_cont'][i] > 3.0:
        varypar = [T,T,T]
    else:
        varypar = [T,F,F]
     
    uvmodelfit(vis       = contvis,
               comptype  = 'P',
               sourcepar = [flx,dx,dy],
               varypar   = varypar,
               niter     = 10)
