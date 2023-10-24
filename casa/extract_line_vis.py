#!/usr/bin/env

import numpy as np
import sys,os,pdb,glob


visfile = 'calibrated_final.ms'

fitspw = '23,29,84,90,114,120,146,152,166,172,218,224,265,271,302,308'
linespw = '25,27,31,86,87,92,116,117,122,148,149,154,168,169,174,200,221,226,267,268,273,304,305,310'


dirs = os.listdir(".")
for i,val in enumerate(dirs):

    if val[0] == 'S':

        field = val[1:]
        
        fieldvis = val+'/S'+field+'.vis'
        linevis  = val+'/S'+field+'.vis.contsub'

        pdb.set_trace()

        if (os.path.isdir(linevis) == False):

            print "\n>>>Extracting "+val

            print "    field visfile"
            split(vis = visfile,
                  outputvis = fieldvis,
                  field = field,
                  datacolumn = 'data')
            
            print "    uvcontsub visfile"
            uvcontsub(vis = fieldvis,
                      spw = linespw,
                      fitspw = fitspw,
                      excludechans=False,
                      combine = 'spw',
                      solint = 'int',
                      fitorder = 1,
                      want_cont = False) 

            os.system('rm -rf '+fieldvis)
            print "    Done!"

        else:

            print "\n!Skipping "+val
