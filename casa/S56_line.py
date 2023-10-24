#!/usr/bin/env


# ================================== Setup =====================================

field   = 56
visfile = '../calibrated_final.ms'
linevis = 'S'+str(field)+'.vis.contsub' 

fitspw  = '0,3,5,8,10,13,15,18'
linespw = '1,2,4,6,7,9,11,12,14,16,17,19'

cell   = '0.03arcsec'
imsize = [640,640]
robust = 0.5

xc,yc = 322,318


# ================ Create continuum subtracted line datasets ===================


fieldvis = 'S'+str(field)+'.vis'

split(vis = visfile,
      outputvis = fieldvis,
      field = field,
      datacolumn = 'data')

    
uvcontsub(vis = fieldvis,
          spw = linespw,
          fitspw = fitspw,
          excludechans=False,
          combine = 'spw',
          solint = 'int',
          fitorder = 1,
          want_cont = False) 

os.system('rm -rf '+fieldvis)


# ============================== 12CO line =====================================


### line info
linename = '12CO'
restfreq = '230.538GHz' 
width    = '1.0km/s'
start    = '-11.0km/s'
nchan    = 42
lineimg  = 'S'+str(field)+'_'+linename

### clear previous stuff
clearcal(vis=linevis)
delmod(vis=linevis)
for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage']:
    rmtables(lineimg+ext)

                        
### clean channels
clean(vis           = linevis,
      imagename     = lineimg,
      mode          = 'velocity',
      start         = start,
      width         = width,
      nchan         = nchan,
      outframe      = 'lsrk',
      veltype       = 'radio',
      restfreq      = restfreq,
      niter         = 2000,
      threshold     = '0.0mJy',
      interactive   = True,
      imsize        = imsize,
      cell          = cell,
      weighting     ='briggs',
      robust        = robust,
      imagermode    = 'csclean')


### maybe something 6-10 km/s?
imview(raster   = [{'file':lineimg+'.image'}],
        contour = [{'file':'S'+str(field)+'_cont.fits'}])

### export cube to fits file    
os.system('rm -f '+lineimg+'.fits')
exportfits(imagename=lineimg+'.image',fitsimage=lineimg+'.fits')
imview(raster   = [{'file':lineimg+'.fits'}],
        contour = [{'file':'S'+str(field)+'_cont.fits'}])


### delete stuff
for ext in ['.flux','.image','.mask','.model','.psf','.residual']: rmtables(lineimg+ext)
os.system('rm -rf cgrid_ft.im *.last *.log')

### make zero-moment map
os.system('rm -rf '+lineimg+'_mom0*')
immoments(imagename  = lineimg+'.fits',
          outfile    = lineimg+'_mom0',
          moments    = [0],
          includepix = [-10.0,1000.0],
          chans      = ('range=[6.0km/s,10.0km/s]'))
exportfits(imagename=lineimg+'_mom0',fitsimage=lineimg+'_mom0.fits')
os.system('rm -rf '+lineimg+'_mom0.image')
os.system('rm -rf '+lineimg+'_mom0')
imview(raster   = [{'file':lineimg+'_mom0.fits'}],
        contour = [{'file':'S'+str(field)+'_cont.fits'}])


### make first-moment map
os.system('rm -rf '+lineimg+'_mom1*')
immoments(imagename  = lineimg+'.fits',
          outfile    = lineimg+'_mom1',
          moments    = [1],
          excludepix = [-100,0.01],
          chans      = ('range=[6km/s,10km/s]'))
exportfits(imagename=lineimg+'_mom1',fitsimage=lineimg+'_mom1.fits')
imview(raster   = [{'file':lineimg+'_mom1.fits'}],
        contour = [{'file':'S'+str(field)+'_cont.fits'}])


### re-center image on source
os.system('rm -rf '+lineimg+'_mom0_crop.fits')
boxwidth = 150.
box = rg.box([xc-boxwidth,yc-boxwidth],[xc+boxwidth,yc+boxwidth])
ia.fromimage(outfile=lineimg+'_mom0_crop.image',infile=lineimg+'_mom0.fits',region=box )
ia.close() 
exportfits(imagename = lineimg+'_mom0_crop.image',fitsimage = lineimg+'_mom0_crop.fits')
os.system('rm -rf '+lineimg+'_mom0_crop.image')


### check image
imview(raster   = [{'file':lineimg+'_mom0_crop.fits'}],
        contour = [{'file':'S'+str(field)+'_cont.fits'}])


### use measure.py to get COG flux
'''
S56/S56_12CO_mom0_crop.fits
Annulus = [2.0, 5.0] arcseconds

pixel scale     0.03 arcsec
image noise    34.092 mJy
beam :  PA   -60    0.31 x  0.23 arcsec   #pix/beam:    90
center pixel    150,  150
Annulus is between 2.0 and 5.0 arcsec
Subtracting offset of -0.658 mJy/beam from image
RMS in annulus is 33.6

#   i  radius    #pix    total      rms     snr   beam_frac
#      (asec)            (mJy)     (mJy)
   0    0.000       1   153.00   34.455     4.4   0.011
   1    0.050       9   158.41   30.144     5.3   0.095
   2    0.100      37   176.25   27.637     6.4   0.335
   3    0.150      81   202.72   30.293     6.7   0.586
   4    0.200     137   232.85   32.079     7.3   0.773
   5    0.250     221   272.80   34.921     7.8   0.905
   6    0.300     317   314.09   41.667     7.5   0.964
   7    0.350     429   359.19   50.008     7.2   0.988
   8    0.400     553   406.74   49.492     8.2   0.996
   9    0.450     709   461.26   50.073     9.2   0.999
  10    0.500     877   500.82   73.512     6.8   1.000
  11    0.550    1049   513.57   58.056     8.8   1.000
  12    0.600    1257   495.27   99.510     5.0   1.000
  13    0.650    1481   447.92   93.831     4.8   1.000
  14    0.700    1709   389.58   92.275     4.2   1.000
  15    0.750    1961   335.35   88.290     3.8   1.000
  16    0.800    2233   306.36   74.216     4.1   1.000
  17    0.850    2537   307.33  114.017     2.7   1.000

PEAK FLUX: radius     total     rms    snr
             0.55    513.57   58.06    8.8
PEAK SNR:  radius     total     rms   snr
             0.45    461.26   50.07    9.2
'''




# ============================== 13CO line =====================================


### line info
linename = '13CO'
restfreq = '220.39868GHz' 
width    = '1.0km/s'
start    = '-11.0km/s'
nchan    = 42
lineimg  = 'S'+str(field)+'_'+linename

### clear previous stuff
clearcal(vis=linevis)
delmod(vis=linevis)
for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage']: rmtables(lineimg+ext)

                        
### clean channels
clean(vis           = linevis,
      imagename     = lineimg,
      mode          = 'velocity',
      start         = start,
      width         = width,
      nchan         = nchan,
      outframe      = 'lsrk',
      veltype       = 'radio',
      restfreq      = restfreq,
      niter         = 2000,
      threshold     = '0.0mJy',
      interactive   = True,
      imsize        = imsize,
      cell          = cell,
      weighting     ='briggs',
      robust        = robust,
      imagermode    = 'csclean')
# could not see line clearly, so just cleaned lightly 
# within continuum region for all channels showing 12CO emission


# can't see emission
imview(raster   = [{'file':lineimg+'.image'}],
        contour = [{'file':'S'+str(field)+'_cont.fits'}])


### export cube to fits file    
os.system('rm -f '+lineimg+'.fits')
exportfits(imagename=lineimg+'.image',fitsimage=lineimg+'.fits')


### make zero-moment map (using same channels as 12CO)
os.system('rm -rf '+lineimg+'_mom0*')
immoments(imagename  = lineimg+'.fits',
          outfile    = lineimg+'_mom0',
          moments    = [0],
          includepix = [-10.0,1000.0],
          chans      = ('range=[6.0km/s,10.0km/s]'))
exportfits(imagename=lineimg+'_mom0',fitsimage=lineimg+'_mom0.fits')
os.system('rm -rf '+lineimg+'_mom0.image')
os.system('rm -rf '+lineimg+'_mom0')

### delete stuff
for ext in ['.flux','.image','.mask','.model','.psf','.residual']: rmtables(lineimg+ext)
os.system('rm -rf cgrid_ft.im *.last *.log')

### re-center image on source
os.system('rm -rf '+lineimg+'_mom0_crop*')
boxwidth = 150.
box = rg.box([xc-boxwidth,yc-boxwidth],[xc+boxwidth,yc+boxwidth])
ia.fromimage(outfile=lineimg+'_mom0_crop.image',infile=lineimg+'_mom0.fits',region=box )
ia.close() 
exportfits(imagename = lineimg+'_mom0_crop.image',fitsimage = lineimg+'_mom0_crop.fits')
os.system('rm -rf '+lineimg+'_mom0_crop.image')


### check image
imview(raster   = [{'file':lineimg+'_mom0_crop.fits'}],
        contour = [{'file':'S'+str(field)+'_cont.fits'}])


### use measure.py to get COG flux
'''
S56/S56_13CO_mom0_crop.fits
Annulus = [2.0, 5.0] arcseconds

pixel scale     0.03 arcsec
image noise    35.403 mJy
beam :  PA   -60    0.32 x  0.25 arcsec   #pix/beam:   101
center pixel    150,  150
Annulus is between 2.0 and 5.0 arcsec
Subtracting offset of -0.676 mJy/beam from image
RMS in annulus is 35.5

#   i  radius    #pix    total      rms     snr   beam_frac
#      (asec)            (mJy)     (mJy)
   0    0.000       1    35.84   24.891     1.4   0.010
   1    0.050       9    35.87   35.516     1.0   0.085
   2    0.100      37    36.41   32.602     1.1   0.305
   3    0.150      81    39.65   37.716     1.1   0.544
   4    0.200     137    45.45   40.929     1.1   0.734
   5    0.250     221    57.64   47.717     1.2   0.878
   6    0.300     317    76.81   48.933     1.6   0.949
   7    0.350     429    97.08   51.754     1.9   0.981
   8    0.400     553   118.86   71.922     1.7   0.993
   9    0.450     709   139.01   68.361     2.0   0.998
  10    0.500     877   157.53   55.260     2.9   1.000
  11    0.550    1049   171.19   74.970     2.3   1.000
  12    0.600    1257   177.18   95.171     1.9   1.000
  13    0.650    1481   187.51   98.105     1.9   1.000
  14    0.700    1709   204.44   91.449     2.2   1.000
  15    0.750    1961   238.38   85.968     2.8   1.000
  16    0.800    2233   276.09  132.490     2.1   1.000
  17    0.850    2537   314.11  121.998     2.6   1.000
  18    0.900    2821   335.16   96.142     3.5   1.000
  19    0.950    3149   353.49  149.723     2.4   1.000

PEAK FLUX: radius     total     rms    snr
             0.95    353.49  149.72    2.4
PEAK SNR:  radius     total     rms   snr
             0.90    335.16   96.14    3.5
'''



# ============================== C18O line =====================================


### line info
linename = 'C18O'
restfreq = '219.56035GHz' 
width    = '1.0km/s'
start    = '-11.0km/s'
nchan    = 42
lineimg  = 'S'+str(field)+'_'+linename


### clear previous stuff
clearcal(vis=linevis)
delmod(vis=linevis)
for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage']: rmtables(lineimg+ext)

                        
### clean channels
clean(vis           = linevis,
      imagename     = lineimg,
      mode          = 'velocity',
      start         = start,
      width         = width,
      nchan         = nchan,
      outframe      = 'lsrk',
      veltype       = 'radio',
      restfreq      = restfreq,
      niter         = 2000,
      threshold     = '0.0mJy',
      interactive   = True,
      imsize        = imsize,
      cell          = cell,
      weighting     ='briggs',
      robust        = robust,
      imagermode    = 'csclean')
# could not see line clearly, so just cleaned lightly 
# within continuum region for all channels with 12CO emission


# can't see emission
imview(raster   = [{'file':lineimg+'.image'}],
        contour = [{'file':'S'+str(field)+'_cont.fits'}])


### export cube to fits file    
os.system('rm -f '+lineimg+'.fits')
exportfits(imagename=lineimg+'.image',fitsimage=lineimg+'.fits')


### make zero-moment map (using channels of 12CO)
os.system('rm -rf '+lineimg+'_mom0*')
immoments(imagename  = lineimg+'.fits',
          outfile    = lineimg+'_mom0',
          moments    = [0],
          includepix = [-10.0,1000.0],
          chans      = ('range=[6.0km/s,10.0km/s]'))
exportfits(imagename=lineimg+'_mom0',fitsimage=lineimg+'_mom0.fits')
os.system('rm -rf '+lineimg+'_mom0.image')
os.system('rm -rf '+lineimg+'_mom0')


### re-center image on source 
boxwidth = 150.
box = rg.box([xc-boxwidth,yc-boxwidth],[xc+boxwidth,yc+boxwidth])
ia.fromimage(outfile=lineimg+'_mom0_crop.image',infile=lineimg+'_mom0.fits',region=box )
ia.close() 
exportfits(imagename = lineimg+'_mom0_crop.image',fitsimage = lineimg+'_mom0_crop.fits')
os.system('rm -rf '+lineimg+'_mom0_crop.image')

### delete stuff
for ext in ['.flux','.image','.mask','.model','.psf','.residual']: rmtables(lineimg+ext)
os.system('rm -rf cgrid_ft.im *.last *.log')

### check image
imview(raster   = [{'file':lineimg+'_mom0_crop.fits'}],
        contour = [{'file':'S'+str(field)+'_cont.fits'}])


### use measure.py to get COG flux
'''

S56/S56_C18O_mom0_crop.fits
Annulus = [2.0, 5.0] arcseconds

pixel scale     0.03 arcsec
image noise    27.152 mJy
beam :  PA   -64    0.32 x  0.25 arcsec   #pix/beam:   103
center pixel    150,  150
Annulus is between 2.0 and 5.0 arcsec
Subtracting offset of 0.616 mJy/beam from image
RMS in annulus is 27.4

#   i  radius    #pix    total      rms     snr   beam_frac
#      (asec)            (mJy)     (mJy)
   0    0.000       1    20.48   26.013     0.8   0.010
   1    0.050       9    21.36   27.897     0.8   0.084
   2    0.100      37    24.01   22.949     1.0   0.300
   3    0.150      81    28.55   26.198     1.1   0.539
   4    0.200     137    32.32   29.679     1.1   0.728
   5    0.250     221    34.64   27.558     1.3   0.875
   6    0.300     317    35.91   37.144     1.0   0.947
   7    0.350     429    30.19   44.818     0.7   0.980
   8    0.400     553    22.15   36.947     0.6   0.993
   9    0.450     709    10.83   60.320     0.2   0.998
  10    0.500     877     1.55   65.123     0.0   1.000
  11    0.550    1049    -6.78   66.964    -0.1   1.000
  12    0.600    1257   -25.05   60.405    -0.4   1.000
  13    0.650    1481   -51.08   65.583    -0.8   1.000
  14    0.700    1709   -84.32   60.904    -1.4   1.000
  15    0.750    1961  -122.05   78.037    -1.6   1.000
  16    0.800    2233  -158.95  106.388    -1.5   1.000
  17    0.850    2537  -189.01  103.652    -1.8   1.000
  18    0.900    2821  -210.83  117.802    -1.8   1.000
  19    0.950    3149  -225.73  105.828    -2.1   1.000

PEAK FLUX: radius     total     rms    snr
             0.30     35.91   37.14    1.0
PEAK SNR:  radius     total     rms   snr
             0.25     34.64   27.56    1.3

'''
