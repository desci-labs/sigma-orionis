#!/usr/bin/env


# ================================== Setup =====================================

field   = 75
visfile = '../calibrated_final.ms'
linevis = 'S'+str(field)+'.vis.contsub' 

fitspw  = '0,3,5,8,10,13,15,18'
linespw = '1,2,4,6,7,9,11,12,14,16,17,19'

cell   = '0.03arcsec'
imsize = [640,640]
robust = 0.5

xc,yc = 321,320


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


### can see emission 10-14
imview(raster   = [{'file':lineimg+'.image'}],
        contour = [{'file':'S'+str(field)+'_cont.fits'}])

### export cube to fits file    
os.system('rm -f '+lineimg+'.fits')
exportfits(imagename=lineimg+'.image',fitsimage=lineimg+'.fits')
imview(raster   = [{'file':lineimg+'.fits'}],
        contour = [{'file':'S'+str(field)+'_cont.fits'}])


### make zero-moment map
os.system('rm -rf '+lineimg+'_mom0*')
immoments(imagename  = lineimg+'.fits',
          outfile    = lineimg+'_mom0',
          moments    = [0],
          includepix = [-10.0,1000.0],
          chans      = ('range=[10km/s,14km/s]'))
exportfits(imagename=lineimg+'_mom0',fitsimage=lineimg+'_mom0.fits')
os.system('rm -rf '+lineimg+'_mom0.image')
imview(raster   = [{'file':lineimg+'_mom0.fits'}],
        contour = [{'file':'S'+str(field)+'_cont.fits'}])

### delete stuff
for ext in ['.flux','.image','.mask','.model','.psf','.residual']: rmtables(lineimg+ext)
os.system('rm -rf cgrid_ft.im *.last *.log')


### make first-moment map within bounds
os.system('rm -rf '+lineimg+'_mom1*')
immoments(imagename  = lineimg+'.fits',
          outfile    = lineimg+'_mom1',
          moments    = [1],
          excludepix = [-100,0.015],
          chans      = ('range=[10km/s,14km/s]'))
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
S75/S75_12CO_mom0_crop.fits
Annulus = [2.0, 5.0] arcseconds

pixel scale     0.03 arcsec
image noise    32.250 mJy
beam :  PA   -60    0.32 x  0.24 arcsec   #pix/beam:    94
center pixel    150,  150
Annulus is between 2.0 and 5.0 arcsec
Subtracting offset of 1.410 mJy/beam from image
RMS in annulus is 31.6

#   i  radius    #pix    total      rms     snr   beam_frac
#      (asec)            (mJy)     (mJy)
   0    0.000       1   148.16   30.491     4.9   0.011
   1    0.050       9   154.56   33.032     4.7   0.091
   2    0.100      37   177.31   30.242     5.9   0.324
   3    0.150      81   216.12   35.924     6.0   0.570
   4    0.200     137   267.19   33.456     8.0   0.758
   5    0.250     221   343.52   36.240     9.5   0.894
   6    0.300     317   420.33   41.021    10.2   0.958
   7    0.350     429   483.80   52.889     9.1   0.985
   8    0.400     553   524.14   47.840    11.0   0.995
   9    0.450     709   552.42   61.570     9.0   0.999
  10    0.500     877   573.05   68.881     8.3   1.000
  11    0.550    1049   597.79   90.530     6.6   1.000
  12    0.600    1257   632.73   82.203     7.7   1.000 <--good
  13    0.650    1481   668.69   90.328     7.4   1.000
  14    0.700    1709   702.59   80.148     8.8   1.000
  15    0.750    1961   739.30   98.495     7.5   1.000
  16    0.800    2233   778.80  115.910     6.7   1.000
  17    0.850    2537   820.06  104.490     7.8   1.000
  18    0.900    2821   851.25  151.670     5.6   1.000

PEAK FLUX: radius     total     rms    snr
             1.35    940.69  257.71    3.7
PEAK SNR:  radius     total     rms   snr
             1.15    935.85   24.68   37.9


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
# within continuum region for all channels with 12 CO emission


# can see emission 12-14 in spectrum
imview(raster   = [{'file':lineimg+'.image'}],
        contour = [{'file':'S'+str(field)+'_cont.fits'}])


### export cube to fits file    
os.system('rm -f '+lineimg+'.fits')
exportfits(imagename=lineimg+'.image',fitsimage=lineimg+'.fits')


### make zero-moment map
os.system('rm -rf '+lineimg+'_mom0*')
immoments(imagename  = lineimg+'.fits',
          outfile    = lineimg+'_mom0',
          moments    = [0],
          includepix = [-10.0,1000.0],
          chans      = ('range=[10km/s,14km/s]'))
exportfits(imagename=lineimg+'_mom0',fitsimage=lineimg+'_mom0.fits')
os.system('rm -rf '+lineimg+'_mom0.image')
## same flux as when using 12CO boundaries

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
S75/S75_13CO_mom0_crop.fits
Annulus = [2.0, 5.0] arcseconds

pixel scale     0.03 arcsec
image noise    35.071 mJy
beam :  PA   -59    0.33 x  0.25 arcsec   #pix/beam:   105
center pixel    150,  150
Annulus is between 2.0 and 5.0 arcsec
Subtracting offset of 0.040 mJy/beam from image
RMS in annulus is 34.8

#   i  radius    #pix    total      rms     snr   beam_frac
#      (asec)            (mJy)     (mJy)
   0    0.000       1    79.75   32.204     2.5   0.010
   1    0.050       9    80.91   32.830     2.5   0.082
   2    0.100      37    85.39   35.164     2.4   0.295
   3    0.150      81    94.88   40.530     2.3   0.531
   4    0.200     137   108.62   43.337     2.5   0.719
   5    0.250     221   132.23   48.346     2.7   0.866
   6    0.300     317   162.08   49.605     3.3   0.942
   7    0.350     429   192.18   57.137     3.4   0.977
   8    0.400     553   220.62   54.343     4.1   0.992
   9    0.450     709   252.48   64.613     3.9   0.998
  10    0.500     877   280.09   66.530     4.2   0.999
  11    0.550    1049   300.24   59.233     5.1   1.000
  12    0.600    1257   309.18   65.453     4.7   1.000
  13    0.650    1481   307.96   80.673     3.8   1.000
  14    0.700    1709   300.54   92.658     3.2   1.000

PEAK FLUX: radius     total     rms    snr
             0.60    309.18   65.45    4.7
PEAK SNR:  radius     total     rms   snr
             0.55    300.24   59.23    5.1
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
# within continuum region for all channels with 12 CO emission


# can't see emission
imview(raster   = [{'file':lineimg+'.image'}],
        contour = [{'file':'S'+str(field)+'_cont.fits'}])


### export cube to fits file    
os.system('rm -f '+lineimg+'.fits')
exportfits(imagename=lineimg+'.image',fitsimage=lineimg+'.fits')


### make zero-moment map (same as 13 CO)
os.system('rm -rf '+lineimg+'_mom0*')
immoments(imagename  = lineimg+'.fits',
          outfile    = lineimg+'_mom0',
          moments    = [0],
          includepix = [-10.0,1000.0],
          chans      = ('range=[12km/s,14km/s]'))
exportfits(imagename=lineimg+'_mom0',fitsimage=lineimg+'_mom0.fits')
os.system('rm -rf '+lineimg+'_mom0.image')

### delete stuff
for ext in ['.flux','.image','.mask','.model','.psf','.residual']: rmtables(lineimg+ext)
os.system('rm -rf cgrid_ft.im *.last *.log')


### re-center image on source 
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

S75/S75_C18O_mom0_crop.fits
Annulus = [2.0, 5.0] arcseconds

pixel scale     0.03 arcsec
image noise    26.108 mJy
beam :  PA   -62    0.33 x  0.26 arcsec   #pix/beam:   107
center pixel    150,  150
Annulus is between 2.0 and 5.0 arcsec
Subtracting offset of 0.402 mJy/beam from image
RMS in annulus is 25.7

#   i  radius    #pix    total      rms     snr   beam_frac
#      (asec)            (mJy)     (mJy)
   0    0.000       1    11.89   20.998     0.6   0.009
   1    0.050       9    14.35   26.411     0.5   0.081
   2    0.100      37    22.26   24.258     0.9   0.292
   3    0.150      81    33.15   23.122     1.4   0.527
   4    0.200     137    43.79   29.235     1.5   0.717
   5    0.250     221    53.05   32.757     1.6   0.866
   6    0.300     317    57.84   33.019     1.8   0.942
   7    0.350     429    56.31   44.547     1.3   0.978
   8    0.400     553    52.31   43.598     1.2   0.992
   9    0.450     709    56.31   54.745     1.0   0.998
  10    0.500     877    62.85   61.587     1.0   0.999
  11    0.550    1049    71.94   56.310     1.3   1.000
  12    0.600    1257    84.76   63.689     1.3   1.000
  13    0.650    1481    91.73  106.790     0.9   1.000
  14    0.700    1709    94.77   97.546     1.0   1.000

PEAK FLUX: radius     total     rms    snr
             0.70     94.77   97.55    1.0
PEAK SNR:  radius     total     rms   snr
             0.30     57.84   33.02    1.8
'''
