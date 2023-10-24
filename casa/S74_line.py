#!/usr/bin/env


# ================================== Setup =====================================

field   = 74
visfile = '../calibrated_final.ms'
linevis = 'S'+str(field)+'.vis.contsub' 

fitspw  = '0,3,5,8,10,13,15,18'
linespw = '1,2,4,6,7,9,11,12,14,16,17,19'

cell   = '0.03arcsec'
imsize = [640,640]
robust = 0.5

xc,yc = 320,319


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

### can see 11-13 km/s
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
 
### make zero-moment map within bounds
### above 3xRMS in spectrum
os.system('rm -rf '+lineimg+'_mom0*')
immoments(imagename  = lineimg+'.fits',
          outfile    = lineimg+'_mom0',
          moments    = [0],
          includepix = [-10.0,1000.0],
          chans      = ('range=[11km/s,14km/s]'))
exportfits(imagename=lineimg+'_mom0',fitsimage=lineimg+'_mom0.fits')
os.system('rm -rf '+lineimg+'_mom0.image')
os.system('rm -rf '+lineimg+'_mom0')
imview(raster   = [{'file':lineimg+'_mom0.fits'}],
        contour = [{'file':'S'+str(field)+'_cont.fits'}])


### make first-moment map within bounds
os.system('rm -rf '+lineimg+'_mom1*')
immoments(imagename  = lineimg+'.fits',
          outfile    = lineimg+'_mom1',
          moments    = [1],
          excludepix = [-100,0.005],
          chans      = ('range=[11km/s,14km/s]'))
exportfits(imagename=lineimg+'_mom1',fitsimage=lineimg+'_mom1.fits')
imview(raster   = [{'file':lineimg+'_mom1.fits'}],
        contour = [{'file':'S'+str(field)+'_cont.fits'}])


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
S74/S74_12CO_mom0_crop.fits
S74/S74_12CO_mom0_crop.fits
Annulus = [2.0, 5.0] arcseconds

pixel scale     0.03 arcsec
image noise    27.100 mJy
beam :  PA   -61    0.32 x  0.23 arcsec   #pix/beam:    94
center pixel    150,  150
Annulus is between 2.0 and 5.0 arcsec
Subtracting offset of 1.575 mJy/beam from image
RMS in annulus is 27

#   i  radius    #pix    total      rms     snr   beam_frac
#      (asec)            (mJy)     (mJy)
   0    0.000       1    83.89   28.140     3.0   0.011
   1    0.050       9    87.26   26.778     3.3   0.092
   2    0.100      37    98.63   25.850     3.8   0.325
   3    0.150      81   116.22   29.789     3.9   0.571
   4    0.200     137   136.14   28.000     4.9   0.758
   5    0.250     221   159.35   33.963     4.7   0.894
   6    0.300     317   170.81   33.411     5.1   0.958
   7    0.350     429   170.19   42.817     4.0   0.985
   8    0.400     553   159.57   47.470     3.4   0.995
   9    0.450     709   146.94   54.953     2.7   0.999
  10    0.500     877   145.25   61.983     2.3   1.000
  11    0.550    1049   152.37   64.717     2.4   1.000
  12    0.600    1257   168.39   74.397     2.3   1.000
  13    0.650    1481   175.45   72.428     2.4   1.000
  14    0.700    1709   171.75   76.383     2.2   1.000
  15    0.750    1961   152.11   76.599     2.0   1.000
  16    0.800    2233   123.16   70.879     1.7   1.000
  17    0.850    2537    81.49  105.157     0.8   1.000
  18    0.900    2821    45.44   80.191     0.6   1.000
  19    0.950    3149    12.68   54.699     0.2   1.000

PEAK FLUX: radius     total     rms    snr
             0.65    175.45   72.43    2.4
PEAK SNR:  radius     total     rms   snr
             0.30    170.81   33.41    5.1
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
# in channels with 12CO emission in continuum region


# some emission seen after cleaning in spectrum
imview(raster   = [{'file':lineimg+'.image'}],
        contour = [{'file':'S'+str(field)+'_cont.fits'}])


### export cube to fits file    
os.system('rm -f '+lineimg+'.fits')
exportfits(imagename=lineimg+'.image',fitsimage=lineimg+'.fits')


### delete stuff
for ext in ['.flux','.image','.mask','.model','.psf','.residual']: rmtables(lineimg+ext)
os.system('rm -rf cgrid_ft.im *.last *.log')


### make zero-moment map (not clear emission)
os.system('rm -rf '+lineimg+'_mom0*')
immoments(imagename  = lineimg+'.fits',
          outfile    = lineimg+'_mom0',
          moments    = [0],
          includepix = [-10.0,1000.0],
          chans      = ('range=[11km/s,14km/s]'))
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


### check image
imview(raster   = [{'file':lineimg+'_mom0_crop.fits'}],
        contour = [{'file':'S'+str(field)+'_cont.fits'}])


### use measure.py to get COG flux
'''
S74/S74_13CO_mom0_crop.fits
Annulus = [2.0, 5.0] arcseconds

pixel scale     0.03 arcsec
image noise    30.660 mJy
beam :  PA   -60    0.33 x  0.25 arcsec   #pix/beam:   105
center pixel    150,  150
Annulus is between 2.0 and 5.0 arcsec
Subtracting offset of 0.840 mJy/beam from image
RMS in annulus is 30.9

#   i  radius    #pix    total      rms     snr   beam_frac
#      (asec)            (mJy)     (mJy)
   0    0.000       1    29.19   30.652     1.0   0.010
   1    0.050       9    29.36   35.420     0.8   0.082
   2    0.100      37    30.11   35.344     0.9   0.296
   3    0.150      81    32.47   27.890     1.2   0.532
   4    0.200     137    36.76   35.747     1.0   0.720
   5    0.250     221    46.75   29.167     1.6   0.867
   6    0.300     317    62.81   36.055     1.7   0.942
   7    0.350     429    81.97   39.147     2.1   0.977
   8    0.400     553    99.81   42.026     2.4   0.992
   9    0.450     709   108.87   46.906     2.3   0.998
  10    0.500     877   107.28   54.960     2.0   0.999
  11    0.550    1049   100.80   59.606     1.7   1.000
  12    0.600    1257    90.00   48.534     1.9   1.000
  13    0.650    1481    87.23   67.284     1.3   1.000
  14    0.700    1709    85.79   82.548     1.0   1.000
  15    0.750    1961    75.05   63.071     1.2   1.000
  16    0.800    2233    55.99  106.793     0.5   1.000
  17    0.850    2537    30.26   91.002     0.3   1.000
  18    0.900    2821     9.68  130.246     0.1   1.000
  19    0.950    3149     0.01   93.621     0.0   1.000

PEAK FLUX: radius     total     rms    snr
             0.45    108.87   46.91    2.3
PEAK SNR:  radius     total     rms   snr
             0.40     99.81   42.03    2.4
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
# in channels with 12CO emission in continuum region


# still cannot see emission
imview(raster   = [{'file':lineimg+'.image'}],
        contour = [{'file':'S'+str(field)+'_cont.fits'}])


### export cube to fits file    
os.system('rm -f '+lineimg+'.fits')
exportfits(imagename=lineimg+'.image',fitsimage=lineimg+'.fits')


### delete stuff
for ext in ['.flux','.image','.mask','.model','.psf','.residual']: rmtables(lineimg+ext)
os.system('rm -rf cgrid_ft.im *.last *.log')

### make zero-moment map
os.system('rm -rf '+lineimg+'_mom0*')
immoments(imagename  = lineimg+'.fits',
          outfile    = lineimg+'_mom0',
          moments    = [0],
          includepix = [-10.0,1000.0],
          chans      = ('range=[11km/s,14km/s]'))
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


### check image
imview(raster   = [{'file':lineimg+'_mom0_crop.fits'}],
        contour = [{'file':'S'+str(field)+'_cont.fits'}])


### use measure.py to get COG flux
'''
S74/S74_C18O_mom0_crop.fits
Annulus = [2.0, 5.0] arcseconds

pixel scale     0.03 arcsec
image noise    22.151 mJy
beam :  PA   -63    0.33 x  0.26 arcsec   #pix/beam:   107
center pixel    150,  150
Annulus is between 2.0 and 5.0 arcsec
Subtracting offset of 0.082 mJy/beam from image
RMS in annulus is 22.1

#   i  radius    #pix    total      rms     snr   beam_frac
#      (asec)            (mJy)     (mJy)
   0    0.000       1    24.23   20.487     1.2   0.009
   1    0.050       9    25.44   24.735     1.0   0.081
   2    0.100      37    29.34   24.594     1.2   0.292
   3    0.150      81    34.81   24.637     1.4   0.527
   4    0.200     137    40.40   23.891     1.7   0.716
   5    0.250     221    46.04   27.934     1.6   0.865
   6    0.300     317    50.90   26.160     1.9   0.941
   7    0.350     429    53.48   37.581     1.4   0.977
   8    0.400     553    55.21   43.497     1.3   0.992
   9    0.450     709    60.12   40.709     1.5   0.998
  10    0.500     877    64.29   54.595     1.2   0.999
  11    0.550    1049    67.71   45.237     1.5   1.000
  12    0.600    1257    65.87   72.293     0.9   1.000
  13    0.650    1481    60.69   50.900     1.2   1.000
  14    0.700    1709    50.69   69.883     0.7   1.000
  15    0.750    1961    39.56   69.397     0.6   1.000
  16    0.800    2233    28.57   63.444     0.5   1.000
  17    0.850    2537    18.39   87.092     0.2   1.000
  18    0.900    2821     8.34  100.177     0.1   1.000
  19    0.950    3149    -6.79  103.218    -0.1   1.000

PEAK FLUX: radius     total     rms    snr
             0.55     67.71   45.24    1.5
PEAK SNR:  radius     total     rms   snr
             0.30     50.90   26.16    1.9
'''
