#!/usr/bin/env


# ================================== Setup =====================================

field   = 29
visfile = '../calibrated_final.ms'
linevis = 'S'+str(field)+'.vis.contsub' 

fitspw  = '0,3,5,8,10,13,15,18'
linespw = '1,2,4,6,7,9,11,12,14,16,17,19'

cell   = '0.03arcsec'
imsize = [640,640]
robust = 0.5

xc,yc = 319,318


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


# emission seen cleary at 5-11 km/s
# spectrum shows 5-12 for dropping below 3xRMS
# unknown RV for this source 
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
### tested different bounds, this gave highest snr
### and same flux as when going out to 13 km/s
os.system('rm -rf '+lineimg+'_mom0*')
immoments(imagename  = lineimg+'.fits',
          outfile    = lineimg+'_mom0',
          moments    = [0],
          includepix = [-10.0,1000.0],
          chans      = ('range=[5km/s,12km/s]'))
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
          excludepix = [-100,0.015],
          chans      = ('range=[5km/s,12km/s]'))
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

S29/S29_12CO_mom0_crop.fits
Annulus = [2.0, 5.0] arcseconds

pixel scale     0.03 arcsec
image noise    39.241 mJy
beam :  PA   -59    0.30 x  0.24 arcsec   #pix/beam:    90
center pixel    150,  150
Annulus is between 2.0 and 5.0 arcsec
Subtracting offset of -0.854 mJy/beam from image
RMS in annulus is 38

#   i  radius    #pix    total      rms     snr   beam_frac
#      (asec)            (mJy)     (mJy)
   0    0.000       1   216.36   37.185     5.8   0.011
   1    0.050       9   222.50   38.650     5.8   0.095
   2    0.100      37   244.03   41.330     5.9   0.337
   3    0.150      81   279.92   32.534     8.6   0.589
   4    0.200     137   328.74   38.804     8.5   0.778
   5    0.250     221   409.50   49.908     8.2   0.909
   6    0.300     317   508.94   51.188     9.9   0.967
   7    0.350     429   621.44   59.646    10.4   0.990
   8    0.400     553   729.71   47.877    15.2   0.997
   9    0.450     709   844.81   88.263     9.6   0.999
  10    0.500     877   932.73   81.465    11.4   1.000
  11    0.550    1049  1003.25   91.076    11.0   1.000
  12    0.600    1257  1078.39   83.048    13.0   1.000
  13    0.650    1481  1142.40   95.455    12.0   1.000
  14    0.700    1709  1186.97   88.086    13.5   1.000
  15    0.750    1961  1202.94  130.006     9.3   1.000
  16    0.800    2233  1204.49   85.179    14.1   1.000
  17    0.850    2537  1194.58  129.376     9.2   1.000

PEAK FLUX: radius     total     rms    snr
             0.80   1204.49   85.18   14.1
PEAK SNR:  radius     total     rms   snr
             0.40    729.71   47.88   15.2

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


# some emission seen after cleaning in spectrum, 6-11 km/s
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
          chans      = ('range=[6km/s,11km/s]'))
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
S29/S29_13CO_mom0_crop.fits
Annulus = [2.0, 5.0] arcseconds

pixel scale     0.03 arcsec
image noise    35.809 mJy
beam :  PA   -63    0.31 x  0.25 arcsec   #pix/beam:    99
center pixel    150,  150
Annulus is between 2.0 and 5.0 arcsec
Subtracting offset of -0.904 mJy/beam from image
RMS in annulus is 35.8

#   i  radius    #pix    total      rms     snr   beam_frac
#      (asec)            (mJy)     (mJy)
   0    0.000       1   111.37   39.157     2.8   0.010
   1    0.050       9   116.74   42.980     2.7   0.087
   2    0.100      37   134.63   38.999     3.5   0.312
   3    0.150      81   160.75   31.921     5.0   0.556
   4    0.200     137   190.04   37.653     5.0   0.746
   5    0.250     221   225.24   47.464     4.7   0.888
   6    0.300     317   249.18   44.363     5.6   0.955
   7    0.350     429   265.74   44.940     5.9   0.984
   8    0.400     553   274.36   66.864     4.1   0.995
   9    0.450     709   275.85   53.656     5.1   0.999
  10    0.500     877   275.36   83.771     3.3   1.000
  11    0.550    1049   271.92   84.166     3.2   1.000
  12    0.600    1257   267.96   71.978     3.7   1.000
  13    0.650    1481   262.79  122.449     2.1   1.000
  14    0.700    1709   258.29  100.753     2.6   1.000

PEAK FLUX: radius     total     rms    snr
             0.45    275.85   53.66    5.1
PEAK SNR:  radius     total     rms   snr
             0.35    265.74   44.94    5.9

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


# can maybe see emission where expected
imview(raster   = [{'file':lineimg+'.image'}],
        contour = [{'file':'S'+str(field)+'_cont.fits'}])


### export cube to fits file    
os.system('rm -f '+lineimg+'.fits')
exportfits(imagename=lineimg+'.image',fitsimage=lineimg+'.fits')


### delete stuff
for ext in ['.flux','.image','.mask','.model','.psf','.residual']: rmtables(lineimg+ext)
os.system('rm -rf cgrid_ft.im *.last *.log')

### make zero-moment map using same limits as 13CO since can't see emission
os.system('rm -rf '+lineimg+'_mom0*')
immoments(imagename  = lineimg+'.fits',
          outfile    = lineimg+'_mom0',
          moments    = [0],
          includepix = [-10.0,1000.0],
          chans      = ('range=[6km/s,11km/s]'))
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
S29/S29_C18O_mom0_crop.fits
Annulus = [2.0, 5.0] arcseconds

pixel scale     0.03 arcsec
image noise    29.546 mJy
beam :  PA   -62    0.32 x  0.26 arcsec   #pix/beam:   102
center pixel    150,  150
Annulus is between 2.0 and 5.0 arcsec
Subtracting offset of -1.122 mJy/beam from image
RMS in annulus is 29.3

#   i  radius    #pix    total      rms     snr   beam_frac
#      (asec)            (mJy)     (mJy)
   0    0.000       1     5.30   30.700     0.2   0.010
   1    0.050       9     5.68   27.592     0.2   0.085
   2    0.100      37     6.79   26.918     0.3   0.304
   3    0.150      81     8.25   31.061     0.3   0.544
   4    0.200     137     9.21   30.060     0.3   0.735
   5    0.250     221     8.96   27.892     0.3   0.880
   6    0.300     317     9.00   39.382     0.2   0.951
   7    0.350     429     7.43   47.018     0.2   0.982
   8    0.400     553     6.65   53.573     0.1   0.994
   9    0.450     709    11.34   58.601     0.2   0.999
  10    0.500     877    18.76   92.350     0.2   1.000
  11    0.550    1049    29.21   74.666     0.4   1.000
  12    0.600    1257    44.67   94.621     0.5   1.000
  13    0.650    1481    61.63  113.212     0.5   1.000
  14    0.700    1709    80.75   84.279     1.0   1.000
  15    0.750    1961   105.72  107.763     1.0   1.000
  16    0.800    2233   131.09  113.431     1.2   1.000
  17    0.850    2537   162.08  119.088     1.4   1.000

PEAK FLUX: radius     total     rms    snr
             0.85    162.08  119.09    1.4
PEAK SNR:  radius     total     rms   snr
             0.85    162.08  119.09    1.4

'''
