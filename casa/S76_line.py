#!/usr/bin/env


# ================================== Setup =====================================

field   = 76
visfile = '../calibrated_final.ms'
linevis = 'S'+str(field)+'.vis.contsub' 

fitspw  = '0,3,5,8,10,13,15,18'
linespw = '1,2,4,6,7,9,11,12,14,16,17,19'

cell   = '0.03arcsec'
imsize = [640,640]
robust = 0.5

xc,yc = 316,317


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

### can see emission ~11-15/16 (confirmed in spectrum)
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
          chans      = ('range=[11.0km/s,15.0km/s]'))
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
          excludepix = [-100,0.012],
          chans      = ('range=[11km/s,15km/s]'))
exportfits(imagename=lineimg+'_mom1',fitsimage=lineimg+'_mom1.fits')
imview(raster   = [{'file':lineimg+'_mom1.fits'}],
        contour = [{'file':'S'+str(field)+'_cont.fits'}])



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
S76/S76_12CO_mom0_crop.fits
Annulus = [2.0, 5.0] arcseconds

pixel scale     0.03 arcsec
image noise    30.395 mJy
beam :  PA   -61    0.32 x  0.24 arcsec   #pix/beam:    94
center pixel    150,  150
Annulus is between 2.0 and 5.0 arcsec
Subtracting offset of -0.870 mJy/beam from image
RMS in annulus is 29.3

#   i  radius    #pix    total      rms     snr   beam_frac
#      (asec)            (mJy)     (mJy)
   0    0.000       1   180.22   32.916     5.5   0.011
   1    0.050       9   185.17   26.320     7.0   0.091
   2    0.100      37   202.37   25.433     8.0   0.323
   3    0.150      81   230.76   31.286     7.4   0.569
   4    0.200     137   267.92   33.553     8.0   0.756
   5    0.250     221   325.83   43.388     7.5   0.892
   6    0.300     317   387.76   39.755     9.8   0.957
   7    0.350     429   449.06   48.554     9.2   0.985
   8    0.400     553   497.35   54.161     9.2   0.995
   9    0.450     709   534.49   69.574     7.7   0.999
  10    0.500     877   551.37   62.271     8.9   1.000
  11    0.550    1049   557.15   68.863     8.1   1.000
  12    0.600    1257   555.93   68.043     8.2   1.000
  13    0.650    1481   552.55  110.277     5.0   1.000
  14    0.700    1709   552.13  100.448     5.5   1.000

PEAK FLUX: radius     total     rms    snr
             0.55    557.15   68.86    8.1
PEAK SNR:  radius     total     rms   snr
             0.30    387.76   39.75    9.8

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

# can't see emission
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
          chans      = ('range=[11.0km/s,15.0km/s]'))
exportfits(imagename=lineimg+'_mom0',fitsimage=lineimg+'_mom0.fits')
os.system('rm -rf '+lineimg+'_mom0.image')
os.system('rm -rf '+lineimg+'_mom0')


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
Annulus = [2.0, 5.0] arcseconds

pixel scale     0.03 arcsec
image noise    33.025 mJy
beam :  PA   -60    0.34 x  0.25 arcsec   #pix/beam:   106
center pixel    150,  150
Annulus is between 2.0 and 5.0 arcsec
Subtracting offset of -0.346 mJy/beam from image
RMS in annulus is 32.8

#   i  radius    #pix    total      rms     snr   beam_frac
#      (asec)            (mJy)     (mJy)
   0    0.000       1    -3.54   35.230    -0.1   0.009
   1    0.050       9    -2.18   30.889    -0.1   0.082
   2    0.100      37     1.47   32.022     0.0   0.294
   3    0.150      81     3.47   35.553     0.1   0.528
   4    0.200     137     0.97   31.768     0.0   0.717
   5    0.250     221   -12.62   38.502    -0.3   0.864
   6    0.300     317   -38.50   41.373    -0.9   0.940
   7    0.350     429   -69.97   49.555    -1.4   0.977
   8    0.400     553   -99.71   54.246    -1.8   0.991
   9    0.450     709  -122.81   66.132    -1.9   0.997
  10    0.500     877  -133.35   52.120    -2.6   0.999
  11    0.550    1049  -130.89   64.313    -2.0   1.000
  12    0.600    1257  -122.71   73.255    -1.7   1.000
  13    0.650    1481  -105.97   72.076    -1.5   1.000
  14    0.700    1709   -90.58   84.153    -1.1   1.000

PEAK FLUX: radius     total     rms    snr
             0.15      3.47   35.55    0.1
PEAK SNR:  radius     total     rms   snr
             0.15      3.47   35.55    0.1
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


# can't see emission
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
          chans      = ('range=[11.0km/s,15.0km/s]'))
exportfits(imagename=lineimg+'_mom0',fitsimage=lineimg+'_mom0.fits')
os.system('rm -rf '+lineimg+'_mom0.image')
os.system('rm -rf '+lineimg+'_mom0')


### re-center image on source
os.system('rm -rf '+lineimg+'_mom0_crop.fits')
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
S76/S76_C18O_mom0_crop.fits
Annulus = [2.0, 5.0] arcseconds

pixel scale     0.03 arcsec
image noise    24.244 mJy
beam :  PA   -63    0.33 x  0.26 arcsec   #pix/beam:   107
center pixel    150,  150
Annulus is between 2.0 and 5.0 arcsec
Subtracting offset of 0.099 mJy/beam from image
RMS in annulus is 24

#   i  radius    #pix    total      rms     snr   beam_frac
#      (asec)            (mJy)     (mJy)
   0    0.000       1     9.13   28.209     0.3   0.009
   1    0.050       9     9.86   20.285     0.5   0.081
   2    0.100      37    12.58   22.881     0.6   0.291
   3    0.150      81    17.27   25.825     0.7   0.526
   4    0.200     137    24.12   20.681     1.2   0.715
   5    0.250     221    35.76   32.896     1.1   0.865
   6    0.300     317    47.09   35.634     1.3   0.941
   7    0.350     429    59.26   35.214     1.7   0.977
   8    0.400     553    68.51   37.665     1.8   0.992
   9    0.450     709    73.68   59.974     1.2   0.998
  10    0.500     877    77.20   58.904     1.3   0.999
  11    0.550    1049    80.15   63.631     1.3   1.000
  12    0.600    1257    86.83   52.638     1.6   1.000
  13    0.650    1481    95.03   66.791     1.4   1.000
  14    0.700    1709    98.28   58.721     1.7   1.000

PEAK FLUX: radius     total     rms    snr
             0.70     98.28   58.72    1.7
PEAK SNR:  radius     total     rms   snr
             0.40     68.51   37.67    1.8
'''
