#!/usr/bin/env


# ================================== Setup =====================================

field   = 88 
visfile = '../calibrated_final.ms'
linevis = 'S'+str(field)+'.vis.contsub' 

fitspw  = '0,3,5,8,10,13,15,18'
linespw = '1,2,4,6,7,9,11,12,14,16,17,19'

cell   = '0.03arcsec'
imsize = [640,640]
robust = 0.5

xc,yc = 290,314


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


# emission clearly seen 12-13 km/s (confirmed in spectrum)
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
          chans      = ('range=[12km/s,13km/s]'))
exportfits(imagename=lineimg+'_mom0',fitsimage=lineimg+'_mom0.fits')
os.system('rm -rf '+lineimg+'_mom0.image')
imview(raster   = [{'file':lineimg+'_mom0.fits'}],
        contour = [{'file':'S'+str(field)+'_cont.fits'}])

### delete stuff
for ext in ['.flux','.image','.mask','.model','.psf','.residual']: rmtables(lineimg+ext)
os.system('rm -rf cgrid_ft.im *.last *.log')


### make first-moment map
os.system('rm -rf '+lineimg+'_mom1*')
immoments(imagename  = lineimg+'.fits',
          outfile    = lineimg+'_mom1',
          moments    = [1],
          excludepix = [-100,0.005],
          chans      = ('range=[12km/s,13km/s]'))
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
S88/S88_12CO_mom0_crop.fits
Annulus = [2.0, 5.0] arcseconds

pixel scale     0.03 arcsec
image noise    20.525 mJy
beam :  PA   -61    0.32 x  0.23 arcsec   #pix/beam:    94
center pixel    150,  150
Annulus is between 2.0 and 5.0 arcsec
Subtracting offset of -0.504 mJy/beam from image
RMS in annulus is 19.4

#   i  radius    #pix    total      rms     snr   beam_frac
#      (asec)            (mJy)     (mJy)
   0    0.000       1   102.92   20.918     4.9   0.011
   1    0.050       9   106.48   18.688     5.7   0.092
   2    0.100      37   119.38   18.548     6.4   0.324
   3    0.150      81   142.06   21.496     6.6   0.570
   4    0.200     137   174.20   20.806     8.4   0.758
   5    0.250     221   229.47   26.630     8.6   0.893
   6    0.300     317   295.33   22.844    12.9   0.957
   7    0.350     429   373.00   32.023    11.6   0.985
   8    0.400     553   452.66   26.591    17.0   0.995
   9    0.450     709   534.05   35.671    15.0   0.999
  10    0.500     877   606.66   44.140    13.7   1.000
  11    0.550    1049   666.24   46.642    14.3   1.000
  12    0.600    1257   717.58   56.139    12.8   1.000
  13    0.650    1481   761.90   48.586    15.7   1.000
  14    0.700    1709   793.68   48.708    16.3   1.000
  15    0.750    1961   821.16   63.819    12.9   1.000
  16    0.800    2233   840.18   68.369    12.3   1.000
  17    0.850    2537   854.47   74.111    11.5   1.000
  18    0.900    2821   860.76   87.780     9.8   1.000 <-- levels off
  19    0.950    3149   862.37  104.811     8.2   1.000
  20    1.000    3505   865.08   58.363    14.8   1.000
  21    1.050    3853   870.53  114.077     7.6   1.000
  22    1.100    4205   879.31  131.622     6.7   1.000

PEAK FLUX: radius     total     rms    snr
             1.10    879.31  131.62    6.7
PEAK SNR:  radius     total     rms   snr
             0.40    452.66   26.59   17.0
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

# emission seen at 12 km/s
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
          chans      = ('range=[12.0km/s,13.0km/s]'))
exportfits(imagename=lineimg+'_mom0',fitsimage=lineimg+'_mom0.fits')
os.system('rm -rf '+lineimg+'_mom0.image')


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

S88/S88_13CO_mom0_crop.fits
Annulus = [2.0, 5.0] arcseconds

pixel scale     0.03 arcsec
image noise    20.434 mJy
beam :  PA   -61    0.34 x  0.25 arcsec   #pix/beam:   106
center pixel    150,  150
Annulus is between 2.0 and 5.0 arcsec
Subtracting offset of -0.864 mJy/beam from image
RMS in annulus is 20.3

#   i  radius    #pix    total      rms     snr   beam_frac
#      (asec)            (mJy)     (mJy)
   0    0.000       1    38.74   19.678     2.0   0.009
   1    0.050       9    39.26   19.025     2.1   0.081
   2    0.100      37    41.25   17.817     2.3   0.293
   3    0.150      81    45.54   22.027     2.1   0.527
   4    0.200     137    52.29   21.117     2.5   0.715
   5    0.250     221    65.88   27.629     2.4   0.863
   6    0.300     317    86.85   26.711     3.3   0.939
   7    0.350     429   112.81   28.902     3.9   0.976
   8    0.400     553   140.94   43.757     3.2   0.991
   9    0.450     709   170.10   33.042     5.1   0.997
  10    0.500     877   191.92   49.924     3.8   0.999
  11    0.550    1049   206.25   44.086     4.7   1.000
  12    0.600    1257   217.71   51.274     4.2   1.000
  13    0.650    1481   232.41   48.257     4.8   1.000
  14    0.700    1709   253.32   53.937     4.7   1.000
  15    0.750    1961   278.31   76.214     3.7   1.000
  16    0.800    2233   304.44   83.549     3.6   1.000
  17    0.850    2537   320.82   59.317     5.4   1.000
  18    0.900    2821   325.77   67.788     4.8   1.000 <-- levels off
  19    0.950    3149   322.68   70.332     4.6   1.000
  20    1.000    3505   320.43   71.576     4.5   1.000
  21    1.050    3853   325.13  131.476     2.5   1.000
  22    1.100    4205   330.56   80.518     4.1   1.000

PEAK FLUX: radius     total     rms    snr
             1.10    330.56   80.52    4.1
PEAK SNR:  radius     total     rms   snr
             0.85    320.82   59.32    5.4


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


# no clear emission seen
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
          chans      = ('range=[12km/s,13km/s]'))
exportfits(imagename=lineimg+'_mom0',fitsimage=lineimg+'_mom0.fits')
os.system('rm -rf '+lineimg+'_mom0.image')


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
S88/S88_C18O_mom0_crop.fits
Annulus = [2.0, 5.0] arcseconds

pixel scale     0.03 arcsec
image noise    16.040 mJy
beam :  PA   -63    0.34 x  0.25 arcsec   #pix/beam:   108
center pixel    150,  150
Annulus is between 2.0 and 5.0 arcsec
Subtracting offset of 0.393 mJy/beam from image
RMS in annulus is 16.2

#   i  radius    #pix    total      rms     snr   beam_frac
#      (asec)            (mJy)     (mJy)
   0    0.000       1     6.15   17.471     0.4   0.009
   1    0.050       9     6.42   14.543     0.4   0.080
   2    0.100      37     7.54   14.696     0.5   0.289
   3    0.150      81     9.42   16.389     0.6   0.522
   4    0.200     137    12.83   15.149     0.8   0.711
   5    0.250     221    19.88   18.137     1.1   0.860
   6    0.300     317    28.69   22.091     1.3   0.938
   7    0.350     429    40.85   26.952     1.5   0.975
   8    0.400     553    52.40   26.183     2.0   0.991
   9    0.450     709    62.65   30.293     2.1   0.997
  10    0.500     877    66.41   31.518     2.1   0.999
  11    0.550    1049    65.50   28.254     2.3   1.000
  12    0.600    1257    61.19   35.409     1.7   1.000
  13    0.650    1481    57.33   39.112     1.5   1.000
  14    0.700    1709    54.35   52.175     1.0   1.000
  15    0.750    1961    49.96   26.166     1.9   1.000
  16    0.800    2233    45.07   38.041     1.2   1.000
  17    0.850    2537    39.27   42.045     0.9   1.000
  18    0.900    2821    35.17   37.008     1.0   1.000
  19    0.950    3149    32.78   41.513     0.8   1.000
  20    1.000    3505    32.98   30.837     1.1   1.000
  21    1.050    3853    39.56   45.144     0.9   1.000
  22    1.100    4205    46.92   69.814     0.7   1.000

PEAK FLUX: radius     total     rms    snr
             0.50     66.41   31.52    2.1
PEAK SNR:  radius     total     rms   snr
             0.55     65.50   28.25    2.3
'''
