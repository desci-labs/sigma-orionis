#!/usr/bin/env


# ======================== Setup ===========================

field   = 29 
visfile = '../calibrated_final.ms'
contspw = '0,3,5,8,10,13,15,18'
contwid = [128,1920,128,1920,128,1920,128,1920]

cell   = '0.03arcsec'
imsize = [640,640]
robust = 0.5


# ======================= Split Off Continuum ========================

### continuum visibilities for this field
contvis = 'S'+str(field)+'_cont.vis'


### average continuum channels for this field
split(vis = visfile,
      outputvis = contvis,
      field = field,
      datacolumn = 'data',
      spw = contspw,
      width = contwid)


### plot uv-distance vs. amplitude
plotms(vis = contvis,
       xaxis = 'uvdist',yaxis='amp',
       coloraxis = 'spw')
# maybe resolved?



                        
# ================== IMAGE CONTINUUM ==================

### continuum image for this field
contimg = 'S'+str(field)+'_cont'


### clear previous stuff
clearcal(vis=contvis)
delmod(vis=contvis)
for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage']: rmtables(contimg+ext)


### perform clean
clean(vis = contvis,
      imagename = contimg,
      field = '0',
      mode = 'mfs',
      psfmode = 'clark',
      imsize = imsize,
      cell = cell,
      weighting = 'briggs',
      robust = robust,
      niter = 1000,
      threshold = '0.0mJy',
      interactive = True,
      imagermode = 'csclean')


### export to fits
os.system('rm -f '+contimg+'.fits')
exportfits(imagename=contimg+'.image', fitsimage=contimg+'.fits')



# ================== MEASURE FLUX WITHIN APERTURE ==================


### check image to get center
imview(raster   = [{'file':contimg+'.fits'}],
        contour = [{'file':contimg+'.fits'}])


### set aperture and source center
xc,yc = 319,318
aper = 1.0


### check stats
img_max = imstat(imagename=contimg+'.image')['max'][0]
img_rms = imstat(imagename=contimg+'.image',region='annulus[['+str(xc)+'pix,'+str(yc)+'pix],['+str(80)+'pix,'+str(120)+'pix]]')['rms'][0]
bmaj = imhead(imagename=contimg+'.image', mode="get", hdkey="beammajor")
bmin = imhead(imagename=contimg+'.image', mode="get", hdkey="beamminor")
print '# Peak = {0:.2f} mJy, rms = {1:.2f} mJy, S/N = {2:.1f}'.format(1000*img_max, 1000*img_rms, img_max/img_rms)
print '# Beam = {0:.2f} x {1:.2f} arcsec'.format(bmaj.get('value'),bmin.get('value'))
# Peak = 5.73 mJy, rms = 0.15 mJy, S/N = 37.1
# Beam = 0.31 x 0.25 arcsec


### measure flux
img_rms = imstat(imagename=contimg+'.image',region='annulus[['+str(xc)+'pix,'+str(yc)+'pix],['+str(80)+'pix,'+str(120)+'pix]]')['rms'][0]
img_flx = imstat(imagename=contimg+'.image',region='circle[['+str(xc)+'pix,'+str(yc)+'pix],'+str(aper)+'arcsec]')['flux'][0]
print '# Flux = {0:.2f} mJy, rms = {1:.2f} mJy, S/N = {2:.1f}'.format(1000*img_flx, 1000*img_rms, img_flx/img_rms)
# Flux = 11.97 mJy, rms = 0.15 mJy, S/N = 77.7



# ======================== Measure flux with UVMODELFIT ==================


### calculate offset from phase center in arcsec
pixscale = 0.03             # must match 'cell'                 
dx = pixscale*(320.0-xc)    # offset to east (left)
dy = pixscale*(yc-320.0)    # offset to north (up)

  
#### measure flux as gaussian (resolved wtih visbin)
uvmodelfit(vis       = contvis,
           comptype  = 'G',
           sourcepar = [img_flx,dx,dy,0.5,0.5,0.0],
           varypar   = [T,T,T,T,T,T],
           niter     = 10)

'''
reduced chi2=2.889
I = 0.0106878 +/- 0.000171925
x = 0.0189835 +/- 0.0020049 arcsec
y = -0.0661702 +/- 0.0025494 arcsec
a = 0.31584 +/- 0.00708159 arcsec
r = 0.503478 +/- 0.0282478
p = 13.8547 +/- 1.87673 deg

'''



# ================== MEASURE FLUX WITH CURVE OF GROWTH ==================


### remove previous image
os.system('rm -f '+contimg+'_crop.fits')

### re-center image on source 
boxwidth = 150.
box = rg.box([xc-boxwidth,yc-boxwidth],[xc+boxwidth,yc+boxwidth])
ia.fromimage(outfile=contimg+'_crop.image',infile=contimg+'.fits',region=box )
ia.close() 
exportfits(imagename = contimg+'_crop.image',fitsimage = contimg+'_crop.fits')
os.system('rm -rf '+contimg+'_crop.image')


### check image at center
imview(raster   = [{'file':contimg+'_crop.fits'}])


### use measure.py to get COG flux
'''
pixel scale     0.03 arcsec
image noise     0.238 mJy
beam :  PA   -62    0.31 x  0.25 arcsec   #pix/beam:    97
center pixel    150,  150
Annulus is between 5.0 and 20.0 arcsec
Subtracting offset of -0.021 mJy/beam from image
RMS in annulus is 0.15

#   i  radius    #pix    total      rms     snr   beam_frac
#      (asec)            (mJy)     (mJy)
   0    0.000       1     5.75    0.166    34.5   0.010
   1    0.050       9     5.87    0.163    36.0   0.089
   2    0.100      37     6.28    0.171    36.8   0.316
   3    0.150      81     6.89    0.144    48.0   0.561
   4    0.200     137     7.58    0.196    38.7   0.751
   5    0.250     221     8.44    0.187    45.1   0.891
   6    0.300     317     9.14    0.168    54.4   0.957
   7    0.350     429     9.70    0.207    46.8   0.985
   8    0.400     553    10.10    0.254    39.8   0.995 <-- looks good in image, matches uvmodelfit
   9    0.450     709    10.43    0.270    38.6   0.999
  10    0.500     877    10.69    0.323    33.1   1.000
  11    0.550    1049    10.93    0.280    39.1   1.000
  12    0.600    1257    11.18    0.290    38.6   1.000
  13    0.650    1481    11.44    0.364    31.4   1.000
  14    0.700    1709    11.68    0.321    36.3   1.000
  15    0.750    1961    11.89    0.301    39.4   1.000
  16    0.800    2233    12.08    0.406    29.8   1.000
  17    0.850    2537    12.26    0.440    27.8   1.000
  18    0.900    2821    12.40    0.717    17.3   1.000
  19    0.950    3149    12.57    0.626    20.1   1.000

PEAK FLUX: radius     total     rms    snr
             0.95     12.57    0.63   20.1
PEAK SNR:  radius     total     rms   snr
             0.30      9.14    0.17   54.4
             
COG not working; flux keeps increasing with radius...

'''



