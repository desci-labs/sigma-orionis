#!/usr/bin/env


# ======================== Setup ===========================

field   = 79
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
# not resolved



                        
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
xc,yc = 319,319
aper = 1.0


### check stats
img_max = imstat(imagename=contimg+'.image')['max'][0]
img_rms = imstat(imagename=contimg+'.image',region='annulus[['+str(xc)+'pix,'+str(yc)+'pix],['+str(80)+'pix,'+str(120)+'pix]]')['rms'][0]
bmaj = imhead(imagename=contimg+'.image', mode="get", hdkey="beammajor")
bmin = imhead(imagename=contimg+'.image', mode="get", hdkey="beamminor")
print '# Peak = {0:.2f} mJy, rms = {1:.2f} mJy, S/N = {2:.1f}'.format(1000*img_max, 1000*img_rms, img_max/img_rms)
print '# Beam = {0:.2f} x {1:.2f} arcsec'.format(bmaj.get('value'),bmin.get('value'))
# Peak = 5.32 mJy, rms = 0.18 mJy, S/N = 29.7
# Beam = 0.32 x 0.25 arcsec


### measure flux
img_rms = imstat(imagename=contimg+'.image',region='annulus[['+str(xc)+'pix,'+str(yc)+'pix],['+str(80)+'pix,'+str(120)+'pix]]')['rms'][0]
img_flx = imstat(imagename=contimg+'.image',region='circle[['+str(xc)+'pix,'+str(yc)+'pix],'+str(aper)+'arcsec]')['flux'][0]
print '# Flux = {0:.2f} mJy, rms = {1:.2f} mJy, S/N = {2:.1f}'.format(1000*img_flx, 1000*img_rms, img_flx/img_rms)
# Flux = 5.69 mJy, rms = 0.18 mJy, S/N = 31.7


# ======================== Measure flux with UVMODELFIT ==================


### calculate offset from phase center in arcsec
pixscale = 0.03             # must match 'cell'                 
dx = pixscale*(320.0-xc)    # offset to east (left)
dy = pixscale*(yc-320.0)    # offset to north (up)

  
#### measure flux as gaussian
uvmodelfit(vis       = contvis,
           comptype  = 'G',
           sourcepar = [img_flx,dx,dy,0.5,0.5,0.0],
           varypar   = [T,T,T,T,T,T],
           niter     = 10)

'''
reduced chi2=2.891
I = 0.00631415 +/- 0.000141095
x = 0.023114 +/- 0.00248323 arcsec
y = -0.0420716 +/- 0.00220721 arcsec
a = 0.0978194 +/- 0.0130321 arcsec <-- bad
r = 1 +/- 0.198694
p = -8.68588 +/- 57.2958 deg <--bad

'''


#### measure flux as point source
uvmodelfit(vis       = contvis,
           comptype  = 'P',
           sourcepar = [img_flx,dx,dy],
           varypar   = [T,T,T],
           niter     = 10)

'''
I = 0.00566196 +/- 9.15907e-05
x = 0.0235379 +/- 0.00233097 arcsec
y = -0.0424842 +/- 0.00205326 arcsec


'''



# ================== MEASURE FLUX WITH CURVE OF GROWTH ==================


### re-center image on source 
os.system('rm -f '+contimg+'_crop.fits')
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
image noise     0.223 mJy
beam :  PA   -61    0.32 x  0.25 arcsec   #pix/beam:   101
center pixel    150,  150
Annulus is between 5.0 and 20.0 arcsec
Subtracting offset of -0.006 mJy/beam from image
RMS in annulus is 0.176

#   i  radius    #pix    total      rms     snr   beam_frac
#      (asec)            (mJy)     (mJy)
   0    0.000       1     5.33    0.183    29.0   0.010
   1    0.100      37     5.45    0.162    33.7   0.305
   2    0.200     137     5.78    0.172    33.5   0.735
   3    0.300     317     6.25    0.188    33.2   0.949
   4    0.400     553     6.64    0.323    20.5   0.994
   5    0.500     877     6.88    0.406    16.9   1.000
   6    0.600    1257     6.92    0.361    19.2   1.000
   7    0.700    1709     6.72    0.478    14.1   1.000
   8    0.800    2233     6.39    0.640    10.0   1.000
   9    0.900    2821     6.08    0.801     7.6   1.000
  10    1.000    3505     5.91    0.490    12.1   1.000
  11    1.100    4205     6.03    0.769     7.8   1.000

PEAK FLUX: radius     total     rms    snr
              0.6      6.92    0.36   19.2
PEAK SNR:  radius     total     rms   snr
              0.1      5.45    0.16   33.7
'''



