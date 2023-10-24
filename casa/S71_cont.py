#!/usr/bin/env


# ======================== Setup ===========================

field   = 71
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
xc,yc = 316,316
aper = 1.0


### check stats
img_max = imstat(imagename=contimg+'.image')['max'][0]
img_rms = imstat(imagename=contimg+'.image',region='annulus[['+str(xc)+'pix,'+str(yc)+'pix],['+str(80)+'pix,'+str(120)+'pix]]')['rms'][0]
bmaj = imhead(imagename=contimg+'.image', mode="get", hdkey="beammajor")
bmin = imhead(imagename=contimg+'.image', mode="get", hdkey="beamminor")
print '# Peak = {0:.2f} mJy, rms = {1:.2f} mJy, S/N = {2:.1f}'.format(1000*img_max, 1000*img_rms, img_max/img_rms)
print '# Beam = {0:.2f} x {1:.2f} arcsec'.format(bmaj.get('value'),bmin.get('value'))
# Peak = 5.83 mJy, rms = 0.17 mJy, S/N = 35.0
# Beam = 0.32 x 0.25 arcsec


### measure flux
img_rms = imstat(imagename=contimg+'.image',region='annulus[['+str(xc)+'pix,'+str(yc)+'pix],['+str(80)+'pix,'+str(120)+'pix]]')['rms'][0]
img_flx = imstat(imagename=contimg+'.image',region='circle[['+str(xc)+'pix,'+str(yc)+'pix],'+str(aper)+'arcsec]')['flux'][0]
print '# Flux = {0:.2f} mJy, rms = {1:.2f} mJy, S/N = {2:.1f}'.format(1000*img_flx, 1000*img_rms, img_flx/img_rms)
# Flux = 7.59 mJy, rms = 0.17 mJy, S/N = 45.6


# ======================== Measure flux with UVMODELFIT ==================


### calculate offset from phase center in arcsec
pixscale = 0.03             # must match 'cell'                 
dx = pixscale*(320.0-xc)    # offset to east (left)
dy = pixscale*(yc-320.0)    # offset to north (up)

  
#### measure flux as gaussian (unresolved with visbin)
uvmodelfit(vis       = contvis,
           comptype  = 'G',
           sourcepar = [img_flx,dx,dy,0.5,0.5,0.0],
           varypar   = [T,T,T,T,T,T],
           niter     = 10)

'''
reduced chi2=2.86
I = 0.00660284 +/- 0.000133871
x = 0.131665 +/- 0.00221619 arcsec
y = -0.130914 +/- 0.00198838 arcsec
a = 0.088763 +/- 0.0122457 arcsec <-- bad
r = 1 +/- 0.213591
p = -2.69873 +/- 57.2958 deg <-- bad

'''

#### measure flux as point source
uvmodelfit(vis       = contvis,
           comptype  = 'P',
           sourcepar = [img_flx,dx,dy],
           varypar   = [T,T,T],
           niter     = 10)

'''
reduced chi2=2.86
I = 0.00606923 +/- 8.83885e-05
x = 0.131167 +/- 0.00208574 arcsec
y = -0.129343 +/- 0.00185912 arcsec
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
image noise     0.225 mJy
beam :  PA   -61    0.32 x  0.25 arcsec   #pix/beam:   101
center pixel    150,  150
Annulus is between 5.0 and 20.0 arcsec
Subtracting offset of -0.007 mJy/beam from image
RMS in annulus is 0.164

#   i  radius    #pix    total      rms     snr   beam_frac
#      (asec)            (mJy)     (mJy)
   0    0.000       1     5.84    0.196    29.8   0.010
   1    0.050       9     5.86    0.153    38.3   0.085
   2    0.100      37     5.93    0.205    28.9   0.305
   3    0.150      81     6.05    0.147    41.0   0.545
   4    0.200     137     6.21    0.200    31.0   0.734 <-- best match
   5    0.250     221     6.47    0.192    33.6   0.879
   6    0.300     317     6.75    0.217    31.2   0.949
   7    0.350     429     7.03    0.285    24.7   0.982
   8    0.400     553     7.25    0.367    19.7   0.994
   9    0.450     709     7.43    0.346    21.5   0.998
  10    0.500     877     7.51    0.371    20.2   1.000
  11    0.550    1049     7.55    0.408    18.5   1.000
  12    0.600    1257     7.60    0.444    17.1   1.000
  13    0.650    1481     7.60    0.382    19.9   1.000
  14    0.700    1709     7.54    0.447    16.9   1.000

PEAK FLUX: radius     total     rms    snr
              0.6      7.60    0.44   17.1
PEAK SNR:  radius     total     rms   snr
              0.2      6.05    0.15   41.0
'''



