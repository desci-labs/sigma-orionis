#!/usr/bin/env


# ======================== Setup ===========================

field   = 75
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
xc,yc = 321,320
aper = 1.0


### check stats
img_max = imstat(imagename=contimg+'.image')['max'][0]
img_rms = imstat(imagename=contimg+'.image',region='annulus[['+str(xc)+'pix,'+str(yc)+'pix],['+str(80)+'pix,'+str(120)+'pix]]')['rms'][0]
bmaj = imhead(imagename=contimg+'.image', mode="get", hdkey="beammajor")
bmin = imhead(imagename=contimg+'.image', mode="get", hdkey="beamminor")
print '# Peak = {0:.2f} mJy, rms = {1:.2f} mJy, S/N = {2:.1f}'.format(1000*img_max, 1000*img_rms, img_max/img_rms)
print '# Beam = {0:.2f} x {1:.2f} arcsec'.format(bmaj.get('value'),bmin.get('value'))
# Peak = 5.31 mJy, rms = 0.18 mJy, S/N = 29.4
# Beam = 0.32 x 0.25 arcsec


### measure flux
img_rms = imstat(imagename=contimg+'.image',region='annulus[['+str(xc)+'pix,'+str(yc)+'pix],['+str(80)+'pix,'+str(120)+'pix]]')['rms'][0]
img_flx = imstat(imagename=contimg+'.image',region='circle[['+str(xc)+'pix,'+str(yc)+'pix],'+str(aper)+'arcsec]')['flux'][0]
print '# Flux = {0:.2f} mJy, rms = {1:.2f} mJy, S/N = {2:.1f}'.format(1000*img_flx, 1000*img_rms, img_flx/img_rms)
# Flux = 9.15 mJy, rms = 0.18 mJy, S/N = 50.7



# ======================== Measure flux with UVMODELFIT ==================


### calculate offset from phase center in arcsec
pixscale = 0.03             # must match 'cell'                 
dx = pixscale*(320.0-xc)    # offset to east (left)
dy = pixscale*(yc-320.0)    # offset to north (up)

  
#### measure flux as gaussian (clearly resolved with visbin)
uvmodelfit(vis       = contvis,
           comptype  = 'G',
           sourcepar = [img_flx,dx,dy,0.5,0.5,0.0],
           varypar   = [T,T,T,T,T,T],
           niter     = 10)

'''
reduced chi2=2.949
I = 0.00856843 +/- 0.000174011
x = -0.0314295 +/- 0.00284998 arcsec
y = -0.00187217 +/- 0.00308034 arcsec
a = 0.312627 +/- 0.0103293 arcsec
r = 0.407306 +/- 0.0336096
p = -32.6559 +/- 2.02278 deg

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
image noise     0.240 mJy
beam :  PA   -60    0.32 x  0.25 arcsec   #pix/beam:   103
center pixel    150,  150
Annulus is between 5.0 and 20.0 arcsec
Subtracting offset of -0.012 mJy/beam from image
RMS in annulus is 0.162

#   i  radius    #pix    total      rms     snr   beam_frac
#      (asec)            (mJy)     (mJy)
   0    0.000       1     5.32    0.176    30.2   0.010
   1    0.050       9     5.40    0.185    29.1   0.084
   2    0.100      37     5.64    0.187    30.1   0.302
   3    0.150      81     6.00    0.180    33.3   0.541
   4    0.200     137     6.42    0.203    31.7   0.730
   5    0.250     221     6.98    0.219    31.9   0.875
   6    0.300     317     7.51    0.232    32.4   0.947
   7    0.350     429     7.98    0.294    27.1   0.980
   8    0.400     553     8.35    0.254    32.9   0.993
   9    0.450     709     8.70    0.335    26.0   0.998
  10    0.500     877     8.93    0.492    18.1   1.000 <-- GOOD
  11    0.550    1049     9.07    0.493    18.4   1.000
  12    0.600    1257     9.18    0.530    17.3   1.000
  13    0.650    1481     9.26    0.590    15.7   1.000
  14    0.700    1709     9.33    0.682    13.7   1.000

PEAK FLUX: radius     total     rms    snr
              0.7      9.33    0.68   13.7
PEAK SNR:  radius     total     rms   snr
              0.2      6.00    0.18   33.3
'''



