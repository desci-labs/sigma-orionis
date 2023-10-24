#!/usr/bin/env


# ======================== Setup ===========================

field   = 86
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
### large offset in dec, right star?
xc,yc = 327,268
aper = 1.0


### check stats
img_max = imstat(imagename=contimg+'.image')['max'][0]
img_rms = imstat(imagename=contimg+'.image',region='annulus[['+str(xc)+'pix,'+str(yc)+'pix],['+str(80)+'pix,'+str(120)+'pix]]')['rms'][0]
bmaj = imhead(imagename=contimg+'.image', mode="get", hdkey="beammajor")
bmin = imhead(imagename=contimg+'.image', mode="get", hdkey="beamminor")
print '# Peak = {0:.2f} mJy, rms = {1:.2f} mJy, S/N = {2:.1f}'.format(1000*img_max, 1000*img_rms, img_max/img_rms)
print '# Beam = {0:.2f} x {1:.2f} arcsec'.format(bmaj.get('value'),bmin.get('value'))
# Peak = 2.04 mJy, rms = 0.17 mJy, S/N = 11.9
# Beam = 0.33 x 0.25 arcsec

### measure flux
img_rms = imstat(imagename=contimg+'.image',region='annulus[['+str(xc)+'pix,'+str(yc)+'pix],['+str(80)+'pix,'+str(120)+'pix]]')['rms'][0]
img_flx = imstat(imagename=contimg+'.image',region='circle[['+str(xc)+'pix,'+str(yc)+'pix],'+str(aper)+'arcsec]')['flux'][0]
print '# Flux = {0:.2f} mJy, rms = {1:.2f} mJy, S/N = {2:.1f}'.format(1000*img_flx, 1000*img_rms, img_flx/img_rms)
# Flux = 2.23 mJy, rms = 0.17 mJy, S/N = 13.0


# ======================== Measure flux with UVMODELFIT ==================


### calculate offset from phase center in arcsec
pixscale = 0.03             # must match 'cell'                 
dx = pixscale*(320.0-xc)    # offset to east (left)
dy = pixscale*(yc-320.0)    # offset to north (up)

  
#### measure flux as gaussian (marginally resolved)
uvmodelfit(vis       = contvis,
           comptype  = 'G',
           sourcepar = [img_flx,dx,dy,0.5,0.5,0.0],
           varypar   = [T,T,T,T,T,T],
           niter     = 10)

'''
reduced chi2=2.86
I = 0.00299639 +/- 0.000162879
x = -0.222304 +/- 0.00661228 arcsec
y = -1.56633 +/- 0.00686869 arcsec
a = 0.20698 +/- 0.023286 arcsec
r = 0.623986 +/- 0.171783
p = 0.662236 +/- 14.0593 deg

'''


#### measure flux as point source
uvmodelfit(vis       = contvis,
           comptype  = 'P',
           sourcepar = [img_flx,dx,dy],
           varypar   = [T,T,T],
           niter     = 10)
'''
I = 0.00227457 +/- 9.17222e-05
x = -0.216963 +/- 0.00601623 arcsec
y = -1.56219 +/- 0.00519628 arcsec

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
image noise     0.180 mJy
beam :  PA   -61    0.33 x  0.25 arcsec   #pix/beam:   104
center pixel    150,  150
Annulus is between 5.0 and 20.0 arcsec
Subtracting offset of 0.006 mJy/beam from image
RMS in annulus is 0.161

#   i  radius    #pix    total      rms     snr   beam_frac
#      (asec)            (mJy)     (mJy)
   0    0.000       1     2.04    0.156    13.1   0.010
   1    0.100      37     2.15    0.176    12.2   0.297
   2    0.200     137     2.43    0.182    13.4   0.722
   3    0.300     317     2.74    0.230    11.9   0.943
   4    0.400     553     2.73    0.283     9.7   0.992
   5    0.500     877     2.44    0.377     6.5   0.999
   6    0.600    1257     2.23    0.497     4.5   1.000
   7    0.700    1709     2.23    0.527     4.2   1.000
   8    0.800    2233     2.31    0.562     4.1   1.000
   9    0.900    2821     2.23    0.808     2.8   1.000

PEAK FLUX: radius     total     rms    snr
              0.3      2.74    0.23   11.9
PEAK SNR:  radius     total     rms   snr
              0.2      2.43    0.18   13.4
'''



