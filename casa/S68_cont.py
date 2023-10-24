#!/usr/bin/env


# ======================== Setup ===========================

field   = 68 
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




                        
# ================== IMAGE CONTINUUM ==================

### continuum image for this field
contimg = 'S'+str(field)+'_cont'


### clear previous stuff
clearcal(vis=contvis)
delmod(vis=contvis)
for ext in ['.flux','.image','.mask','.model','.pbcor',
            '.psf','.residual','.flux.pbcoverage']: rmtables(contimg+ext)


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
      niter = 500,
      threshold = '0.0mJy',
      interactive = False,
      imagermode = 'csclean')


### export to fits
os.system('rm -f '+contimg+'.fits')
exportfits(imagename=contimg+'.image', fitsimage=contimg+'.fits')



# ================== MEASURE FLUX WITHIN APERTURE ==================


### check image to get center
imview(raster   = [{'file':contimg+'.fits'}],
        contour = [{'file':contimg+'.fits'}])


### set aperture and source center
xc,yc = 320,320
aper = 1.0


### check stats
img_max = imstat(imagename=contimg+'.image')['max'][0]
img_rms = imstat(imagename=contimg+'.image',region='annulus[['+str(xc)+'pix,'+str(yc)+'pix],['+str(80)+'pix,'+str(120)+'pix]]')['rms'][0]
bmaj = imhead(imagename=contimg+'.image', mode="get", hdkey="beammajor")
bmin = imhead(imagename=contimg+'.image', mode="get", hdkey="beamminor")
print imhead(imagename=contimg+'.image', mode="get", hdkey="object")
print '# Peak = {0:.2f} mJy, rms = {1:.2f} mJy, S/N = {2:.1f}'.format(1000*img_max, 1000*img_rms, img_max/img_rms)
print '# Beam = {0:.2f} x {1:.2f} arcsec'.format(bmaj.get('value'),bmin.get('value'))
# Peak = 0.58 mJy, rms = 0.14 mJy, S/N = 4.0
# Beam = 0.32 x 0.25 arcsec


### measure flux
img_rms = imstat(imagename=contimg+'.image',region='annulus[['+str(xc)+'pix,'+str(yc)+'pix],['+str(80)+'pix,'+str(120)+'pix]]')['rms'][0]
img_flx = imstat(imagename=contimg+'.image',region='circle[['+str(xc)+'pix,'+str(yc)+'pix],'+str(aper)+'arcsec]')['flux'][0]
print '# Flux = {0:.2f} mJy, rms = {1:.2f} mJy, S/N = {2:.1f}'.format(1000*img_flx, 1000*img_rms, img_flx/img_rms)
# Flux = -0.67 mJy, rms = 0.14 mJy, S/N = -4.7


# ======================== Measure flux with UVMODELFIT ==================


### calculate offset from phase center in arcsec
pixscale = 0.03             # must match 'cell'                 
dx = pixscale*(320.0-xc)    # offset to east (left)
dy = pixscale*(yc-320.0)    # offset to north (up)

#### measure flux as point source
uvmodelfit(vis       = contvis,
           comptype  = 'P',
           sourcepar = [img_flx,dx,dy],
           varypar   = [T,F,F],
           niter     = 10)
'''
I = 8.64826e-05 +/- 8.69971e-05
x = 0 +/- 0 arcsec
y = 0 +/- 0 arcsec
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
image noise     0.146 mJy
beam :  PA   -62    0.32 x  0.25 arcsec   #pix/beam:   101
center pixel    150,  150
Annulus is between 5.0 and 20.0 arcsec
Subtracting offset of -0.005 mJy/beam from image
RMS in annulus is 0.158

#   i  radius    #pix    total      rms     snr   beam_frac
#      (asec)            (mJy)     (mJy)
   0    0.000       1    -0.02    0.160    -0.1   0.010
   1    0.050       9    -0.03    0.129    -0.2   0.085
   2    0.100      37    -0.07    0.161    -0.4   0.305
   3    0.150      81    -0.12    0.147    -0.8   0.545
   4    0.200     137    -0.18    0.138    -1.3   0.735
   5    0.250     221    -0.25    0.182    -1.4   0.879
   6    0.300     317    -0.29    0.210    -1.4   0.950
   7    0.350     429    -0.29    0.207    -1.4   0.982
   8    0.400     553    -0.26    0.219    -1.2   0.994
   9    0.450     709    -0.23    0.279    -0.8   0.998
  10    0.500     877    -0.22    0.333    -0.7   1.000
  11    0.550    1049    -0.26    0.444    -0.6   1.000
  12    0.600    1257    -0.37    0.405    -0.9   1.000
  13    0.650    1481    -0.49    0.444    -1.1   1.000
  14    0.700    1709    -0.56    0.407    -1.4   1.000
  15    0.750    1961    -0.56    0.481    -1.2   1.000
  16    0.800    2233    -0.49    0.391    -1.3   1.000
  17    0.850    2537    -0.42    0.445    -1.0   1.000
  18    0.900    2821    -0.40    0.428    -0.9   1.000
  19    0.950    3149    -0.43    0.501    -0.9   1.000

PEAK FLUX: radius     total     rms    snr
              0.0     -0.02    0.16   -0.1
PEAK SNR:  radius     total     rms   snr
              0.0     -0.02    0.16   -0.1

'''



