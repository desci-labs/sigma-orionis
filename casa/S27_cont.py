#!/usr/bin/env


# ======================== Setup ===========================

field   = 27 
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
      niter = 100,
      threshold = '0.0mJy',
      interactive = False,
      imagermode = 'csclean')
### couldn't see clear emission so just cleaned lightly

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
print imhead(imagename=contimg+'.image', mode="get", hdkey="object")
print '# Peak = {0:.2f} mJy, rms = {1:.2f} mJy, S/N = {2:.1f}'.format(1000*img_max, 1000*img_rms, img_max/img_rms)
print '# Beam = {0:.2f} x {1:.2f} arcsec'.format(bmaj.get('value'),bmin.get('value'))
# Peak = 0.60 mJy, rms = 0.15 mJy, S/N = 4.0
# Beam = 0.31 x 0.25 arcsec


### measure flux
img_rms = imstat(imagename=contimg+'.image',region='annulus[['+str(xc)+'pix,'+str(yc)+'pix],['+str(80)+'pix,'+str(120)+'pix]]')['rms'][0]
img_flx = imstat(imagename=contimg+'.image',region='circle[['+str(xc)+'pix,'+str(yc)+'pix],'+str(aper)+'arcsec]')['flux'][0]
print '# Flux = {0:.2f} mJy, rms = {1:.2f} mJy, S/N = {2:.1f}'.format(1000*img_flx, 1000*img_rms, img_flx/img_rms)
# Flux = 0.20 mJy, rms = 0.15 mJy, S/N = 1.4



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
I = 0.000522495 +/- 8.08082e-05
x = 0.0396013 +/- 0.0212211 arcsec
y = -0.0806748 +/- 0.0193361 arcsec

'''

### Note: still detection if using fixed position
###       or using [320,320] as source location
### Note: doesn't matter if you clean whole image
###       or select region that might be emission
### Note: confirmed by COG method

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
image noise     0.150 mJy
beam :  PA   -63    0.31 x  0.25 arcsec   #pix/beam:    96
center pixel    150,  150
Annulus is between 5.0 and 20.0 arcsec
Subtracting offset of 0.012 mJy/beam from image
RMS in annulus is 0.147

#   i  radius    #pix    total      rms     snr   beam_frac
#      (asec)            (mJy)     (mJy)
   0    0.000       1     0.45    0.168     2.7   0.010
   1    0.050       9     0.45    0.159     2.8   0.089
   2    0.100      37     0.46    0.160     2.9   0.319
   3    0.150      81     0.46    0.133     3.5   0.565 <-- good
   4    0.200     137     0.44    0.193     2.3   0.755
   5    0.250     221     0.40    0.179     2.2   0.894
   6    0.300     317     0.36    0.224     1.6   0.959
   7    0.350     429     0.33    0.230     1.4   0.986
   8    0.400     553     0.33    0.243     1.3   0.996
   9    0.450     709     0.32    0.304     1.0   0.999
  10    0.500     877     0.28    0.373     0.8   1.000
  11    0.550    1049     0.23    0.332     0.7   1.000
  12    0.600    1257     0.20    0.377     0.5   1.000
  13    0.650    1481     0.20    0.307     0.7   1.000
  14    0.700    1709     0.21    0.358     0.6   1.000
  15    0.750    1961     0.14    0.363     0.4   1.000
  16    0.800    2233     0.02    0.409     0.1   1.000
  17    0.850    2537    -0.12    0.391    -0.3   1.000
  18    0.900    2821    -0.20    0.612    -0.3   1.000
  19    0.950    3149    -0.20    0.357    -0.6   1.000

PEAK FLUX: radius     total     rms    snr
             0.15      0.46    0.13    3.5
PEAK SNR:  radius     total     rms   snr
             0.15      0.46    0.13    3.5

'''



