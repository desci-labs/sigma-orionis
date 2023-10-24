#!/usr/bin/env


# ======================== Setup ===========================

field   = 4 
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
xc,yc = 317,321
aper = 1.0


### check stats
img_max = imstat(imagename=contimg+'.image')['max'][0]
img_rms = imstat(imagename=contimg+'.image',region='annulus[['+str(xc)+'pix,'+str(yc)+'pix],['+str(80)+'pix,'+str(120)+'pix]]')['rms'][0]
bmaj = imhead(imagename=contimg+'.image', mode="get", hdkey="beammajor")
bmin = imhead(imagename=contimg+'.image', mode="get", hdkey="beamminor")
print imhead(imagename=contimg+'.image', mode="get", hdkey="object")
print '# Peak = {0:.2f} mJy, rms = {1:.2f} mJy, S/N = {2:.1f}'.format(1000*img_max, 1000*img_rms, img_max/img_rms)
print '# Beam = {0:.2f} x {1:.2f} arcsec'.format(bmaj.get('value'),bmin.get('value'))
# Peak = 0.64 mJy, rms = 0.14 mJy, S/N = 4.6
# Beam = 0.30 x 0.25 arcsec


### measure flux
img_rms = imstat(imagename=contimg+'.image',region='annulus[['+str(xc)+'pix,'+str(yc)+'pix],['+str(80)+'pix,'+str(120)+'pix]]')['rms'][0]
img_flx = imstat(imagename=contimg+'.image',region='circle[['+str(xc)+'pix,'+str(yc)+'pix],'+str(aper)+'arcsec]')['flux'][0]
print '# Flux = {0:.2f} mJy, rms = {1:.2f} mJy, S/N = {2:.1f}'.format(1000*img_flx, 1000*img_rms, img_flx/img_rms)
# Flux = 0.25 mJy, rms = 0.14 mJy, S/N = 1.8


# ======================== Measure flux with UVMODELFIT ==================


### calculate offset from phase center in arcsec
pixscale = 0.03             # must match 'cell'                 
dx = pixscale*(320.0-xc)    # offset to east (left)
dy = pixscale*(yc-320.0)    # offset to north (up)


#### measure flux as point source
uvmodelfit(vis       = contvis,
           comptype  = 'P',
           sourcepar = [img_flx,dx,dy],
           varypar   = [T,T,T],
           niter     = 10)
'''
reduced chi2=2.8237
I = 0.000532189 +/- 7.77046e-05
x = 0.0949058 +/- 0.0197072 arcsec
y = 0.0300131 +/- 0.018265 arcsec

'''


#### measure flux as gaussian
uvmodelfit(vis       = contvis,
           comptype  = 'G',
           sourcepar = [img_flx,dx,dy,0.5,0.5,0.0],
           varypar   = [T,T,T,T,T,T],
           niter     = 10)

'''
reduced chi2=2.82376
I = 0.000697503 +/- 0.000136838
x = 0.0714228 +/- 0.0233172 arcsec
y = 0.0173061 +/- 0.0217583 arcsec
a = 0.160133 +/- 0.0955601 arcsec <-- bad
r = 1 +/- 0.782773 <-- bad
p = -76.3977 +/- 57.2958 deg <-- bad
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
image noise     0.139 mJy
beam :  PA   -62    0.30 x  0.25 arcsec   #pix/beam:    95
center pixel    150,  150
Annulus is between 5.0 and 20.0 arcsec
Subtracting offset of 0.016 mJy/beam from image
RMS in annulus is 0.135

#   i  radius    #pix    total      rms     snr   beam_frac
#      (asec)            (mJy)     (mJy)
   0    0.000       1     0.47    0.131     3.6   0.011
   1    0.050       9     0.48    0.126     3.8   0.091
   2    0.100      37     0.50    0.123     4.1   0.323 <-- good
   3    0.150      81     0.52    0.143     3.6   0.571
   4    0.200     137     0.52    0.134     3.9   0.761 <-- good
   5    0.250     221     0.51    0.156     3.3   0.899
   6    0.300     317     0.47    0.223     2.1   0.961
   7    0.350     429     0.41    0.212     1.9   0.987
   8    0.400     553     0.36    0.251     1.4   0.996
   9    0.450     709     0.34    0.274     1.2   0.999
  10    0.500     877     0.33    0.290     1.1   1.000
  11    0.550    1049     0.34    0.317     1.1   1.000
  12    0.600    1257     0.37    0.393     0.9   1.000
  13    0.650    1481     0.36    0.351     1.0   1.000
  14    0.700    1709     0.30    0.400     0.7   1.000

PEAK FLUX: radius     total     rms    snr
              0.2      0.52    0.13    3.9
PEAK SNR:  radius     total     rms   snr
              0.1      0.50    0.12    4.1

'''



