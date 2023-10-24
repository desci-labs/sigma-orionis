#!/usr/bin/env


# ======================== Setup ===========================

field   = 90
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
xc,yc = 321,319
aper = 1.0


### check stats
img_max = imstat(imagename=contimg+'.image')['max'][0]
img_rms = imstat(imagename=contimg+'.image',region='annulus[['+str(xc)+'pix,'+str(yc)+'pix],['+str(80)+'pix,'+str(120)+'pix]]')['rms'][0]
bmaj = imhead(imagename=contimg+'.image', mode="get", hdkey="beammajor")
bmin = imhead(imagename=contimg+'.image', mode="get", hdkey="beamminor")
print imhead(imagename=contimg+'.image', mode="get", hdkey="object")
print '# Peak = {0:.2f} mJy, rms = {1:.2f} mJy, S/N = {2:.1f}'.format(1000*img_max, 1000*img_rms, img_max/img_rms)
print '# Beam = {0:.2f} x {1:.2f} arcsec'.format(bmaj.get('value'),bmin.get('value'))
# Peak = 1.64 mJy, rms = 0.16 mJy, S/N = 9.9
# Beam = 0.33 x 0.25 arcsec


### measure flux
img_rms = imstat(imagename=contimg+'.image',region='annulus[['+str(xc)+'pix,'+str(yc)+'pix],['+str(80)+'pix,'+str(120)+'pix]]')['rms'][0]
img_flx = imstat(imagename=contimg+'.image',region='circle[['+str(xc)+'pix,'+str(yc)+'pix],'+str(aper)+'arcsec]')['flux'][0]
print '# Flux = {0:.2f} mJy, rms = {1:.2f} mJy, S/N = {2:.1f}'.format(1000*img_flx, 1000*img_rms, img_flx/img_rms)
# Flux = 2.48 mJy, rms = 0.16 mJy, S/N = 15.0


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
I = 0.00226803 +/- 0.000177172
x = -0.0159726 +/- 0.0112931 arcsec
y = -0.013262 +/- 0.00963527 arcsec
a = 0.248182 +/- 0.040459 arcsec <-- bad
r = 0.649011 +/- 0.163167
p = -52.6514 +/- 16.3033 deg
'''


#### measure flux as point source
uvmodelfit(vis       = contvis,
           comptype  = 'P',
           sourcepar = [img_flx,dx,dy],
           varypar   = [T,T,T],
           niter     = 10)
'''
I = 0.00162877 +/- 9.42577e-05
x = -0.0182142 +/- 0.00849569 arcsec
y = -0.0180444 +/- 0.00735314 arcsec
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

'''
pixel scale     0.03 arcsec
image noise     0.179 mJy
beam :  PA   -62    0.33 x  0.25 arcsec   #pix/beam:   103
center pixel    150,  150
Annulus is between 5.0 and 20.0 arcsec
Subtracting offset of -0.015 mJy/beam from image
RMS in annulus is 0.176

#   i  radius    #pix    total      rms     snr   beam_frac
#      (asec)            (mJy)     (mJy)
   0    0.000       1     1.65    0.173     9.5   0.010
   1    0.050       9     1.67    0.174     9.6   0.084
   2    0.100      37     1.72    0.191     9.0   0.302
   3    0.150      81     1.82    0.116    15.7   0.540
   4    0.200     137     1.93    0.189    10.2   0.729
   5    0.250     221     2.10    0.171    12.3   0.874
   6    0.300     317     2.26    0.212    10.7   0.947
   7    0.350     429     2.41    0.257     9.4   0.980
   8    0.400     553     2.52    0.307     8.2   0.993
   9    0.450     709     2.59    0.323     8.0   0.998
  10    0.500     877     2.64    0.282     9.3   1.000
  11    0.550    1049     2.68    0.351     7.7   1.000
  12    0.600    1257     2.70    0.410     6.6   1.000
  13    0.650    1481     2.74    0.407     6.7   1.000
  14    0.700    1709     2.76    0.367     7.5   1.000
  15    0.750    1961     2.75    0.419     6.6   1.000
  16    0.800    2233     2.74    0.344     8.0   1.000
  17    0.850    2537     2.78    0.291     9.5   1.000
  18    0.900    2821     2.84    0.474     6.0   1.000
  19    0.950    3149     2.95    0.250    11.8   1.000

PEAK FLUX: radius     total     rms    snr
             0.95      2.95    0.25   11.8
PEAK SNR:  radius     total     rms   snr
             0.15      1.82    0.12   15.7

'''



