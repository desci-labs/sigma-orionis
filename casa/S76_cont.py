#!/usr/bin/env


# ======================== Setup ===========================

field   = 76
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
# not resolved?



                        
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
xc,yc = 316,317
aper = 1.0


### check stats
img_max = imstat(imagename=contimg+'.image')['max'][0]
img_rms = imstat(imagename=contimg+'.image',region='annulus[['+str(xc)+'pix,'+str(yc)+'pix],['+str(80)+'pix,'+str(120)+'pix]]')['rms'][0]
bmaj = imhead(imagename=contimg+'.image', mode="get", hdkey="beammajor")
bmin = imhead(imagename=contimg+'.image', mode="get", hdkey="beamminor")
print '# Peak = {0:.2f} mJy, rms = {1:.2f} mJy, S/N = {2:.1f}'.format(1000*img_max, 1000*img_rms, img_max/img_rms)
print '# Beam = {0:.2f} x {1:.2f} arcsec'.format(bmaj.get('value'),bmin.get('value'))
# Peak = 8.59 mJy, rms = 0.18 mJy, S/N = 47.9
# Beam = 0.32 x 0.25 arcsec


### measure flux
img_rms = imstat(imagename=contimg+'.image',region='annulus[['+str(xc)+'pix,'+str(yc)+'pix],['+str(80)+'pix,'+str(120)+'pix]]')['rms'][0]
img_flx = imstat(imagename=contimg+'.image',region='circle[['+str(xc)+'pix,'+str(yc)+'pix],'+str(aper)+'arcsec]')['flux'][0]
print '# Flux = {0:.2f} mJy, rms = {1:.2f} mJy, S/N = {2:.1f}'.format(1000*img_flx, 1000*img_rms, img_flx/img_rms)
# Flux = 14.21 mJy, rms = 0.18 mJy, S/N = 79.2



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
reduced chi2=2.94955:
I = 0.0135896 +/- 0.000161026
x = 0.1304 +/- 0.00162508 arcsec
y = -0.0961022 +/- 0.00146854 arcsec
a = 0.196042 +/- 0.00503134 arcsec
r = 1 +/- 0.0398303
p = 62.4766 +/- 57.2958 deg
'''

#### measure flux as point source
uvmodelfit(vis       = contvis,
           comptype  = 'P',
           sourcepar = [img_flx,dx,dy],
           varypar   = [T,T,T],
           niter     = 10)

'''
I = 0.00959498 +/- 8.65311e-05
x = 0.126494 +/- 0.00131396 arcsec
y = -0.0934525 +/- 0.00115512 arcsec
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
image noise     0.314 mJy
beam :  PA   -61    0.32 x  0.25 arcsec   #pix/beam:   102
center pixel    150,  150
Annulus is between 5.0 and 20.0 arcsec
Subtracting offset of -0.020 mJy/beam from image
RMS in annulus is 0.189

#   i  radius    #pix    total      rms     snr   beam_frac
#      (asec)            (mJy)     (mJy)
   0    0.000       1     8.61    0.153    56.3   0.010
   1    0.050       9     8.72    0.179    48.6   0.085
   2    0.100      37     9.09    0.138    66.0   0.303
   3    0.150      81     9.67    0.176    54.8   0.543
   4    0.200     137    10.34    0.186    55.7   0.732
   5    0.250     221    11.25    0.219    51.4   0.877
   6    0.300     317    12.12    0.270    44.9   0.948
   7    0.350     429    12.89    0.267    48.4   0.981
   8    0.400     553    13.51    0.312    43.3   0.993
   9    0.450     709    14.02    0.368    38.1   0.998
  10    0.500     877    14.30    0.271    52.8   1.000
  11    0.550    1049    14.41    0.481    30.0   1.000
  12    0.600    1257    14.44    0.489    29.5   1.000
  13    0.650    1481    14.43    0.506    28.5   1.000
  14    0.700    1709    14.44    0.379    38.1   1.000

PEAK FLUX: radius     total     rms    snr
              0.7     14.44    0.38   38.1
PEAK SNR:  radius     total     rms   snr
              0.1      9.09    0.14   66.0
'''



