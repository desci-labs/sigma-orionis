#!/usr/bin/env


# ======================== Setup ===========================

field   = 16 
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
# resolved



                        
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
xc,yc = 318,320
aper = 1.0


### check stats
img_max = imstat(imagename=contimg+'.image')['max'][0]
img_rms = imstat(imagename=contimg+'.image',region='annulus[['+str(xc)+'pix,'+str(yc)+'pix],['+str(80)+'pix,'+str(120)+'pix]]')['rms'][0]
bmaj = imhead(imagename=contimg+'.image', mode="get", hdkey="beammajor")
bmin = imhead(imagename=contimg+'.image', mode="get", hdkey="beamminor")
print imhead(imagename=contimg+'.image', mode="get", hdkey="object")
print '# Peak = {0:.2f} mJy, rms = {1:.2f} mJy, S/N = {2:.1f}'.format(1000*img_max, 1000*img_rms, img_max/img_rms)
print '# Beam = {0:.2f} x {1:.2f} arcsec'.format(bmaj.get('value'),bmin.get('value'))
# Peak = 5.06 mJy, rms = 0.15 mJy, S/N = 33.5
# Beam = 0.30 x 0.25 arcsec


### measure flux
img_rms = imstat(imagename=contimg+'.image',region='annulus[['+str(xc)+'pix,'+str(yc)+'pix],['+str(80)+'pix,'+str(120)+'pix]]')['rms'][0]
img_flx = imstat(imagename=contimg+'.image',region='circle[['+str(xc)+'pix,'+str(yc)+'pix],'+str(aper)+'arcsec]')['flux'][0]
print '# Flux = {0:.2f} mJy, rms = {1:.2f} mJy, S/N = {2:.1f}'.format(1000*img_flx, 1000*img_rms, img_flx/img_rms)
# Flux = 6.98 mJy, rms = 0.15 mJy, S/N = 46.2


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
I = 0.00567564 +/- 0.00012283
x = 0.0542803 +/- 0.00216842 arcsec
y = 0.00376441 +/- 0.00217767 arcsec
a = 0.117163 +/- 0.011083 arcsec <-- good
r = 0.568682 +/- 0.15627
p = -5.41797 +/- 10.5141 deg

'''

### measure as point source (since unresolved)
uvmodelfit(vis       = contvis,
           comptype  = 'P',
           sourcepar = [img_flx,dx,dy],
           varypar   = [T,T,T],
           niter     = 10)
'''
I = 0.00515518 +/- 7.9091e-05
x = 0.0540966 +/- 0.00208561 arcsec
y = 0.00440632 +/- 0.00193195 arcsec
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
image noise     0.195 mJy
beam :  PA   -62    0.30 x  0.25 arcsec   #pix/beam:    96
center pixel    150,  150
Annulus is between 5.0 and 20.0 arcsec
Subtracting offset of 0.006 mJy/beam from image
RMS in annulus is 0.154

#   i  radius    #pix    total      rms     snr   beam_frac
#      (asec)            (mJy)     (mJy)
   0    0.000       1     5.05    0.151    33.4   0.010
   1    0.050       9     5.07    0.148    34.3   0.089
   2    0.100      37     5.14    0.134    38.4   0.318
   3    0.150      81     5.24    0.156    33.6   0.565
   4    0.200     137     5.34    0.160    33.5   0.755
   5    0.250     221     5.45    0.167    32.6   0.895
   6    0.300     317     5.57    0.195    28.5   0.959 <-- matches image, matches UV
   7    0.350     429     5.67    0.201    28.2   0.986
   8    0.400     553     5.79    0.278    20.8   0.996
   9    0.450     709     5.95    0.308    19.3   0.999
  10    0.500     877     6.12    0.250    24.5   1.000
  11    0.550    1049     6.30    0.367    17.2   1.000
  12    0.600    1257     6.52    0.424    15.4   1.000
  13    0.650    1481     6.72    0.292    23.0   1.000
  14    0.700    1709     6.83    0.442    15.5   1.000 <-- high
  15    0.750    1961     6.83    0.279    24.4   1.000
  16    0.800    2233     6.77    0.443    15.3   1.000
  17    0.850    2537     6.70    0.332    20.2   1.000
  18    0.900    2821     6.69    0.244    27.4   1.000
  19    0.950    3149     6.72    0.401    16.8   1.000

PEAK FLUX: radius     total     rms    snr
             0.70      6.83    0.44   15.5
PEAK SNR:  radius     total     rms   snr
             0.10      5.14    0.13   38.4

'''



