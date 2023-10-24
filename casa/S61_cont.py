#!/usr/bin/env


# ======================== Setup ===========================

field   = 61
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

### BINARY?


# ================== MEASURE FLUX WITHIN APERTURE ==================


### check image to get center
imview(raster   = [{'file':contimg+'.fits'}],
        contour = [{'file':contimg+'.fits'}])


### set aperture and source center
### large offset in dec, right star?
xc,yc = 321,322
aper = 1.0


### check stats
img_max = imstat(imagename=contimg+'.image')['max'][0]
img_rms = imstat(imagename=contimg+'.image',region='annulus[['+str(xc)+'pix,'+str(yc)+'pix],['+str(80)+'pix,'+str(120)+'pix]]')['rms'][0]
bmaj = imhead(imagename=contimg+'.image', mode="get", hdkey="beammajor")
bmin = imhead(imagename=contimg+'.image', mode="get", hdkey="beamminor")
print '# Peak = {0:.2f} mJy, rms = {1:.2f} mJy, S/N = {2:.1f}'.format(1000*img_max, 1000*img_rms, img_max/img_rms)
print '# Beam = {0:.2f} x {1:.2f} arcsec'.format(bmaj.get('value'),bmin.get('value'))
# Peak = 2.28 mJy, rms = 0.16 mJy, S/N = 14.6
# Beam = 0.32 x 0.25 arcsec

### measure flux
img_rms = imstat(imagename=contimg+'.image',region='annulus[['+str(xc)+'pix,'+str(yc)+'pix],['+str(80)+'pix,'+str(120)+'pix]]')['rms'][0]
img_flx = imstat(imagename=contimg+'.image',region='circle[['+str(xc)+'pix,'+str(yc)+'pix],'+str(aper)+'arcsec]')['flux'][0]
print '# Flux = {0:.2f} mJy, rms = {1:.2f} mJy, S/N = {2:.1f}'.format(1000*img_flx, 1000*img_rms, img_flx/img_rms)
# Flux = 3.34 mJy, rms = 0.16 mJy, S/N = 21.3


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
reduced chi2=2.82
I = 0.00318182 +/- 0.000149884
x = -0.0449079 +/- 0.00576564 arcsec
y = 0.0539021 +/- 0.00525231 arcsec
a = 0.158712 +/- 0.0251245 arcsec <-- bad
r = 1 +/- 0.196639
p = -50.5188 +/- 57.2958 deg <---bad
'''


#### measure flux as point source
uvmodelfit(vis       = contvis,
           comptype  = 'P',
           sourcepar = [img_flx,dx,dy],
           varypar   = [T,T,T],
           niter     = 10)

'''
I = 0.00248845 +/- 8.60239e-05
x = -0.0423129 +/- 0.00486782 arcsec
y = 0.051797 +/- 0.00437196 arcsec

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
image noise     0.170 mJy
beam :  PA   -61    0.32 x  0.25 arcsec   #pix/beam:    99
center pixel    150,  150
Annulus is between 5.0 and 20.0 arcsec
Subtracting offset of 0.006 mJy/beam from image
RMS in annulus is 0.171

#   i  radius    #pix    total      rms     snr   beam_frac
#      (asec)            (mJy)     (mJy)
   0    0.000       1     2.28    0.141    16.1   0.010
   1    0.050       9     2.30    0.173    13.3   0.087
   2    0.100      37     2.38    0.159    15.0   0.310
   3    0.150      81     2.50    0.175    14.3   0.553
   4    0.200     137     2.62    0.144    18.2   0.742
   5    0.250     221     2.76    0.176    15.7   0.885
   6    0.300     317     2.88    0.203    14.2   0.953
   7    0.350     429     2.97    0.240    12.4   0.983
   8    0.400     553     3.04    0.286    10.7   0.995
   9    0.450     709     3.11    0.328     9.5   0.999
  10    0.500     877     3.15    0.364     8.7   1.000
  11    0.550    1049     3.17    0.431     7.3   1.000
  12    0.600    1257     3.21    0.346     9.3   1.000
  13    0.650    1481     3.26    0.257    12.7   1.000

PEAK FLUX: radius     total     rms    snr
              0.7      3.26    0.26   12.7
PEAK SNR:  radius     total     rms   snr
              0.2      2.62    0.14   18.2
'''



