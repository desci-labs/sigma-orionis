#!/usr/bin/env


# ======================== Setup ===========================

field   = 67 
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
xc,yc = 315,318
aper = 1.0


### check stats
img_max = imstat(imagename=contimg+'.image')['max'][0]
img_rms = imstat(imagename=contimg+'.image',region='annulus[['+str(xc)+'pix,'+str(yc)+'pix],['+str(80)+'pix,'+str(120)+'pix]]')['rms'][0]
bmaj = imhead(imagename=contimg+'.image', mode="get", hdkey="beammajor")
bmin = imhead(imagename=contimg+'.image', mode="get", hdkey="beamminor")
print imhead(imagename=contimg+'.image', mode="get", hdkey="object")
print '# Peak = {0:.2f} mJy, rms = {1:.2f} mJy, S/N = {2:.1f}'.format(1000*img_max, 1000*img_rms, img_max/img_rms)
print '# Beam = {0:.2f} x {1:.2f} arcsec'.format(bmaj.get('value'),bmin.get('value'))
# Peak = 1.40 mJy, rms = 0.15 mJy, S/N = 9.1
# Beam = 0.32 x 0.25 arcsec


### measure flux
img_rms = imstat(imagename=contimg+'.image',region='annulus[['+str(xc)+'pix,'+str(yc)+'pix],['+str(80)+'pix,'+str(120)+'pix]]')['rms'][0]
img_flx = imstat(imagename=contimg+'.image',region='circle[['+str(xc)+'pix,'+str(yc)+'pix],'+str(aper)+'arcsec]')['flux'][0]
print '# Flux = {0:.2f} mJy, rms = {1:.2f} mJy, S/N = {2:.1f}'.format(1000*img_flx, 1000*img_rms, img_flx/img_rms)
# Flux = 1.18 mJy, rms = 0.15 mJy, S/N = 7.7



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
reduced chi2=2.86
I = 0.00139914 +/- 0.000126686
x = 0.13783 +/- 0.00905661 arcsec
y = -0.0561021 +/- 0.00807948 arcsec
a = 2.47696e-07 +/- 17456.3 arcsec <--bad
r = 1 +/- 1.0404e+11 <--bad
p = -8.70779 +/- 57.2958 deg <--bad
'''


#### measure flux as point source
uvmodelfit(vis       = contvis,
           comptype  = 'P',
           sourcepar = [img_flx,dx,dy],
           varypar   = [T,T,T],
           niter     = 10)
'''
I = 0.00140819 +/- 8.94273e-05
x = 0.137843 +/- 0.00899837 arcsec
y = -0.0561064 +/- 0.00802753 arcsec
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
image noise     0.164 mJy
beam :  PA   -62    0.32 x  0.25 arcsec   #pix/beam:   100
center pixel    150,  150
Annulus is between 5.0 and 20.0 arcsec
Subtracting offset of -0.005 mJy/beam from image
RMS in annulus is 0.151

#   i  radius    #pix    total      rms     snr   beam_frac
#      (asec)            (mJy)     (mJy)
   0    0.000       1     1.40    0.156     9.0   0.010
   1    0.050       9     1.40    0.179     7.8   0.086
   2    0.100      37     1.39    0.159     8.8   0.307
   3    0.150      81     1.39    0.154     9.0   0.548
   4    0.200     137     1.38    0.191     7.2   0.738
   5    0.250     221     1.39    0.215     6.4   0.881
   6    0.300     317     1.41    0.216     6.5   0.951
   7    0.350     429     1.46    0.248     5.9   0.982
   8    0.400     553     1.52    0.265     5.7   0.994
   9    0.450     709     1.62    0.267     6.1   0.998
  10    0.500     877     1.70    0.273     6.2   1.000
  11    0.550    1049     1.75    0.367     4.8   1.000
  12    0.600    1257     1.79    0.388     4.6   1.000
  13    0.650    1481     1.79    0.518     3.5   1.000
  14    0.700    1709     1.80    0.630     2.9   1.000
  15    0.750    1961     1.80    0.362     5.0   1.000
  16    0.800    2233     1.81    0.564     3.2   1.000
  17    0.850    2537     1.77    0.549     3.2   1.000
  18    0.900    2821     1.68    0.561     3.0   1.000
  19    0.950    3149     1.55    0.545     2.8   1.000

PEAK FLUX: radius     total     rms    snr
              0.8      1.81    0.56    3.2
PEAK SNR:  radius     total     rms   snr
              0.2      1.39    0.15    9.0

'''



