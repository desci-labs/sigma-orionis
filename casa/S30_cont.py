#!/usr/bin/env


# ======================== Setup ===========================

field   = 30
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
xc,yc = 315,317
aper = 1.0


### check stats
img_max = imstat(imagename=contimg+'.image')['max'][0]
img_rms = imstat(imagename=contimg+'.image',region='annulus[['+str(xc)+'pix,'+str(yc)+'pix],['+str(80)+'pix,'+str(120)+'pix]]')['rms'][0]
bmaj = imhead(imagename=contimg+'.image', mode="get", hdkey="beammajor")
bmin = imhead(imagename=contimg+'.image', mode="get", hdkey="beamminor")
print imhead(imagename=contimg+'.image', mode="get", hdkey="object")
print '# Peak = {0:.2f} mJy, rms = {1:.2f} mJy, S/N = {2:.1f}'.format(1000*img_max, 1000*img_rms, img_max/img_rms)
print '# Beam = {0:.2f} x {1:.2f} arcsec'.format(bmaj.get('value'),bmin.get('value'))
# Peak = 0.67 mJy, rms = 0.15 mJy, S/N = 4.5
# Beam = 0.31 x 0.25 arcsec

### measure flux
img_rms = imstat(imagename=contimg+'.image',region='annulus[['+str(xc)+'pix,'+str(yc)+'pix],['+str(80)+'pix,'+str(120)+'pix]]')['rms'][0]
img_flx = imstat(imagename=contimg+'.image',region='circle[['+str(xc)+'pix,'+str(yc)+'pix],'+str(aper)+'arcsec]')['flux'][0]
print '# Flux = {0:.2f} mJy, rms = {1:.2f} mJy, S/N = {2:.1f}'.format(1000*img_flx, 1000*img_rms, img_flx/img_rms)
# Flux = 1.41 mJy, rms = 0.15 mJy, S/N = 9.5


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
I = 0.00125083 +/- 0.000175662
x = 0.128864 +/- 0.0221255 arcsec
y = -0.0608445 +/- 0.0206058 arcsec
a = 0.268483 +/- 0.0585279 arcsec <-- bad
r = 1 +/- 0.337656
p = 15.6453 +/- 57.2958 deg <-- bad

'''


#### measure flux as point source
uvmodelfit(vis       = contvis,
           comptype  = 'P',
           sourcepar = [img_flx,dx,dy],
           varypar   = [T,T,T],
           niter     = 10)
'''
I = 0.000711 +/- 7.95105e-05
x = 0.152952 +/- 0.0154072 arcsec
y = -0.0971363 +/- 0.0140592 arcsec

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
image noise     0.153 mJy
beam :  PA   -62    0.31 x  0.25 arcsec   #pix/beam:    97
center pixel    150,  150
Annulus is between 5.0 and 20.0 arcsec
Subtracting offset of 0.011 mJy/beam from image
RMS in annulus is 0.154

#   i  radius    #pix    total      rms     snr   beam_frac
#      (asec)            (mJy)     (mJy)
   0    0.000       1     0.65    0.141     4.6   0.010
   1    0.050       9     0.66    0.146     4.5   0.089
   2    0.100      37     0.68    0.130     5.2   0.317
   3    0.150      81     0.72    0.113     6.3   0.563 <-- good
   4    0.200     137     0.76    0.165     4.6   0.753
   5    0.250     221     0.82    0.222     3.7   0.893
   6    0.300     317     0.85    0.167     5.1   0.958
   7    0.350     429     0.87    0.206     4.2   0.986
   8    0.400     553     0.86    0.227     3.8   0.996
   9    0.450     709     0.83    0.308     2.7   0.999
  10    0.500     877     0.83    0.316     2.6   1.000
  11    0.550    1049     0.85    0.381     2.2   1.000
  12    0.600    1257     0.92    0.436     2.1   1.000
  13    0.650    1481     1.01    0.512     2.0   1.000
  14    0.700    1709     1.10    0.467     2.4   1.000

PEAK FLUX: radius     total     rms    snr
             0.70      1.10    0.47    2.4
PEAK SNR:  radius     total     rms   snr
             0.15      0.72    0.11    6.3

'''



