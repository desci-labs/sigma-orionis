#!/usr/bin/env


# ======================== Setup ===========================

field   = 65
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
xc,yc = 321,316
aper = 1.0


### check stats
img_max = imstat(imagename=contimg+'.image')['max'][0]
img_rms = imstat(imagename=contimg+'.image',region='annulus[['+str(xc)+'pix,'+str(yc)+'pix],['+str(80)+'pix,'+str(120)+'pix]]')['rms'][0]
bmaj = imhead(imagename=contimg+'.image', mode="get", hdkey="beammajor")
bmin = imhead(imagename=contimg+'.image', mode="get", hdkey="beamminor")
print imhead(imagename=contimg+'.image', mode="get", hdkey="object")
print '# Peak = {0:.2f} mJy, rms = {1:.2f} mJy, S/N = {2:.1f}'.format(1000*img_max, 1000*img_rms, img_max/img_rms)
print '# Beam = {0:.2f} x {1:.2f} arcsec'.format(bmaj.get('value'),bmin.get('value'))
# Peak = 1.68 mJy, rms = 0.15 mJy, S/N = 10.9
# Beam = 0.32 x 0.25 arcsec


### measure flux
img_rms = imstat(imagename=contimg+'.image',region='annulus[['+str(xc)+'pix,'+str(yc)+'pix],['+str(80)+'pix,'+str(120)+'pix]]')['rms'][0]
img_flx = imstat(imagename=contimg+'.image',region='circle[['+str(xc)+'pix,'+str(yc)+'pix],'+str(aper)+'arcsec]')['flux'][0]
print '# Flux = {0:.2f} mJy, rms = {1:.2f} mJy, S/N = {2:.1f}'.format(1000*img_flx, 1000*img_rms, img_flx/img_rms)
# Flux = 2.38 mJy, rms = 0.15 mJy, S/N = 15.5


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
           niter     = 5)

'''
I = 0.00145928 +/- 0.000121787
x = -0.0497924 +/- 0.00838384 arcsec
y = -0.110417 +/- 0.00744949 arcsec
a = 0.000184914 +/- 21.8111 arcsec
r = 6.89519e-06 +/- 1.90145e+10
p = -3.39715 +/- 9.99102e+06 deg

'''


#### measure flux as point source
uvmodelfit(vis       = contvis,
           comptype  = 'P',
           sourcepar = [img_flx,dx,dy],
           varypar   = [T,T,T],
           niter     = 10)
'''
I = 0.00171309 +/- 8.54171e-05
x = -0.0478716 +/- 0.0071417 arcsec
y = -0.109489 +/- 0.00634578 arcsec
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
image noise     0.163 mJy
beam :  PA   -61    0.32 x  0.25 arcsec   #pix/beam:   101
center pixel    150,  150
Annulus is between 5.0 and 20.0 arcsec
Subtracting offset of -0.004 mJy/beam from image
RMS in annulus is 0.163

#   i  radius    #pix    total      rms     snr   beam_frac
#      (asec)            (mJy)     (mJy)
   0    0.000       1     1.68    0.140    12.0   0.010
   1    0.050       9     1.67    0.136    12.2   0.085
   2    0.100      37     1.63    0.142    11.4   0.306
   3    0.150      81     1.57    0.184     8.6   0.546
   4    0.200     137     1.51    0.181     8.4   0.736
   5    0.250     221     1.47    0.174     8.4   0.880
   6    0.300     317     1.45    0.221     6.6   0.950
   7    0.350     429     1.48    0.233     6.4   0.982
   8    0.400     553     1.54    0.284     5.4   0.994
   9    0.450     709     1.64    0.332     4.9   0.998

PEAK FLUX: radius     total     rms    snr
              0.0      1.68    0.14   12.0
PEAK SNR:  radius     total     rms   snr
              0.1      1.67    0.14   12.2

'''



