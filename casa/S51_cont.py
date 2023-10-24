#!/usr/bin/env


# ======================== Setup ===========================

field   = 51
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
xc,yc = 306,306
aper = 1.0


### check stats
img_max = imstat(imagename=contimg+'.image')['max'][0]
img_rms = imstat(imagename=contimg+'.image',region='annulus[['+str(xc)+'pix,'+str(yc)+'pix],['+str(80)+'pix,'+str(120)+'pix]]')['rms'][0]
bmaj = imhead(imagename=contimg+'.image', mode="get", hdkey="beammajor")
bmin = imhead(imagename=contimg+'.image', mode="get", hdkey="beamminor")
print imhead(imagename=contimg+'.image', mode="get", hdkey="object")
print '# Peak = {0:.2f} mJy, rms = {1:.2f} mJy, S/N = {2:.1f}'.format(1000*img_max, 1000*img_rms, img_max/img_rms)
print '# Beam = {0:.2f} x {1:.2f} arcsec'.format(bmaj.get('value'),bmin.get('value'))
# Peak = 0.65 mJy, rms = 0.15 mJy, S/N = 4.3
# Beam = 0.30 x 0.25 arcsec


### measure flux
img_rms = imstat(imagename=contimg+'.image',region='annulus[['+str(xc)+'pix,'+str(yc)+'pix],['+str(80)+'pix,'+str(120)+'pix]]')['rms'][0]
img_flx = imstat(imagename=contimg+'.image',region='circle[['+str(xc)+'pix,'+str(yc)+'pix],'+str(aper)+'arcsec]')['flux'][0]
print '# Flux = {0:.2f} mJy, rms = {1:.2f} mJy, S/N = {2:.1f}'.format(1000*img_flx, 1000*img_rms, img_flx/img_rms)
# Flux = -0.09 mJy, rms = 0.15 mJy, S/N = -0.6


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
I = 0.000716489 +/- 0.000137354
x = 0.0753879 +/- 0.0228115 arcsec
y = 0.0219298 +/- 0.021288 arcsec
a = 0.161232 +/- 0.0953772 arcsec <--bad
r = 1 +/- 0.752462
p = -52.5902 +/- 57.2958 deg <-- bad


'''


#### measure flux as point source
uvmodelfit(vis       = contvis,
           comptype  = 'P',
           sourcepar = [img_flx,dx,dy],
           varypar   = [T,T,T],
           niter     = 10)
'''
I = 0.000519734 +/- 8.40338e-05
x = 0.425928 +/- 0.0225313 arcsec
y = -0.430485 +/- 0.0203055 arcsec
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
image noise     0.153 mJy
beam :  PA   -61    0.31 x  0.25 arcsec   #pix/beam:    98
center pixel    150,  150
Annulus is between 5.0 and 20.0 arcsec
Subtracting offset of -0.015 mJy/beam from image
RMS in annulus is 0.16

#   i  radius    #pix    total      rms     snr   beam_frac
#      (asec)            (mJy)     (mJy)
   0    0.000       1     0.61    0.133     4.6   0.010
   1    0.050       9     0.60    0.136     4.4   0.088
   2    0.100      37     0.58    0.112     5.2   0.314 <-- good
   3    0.150      81     0.55    0.156     3.5   0.558
   4    0.200     137     0.53    0.150     3.5   0.748
   5    0.250     221     0.53    0.175     3.0   0.889
   6    0.300     317     0.55    0.199     2.8   0.956
   7    0.350     429     0.60    0.227     2.6   0.985
   8    0.400     553     0.63    0.243     2.6   0.995
   9    0.450     709     0.65    0.351     1.8   0.999
  10    0.500     877     0.63    0.301     2.1   1.000
  11    0.550    1049     0.58    0.336     1.7   1.000
  12    0.600    1257     0.55    0.296     1.9   1.000
  13    0.650    1481     0.52    0.393     1.3   1.000
  14    0.700    1709     0.52    0.489     1.1   1.000

PEAK FLUX: radius     total     rms    snr
             0.45      0.65    0.35    1.8
PEAK SNR:  radius     total     rms   snr
             0.10      0.58    0.11    5.2
'''
