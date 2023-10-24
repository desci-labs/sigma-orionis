#!/usr/bin/env


# ======================== Setup ===========================

field   = 32
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



# ================== MEASURE FLUX WITHIN APERTURE ==================


### check image to get center
imview(raster   = [{'file':contimg+'.fits'}],
        contour = [{'file':contimg+'.fits'}])


### set aperture and source center
xc,yc = 314,320
aper = 1.0


### check stats
img_max = imstat(imagename=contimg+'.image')['max'][0]
img_rms = imstat(imagename=contimg+'.image',region='annulus[['+str(xc)+'pix,'+str(yc)+'pix],['+str(80)+'pix,'+str(120)+'pix]]')['rms'][0]
bmaj = imhead(imagename=contimg+'.image', mode="get", hdkey="beammajor")
bmin = imhead(imagename=contimg+'.image', mode="get", hdkey="beamminor")
print '# Peak = {0:.2f} mJy, rms = {1:.2f} mJy, S/N = {2:.1f}'.format(1000*img_max, 1000*img_rms, img_max/img_rms)
print '# Beam = {0:.2f} x {1:.2f} arcsec'.format(bmaj.get('value'),bmin.get('value'))
# Peak = 1.88 mJy, rms = 0.14 mJy, S/N = 13.5
# Beam = 0.31 x 0.25 arcsec


### measure flux
img_rms = imstat(imagename=contimg+'.image',region='annulus[['+str(xc)+'pix,'+str(yc)+'pix],['+str(80)+'pix,'+str(120)+'pix]]')['rms'][0]
img_flx = imstat(imagename=contimg+'.image',region='circle[['+str(xc)+'pix,'+str(yc)+'pix],'+str(aper)+'arcsec]')['flux'][0]
print '# Flux = {0:.2f} mJy, rms = {1:.2f} mJy, S/N = {2:.1f}'.format(1000*img_flx, 1000*img_rms, img_flx/img_rms)
# Flux = 2.28 mJy, rms = 0.14 mJy, S/N = 16.4


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
reduced chi2=2.949
I = 0.00193518 +/- 0.000119158
x = 0.171186 +/- 0.00639525 arcsec
y = 0.0156325 +/- 0.00549957 arcsec
a = 0.0949431 +/- 0.0345049 arcsec <--bad
r = 1e-08 +/- 3.58135e+07 <--- bad
p = 68.0728 +/- 29.0083 deg

'''



#### measure flux as point source
uvmodelfit(vis       = contvis,
           comptype  = 'P',
           sourcepar = [img_flx,dx,dy],
           varypar   = [T,T,T],
           niter     = 10)

'''
I = 0.00190249 +/- 7.97899e-05
x = 0.171043 +/- 0.00578967 arcsec
y = 0.015111 +/- 0.0052669 arcsec
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
image noise     0.150 mJy
beam :  PA   -62    0.31 x  0.25 arcsec   #pix/beam:    96
center pixel    150,  150
Annulus is between 5.0 and 20.0 arcsec
Subtracting offset of -0.007 mJy/beam from image
RMS in annulus is 0.155

#   i  radius    #pix    total      rms     snr   beam_frac
#      (asec)            (mJy)     (mJy)
   0    0.000       1     1.89    0.151    12.5   0.010
   1    0.050       9     1.88    0.142    13.3   0.089
   2    0.100      37     1.88    0.160    11.7   0.318
   3    0.150      81     1.86    0.107    17.4   0.564 <-- good
   4    0.200     137     1.82    0.143    12.8   0.754
   5    0.250     221     1.76    0.145    12.1   0.894
   6    0.300     317     1.69    0.182     9.3   0.958
   7    0.350     429     1.62    0.193     8.4   0.986
   8    0.400     553     1.60    0.202     7.9   0.996
   9    0.450     709     1.63    0.266     6.1   0.999

PEAK FLUX: radius     total     rms    snr
             0.00      1.89    0.15   12.5
PEAK SNR:  radius     total     rms   snr
             0.15      1.86    0.11   17.4
'''



