#!/usr/bin/env


# ======================== Setup ===========================

field   = 56
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
### large offset in dec, right star?
xc,yc = 322,318
aper = 1.0


### check stats
img_max = imstat(imagename=contimg+'.image')['max'][0]
img_rms = imstat(imagename=contimg+'.image',region='annulus[['+str(xc)+'pix,'+str(yc)+'pix],['+str(80)+'pix,'+str(120)+'pix]]')['rms'][0]
bmaj = imhead(imagename=contimg+'.image', mode="get", hdkey="beammajor")
bmin = imhead(imagename=contimg+'.image', mode="get", hdkey="beamminor")
print '# Peak = {0:.2f} mJy, rms = {1:.2f} mJy, S/N = {2:.1f}'.format(1000*img_max, 1000*img_rms, img_max/img_rms)
print '# Beam = {0:.2f} x {1:.2f} arcsec'.format(bmaj.get('value'),bmin.get('value'))
# Peak = 1.80 mJy, rms = 0.17 mJy, S/N = 10.6
# Beam = 0.31 x 0.25 arcsec

### measure flux
img_rms = imstat(imagename=contimg+'.image',region='annulus[['+str(xc)+'pix,'+str(yc)+'pix],['+str(80)+'pix,'+str(120)+'pix]]')['rms'][0]
img_flx = imstat(imagename=contimg+'.image',region='circle[['+str(xc)+'pix,'+str(yc)+'pix],'+str(aper)+'arcsec]')['flux'][0]
print '# Flux = {0:.2f} mJy, rms = {1:.2f} mJy, S/N = {2:.1f}'.format(1000*img_flx, 1000*img_rms, img_flx/img_rms)
# Flux = 1.84 mJy, rms = 0.17 mJy, S/N = 10.8


# ======================== Measure flux with UVMODELFIT ==================


### calculate offset from phase center in arcsec
pixscale = 0.03             # must match 'cell'                 
dx = pixscale*(320.0-xc)    # offset to east (left)
dy = pixscale*(yc-320.0)    # offset to north (up)

  
#### measure flux as gaussian (not resolved with visbin)
uvmodelfit(vis       = contvis,
           comptype  = 'G',
           sourcepar = [img_flx,dx,dy,0.5,0.5,0.0],
           varypar   = [T,T,T,T,T,T],
           niter     = 10)

'''
reduced chi2=2.80
I = 0.00250435 +/- 0.000157932
x = -0.0469227 +/- 0.00668091 arcsec
y = -0.0540746 +/- 0.00761597 arcsec
a = 0.200028 +/- 0.0254074 arcsec <--meh ok
r = 0.415858 +/- 0.239062
p = 8.08264 +/- 10.933 deg
'''


#### measure flux as point source
uvmodelfit(vis       = contvis,
           comptype  = 'P',
           sourcepar = [img_flx,dx,dy],
           varypar   = [T,T,T],
           niter     = 10)
'''
I = 0.00196571 +/- 9.10857e-05
x = -0.0481259 +/- 0.00645832 arcsec
y = -0.0563281 +/- 0.00582532 arcsec
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
image noise     0.173 mJy
beam :  PA   -61    0.31 x  0.25 arcsec   #pix/beam:    99
center pixel    150,  150
Annulus is between 5.0 and 20.0 arcsec
Subtracting offset of -0.013 mJy/beam from image
RMS in annulus is 0.166

#   i  radius    #pix    total      rms     snr   beam_frac
#      (asec)            (mJy)     (mJy)
   0    0.000       1     1.81    0.161    11.3   0.010
   1    0.100      37     1.92    0.164    11.7   0.311
   2    0.200     137     2.13    0.163    13.1   0.745
   3    0.300     317     2.28    0.248     9.2   0.954
   4    0.400     553     2.30    0.267     8.6   0.995 <-- goood
   5    0.500     877     2.24    0.400     5.6   1.000
   6    0.600    1257     2.03    0.399     5.1   1.000
   7    0.700    1709     2.03    0.492     4.1   1.000
   8    0.800    2233     2.34    0.365     6.4   1.000
   9    0.900    2821     2.40    0.563     4.3   1.000

PEAK FLUX: radius     total     rms    snr
              0.9      2.40    0.56    4.3
PEAK SNR:  radius     total     rms   snr
              0.2      2.13    0.16   13.1
'''



