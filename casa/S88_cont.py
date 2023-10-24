#!/usr/bin/env


# ======================== Setup ===========================

field   = 88 
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
# source is resolved



                        
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


# set aperture and source center
xc,yc = 290,314
aper = 1.0


### check stats
img_max = imstat(imagename=contimg+'.image')['max'][0]
img_rms = imstat(imagename=contimg+'.image',region='annulus[['+str(xc)+'pix,'+str(yc)+'pix],['+str(80)+'pix,'+str(120)+'pix]]')['rms'][0]
bmaj = imhead(imagename=contimg+'.image', mode="get", hdkey="beammajor")
bmin = imhead(imagename=contimg+'.image', mode="get", hdkey="beamminor")
print '# Peak = {0:.2f} mJy, rms = {1:.2f} mJy, S/N = {2:.1f}'.format(1000*img_max, 1000*img_rms, img_max/img_rms)
print '# Beam = {0:.2f} x {1:.2f} arcsec'.format(bmaj.get('value'),bmin.get('value'))
# Peak = 4.87 mJy, rms = 0.16 mJy, S/N = 29.7
# Beam = 0.33 x 0.25 arcsec


### measure flux
img_rms = imstat(imagename=contimg+'.image',region='annulus[['+str(xc)+'pix,'+str(yc)+'pix],['+str(80)+'pix,'+str(120)+'pix]]')['rms'][0]
img_flx = imstat(imagename=contimg+'.image',region='circle[['+str(xc)+'pix,'+str(yc)+'pix],'+str(aper)+'arcsec]')['flux'][0]
print '# Flux = {0:.2f} mJy, rms = {1:.2f} mJy, S/N = {2:.1f}'.format(1000*img_flx, 1000*img_rms, img_flx/img_rms)
# Flux = 15.09 mJy, rms = 0.16 mJy, S/N = 92.0




# ======================== Measure flux with UVMODELFIT ==================


### calculate offset from phase center in arcsec
pixscale = 0.03             # must match 'cell'                 
dx = pixscale*(320.0-xc)    # offset to east (left)
dy = pixscale*(yc-320.0)    # offset to north (up)

  
#### measure flux as gaussian (clearly resolved)
uvmodelfit(vis       = contvis,
           comptype  = 'G',
           sourcepar = [img_flx,dx,dy,0.5,0.5,0.0],
           varypar   = [T,T,T,T,T,T],
           niter     = 10)

'''
reduced chi2=2.8845
I = 0.0153789 +/- 0.00025158
x = 0.894028 +/- 0.00373495 arcsec
y = -0.168394 +/- 0.00362202 arcsec
a = 0.45679 +/- 0.0102281 arcsec <-- good
r = 0.845369 +/- 0.0272296
p = -33.6203 +/- 5.54293 deg

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



### use measure.py to get COG flux
'''
pixel scale     0.03 arcsec
image noise     0.260 mJy
beam :  PA   -61    0.33 x  0.25 arcsec   #pix/beam:   102
center pixel    150,  150
Annulus is between 5.0 and 20.0 arcsec
Subtracting offset of -0.020 mJy/beam from image
RMS in annulus is 0.172

#   i  radius    #pix    total      rms     snr   beam_frac
#      (asec)            (mJy)     (mJy)
   0    0.000       1     4.88    0.150    32.7   0.010
   1    0.100      37     5.42    0.155    35.0   0.302
   2    0.200     137     6.98    0.163    42.8   0.730
   3    0.300     317     9.68    0.211    45.8   0.947
   4    0.400     553    12.19    0.280    43.5   0.993
   5    0.500     877    13.85    0.314    44.2   1.000
   6    0.600    1257    14.65    0.497    29.5   1.000
   7    0.700    1709    15.12    0.431    35.1   1.000
   8    0.800    2233    15.31    0.505    30.3   1.000
   9    0.900    2821    15.47    0.718    21.5   1.000

PEAK FLUX: radius     total     rms    snr
              0.9     15.47    0.72   21.5 <-- good
PEAK SNR:  radius     total     rms   snr
              0.3      9.68    0.21   45.8

'''


